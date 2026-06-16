/* Round-1 — per-round process aliases for NF 25.x uniqueness constraint */
include { MERGE_HYBRID        as MERGE_HYBRID_R1          } from '../modules/deepvariant'
include { PEPPER_MARGIN_DV    as PEPPER_MARGIN_DV_R1      } from '../modules/deepvariant'
include { DV_MERGE_CHR_VCFS   as DV_MERGE_CHR_VCFS_R1    } from '../modules/deepvariant'
include { DV_MAKE_EXAMPLES    as DV_MAKE_EXAMPLES_HYB_R1  } from '../modules/deepvariant'
include { DV_CALL_VARIANTS    as DV_CALL_VARIANTS_HYB_R1  } from '../modules/deepvariant'
include { DV_POSTPROCESS      as DV_POSTPROCESS_HYB_R1    } from '../modules/deepvariant'
include { DV_MAKE_EXAMPLES    as DV_MAKE_EXAMPLES_ONT_R1  } from '../modules/deepvariant'
include { DV_CALL_VARIANTS    as DV_CALL_VARIANTS_ONT_R1  } from '../modules/deepvariant'
include { DV_POSTPROCESS      as DV_POSTPROCESS_ONT_R1    } from '../modules/deepvariant'

workflow DEEPVARIANT_R1 {
    take:
    refs
    wm_bams
    bwa_bams

    main:
    def longPlats  = params.dv_long_platforms .tokenize(',').collect { it.trim() }
    def shortPlats = params.dv_short_platforms.tokenize(',').collect { it.trim() }

    // ---- Track A inputs ----
    def hybrid_input = wm_bams
        .filter { hap, ver, plat, bam, bai -> longPlats.contains(plat) }
        .combine(
            bwa_bams.filter { hap, ver, plat, bam, idx -> shortPlats.contains(plat) },
            by: [0, 1]
        )

    // hyb_merged: [ hap, ver_from, combo, bam, bai ]
    // If mapping_dir is set, load hybrid BAMs from there and run MERGE_HYBRID
    // only for hap/combo items absent from the directory (per-item re-entry).
    def hyb_merged
    if ( params.mapping_dir ) {
        def mdir = params.mapping_dir.replaceAll('/$', '')
        def _pfx = "${params.asm_name}_${params.asm_ver}"
        def hyb_existing = Channel.fromPath("${mdir}/${_pfx}.*.*/${_pfx}.*.*.bam")
            .filter { f -> !(f.name =~ /\.(hifi|ont|illumina|element)\.pri\.bam$/) &&
                           !(f.name =~ /\.dedup\.pri\.bam$/) &&
                           (f.name =~ /\.[^.]+_[^.]+\.bam$/) }   // combo must contain '_' (e.g. hifi_illumina)
            .map { f ->
                def m = f.name =~ /\.([^.]+)\.([^.]+)\.bam$/
                if (!m) error "mapping_dir: cannot parse hap/combo from filename: ${f.name}"
                tuple(m[0][1], params.asm_ver, m[0][2], f, file("${f}.bai"))
            }
        def hybrid_input_keyed = hybrid_input
            .map { hap, ver, long_plat, long_bam, long_bai, short_plat, short_bam, short_bai ->
                   tuple(hap, ver, "${long_plat}_${short_plat}",
                         long_plat, long_bam, long_bai, short_plat, short_bam, short_bai) }
        def hyb_missing_input = hybrid_input_keyed
            .join(hyb_existing, by: [0,1,2], remainder: true)
            .filter { it[-1] == null }
            .map    { hap,ver,combo,long_plat,long_bam,long_bai,short_plat,short_bam,short_bai,_null ->
                      tuple(hap, ver, long_plat, long_bam, long_bai, short_plat, short_bam, short_bai) }
        hyb_merged = hyb_existing.mix( MERGE_HYBRID_R1(hyb_missing_input) )
    } else {
        hyb_merged = MERGE_HYBRID_R1(hybrid_input)
    }

    def hyb_merged_ref = hyb_merged.combine(refs, by: [0, 1])
    def hyb_input      = hyb_merged_ref
        .map { hap, ver, combo, bam, bai, ref, gzi, fai ->
               tuple(hap, ver, combo, 'HYBRID_PACBIO_ILLUMINA',
                     (hap == 'dip') ? params.dv_mq_dip.toInteger() : params.dv_mq_hap.toInteger(),
                     ref, gzi, fai, bam, bai) }

    // ---- Track B inputs ----
    def dip_refs         = refs.filter    { hap, ver, ref, gzi, fai -> hap == 'dip' }
    def ont_dip_bams     = wm_bams.filter { hap, ver, plat, bam, bai -> hap == 'dip' && plat == 'ont' }
    def ont_dip_with_ref = ont_dip_bams.combine(dip_refs, by: [0, 1])

    // R9 path
    def r9_input = ( params.ont_chemistry == 'r9' )
        ? ont_dip_with_ref.flatMap { hap, ver, plat, bam, bai, ref, gzi, fai ->
              fai.text.readLines().collect { line -> line.split('\t')[0] }
                  .collect { region -> tuple(hap, ver, region, params.dv_mq_dip.toInteger(), ref, gzi, fai, bam, bai) } }
        : Channel.empty()
    PEPPER_MARGIN_DV_R1(r9_input)
    def r9_merged = DV_MERGE_CHR_VCFS_R1( PEPPER_MARGIN_DV_R1.out.vcfs.groupTuple(by: [0, 1, 2]) )

    // R10 input — defined here so it's in scope for the examples re-entry block below
    def ont_input = ( params.ont_chemistry != 'r9' )
        ? ont_dip_with_ref.map { hap, ver, plat, bam, bai, ref, gzi, fai ->
              tuple(hap, ver, plat, 'ONT_R104', params.dv_mq_dip.toInteger(), ref, gzi, fai, bam, bai) }
        : Channel.empty()

    // ---- DV step 1: make_examples — with per-item re-entry ----
    // If deepvariant_dir is set, load existing examples/ dirs from there and
    // run DV_MAKE_EXAMPLES only for hap/combo/mq items absent from the directory.
    def hyb_ex
    def ont_ex
    if ( params.deepvariant_dir ) {
        def dvdir = params.deepvariant_dir.replaceAll('/$', '')
        def _pfx  = "${params.asm_name}_${params.asm_ver}"
        def ex_ch = Channel.fromPath("${dvdir}/${_pfx}.*.*.MQ*/examples", type: 'dir')
            .map { d ->
                // parent dir name: {asm_name}_{ver}.{hap}.{combo}.MQ{mq}
                def m = d.parent.name =~ /\.([^.]+)\.([^.]+)\.MQ(\d+)$/
                if (!m) error "deepvariant_dir: cannot parse hap/combo/mq from dir: ${d.parent.name}"
                tuple(m[0][1], params.asm_ver, m[0][2], m[0][3].toInteger(), d)
            }
        def hyb_ex_existing = ex_ch
            .filter { hap,ver,combo,mq,d -> combo.contains('illumina') || combo.contains('element') }
            .map    { hap,ver,combo,mq,d -> tuple(hap, ver, combo, 'HYBRID_PACBIO_ILLUMINA', mq, d) }
        def ont_ex_existing = ex_ch
            .filter { hap,ver,combo,mq,d -> !(combo.contains('illumina') || combo.contains('element')) }
            .map    { hap,ver,combo,mq,d -> tuple(hap, ver, combo, 'ONT_R104', mq, d) }

        // Left-outer join: items with no existing examples dir get null as last element
        // → feed those to DV_MAKE_EXAMPLES; mix results with what was already available.
        def hyb_missing = hyb_input
            .join(hyb_ex_existing, by: [0,1,2,3,4], remainder: true)
            .filter { it.last() == null }
            .map    { hap,ver,combo,mode,mq,ref,gzi,fai,bam,bai,_null ->
                      tuple(hap,ver,combo,mode,mq,ref,gzi,fai,bam,bai) }
        def ont_missing = ont_input
            .join(ont_ex_existing, by: [0,1,2,3,4], remainder: true)
            .filter { it.last() == null }
            .map    { hap,ver,combo,mode,mq,ref,gzi,fai,bam,bai,_null ->
                      tuple(hap,ver,combo,mode,mq,ref,gzi,fai,bam,bai) }

        hyb_ex = hyb_ex_existing.mix( DV_MAKE_EXAMPLES_HYB_R1(hyb_missing) )
        ont_ex = ont_ex_existing.mix( DV_MAKE_EXAMPLES_ONT_R1(ont_missing)  )
    } else {
        hyb_ex = DV_MAKE_EXAMPLES_HYB_R1(hyb_input)
        ont_ex = DV_MAKE_EXAMPLES_ONT_R1(ont_input)
    }

    // ---- DV steps 2 + 3 ----
    // hyb_exref / ont_exref provide the ref context that POSTPROCESS needs.
    // When deepvariant_dir is set the BAM channels may be empty (no mapping_dir
    // provided), so build exref directly from examples dirs combined with refs.
    // Also check for existing call_variants_output shards to skip step 2 per item.
    def hyb_exref
    def ont_exref
    def hyb_calls
    def ont_calls
    if ( params.deepvariant_dir ) {
        def dvdir = params.deepvariant_dir.replaceAll('/$', '')
        def _pfx  = "${params.asm_name}_${params.asm_ver}"

        // Discover existing call_variants_output shards: group all shards per hap/combo/mq.
        // Parent dir: {asm_name}_{ver}.{hap}.{combo}.MQ{mq}
        def cvo_ch = Channel.fromPath("${dvdir}/${_pfx}.*.*.MQ*/call_variants_output-0000*-of-0000*.tfrecord.gz")
            .map { f ->
                def m = f.parent.name =~ /\.([^.]+)\.([^.]+)\.MQ(\d+)$/
                if (!m) error "deepvariant_dir: cannot parse hap/combo/mq from dir: ${f.parent.name}"
                tuple(m[0][1], params.asm_ver, m[0][2], m[0][3].toInteger(), f)
            }
            .groupTuple(by: [0,1,2,3])  // [ hap, ver, combo, mq, [shards…] ]

        def hyb_cvo_existing = cvo_ch
            .filter { hap,ver,combo,mq,shards -> combo.contains('illumina') || combo.contains('element') }
            .map    { hap,ver,combo,mq,shards -> tuple(hap, ver, combo, 'HYBRID_PACBIO_ILLUMINA', mq, shards) }
        def ont_cvo_existing = cvo_ch
            .filter { hap,ver,combo,mq,shards -> !(combo.contains('illumina') || combo.contains('element')) }
            .map    { hap,ver,combo,mq,shards -> tuple(hap, ver, combo, 'ONT_R104', mq, shards) }

        // Left-outer join on [hap,ver,combo,mode,mq]: items with no existing shards go to DV_CALL_VARIANTS.
        def hyb_call_missing = hyb_ex
            .join(hyb_cvo_existing.map { hap,ver,combo,mode,mq,shards -> tuple(hap,ver,combo,mode,mq) },
                  by: [0,1,2,3,4], remainder: true)
            .filter { it.last() == null }
            .map    { hap,ver,combo,mode,mq,d,_null -> tuple(hap,ver,combo,mode,mq,d) }
        def ont_call_missing = ont_ex
            .join(ont_cvo_existing.map { hap,ver,combo,mode,mq,shards -> tuple(hap,ver,combo,mode,mq) },
                  by: [0,1,2,3,4], remainder: true)
            .filter { it.last() == null }
            .map    { hap,ver,combo,mode,mq,d,_null -> tuple(hap,ver,combo,mode,mq,d) }

        // For re-entry shards, join with hyb_ex/ont_ex to carry examples alongside —
        // DV_POSTPROCESS needs the examples dir staged in the work directory.
        def hyb_cvo_with_ex = hyb_cvo_existing
            .join(hyb_ex.map { hap,ver,combo,mode,mq,d -> tuple(hap,ver,combo,mode,mq,d) },
                  by: [0,1,2,3,4])
            .map { hap,ver,combo,mode,mq,shards,d -> tuple(hap,ver,combo,mode,mq,shards) }
        def ont_cvo_with_ex = ont_cvo_existing
            .join(ont_ex.map { hap,ver,combo,mode,mq,d -> tuple(hap,ver,combo,mode,mq,d) },
                  by: [0,1,2,3,4])
            .map { hap,ver,combo,mode,mq,shards,d -> tuple(hap,ver,combo,mode,mq,shards) }

        hyb_calls = hyb_cvo_with_ex.mix( DV_CALL_VARIANTS_HYB_R1(hyb_call_missing) )
        ont_calls = ont_cvo_with_ex.mix( DV_CALL_VARIANTS_ONT_R1(ont_call_missing)  )

        hyb_exref = hyb_ex.combine(refs, by: [0, 1])
            .map { hap,ver,combo,mode,mq,examples_dir,ref,gzi,fai ->
                   tuple(hap,ver,combo,mode,mq,ref,gzi,fai,examples_dir) }
        ont_exref = ont_ex.combine(dip_refs, by: [0, 1])
            .map { hap,ver,combo,mode,mq,examples_dir,ref,gzi,fai ->
                   tuple(hap,ver,combo,mode,mq,ref,gzi,fai,examples_dir) }
    } else {
        hyb_calls = DV_CALL_VARIANTS_HYB_R1(hyb_ex)
        ont_calls = DV_CALL_VARIANTS_ONT_R1(ont_ex)
        hyb_exref = hyb_input.join(hyb_ex, by: [0,1,2,3,4])
            .map { hap,ver,combo,mode,mq,ref,gzi,fai,bam,bai,examples_dir ->
                   tuple(hap,ver,combo,mode,mq,ref,gzi,fai,examples_dir) }
        ont_exref = ont_input.join(ont_ex, by: [0,1,2,3,4])
            .map { hap,ver,combo,mode,mq,ref,gzi,fai,bam,bai,examples_dir ->
                   tuple(hap,ver,combo,mode,mq,ref,gzi,fai,examples_dir) }
    }

    // ---- DV step 3: postprocess — with per-item re-entry on completed VCFs ----
    // If deepvariant_dir is set, scan for already-produced VCFs (vcf + tbi + gvcf + gvcf.tbi)
    // and skip DV_POSTPROCESS for those hap/combo/mq items.
    def hyb_post_input = hyb_calls.join(hyb_exref, by:[0,1,2,3,4])
        .map { hap,ver,combo,mode,mq,cvo,ref,gzi,fai,examples_dir ->
               tuple(hap,ver,combo,mode,mq,ref,gzi,fai,cvo,examples_dir) }
    def ont_post_input = ont_calls.join(ont_exref, by:[0,1,2,3,4])
        .map { hap,ver,combo,mode,mq,cvo,ref,gzi,fai,examples_dir ->
               tuple(hap,ver,combo,mode,mq,ref,gzi,fai,cvo,examples_dir) }

    def hyb_post_existing
    def ont_post_existing
    if ( params.deepvariant_dir ) {
        def dvdir = params.deepvariant_dir.replaceAll('/$', '')
        def _pfx  = "${params.asm_name}_${params.asm_ver}"
        // Parent dir: {asm_name}_{ver}.{hap}.{combo}.MQ{mq}
        // VCF file:   dv_{mode}_MQ{mq}.{hap}.vcf.gz  (mode ∈ {HYBRID_PACBIO_ILLUMINA, ONT_R104})
        def vcf_ch = Channel.fromPath("${dvdir}/${_pfx}.*.*.MQ*/dv_*_MQ*.*.vcf.gz")
            .filter { f -> !f.name.contains('.gvcf.') }
            .map { f ->
                def pm = f.parent.name =~ /\.([^.]+)\.([^.]+)\.MQ(\d+)$/
                def fm = f.name =~ /^dv_(.+)_MQ\d+\.[^.]+\.vcf\.gz$/
                if (!pm || !fm) error "deepvariant_dir: cannot parse VCF path: ${f}"
                def hap   = pm[0][1]
                def combo = pm[0][2]
                def mq    = pm[0][3].toInteger()
                def mode  = fm[0][1]
                def base  = "${f.parent}/dv_${mode}_MQ${mq}.${hap}"
                def vcf   = file("${base}.vcf.gz")
                def tbi   = file("${base}.vcf.gz.tbi")
                def gvcf  = file("${base}.gvcf.gz")
                def gtbi  = file("${base}.gvcf.gz.tbi")
                tuple(hap, params.asm_ver, combo, mode, mq, vcf, tbi, gvcf, gtbi)
            }
            .filter { hap,ver,combo,mode,mq,vcf,tbi,gvcf,gtbi ->
                      tbi.exists() && gvcf.exists() && gtbi.exists() }
        hyb_post_existing = vcf_ch.filter { it[3] == 'HYBRID_PACBIO_ILLUMINA' }
        ont_post_existing = vcf_ch.filter { it[3] == 'ONT_R104' }
    } else {
        hyb_post_existing = Channel.empty()
        ont_post_existing = Channel.empty()
    }

    def hyb_post_missing = hyb_post_input
        .join(hyb_post_existing.map { hap,ver,combo,mode,mq,vcf,tbi,gvcf,gtbi -> tuple(hap,ver,combo,mode,mq) },
              by:[0,1,2,3,4], remainder: true)
        .filter { it.last() == null }
        .map    { hap,ver,combo,mode,mq,ref,gzi,fai,cvo,examples_dir,_null ->
                  tuple(hap,ver,combo,mode,mq,ref,gzi,fai,cvo,examples_dir) }
    def ont_post_missing = ont_post_input
        .join(ont_post_existing.map { hap,ver,combo,mode,mq,vcf,tbi,gvcf,gtbi -> tuple(hap,ver,combo,mode,mq) },
              by:[0,1,2,3,4], remainder: true)
        .filter { it.last() == null }
        .map    { hap,ver,combo,mode,mq,ref,gzi,fai,cvo,examples_dir,_null ->
                  tuple(hap,ver,combo,mode,mq,ref,gzi,fai,cvo,examples_dir) }

    DV_POSTPROCESS_HYB_R1(hyb_post_missing)
    DV_POSTPROCESS_ONT_R1(ont_post_missing)

    def dv_vcfs_ch = hyb_post_existing
        .mix( DV_POSTPROCESS_HYB_R1.out[0] )
        .mix( ont_post_existing )
        .mix( DV_POSTPROCESS_ONT_R1.out[0] )
        .mix( r9_merged )

    emit:
    dv_vcfs = dv_vcfs_ch
}
