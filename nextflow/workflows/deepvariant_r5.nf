/* Round-5 — per-round process aliases for NF 25.x uniqueness constraint */
include { MERGE_HYBRID        as MERGE_HYBRID_R5          } from '../modules/deepvariant'
include { PEPPER_MARGIN_DV    as PEPPER_MARGIN_DV_R5      } from '../modules/deepvariant'
include { DV_MERGE_CHR_VCFS   as DV_MERGE_CHR_VCFS_R5    } from '../modules/deepvariant'
include { DV_MAKE_EXAMPLES    as DV_MAKE_EXAMPLES_HYB_R5  } from '../modules/deepvariant'
include { DV_CALL_VARIANTS    as DV_CALL_VARIANTS_HYB_R5  } from '../modules/deepvariant'
include { DV_POSTPROCESS      as DV_POSTPROCESS_HYB_R5    } from '../modules/deepvariant'
include { DV_MAKE_EXAMPLES    as DV_MAKE_EXAMPLES_ONT_R5  } from '../modules/deepvariant'
include { DV_CALL_VARIANTS    as DV_CALL_VARIANTS_ONT_R5  } from '../modules/deepvariant'
include { DV_POSTPROCESS      as DV_POSTPROCESS_ONT_R5    } from '../modules/deepvariant'

workflow DEEPVARIANT_R5 {
    take:
    refs
    wm_bams
    bwa_bams

    main:
    def longPlats  = params.dv_long_platforms .tokenize(',').collect { it.trim() }
    def shortPlats = params.dv_short_platforms.tokenize(',').collect { it.trim() }

    // ---- Track A: Hybrid HiFi + short-read ----
    def hybrid_input = wm_bams
        .filter { hap, ver, plat, bam, bai -> longPlats.contains(plat) }
        .combine(
            bwa_bams.filter { hap, ver, plat, bam, idx -> shortPlats.contains(plat) },
            by: [0, 1]
        )
    def hyb_merged       = MERGE_HYBRID_R5(hybrid_input)
    def hyb_merged_ref   = hyb_merged.combine(refs, by: [0, 1])
    def hyb_input        = hyb_merged_ref
        .map { hap, ver, combo, bam, bai, ref, gzi, fai ->
               tuple(hap, ver, combo, 'HYBRID_PACBIO_ILLUMINA',
                     (hap == 'dip') ? params.dv_mq_dip.toInteger() : params.dv_mq_hap.toInteger(),
                     ref, gzi, fai, bam, bai) }
    def hyb_ex    = DV_MAKE_EXAMPLES_HYB_R5(hyb_input)
    def hyb_exref = hyb_input.join(hyb_ex, by: [0,1,2,3,4])
        .map { hap,ver,combo,mode,mq,ref,gzi,fai,bam,bai,examples_dir ->
               tuple(hap,ver,combo,mode,mq,ref,gzi,fai,examples_dir) }
    def hyb_calls = DV_CALL_VARIANTS_HYB_R5(hyb_ex)
    DV_POSTPROCESS_HYB_R5(
        hyb_calls.join(hyb_exref, by:[0,1,2,3,4])
            .map { hap,ver,combo,mode,mq,cvo,ref,gzi,fai,examples_dir ->
                   tuple(hap,ver,combo,mode,mq,ref,gzi,fai,cvo,examples_dir) }
    )

    // ---- Track B: ONT dip-only ----
    def dip_refs         = refs.filter    { hap, ver, ref, gzi, fai -> hap == 'dip' }
    def ont_dip_bams     = wm_bams.filter { hap, ver, plat, bam, bai -> hap == 'dip' && plat == 'ont' }
    def ont_dip_with_ref = ont_dip_bams.combine(dip_refs, by: [0, 1])

    // R9 path — processes always declared; input channel is empty when ont_chemistry != 'r9'
    def r9_input = ( params.ont_chemistry == 'r9' )
        ? ont_dip_with_ref.flatMap { hap, ver, plat, bam, bai, ref, gzi, fai ->
              fai.text.readLines().collect { line -> line.split('\t')[0] }
                  .collect { region -> tuple(hap, ver, region, params.dv_mq_dip.toInteger(), ref, gzi, fai, bam, bai) } }
        : Channel.empty()
    PEPPER_MARGIN_DV_R5(r9_input)
    def r9_merged = DV_MERGE_CHR_VCFS_R5( PEPPER_MARGIN_DV_R5.out.vcfs.groupTuple(by: [0, 1, 2]) )

    // R10 path — input channel is empty when ont_chemistry == 'r9'
    def ont_input = ( params.ont_chemistry != 'r9' )
        ? ont_dip_with_ref.map { hap, ver, plat, bam, bai, ref, gzi, fai ->
              tuple(hap, ver, plat, 'ONT_R104', params.dv_mq_dip.toInteger(), ref, gzi, fai, bam, bai) }
        : Channel.empty()
    def ont_ex    = DV_MAKE_EXAMPLES_ONT_R5(ont_input)
    def ont_exref = ont_input.join(ont_ex, by: [0,1,2,3,4])
        .map { hap,ver,combo,mode,mq,ref,gzi,fai,bam,bai,examples_dir ->
               tuple(hap,ver,combo,mode,mq,ref,gzi,fai,examples_dir) }
    def ont_calls = DV_CALL_VARIANTS_ONT_R5(ont_ex)
    DV_POSTPROCESS_ONT_R5(
        ont_calls.join(ont_exref, by:[0,1,2,3,4])
            .map { hap,ver,combo,mode,mq,cvo,ref,gzi,fai,examples_dir ->
                   tuple(hap,ver,combo,mode,mq,ref,gzi,fai,cvo,examples_dir) }
    )

    def dv_vcfs_ch = DV_POSTPROCESS_HYB_R5.out[0]
        .mix( DV_POSTPROCESS_ONT_R5.out[0] )
        .mix( r9_merged )

    emit:
    dv_vcfs = dv_vcfs_ch
}
