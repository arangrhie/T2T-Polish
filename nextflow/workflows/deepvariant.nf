/*
 * workflows/deepvariant.nf
 *
 * DEEPVARIANT sub-workflow.
 *
 * Runs two independent variant-calling tracks in parallel:
 *
 *  Track A — Hybrid (HiFi + short-read), all three haps (hap1, hap2, dip):
 *    MERGE_HYBRID  →  DV_THREE_STEP_HYBRID (HYBRID_PACBIO_ILLUMINA)
 *    MQ: hap1/hap2 → params.dv_mq_hap (0), dip → params.dv_mq_dip (5)
 *
 *  Track B — ONT single-platform, dip only:
 *    R10 (params.ont_chemistry = 'r10', default):
 *      DV_THREE_STEP_ONT (ONT_R104)
 *      MQ: params.dv_mq_dip (5)
 *    R9  (params.ont_chemistry = 'r9'):
 *      PEPPER_MARGIN_DV (per-chromosome GPU scatter)  →  DV_MERGE_CHR_VCFS
 *      MQ: params.dv_mq_dip (5)
 *
 * DV_THREE_STEP_HYBRID and DV_THREE_STEP_ONT are aliases of the same
 * DV_THREE_STEP sub-workflow (workflows/dv_three_step.nf).  Two aliases are
 * required because NF 25.x enforces that each process module is bound at most
 * once globally — aliased workflow includes create separate process bindings.
 *
 * Inputs:
 *   refs     — [ hap, ver_from, ref_fa_gz, ref_fai ]
 *   wm_bams  — [ hap, ver_from, platform, bam, bai ]
 *   bwa_bams — [ hap, ver_from, platform, bam, csi ]
 */

include { MERGE_HYBRID;
          PEPPER_MARGIN_DV;
          DV_MERGE_CHR_VCFS } from '../modules/deepvariant'

include { DV_THREE_STEP_HYBRID } from './dv_three_step_hybrid'
include { DV_THREE_STEP_ONT    } from './dv_three_step_ont'

// ---------------------------------------------------------------------------
workflow DEEPVARIANT {
// ---------------------------------------------------------------------------
    take:
    refs     // [ hap, ver_from, ref_fa_gz, ref_fai ]    — one per hap
    wm_bams  // [ hap, ver_from, platform, bam, bai ]    — long-read primary BAMs
    bwa_bams // [ hap, ver_from, platform, bam, csi ]    — short-read merged BAMs

    main:

    def longPlats  = params.dv_long_platforms .tokenize(',').collect { it.trim() }
    def shortPlats = params.dv_short_platforms.tokenize(',').collect { it.trim() }

    // Accumulated DV output channel.
    // Shape: [ hap, ver_from, combo, mode, mq, vcf_gz, vcf_gz_tbi, gvcf_gz, gvcf_gz_tbi ]
    dv_vcfs = Channel.empty()

    // =========================================================================
    // TRACK A — Hybrid HiFi + short-read → HYBRID_PACBIO_ILLUMINA DV
    // Runs on hap1, hap2, and dip in parallel.
    // =========================================================================
    def dv_hybrid_mode = 'HYBRID_PACBIO_ILLUMINA'

    def hybrid_input = wm_bams
        .filter { hap, ver, plat, bam, bai -> longPlats.contains(plat) }
        .combine(
            bwa_bams.filter { hap, ver, plat, bam, idx -> shortPlats.contains(plat) },
            by: [0, 1]   // join on (hap, ver_from)
        )
    // [ hap, ver_from, long_plat, long_bam, long_bai, short_plat, short_bam, short_idx ]

    def merged = MERGE_HYBRID(hybrid_input)
    // [ hap, ver_from, combo, bam, bai ]

    // Attach the per-hap reference joined on (hap, ver_from)
    def merged_with_ref = merged.combine(refs, by: [0, 1])
    // [ hap, ver_from, combo, bam, bai, ref_fa_gz, ref_fai ]

    // Assign MQ by hap: dip → dv_mq_dip (5), hap1/hap2 → dv_mq_hap (0)
    def hybrid_make_ex_input = merged_with_ref
        .map { hap, ver, combo, bam, bai, ref, fai ->
               def mq = (hap == 'dip') ? params.dv_mq_dip.toInteger()
                                       : params.dv_mq_hap.toInteger()
               tuple(hap, ver, combo, dv_hybrid_mode, mq, ref, fai, bam, bai) }
    // [ hap, ver_from, combo, mode, mq, ref_fa_gz, ref_fai, bam, bai ]

    DV_THREE_STEP_HYBRID(hybrid_make_ex_input)

    dv_vcfs = dv_vcfs.mix( DV_THREE_STEP_HYBRID.out.dv_vcfs )

    // =========================================================================
    // TRACK B — ONT single-platform, dip only
    // =========================================================================

    // dip reference only (one per round)
    def dip_refs = refs.filter { hap, ver, ref, fai -> hap == 'dip' }

    // ONT BAMs for dip only
    def ont_dip_bams = wm_bams
        .filter { hap, ver, plat, bam, bai -> hap == 'dip' && plat == 'ont' }
    // [ 'dip', ver_from, 'ont', bam, bai ]

    def ont_dip_with_ref = ont_dip_bams.combine(dip_refs, by: [0, 1])
    // [ 'dip', ver_from, 'ont', bam, bai, ref_fa_gz, ref_fai ]

    if ( params.ont_chemistry == 'r9' ) {

        // -------------------------------------------------------------------
        // R9 path: scatter PEPPER_MARGIN_DV one task per chromosome,
        //          then merge all per-chr VCFs.
        // -------------------------------------------------------------------
        def ont_r9_input = ont_dip_with_ref
            .flatMap { hap, ver, plat, bam, bai, ref, fai ->
                       def mq = params.dv_mq_dip.toInteger()
                       fai.text.readLines()
                           .collect { line -> line.split('\t')[0] }
                           .collect { region ->
                               tuple(hap, region, mq, ref, fai, bam, bai) }
                     }
        // [ 'dip', region, mq, ref_fa_gz, ref_fai, bam, bai ]

        def pepper_r9 = PEPPER_MARGIN_DV(ont_r9_input)
        def per_chr_vcfs = pepper_r9.out.vcfs
        // [ 'dip', mq, vcf.gz, vcf.gz.tbi ]  — one item per chromosome

        pepper_r9.out.status
            .map { it.text.trim() }
            .filter { it }
            .view { msg -> "[PEPPER_MARGIN_DV] ${msg}" }

        def merged_vcf_input = per_chr_vcfs
            .groupTuple(by: [0, 1])
            // [ 'dip', mq, [vcf1, vcf2, ...], [tbi1, tbi2, ...] ]

        DV_MERGE_CHR_VCFS(merged_vcf_input)

    } else {

        // -------------------------------------------------------------------
        // R10 path (default): standard DeepVariant three-step with ONT_R104.
        // -------------------------------------------------------------------
        def ont_r10_make_ex_input = ont_dip_with_ref
            .map { hap, ver, plat, bam, bai, ref, fai ->
                   def mq = params.dv_mq_dip.toInteger()
                   tuple(hap, ver, plat, 'ONT_R104', mq, ref, fai, bam, bai) }
        // [ hap, ver_from, combo, mode, mq, ref_fa_gz, ref_fai, bam, bai ]

        DV_THREE_STEP_ONT(ont_r10_make_ex_input)
        dv_vcfs = dv_vcfs.mix( DV_THREE_STEP_ONT.out.dv_vcfs )
    }

    emit:
    dv_vcfs  // [ hap, ver_from, combo, mode, mq, vcf_gz, vcf_gz_tbi, gvcf_gz, gvcf_gz_tbi ]
}
