/*
 * workflows/dv_three_step_ont.nf
 *
 * Three-step DeepVariant pipeline for the ONT R10 track (Track B).
 * Processes are included under _ONT aliases so NF 25.x sees them as
 * distinct from the Hybrid track's process bindings.
 *
 * Input:  [ hap, ver_from, combo, mode, mq, ref_fa_gz, ref_fai, bam, bai ]
 * Emits:  dv_vcfs — [ hap, ver_from, combo, mode, mq, vcf_gz, vcf_gz_tbi, gvcf_gz, gvcf_gz_tbi ]
 */

include { DV_MAKE_EXAMPLES  as DV_MAKE_EXAMPLES_ONT  } from '../modules/deepvariant'
include { DV_CALL_VARIANTS  as DV_CALL_VARIANTS_ONT  } from '../modules/deepvariant'
include { DV_POSTPROCESS    as DV_POSTPROCESS_ONT    } from '../modules/deepvariant'

workflow DV_THREE_STEP_ONT {

    take:
    make_ex_input  // [ hap, ver_from, combo, mode, mq, ref_fa_gz, ref_fai, bam, bai ]

    main:

    def examples = DV_MAKE_EXAMPLES_ONT(make_ex_input)
    // [ hap, ver_from, combo, mode, mq, tf_records, gvcf_tf_records ]

    def examples_with_ref = make_ex_input
        .join( examples, by: [0, 1, 2, 3, 4] )
        .map { hap, ver, combo, mode, mq, ref, fai, bam, bai, tf, gvcf_tf ->
               tuple(hap, ver, combo, mode, mq, ref, fai, tf, gvcf_tf) }

    def calls = DV_CALL_VARIANTS_ONT(examples)
    // [ hap, ver_from, combo, mode, mq, call_variants_output ]

    def postprocess_input = calls
        .join( examples_with_ref, by: [0, 1, 2, 3, 4] )
        .map { hap, ver, combo, mode, mq, cvo, ref, fai, _tf, gvcf_tf ->
               tuple(hap, ver, combo, mode, mq, ref, fai, cvo, gvcf_tf) }

    DV_POSTPROCESS_ONT(postprocess_input)

    def dv_vcfs_ch = DV_POSTPROCESS_ONT.out[0]

    emit:
    dv_vcfs = dv_vcfs_ch
    // [ hap, ver_from, combo, mode, mq, vcf_gz, vcf_gz_tbi, gvcf_gz, gvcf_gz_tbi ]
}
