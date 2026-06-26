/*
 * workflows/evaluation.nf
 *
 * Post-assembly evaluation workflow: Merqury QV, pattern analysis, and issue detection.
 * Runs when params.evaluate = true to provide assembly quality metrics after mapping.
 * 
 * All steps consolidated into a single process to minimize scheduling overhead.
 */

include { FULL_EVALUATION } from '../modules/evaluation'

workflow EVALUATION {
    take:
    asm_dip           // [ name, asm.dip.fa.gz, asm.dip.fa.gz.fai, ver ]
    wm_pri_pafs       // [ hap, ver_from, platform, paf ]
    
    main:
    // Access meryl_db_path directly from params to avoid channel/scalar type issues
    def meryl_hybrid_path = params.hybrid_meryl

    // Extract HiFi and ONT PAFs for the dip haplotype directly from the
    // SAM2PAF output channel. Using the real channel makes FULL_EVALUATION
    // depend on SAM2PAF actually completing — building file() paths from
    // BAM tuples instead would create no dataflow edge and let the process
    // fire before (or without) the PAFs being produced.
    def dip_hifi_paf = wm_pri_pafs
        .filter { hap, v, plat, paf -> hap == 'dip' && plat == 'hifi' }
        .map    { hap, v, plat, paf -> paf }

    def dip_ont_paf = wm_pri_pafs
        .filter { hap, v, plat, paf -> hap == 'dip' && plat == 'ont' }
        .map    { hap, v, plat, paf -> paf }

    // Combine HiFi and ONT PAFs with assembly info into single input tuple for unified evaluation process.
    // NOTE: Nextflow's .combine() FLATTENS tuples, so the resulting items are
    // [name, asm_fa_gz, asm_fai, ver, hifi_paf, ont_paf] — six separate args, not nested.
    def eval_input = asm_dip
        .combine(dip_hifi_paf)
        .combine(dip_ont_paf)
        .map { name, asm_fa_gz, asm_fai, ver, hifi_paf, ont_paf ->
            tuple(name, ver, asm_fa_gz, asm_fai, file(meryl_hybrid_path), hifi_paf, ont_paf)
        }
    
    FULL_EVALUATION(eval_input)

    emit:
    results = FULL_EVALUATION.out
}
