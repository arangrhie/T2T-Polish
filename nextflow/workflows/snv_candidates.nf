/*
 * workflows/snv_candidates.nf
 *
 * SNV_CANDIDATES sub-workflow
 *
 * Collects SNV polishing candidates from four DeepVariant VCFs, runs bcftools
 * filtering / intersection, validates with Merfin, applies the final candidate
 * set to all three haplotype references (hap1, hap2, dip), and prepares the
 * polished FASTAs as inputs for the next polishing round.
 *
 * Takes:
 *   dv_vcfs   — channel from DEEPVARIANT, items:
 *                 [hap, combo, mode, mq, vcf_gz, vcf_gz_tbi, gvcf_gz, gvcf_gz_tbi]
 *   refs      — channel of [ hap, ref_fa_gz, ref_gzi, ref_fai ] (hap1, hap2, dip)
 *   ver_from  — String: current assembly version tag (e.g. 'v0.1')
 *   ver_to    — String: polished assembly version tag (e.g. 'v0.2')
 *
 * Emits:
 *   merfin_loose_vcf  — snv_candidates.merfin-loose.vcf.gz
 *   candidates_vcf    — snv_candidates.vcf.gz (pre-Merfin-loose concat)
 *   next_refs         — [ hap, fa.gz, gzi, fai ] for hap1, hap2, dip — ready for round N+1
 *   chain             — per-hap liftover chain files
 *
 * Required params:
 *   params.hybrid_meryl   — path to hybrid read k-mer meryl DB
 *   params.merfin_peak    — integer peak coverage for Merfin
 *   params.snv_outdir     — output directory
 */

include { SNV_FILTER_INTERSECT  } from '../modules/snv_candidates'
include { SNV_MERFIN            } from '../modules/snv_candidates'
include { SNV_APPLY_CONSENSUS   } from '../modules/snv_candidates'
include { PREPARE_NEXT_ROUND    } from '../modules/snv_candidates'

workflow SNV_CANDIDATES {
    take:
    dv_vcfs   // [hap, ver_from, combo, mode, mq, vcf_gz, tbi, gvcf_gz, gvcf_gz_tbi]
    refs      // [hap, ver_from, ref_fa_gz, ref_fai]
    ver_from  // String — e.g. 'assembly_v0.1'
    ver_to    // String — e.g. 'assembly_v0.2'

    main:

    // -------------------------------------------------------------------------
    // Identify the four VCFs we need from the DV output channel.
    //
    //   hybrid_to_dip   — mode=HYBRID_PACBIO_ILLUMINA, hap=dip,  MQ=dv_mq_dip
    //   hybrid_to_hap1  — mode=HYBRID_PACBIO_ILLUMINA, hap=hap1, MQ=dv_mq_hap
    //   hybrid_to_hap2  — mode=HYBRID_PACBIO_ILLUMINA, hap=hap2, MQ=dv_mq_hap
    //   ont_to_dip      — mode=ONT_R* (any),           hap=dip,  MQ=dv_mq_dip
    // -------------------------------------------------------------------------

    hyb_dip_vcf = dv_vcfs
        .filter { hap, ver, combo, mode, mq, vcf, tbi, gvcf, gvcf_tbi ->
            hap  == 'dip' &&
            mode ==~ /HYBRID_PACBIO.*/ &&
            mq   == params.dv_mq_dip }
        .map { hap, ver, combo, mode, mq, vcf, tbi, gvcf, gvcf_tbi -> [vcf, tbi] }

    hyb_hap1_vcf = dv_vcfs
        .filter { hap, ver, combo, mode, mq, vcf, tbi, gvcf, gvcf_tbi ->
            hap  == 'hap1' &&
            mode ==~ /HYBRID_PACBIO.*/ &&
            mq   == params.dv_mq_hap }
        .map { hap, ver, combo, mode, mq, vcf, tbi, gvcf, gvcf_tbi -> [vcf, tbi] }

    hyb_hap2_vcf = dv_vcfs
        .filter { hap, ver, combo, mode, mq, vcf, tbi, gvcf, gvcf_tbi ->
            hap  == 'hap2' &&
            mode ==~ /HYBRID_PACBIO.*/ &&
            mq   == params.dv_mq_hap }
        .map { hap, ver, combo, mode, mq, vcf, tbi, gvcf, gvcf_tbi -> [vcf, tbi] }

    ont_dip_vcf = dv_vcfs
        .filter { hap, ver, combo, mode, mq, vcf, tbi, gvcf, gvcf_tbi ->
            hap  == 'dip' &&
            mode ==~ /ONT.*/ &&
            mq   == params.dv_mq_dip }
        .map { hap, ver, combo, mode, mq, vcf, tbi, gvcf, gvcf_tbi -> [vcf, tbi] }

    // Combine the four VCF pairs into a single tuple for SNV_FILTER_INTERSECT
    vcf_inputs = ont_dip_vcf
        .combine(hyb_dip_vcf)
        .combine(hyb_hap1_vcf)
        .combine(hyb_hap2_vcf)
        .map { ont_vcf, ont_tbi,
               hyb_dip_vcf_, hyb_dip_tbi_,
               hap1_vcf, hap1_tbi,
               hap2_vcf, hap2_tbi ->
            tuple(ont_vcf, ont_tbi,
                  hyb_dip_vcf_, hyb_dip_tbi_,
                  hap1_vcf, hap1_tbi,
                  hap2_vcf, hap2_tbi)
        }

    // Extract the dip reference for Merfin (heavy node needs the dip ref)
    dip_ref = refs
        .filter { hap, ver, fa, fai -> hap == 'dip' }
        .map    { hap, ver, fa, fai -> [fa, fai] }
        .first()

    // -------------------------------------------------------------------------
    // Run bcftools filtering / intersection (light node)
    // -------------------------------------------------------------------------
    SNV_FILTER_INTERSECT(ver_from, ver_to, vcf_inputs)

    // -------------------------------------------------------------------------
    // Run Merfin -strict and -loose (heavy node)
    // -------------------------------------------------------------------------
    SNV_MERFIN(
        ver_from,
        ver_to,
        SNV_FILTER_INTERSECT.out,
        dip_ref.map { fa, fai -> fa },
        dip_ref.map { fa, fai -> fai },
        file(params.hybrid_meryl)
    )

    // -------------------------------------------------------------------------
    // Apply the validated candidates to the dip reference (light node).
    // hap1 and hap2 are extracted from the polished dip in PREPARE_NEXT_ROUND.
    // -------------------------------------------------------------------------
    SNV_APPLY_CONSENSUS(
        ver_from,
        ver_to,
        dip_ref.map { fa, fai -> fa },
        dip_ref.map { fa, fai -> fai },
        SNV_MERFIN.out.loose_vcf,
        SNV_MERFIN.out.loose_vcf_tbi
    )

    // -------------------------------------------------------------------------
    // Extract hap1 and hap2 from the polished dip by sequence name, using
    // the FAI files of the current-round hap1/hap2 references.
    // -------------------------------------------------------------------------
    hap1_fai_orig = refs
        .filter { hap, ver, fa, fai -> hap == 'hap1' }
        .map    { hap, ver, fa, fai -> fai }
        .first()

    hap2_fai_orig = refs
        .filter { hap, ver, fa, fai -> hap == 'hap2' }
        .map    { hap, ver, fa, fai -> fai }
        .first()

    PREPARE_NEXT_ROUND(
        ver_from,
        ver_to,
        SNV_APPLY_CONSENSUS.out.polished_dip_fa.map { ver, fa_gz, gzi, fai -> fa_gz },
        SNV_APPLY_CONSENSUS.out.polished_dip_fa.map { ver, fa_gz, gzi, fai -> fai },
        hap1_fai_orig,
        hap2_fai_orig
    )

    // Merge the three per-hap next-round ref tuples into one channel
    next_refs_ch = PREPARE_NEXT_ROUND.out.hap1_ref
        .mix(PREPARE_NEXT_ROUND.out.hap2_ref)
        .mix(SNV_APPLY_CONSENSUS.out.polished_dip_fa
            .map { ver, fa_gz, gzi, fai -> tuple('dip', fa_gz, gzi, fai) })

    def merfin_loose_vcf_ch  = SNV_MERFIN.out.loose_vcf
    def candidates_vcf_ch    = SNV_MERFIN.out.candidates_vcf
    def chain_ch             = SNV_APPLY_CONSENSUS.out.chain

    emit:
    merfin_loose_vcf = merfin_loose_vcf_ch
    candidates_vcf   = candidates_vcf_ch
    next_refs        = next_refs_ch   // [hap, fa.gz, gzi, fai] — feed into next round's BUILD_REFS_FROM_FILES
    chain            = chain_ch
}
