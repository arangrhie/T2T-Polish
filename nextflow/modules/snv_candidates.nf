/*
 * modules/snv_candidates.nf
 *
 * Processes for collecting SNV candidates from DeepVariant VCFs, filtering,
 * intersecting, and validating with Merfin.
 *
 * Mirrors variant_call/snv_candidates.sh, split into two resource tiers:
 *
 *   SNV_FILTER_INTERSECT — reheader, PASS-filter, bcftools isec/concat
 *                          Light node: 12 CPUs / 12 GB
 *
 *   SNV_MERFIN           — merfin -strict on hybrid_to_dip VCF;
 *                          merfin -loose on final snv_candidates
 *                          Heavy node: 12 CPUs / 120 GB
 *
 * Replaces: variant_call/snv_candidates.sh
 *
 * Inputs (all VCFs expected as *.vcf.gz with co-located *.vcf.gz.tbi):
 *   ont_dip_vcf     — ONT-to-dip VCF (from DEEPVARIANT Track B)
 *   hyb_dip_vcf     — Hybrid-to-dip VCF (from DEEPVARIANT Track A, hap=dip)
 *   hyb_hap1_vcf    — Hybrid-to-hap1 VCF (from DEEPVARIANT Track A, hap=hap1)
 *   hyb_hap2_vcf    — Hybrid-to-hap2 VCF (from DEEPVARIANT Track A, hap=hap2)
 *   dip_ref_fa_gz   — Dip reference FASTA (with .fai alongside)
 *   hybrid_meryl    — Pre-built hybrid read k-mer database (params.hybrid_meryl)
 *   peak            — Merfin peak value (params.merfin_peak)
 */

/*
 * Step 1: Reheader all four VCFs to a common sample name, filter for PASS,
 * and run all bcftools concat / isec intersections.
 *
 * Output: a single channel item containing the pre-Merfin intermediate VCFs
 * needed by SNV_MERFIN.
 */
process SNV_FILTER_INTERSECT {
    label 'quick_snv_filter'
    tag "snv_filter:${ver_from}"
    publishDir "${params.snv_outdir}/${ver_from}_to_${ver_to}/intermediates", mode: 'link', overwrite: true,
               enabled: params.keep_snv_intermediates

    input:
    val(ver_from)
    val(ver_to)
    tuple path(ont_dip_vcf),   path(ont_dip_tbi),
          path(hyb_dip_vcf),   path(hyb_dip_tbi),
          path(hyb_hap1_vcf),  path(hyb_hap1_tbi),
          path(hyb_hap2_vcf),  path(hyb_hap2_tbi)

    output:
    tuple path("hybrid_to_dip.PASS.vcf.gz"),      path("hybrid_to_dip.PASS.vcf.gz.tbi"),
          path("hybrid_to_hap.PASS.vcf.gz"),       path("hybrid_to_hap.PASS.vcf.gz.tbi"),
          path("ont_to_dip.PASS.vcf.gz"),          path("ont_to_dip.PASS.vcf.gz.tbi"),
          path("hybrid_to_dip.PASS.merfin_input.vcf.gz"),
          path("hybrid_to_dip.PASS.merfin_input.vcf.gz.tbi"),
          path("snv_pre_merfin.vcf.gz"),           path("snv_pre_merfin.vcf.gz.tbi")

    script:
    """
    set -euo pipefail
    module load bcftools || true

    cpus=${task.cpus}

    # -------------------------------------------------------------------------
    # 1. Reheader all VCFs to a common sample name so intersections are clean.
    # -------------------------------------------------------------------------
    echo "sample" > name.txt
    bcftools reheader -s name.txt ${ont_dip_vcf}  -o ont_to_dip.vcf.gz
    bcftools reheader -s name.txt ${hyb_dip_vcf}  -o hybrid_to_dip.vcf.gz
    bcftools reheader -s name.txt ${hyb_hap1_vcf} -o hybrid_to_hap1.vcf.gz
    bcftools reheader -s name.txt ${hyb_hap2_vcf} -o hybrid_to_hap2.vcf.gz

    tabix -p vcf ont_to_dip.vcf.gz
    tabix -p vcf hybrid_to_dip.vcf.gz
    tabix -p vcf hybrid_to_hap1.vcf.gz
    tabix -p vcf hybrid_to_hap2.vcf.gz

    # -------------------------------------------------------------------------
    # 2. Merge hap1 + hap2 into a single hybrid_to_hap VCF.
    # -------------------------------------------------------------------------
    bcftools concat -a -D \\
        -Oz --threads=\$cpus --no-version \\
        hybrid_to_hap1.vcf.gz hybrid_to_hap2.vcf.gz > hybrid_to_hap.vcf.gz
    bcftools index -t hybrid_to_hap.vcf.gz

    # -------------------------------------------------------------------------
    # 3. Filter all VCFs for PASS calls.
    # -------------------------------------------------------------------------
    filter_PASS() {
        local IN=\$1
        local OUT=\${IN/.vcf.gz/.PASS.vcf.gz}
        bcftools view -f "PASS" -Oz --threads=\$cpus --no-version \$IN > \$OUT
        bcftools index -t -f \$OUT
    }
    filter_PASS hybrid_to_hap.vcf.gz
    filter_PASS hybrid_to_dip.vcf.gz
    filter_PASS ont_to_dip.vcf.gz

    # -------------------------------------------------------------------------
    # 4. Consensus error: sites in hybrid_to_hap AND hybrid_to_dip that are
    #    not heterozygous.  (isec 0002 = shared by both inputs, after -e het)
    # -------------------------------------------------------------------------
    bcftools isec -e 'GT="het"' \\
        -Oz --threads=\$cpus --no-version \\
        -p isec_hap_dip_no_het \\
        hybrid_to_hap.PASS.vcf.gz hybrid_to_dip.PASS.vcf.gz
    cp isec_hap_dip_no_het/0002.vcf.gz     hybrid_and_dip.PASS.no_het.vcf.gz
    cp isec_hap_dip_no_het/0002.vcf.gz.tbi hybrid_and_dip.PASS.no_het.vcf.gz.tbi

    # -------------------------------------------------------------------------
    # 5. Dip homozygous high-confidence sites (GQ>10, AF>0.8).
    # -------------------------------------------------------------------------
    bcftools view -i 'GT="hom" & GQ>10 & AF>0.8' \\
        -Oz --threads=\$cpus --no-version \\
        hybrid_to_dip.PASS.vcf.gz > hybrid_to_dip.PASS.hom_GQ10_AF0.8.vcf.gz
    bcftools index -t hybrid_to_dip.PASS.hom_GQ10_AF0.8.vcf.gz

    # -------------------------------------------------------------------------
    # 6. Hap het (GQ>20) ∩ ONT hom (GQ>20, AF>0.8) ∩ dip alt
    #    — phasing errors fixable with ONT support.
    # -------------------------------------------------------------------------
    bcftools view -i 'GT="het" & GQ>20' \\
        -Oz --threads=\$cpus --no-version \\
        hybrid_to_hap.PASS.vcf.gz > hybrid_to_hap.PASS.het_GQ20.vcf.gz
    bcftools index -t hybrid_to_hap.PASS.het_GQ20.vcf.gz

    bcftools view -i 'GT="hom" & GQ>20 & AF>0.8' \\
        -Oz --threads=\$cpus --no-version \\
        ont_to_dip.PASS.vcf.gz > ont_to_dip.PASS.hom_GQ20_AF80.vcf.gz
    bcftools index -t ont_to_dip.PASS.hom_GQ20_AF80.vcf.gz

    # hap het ∩ ONT hom
    bcftools isec \\
        -p isec_hap_het_ont_hom \\
        -Oz --threads=\$cpus --no-version \\
        hybrid_to_hap.PASS.het_GQ20.vcf.gz ont_to_dip.PASS.hom_GQ20_AF80.vcf.gz
    cp isec_hap_het_ont_hom/0002.vcf.gz     hap_het_and_ont_hom.vcf.gz
    cp isec_hap_het_ont_hom/0002.vcf.gz.tbi hap_het_and_ont_hom.vcf.gz.tbi

    # require dip "alt" support
    bcftools view -i 'GT="alt"' \\
        -Oz --threads=\$cpus --no-version \\
        hybrid_to_dip.PASS.vcf.gz > hybrid_to_dip.PASS.alt.vcf.gz
    bcftools index -t hybrid_to_dip.PASS.alt.vcf.gz

    bcftools isec \\
        -p isec_hap_het_ont_hom_dip_alt \\
        -Oz --threads=\$cpus --no-version \\
        hap_het_and_ont_hom.vcf.gz hybrid_to_dip.PASS.alt.vcf.gz
    cp isec_hap_het_ont_hom_dip_alt/0002.vcf.gz     hap_het_and_ont_hom_and_dip_alt.vcf.gz
    cp isec_hap_het_ont_hom_dip_alt/0002.vcf.gz.tbi hap_het_and_ont_hom_and_dip_alt.vcf.gz.tbi

    # -------------------------------------------------------------------------
    # 7. Merge the three candidate sets (pre-Merfin).
    # -------------------------------------------------------------------------
    bcftools concat -a -D \\
        -Oz --threads=\$cpus --no-version \\
        hybrid_and_dip.PASS.no_het.vcf.gz \\
        hybrid_to_dip.PASS.hom_GQ10_AF0.8.vcf.gz \\
        hap_het_and_ont_hom_and_dip_alt.vcf.gz > snv_pre_merfin.vcf.gz
    bcftools index -t snv_pre_merfin.vcf.gz

    # Force GT=1/1 for Merfin -loose input
    { bcftools view -h snv_pre_merfin.vcf.gz;
      bcftools view -H snv_pre_merfin.vcf.gz \\
        | cut -f1-8 | awk '{print \$0"\\tGT\\t1/1"}'; } \\
        | bcftools view -Oz --threads=\$cpus --no-version > snv_pre_merfin.hom.vcf.gz
    bcftools index -t snv_pre_merfin.hom.vcf.gz

    # Also expose the PASS hybrid_to_dip for merfin -strict step
    cp hybrid_to_dip.PASS.vcf.gz     hybrid_to_dip.PASS.merfin_input.vcf.gz
    cp hybrid_to_dip.PASS.vcf.gz.tbi hybrid_to_dip.PASS.merfin_input.vcf.gz.tbi
    """

    stub:
    """
    touch hybrid_to_dip.PASS.vcf.gz        hybrid_to_dip.PASS.vcf.gz.tbi
    touch hybrid_to_hap.PASS.vcf.gz        hybrid_to_hap.PASS.vcf.gz.tbi
    touch ont_to_dip.PASS.vcf.gz           ont_to_dip.PASS.vcf.gz.tbi
    touch hybrid_to_dip.PASS.merfin_input.vcf.gz hybrid_to_dip.PASS.merfin_input.vcf.gz.tbi
    touch snv_pre_merfin.vcf.gz            snv_pre_merfin.vcf.gz.tbi
    """
}

/*
 * Step 2: Run Merfin.
 *   a. merfin -strict on hybrid_to_dip PASS VCF  (consensus-error filter)
 *   b. merfin -loose  on the merged snv_pre_merfin candidates
 *
 * This step needs 120 GB memory (meryl k-mer lookup) and is submitted to a
 * high-memory Slurm node via the 'norm_snv_merfin' resource label.
 */
process SNV_MERFIN {
    label 'norm_snv_merfin'
    tag "snv_merfin:${ver_from}_to_${ver_to}"
    publishDir "${params.snv_outdir}/${ver_from}_to_${ver_to}", mode: 'link', overwrite: true,
               saveAs: { fn -> fn.startsWith("snv_candidates.merfin-loose") ? fn : null }
    publishDir "${params.snv_outdir}/${ver_from}_to_${ver_to}/intermediates", mode: 'link', overwrite: true,
               saveAs: { fn -> fn.startsWith("snv_candidates.merfin-loose") ? null : fn }

    input:
    val(ver_from)
    val(ver_to)
    tuple path(hyb_dip_pass_vcf),    path(hyb_dip_pass_tbi),
          path(hyb_hap_pass_vcf),    path(hyb_hap_pass_tbi),
          path(ont_dip_pass_vcf),    path(ont_dip_pass_tbi),
          path(merfin_input_vcf),    path(merfin_input_tbi),
          path(snv_pre_merfin_vcf),  path(snv_pre_merfin_tbi)
    path(dip_ref_fa_gz)
    path(dip_ref_fai)
    path(hybrid_meryl)

    output:
    path("hybrid_to_dip.PASS.merfin-strict.vcf.gz"),   emit: strict_vcf
    path("hybrid_to_dip.PASS.merfin-strict.vcf.gz.tbi")
    path("snv_candidates.merfin-loose.vcf.gz"),         emit: loose_vcf
    path("snv_candidates.merfin-loose.vcf.gz.tbi"),     emit: loose_vcf_tbi
    path("snv_candidates.vcf.gz"),                      emit: candidates_vcf
    path("snv_candidates.vcf.gz.tbi")

    script:
    def peak   = params.merfin_peak
    def merfin = params.merfin ?: (params.tools ? "${params.tools}/merfin/1.1/bin/merfin" : 'merfin')
    """
    set -euo pipefail
    module load meryl  || true
    module load bcftools || true

    cpus=${task.cpus}

    # Build sequence k-mers from the dip reference
    seqmer=dip_ref.meryl
    meryl count k=31 ${dip_ref_fa_gz} output \$seqmer

    # -------------------------------------------------------------------------
    # a. merfin -strict: filter hybrid_to_dip PASS VCF
    # -------------------------------------------------------------------------
    ${merfin} -strict -peak ${peak} \\
        -sequence  ${dip_ref_fa_gz} \\
        -seqmers   \$seqmer \\
        -output    hybrid_to_dip.PASS.merfin \\
        -readmers  ${hybrid_meryl} \\
        -vcf       ${merfin_input_vcf}

    mv hybrid_to_dip.PASS.merfin.filter.vcf hybrid_to_dip.PASS.merfin-strict.vcf
    bcftools view -Oz --threads=\$cpus --no-version \\
        hybrid_to_dip.PASS.merfin-strict.vcf > hybrid_to_dip.PASS.merfin-strict.vcf.gz
    bcftools index -t hybrid_to_dip.PASS.merfin-strict.vcf.gz

    # -------------------------------------------------------------------------
    # b. Final candidate merge: add the strict-filtered VCF alongside the
    #    other two candidate sets.
    # -------------------------------------------------------------------------
    bcftools concat -a -D \\
        -Oz --threads=\$cpus --no-version \\
        ${snv_pre_merfin_vcf} \\
        hybrid_to_dip.PASS.merfin-strict.vcf.gz > snv_candidates.vcf.gz
    bcftools index -t snv_candidates.vcf.gz

    # GT → 1/1 for merfin -loose
    { bcftools view -h snv_candidates.vcf.gz;
      bcftools view -H snv_candidates.vcf.gz \\
        | cut -f1-8 | awk '{print \$0"\\tGT\\t1/1"}'; } \\
        | bcftools view -Oz --threads=\$cpus --no-version > snv_candidates.hom.vcf.gz
    bcftools index -t snv_candidates.hom.vcf.gz

    # -------------------------------------------------------------------------
    # c. merfin -loose: final SNV candidate filter
    # -------------------------------------------------------------------------
    ${merfin} -loose -peak ${peak} \\
        -sequence  ${dip_ref_fa_gz} \\
        -seqmers   \$seqmer \\
        -output    snv_candidates.merfin \\
        -readmers  ${hybrid_meryl} \\
        -vcf       snv_candidates.hom.vcf.gz

    mv snv_candidates.merfin.filter.vcf snv_candidates.merfin-loose.vcf
    bcftools view -Oz --threads=\$cpus --no-version \\
        snv_candidates.merfin-loose.vcf > snv_candidates.merfin-loose.vcf.gz
    bcftools index -t snv_candidates.merfin-loose.vcf.gz
    """

    stub:
    """
    touch hybrid_to_dip.PASS.merfin-strict.vcf.gz hybrid_to_dip.PASS.merfin-strict.vcf.gz.tbi
    touch snv_candidates.merfin-loose.vcf.gz       snv_candidates.merfin-loose.vcf.gz.tbi
    touch snv_candidates.vcf.gz                    snv_candidates.vcf.gz.tbi
    """
}

/*
 * Step 3: Apply the validated SNV candidates to the dip reference to produce
 * a polished dip assembly FASTA and a liftover chain file.
 *
 * The VCF is applied only to the dip reference; hap1 and hap2 are then
 * extracted from the polished dip in PREPARE_NEXT_ROUND using the sequence
 * names recorded in the original hap1/hap2 FAI files.  This guarantees that
 * MT, EBV, rDNA, and any other shared sequences appear exactly once in the
 * dip and are consistently carried into hap1/hap2 without duplication.
 *
 * Inputs:
 *   ver_from         — current version string, e.g. 'v0.1'
 *   ver_to           — polished version string, e.g. 'v0.2'
 *   dip_ref_fa_gz    — bgzipped FASTA of the current dip reference
 *   dip_ref_fai      — .fai index of dip_ref_fa_gz
 *   merfin_loose_vcf — validated SNV candidates (GT=1/1)
 *   merfin_loose_tbi — tabix index for the VCF
 *
 * Outputs:
 *   polished_dip_fa  — [ ver_to, <ver_to>.dip.fa.gz, <ver_to>.dip.fa.gz.fai ]
 *   chain            — liftover chain file (dip)
 *
 * Light node: consensus application is fast and memory-light.
 */
process SNV_APPLY_CONSENSUS {
    label 'quick_snv_filter'
    tag "consensus:${ver_from}→${ver_to}"
    publishDir "${params.snv_outdir}/${ver_from}_to_${ver_to}", mode: 'link', overwrite: true

    input:
    val(ver_from)
    val(ver_to)
    path(dip_ref_fa_gz)
    path(dip_ref_fai)
    path(merfin_loose_vcf)
    path(merfin_loose_tbi)

    output:
    tuple val(ver_to),
          path("${params.asm_name}_${ver_to}.dip.fa.gz"),
          path("${params.asm_name}_${ver_to}.dip.fa.gz.gzi"),
          path("${params.asm_name}_${ver_to}.dip.fa.gz.fai"),
          emit: polished_dip_fa
    path("${ver_from}_to_${ver_to}.dip.chain"), emit: chain

    script:
    def out = "${params.asm_name}_${ver_to}.dip.fa.gz"
    """
    set -euo pipefail
    module load bcftools samtools || true

    bcftools consensus -H 1 \\
        --chain ${ver_from}_to_${ver_to}.dip.chain \\
        -f ${dip_ref_fa_gz} \\
        ${merfin_loose_vcf} \\
        | bgzip --index -I ${out}.gzi -@ ${task.cpus} > ${out}
    ${params.samtools} faidx -@ ${task.cpus} ${out}
    """

    stub:
    def out = "${params.asm_name}_${ver_to}.dip.fa.gz"
    """
    touch ${out} ${out}.gzi ${out}.fai
    touch ${ver_from}_to_${ver_to}.dip.chain
    """
}

/*
 * Step 4: Extract hap1 and hap2 from the polished dip FASTA using the
 * sequence names recorded in the original hap1 and hap2 FAI files.
 *
 * This is the correct way to split the dip back into its haplotypes: the
 * per-hap sequence names are stable across rounds (only the bases change),
 * so `samtools faidx -r <(cut -f1 hap.fai)` faithfully recovers each
 * haplotype without duplicating shared sequences (MT, EBV, rDNA, etc.).
 *
 * Inputs:
 *   ver_to        — version string for output filenames (e.g. 'v0.2')
 *   dip_fa_gz     — bgzipped polished dip FASTA from SNV_APPLY_CONSENSUS (.fa.gz)
 *   dip_fai       — .fai index for dip_fa_gz
 *   hap1_fai_orig — FAI of the *current-round* hap1 reference (sequence names)
 *   hap2_fai_orig — FAI of the *current-round* hap2 reference (sequence names)
 *
 * Outputs (each as a separate named emit):
 *   hap1_ref — [ 'hap1', <ver_to>.hap1.fa.gz, <ver_to>.hap1.fa.gz.fai ]
 *   hap2_ref — [ 'hap2', <ver_to>.hap2.fa.gz, <ver_to>.hap2.fa.gz.fai ]
 *   dip_ref  — [ 'dip',  <ver_to>.dip.fa.gz,  <ver_to>.dip.fa.gz.fai  ]
 *              (copied from the SNV_APPLY_CONSENSUS output for co-location)
 */
process PREPARE_NEXT_ROUND {
    label 'quick_snv_filter'
    tag "prepare:${ver_from}_to_${ver_to}"
    publishDir "${params.outdir}/assemblies", mode: 'link', overwrite: true

    input:
    val(ver_from)
    val(ver_to)
    path(dip_fa_gz)
    path(dip_fai)
    path(hap1_fai_orig)
    path(hap2_fai_orig)

    output:
    tuple val('hap1'), path("${params.asm_name}_${ver_to}.hap1.fa.gz"),
                       path("${params.asm_name}_${ver_to}.hap1.fa.gz.gzi"),
                       path("${params.asm_name}_${ver_to}.hap1.fa.gz.fai"), emit: hap1_ref
    tuple val('hap2'), path("${params.asm_name}_${ver_to}.hap2.fa.gz"),
                       path("${params.asm_name}_${ver_to}.hap2.fa.gz.gzi"),
                       path("${params.asm_name}_${ver_to}.hap2.fa.gz.fai"), emit: hap2_ref

    script:
    def h1 = "${params.asm_name}_${ver_to}.hap1.fa.gz"
    def h2 = "${params.asm_name}_${ver_to}.hap2.fa.gz"
    """
    set -euo pipefail
    module load samtools || true

    # Extract hap1 sequences from the polished dip by name
    cut -f1 ${hap1_fai_orig} > hap1_names.txt
    ${params.samtools} faidx -@ ${task.cpus} -r hap1_names.txt ${dip_fa_gz} \\
        | bgzip --index -I ${h1}.gzi -@ ${task.cpus} > ${h1}
    ${params.samtools} faidx -@ ${task.cpus} ${h1}

    # Extract hap2 sequences from the polished dip by name
    cut -f1 ${hap2_fai_orig} > hap2_names.txt
    ${params.samtools} faidx -@ ${task.cpus} -r hap2_names.txt ${dip_fa_gz} \\
        | bgzip --index -I ${h2}.gzi -@ ${task.cpus} > ${h2}
    ${params.samtools} faidx -@ ${task.cpus} ${h2}
    """

    stub:
    def h1 = "${params.asm_name}_${ver_to}.hap1.fa.gz"
    def h2 = "${params.asm_name}_${ver_to}.hap2.fa.gz"
    """
    touch ${h1} ${h1}.gzi ${h1}.fai
    touch ${h2} ${h2}.gzi ${h2}.fai
    """
}
