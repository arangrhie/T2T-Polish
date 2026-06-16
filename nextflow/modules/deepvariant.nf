/*
 * modules/deepvariant.nf
 *
 * Processes for hybrid BAM merging and DeepVariant variant calling.
 *
 * Pipeline order (hybrid HiFi+short-read mode):
 *   MERGE_HYBRID  →  DV_MAKE_EXAMPLES  →  DV_CALL_VARIANTS  →  DV_POSTPROCESS
 *
 * Pipeline order (ONT R10 mode, dip only):
 *   DV_MAKE_EXAMPLES (ONT_R104)  →  DV_CALL_VARIANTS  →  DV_POSTPROCESS
 *
 * Pipeline order (ONT R9 mode, dip only):
 *   PEPPER_MARGIN_DV (per chromosome, GPU scatter)  →  DV_MERGE_CHR_VCFS
 *
 * Replaces:  deepvariant/_submit_mrg_hybrid_dv.sh,
 *            deepvariant/_submit_deepvariant_with_minqual.sh,
 *            deepvariant/_submit_ont_r9_pepper_margin_dv.sh,
 *            deepvariant/_submit_winnowmap_dv_ont_r9.sh,
 *            deepvariant/_submit_winnowmap_dv_ont_r10.sh,
 *            deepvariant/merge_hybrid.sh,
 *            deepvariant/ont_r9_pepper_margin_dv.sh,
 *            deepvariant/merge_per_chr_vcfs.sh,
 *            deepvariant/step1_with_minqual.sh,
 *            deepvariant/step2_with_minqual.sh,
 *            deepvariant/step3_with_minqual.sh
 */

/*
 * Merge one HiFi BAM and one short-read (Illumina/Element) BAM into a
 * single hybrid BAM used by DeepVariant HYBRID_PACBIO_ILLUMINA mode.
 * Replaces merge_hybrid.sh.
 */
process MERGE_HYBRID {
    label 'merge_hybrid'
    tag "${hap}:${ver_from}:${long_plat}+${short_plat}"
    publishDir "${params.mapping_outdir}/${params.asm_name}_${ver_from}.${hap}.${long_plat}_${short_plat}",
               mode: 'link', overwrite: true

    input:
    tuple val(hap), val(ver_from),
          val(long_plat),  path(long_bam),  path(long_bai),
          val(short_plat), path(short_bam), path(short_bai)

    output:
    tuple val(hap), val(ver_from), val("${long_plat}_${short_plat}"),
          path("${params.asm_name}_${ver_from}.${hap}.${long_plat}_${short_plat}.bam"),
          path("${params.asm_name}_${ver_from}.${hap}.${long_plat}_${short_plat}.bam.bai")

    script:
    def out = "${params.asm_name}_${ver_from}.${hap}.${long_plat}_${short_plat}.bam"
    """
    set -euo pipefail
    module load samtools/1.21 || true
    ${params.samtools} merge -@${task.cpus} -O bam -o ${out} ${long_bam} ${short_bam}
    ${params.samtools} index -@${task.cpus} ${out}
    """

    stub:
    def out = "${params.asm_name}_${ver_from}.${hap}.${long_plat}_${short_plat}.bam"
    """
    touch ${out} ${out}.bai
    """
}

/*
 * DeepVariant step 1: make_examples (scatter across shards with GNU parallel).
 * Replaces step1_with_minqual.sh.
 *
 * Inputs:
 *   hap        — haplotype tag (hap1 / hap2 / dip)
 *   combo      — platform combo tag used in output names (e.g. hifi_illumina)
 *   mode       — DV mode string: PACBIO | ONT_R104 | WGS | HYBRID_PACBIO_ILLUMINA
 *   minqual    — min mapping quality (use 0 for dip, 5 for hap; -1 means DV default)
 *   ref_fa_gz  — reference FASTA (.gz, with .fai alongside)
 *   bam        — input BAM (merged hybrid or single-platform)
 *   bai        — BAM index
 */
process DV_MAKE_EXAMPLES {
    label 'dv_make_examples'
    tag "${hap}:${ver_from}:${combo}:MQ${minqual}"
    publishDir "${params.dv_outdir}/${params.asm_name}_${ver_from}.${hap}.${combo}.MQ${minqual}",
               mode: 'link', overwrite: true

    input:
    tuple val(hap), val(ver_from), val(combo), val(mode), val(minqual),
          path(ref_fa_gz), path(ref_gzi), path(ref_fai),
          path(bam), path(bai)

    output:
    tuple val(hap), val(ver_from), val(combo), val(mode), val(minqual),
          path("examples")   // entire directory — contents depend on DV version

    script:
    def mq_arg = (minqual.toInteger() < 0) ? '' : "--min_mapping_quality ${minqual}"
    def extra_args
    switch (mode) {
        case 'WGS':
            extra_args = "--channels insert_size ${mq_arg}"
            break
        case 'PACBIO':
            extra_args = "--add_hp_channel --alt_aligned_pileup diff_channels --max_reads_per_partition 600 ${mq_arg} --parse_sam_aux_fields --partition_size 25000 --phase_reads --pileup_image_width 199 --norealign_reads --sort_by_haplotypes --track_ref_reads --vsc_min_fraction_indels 0.12"
            break
        case 'ONT_R104':
            extra_args = "--add_hp_channel --alt_aligned_pileup diff_channels --max_reads_per_partition 600 ${mq_arg} --parse_sam_aux_fields --partition_size 25000 --phase_reads --pileup_image_width 199 --norealign_reads --sort_by_haplotypes --track_ref_reads --vsc_min_fraction_indels 0.12 --vsc_min_fraction_snps 0.08"
            break
        case 'HYBRID_PACBIO_ILLUMINA':
            extra_args = "${mq_arg}"
            break
        default:
            error "DV_MAKE_EXAMPLES: unknown mode '${mode}'"
    }
    def n = params.dv_n_shard
    """
    set -euo pipefail
    module load deepvariant/1.6.1 || true
    module load parallel          || true

    mkdir -p examples logs-parallel

    seq 0 \$(( ${n} - 1 )) \
      | parallel -j ${task.cpus} --halt 2 \
          --joblog logs-parallel/log \
          --res    logs-parallel \
        make_examples        \
          --mode calling     \
          --ref  ${ref_fa_gz} \
          --reads ${bam}     \
          --examples examples/examples.tfrecord@${n}.gz \
          --gvcf  examples/examples.gvcf.tfrecord@${n}.gz \
          --sample_name ${hap} \
          ${extra_args}      \
          --task {}
    """

    stub:
    """
    mkdir -p examples
    touch examples/examples.tfrecord-00000-of-00001.gz
    touch examples/examples.gvcf.tfrecord-00000-of-00001.gz
    """
}

/*
 * DeepVariant step 2: call_variants (GPU).
 * Replaces step2_with_minqual.sh (including its retry-on-zero-size logic,
 * now handled by Nextflow's built-in retry mechanism).
 */
process DV_CALL_VARIANTS {
    label 'dv_call_variants'
    tag "${hap}:${ver_from}:${combo}:MQ${minqual}"
    // call_variants_output shards are always published — they are required by
    // DV_POSTPROCESS (step 3) and also serve as a re-entry checkpoint.
    publishDir "${params.dv_outdir}/${params.asm_name}_${ver_from}.${hap}.${combo}.MQ${minqual}",
               mode: 'link', overwrite: true

    input:
    tuple val(hap), val(ver_from), val(combo), val(mode), val(minqual),
          path("examples")   // directory from DV_MAKE_EXAMPLES

    output:
    tuple val(hap), val(ver_from), val(combo), val(mode), val(minqual),
          path('call_variants_output-*-of-*.tfrecord.gz')

    script:
    def mode_lc = mode.toLowerCase()
    def n       = params.dv_n_shard
    """
    set -euo pipefail
    module load deepvariant/1.6.1 || true

    call_variants \
        --outfile call_variants_output.tfrecord.gz \
        --writer_threads 4 \
        --examples examples/examples.tfrecord@${n}.gz \
        --checkpoint /opt/models/${mode_lc}

    # Abort if any shard is zero-size (silently-dying GPU threads)
    for cvo in call_variants_output-*-of-*.tfrecord.gz; do
        if [[ ! -s \$cvo ]]; then
            echo "call_variants produced zero-size output: \$cvo. Failing task for retry."
            exit 1
        fi
    done
    """

    stub:
    def nw = 4
    """
    for i in \$(seq -w 0 \$(( ${nw} - 1 ))); do
        touch call_variants_output-\$(printf '%05d' \$i)-of-\$(printf '%05d' ${nw}).tfrecord.gz
    done
    """
}

/*
 * DeepVariant step 3: postprocess_variants.
 * Replaces step3_with_minqual.sh.
 * Emits the final VCF and gVCF pair.
 */
process DV_POSTPROCESS {
    label 'dv_postprocess'
    tag "${hap}:${ver_from}:${combo}:MQ${minqual}"
    publishDir "${params.dv_outdir}/${params.asm_name}_${ver_from}.${hap}.${combo}.MQ${minqual}",
               mode: 'link', overwrite: true

    input:
    tuple val(hap), val(ver_from), val(combo), val(mode), val(minqual),
          path(ref_fa_gz), path(ref_gzi), path(ref_fai),
          path(call_variants_shards),   // list of call_variants_output-*-of-*.tfrecord.gz
          path("examples")   // directory from DV_MAKE_EXAMPLES (needed for gvcf tfrecords)

    output:
    tuple val(hap), val(ver_from), val(combo), val(mode), val(minqual),
          path("dv_${mode}_MQ${minqual}.${hap}.vcf.gz"),
          path("dv_${mode}_MQ${minqual}.${hap}.vcf.gz.tbi"),
          path("dv_${mode}_MQ${minqual}.${hap}.gvcf.gz"),
          path("dv_${mode}_MQ${minqual}.${hap}.gvcf.gz.tbi")

    script:
    def n      = params.dv_n_shard
    def sample = hap
    def vcf    = "dv_${mode}_MQ${minqual}.${hap}.vcf.gz"
    def gvcf   = "dv_${mode}_MQ${minqual}.${hap}.gvcf.gz"
    """
    set -euo pipefail
    module load deepvariant/1.6.1 || true

    postprocess_variants \
        --ref       ${ref_fa_gz} \
        --infile    call_variants_output.tfrecord.gz \
        --outfile   ${vcf} \
        --gvcf_outfile ${gvcf} \
        --nonvariant_site_tfrecord_path examples/examples.gvcf.tfrecord@${n}.gz \
        --cpus      ${task.cpus} \
        --sample_name ${params.dv_sample ?: sample}
    """

    stub:
    def vcf  = "dv_${mode}_MQ${minqual}.${hap}.vcf.gz"
    def gvcf = "dv_${mode}_MQ${minqual}.${hap}.gvcf.gz"
    """
    touch ${vcf} ${vcf}.tbi ${gvcf} ${gvcf}.tbi
    """
}

/*
 * ONT R9: run pepper_margin_deepvariant for a single chromosome region.
 * This process is scattered — one task per chromosome from the .fai.
 * Replaces ont_r9_pepper_margin_dv.sh.
 * Note that pepper_margin_deepvariant is not able to parse out the region
 * as expected when a - is in the sequence name.
 * Instead, we are dumping a sub-bam and bai to avoid using -r.
 * pepper_margin_deepvariant is also does not accept .csi indices.
 *
 * Inputs:
 *   hap       — haplotype tag (dip only for this process)
 *   ver_from  — assembly version tag
 *   region    — chromosome / contig name (from ref.fai)
 *   mq        — min mapping quality (0 for dip)
 *   ref_fa_gz — reference FASTA (.gz, with .fai alongside)
 *   ref_fai   — .fai index
 *   bam       — input BAM
 *   bai       — BAM index
 */
process PEPPER_MARGIN_DV {
    label 'dv_pepper_margin'
    tag "${hap}:${region}:MQ${mq}"
    // Per-chromosome output dirs are named dv_ONT_R9_MQ{mq}_{region}/
    publishDir "${params.dv_outdir}/${params.asm_name}.${hap}.ont.MQ${mq}/region",
               mode: 'link', overwrite: true

    input:
    tuple val(hap), val(ver_from), val(region), val(mq),
          path(ref_fa_gz), path(ref_gzi), path(ref_fai),
          path(bam), path(bai)

    output:
    tuple val(hap), val(ver_from), val(mq),
          path("dv_ONT_R9_MQ${mq}_${region}/dv_ONT_R9_MQ${mq}_${region}.vcf.gz"),
          path("dv_ONT_R9_MQ${mq}_${region}/dv_ONT_R9_MQ${mq}_${region}.vcf.gz.tbi"),
            optional: true,
            emit: vcfs
        path("pepper_status.log"),
           optional: true,
           emit: status

    script:
    def mq_opts = (mq.toInteger() >= 0)
                      ? "--pepper_min_mapq ${mq} --dv_min_mapping_quality ${mq}"
                      : ''
    """
    set -euo pipefail
    module load pepper_deepvariant/0.8 || true
    module load samtools    || true

    samtools view -hb -@ ${task.cpus} --write-index -o ${region}.bam##idx##${region}.bam.bai ${bam} ${region}
    set +e
    run_pepper_margin_deepvariant call_variant \
        -b ${region}.bam \
        -f ${ref_fa_gz} \
        -p dv_ONT_R9_MQ${mq}_${region} \
        -o dv_ONT_R9_MQ${mq}_${region} \
        ${mq_opts} \
        -t ${task.cpus} \
        --ont_r9_guppy5_sup \
        --gpu
    pm_rc=\$?
    set -e

    if [[ \${pm_rc} -ne 0 ]]; then
        log_file=dv_ONT_R9_MQ${mq}_${region}/logs/2_margin_haplotag.log
        if [[ -s "\${log_file}" ]]; then
            last_line=\$(tail -n 1 "\${log_file}" 2>/dev/null || true)
            if [[ "\${last_line}" == *"No valid VCF entries found!"* ]]; then
                echo "WARNING :: PEPPER_MARGIN_DV returned no variants for ${region} (expected for some small chromosomes). Marking task complete." > pepper_status.log
                exit 0
            fi
            echo "ERROR :: run_pepper_margin_deepvariant failed (exit \${pm_rc}). Last log line: \${last_line}" >&2
        else
            echo "ERROR :: run_pepper_margin_deepvariant failed (exit \${pm_rc}) and no 2_margin_haplotag.log was found." >&2
        fi
        exit \${pm_rc}
    fi

    """

    stub:
    """
    mkdir -p dv_ONT_R9_MQ${mq}_${region}
    touch dv_ONT_R9_MQ${mq}_${region}/dv_ONT_R9_MQ${mq}_${region}.vcf.gz
    touch dv_ONT_R9_MQ${mq}_${region}/dv_ONT_R9_MQ${mq}_${region}.vcf.gz.tbi
    """
}

/*
 * Merge per-chromosome VCFs produced by PEPPER_MARGIN_DV into a single
 * genome-wide VCF.  Replaces merge_per_chr_vcfs.sh.
 *
 * Inputs:
 *   hap      — haplotype tag (dip)
 *   ver_from — assembly version tag
 *   mq       — MQ threshold used in scatter step
 *   vcf_list — list of per-chromosome VCF.gz paths (from PEPPER_MARGIN_DV)
 *   tbi_list — matching .tbi index paths
 */
process DV_MERGE_CHR_VCFS {
    label 'quick_small'
    tag "${hap}:ONT_R9:MQ${mq}"
    publishDir "${params.dv_outdir}/${params.asm_name}_${ver_from}.${hap}.ont.MQ${mq}",
               mode: 'link', overwrite: true

    input:
    tuple val(hap), val(ver_from), val(mq),
          path(vcfs),   // staged as an array; bcftools will sort/concat
          path(tbis)

    output:
    tuple val(hap), val(ver_from), val("ont"), val("ONT_R9"), val(mq),
          path("dv_ONT_R9_MQ${mq}.${hap}.vcf.gz"),
          path("dv_ONT_R9_MQ${mq}.${hap}.vcf.gz.tbi"),
          path("dv_ONT_R9_MQ${mq}.${hap}.gvcf.gz"),
          path("dv_ONT_R9_MQ${mq}.${hap}.gvcf.gz.tbi")

    script:
    """
    set -euo pipefail
    module load bcftools     || true

    # Write a sorted file list (bcftools concat needs correct contig order)
    ls -v dv_ONT_R9_MQ${mq}_*.vcf.gz \
        > r9_files_to_mrg.list

    bcftools concat -D -a \
        --threads ${task.cpus} \
        --no-version \
        -Oz -o dv_ONT_R9_MQ${mq}.${hap}.vcf.gz \
        -f r9_files_to_mrg.list

    bcftools index --tbi dv_ONT_R9_MQ${mq}.${hap}.vcf.gz

    # Create empty gVCF files (pepper_margin_deepvariant does not produce gVCFs)
    touch dv_ONT_R9_MQ${mq}.${hap}.gvcf.gz
    touch dv_ONT_R9_MQ${mq}.${hap}.gvcf.gz.tbi
    """

    stub:
    """
    touch dv_ONT_R9_MQ${mq}.${hap}.vcf.gz dv_ONT_R9_MQ${mq}.${hap}.vcf.gz.tbi
    touch dv_ONT_R9_MQ${mq}.${hap}.gvcf.gz dv_ONT_R9_MQ${mq}.${hap}.gvcf.gz.tbi
    """
}
