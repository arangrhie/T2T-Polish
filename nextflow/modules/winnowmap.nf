/*
 * modules/winnowmap.nf
 *
 * Processes for long-read (HiFi / ONT) alignment with Winnowmap.
 *
 * Pipeline order:
 *   MERYL_REPETITIVE  →  WINNOWMAP_MAP (scatter)  →  WINNOWMAP_MERGE
 *     →  WINNOWMAP_FILTER  →  SAM2PAF
 *
 * Replaces:  winnowmap/_submit.sh, init.sh, map.sh, merge.sh, filt.sh,
 *            coverage/sam2paf.sh
 */

/*
 * Build repetitive_k15.txt with meryl.
 * Replaces init.sh.  Runs once per (hap, ref) on the quick partition.
 */
process MERYL_REPETITIVE {
    label 'quick_meryl'
    tag "${hap}:${ver_from}"
    publishDir "${params.outdir}/assemblies", mode: 'link', overwrite: false

    input:
    tuple val(hap), val(ver_from), path(ref_fa_gz), path(ref_fai)

    output:
    tuple val(hap), val(ver_from), path(ref_fa_gz), path(ref_fai),
          path("${params.asm_name}_${ver_from}.${hap}.repetitive_k15.txt")

    script:
    def out = "${params.asm_name}_${ver_from}.${hap}.repetitive_k15.txt"
    """
    set -euo pipefail
    module load meryl || true
    meryl count k=15 ${ref_fa_gz} output merylDB
    meryl print greater-than distinct=0.9998 merylDB > ${out}
    """

    stub:
    def out = "${params.asm_name}_${ver_from}.${hap}.repetitive_k15.txt"
    """
    touch ${out}
    """
}

/*
 * Align one read file to the reference with Winnowmap.
 * Replaces map.sh.  One process instance per read file (Nextflow scatter).
 * winnowmap → samtools sort → samtools index
 *
 * Per-read BAMs are intermediate files.  They are only published to the
 * results directory when params.keep_intermediates = true (default: false).
 */
process WINNOWMAP_MAP {
    label 'norm_map'
    tag "${hap}:${ver_from}:${platform}:${reads.name}"
    publishDir "${params.mapping_outdir}/${params.asm_name}_${ver_from}.${hap}.${platform}",
               mode: 'link', overwrite: false,
               enabled: params.keep_intermediates

    input:
    tuple val(hap), val(ver_from), path(ref_fa_gz), path(ref_fai), path(rep_txt),
          val(platform), path(reads)

    output:
    tuple val(hap), val(ver_from), val(platform),
          path("${reads.simpleName}.sort.bam"),
          path("${reads.simpleName}.sort.bam.bai")

    script:
    def preset = (platform == 'hifi') ? 'map-pb' : 'map-ont'
    def cpus   = task.cpus - 2      // leave 2 threads for samtools
    def tmp    = '/lscratch/$SLURM_JOB_ID'
    def out    = "${reads.simpleName}.sort"
    """
    set -euo pipefail
    module load winnowmap/2.03 || true
    module load samtools       || true

    winnowmap --MD -W ${rep_txt} -ax ${preset} -I12g -t${cpus} -y \
        ${ref_fa_gz} ${reads} > ${tmp}/${out}.sam

    ${params.samtools} sort -@${cpus} -m2G \
        -T ${tmp}/${out}.tmp -O bam \
        -o ${tmp}/${out}.bam ${tmp}/${out}.sam
    rm ${tmp}/${out}.sam

    mv ${tmp}/${out}.bam ./
    ${params.samtools} index -@${task.cpus} ${out}.bam
    """

    stub:
    """
    touch ${reads.simpleName}.sort.bam ${reads.simpleName}.sort.bam.bai
    """
}

/*
 * Merge all per-read BAMs for one (hap, platform) into a single BAM.
 * Replaces winnowmap/merge.sh.
 */
process WINNOWMAP_MERGE {
    label 'norm_merge_wm'
    tag "${hap}:${ver_from}:${platform}"
    publishDir "${params.mapping_outdir}/${params.asm_name}_${ver_from}.${hap}.${platform}", mode: 'link', overwrite: true

    input:
    tuple val(hap), val(ver_from), val(platform), path(bams), path(bais)

    output:
    tuple val(hap), val(ver_from), val(platform),
          path("${params.asm_name}_${ver_from}.${hap}.${platform}.bam"),
          path("${params.asm_name}_${ver_from}.${hap}.${platform}.bam.bai")

    script:
    def out = "${params.asm_name}_${ver_from}.${hap}.${platform}.bam"
    """
    set -euo pipefail
    module load samtools/1.21 || true
    ${params.samtools} merge -O bam -@${task.cpus} ${out} ${bams}
    ${params.samtools} index -@${task.cpus} ${out}
    """

    stub:
    def out = "${params.asm_name}_${ver_from}.${hap}.${platform}.bam"
    """
    touch ${out} ${out}.bai
    """
}

/*
 * Filter out unmapped + supplementary reads (-F0x104).
 * Replaces filt.sh.
 */
process WINNOWMAP_FILTER {
    label 'norm_filter'
    tag "${hap}:${ver_from}:${platform}"
    publishDir "${params.mapping_outdir}/${params.asm_name}_${ver_from}.${hap}.${platform}", mode: 'link', overwrite: true

    input:
    tuple val(hap), val(ver_from), val(platform), path(bam), path(bai)

    output:
    tuple val(hap), val(ver_from), val(platform),
          path("${params.asm_name}_${ver_from}.${hap}.${platform}.pri.bam"),
          path("${params.asm_name}_${ver_from}.${hap}.${platform}.pri.bam.bai")

    script:
    def out = "${params.asm_name}_${ver_from}.${hap}.${platform}.pri.bam"
    """
    set -euo pipefail
    module load samtools || true
    ${params.samtools} view -F0x104 -@${task.cpus} -hb ${bam} > ${out}
    ${params.samtools} index -@${task.cpus} ${out}
    """

    stub:
    def out = "${params.asm_name}_${ver_from}.${hap}.${platform}.pri.bam"
    """
    touch ${out} ${out}.bai
    """
}

/*
 * Convert filtered BAM → PAF using paftools.js (k8 + minimap2 misc).
 * Replaces coverage/sam2paf.sh.
 */
process SAM2PAF {
    label 'norm_filter'
    tag "${hap}:${ver_from}:${platform}"
    publishDir "${params.mapping_outdir}/${params.asm_name}_${ver_from}.${hap}.${platform}", mode: 'link', overwrite: true

    input:
    tuple val(hap), val(ver_from), val(platform), path(bam), path(bai)

    output:
    tuple val(hap), val(ver_from), val(platform),
          path("${params.asm_name}_${ver_from}.${hap}.${platform}.pri.paf")

    script:
    def out = "${params.asm_name}_${ver_from}.${hap}.${platform}.pri.paf"
    """
    set -euo pipefail
    module load minimap2/2.28 || true
    module load samtools       || true
    mm2_misc=\$(dirname \$(which minimap2))/misc
    ${params.samtools} view -h -@${task.cpus} ${bam} \
        | ${params.k8} \${mm2_misc}/paftools.js sam2paf - \
        | cut -f1-16 > ${out}
    """

    stub:
    def out = "${params.asm_name}_${ver_from}.${hap}.${platform}.pri.paf"
    """
    touch ${out}
    """
}
