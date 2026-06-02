/*
 * modules/bwa.nf
 *
 * Processes for short-read (Illumina / Element) alignment with BWA-MEM.
 *
 * Pipeline order:
 *   BWA_INDEX  →  BWA_MAP (scatter)  →  BWA_MERGE
 *
 * Replaces:  bwa/_submit_bwa.sh, bwa_index.sh, bwa.sh, bwa/merge.sh
 */

/*
 * Build BWA index for the reference.
 * Replaces bwa_index.sh.  Runs once per hap ref on the quick partition.
 */
process BWA_INDEX {
    label 'quick_small'
    tag "${hap}"
    publishDir "${params.outdir}/assemblies", mode: 'link', overwrite: false,
               saveAs: { fn -> fn == ref_fa_gz.name || fn == ref_fai.name ? null : fn }

    input:
    tuple val(hap), val(ver_from), path(ref_fa_gz), path(ref_fai)

    output:
    tuple val(hap), val(ver_from), path(ref_fa_gz), path(ref_fai),
          path("${ref_fa_gz}.amb"), path("${ref_fa_gz}.ann"),
          path("${ref_fa_gz}.bwt"), path("${ref_fa_gz}.pac"),
          path("${ref_fa_gz}.sa")

    script:
    """
    set -euo pipefail
    module load bwa      || true
    module load samtools || true
    bwa index ${ref_fa_gz}
    """

    stub:
    """
    touch ${ref_fa_gz}.amb ${ref_fa_gz}.ann ${ref_fa_gz}.bwt ${ref_fa_gz}.pac ${ref_fa_gz}.sa
    """
}

/*
 * Align one R1/R2 pair with bwa mem, then fixmate → sort → markdup → filter secondary.
 * Replaces bwa.sh.  One process instance per read pair (Nextflow scatter).
 *
 * Per-pair BAMs are intermediate files.  They are only published to the
 * results directory when params.keep_intermediates = true (default: false).
 */
process BWA_MAP {
    label 'norm_bwa_map'
    tag "${hap}:${ver_from}:${platform}:${r1.name}"
    publishDir "${params.mapping_outdir}/${params.asm_name}_${ver_from}.${hap}.${platform}",
               mode: 'link', overwrite: true,
               enabled: params.keep_intermediates

    input:
    tuple val(hap), val(ver_from), path(ref_fa_gz), path(ref_fai),
          path(amb), path(ann), path(bwt), path(pac), path(sa),
          val(platform), path(r1), path(r2)

    output:
    tuple val(hap), val(ver_from), val(platform),
          path("${r1.simpleName}.dedup.pri.bam"),
          path("${r1.simpleName}.dedup.pri.bam.csi")

    script:
    def tmp = '/lscratch/$SLURM_JOB_ID'
    def out = r1.simpleName
    """
    set -euo pipefail
    module load bwa      || true
    module load samtools || true

    bwa mem -t ${task.cpus} ${ref_fa_gz} ${r1} ${r2} \
        > ${tmp}/${out}.sam

    ${params.samtools} fixmate -m -@${task.cpus} \
        ${tmp}/${out}.sam ${tmp}/${out}.fix.bam
    rm ${tmp}/${out}.sam

    ${params.samtools} sort -@${task.cpus} -O bam \
        -T ${tmp}/${out}.tmp \
        -o ${tmp}/${out}.bam ${tmp}/${out}.fix.bam
    rm ${tmp}/${out}.fix.bam

    ${params.samtools} index -@${task.cpus} ${tmp}/${out}.bam

    ${params.samtools} markdup -r -@${task.cpus} \
        ${tmp}/${out}.bam ${tmp}/${out}.dedup.bam
    rm ${tmp}/${out}.bam

    ${params.samtools} view -@${task.cpus} -F0x100 -hb --write-index \
        -o ${tmp}/${out}.dedup.pri.bam \
        ${tmp}/${out}.dedup.bam

    mv ${tmp}/${out}.dedup.pri.bam     ./
    mv ${tmp}/${out}.dedup.pri.bam.csi ./
    """

    stub:
    """
    touch ${r1.simpleName}.dedup.pri.bam ${r1.simpleName}.dedup.pri.bam.csi
    """
}

/*
 * Merge all per-read-pair dedup.pri BAMs for one (hap, platform).
 * Replaces bwa/merge.sh.
 */
process BWA_MERGE {
    label 'norm_merge_bwa'
    tag "${hap}:${ver_from}:${platform}"
    publishDir "${params.mapping_outdir}/${params.asm_name}_${ver_from}.${hap}.${platform}",
                mode: 'link', overwrite: true

    input:
    tuple val(hap), val(ver_from), val(platform), path(bams), path(csis)

    output:
    tuple val(hap), val(ver_from), val(platform),
          path("${params.asm_name}_${ver_from}.${hap}.${platform}.dedup.pri.bam"),
          path("${params.asm_name}_${ver_from}.${hap}.${platform}.dedup.pri.bam.csi")

    script:
    def out = "${params.asm_name}_${ver_from}.${hap}.${platform}.dedup.pri.bam"
    """
    set -euo pipefail
    module load samtools || true
    num=\$(echo "${bams}" | wc -w)
    if [[ "\$num" -eq 1 ]]; then
        ln -L ${bams} ${out}
        ln -L ${csis} ${out}.csi
    else
        ${params.samtools} merge --write-index -@${task.cpus} -O bam -b <(ls *.dedup.pri.bam) ${out}
    fi
    """

    stub:
    def out = "${params.asm_name}_${ver_from}.${hap}.${platform}.dedup.pri.bam"
    """
    touch ${out} ${out}.csi
    """
}
