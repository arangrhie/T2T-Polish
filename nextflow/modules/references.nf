/*
 * modules/references.nf
 *
 * Processes for building hap1 / hap2 / dip reference FA.GZ files.
 *
 * params.ebv_fasta_gz is optional.  Set it to null (or omit it) if the
 * assembly does not contain an EBV sequence.
 */

/*
 * Concatenate hap_fa + mito + [ebv] + [rdna] into a single reference and index it.
 * Used for hap1 and hap2.
 */
process BUILD_HAP_REFERENCES {
    label 'norm_build_ref'
    tag "${hap}:${params.asm_ver}"
    publishDir "${params.outdir}/assemblies", mode: 'link', overwrite: true

    input:
        tuple val(hap),
            path(hap_fa, stageAs: 'in_hap/*'),
            path(mito,   stageAs: 'in_mito/*'),
            path(ebv,    stageAs: 'in_ebv/*'),
            path(rdna,   stageAs: 'in_rdna/*')

    output:
    tuple val(hap), path("${params.asm_name}_${params.asm_ver}.${hap}.fa.gz"),
                    path("${params.asm_name}_${params.asm_ver}.${hap}.fa.gz.gzi"),
                    path("${params.asm_name}_${params.asm_ver}.${hap}.fa.gz.fai")

    script:
    def out      = "${params.asm_name}_${params.asm_ver}.${hap}.fa.gz"
    def ebv_arg  = ebv.toString().contains('NO_FILE_') ? '' : "${ebv}"
    def rdna_arg = rdna.toString().contains('NO_FILE_') ? '' : "${rdna}"
    def cpus     = task.cpus
    """
    set -euo pipefail
    module load samtools || true

    # cat if already bgzipped, otherwise compress on the fly
    { for f in ${hap_fa} ${mito} ${ebv_arg} ${rdna_arg}; do
          case "\$f" in *.gz|*.bgz) cat "\$f" ;; *) bgzip -@${cpus} -c "\$f" ;; esac
      done; } > ${out}
    bgzip --reindex -@${cpus} ${out}
    ${params.samtools} faidx -@${cpus} ${out}
    """

    stub:
    def out = "${params.asm_name}_${params.asm_ver}.${hap}.fa.gz"
    """
    touch ${out} ${out}.gzi ${out}.fai
    """
}

/*
 * Concatenate hap1 + hap2 + mito + [ebv] + [rdna] into the dip reference and index it.
 * Replaces the former two-step MAKE_DIP_HAPS → BUILD_HAP_REFERENCES(dip) chain.
 */
process BUILD_DIP_REFERENCE {
    label 'norm_build_ref'
    tag "dip:${params.asm_ver}.dip"
    publishDir "${params.outdir}/assemblies", mode: 'link', overwrite: true

    input:
    path(h1,   stageAs: 'in_h1/*')
    path(h2,   stageAs: 'in_h2/*')
    path(mito, stageAs: 'in_mito/*')
    path(ebv,  stageAs: 'in_ebv/*')
    path(rdna, stageAs: 'in_rdna/*')

    output:
    tuple val('dip'), path("${params.asm_name}_${params.asm_ver}.dip.fa.gz"),
                      path("${params.asm_name}_${params.asm_ver}.dip.fa.gz.gzi"),
                      path("${params.asm_name}_${params.asm_ver}.dip.fa.gz.fai")

    script:
    def out      = "${params.asm_name}_${params.asm_ver}.dip.fa.gz"
    def mito_arg  = mito.toString().contains('NO_FILE_') ? '' : "${mito}"
    def ebv_arg  = ebv.toString().contains('NO_FILE_') ? '' : "${ebv}"
    def rdna_arg = rdna.toString().contains('NO_FILE_') ? '' : "${rdna}"
    def cpus     = task.cpus
    """
    set -euo pipefail
    module load samtools || true

    # cat if already bgzipped, otherwise compress on the fly
    { for f in ${h1} ${h2} ${mito_arg} ${ebv_arg} ${rdna_arg}; do
          case "\$f" in *.gz|*.bgz) cat "\$f" ;; *) bgzip -@${cpus} -c "\$f" ;; esac
      done; } > ${out}
    bgzip --reindex -@${cpus} ${out}
    ${params.samtools} faidx -@${cpus} ${out}
    """

    stub:
    def out = "${params.asm_name}_${params.asm_ver}.dip.fa.gz"
    """
    touch ${out} ${out}.gzi ${out}.fai
    """
}
