#! /bin/bash
##############################################
### Authors: Ivan Sovic and Ann Mc Cartney ###
### Date: September 2021                   ###
##############################################

### Dependencies: meryl, minimap2, merfin, bcftools, racon, winnowmap, falconc (pbipa)

if [ "$#" -ne 6 ]; then
    echo "Automated polishing of draft genomes."
    echo "Wrong number of arguments."
    echo ""
    echo "Usage:"
    echo "  $0 num_threads iterations <in_draft.fasta/fastq> <in_reads.fasta/fastq> <in_readmers> <out_prefix>"
    echo ""
    echo "num_threads - number of threads to use."
    echo "iterations - number of polishing iterations to perform"
    echo "in_draft_fasta - path to the input FASTA/FASTQ containing draft sequences for polishing."
    echo "in_reads - path to the input reads file, in FASTA/FASTQ format."
    echo "in_readmers - path to a Meryl database of Illumina mers."
    echo "out_prefix - prefix of the output files. A folder will automatically be created if one does not exist."
    echo ""
    exit
fi

set -vex

# Input parameters.
CMD_OPT_NUM_THREADS=$1
CMD_OPT_ITERATIONS=$2
CMD_OPT_IN_DRAFT_FASTA=$3
CMD_OPT_IN_READS=$4
CMD_OPT_IN_READMERS=$5
CMD_OPT_OUT_PREFIX=$6

# Dependencies.
RACON=racon
WINNOWMAP=winnowmap
FALCONC=falconc
MERYL=meryl
MERFIN=merfin
BCFTOOLS=bcftools

function run_one_iteration {
    local out_prefix=$1; shift
    local in_draft=$1; shift
    local threads=$1; shift
    local in_dataset=$1; shift
    local in_readmers=$1; shift         # IlluminaPCRfree.k21.gt1.meryl

    # Get the absolute paths.
    mkdir -p $(dirname ${out_prefix})
    out_prefix=$(realpath ${out_prefix})

    # Generate repeitive 15-mers to downweight.
    local out_winnowmap_bam=${out_prefix}.winnowmap.bam
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x" -o ${out_winnowmap_bam}.count15k.memtime \
    ${MERYL} count k=15 ${in_draft} output merylDB
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x" -o ${out_winnowmap_bam}.print15k.memtime \
    ${MERYL} print greater-than distinct=0.9998 merylDB > ${out_winnowmap_bam}.repetitive_k15.txt

    # Map the reads using Winnowmap.
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x" -o ${out_winnowmap_bam}.memtime \
    ${WINNOWMAP} --MD -W ${out_winnowmap_bam}.repetitive_k15.txt -ax map-pb -t ${threads} ${in_draft} ${in_dataset} | samtools view -Sb > ${out_winnowmap_bam}

    # Sort the BAM file.
    local out_winnowmap_sorted_bam=${out_prefix}.winnowmap.sorted.bam
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x" -o ${out_winnowmap_sorted_bam}.memtime \
    samtools sort --threads ${threads} -o ${out_winnowmap_sorted_bam} ${out_winnowmap_bam}

    # Filtering the BAM file.
    local out_falconc_sam=${out_prefix}.falconc.sam
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x" -o ${out_falconc_sam}.memtime \
    falconc bam-filter-clipped -t -F 0x104 --input-fn ${out_winnowmap_sorted_bam} --output-fn ${out_falconc_sam} --output-count-fn ${out_falconc_sam}.filtered_aln_count.txt 2>&1 | tee ${out_falconc_sam}.falconc.log

    # Polish using Racon.
    local out_racon_fasta=${out_prefix}.racon.fasta
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x" -o ${out_racon_fasta}.memtime \
    ${RACON} -t ${threads} ${in_dataset} ${out_falconc_sam} ${in_draft} -L ${out_racon_fasta} -S > ${out_racon_fasta}

    # Generate the Meryl database.
    local out_meryl=${out_prefix}.racon.meryl
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x" -o ${out_meryl}.memtime \
    ${MERYL} count k=21 ${in_draft} output ${out_meryl}

    # Run Merfin.
    local out_merfin=${out_prefix}.racon.merfin
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x" -o ${out_merfin}.memtime \
    ${MERFIN} -polish -sequence ${in_draft} -seqmers ${out_meryl} -readmers ${in_readmers} -peak 106.7 -vcf ${out_racon_fasta}.vcf -output ${out_merfin} -threads ${threads}

    # Call Consensus
    local out_consensus=${out_prefix}.consensus.fasta
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x" -o ${out_consensus}.bcftools_view.memtime \
    ${BCFTOOLS} view -Oz ${out_merfin}.polish.vcf > ${out_merfin}.polish.vcf.gz
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x" -o ${out_consensus}.bcftools_index.memtime \
    ${BCFTOOLS} index ${out_merfin}.polish.vcf.gz
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x" -o ${out_consensus}.bcftools_consensus.memtime \
    ${BCFTOOLS} consensus ${out_merfin}.polish.vcf.gz -f ${in_draft} -H 1 > ${out_consensus}
}

function run_all {
    local out_prefix=$1; shift
    local threads=$1; shift
    local iterations=$1; shift
    local in_draft=$1; shift
    local in_reads=$1; shift
    local in_readmers=$1; shift         # IlluminaPCRfree.k21.gt1.meryl

    mkdir -p $(dirname ${out_prefix})
    cp ${in_draft} ${out_prefix}.iter_0.consensus.fasta
    for ((i = 0; i < ${iterations} ; i++)); do
        next_i=$(( i + 1 ))
        run_one_iteration ${out_prefix}.iter_${next_i} ${out_prefix}.iter_${i}.consensus.fasta ${threads} ${in_reads} ${in_readmers}
    done
}

run_all ${CMD_OPT_OUT_PREFIX} ${CMD_OPT_NUM_THREADS} ${CMD_OPT_ITERATIONS} ${CMD_OPT_IN_DRAFT_FASTA} ${CMD_OPT_IN_READS} ${CMD_OPT_IN_READMERS}
