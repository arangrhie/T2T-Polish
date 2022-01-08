import argparse
import pysam
from pysam import VariantFile, FastaFile


def generate_homopolymer_plots(bed_file, fasta_file, bam_file):
    bed_file_records = open(bed_file, 'r')
    for line in bed_file_records:
        contig, start_pos, end_pos = line.rstrip().split('\t')
        start_pos = int(start_pos)
        end_pos = int(end_pos)
        if start_pos < 1000:
            continue
        if end_pos - start_pos > 50:
            continue

        samfile = pysam.AlignmentFile(bam_file, "rb")

        assembly_fasta_file = FastaFile(fasta_file)
        reference_sequence = assembly_fasta_file.fetch(reference=contig, start=start_pos, end=start_pos + 200)

        reference_homopolymer_index_start = 1
        reference_homopolymer_index_end = 1
        homopolymer_base = reference_sequence[reference_homopolymer_index_start]
        # print(homopolymer_base)
        while reference_homopolymer_index_end < len(reference_sequence) and reference_sequence[reference_homopolymer_index_end] == homopolymer_base:
            reference_homopolymer_index_end += 1

        # print(reference_sequence[reference_homopolymer_index_start:reference_homopolymer_index_end])
        reference_homopolymer_length = reference_homopolymer_index_end - reference_homopolymer_index_start

        all_reads = samfile.fetch(contig, start_pos - 1, end_pos)

        read_homopolymers = []
        for read in all_reads:
            aligned_pairs = read.get_aligned_pairs()

            start_index = 0
            for index, position in aligned_pairs:
                if index is None:
                    continue
                if position == start_pos:
                    start_index = index + 1
                    break
            if read.query_sequence is None:
                continue
            if start_index == len(read.query_sequence):
                continue
            homopolymer_base = read.query_sequence[start_index]
            # print(homopolymer_base)
            end_index = start_index
            while end_index < len(read.query_sequence) and read.query_sequence[end_index] == homopolymer_base:
                end_index += 1
            read_homopolymer_length = end_index - start_index
            read_homopolymers.append(read_homopolymer_length)

        print(contig + "\t" + str(start_pos) + "\t" + str(end_pos) + "\t" + str(reference_homopolymer_length) + "\t" + str(','.join([str(x) for x in read_homopolymers])))


if __name__ == '__main__':
    '''
    Processes arguments and performs tasks.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bed",
        type=str,
        required=True,
        help="Path to a BED file."
    )
    parser.add_argument(
        "--fasta",
        type=str,
        required=False,
        help="Path to a Fasta file."
    )
    parser.add_argument(
        "--bam",
        type=str,
        required=False,
        help="Path to a BAM file."
    )
    FLAGS, unparsed = parser.parse_known_args()

    generate_homopolymer_plots(FLAGS.bed, FLAGS.fasta, FLAGS.bam)

