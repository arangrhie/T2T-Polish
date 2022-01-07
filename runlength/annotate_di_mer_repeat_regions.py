import argparse
import pysam
from pysam import VariantFile, FastaFile


def extend_dimers(reference_sequence, di_mer_base, start_index):
    end_index = start_index + 2
    while end_index < len(reference_sequence) - 1:
        di_mer = reference_sequence[end_index] + reference_sequence[end_index + 1]
        if di_mer == di_mer_base:
            end_index += 2
        else:
            break

    return end_index


def find_dimer_repeats(bam_file, vcf_file, fasta_file):
    assembly_fasta_file = FastaFile(fasta_file)
    small_variant_vcf = VariantFile(vcf_file)
    samfile = pysam.AlignmentFile(bam_file, "rb")

    for rec in small_variant_vcf.fetch():
        alternate_allele = rec.alleles[1]
        if len(alternate_allele) > 50:
            continue
        rec_len = rec.stop - rec.start
        if rec_len > 20:
            continue

        # if rec.contig != 'chr22':
        #     continue
        reference_start = rec.start-200
        reference_end = rec.stop+200
        reference_sequence = assembly_fasta_file.fetch(reference=rec.contig, start=rec.start, end=rec.stop+200)

        in_dimer = False
        end_index = 1
        dimer_base = '**'
        reference_dimer_length = 0
        for i in range(len(reference_sequence)-1):
            if reference_sequence[i] != reference_sequence[i+1]:
                dimer_base = reference_sequence[i] + reference_sequence[i+1]
                end_index = extend_dimers(reference_sequence, dimer_base, i)
                reference_dimer_length = int((end_index - i)/2)
                if i == 1 and reference_dimer_length > 1:
                    # print("----------------------")
                    # print(rec, end='')
                    # print(reference_sequence[i:end_index])
                    # print("FOUND", i, end_index, reference_dimer_length)
                    # print("######################")
                    in_dimer = True
                    break

        if not in_dimer:
            continue

        all_reads = samfile.fetch(rec.contig, rec.start - 10, rec.start+end_index)

        read_dimers = []
        for read in all_reads:
            aligned_pairs = read.get_aligned_pairs()

            read_start_index = -1
            for index, position in aligned_pairs:
                if index is None:
                    continue
                if position == rec.start - 1:
                    read_start_index = index + 2
                    break

            if read_start_index < 0:
                continue
            if read.query_sequence is None:
                continue
            read_end_index = read_start_index + end_index
            if read_start_index >= len(read.query_sequence) or read_end_index >= len(read.query_sequence):
                continue

            read_sequence = read.query_sequence[read_start_index:]

            read_dimer_base = read_sequence[0:2]
            if read_dimer_base != dimer_base:
                continue

            read_end_index_late = extend_dimers(read_sequence, dimer_base, 0)
            read_dimer_length = int((read_end_index_late - 0)/2)
            read_dimers.append(read_dimer_length)
            # print(read.query_sequence[read_start_index:read_end_index], read_start_index, read_end_index, read_dimer_length)

        if len(read_dimers) == 0:
            continue
        print(str(rec.contig) + "\t" + str(rec.start) + "\t" + str(rec.start+end_index) + "\t" + str(reference_dimer_length) + "\t" + str(','.join([str(x) for x in read_dimers])))





if __name__ == '__main__':
    '''
    Processes arguments and performs tasks.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bam",
        "-b",
        type=str,
        required=True,
        help="Path to a VCF file."
    )
    parser.add_argument(
        "--vcf",
        "-v",
        type=str,
        required=True,
        help="Path to a VCF file."
    )
    parser.add_argument(
        "--fasta",
        "-f",
        type=str,
        required=False,
        help="Path to the assembly file."
    )
    FLAGS, unparsed = parser.parse_known_args()

    find_dimer_repeats(FLAGS.bam, FLAGS.vcf, FLAGS.fasta)

