import argparse
import pysam
from pysam import VariantFile, FastaFile


def find_homopolymer_cases(vcf_file, fasta_file):
    small_variant_vcf = VariantFile(vcf_file)
    assembly_fasta_file = FastaFile(fasta_file)

    homopolymer_changes = 0
    total_records = 0
    for rec in small_variant_vcf.fetch():
        alternate_allele = rec.alleles[1]
        if len(alternate_allele) > 50:
            continue
        rec_len = rec.stop - rec.start
        if rec_len > 50:
            continue
        total_records += 1
        reference_start = rec.start-200
        reference_end = rec.stop+200
        reference_sequence = assembly_fasta_file.fetch(reference=rec.contig, start=rec.start-200, end=rec.stop+200)

        homopolymer_positions = [0] * len(reference_sequence)
        for i in range(0, len(reference_sequence)):
            if i == 0:
                homopolymer_positions[i] = reference_start
            elif reference_sequence[i] == reference_sequence[i-1]:
                homopolymer_positions[i] = homopolymer_positions[i-1]
            else:
                homopolymer_positions[i] = i + reference_start

        homopolymer_start = 0
        homopolymer_end = 0
        # print(rec, end='')
        for i in range(0, len(reference_sequence)):
            if i + reference_start == rec.start:
                homopolymer_start = homopolymer_positions[i]

            if i + reference_start > rec.stop and homopolymer_positions[i] != homopolymer_start:
                homopolymer_end = max(homopolymer_positions[i], rec.stop + 1)
                break

        sequence_in_assembly = assembly_fasta_file.fetch(reference=rec.contig, start=rec.start, end=rec.stop + 1)

        polished_homopolymer = assembly_fasta_file.fetch(reference=rec.contig, start=homopolymer_start, end=rec.start) + alternate_allele + assembly_fasta_file.fetch(reference=rec.contig, start=rec.stop, end=homopolymer_end)
        sequence_in_polished = assembly_fasta_file.fetch(reference=rec.contig, start=rec.start, end=rec.start) + alternate_allele + assembly_fasta_file.fetch(reference=rec.contig, start=rec.stop, end=rec.stop+1)

        homopolymer_record_end = homopolymer_start
        while reference_sequence[homopolymer_record_end - reference_start] == reference_sequence[homopolymer_start - reference_start]:
            homopolymer_record_end += 1
        # print(assembly_fasta_file.fetch(reference=rec.contig, start=rec.start-1, end=rec.start))
        # print(alternate_allele)
        # print(assembly_fasta_file.fetch(reference=rec.contig, start=rec.stop, end=rec.stop+10))
        # print("Assembly", sequence_in_assembly)
        # print("Polish", sequence_in_polished)
        # if rec.contig != 'chr22':
        #     continue

        # print(rec, end='')
        # print(sequence_in_assembly)
        # print(sequence_in_polished)
        true_homopolymer = True
        if len(sequence_in_assembly) > 1:
            start_index = 2
            while start_index < len(sequence_in_assembly):
                if sequence_in_assembly[start_index] != sequence_in_assembly[start_index - 1]:
                    true_homopolymer = False
                    break
                start_index += 1

        if len(sequence_in_polished) > 1:
            start_index = 2
            while start_index < len(sequence_in_polished):
                if sequence_in_polished[start_index] != sequence_in_polished[start_index - 1]:
                    true_homopolymer = False
                    break
                start_index += 1

        if not true_homopolymer:
            pass
        else:
            print(rec.contig + "\t" + str(rec.start) + "\t" + str(rec.stop))
        # print("#############")

        # total_records += 1

    # print("TOTAL RECORDS: ", total_records)
    # print("TOTAL HOM_CALLS: ", homopolymer_changes)

if __name__ == '__main__':
    '''
    Processes arguments and performs tasks.
    '''
    parser = argparse.ArgumentParser()
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

    find_homopolymer_cases(FLAGS.vcf, FLAGS.fasta)

