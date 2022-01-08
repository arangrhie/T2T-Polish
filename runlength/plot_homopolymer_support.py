import argparse
import pysam
from pysam import VariantFile, FastaFile
import matplotlib.pyplot as plt
from matplotlib import pyplot
pyplot.rcParams['pdf.fonttype'] = 42
pyplot.switch_backend('agg')
plt.rcParams['ytick.labelsize'] = 5

def plot_homopolmer_support(bed_file, axs):
    bed_file_records = open(bed_file, 'r')

    read_supports = []
    for line in bed_file_records:
        line_split = line.rstrip().split('\t')
        if len(line_split) < 5:
            continue
        contig, start_pos, end_pos, ref_hom_length, read_hom_observed = line_split

        R = int(ref_hom_length)
        read_hom_observed = [int(x)-R for x in read_hom_observed.split(',')]

        if R < 2:
            continue
        if len(read_hom_observed) < 10:
            continue

        support = [r for r in read_hom_observed]
        read_supports.extend(support)

    total_bins = max(read_supports) - min(read_supports)
    n, bins, patches = axs.hist(read_supports, bins=total_bins, width=0.98)
    axs.set_xlim(-5, 5)
    axs.set_xticks([i + 0.5 for i in range(-5, 5, 1)])
    axs.set_xticklabels([str(i) for i in range(-5, 5, 1)], fontsize=5)
    axs.set_xlabel("Observed difference (read - asm)", fontsize=5)
    axs.set_ylabel("Count", fontsize=5)
    axs.set_title(bed_file, fontsize=5)


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
        "--bed2",
        type=str,
        required=True,
        help="Path to a BED2 files for v1.0."
    )
    FLAGS, unparsed = parser.parse_known_args()

    all_beds_v09 = FLAGS.bed.rstrip().split(',')
    all_beds_v10 = FLAGS.bed2.rstrip().split(',')

    fig, axes = plt.subplots(2, len(all_beds_v09), figsize=(6, 3))

    for i, bed in enumerate(all_beds_v09):
        plot_homopolmer_support(bed, axes[0][i])

    for i, bed in enumerate(all_beds_v10):
        plot_homopolmer_support(bed, axes[1][i])
    plt.savefig('./homopolymer_support.png', format='png', dpi=300)
    # plt.show()

