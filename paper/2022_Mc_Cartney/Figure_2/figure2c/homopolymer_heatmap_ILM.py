import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

files = ["./inputs/t2t_chm13_20200921_illumina_pcrFree_runlength_matrix/t2t_chm13_v1_Illumina_PCRFree_RLE_matrix.run_lengths.A.tsv",
         "./inputs/t2t_chm13_20200921_illumina_pcrFree_runlength_matrix/t2t_chm13_v1_Illumina_PCRFree_RLE_matrix.run_lengths.C.tsv",
         "./inputs/t2t_chm13_20200921_illumina_pcrFree_runlength_matrix/t2t_chm13_v1_Illumina_PCRFree_RLE_matrix.run_lengths.G.tsv",
         "./inputs/t2t_chm13_20200921_illumina_pcrFree_runlength_matrix/t2t_chm13_v1_Illumina_PCRFree_RLE_matrix.run_lengths.T.tsv"]
names = ["A", "C", "G", "T"]
fig, axes = plt.subplots(nrows=1, ncols=len(files), figsize=(16, 8))
RLE_min = 0
RLE_max = 50
interval = 10
for i in range(len(files)):
    homopolymer_mat = np.loadtxt(files[i], skiprows=1, dtype=np.uint64)
    homopolymer_mat = np.delete(homopolymer_mat, 0, 1)
    # homopolymer_mat = homopolymer_mat + 0.00000001
    homopolymer_mat = homopolymer_mat[RLE_min:RLE_max,RLE_min:RLE_max]
    print(homopolymer_mat)
    exit()

    sum_of_rows = homopolymer_mat.sum(axis=1)
    homopolymer_mat_normalized = homopolymer_mat / sum_of_rows[:, np.newaxis]

    im = axes[i].imshow(homopolymer_mat_normalized, cmap="Blues")

    # We want to show all ticks...
    axes[i].set_xticks([0] + [x for x in np.arange(interval, len(homopolymer_mat_normalized), interval)])
    axes[i].set_yticks([0] + [x for x in np.arange(interval, len(homopolymer_mat_normalized), interval)])
    # ... and label them with the respective list entries
    xy_ticklables = [str(max(1, RLE_min-interval))] + [str(x) for x in np.arange(RLE_min+interval, RLE_max, interval)]

    axes[i].set_xticklabels(xy_ticklables)
    axes[i].set_yticklabels(xy_ticklables)

    axes[i].set_title(names[i], fontsize=16)

    axes[i].set_xlabel("Observed in Read alignments")
    axes[i].set_ylabel("Observed in Assembly")

fig.tight_layout()
plt.suptitle("Run-length distribution observed in Illumina PCRfree read alignments (Primary - q5+).", fontsize=16, y=0.8)
# plt.show()
plt.savefig('./Homopolymer_distribution_ILM_'+ str(RLE_max) +'.png', format='png', dpi=300)