library(data.table)
library(ggplot2)
library(scales)

# Color
BLUE="#377EB8"
GREEN="#4DAF4A"
ORANGE="#FF7F00"

plot_bars <- function(dat=NULL, title="", palette=NULL) {
    p <- ggplot(dat, aes(x = CHR, y = MISSING, fill = TYPE)) +
        geom_col(position = position_stack(reverse = TRUE)) +
        scale_fill_manual(values = palette, labels = c("HiFi Consensus", "Patches")) +
        theme_bw() +
        theme(
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold", size = 7),
            axis.text = element_text(size = 6),
            axis.title = element_blank(),
            legend.title = element_text(size = 5, face = "bold"),
            legend.text = element_text(size = 5),
            legend.key.size = unit(2, "mm"),
            legend.key.height = unit(2, "mm"), #change legend key height
            legend.position = c(0.8, 0.83),
            legend.margin=margin(c(0,0,0,0))
        ) +
        scale_x_discrete(drop = FALSE) +
        scale_y_continuous(label = comma, limits = c(0, 2000)) +
        ggtitle(title) +
        labs(fill = "Found in")
    return(p)
}

Patches="black"

hifiPalette = c(BLUE, Patches)
illmPalette = c(GREEN, Patches)
hybrPalette = c(ORANGE, Patches)

dat=fread("ED_Fig1b-d_missing_kmers.tab", header = T)
head(dat)

# Fix order
dat$CHR = factor(dat$CHR, levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                                                 "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
                                                 "21", "22", "X", "M"))
dat$TYPE = factor(dat$TYPE, levels=c("HiFi_Consensus", "Patches"))
dat_hifi=dat[dat$PLATFORM=="HiFi",]
dat_illm=dat[dat$PLATFORM=="Illumina",]
dat_hybr=dat[dat$PLATFORM=="Hybrid",]

plot_bars(dat = dat_hifi, title = "k-mers missing in HiFi", palette = hifiPalette)
ggsave("ED_Fig1b_missing_kmer_hifi.jpg", width = 3.5, height = 1.7)

plot_bars(dat = dat_illm, title = "k-mers missing in Illumina", palette = illmPalette)
ggsave("ED_Fig1c_missing_kmer_illumina.jpg", width = 3.5, height = 1.7)

plot_bars(dat = dat_hybr, title = "k-mers missing in Hybrid", palette = hybrPalette)
ggsave("ED_Fig1d_missing_kmer_hybrid.jpg", width = 3.5, height = 1.7)
