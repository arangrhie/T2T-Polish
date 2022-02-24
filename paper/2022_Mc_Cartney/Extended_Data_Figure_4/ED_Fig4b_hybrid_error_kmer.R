library(data.table)
library(ggplot2)
library(scales)

# Color
v09="darkgray"
v10="black"
asmPalette=c(v09, v10)

dat=fread("ED_Fig4b_error_kmer_v0.9_v1.0.tab", header = T)
head(dat)

# Fix order
dat$Chromosome = factor(dat$Chromosome, levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                                                 "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
                                                 "21", "22", "X", "M"))

ggplot(dat, aes(x = Chromosome, y = Errors, fill = Assembly)) +
    geom_bar(stat="identity", position="dodge") +
    scale_fill_manual(values = asmPalette) +
    theme_bw() +
    theme(
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 6, face = "bold"),
        legend.title = element_text(size = 6, face = "bold"),
        legend.text = element_text(size = 6),
        legend.key.size = unit(3, "mm"),
        legend.key.height = unit(3, "mm"), #change legend key height
    )
ggsave("ED_Fig4b_error_kmer_v0.9_v1.0.jpg", width = 5, height = 1.5)
