library(ggplot2)
library(scales)
library(data.table)
library(cowplot)

setwd("T2T-Polishing/")

BLUE="#377EB8"
GREEN="#4DAF4A"

plot_hist <- function(dat = NULL, title = "", x_max = 0, y_max, count_col) {
    p <- ggplot(data=dat, aes(x=V1, y=V2, group=count_col)) +
        geom_line(color = count_col, size = 0.2) +
        theme_bw() +
        ggtitle(title) +
        scale_x_continuous("k-mer multiplicity") +
        scale_y_continuous("Count") +
        coord_cartesian(xlim=c(0,x_max), ylim=c(0,y_max)) +
        theme(text = element_text(size = 6),
              plot.title = element_text(size=6, face = "bold"),
              axis.title = element_text(size=6, face = "bold"),
              axis.text  = element_text(size=5))
    return(p)
}

plot_heatmap <- function(dat = NULL, cnt, title = "", x_min = 0, x_max, y_max, pattern, counts, low_col) {
    p <- ggplot(dat, aes(x = V2, y = V1, fill = cnt)) +
        geom_tile() +
        scale_fill_gradient(low = low_col,
                            high = "white",
                            guide = "colorbar") +
        theme_bw() +
        ggtitle(title) +
        labs(fill = counts) +
        scale_x_continuous("k-mer multiplicity", expand = c(0, 0), limits=c(x_min, x_max)) +
        scale_y_continuous(paste(pattern, "in 21-mers", sep=" "), expand = c(0, 0), limits=c(0, y_max)) +
        theme(text = element_text(size = 6),
              plot.title = element_text(size=6, face = "bold"),
              axis.title = element_text(size=6, face = "bold"),
              axis.text  = element_text(size=5),
              legend.key.size = unit(2, "mm"),
              legend.key.height = unit(2, "mm"), #change legend key height
              legend.text = element_text(size=5),
              legend.title = element_text(size=5),
              legend.margin=margin(c(0,0,0,0)),
              legend.position = c(0.85, 0.2))
    return(p)
}

# peak for normalizing copy numbers
peak1=31
peak2=105
cp=3

x_max1=peak1 * cp
x_max2=peak2 * cp

dat1=fread("illumina.0.hifi.hist", header=F)
p1=plot_hist(dat1, "Missing in Illumina, found in HiFi", x_max = x_max1, y_max = 500, count_col = BLUE)

dat2=fread("hifi.0.illm.hist", header=F)
p2=plot_hist(dat2, "Missing in HiFi, found in Illumina", x_max = x_max2, y_max = 150, count_col = GREEN)

p1
p2

dat1=fread("illumina.0.hifi.meryl.GC.hist", header=F)
dat1$LogCount=log(dat1$V3, base = 10)
summary(dat1)

dat2=fread("hifi.0.illm.meryl.GC.hist", header=F)
dat2$LogCount=log(dat2$V3, base = 10)
summary(dat2)
y_max=max(dat1$V1)

p3=plot_heatmap(dat1, dat1$LogCount,
             "Missing in Illumina, found in HiFi",
             x_min = 0, x_max = x_max1, y_max = y_max,
             pattern = "GC", counts = "log(Count)",
             low_col = BLUE)

p4=plot_heatmap(dat2, dat2$LogCount,
             "Missing in HiFi, found in Illumina",
             x_min = 0, x_max = x_max2, y_max = y_max,
             pattern = "GC", counts = "log(Count)",
             low_col = GREEN)
p3
p4

plot_grid(p1, p2, p3, p4, nrow = 2, rel_widths = c(0.5, 0.5), rel_heights = c(0.35, 0.65))
ggsave(paste("Fig2ab.pdf", sep=""), width = 4, height = 3)
