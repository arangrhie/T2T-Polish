library(ggplot2)
library(scales)
library(data.table)
library(cowplot)

BLUE="#377EB8"
GREEN="#4DAF4A"

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
              axis.text  = element_text(size=6),
              legend.key.size = unit(3, "mm"),
              legend.key.height = unit(3, "mm"), #change legend key height
              legend.text = element_text(size=6),
              legend.title = element_text(size=6, face = "bold")
              )
    return(p)
}

# peak for normalizing copy numbers
peak1=31
peak2=105
cp=3

x_max1=peak1 * cp
x_max2=peak2 * cp

dat1=fread("ED_Fig1a_illumina.0.hifi.meryl.GA_TC.hist", header=F)
dat1$LogCount=log(dat1$V3, base = 10)
summary(dat1)

dat2=fread("ED_Fig1a_hifi.0.illm.meryl.GA_TC.hist", header=F)
dat2$LogCount=log(dat2$V3, base = 10)
summary(dat2)
y_max=max(dat1$V1)

p1=plot_heatmap(dat1, dat1$LogCount,
                "Missing in Illumina, found in HiFi",
                x_min = 0, x_max = x_max1, y_max = y_max,
                pattern = "GA", counts = "log(Count)",
                low_col = BLUE)

p2=plot_heatmap(dat2, dat2$LogCount,
                "Missing in HiFi, found in Illumina",
                x_min = 0, x_max = x_max2, y_max = y_max,
                pattern = "GA", counts = "log(Count)",
                low_col = GREEN)

plot_grid(p1, p2, nrow = 1, rel_widths = c(0.5, 0.5))
ggsave(paste("ED_Fig1a.jpg", sep=""), width = 5, height = 1.8)
