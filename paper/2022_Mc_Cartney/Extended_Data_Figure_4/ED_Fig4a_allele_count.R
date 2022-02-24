library(ggplot2)
library(scales)
library(cowplot)
library(grid)
library(gridExtra)

plot_points_by_qual <- function(dat, title="", QUAL, category, xmax=0, ymax=0) {
    if (xmax==0) {
        xmax = max(dat$REF_CNT)
    }
    if (ymax==0) {
        ymax = max(dat$ALT_CNT)
    }
    p <- ggplot(dat, mapping = aes(x=REF_CNT, y=ALT_CNT, color=QUAL)) +
        theme_classic() +
        geom_point(alpha=0.3, size = 0.2) +
        scale_colour_gradient(low = "white", high = "black") +
        ggtitle(title) +
        scale_x_continuous(labels = comma, limits = c(0, xmax)) +
        scale_y_continuous(labels = comma, limits = c(0, ymax)) +
        labs(color = category) +
        xlab("Ref. Count") + ylab("Alt. Count") +
        theme(
            plot.title = element_text(face="bold", size=7),
            axis.title = element_text(face="bold", size=6),
            axis.text = element_text(face="bold", size=5),
            legend.position = c(0.95,0.80),
            legend.key = element_blank(),
            legend.text = element_text(size = 5),
            legend.title = element_text(face="bold", size = 5),
            legend.key.size = unit(2, "mm"),
            legend.margin=margin(c(0,0,0,0)))
    return(p)
}

# v0.9
dat=read.table("ED_Fig4a_v0.9_combined.cnt", header=F)
names(dat) <- c("Chr", "Pos", "Ref", "Alt", "QUAL", "FILTER", "GQ", "GL", "GT", "REF_CNT", "ALT_CNT")
dat_filt=dat[dat$REF_CNT < 200.0 & dat$ALT_CNT < 200.0 & dat$GT == "HOM", ]
summary(dat_filt)

p1 = plot_points_by_qual(dat_filt, "v0.9", dat$GQ, "GQ", xmax = 200, ymax = 200)
p1

# v1.0
dat=read.table("ED_Fig4/ED_Fig4a_v1.0_combined.cnt", header=F)
names(dat) <- c("Chr", "Pos", "Ref", "Alt", "QUAL", "FILTER", "GQ", "GL", "GT", "REF_CNT", "ALT_CNT")
dat_filt=dat[dat$REF_CNT < 200.0 & dat$ALT_CNT < 200.0 & dat$GT == "HOM", ]
summary(dat_filt)

p2 = plot_points_by_qual(dat_filt, "v1.0", dat$GQ, "GQ", xmax = 200, ymax = 200)
p2

plot_grid(p1, p2, nrow = 1)
ggsave(file = "ED_Fig4a_hom_gq.200x.jpg", width = 2.8, height = 1.5)

