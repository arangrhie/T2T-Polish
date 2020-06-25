library(ggplot2)
library(scales)
library(cowplot)
library(grid)
library(gridExtra)

getwd()
setwd("T2T-Polishing")

plot_points <- function(dat, title="", xmax=0, ymax=0) {
    if (xmax==0) {
        xmax = max(dat$REF_CNT)
    }
    if (ymax==0) {
        ymax = max(dat$ALT_CNT)
    }
    p <- ggplot(dat, mapping = aes(x=REF_CNT, y=ALT_CNT, fill=GT, color=GT)) +
        theme_classic() +
        geom_point(shape=21, alpha=0.3) +
        scale_color_brewer(palette = "Set1") +
        scale_fill_brewer(palette = "Set1") +
        scale_x_continuous(labels = comma, limits = c(0, xmax)) +
        scale_y_continuous(labels = comma, limits = c(0, ymax)) +
        ggtitle(title) +
        theme(
            plot.title = element_text(face="bold", family = "Arial", size=7),
            axis.title = element_text(face="bold", family = "Arial", size=7),
            axis.text = element_text(face="bold", family = "Arial", size=6),
            legend.position = "none",
            legend.key = element_blank())
    return(p)
}

plot_points_by_qual <- function(dat, title="", mid, QUAL, xmax=0, ymax=0) {
    if (xmax==0) {
        xmax = max(dat$REF_CNT)
    }
    if (ymax==0) {
        ymax = max(dat$ALT_CNT)
    }
    p <- ggplot(dat, mapping = aes(x=REF_CNT, y=ALT_CNT, color=QUAL)) +
        theme_classic() +
        geom_point(alpha=0.3) +
        scale_color_gradient2(midpoint=mid, low="blue", mid="white",
                              high="red", space ="Lab" ) +
        ggtitle(title) +
        scale_x_continuous(labels = comma, limits = c(0, xmax)) +
        scale_y_continuous(labels = comma, limits = c(0, ymax)) +
        theme(
            plot.title = element_text(face="bold", family = "Arial", size=7),
            axis.title = element_text(face="bold", family = "Arial", size=7),
            axis.text = element_text(face="bold", family = "Arial", size=6),
            legend.position = "none",
            legend.key = element_blank())
    return(p)
}

# Get legends
legend_gt <- function(dat) {
    p<-ggplot(dat, mapping = aes(x=REF_CNT, y=ALT_CNT, fill=GT, color=GT)) +
        theme_classic() +
        geom_point(shape=21, alpha=0.3) +
        scale_color_brewer(palette = "Set1") +
        scale_fill_brewer(palette = "Set1") +
        theme(
            legend.position = "left",
            legend.key = element_blank(),
            legend.text = element_text(family = "Arial", size = 6),
            legend.title = element_text(face="bold", family = "Arial", size = 7))
    return(get_legend(p))
}

legend_qual <- function(dat, mid, QUAL, category) {
    p<-ggplot(dat, mapping = aes(x=REF_CNT, y=ALT_CNT, color=QUAL)) +
        theme_classic() +
        geom_point(alpha=0.3) +
        scale_color_gradient2(midpoint=mid, low="blue", mid="white",
                              high="red", space ="Lab" ) +
        labs(color = category) +
        theme(
            legend.position = "left",
            legend.key = element_blank(),
            legend.text = element_text(family = "Arial", size = 6),
            legend.title = element_text(face="bold", family = "Arial", size = 7))
    return(get_legend(p))
}


# Missing
dat=read.table("input/10x_deepvariant.chr.PASS.chr_only.allele.cnt", header=F)
dat=read.table("input/10x_deepvariant.realign.PASS.missing.cnt", header=F)
dat=read.table("input/10x_deepvariant.mq60_bq60_realign.PASS.missing.cnt", header=F)
names(dat) <- c("Chr", "Pos", "Ref", "Alt", "QUAL", "FILTER", "GQ", "GT", "REF_CNT", "ALT_CNT")
summary(dat$GT)
chr_only.noMulti=dat[dat$GT != "MULTI",]
head(chr_only.noMulti)
chr_only.noMulti$ALT_CNT=as.numeric(chr_only.noMulti$ALT_CNT)

# Non-Missing
dat=read.table("input/10x_deepvariant.chr.PASS.non_chr_only.allele.cnt", header=F)
dat=read.table("input/10x_deepvariant.realign.PASS.non_missing.cnt", header=F)
dat=read.table("input/10x_deepvariant.mq60_bq60_realign.PASS.non_missing.cnt", header=F)
names(dat) <- c("Chr", "Pos", "Ref", "Alt", "QUAL", "FILTER", "GQ", "GT", "REF_CNT", "ALT_CNT")
summary(dat$GT)
non_chr_only.noMulti=dat[dat$GT != "MULTI",]
head(non_chr_only.noMulti)
non_chr_only.noMulti$ALT_CNT=as.numeric(non_chr_only.noMulti$ALT_CNT)

# Legends
data_sum <- rbind(chr_only.noMulti, non_chr_only.noMulti)
mid_qual<-mean(data_sum$QUAL)
mid_gq<-mean(data_sum$GQ)
mid_qual
mid_gq
legend1=legend_gt(chr_only.noMulti)
legend2=legend_qual(data_sum, mid_qual, data_sum$QUAL, "QUAL")
legend3=legend_qual(data_sum, mid_gq, data_sum$GQ, "GQ")

# Subset to max 70x
non_chr_only.noMulti.sub=non_chr_only.noMulti[non_chr_only.noMulti$REF_CNT < 70 & non_chr_only.noMulti$ALT_CNT < 70,]

p1 = plot_points(chr_only.noMulti, "Merqury \"Missing\"")
p2 = plot_points(non_chr_only.noMulti, "\"Non-Missing\"", xmax = 500, ymax = 500)
p3 = plot_points(non_chr_only.noMulti.sub, "\"Non-Missing, <70x\"", xmax = 70, ymax = 70)

p4 = plot_points_by_qual(chr_only.noMulti, "Merqury \"Missing\"", mid_qual, chr_only.noMulti$QUAL)
p5 = plot_points_by_qual(non_chr_only.noMulti, "\"Non-Missing\"", mid_qual, non_chr_only.noMulti$QUAL, xmax = 500, ymax = 500)
p6 = plot_points_by_qual(non_chr_only.noMulti.sub, "\"Non-Missing, <70x\"", mid_qual, non_chr_only.noMulti.sub$QUAL, xmax = 70, ymax = 70)

p7 = plot_points_by_qual(chr_only.noMulti, "Merqury \"Missing\"", mid_gq, chr_only.noMulti$GQ)
p8 = plot_points_by_qual(non_chr_only.noMulti, "\"Non-Missing\"", mid_gq, non_chr_only.noMulti$GQ, xmax = 500, ymax = 500)
p9 = plot_points_by_qual(non_chr_only.noMulti.sub, "\"Non-Missing, <70x\"", mid_gq, non_chr_only.noMulti.sub$GQ, xmax = 70, ymax = 70)

plot_grid(legend1, p1, p2, p3, legend2, p4, p5, p6, legend3, p7, p8, p9, nrow = 3, rel_widths = c(0.1, 0.3, 0.3, 0.3))
ggsave(file = "output/10x_deepvariant_pass_allele_count.png", width = 8, height = 8)
ggsave(file = "output/10x_deepvariant_pass_realign.png", width = 8, height = 8)

# QUAL > mid_qual
mid_qual=30
chr_only.noMulti.qual=chr_only.noMulti[chr_only.noMulti$QUAL>mid_qual,]
non_chr_only.noMulti.qual=non_chr_only.noMulti[non_chr_only.noMulti$QUAL>mid_qual,]
non_chr_only.noMulti.sub.qual=non_chr_only.noMulti.sub[non_chr_only.noMulti.sub$QUAL>mid_qual,]

# GQ > mid_gq
mid_gq=30
chr_only.noMulti.gq=chr_only.noMulti[chr_only.noMulti$GQ>mid_gq,]
non_chr_only.noMulti.gq=non_chr_only.noMulti[non_chr_only.noMulti$GQ>mid_gq,]
non_chr_only.noMulti.sub.gq=non_chr_only.noMulti.sub[non_chr_only.noMulti.sub$GQ>mid_gq,]

# QUAL > mid_qual & GQ > mid_gq
summary(chr_only.noMulti)
summary(chr_only.noMulti.qual)
summary(chr_only.noMulti.gq)

summary(non_chr_only.noMulti)
summary(non_chr_only.noMulti.qual)
summary(non_chr_only.noMulti.gq)

p1 = plot_points_by_qual(chr_only.noMulti.qual, "Merqury \"Missing\"", mid_qual, chr_only.noMulti.qual$QUAL, xmax = 70, ymax = 70)
p2 = plot_points_by_qual(non_chr_only.noMulti.qual, "\"Non-Missing\"", mid_qual, non_chr_only.noMulti.qual$QUAL, xmax = 500, ymax = 500)
p3 = plot_points_by_qual(non_chr_only.noMulti.sub.qual, "\"Non-Missing, <70x\"", mid_qual, non_chr_only.noMulti.sub.qual$QUAL, xmax = 70, ymax = 70)

p4 = plot_points_by_qual(chr_only.noMulti.gq, "Merqury \"Missing\"", mid_gq, chr_only.noMulti.gq$GQ, xmax = 70, ymax = 70)
p5 = plot_points_by_qual(non_chr_only.noMulti.gq, "\"Non-Missing\"", mid_gq, non_chr_only.noMulti.gq$GQ, xmax = 500, ymax = 500)
p6 = plot_points_by_qual(non_chr_only.noMulti.sub.gq, "\"Non-Missing, <70x\"", mid_gq, non_chr_only.noMulti.sub.gq$GQ, xmax = 70, ymax = 70)

plot_grid(legend2, p1, p2, p3, legend3, p4, p5, p6, nrow = 2, rel_widths = c(0.1, 0.3, 0.3, 0.3))
ggsave(file = "output/10x_deepvariant_pass_allele_count.filt.30.png", width = 8, height = 5)
ggsave(file = "output/10x_deepvariant_pass_realign.filt.30.png", width = 8, height = 5)

