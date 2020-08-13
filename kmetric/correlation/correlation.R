setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(ggplot2)
library(dplyr)
library(ggpubr)

kmers <- read.csv("illumina_vs_hifi20kb.txt", header = FALSE, sep = "\t", na.strings="NA")

kmers <- kmers %>% mutate(missing = if_else(is.na(V2) | is.na(V3), 'missing', 'non_missing'))

png("correlation_readK.png", width = 2000, height = 2000)

# Add regression line
ggplot(sampled_kmers, aes(x=V3, y=V6)) + geom_point() + geom_smooth(method = lm)+
  theme(
    plot.margin = margin(t = 20, b = 20),
    axis.title.x = element_text(family = "Arial",
                                size = rel(4),colour = "black"), 
    axis.title.y = element_text(family = "Arial",
                                size = rel(4),colour = "black"), 
    axis.text = element_text(family = "Arial",
                             size = rel(4), colour = "black"),
    axis.ticks.length = unit(.50, "cm"),
    axis.ticks = element_line(colour = "black", size=(0.6)),
    legend.text=element_text(size=30),
    legend.title=element_blank()
  )+
  xlab("Illumina PCR free readK")+
  ylab("20kb readK")+
  stat_cor(method = "pearson")+
  xlim(c(0,250000))+
  ylim(c(0,250000))+
  geom_jitter()


dev.off()

kmers <- read.csv("illumina_vs_hifi20kb.txt", header = FALSE, sep = "\t", na.strings="NA")

kmers <- kmers %>% mutate(missing = if_else(is.na(V2) | is.na(V3), 'missing', 'non_missing'))

kmers[is.na(kmers)] <- 0

png("correlation_Kstar.png", width = 2000, height = 2000)

# Add regression line
ggplot(kmers%>%arrange(desc(missing)), aes(x=V2, y=V3)) + geom_point(aes(color=missing, size = V1)) +
  geom_smooth(method = lm)+
  theme(
    plot.margin = margin(t = 20, b = 20),
    axis.title.x = element_text(family = "Arial",
                                size = rel(4),colour = "black"), 
    axis.title.y = element_text(family = "Arial",
                                size = rel(4),colour = "black"), 
    axis.text = element_text(family = "Arial",
                             size = rel(4), colour = "black"),
    axis.ticks.length = unit(.50, "cm"),
    axis.ticks = element_line(colour = "black", size=(0.6)),
    legend.text=element_text(size=30),
    legend.title=element_blank(),
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                    colour = "grey"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                    colour = "grey")
  )+
  xlab("Illumina PCR free K*")+
  ylab("20kb K*")+
  xlim(c(-50,50))+
  ylim(c(-50,50))+
  scale_color_manual(breaks = c("missing", "non_missing"),
                     values=c("red", "black"))

dev.off()

# Point + regression line
# Remove the confidence interval 
b + geom_point() + 
  geom_smooth(method = lm, se = FALSE)

# loess method: local regression fitting
b + geom_point() + geom_smooth()