library(ggplot2)

dat = read.table("chm13v1.1.perfect.lens")
new=data.frame(idy=dat[,1])
new$Assembly="CHM13v1.1"
dat = read.table("chm13v1.0.perfect.lens")
old=data.frame(idy=dat[,1])
old$Assembly="CHM13v1.0"
vegLengths=rbind(old, new)

ggplot(vegLengths, aes(idy, color = Assembly, weight=idy)) +
    geom_density(alpha = 0.2, position="identity", aes(y = ..density..)) +
    scale_color_brewer(palette = "Set2") +
    scale_x_continuous(breaks=c(1,10,100,1000,10000), trans="log1p", expand=c(0,0)) +
    theme_bw() +
    theme(text = element_text(size = 6),
          plot.title = element_text(size=6, face = "bold"),
          axis.title = element_text(size=6, face = "bold"),
          axis.text  = element_text(size=6),
          legend.key.size = unit(3, "mm"),
          legend.key.height = unit(3, "mm"), #change legend key height
          legend.title = element_text(size=6, face = "bold"),
          legend.text = element_text(size=6),
          legend.margin=margin(c(0,0,0,0)),
          legend.position = c(0.85, 0.2)) +
    xlab("Perfect Telomere (bp)")

ggsave("ED_Fig6b.pdf", width=160, height=60, units = "mm")
