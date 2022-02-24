require("ggplot2")
require("scales")

dat=read.table("ED_Fig3cd_single_cum_10k.hist", header=FALSE)
head(dat)

fancy_scientific <- function(d) {
    # turn in to character string in scientific notation
    d <- format(d, scientific = TRUE)
    # quote the part before the exponent to keep all the digits and turn the 'e+' into 10^ format
    d <- gsub("^(.*)e\\+", "'\\1'%*%10^", d)
    # convert 0x10^00 to 0
    d <- gsub("\\'0[\\.0]*\\'(.*)", "'0'", d)
    # return this as an expression
    parse(text=d)
}

x_max=max(dat$V1)
y_max=max(dat$V2)*1.1

# Entire plot
ggplot(data=dat, aes(x=V1, y=V2)) +
    geom_line() +
    theme_bw() +
    xlab("Num. markers in 10kb") +
    ylab("Count") +
    scale_x_continuous(labels = comma) +
    scale_y_continuous(labels=fancy_scientific) +
    coord_cartesian(xlim=c(0,x_max), ylim=c(0,y_max))

# Subset to y max 150
y_max=150
ggplot(data=dat, aes(x=V1, y=V2)) +
    geom_line(size = 0.2) +
    theme_bw() +
    theme(
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 6, face = "bold")) +
    xlab("Num. markers in 10kb") +
    ylab("Count") +
    scale_x_continuous(labels = comma) +
    scale_y_continuous(labels=fancy_scientific) +
    coord_cartesian(xlim=c(0,x_max), ylim=c(0,y_max))
ggsave("ED_Fig3c_single_count_max150.jpg", width=3, height = 2.5)

y_max=max(dat$V3)*1.1
ggplot(data=dat, aes(x=V1, y=V3)) +
    geom_line(size = 0.2) +
    theme_bw() +
    theme(
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 6, face = "bold")) +
    xlab("Num. markers in 10kb") +
    ylab("Cumulative Bases") +
    scale_x_continuous(labels = comma) +
    scale_y_continuous(labels=fancy_scientific) +
    coord_cartesian(xlim=c(0,x_max), ylim=c(0,y_max))
ggsave("ED_Fig3d_single_cum.jpg", width=3, height = 2.5)


