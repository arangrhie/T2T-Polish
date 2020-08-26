path<-dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)

na.zero <- function (x) {
  x[is.na(x)] <- 0
  return(x)
}

histfile <- "mBalMus1.k21.hist"
kmer_prof <- read.csv(file=histfile,sep="\t", header=FALSE)

x = kmer_prof[[1]]
y = kmer_prof[[2]]

model_sum=summary(model_4peaks[[1]])

kcovfloor = floor(min_max(model_sum$coefficients['kmercov',])[[1]])
het  = min_max(model_sum$coefficients['r',])
dups = min_max(model_sum$coefficients['bias',])
kcov = min_max(model_sum$coefficients['kmercov',])
mlen = min_max(model_sum$coefficients['length',])
md   = min_max(model_sum$coefficients['d',])

amlen = (mlen[1] + mlen[2]) / 2
ahet  = (het[1]  + het[2])  / 2
amd   = (md[1]   + md[2])   / 2
akcov = (kcov[1] + kcov[2]) / 2
adups = (dups[1] + dups[2]) / 2

unique_hist <- (2 * (1 - amd) * (1 - (1 - ahet)^k))                         * dnbinom(x, size = akcov     / adups, mu = akcov)     * amlen +
  ((amd * (1 - (1 - ahet)^k)^2) + (1 - amd) * ((1 - ahet)^k))  * dnbinom(x, size = akcov * 2 / adups, mu = akcov * 2) * amlen 

one_hist <- ((2*(1-amd)*(1-(1-ahet)^k)) + (2*amd*(1-(1-ahet)^k)^2) + (2*amd*((1-ahet)^k)*(1-(1-ahet)^k))) * dnbinom(x, size = akcov   / adups, mu = akcov)     * amlen   
two_hist <- (((1-amd)*((1-ahet)^k)) + (amd*(1-(1-ahet)^k)^2))                                             * dnbinom(x, size = akcov*2 / adups, mu = akcov * 2) * amlen
thr_hist <- (2*amd*((1-ahet)^k)*(1-(1-ahet)^k))                                                           * dnbinom(x, size = akcov*3 / adups, mu = akcov * 3) * amlen
fou_hist <- (amd*(1-ahet)^(2*k))                                                                          * dnbinom(x, size = akcov*4 / adups, mu = akcov * 4) * amlen

total_kmers = sum(as.numeric(x*y))
unique_kmers = sum(as.numeric(x*unique_hist))
total_error_kmers = sum(as.numeric(error_kmers * x[1:error_xcutoff_ind]))
repeat_kmers = total_kmers - unique_kmers - total_error_kmers
repeat_len=repeat_kmers/(2*kcov)

pred=predict(model_4peaks[[1]], newdata=data.frame(x))

## Compute error rate, by counting kmers unexplained by model through first peak
## truncate errors as soon as it goes to zero, dont allow it to go back up
error_xcutoff = kcovfloor
error_xcutoff_ind = which(x==error_xcutoff)

error_kmers = y[1:error_xcutoff_ind] - pred[1:error_xcutoff_ind]
error_kmers = pmax(error_kmers, 1e-10)

repeat_kmers = y[1:length(y)] - pred[1:length(y)]
repeat_kmers <- c(rep(0,akcov*5-1),repeat_kmers[ceiling(akcov*5)[[1]]:length(repeat_kmers)])

zer_hist_PMF <- error_kmers/sum(y)
zer_hist_PMF[zer_hist_PMF<0] <- 0
one_hist_PMF <- one_hist/sum(y)
two_hist_PMF <- two_hist/sum(y)
thr_hist_PMF <- thr_hist/sum(y)
fou_hist_PMF <- fou_hist/sum(y)
rep_hist_PMF <- repeat_kmers/sum(y)

dev.off()
plot.new()

plot(kmer_prof,type="n", main="GenomeScope Profile\n", xlab="Coverage", ylab="Probability",xlim=c(0,200), ylim=c(0,0.001))

lines(zer_hist_PMF, lwd=1, col="black")
lines(one_hist_PMF, col="red")
lines(two_hist_PMF, col="green")
lines(thr_hist_PMF, col="purple")
lines(fou_hist_PMF, col="blue")
lines(rep_hist_PMF, col="orange")

complete <- rep(0,length(x))

peaks<-akcov * 1:60000

abline(v=peaks, col="black", lty=2, lwd=0.1)

start <- 0
step <- 10000

for (i in peaks[peaks >= start & peaks <= start+step]) {
  
  cat("\nfrequency:",y[i])
  cat("\ncopy number:",i)
  cat("\nmean:",i)
  cat("\nstandard deviation:",adups*i/sqrt(i))
  
  #lines(dnorm(x, mean = i, sd = adups*i/sqrt(i), log = FALSE) * (y[i]*i / repeat_kmers) * repeat_len[[1]], col=cols[i], lwd=1, lty=1)
  
  new <- dnorm(x, mean = i, sd = adups*i/sqrt(i), log = FALSE) * (y[i]*i / repeat_kmers)
  
  cat("\nmax:",max(new), ", index:", which.max(new))
  
  #compute the sum distribution 
  complete <- complete + new
  
}

lines(complete, col="red", lwd=1, lty=1)

lookup_table <- NULL

for (i in 1:200){
  
  totalP<-sum(na.zero(zer_hist_PMF[i]),
              na.zero(one_hist_PMF[i]),
              na.zero(two_hist_PMF[i]),
              na.zero(thr_hist_PMF[i]),
              na.zero(fou_hist_PMF[i]))
  
  cat (i," ",na.zero(zer_hist_PMF[i])/totalP,
             na.zero(one_hist_PMF[i])/totalP,
             na.zero(two_hist_PMF[i])/totalP,
             na.zero(thr_hist_PMF[i])/totalP,
             na.zero(fou_hist_PMF[i])/totalP,"\n")
  
  lookup_table<-rbind(lookup_table,c(i,na.zero(zer_hist_PMF[i])/totalP,
       na.zero(one_hist_PMF[i])/totalP,
       na.zero(two_hist_PMF[i])/totalP,
       na.zero(thr_hist_PMF[i])/totalP,
       na.zero(fou_hist_PMF[i])/totalP))
  
}

lookup_table<-data.frame(lookup_table)

colors <- c("black","red","green","purple","blue")

dev.off()

plot(lookup_table$X1,type="n", xlab="Coverage", ylab="Probability",xlim=c(1,200), ylim=c(0,1))

a<-0
for (i in 2:6) {
  a<-a+1
  lines(lookup_table$X1, lookup_table[,i], type="l", lwd=1.5,
        col=colors[a])
}

#homozygous diploid case d=1, r=0
(d*(1-r)^(2*k))  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov * 4) * length

#haploid case d=0, r=0
((1-d)*((1-r)^k)) * dnbinom(x, size = kmercov*2 / bias, mu = kmercov * 2) * length
