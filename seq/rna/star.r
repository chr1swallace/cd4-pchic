library(ggplot2)
library(reshape2)

alnplot <- function(d, file, ylab) {
  d$id <- rownames(d)
  d$id <- factor(c("activated-1", "nonactivated-1", "activated-2", "nonactivated-2"), levels=c("nonactivated-1","activated-1","nonactivated-2","activated-2"))
  m <- melt(d)
  pdf(file)
  print(ggplot(m, aes(x=id, y=value, fill=variable)) + xlab('Sample') + ylab(ylab) + guides(fill=guide_legend(title="Mapping class")) + geom_bar(stat = "identity", position="fill"))
  dev.off()
}

d <- read.csv('data_aln.csv')
rownames(d) <- factor(c("activated-1", "nonactivated-1", "activated-2", "nonactivated-2"), levels=c("nonactivated-1","activated-1","nonactivated-2","activated-2"))
alnplot(as.data.frame(data.matrix(d)[,1] * data.matrix(d)[,2:7]), 'data_aln.cnt.pdf', 'Read ratio')
