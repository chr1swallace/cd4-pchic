library(ggplot2)
library(plyr)

S <- 0.005
H <- 5
E <- 100

hic <- read.delim("hic.csv")
oeg <- read.delim("oeg.csv")
rna <- read.delim("rna.csv")
rna$padj[is.na(rna$padj)] <- 1

print(nrow(rna))
print(nrow(rna[rna$padj < S,]))
print(nrow(rna[rna$padj < S & rna$log2FoldChange > 0,]))
print(nrow(rna[rna$padj < S & rna$log2FoldChange < 0,]))
print(table(rna[rna$padj < S, 'biotype']))

a <- (rna$activated.1 + rna$activated.2) / 2 >= E
b <- (rna$nonactivated.1 + rna$nonactivated.2) / 2 >= E
print(sum(a)) # expressed act
print(sum(b)) # expressed nonact
print(sum(a | b)) # expressed union

print(nrow(hic))
c <- hic$activated >= H & hic$nonactivated < H
d <- hic$activated < H & hic$nonactivated >= H
print(nrow(hic[c | d,])) # diff contact
print(nrow(hic[c,])) # upreg
print(nrow(hic[d,])) # downreg

e <- hic$activated >= H
f <- hic$nonactivated >= H
print(sum(e))
print(sum(f))
print(sum(e | f))

png("1.png")
ggplot(hic, aes(log2(nonactivated + 0.1), log2(activated + 0.1))) + geom_point() + geom_rug(alpha=0.1)
dev.off()

all <- merge(merge(hic, oeg), rna)
sig <- all[all$padj < S,]
sgu <- unique(sig[,colnames(hic)])

png("2.png")
ggplot(sgu, aes(log2(nonactivated + 0.1), log2(activated + 0.1))) + geom_point() + geom_rug(alpha=0.1)
dev.off()

act <- ddply(hic[hic$activated >= H,], 'oeID', function(x) data.frame(nactivated=nrow(x)))
nct <- ddply(hic[hic$nonactivated >= H,], 'oeID', function(x) data.frame(nnonactivated=nrow(x)))
nin <- merge(act, nct, all=T)
nin$nactivated[is.na(nin$nactivated)] <- 0
nin$nnonactivated[is.na(nin$nnonactivated)] <- 0

png("3.png")
ggplot(nin, aes(nnonactivated, nactivated)) + geom_point() + geom_jitter() + geom_rug(alpha=0.1)
dev.off()

nnn <- merge(merge(nin, oeg), rna)
snn <- nnn[nnn$padj < S,]
unn <- unique(snn[,colnames(nin)])

png("4.png")
ggplot(unn, aes(nnonactivated, nactivated)) + geom_point() + geom_jitter() + geom_rug(alpha=0.1)
dev.off()
