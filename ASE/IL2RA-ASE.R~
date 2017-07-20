## self contained

## create right hand panels on IL2RA ASE figure given supplementary data items
## supp-data-IL2RA-ASE.csv
## supp-data-IL2RA-quant.csv

library(ggplot2)
library(cowplot); theme_set(theme_cowplot())
library(data.table)
library(magrittr)

setwd(CD4CHIC.EXTDATA)
x <- fread("supp-data-10-ASE.csv")

## drop failed sample with v low reads
x[(A1+A2)<=2000,]
x <- x[(A1+A2)>100,]
x <- x[Expt!="D181" | (A1+A2>2000),] # one outlier for D181 from other replicates with low reads in each time/stim condition.
w <- 10
x[,stim:=factor(stim,levels=c("genomic","time0","non","act"))]
levels(x$stim) <- c("genomic","time0","non-activated","activated")
x[,snp:=factor(snp,levels=c("rs61839660","rs12722495","rs12244380"))]
levels(x$snp) <- c("rs61839660 (intron 1)","rs12722495 (intron 1)","rs12244380 (3' UTR)")

library(RColorBrewer)
cols <- brewer.pal(4,"Set1")
cols[1] <- "black"
x[,logprop:=log(A2/(A1+A2))]
x[,logratio:=log(A2/A1)]
x[,mean.logprop:=mean(logprop), by=c("time","stim","snp","Expt")]
x[,mean.logratio:=mean(logratio), by=c("time","stim","snp","Expt")]
x <- x[order(Expt,time,snp),]

## basic qc
x[stim=="genomic",exp(mean(logprop)),by="snp"]
x[stim=="genomic",exp(mean(logprop)),by=c("snp","Expt")]
ggplot(x, aes(x=log2(A1+A2))) + geom_histogram(binwidth=0.1)

qc <- x[,.(totalreads=sum(A1+A2),
           minreads=min(A1+A2),
           mean.logprop=mean(logprop),
           sd.logprop=sd(logprop)),
        by=c("time","stim","snp","Expt")]
qc[,cov:=sd.logprop/abs(mean.logprop)]
ggplot(qc, aes(x=log10(totalreads),y=cov)) + geom_point() + facet_wrap(~snp) + geom_hline(yintercept=c(0.2,0.25,0.3),col="violet")
ggplot(qc, aes(x=log10(minreads),y=cov)) + geom_point() + facet_wrap(~snp) + geom_hline(yintercept=c(0.2,0.25,0.3),col="violet")

qc[cov>0.2,]

x <- x[!(stim=="non-activated" & snp=="rs61839660 (intron 1)" & Expt=="D135" & time==120),]
## recalculate means
x[,mean.logprop:=mean(logprop), by=c("time","stim","snp","Expt")]
x[,mean.logratio:=mean(logratio), by=c("time","stim","snp","Expt")]
x <- x[!(time %in% c(30,60)),]

################################################################################
## PAIRED TESTS WITHIN INDIVIDUALS (WHERE gDNA exists for same indiv)
RESULTS <- unique(x[cell=="Total CD4+",.(snp,stim,time,Expt)])
RESULTS[,c("estimate","stat","p"):=list(as.numeric(NA),as.numeric(NA),as.numeric(NA))]
myf <- function(y,by) {
    yg<-y[by=="genomic"]
    yc<-y[by!="genomic"]
    if(length(yg)<2 | length(yc)<2)
        return(list(estimate=NA,stat=NA,p=NA))
    T <- t.test(yg,yc)
    ## list(estimate=diff(T$estimate),stat=T$statistic,p=T$p.value) 
    W <- wilcox.test(yg,yc)
    list(estimate=diff(T$estimate),stat=W$statistic,p=W$p.value)
}
for(i in 1:nrow(RESULTS)) {
    s <- RESULTS[i,]$stim %>% as.character()
    tm <- RESULTS[i,]$time
    ex <- RESULTS[i,]$Expt
    sn <- RESULTS[i,]$snp
    ret <- with(x[stim %in% c("genomic",s) & time %in% c(-30,tm) & Expt==ex & snp==sn,],
                myf(logprop,stim))
    RESULTS[i, c("estimate","stat","p"):=list(ret$estimate,ret$stat,ret$p)]
}
RESULTS[stim=="time0",]
RESULTS[stim=="non-activated",]
RESULTS[stim=="activated",]


################################################################################

## UNPAIRED TESTS ACROSS INDIVIDUALS

ux <- unique(x,by=c("Expt","stim","snp","time"))
R <- unique(x[stim!="genomic",.(snp,stim,time)])
R$p <- as.numeric(NA)
for(i in 1:nrow(R)) {
    tmp0 <- ux[ux$stim=="genomic" & ux$snp==R$snp[i],]
    tmp <- ux[ux$stim==R$stim[i] & ux$snp==R$snp[i] & ux$time==R$time[i],]
    tt <- wilcox.test(tmp0$mean.logratio, tmp$mean.logratio)
    R[i,p:=tt$p.value]
}
R
    
ux[,mean.logprop.global:=mean(mean.logprop), by=c("time","stim","snp")]
ux[,mean.logratio.global:=mean(mean.logratio), by=c("time","stim","snp")]
ux <- merge(ux, R,by=c("stim","time","snp"),all=TRUE)
tail(ux)

ux[,plab:=paste0("p=",format.pval(round(p,3),digits=3))]
ux[,mean.ratio:=exp(mean.logratio)]
ux[,mean.ratio.global:=exp(mean.logratio.global)]
ux[,x:=paste(time,stim,sep="\n")]
dput(levels(as.factor(ux$x)))
ux[,x:=factor(x, levels=c("-30\ngenomic", "0\ntime0", "120\nnon-activated", "240\nnon-activated",
"120\nactivated", "240\nactivated"))]
levels(ux$x) <- c("gDNA", "time0", "non\n2 hrs", "non\n4 hrs",
"act\n2 hrs", "act\n4 hrs")
ux[,xn:=as.numeric(x)]

ux2 <- unique(ux,by=c("stim","snp","time"))
ux2[,cell:=paste(cell,"(mean)")]
ux2[,pstar:=""]
ux2[p<0.01,pstar:="*"]
ux2[p<0.001,pstar:="**"]
ux2[p<0.0001,pstar:="***"]
ux2[,ylab:=max(ux$mean.ratio)+0.1]
plotsnp <- function(s) {
    ux <- ux[snp==s,]
    ## ux2 <- ux2[snp==s,]
ggplot(ux, aes(x=xn,y=mean.ratio,col=stim)) + 
geom_hline(aes(yintercept=mean.ratio.global),data=ux[stim=="genomic",],
           col="grey50",linetype="dashed") +
geom_point(aes(y=mean.ratio.global),data=ux2[ux2$snp==s,],size=3,pch=3)  +
geom_point(pch=4)  +
geom_text(aes(x=xn,y=ylab,label=pstar), data=ux2[!is.na(p) & ux2$snp==s,],show.legend=FALSE) +
geom_text(aes(y=mean.ratio.global,label=plab),
          data=ux2[!is.na(p) & ux2$snp==s,],nudge_x=-0.4,nudge_y=0,size=4,
          show.legend=FALSE) +
background_grid() +
scale_colour_manual("Condition",values=cols,labels=c("gDNA","cDNA","cDNA (n)","cDNA (a)"),guide=FALSE) +
scale_shape_manual("Cell type",values=c(1,19,22,15)) +
    ## "1"="Central Memory CD4+","19"="Central Memory CD4+ (mean)",
    ##                           "22"="Total CD4+","15"="Total CD4+ (mean)")) +
scale_y_log10("Allelic ratio (log scale)",breaks=c(1,1.25,1.5,2)) +
scale_x_continuous("Condition",breaks=1:6,labels=levels(ux$x)) +
theme(legend.position="none") #c(0.9,0.6)
}
tplotsnp <- function(s) {
    ux <- ux[snp==s,]
    ## ux2 <- ux2[snp==s,]
ggplot(ux, aes(y=xn,x=mean.ratio,col=stim)) + 
geom_vline(aes(xintercept=mean.ratio.global),data=ux[stim=="genomic",],
           col="grey50",linetype="dashed") +
geom_point(aes(x=mean.ratio.global),data=ux2[ux2$snp==s,],size=3,pch=3)  +
geom_point(pch=4)  +
geom_text(aes(y=xn,x=ylab,label=pstar), data=ux2[!is.na(p) & ux2$snp==s,],show.legend=FALSE) +
geom_text(aes(x=mean.ratio.global,label=plab),
          data=ux2[!is.na(p) & ux2$snp==s,],nudge_y=0.2,nudge_x=0,size=4,
          show.legend=FALSE) +
background_grid() +
scale_colour_manual("Condition",values=cols,labels=c("gDNA","cDNA","cDNA (n)","cDNA (a)"),guide=FALSE) +
scale_shape_manual("Cell type",values=c(1,19,22,15)) +
    ## "1"="Central Memory CD4+","19"="Central Memory CD4+ (mean)",
    ##                           "22"="Total CD4+","15"="Total CD4+ (mean)")) +
scale_x_log10("Allelic ratio (log scale)",breaks=c(1,1.25,1.5,2),expand=c(0.05,0.05)) +
scale_y_reverse("Condition",breaks=1:6,labels=levels(ux$x)) +
theme(legend.position="none") #c(0.9,0.6)
}

snps <- unique(ux$snp)
tplotsnp(snps[[1]]) + ggtitle(snps[[1]])
tplotsnp(snps[[2]]) + ggtitle(snps[[2]])
tplotsnp(snps[[3]]) + ggtitle(snps[[3]])

w=3/7
for(i in 1:3) {
    tplotsnp(snps[[i]])
    f <- file.path(CD4CHIC.OUT,"paper",paste0("figure-ASE-",sub(" .*","",snps[[i]]),".pdf"))
    #ggsave(f,width=8.5*w,height=3.5*w)
    ggsave(f,width=7.8*w,height=5)
#    system(paste("display",f))
}

unique(x[,.(snp,alleles)])
## for Olly's poster
i=2; 
plotsnp(snps[[i]])
f <- file.path("/scratch/wallace",paste0("figure-ASE-disp-",sub(" .*","",snps[[i]]),".pdf"))
ggsave(f,width=27,height=8,units="cm")
system(paste("display",f))
f

ux[time<0,time:=NA]
tmp <- ux[,.(N=.N),by=c("time","stim","snp")]
write.table(tmp,file=file.path(CD4CHIC.EXTDATA,"supp-table-12-ASE-sample-counts.csv"))   
