## self contained

## create right hand panels on IL2RA ASE figure given supplementary data items
## supp-data-IL2RA-ASE.csv
## supp-data-IL2RA-quant.csv

library(ggplot2)
library(cowplot); theme_set(theme_cowplot())
library(data.table)
library(magrittr)

setwd(CD4CHIC.EXTDATA)
x <- fread("supp-data-IL2RA-ASE-sep2016.csv")

## drop failed sample with v low reads
x <- x[(A1+A2)>100,]
w <- 10
x[,stim:=factor(stim,levels=c("genomic","time0","non","act"))]
levels(x$stim) <- c("genomic","time0","non-activated","activated")
x$rev <- !(x$snp %in% c("rs12244380","rs61839660"))
x[,snp:=factor(snp,levels=c("rs61839660","rs12722495","rs12244380"))]
levels(x$snp) <- c("rs61839660 (intron 1)","rs12722495 (intron 1)","rs12244380 (3' UTR)")
library(RColorBrewer)
cols <- brewer.pal(4,"Set1")
cols[1] <- "grey30"
x[,logprop:=log(A2/(A1+A2))]
x[,logratio:=log(A2/A1)]
## x[snp %in% c("rs12722495","rs7909519","rs12722522","rs2104286"),logprop:=log(A2/(A1+A2))]
x[,mean.logprop:=mean(logprop), by=c("time","stim","snp","Expt")]
x[,mean.logratio:=mean(logratio), by=c("time","stim","snp","Expt")]
x <- x[order(Expt,time,snp),]
x$cell <- sub("central memory","Central Memory CD4+",x$cell)
## basic qc
x[stim=="genomic",exp(mean(logprop)),by="snp"]
## should be close to 0.5
## rs7909519 is 0.46 - ok?
x[stim=="genomic",exp(mean(logprop)),by=c("snp","Expt")]
ggplot(x, aes(x=log2(A1+A2))) + geom_histogram(binwidth=0.1)
#x <- x[log2(A1+A2)>9,]

qc <- x[,.(totalreads=sum(A1+A2),
           mean.logprop=mean(logprop),
           sd.logprop=sd(logprop)),
        by=c("time","stim","snp","Expt")]
qc[,cov:=sd.logprop/abs(mean.logprop)]
ggplot(qc, aes(x=log2(totalreads),y=cov)) + geom_point() + facet_wrap(~snp)

## one outlier by cov
qc[cov>0.25,]
x <- x[!(stim=="non-activated" & snp=="rs61839660 (intron 1)" & Expt=="D135" & time==120),]
## recalculate means
x[,mean.logprop:=mean(logprop), by=c("time","stim","snp","Expt")]
x[,mean.logratio:=mean(logratio), by=c("time","stim","snp","Expt")]


################################################################################
## PAIRED TESTS WITHIN INDIVIDUALS, Total CD4 only
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

## PAIRED TESTS WITHIN INDIVIDUALS, Memory CD4 only
MRESULTS <- unique(x[cell=="Central Memory CD4+",.(snp,stim,time,Expt,rev)])
MRESULTS[,c("estimate","stat","p"):=list(as.numeric(NA),as.numeric(NA),as.numeric(NA))]
## average genomic ratios within individuals
myf <- function(y0,y) {
    if(length(y0)<2 | length(y)<2)
        return(list(estimate=NA,stat=NA,p=NA))
    T <- t.test(y0,y)
    ## list(estimate=diff(T$estimate),stat=T$statistic,p=T$p.value) 
    W <- wilcox.test(y0,y)
    list(estimate=diff(T$estimate),stat=W$statistic,p=W$p.value)
}
for(i in 1:nrow(MRESULTS)) {
    s <- MRESULTS[i,]$stim %>% as.character()
    tm <- MRESULTS[i,]$time
    ex <- MRESULTS[i,]$Expt
    sn <- MRESULTS[i,]$snp
    ctl <- unique(with(x[stim =="genomic" & snp==sn & cell=="Central Memory CD4+",],mean.logprop))
    ret <- with(x[stim==s & time==tm & Expt==ex & snp==sn,],
                myf(ctl,logprop))
    MRESULTS[i, c("estimate","stat","p"):=list(ret$estimate,ret$stat,ret$p)]
}
MRESULTS[stim=="time0",]
MRESULTS[stim=="activated",]

## PAIRED TESTS WITHIN INDIVIDUALS, Memory CD4 only, time 0 baseline
tmp1 <- unique(x[cell=="Central Memory CD4+" & stim=="time0",.(snp,stim,time,Expt,rev,mean.logratio)])
tmp2 <- unique(x[cell=="Central Memory CD4+" & stim=="activated",.(snp,stim,time,Expt,rev,mean.logratio)])
tmp <- merge(tmp1,tmp2,by=c("snp","Expt","rev"))
for(s in unique(tmp$snp)) {
    print(with(tmp[snp==s,], t.test(mean.logratio.x, mean.logratio.y,paired=TRUE)))
    print(with(tmp[snp==s,], wilcox.test(mean.logratio.x, mean.logratio.y,paired=TRUE)))
}

################################################################################

## UNPAIRED TESTS ACROSS INDIVIDUALS

ux <- unique(x,by=c("Expt","stim","snp","time"))
ux[snp=="rs12722495 (intron 1)",mean.logratio:=-mean.logratio]
R <- unique(x[stim!="genomic",.(snp,stim,time,cell)])
R$p <- as.numeric(NA)

for(i in 1:nrow(R)) {
    tmp0 <- ux[ux$stim=="genomic" & ux$snp==R$snp[i],]
    ctl <- with(tmp0, logratio)
    tmp <- ux[ux$cell==R[i]$cell & ux$stim==R$stim[i] & ux$snp==R$snp[i] & ux$time==R$time[i],]
    test <- with(tmp, ifelse(rev,mean(tmp0$logratio)-logratio,logratio))
#    test <- with(tmp, logratio)
    tt <- wilcox.test(ctl,test)
    R[i,p:=tt$p.value]
}
R
    
ux[,mean.logprop.global:=mean(ifelse(rev,-mean.logprop,mean.logprop)), by=c("time","stim","snp")]
ux[,mean.logratio.global:=mean(ifelse(rev,-mean.logratio,mean.logratio)), by=c("time","stim","snp")]
## ux[,mean.logprop.global:=mean(mean.logprop), by=c("time","stim","snp")]
## ux[,mean.logratio.global:=mean(mean.logratio), by=c("time","stim","snp")]
ux <- merge(ux, R,by=c("stim","time","snp","cell"),all=TRUE)
ux

ux[,plab:=paste0("p=",format.pval(round(p,3),digits=3))]
ux[,mean.ratio:=exp(ifelse(rev,-mean.logratio,mean.logratio))]
#ux[,mean.ratio:=exp(mean.logratio)]
##here
ux[,mean.ratio.global:=exp(mean.logratio.global)]
ux[,x:=paste(time,stim,sep="\n")]
dput(levels(as.factor(ux$x)))
ux[,x:=factor(x, levels=c("-30\ngenomic", "0\ntime0", "120\nnon-activated", "240\nnon-activated",
"120\nactivated", "240\nactivated"))]
levels(ux$x) <- c("genomic", "time0", "non-act'd\n2 hours", "non-act'd\n4 hours",
"act'd\n2 hours", "act'd\n4 hours")
ux[,xn:=as.numeric(x)]
ux[ux$cell=="Central Memory CD4+",xn:=xn+0.1]

ux2 <- unique(ux,by=c("stim","snp","time","cell"))
ux2[,cell:=paste(cell,"(mean)")]
ux2[,pstar:=""]
ux2[p<0.01,pstar:="*"]
ux2[p<0.001,pstar:="**"]
ux2[p<0.0001,pstar:="***"]
ux2[,ylab:=max(ux$mean.ratio)+0.1]
plotsnp <- function(s) {
    ux <- ux[snp==s,]
    ## ux2 <- ux2[snp==s,]
ggplot(ux, aes(x=xn,y=mean.ratio,col=stim,pch=cell)) + 
geom_hline(aes(yintercept=mean.ratio.global),data=ux[stim=="genomic",],
           col=cols[1],linetype="dashed") +
geom_point()  +
geom_point(aes(y=mean.ratio.global),data=ux2[ux2$snp==s,],size=3)  +
geom_text(aes(x=xn,y=ylab,label=pstar), data=ux2[!is.na(p) & ux2$snp==s,],show.legend=FALSE) +
## geom_text(aes(y=mean.ratio.global,label=plab),
##           data=ux2[!is.na(p) & ux2$snp==s,],nudge_x=-0.4,nudge_y=0,size=4,
##           show.legend=FALSE) +
background_grid() +
scale_colour_manual("Condition",values=cols,labels=c("gDNA","cDNA","cDNA (n)","cDNA (a)"),guide=FALSE) +
scale_shape_manual("Cell type",values=c(1,19,22,15)) +
    ## "1"="Central Memory CD4+","19"="Central Memory CD4+ (mean)",
    ##                           "22"="Total CD4+","15"="Total CD4+ (mean)")) +
scale_y_log10("Allelic ratio (log scale)",breaks=c(1,1.25,1.5,2)) +
scale_x_continuous("Condition",breaks=1:6,labels=levels(ux$x)) +
ggtitle(s) + 
theme(legend.position=c(0.9,0.6))
}
snps <- unique(ux$snp)
plotsnp(snps[[1]])
plotsnp(snps[[2]])
plotsnp(snps[[3]])


plotcmem <- function(s) {
ggplot(ux[snp==s & cell=="Central Memory CD4+" & x!="genomic",], aes(x=x, y=mean.ratio,group=Expt)) + geom_point() + geom_path(colour="gray60") + geom_point(aes(y=mean.ratio.global),data=ux2[ux2$snp==s & cell=="Central Memory CD4+ (mean)" & x!="genomic",],size=30,pch="-") +
scale_y_log10("Allelic ratio (log scale)",breaks=c(0.6,0.8,1.0,1.2,1.4)) +
geom_hline(aes(yintercept=mean.ratio.global),data=ux[stim=="genomic" & snp==s,],
           col=cols[1],linetype="dashed") +
scale_x_discrete("Condition")
}
plotcmem(snps[[2]])
plotcmem(snps[[3]])

library(cowplot)
ggplot(ux[cell=="Central Memory CD4+" & x!="genomic",], aes(x=x, y=mean.ratio,group=Expt)) +
geom_path(colour="gray60") +
geom_point() +
#geom_point(aes(y=mean.ratio.global),data=ux2[cell=="Central Memory CD4+ (mean)" & x!="genomic",],size=30,pch="-") +
scale_y_log10("Allelic ratio (log scale)",breaks=c(0.6,0.8,1.0,1.2,1.4)) +
geom_hline(aes(yintercept=mean.ratio.global),data=ux[stim=="genomic" & snp %in% snps[2:3],],
           col=cols[1],linetype="dashed") +
scale_x_discrete("Condition",labels=c("not act'd", "act'd (4h)")) +
theme(strip.background=element_blank()) +
facet_grid(.~snp)

f <- file.path(CD4CHIC.OUT,"paper","figure-ase-centralmem.png")
ggsave(f,width=5.5,height=4)
system(paste("display ",f))

scale_shape_manual("Cell type",values=c(1,19,22,15))

geom_path(aes(x=time,y=exp(mean.logratio)),col=cols[3],data=ux2[stim %in% c("time0","non-activated"),],linetype=3) +
geom_path(aes(x=time,y=exp(mean.logratio)),col=cols[4],data=ux2[stim %in% c("time0","activated"),],linetype=3)+
geom_segment(aes(x=time-w,xend=time+w,y=exp(mean.logratio),yend=exp(mean.logratio)),size=3,data=ux) +
geom_point() + geom_smooth() +
geom_text(aes(x=time-w*4,y=exp(mean.logratio),
              label=plab),
          data=ux[same.dir==TRUE,],nudge_x=-2,nudge_y=0.01,size=4) +
scale_colour_manual("Condition",values=cols,labels=c("gDNA","cDNA","cDNA (n)","cDNA (a)")) +
scale_x_continuous("Time (hours) for cDNA samples",breaks=c(0,120,240),labels=c("0","2","4"),limits=c(-10,250)
                   ,expand=c(0.2,0.05)) +
facet_grid(snp~.) +
scale_y_continuous("Proportion of reads carrying allele A2") +
theme(legend.position="bottom") +
background_grid()


ggplot(ux, aes(x=time,y=exp(logratio),col=stim)) + 
geom_hline(aes(yintercept=exp(mean.logprop)),data=ux2[stim=="genomic",],col=cols[1],linetype="dashed") +
geom_path(aes(x=time,y=exp(mean.logprop)),col=cols[3],data=ux2[stim %in% c("time0","non-activated"),],linetype=3) +
geom_path(aes(x=time,y=exp(mean.logprop)),col=cols[4],data=ux2[stim %in% c("time0","activated"),],linetype=3)+
geom_segment(aes(x=time-w,xend=time+w,y=exp(mean.logprop),yend=exp(mean.logprop)),size=3,data=ux) +
geom_point() + #geom_smooth() +
geom_text(aes(x=time-w*4,y=exp(mean.logprop),
              label=plab),
          data=ux[same.dir==TRUE,],nudge_x=-2,nudge_y=0.01,size=4) +
scale_colour_manual("Condition",values=cols,labels=c("gDNA","cDNA","cDNA (n)","cDNA (a)")) +
scale_x_continuous("Time (hours) for cDNA samples",breaks=c(0,120,240),labels=c("0","2","4")## ,limits=c(-10,250)
                   ,expand=c(0.2,0.05)) +
facet_grid(snp~.) +
scale_y_continuous("Proportion of reads carrying allele A2") +
theme(legend.position="bottom") +
background_grid()

f <- file.path(CD4CHIC.OUT,"paper","figure-ASE-grouped.png")
ggsave(f,width=5.5,height=8.5)
system(paste("display ",f))

################################################################################

ux2[,x:=paste(time,stim,sep="\n")]
dput(levels(as.factor(ux2$x)))
ux2[,x:=factor(x, levels=c("-30\ngenomic", "0\ntime0", "120\nnon-activated", "240\nnon-activated",
"120\nactivated", "240\nactivated"))]
levels(ux2$x) <- c("genomic", "time0", "non-act'd\n2 hours", "non-act'd\n4 hours",
"act'd\n2 hours", "act'd\n4 hours")
ux2 <- merge(ux2,ux[,.(stim,time,snp,t.p,same.dir,plab)],
             by=c("stim","time","snp"))


library(colorspace)

cols <- c("#333333", "#276EA8", "#3D9F3A", "#883E93")
plotsnp <- function(ux2) {

    ux2 <- ux2[order(snp,x,Expt),]
lp2r <- function(lp) {
    p <- exp(lp)
    p/(1-p)
}
ux2[,mean.logratio:=log2(lp2r(mean.logprop))]
ux2[,mean.logratio.byindiv:=log2(lp2r(mean.logprop.byindiv))]
ux2[,mean.ratio:=lp2r(mean.logprop)]
ux2[,mean.ratio.byindiv:=lp2r(mean.logprop.byindiv)]

ggplot(ux2, aes(x=x,col=stim)) +
geom_hline(aes(yintercept=mean.ratio),data=ux2[stim=="genomic",],col="grey") +
geom_point(aes(y=mean.ratio),pch=19) +
geom_point(aes(y=mean.ratio.byindiv),pch=4) +
## geom_path(aes(y=mean.ratio.byindiv,group=Expt),linetype="dashed",col="grey") +
## geom_path(aes(y=mean.ratio2,group=snp),linetype="dashed") +
geom_text(aes(y=mean.ratio,label=plab),
          data=ux2[!is.na(t.p),],nudge_x=-0.4,nudge_y=0,size=4,
          show.legend=FALSE) +
background_grid() +
scale_colour_manual("Condition",values=cols,labels=c("gDNA","cDNA","cDNA (n)","cDNA (a)")) +
scale_y_log10("Allelic ratio (log scale)",breaks=c(1,1.25,1.5,2)) +
scale_x_discrete("Condition") +
theme(legend.position="none") 
#facet_grid(snp ~ .) 

}



bak <- ux2

for(isnp in levels(bak$snp)) {
    plotsnp(ux2[ux2$snp==isnp,])
    f <- paste0("figure-",sub(" .*","",isnp),".png")
    f <- file.path(CD4CHIC.OUT,"paper",f)
    ggsave(f,width=8.5,height=3.5)
    f <- sub("fig","dispfig",f)
    ggsave(f,width=6.5,height=3.5)
#    system(paste("display ",f))
}

