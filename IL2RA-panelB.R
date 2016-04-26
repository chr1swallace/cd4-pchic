## self contained

## create right hand panels on IL2RA ASE figure given supplementary data items
## supp-data-IL2RA-ASE.csv
## supp-data-IL2RA-quant.csv

library(ggplot2)
library(cowplot)
library(data.table)
library(magrittr)

x <- fread(file.path(CD4CHIC.OUT,"paper/supp-data-IL2RA-ASE.csv"))

## drop failed sample with v low reads
x <- x[(T+C)>100,]

x[,N:=C+T]
x[,Y:=T/N]

w <- 10


x[,stim:=factor(stim,levels=c("genomic","time0","non","act"))]
levels(x$stim) <- c("genomic","time0","non-activated","activated")

theme_set(theme_cowplot())

library(RColorBrewer)
cols <- brewer.pal(4,"Set1")
cols[1] <- "grey30"
x[,mean.logratio:=mean(log(T/(T+C))), by=c("time","stim","Expt")]
x <- x[order(Expt,time),]

################################################################################

## TESTS WITHIN INDIVIDUALS

myf <- function(y,by) {
    yg<-y[by=="genomic"]
    yc<-y[by!="genomic"]
    if(length(yg)<2 | length(yc)<2)
        return(list(estimate=NA,stat=NA,p=NA))
    T <- t.test(log(yg),log(yc))
    list(estimate=diff(T$estimate),stat=T$statistic,p=T$p.value) 
}

RESULTS <- unique(x[stim!="genomic",.(stim,time,Expt)])
RESULTS[,c("estimate","stat","p"):=list(as.numeric(NA),as.numeric(NA),as.numeric(NA))]
for(i in 1:nrow(RESULTS)) {
    s <- RESULTS[i,]$stim %>% as.character()
    tm <- RESULTS[i,]$time
    ex <- RESULTS[i,]$Expt
    ret <- with(x[stim %in% c("genomic",s) & time %in% c(-30,tm) & Expt==ex,],
                myf(Y,stim))
    RESULTS[i, c("estimate","stat","p"):=list(ret$estimate,ret$stat,ret$p)]
}
RESULTS[stim=="time0",]
RESULTS[stim=="non-activated",]
RESULTS[stim=="activated",]

## SUPPLEMENTARY FIGURE

p <- ggplot(x, aes(x=time,y=T/(T+C),col=stim)) +
geom_hline(aes(yintercept=exp(mean.logratio)),data=x[time<0,],col=cols[1],linetype="dashed") +
geom_point() + facet_wrap(~ Expt) + #geom_smooth() +
geom_path(aes(x=time,y=exp(mean.logratio)),col=cols[3],data=x[stim %in% c("time0","non-activated"),]) +
geom_path(aes(x=time,y=exp(mean.logratio)),col=cols[4],data=x[stim %in% c("time0","activated"),]) +
geom_segment(aes(x=time-w,xend=time+w,y=exp(mean.logratio),yend=exp(mean.logratio)),size=3) +
scale_colour_manual("Condition",values=cols,labels=c("gDNA","cDNA","cDNA (non)","cDNA (act)")) +
scale_x_continuous("Time (hours) for cDNA samples",breaks=c(0,120,240),labels=c("0","2","4")) +
scale_y_continuous("Proportion of reads carrying T") + 
theme(legend.position="bottom")
p

ggsave(file.path(CD4CHIC.OUT,"paper/figure-ASE-individuals.png"),width=6,height=8)

################################################################################

## TESTS ACROSS INDIVIDUALS

y <- x
y[,medy:=median(T/(T+C)),by=c("Expt","time","stim")]
y <- unique(y,by=c("Expt","time","stim"))
setkey(y,Expt,time,stim)

wilcox.test(y[stim=="genomic",][["medy"]],
            y[stim=="time0",][["medy"]])
wilcox.test(y[stim=="genomic",][["medy"]],
            y[time==120 & stim=="activated",][["medy"]])
wilcox.test(y[stim=="genomic",][["medy"]],
            y[time==120 & stim=="non-activated",][["medy"]])
wilcox.test(y[stim=="genomic",][["medy"]],
            y[time==240 & stim=="activated",][["medy"]])
wilcox.test(y[stim=="genomic",][["medy"]],
            y[time==240 & stim=="non-activated",][["medy"]])

x <- x[order(time),]
x[,mean.logratio:=mean(log(T/(T+C))), by=c("time","stim")]
ggplot(x, aes(x=time,y=T/(T+C),col=stim)) + geom_point() + #geom_smooth() +
geom_hline(aes(yintercept=exp(mean.logratio)),data=x[stim=="genomic",],col=cols[1],linetype="dashed") +
geom_path(aes(x=time,y=exp(mean.logratio)),col=cols[3],data=x[stim %in% c("time0","non-activated"),]) +
geom_path(aes(x=time,y=exp(mean.logratio)),col=cols[4],data=x[stim %in% c("time0","activated"),])+
geom_segment(aes(x=time-w,xend=time+w,y=exp(mean.logratio),yend=exp(mean.logratio)),size=3,data=x) +
scale_colour_manual("Condition",values=cols,labels=c("gDNA","cDNA","cDNA (n)","cDNA (a)")) +
scale_x_continuous("Time (hours) for cDNA samples",breaks=c(0,120,240),labels=c("0","2","4")## ,limits=c(-10,250)
                   ) +
scale_y_continuous("Proportion of reads carrying T") + 
theme(legend.position="bottom")

f <- file.path(CD4CHIC.OUT,"paper","figure-ASE-grouped.png")
ggsave(f,width=5,height=5)
system(paste("display ",f))

################################################################################

## QUANTIFICATION

z <- fread(file.path(CD4CHIC.OUT,"paper/supp-data-IL2RA-quant.csv"))
z <- merge(z,y[,.(stim,Expt,time,medy)], by=c("stim","Expt","time"))
z[,Tconc:=log2(conc * medy)]
z[,Cconc:=log2(conc * (1-medy))]
z[,c("meanC","meanT"):=list(mean(Cconc), mean(Tconc)),
  by=c("time","stim")]
z <- melt(z,measure.vars=list(c("Tconc","Cconc"), logconc=c("meanT","meanC")))
z[,c("Tconc","Cconc"):=list(exp(value1),exp(value2))]

## duplicate time 0 for graphing
z <- rbind(z,z[time==0,])
z[time==0,stim:=ifelse(duplicated(Expt),"non-activated","activated")]

z <- z[order(time),]
ggplot(z, aes(x=time,colour=stim,pch=variable,linetype=variable)) +
geom_point(aes(y=Tconc),size=2) +
geom_point(aes(y=Cconc),size=2) +
geom_path(aes(y=Cconc,group=paste(stim,variable))) +
geom_path(aes(y=Tconc,group=paste(Expt,stim,variable)),size=0.4,alpha=0.5) +
scale_y_log10("Concentration by allele (ng/µL, log scale)") +
scale_x_continuous("Time (hours)",breaks=c(0,120,240),labels=c("0","2","4"),expand=c(0.05,0.05)) +
scale_colour_manual("Condition",values=cols[-c(1:2)],labels=c("cDNA (act)","cDNA (non)")) +
scale_shape_discrete("Allele",labels=c("T","C")) +
scale_linetype_discrete("Allele",labels=c("T","C")) +
theme(legend.position="bottom")

f <- file.path(CD4CHIC.OUT,"paper","figure-ASE-quant.png")
ggsave(file=f,width=5,height=5)
system(paste("display",f))
