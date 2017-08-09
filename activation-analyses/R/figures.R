## * Set up
library(WGCNA)
library(data.table)
library(parallel); options(mc.cores=2)
library(ggplot2); theme_set(theme_bw())
library(magrittr)
library(gridExtra)
library(reshape)
setwd(file.path(CD4CHIC.ROOT,"activation-analyses/R"))

refactor.mods <- function(x) {
    factor(sub("grey","gray",as.character(x)),
           levels=c("black","green","yellow","red","gray","blue","brown","turquoise"))
}

options(stringsAsFactors = FALSE);
enableWGCNAThreads(3)
source("common.R")

mycols <- unique(as.character(modules$module))
names(mycols) <- mycols
mycols[mycols=="NA"] <- "grey20"
mycols["grey"] <- "grey40"
mycols["gray"] <- "grey40"
mycols["yellow"] <- "goldenrod"
mycols["brown"] <- "SaddleBrown"
mycols["grey"] <- "grey40"
mycols["green"] <- "DarkGreen"
mycols["blue"] <- "DarkBlue"
colScale <- scale_colour_manual(name = "module",values = mycols, na.value="purple",guide=FALSE)

#source(file.path(CD4CHIC.ROOT,"activation-analyses/R/collect-data.R"))

## * Load microarray data

(load(file.expression))

## ## modules
spower <- "combined"
geneSummary.file<-out.file('geneSummary')
(load(geneSummary.file))
if("colors" %in% colnames(geneInfo)) {
    geneInfo$module <- geneInfo$colors
} else {
    geneInfo$module <- geneInfo$all_12
}
modules <- geneInfo[,-grep("^all_|^treated_|^colors", colnames(geneInfo))]

## ## get rna-seq limma data @ 4hrs
## rnaseq <- get.rnaseq()
## ens2gene <- fread(file.path(CD4CHIC.DATA,"merged_samples_12Apr2015_full_denorm_bait2baits_e75.tab"),
##                   select=1:3) %>% unique()
## names(ens2gene) <- c("id","gene","biotype")
## setkey(ens2gene,id)
## setkey(rnaseq,id)
## rnaseq <- merge(rnaseq,ens2gene,all.x=TRUE,by="id")

## microarray data
d <- file.path(CD4CHIC.OUT,"modules","limma")
load(file.path(d,"microarray-vs-time0.RData")) # RESULTS
for(nm in names(RESULTS)) {
  RESULTS[[nm]]$time <- nm
}
O <- do.call("rbind",RESULTS)
O <- as.data.table(O)
write.table(O[,.(probeset.id,ens.gene.id,gene.symbol,time,logFC,adj.P.Val)],file=file.path(CD4CHIC.EXTDATA,"microarray-diffexpr.csv"),row.names=FALSE)

load(file=file.path(d,"microarray-vs-matched-unstim.RData")) # PRESULTS

################################################################################

## * Fig 1: modules - timecourse of expression

################################################################################

## NB WGCNA not installed on darwin
library(WGCNA)
{spower <- 12
 f <- out.file("bwnet")
 (load(f))
 bwnet <- bwnet$all
 s.exp <- s.exp$all
}
colors = labels2colors(bwnet$colors)

## gene counts per module
message("gene counts per module")
table(colors)

## dendrogram
plotter <- function(bwnet,...) {
    mergedColors = labels2colors(bwnet$colors)
    ## Plot the dendrogram and the module colors underneath
    plotDendroAndColors(bwnet$dendrograms[[1]], mergedColors[bwnet$blockGenes[[1]]],
                        "Module colors",
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05,...)
}
plotter(bwnet,main=paste(network.type,"all samples"))

nGenes = ncol(s.exp)
nSamples = nrow(s.exp)
pheno.cols <- c('uniqueID','pair','sex','time','treatment')

## moduleEigengenes
## averageExpr
##     If align == "along average", a dataframe containing average normalized expression in each module. The columns are named by the corresponding color with an "AE" prepended, e.g., AEturquoise etc.


getMEs <- function(exp,net) {
    cols <- labels2colors(net$colors)
    MEs <- moduleEigengenes(exp, cols, align="along average")$averageExpr #$averageExpr #$eigengenes
    pheno <- pheno.data[match(rownames(exp), pheno.data$sname), pheno.cols ]
    df <- cbind(pheno,MEs)
}
## Module eigengenes
MEs <- getMEs(s.exp,bwnet)
MEs$samples <- "all"

## module counts
cols <- labels2colors(bwnet$colors)
tt <- table(cols)
count.df <- data.frame(variable=names(tt),n=tt)

## Module patterns
    
df <- melt(MEs, c(pheno.cols,"samples"))
    df <- within(df, {
        id.treat <- paste0(df$uniqueID,df$treatment)
        treatment <- factor(df$treatment, levels=c("US","S"))
    })
    df0 <- df[,setdiff(colnames(df), c("time","value"))] %>% unique()
    df0 <- within(df0, {
        time <- 0
        value <- 0
    })
    df <- rbind(df,df0)
    library(dplyr)
    df$time <- as.factor(as.numeric(df$time))
    df <- df[order(df$treat,df$time),]
    df <- group_by(df, time, treatment, variable, samples)
dfs <- summarise(df, value=mean(value))


df$variable %<>% as.character() %>% sub("ME|AE","",.) %>%
factor(., levels=c("black","green","yellow","red","grey","blue","turquoise", "brown"  ))
dfs$variable %<>% as.character() %>% sub("ME|AE","",.) %>%
factor(., levels=c("black","green","yellow","red","grey","blue","turquoise", "brown"  ))
df$treatment <- relevel(df$treatment,"S")
levels(df$treatment) <- c("activated","non-act'd")
dfs$treatment <- relevel(dfs$treatment,"S")
levels(dfs$treatment) <- c("activated","non-act'd")

(load(file.path(CD4CHIC.DATA,"modules","hallmark-enrichment-results.RData")))
top3 <- lapply(top3, function(x) {
    x <- x[x>0]
    x <- x[1:min(length(x),3)]
    names(x) <- tolower(names(x)) %>% sub("hallmark_","",.) %>% gsub("_"," ",.) %>%
        sub("mtorc1","mTORC1",.) %>% sub("myc","Myc",.) 
    up <- c("il2","tnf","nfkb","stat","e2f","dna","tgf")
    for(uu in up)
        names(x) %<>% sub(uu,toupper(uu),.)
    names(x) %<>% sub("TNFa","TNF*alpha",.)
    names(x) %<>% sub("TGF beta ","TGF*beta~",.)
    names(x) %<>% sub("NFKB","NF*kappa*B",.)
    names(x) %<>% sub("interferon alpha ","interferon~alpha~",.)
    names(x) %<>% gsub(" ","~",.)
    x
})
top3
top3.df <- lapply(top3,function(x) {
    if(!any(is.na(x))){
        names(x)
#        paste(names(x),collapse="\n")
    } else {
        NULL
    }
})
top3.df <- top3.df[!sapply(top3.df,is.null)]
top3.df

top3.df <- data.frame(variable=rep(names(top3.df), sapply(top3.df,length)),
                      y=unlist(lapply(top3.df, seq_along)),
                      label=unlist(top3.df),
                      stringsAsFactors=FALSE)
#top3.df$variable <- factor(top3.df$variable,levels=levels(df$variable))
count.df$n.Freq <- paste0("n=",gsub("n=","",count.df$n.Freq))
count.df$variable <- factor(as.character(count.df$variable),levels=levels(df$variable))

df$variable %<>% refactor.mods()
dfs$variable %<>% refactor.mods()
count.df$variable %<>% refactor.mods()
top3.df$variable %<>% refactor.mods()

library(cowplot)
p <- ggplot(df,aes(x=time,y=value,col=variable)) +
geom_path(mapping=aes(group=id.treat,linetype=treatment),size=0.1) +
geom_hline(yintercept=0,col="grey") +
geom_path(data=dfs,aes(x=time,y=value,colour=variable,group=treatment,linetype=treatment),size=1) +
geom_text(mapping=aes(label=sub("~v1","",label),y=2 - (y-1)*0.25),x=1,data=top3.df,hjust=0,vjust=1,size=5,parse=TRUE) +
geom_text(mapping=aes(label=n.Freq),x=5,y=-2,data=count.df,hjust=1,vjust=1,size=5,fontface="bold") +
colScale +
scale_linetype_discrete("Treatment") +
labs(x="Time (hours)",y="Average change in expression from time 0 (log2 scale)") +
    facet_wrap(~variable,nrow=2) +
theme(legend.justification=c(0,0), legend.position=c(0,0),
      legend.background = element_rect(size=0.5,colour="black",linetype="solid"),
      text = element_text(size=16),
      legend.title = element_text(size=16))
                                        #    scale_colour_manual("Treated",values=c("grey20","darkblue")) +
#ggtitle(paste("Eigengenes in each individual and averaged (thick line) by module,",network.type,"power =",spower))
    p

ff <- file.path(CD4CHIC.OUT,"paper","figure-modules.pdf")
pdf(ff,width=12,height=8)
p
dev.off()
tiff(sub("pdf","tiff",ff),width=12,height=8,units="in",res=600)
p
dev.off()
png(sub("pdf","png",ff),width=12,height=8,units="in",res=600)
p
dev.off()

#system(paste("display ",ff))

######################################################################

## * Fig 2

######################################################################

expr$module %<>% refactor.mods()
m <- merge(ints,expr,by="id")

## ** Violin plots

geom_my <- function(...) { geom_violin(...) }
plotter <- function(d, y) {
    dd <- d[,.N,by=module]
    ggplot(d, aes_string(x="module",y=y,col="module")) +
    geom_my() +
    geom_hline(yintercept=0,col="grey",linetype="dashed") +
    geom_text(aes(x=module,label=paste0("n=",N)),y=-10,data=dd) +
    colScale +
    ylim(-12,12)
}

D <- list(expr[FDR.expr<0.01,.(module,logFC.expr)],
          m[pmax(Total_CD4_Activated, Total_CD4_NonActivated)>5 & FDR.expr<0.01,.(module,logFC.chic)],
          m[pmax(Total_CD4_Activated, Total_CD4_NonActivated)>5 & FDR.expr<0.01 & FDR.chic<0.1,
            .(module,logFC.chic)])
D <- lapply(D, function(d) setnames(d,c("module","y")))
D[[1]][,what:="Differential expression"]
D[[2]][,what:="All interactions"]
D[[3]][,what:="Differential interactions"]
D <- do.call("rbind",D[c(1,3)])
D[,what:=factor(what,levels=c("Differential expression",## "All interactions",
                              "Differential interactions"))]
dd <- D[,.N,by=c("what","module")]
p <- ggplot(D, aes_string(x="module",y="y",col="module")) +
geom_hline(yintercept=0,col="grey",linetype="dashed") +
geom_my(fill="#ffffff00") +
geom_text(aes(x=module,label=paste0("n=",N)),y=-12,data=dd) +
colScale +
facet_grid(what ~ .) +
ylim(-12,12) +
theme(strip.text=element_text(size=12)) +
labs(y="Log2 fold change",x="Gene module")
p + background_grid()

f <- file.path(CD4CHIC.OUT,"paper/violin-by-module2.jpg")
jpeg(f,height=5,width=6,units="in",res=600)
p + background_grid()
dev.off()

## system(paste("display",f))

## ** Dot plots
Dm <- D[,.(y=median(y)),by=c("module","what")]
Dm
tt <- with(Dm,cor.test(y.expr,y.int))
str(tt)
Dm <- merge(Dm[what=="Differential expression",],Dm[what=="Differential interactions",],by=c("module"),suffixes=c(".expr",".int"))
p2 <- ggplot(Dm,aes(x=y.expr,y=y.int)) + geom_point(aes(col=module),size=2) + geom_smooth(method="lm",se=FALSE,col="black") + colScale +
draw_text(sprintf("rho=%2.1f p=%4.3f",tt$estimate,tt$p.value),x=1,y=-0.7) +
labs(x="Median log2 fold change in expression",y="Median log2 fold change in interaction")


f <- file.path(CD4CHIC.OUT,"paper/violin-by-module-spearman2.jpg")
jpeg(f,height=5,width=6,units="in",res=600)
p2 + background_grid()
dev.off()
#system(paste("display",f))


## ** Beeplots
library(ggbeeswarm) # github.com/eclarke/ggbeeswarm
D <- list(expr[FDR.expr<0.01,.(module,logFC.expr)],
          m[pmax(Total_CD4_Activated, Total_CD4_NonActivated)>5 & FDR.expr<0.01,.(module,logFC.chic)],
          m[pmax(Total_CD4_Activated, Total_CD4_NonActivated)>5 & FDR.expr<0.01 & FDR.chic<0.1,
            .(module,logFC.chic)])
D <- lapply(D, function(d) setnames(d,c("module","y")))
D[[1]][,what:="Expression"]
D[[2]][,what:="All interactions"]
D[[3]][,what:="Differential interactions"]
D <- do.call("rbind",D)
D[,what:=factor(what,levels=c("Expression","All interactions","Differential interactions"))]
dd <- D[,.N,by=c("what","module")]
im <- unique(D[,im:=sum(y>0)/sum(y<0),by=c("what","module")],by=c("what","module"))
im[,im:=round(log(im),2)]

wh <- which(D$what=="All interactions")
D <- rbind(D[-wh,], D[sample(wh,10000),])

p <- ggplot(D, aes_string(x="module",y="y")) +
geom_hline(yintercept=0,col="grey",linetype="dashed") +
geom_quasirandom(size=0.5,alpha=0.1,col="grey") +
#geom_violin(fill="#ffffff00") +
geom_text(aes(x=module,label=paste0("n=",N),col=module),y=-12,data=dd) +
geom_pointrange(aes(x=module,y=im,ymin=0,ymax=im,col=module),data=im) +
colScale +
facet_grid(what ~ .) +
ylim(-13,12) +
labs(y="Log2 fold change")

f <- file.path(CD4CHIC.OUT,"paper/bees-by-module.png")
png(f,height=8,width=8,units="in",res=600)
p
dev.off()
## system(paste("display",f))



######################################################################

## *                   microarray vs rnaseq                          

######################################################################

## microarray data
d <- file.path(CD4CHIC.OUT,"modules","limma")
load(file.path(d,"microarray-vs-time0.RData")) # RESULTS
MA <- lapply(RESULTS,as.data.table)
for(nm in names(MA))
    MA[[nm]][,time:=as.numeric(nm)]
MA <- do.call("rbind",MA)
setnames(MA,"ens.gene.id","id")
MA <- MA[,.(id,gene.symbol,time,logFC,adj.P.Val)]
MA$method <- "ma"

## rnaseq
RNA <- copy(rnaseq)
RNA$time <- 4
RNA$method <- "RNAseq"

setnames(RNA,"gene","gene.symbol")
df <- rbind(MA,RNA[,.(id,gene.symbol,time,logFC,adj.P.Val,method)])

setnames(MA,c("logFC","adj.P.Val"),c("logFC.MA","FDR.MA"))
setnames(RNA,c("logFC","adj.P.Val"),c("logFC.RNAseq","FDR.RNAseq"))
setkey(MA,id,gene.symbol)
setkey(RNA,id,gene.symbol)
df4 <- merge(MA[time==4,.(id,gene.symbol,logFC.MA,FDR.MA)],
             RNA[biotype=="protein_coding" & !is.na(logFC.RNAseq),.(id,gene.symbol,logFC.RNAseq,FDR.RNAseq)],
             all=TRUE)
ggplot(df4,aes(x=logFC.MA,y=logFC.RNAseq)) + geom_point() + geom_abline(col="red") + geom_smooth(method="lm")
cor <- with(df4,cor.test(logFC.MA,logFC.RNAseq,use="pair"))

lm(logFC.RNAseq ~ logFC.MA, data=df4) %>% summary()
##              Estimate Std. Error t value Pr(>|t|)    
## (Intercept) -0.071611   0.006873  -10.42   <2e-16 ***
## logFC.MA     1.145505   0.008799  130.18   <2e-16 ***


with(df4,table(ma=is.na(logFC.MA),rnaseq=is.na(logFC.RNAseq)))
df4$logFC.MA[is.na(df4$logFC.MA)] <- 0
df4$logFC.RNAseq[is.na(df4$logFC.RNAseq)] <- 0

library(cowplot)
theme_set(theme_cowplot())
ggplot(df4,aes(x=logFC.MA,y=logFC.RNAseq)) +
geom_hline(yintercept=c(-4,0,4,8),colour="grey",size=0.1) +
geom_vline(xintercept=c(-4,0,4,8),colour="grey",size=0.1) +
geom_point(aes(col=is.na(FDR.MA)|is.na(FDR.RNAseq)),size=1) +
geom_abline(col="black",linetype="dashed") +
geom_smooth(formula = y ~ x , method="lm", data=subset(df4,!is.na(df4$FDR.MA) & !is.na(df4$FDR.RNAseq)),col="blue",se=FALSE) +
scale_colour_manual("Missing",values=c("black","gray40")) +
scale_y_continuous(breaks=c(-4,0,4,8)) +
scale_x_continuous(breaks=c(-4,0,4,8)) +
theme(legend.position="none") +
annotate("text",x=7.5,y=-4.5,label=paste0("rho == ",round(cor$estimate,2)),parse=TRUE,hjust=1,vjust=0) +
labs(x="Microarray: log2 fold change", y="RNAseq: log2 fold change")

ggsave(file.path(CD4CHIC.OUT,"paper","figure-microarray-vs-rnaseq.pdf"),height=6,width=8)

## ggplot(df4,aes(x=logFC.MA,y=-log10(FDR.MA),col=-log10(FDR.RNAseq))) + geom_point()
## ggplot(df4,aes(x=logFC.RNAseq,y=-log10(FDR.RNAseq),col=-log10(FDR.MA))) + geom_point()

######################################################################

## *                PIRs per bait & vice versa                      ##

######################################################################

f <- function(x) {
    table(x$biotype)
}
message("N int in act only\t");print(f(ints[Total_CD4_Activated>5 & Total_CD4_NonActivated<=5,]))
message("N int in non only\t",nrow(ints[Total_CD4_Activated<=5 & Total_CD4_NonActivated>5,]))
message("N int in act+non\t",nrow(ints[Total_CD4_Activated>5 & Total_CD4_NonActivated>5,]))

f <- function(x) {
    unique(x$baitID)
}
g.both <- f(ints[Total_CD4_Activated>5 & Total_CD4_NonActivated>5,])
g.act <- f(ints[Total_CD4_Activated>5 & Total_CD4_NonActivated<=5,]) %>% setdiff(., g.both)
g.non <- f(ints[Total_CD4_Activated<5 & Total_CD4_NonActivated>=5,]) %>% setdiff(., g.both)
message("Genes with int in act only\t", length(g.act))
message("Genes with in non only\t",length(g.non))
message("Genes with in act+non\t",length(g.both))

require(scales)
mylog_trans <- function(base=exp(1), from=0) 
{
  trans <- function(x) log(x, base)-from
  inv <- function(x) base^(x+from)
  trans_new("mylog", trans, inv, log_breaks(base=base), 
            domain = c(base^from, Inf))
}
plotter <- function(dt,x,ti="activated") {
    if(x=="baitsperprey") {
        xl <- "Baits per PIR"
        xm <- 80
        bw <- 1
        col <- "grey"
    } else {
        xl <- "PIRs per bait"
        col <- "black"
        xm <- 450
        bw <- 1
    }
    p <- ggplot(dt,aes_string(x=x)) +
    geom_histogram(col=col,binwidth=bw) +
    xlim(0,xm) +
    scale_y_continuous("Counts (log scale)",
                       trans = mylog_trans(base=10, from=-0.1),
                       breaks=10^seq(0,6),
                       expand=c(0,0)) +
    labs(x="Number") +
    ggtitle(paste0(xl, " (", tolower(substring(ti,1,1)), ")"))
}

tmp1 <- ints[Total_CD4_NonActivated>5,]
tmp1[,baitsperprey:=length(unique(baitID)),by="oeID"]
tmp1[,preysperbait:=length(unique(oeID)),by="baitID"]
p1 <- plotter(unique(tmp1,by="oeID"),"baitsperprey","non-activated")
p2 <- plotter(unique(tmp1,by="baitID"),"preysperbait","non-activated")
d1 <- ggplot(tmp1, aes(x=abs( (baitStart+baitLength/2) - (oeStart+oeLength/2) ))) + geom_histogram(col="grey") + scale_x_log10("Distance (log scale)",
                                                                                                                         breaks=10^c(4:8),labels=c("10kb","","1mb","","100mb")) +
ggtitle("Int. distance (n)")

tmp2 <- ints[Total_CD4_Activated>5,]
tmp2[,baitsperprey:=length(unique(baitID)),by="oeID"]
tmp2[,preysperbait:=length(unique(oeID)),by="baitID"]
p3 <- plotter(unique(tmp2,by="oeID"),"baitsperprey","Activated")
p4 <- plotter(unique(tmp2,by="baitID"),"preysperbait","Activated")
d2 <- ggplot(tmp2, aes(x=abs( (baitStart+baitLength/2) - (oeStart+oeLength/2) ))) + geom_histogram(col="grey") + scale_x_log10("Distance (log scale)",
                                                                                                                         breaks=10^c(4:8),labels=c("10kb","","1mb","","100mb")) +
ggtitle("Int. distance (a)")

plot_grid(p2,p1,d1,p4,p3,d2,ncol=3,labels=c("a","b","c","d","e","f"))

summary(tmp1[,abs( (baitStart+baitLength/2) - (oeStart+oeLength/2) )])
summary(tmp2[,abs( (baitStart+baitLength/2) - (oeStart+oeLength/2) )])

f <- file.path(CD4CHIC.OUT,"paper/bait-prey-counts.pdf")
pdf(f,height=6,width=8)
plot_grid(p2,p1,d1,p4,p3,d2,ncol=3,labels=c("a","b","c","d","e","f"))
dev.off()
system(paste("display ",f))

message("PAPER: each interacting baited fragment connected to a median of ...")
L <- list(n.bpp=unique(tmp1,by="oeID")$baitsperprey,
          n.ppb=unique(tmp1,by="baitID")$preysperbait,
          a.bpp=unique(tmp2,by="oeID")$baitsperprey,
          a.ppb=unique(tmp2,by="baitID")$preysperbait)
lapply(L,summary)

message("PAPER: median distance between interacting fragments of ",
with(ints[Total_CD4_NonActivated>5 | Total_CD4_Activated>5, ],
     median(abs( (baitStart+baitLength/2) - (oeStart+oeLength/2) )))/1000,
"kb")

######################################################################

## *                Interactions vs expression                        ##

######################################################################

m <- merge(ints,rnaseq[,.(id,logFC.expr,FDR.expr,activ.1,activ.2,unstm.1,unstm.2)],by="id")
m[,y.act:=log2(pmean(activ.1,activ.2)+1)]
m[,y.non:=log2(pmean(unstm.1,unstm.2)+1)]
m[,nfrag.act:=sum(Total_CD4_Activated>5),by=c("baitID","id")]
m[,nfrag.non:=sum(Total_CD4_NonActivated>5),by=c("baitID","id")]
x <- unique(m,by=c("baitID","id"))

## SUPPLEMENTARY FIGURE
plot.expr.act <- ggplot(x, aes(x=log10(nfrag.act+1),y=y.act)) + geom_point() + geom_smooth(method="lm") +
  scale_x_continuous("Number of PIR fragments (log scale)",breaks=log10(c(1,2,11,101)),labels=c(0,1,10,100)) +
  labs(y="log2 cpm") + ggtitle("Activated")
plot.expr.non <- ggplot(x, aes(x=log10(nfrag.non+1),y=y.non)) + geom_point() + geom_smooth(method="lm") +
  scale_x_continuous("Number of PIR fragments (log scale)",breaks=log10(c(1,2,11,101)),labels=c(0,1,10,100)) +
  labs(y="log2 cpm") + ggtitle("Non-Activated")
plot_grid(plot.expr.act+ panel_border(),plot.expr.non+ panel_border()) 
ggsave(file.path(CD4CHIC.OUT, "paper/figure-expr-vs-nfrag.jpg"), height=5,width=8,units="in")
m$class <- "gene"
m[gene %in% TF[TF$TFClass_human=="TFclass",]$gene_symbol, class:="TF"]
m[gene %in% unlist(cytokines), class:="cytokine"]
x <- unique(m,by=c("baitID","id"))

plot.expr.act <- ggplot(x, aes(x=log10(nfrag.act+1),y=y.act)) + geom_point() + geom_smooth(method="lm") +
  scale_x_continuous("Number of PIR fragments (log scale)",breaks=log10(c(1,2,11,101)),labels=c(0,1,10,100)) +
facet_wrap(~class) + labs(y="log2 cpm") + ggtitle("Activated")
plot.expr.non <- ggplot(x, aes(x=log10(nfrag.non+1),y=y.non)) + geom_point() + geom_smooth(method="lm") +
  scale_x_continuous("Number of PIR fragments (log scale)",breaks=log10(c(1,2,11,101)),labels=c(0,1,10,100)) +
facet_wrap(~class) +
  labs(y="log2 cpm") + ggtitle("Non-Activated")
plot_grid(plot.expr.act+ panel_border(),plot.expr.non+ panel_border()) 

