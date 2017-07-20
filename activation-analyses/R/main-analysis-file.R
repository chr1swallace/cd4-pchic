## * Startup and collect data

library(data.table)
library(magrittr)
#library(WGCNA)
# library(dplyr)
library(ggplot2); theme_set(theme_minimal())
#library(NetIndices)
library(cowplot) 
setwd(file.path(CD4CHIC.ROOT,"activation-analyses"))
source("R/common.R")
library(boot)
library(gridExtra)
source(file.path(CD4CHIC.ROOT,"activation-analyses/R/collect-data.R"))

## * Summary stats

ints[,class:=paste0(ifelse(Total_CD4_NonActivated>5,"N","-"),
                    ifelse(Total_CD4_Activated>5,"A","-"))]
ints[,class:=factor(class, levels=c("--","N-","-A","NA"))]
ints[,gene.class:=levels(class)[ max(as.numeric(class)) ],by="id"]
ints[,gene.class:=factor(gene.class, levels=levels(class))]
with(ints[FDR.chic<0.01],table(sign(logFC.chic)))

baitoe <- unique(ints, by=c("baitID","oeID"))
with(baitoe, table(class))[-1]

baitoe <- unique(ints[biotype=="protein_coding",], by=c("baitID","oeID"))
with(baitoe, table(class))[-1]

b2g <- get.b2gene()
message("n. baited frags")
length(unique(b2g$baitID))
message("n. baited with int in CD4")
length(unique(ints[class!="--",]$baitID))

with(ints,table(biotype,class))[,-1]
with(unique(ints,by="id"),table(biotype,gene.class))[,-1]

## * Restrict to protein coding from here on

ints.bak <- copy(ints)
ints <- ints[biotype=="protein_coding",]

## * basic stats: baits per prey and vice versa

## *** counts of numbers of interactions/genes expressed > 1 trans/mill in each state + overlap
##' interactions:
n1 <- (nrow(ints[Total_CD4_NonActivated>5 & Total_CD4_Activated<5,])) %>% floor()
n2 <- (nrow(ints[Total_CD4_NonActivated>5 & Total_CD4_Activated>5,])) %>% floor()
n3 <- (nrow(ints[Total_CD4_NonActivated<5 & Total_CD4_Activated>5,])) %>% floor()
c(n1,n2,n3)
sum(c(n1,n2,n3))
c(n1,n2,n3)/sum(c(n1,n2,n3))

##' genes expressed:
n1 <- (nrow(expr[expr.reads.non>1 & expr.reads.act<1,])) %>% floor()
n2 <- (nrow(expr[expr.reads.non>1 & expr.reads.act>1,])) %>% floor()
n3 <- (nrow(expr[expr.reads.non<1 & expr.reads.act>1,])) %>% floor()
c(n1,n2,n3)
sum(c(n1,n2,n3))
c(n1,n2,n3)/sum(c(n1,n2,n3))

##' interacting baits
int.baits <- ints[Total_CD4_NonActivated>5 | Total_CD4_Activated>5,]$baitID %>% unique()
non.baits <- ints[Total_CD4_NonActivated>5,]$baitID %>% unique()
act.baits <- ints[Total_CD4_Activated>5,]$baitID %>% unique()

message("PAPER: ",length(int.baits)," had at least one interaction")
setdiff(non.baits,act.baits) %>% length()
setdiff(act.baits,non.baits) %>% length()
intersect(non.baits,act.baits) %>% length()

##' interacting preys
int.preys <- ints[Total_CD4_NonActivated>5 | Total_CD4_Activated>5,]$oeID %>% unique()
non.preys <- ints[Total_CD4_NonActivated>5,]$oeID %>% unique()
act.preys <- ints[Total_CD4_Activated>5,]$oeID %>% unique()
length(int.preys)
setdiff(non.preys,act.preys) %>% length()
setdiff(act.preys,non.preys) %>% length()
intersect(non.preys,act.preys) %>% length()

setdiff(int.preys,int.baits) %>% length()
setdiff(non.preys,c(int.baits,act.preys)) %>% length()
setdiff(act.preys,c(int.baits,non.preys)) %>% length()
intersect(non.preys,act.preys) %>% setdiff(., int.baits) %>% length()


##' genes expressed, restricted to non prom:
tmp <- rnaseq[!is.na(FDR.expr),]
with(tmp,table(type,class))
with(tmp[FDR.expr<0.01,],table(type,class))


n1 <- (nrow(tmp[id %in% ints$id & expr.reads.non>1 & expr.reads.act<1,])) %>% floor()
n2 <- (nrow(expr[id %in% ints$id & expr.reads.non>1 & expr.reads.act>1,])) %>% floor()
n3 <- (nrow(expr[id %in% ints$id & expr.reads.non<1 & expr.reads.act>1,])) %>% floor()
c(n1,n2,n3)
sum(c(n1,n2,n3))
c(n1,n2,n3)/sum(c(n1,n2,n3))

summary(expr$expr.reads.non)
summary(expr$expr.reads.act)
ks.test(expr$expr.reads.non,expr$expr.reads.act)
ks.test(expr$expr.reads.non/sum(expr$expr.reads.non),expr$expr.reads.act/sum(expr$expr.reads.act))
qqplot(expr$expr.reads.non,expr$expr.reads.act)
qqplot(log(expr$expr.reads.non),log(expr$expr.reads.act))

xx <- unique(expr[pmax(expr.reads.non,expr.reads.act)>1,],by="id")
nrow(xx) # genes expressed
with(xx, table(FC=sign(logFC.expr),FDR=FDR.expr<0.01))

## * Restrict to non-promiscuous 
ints[,ngenesperbait:=length(unique(id)),by="baitID"]

## write.table(ints[pmax(Total_CD4_Activated,Total_CD4_NonActivated)>5,
##                  .(baitID,oeID,
##                    baitChr,baitStart,baitLength,
##                    oeChr,oeStart,oeLength,
##                    id,gene,Total_CD4_Activated, Total_CD4_NonActivated,
##                    b2b,ngenesperbait, P1.non,P2.non,P3.non,P1.act,P2.act,P3.act,logFC,FDR)],
##             file=file.path(CD4CHIC.OUT,"paper/FINAL/supp-table-interactions.csv"))

table(unique(ints,by="baitID")$ngenesperbait)
table(unique(ints,by="baitID")$ngenesperbait==1)
ints <- ints[ngenesperbait==1,]    


## * SOME FIGURES AND SUPPLEMENTAL FIGURES
## expression modules
## volcano/violin plots
## preys per bait etc
## microarray vs rnaseq
## number of prey vs expr

source(file.path(CD4CHIC.ROOT,"activation-analyses/R/figures.R"))

## erna vs prot genes
source(file.path(CD4CHIC.ROOT,"activation-analyses/R/eRNA.R"))


## * Count contiguous blocks, make intb (needed for following sections)

setkey(ints, baitID, oeID)
ints <- ints[order(ints$baitID,ints$oeID),]
#act[gene=="AHR",c("baitID","oeID",v,"oenext","oeprev","runid","runlength","peak"),with=FALSE]
addblocks <- function(x,v,l=NULL) {
    i <- which(x[[v]]>=5)
    act <- unique(x[i,],by=c("baitID","oeID"))
    act[,oenext:=c(0,diff(oeID)),by=baitID]
    act[,oeprev:=c(diff(oeID),0),by=baitID]
    act[,runid:=cumsum(as.numeric(oenext!=1)),by=baitID]
    act[,runlength:=sum(oenext==1)+1,by=c("baitID","runid")]
    head(act[gene=="AHR",])
    ##
    ## peak or not
    act[,score:=act[[v]]]
    act[,score.prev:=c(0,score[-length(score)]),by=c("baitID","runid")]
    act[,score.next:=c(score[-1],0),by=c("baitID","runid")]
    act[,peak:=score>score.prev & score>score.next, by=c("baitID","runid")]
    act[is.na(peak),peak:=FALSE]
    act[,npeak:=sum(peak),by=c("baitID")]
    act[,nrun:=length(unique(runid)),by=c("baitID")]
    act[,nfrag:=.N,by=c("baitID")]
    ##
    ## print some summaries
    tmp <- unique(act,by=c("baitID"))
    for(v in c("nrun","nfrag")) {
        message(v," ",l)
        quantile(tmp[[v]]) %>% print()
    }
    tmp <- unique(act,by=c("baitID","runid"))
    for(v in c("runlength","npeak")) {
        message(v," ",l)
        quantile(tmp[[v]]) %>% print()
    }
    ##
    ## return
    nm <- c("runid","runlength","peak","npeak","nrun","nfrag")
    if(!is.null(l)) {
        setnames(act, nm, paste(nm,l,sep="."))
        nm <- paste(nm,l,sep=".")
    }
    merge(x,act[,c("baitID","oeID",nm),with=FALSE],by=c("baitID","oeID"),all.x=TRUE)
}

ints[,CD4:=pmax(Total_CD4_Activated,Total_CD4_NonActivated)]
intb <- addblocks(ints,"Total_CD4_Activated","act")
intb <- addblocks(intb,"Total_CD4_NonActivated","non")
intb <- addblocks(intb,"CD4","cd4")
head(intb[gene=="AHR",])  ## NB to label a run, need bait AND run.act or run.non


## * Chromatin peaks at baits/oe

# self contained - run externally
#source(file.path(CD4CHIC.ROOT,"activation-analyses/R/chromatin-v2.R"))


## * Gene expression and mRNA half life (can skip)

## ** Expression differences vs HL 

##'
##' LOG FC is significantly less for genes with shorter half lives in
##' activated cells and slightly higher for genes with shorter half lives in
##' non-activated
lm(logFC.expr ~ log(HL.328) + log(HL.M), data=expr) %>% summary()
lm(log(HL.328) ~ module, data=expr) %>% summary()
lm(log(HL.M) ~ module, data=expr) %>% summary()

colScale <- make.colScale(unique(expr$module))

library(cowplot)
tmp <- expr[complete.cases(HL.328,HL.M,module),]
## library(preprocessCore)
## tmpn <- normalize.quantiles(as.matrix(tmp[,.(HL.328,HL.M)]))
## tmp[,NHL.328:=tmpn[,1]]
## tmp[,NHL.M:=tmpn[,2]]
tmp2 <- melt(tmp[,.(HL.328,HL.M,module,gene)],c("gene","module"))
    tmp2[,variable:=sub("HL.M","Non-activated",variable)]
    tmp2[,variable:=sub("HL.328","Activated",variable)]
tmp2[,variable:=factor(variable,levels=c("Non-activated","Activated"))]

## p1 <- ggplot(tmp, aes(y=logFC.expr,x=log(HL.M))) + geom_point() + geom_smooth(method="lm") + labs(x="Half life (resting cells)",y="log2 fold change expression")
## p2 <-
##     ggplot(tmp, aes(y=logFC.expr,x=log(HL.328),col=module)) + geom_point() + geom_smooth(method="lm") + labs(x="Half life (activated cells)",y="log2 fold change expression") + facet_wrap(~module) + colScale
## p4 <- ggplot(tmp,
##        aes(x=module,y=log(HL.328),col=module)) + geom_boxplot() + geom_point(alpha=0.3) + colScale + labs(x="Module",y="Half life (activated cells)")

ggplot(tmp2,aes(x=module,y=value/60,col=module)) + geom_boxplot() + geom_point(alpha=0.3) + #geom_path(alpha=0.2,aes(group=gene)) +
facet_grid(.~variable) + colScale +
  scale_y_log10("half life, hours (log scale)",breaks=c(1,2,5,10,20,50,100)) +
  labs(x="Timecourse gene module") +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
ggsave(file=file.path(CD4CHIC.OUT,"paper/halflife-by-module.pdf"),height=6,width=8)
ggsave(file=file.path(CD4CHIC.OUT,"paper/halflife-by-module.tiff"),height=6,width=8)
write.table(tmp[,.(gene,module,HL.M,HL.328)],file=file.path(CD4CHIC.OUT,"paper/halflife-by-module.tsv"))

models <- list(act=lm(HL.328 ~ module, data=tmp),
               non=lm(HL.M ~ module , data=tmp))
lapply(models, summary)
lapply(models,anova)
##' HL in the two states are correlated, and on average lower in activated (latter probably bias)
with(expr,cor.test(log(HL.M),log(HL.328),use="pair"))
ggplot(subset(expr,!is.na(logFC.expr)),
       aes(x=log(HL.M),y=log(HL.328))) + geom_point() + geom_smooth(method="lm") + geom_abline(col="blue") + facet_wrap(~module)

## * Number of interacting fragments vs expression

## ** Absolute expression levels
chromatin <- get.chromatin()
m <- merge(intb,expr[,.(id,logFC.expr,FDR.expr,HL.328,HL.M,activ.1,activ.2,unstm.1,unstm.2)],by="id")
m <- merge(m,chromatin,by.x="oeID",by.y="name")
m <- m[!is.na(FDR.expr) & !b2b,]
m[,y.act:=log2(activ.1+activ.2+1)]
m[,y.non:=log2(unstm.1+unstm.2+1)]
m[,nfrag.act:=sum(Total_CD4_Activated>5),by=c("baitID","id")]
m[,nfrag.non:=sum(Total_CD4_NonActivated>5),by=c("baitID","id")]

## adjust for transcript length
tl <- get.translength()
tmp <- unique(m,by="id")
tmp <- merge(tmp,tl,by="id")
tmp <- tmp[is.finite(y.act) & is.finite(y.non) & y.non>0,]
tmp$res.act <- loess(log(y.act) ~ max.trans.length, data=tmp) %>% residuals() %>% exp()
tmp$res.non <- loess(log(y.non) ~ max.trans.length, data=tmp) %>% residuals() %>% exp()
m <- merge(m,tmp[,.(id,res.act,res.non)],by="id")

## adjust for half life
tmp <- unique(m,by="id")
tmp <- tmp[complete.cases(HL.328,HL.M),]
tmp$hlres.act <- loess(log(res.act) ~ HL.328, data=tmp) %>% residuals() %>% exp()
tmp$hlres.non <- loess(log(res.non) ~ HL.M, data=tmp) %>% residuals() %>% exp()
m <- merge(m,tmp[,.(id,hlres.act,hlres.non)],by="id")

x <- unique(m,by=c("baitID","id"))

## models really v similar with/without half life adjustment, just sample size vs signal/noise ratio
glm(res.act ~ nfrag.act, data=x) %>% summary()
glm(hlres.act ~ nfrag.act, data=x) %>% summary()

glm(res.non ~ nfrag.non, data=x) %>% summary()
glm(hlres.non ~ nfrag.non, data=x) %>% summary()

## SUPPLEMENTARY FIGURE - see figures.R

mg <- glm(res.non ~ log(nfrag.non+1), family=Gamma(link="log"), data=x)
ml <- glm(res.non ~ log(nfrag.non+1), data=x)

## par(mfrow=c(2,2))
## plot(mg)
## plot(ml)

##' use linear model!  Hip hip hooray!!!

##' NB, half life lower, in activated, for genes with large number of interactions
glm(log(HL.M) ~ log(nfrag.non+1), data=x) %>% summary()
glm(log(HL.328) ~ log(nfrag.act+1), data=x) %>% summary()
glm(log(HL.328)-log(HL.M) ~ log(nfrag.act+1) + log(nfrag.non+1), data=x) %>% summary()
glm(log(HL.328) ~ logFC.expr, data=x) %>% summary()
glm(log(HL.M) ~ logFC.expr, data=x) %>% summary()
glm(logFC.expr ~ log(HL.328) + log(HL.M),data=x) %>% summary()
glm(logFC.expr ~ log(HL.328) + log(HL.M) + log(nfrag.non+1) + log(nfrag.act+1),data=x) %>% summary()


## ** Change in # interacting fragments vs change in expression

plot.expr.delta <- ggplot(x, aes(x=log10(nfrag.act+1) - log10(nfrag.non+1),y=logFC.expr)) + geom_point() + geom_smooth() +
  scale_x_continuous("Ratio of  prey fragments in activated vs non-activated cells (log scale)",breaks=log10(c(0.1,0.5,1,2,10)),labels=c("1/10","1/2",1,2,10)) +
  labs(y="Log2 fold change")
#plot.expr.delta

library(rms)
m$HL.M <- m$HL.M/60 # convert to hours
m$HL.328 <- m$HL.328/60 # convert to hours
m$up <- ifelse(m$logFC.chic>0,1,0)
m <- m[pmax(Total_CD4_Activated,Total_CD4_NonActivated)>=5,]

tmp <- unique(m,by="id")
tmp <- tmp[complete.cases(HL.M,HL.328),]
tmp$logFC.expr.adj <- residuals(lm(logFC.expr ~ HL.M + HL.328, data=tmp))
m <- merge(m,tmp[,.(id,logFC.expr.adj)],by="id")

m[,c("sig.up","sig.down","chic.up","chic.down"):=
     list(logFC.chic>0 & FDR.chic<0.1 & Total_CD4_Activated>=5 & Total_CD4_NonActivated<5,
          logFC.chic<0 & FDR.chic<0.1 & Total_CD4_Activated<5 & Total_CD4_NonActivated>=5,
          ifelse(logFC.chic>0 & FDR.chic<0.5 & Total_CD4_Activated>=5 & Total_CD4_NonActivated<5,1,0),
          ifelse(logFC.chic<0 & FDR.chic<0.5 & Total_CD4_Activated<5 & Total_CD4_NonActivated>=5,1,0))]

 with(m,table(Total_CD4_Activated>5,Total_CD4_NonActivated>5,sig.up,sig.down))
 with(m,table(chic.up,chic.down,sig.up,sig.down))
 with(m,table(Total_CD4_Activated>5,Total_CD4_NonActivated>5,chic.up,chic.down))

source(file.path(CD4CHIC.ROOT, "activation-analyses/R/robust-models.R"))

## *** witHout HL
##' basic effect of differential chic
modelcmp(list("logFC.expr ~ chic.up",
              "logFC.expr ~ chic.down",
              "logFC.expr ~ chic.up + chic.down",
              "logFC.expr ~ sig.up",
              "logFC.expr ~ sig.down",
              "logFC.expr ~ sig.up + sig.down",
              "logFC.expr ~ chic.up + sig.up + chic.down + sig.down"),
         subset(m,complete.cases(sig.up,sig.down)))
##' final model
mod.chicdiff <- model("logFC.expr ~ sig.up + sig.down", m)
model("logFC.expr ~ sig.up + sig.down", m)

##' compare to a model which assumes up and down are equivalent
m[,sig.num:=sig.up - sig.down]
m[,sig.dir:=sig.down>0]
modelcmp(list("logFC.expr ~ sig.up + sig.down",
              "logFC.expr ~ sig.num",
              "logFC.expr ~ sig.num + sig.up"),
         subset(m,complete.cases(sig.up,sig.down)))
##              Coef   Lower.95  Upper.95            P
## sig.num 0.3303420 0.05782938 0.6028547 0.0175048883
## sig.up  0.9738593 0.46887399 1.4788445 0.0001569284

## *** with HL - very little difference
modelcmp(list("logFC.expr.adj ~ chic.up",
              "logFC.expr.adj ~ chic.down",
              "logFC.expr.adj ~ chic.up + chic.down",
              "logFC.expr.adj ~ sig.up",
              "logFC.expr.adj ~ sig.down",
              "logFC.expr.adj ~ sig.up + sig.down",
              "logFC.expr.adj ~ chic.up + sig.up + chic.down + sig.down"),
         subset(m,complete.cases(sig.up,sig.down)))
##' final model
model("logFC.expr.adj ~ sig.up + sig.down", m)

modelcmp(list("logFC.expr ~ sig.up + sig.down",
              "logFC.expr ~ sig.up + sig.down + log(HL.M)",
              "logFC.expr ~ sig.up + sig.down + log(HL.328)",
              "logFC.expr ~ sig.up + sig.down + log(HL.328) + log(HL.M)"),
         subset(m,complete.cases(sig.up,sig.down,HL.M,HL.328)))
##' final model
mod.chicdiff.HL <- model("logFC.expr ~ sig.up + sig.down + log(HL.328) + log(HL.M)", m)
mod.chicdiff.HL

## *** write chosen model
library(xtable)
writer <- function(o, caption="", row.names=character(nrow(o)),...) {
    print(o)
    o <- as.data.frame(o)
    if(length(wh <- which(colnames(o)=="P")))
        o[,wh] <- mysan(o[,wh])
    ot <- xtable(o,caption=caption)
rownames(ot) <- row.names
print(ot,booktabs=TRUE,sanitize.text.function=mysan,caption.placement="top",...)
}

writer(mod.chicdiff,
       caption="Effect upon change in gene expression (log$_2$ fold change) of gaining or losing a PIR upon activation of CD4$^+$ cells",
       row.names=c("Gain","Loss"),
       file=file.path(CD4CHIC.OUT,"paper/table-gainloss-vs-expr.tex"),
       append=FALSE # FIRST TABLE
       )
writer(mod.chicdiff.HL,
       caption="Effect upon change in gene expression (log$_2$ fold change) of gaining or losing a PIR upon activation of CD4$^+$ cells, including information on mRNA half life",
       row.names=c("Gain","Loss","Half life (activated)","Half life (non-activated)"),
       file=file.path(CD4CHIC.OUT,"paper/table-gainloss-vs-expr.tex"),
       append=TRUE
       )

## *** FIGURE plot chosen model

## functions
plotdata <- NULL
plotdata <- function(sf,sc) {
    s2 <- data.frame(Estimate=coefficients(sf),
                     se=vcov(sf) %>% diag() %>% sqrt())
    s2$Estimate[-1] <- s2$Estimate[-1] + s2$Estimate[1]
    s2$x <- rownames(s2) %>% as.character()
    s2$group <- ifelse(grepl("down",s2$x),"loss",
                ifelse(grepl("ac27",s2$x),"gain.ac27","gain"))
    s2$x %<>% sub("as.factor(up.ac27)","",.,fixed=TRUE)
    s2$x %<>% sub("as.factor(up)","",.,fixed=TRUE)
    s2$x %<>% sub("as.factor(down)","-",.,fixed=TRUE)
    s2$x %<>% gsub("fn.*up=","",.)
    s2$x %<>% gsub("fn.*down=","",.)
    s2$x[1] <- 0
    s2$x <- as.numeric(s2$x)
    s2 <- s2[c(1,1:nrow(s2)),]
    s2$group[1:2] <- c("gain","loss")
    df <- with(s2, data.frame(nsig.up=ifelse(group=="gain",x,0),
                              nsig.down=ifelse(group=="gain",0,x)))
    colnames(df) <- names(coef(sc))[-1]
    p <- predict(sc, newdata=df,se=TRUE)
    s2$Pred <- p$linear.predictors
    s2$lci <- p$linear.predictors - 1.96 * p$se.fit
    s2$uci <- p$linear.predictors + 1.96 * p$se.fit
    s2$x <- ifelse(s2$group=="gain", s2$x, -1 * s2$x)
    return(s2)
}
plotter <- function(s2,N=5) {
s2 <- s2[abs(s2$x)<=N,]
    ggplot(s2,aes(x=x,y=Estimate,ymax=Estimate+1.96*se,ymin=Estimate-1.96*se,group=group)) +
geom_smooth(method="lm",aes(y=Pred,ymin=lci,ymax=uci),data=subset(s2,group=="gain"),stat="identity",size=0.7) +
geom_smooth(method="lm",aes(y=Pred,ymin=lci,ymax=uci),data=subset(s2,group=="loss"),stat="identity",size=0.7) +
  geom_hline(yintercept=0,linetype="dashed",size=0.5) +
geom_pointrange(size=0.8) +
scale_x_continuous(breaks=c(-N:N),limits=c(-N-0.5,N+0.5)) +
labs(x="<- PIR lost                 PIR gained ->",
     y="log2 fold change in gene expression") +
background_grid(size.major=0.3,colour.major="grey80") }

N <- 5
mm <- m[,.(nsig.up=sum(sig.up,na.rm=TRUE),
           nsig.down=sum(sig.down,na.rm=TRUE)),
        by=c("baitID","id","logFC.expr","HL.328","HL.M")]
mm <- mm[nsig.up<=N & nsig.down<=N,]

mm[,nsig.change := nsig.up - nsig.down]
model(logFC.expr ~ nsig.change, mm, TRUE)
model(logFC.expr ~ nsig.change + nsig.up, mm, TRUE)
model(logFC.expr ~ nsig.up + nsig.down, mm, TRUE)
model(logFC.expr ~ nsig.up + nsig.down + HL.328 + HL.M, mm, TRUE)

mm[,c("fnsig.up","fnsig.down"):=
    list(factor(nsig.up), factor(nsig.down))]
sf <- model(logFC.expr ~ fnsig.up + fnsig.down -1, mm, FALSE)
sc <- model(logFC.expr ~ nsig.up + nsig.down, mm, FALSE)
data <- plotdata(sf,sc)

plotter(data) + background_grid()
ggsave(file.path(CD4CHIC.OUT,"paper/figure-change-expr-number-prey.pdf"),height=4,width=6)
ggsave(file.path(CD4CHIC.OUT,"paper/figure-change-expr-number-prey.tiff"),height=4,width=6)
write.table(data,file=file.path(CD4CHIC.OUT,"paper/figure-change-expr-number-prey.csv"))


## ** add erna
erna <- get.erna()
erna <- erna[type %in% c( "regulatory",
                         "intergenic.reg") & !is.na(P.Value),]
gr1 <- with(m, GRanges(Rle(oeChr), IRanges(oeStart, width=oeLength)))
addgr <- function(gr2) {
    tmp <- rep(FALSE,length(gr1))
    o <- findOverlaps(gr1,gr2)
    tmp[unique(o@from)] <- TRUE
    return(tmp)
}

## in any erna
gr2 <- with(erna, GRanges(Rle(paste0("chr",chr)), IRanges(start,end))); length(gr2)
m[,in.erna:=addgr(gr2)]
## in regulatory erna
gr2 <- with(erna, GRanges(Rle(paste0("chr",chr)), IRanges(start,end))); length(gr2)
m[,in.erna.reg:=addgr(gr2)]
## in upregulated erna
gr2 <- with(erna[FDR.erna<0.01 & logFC.erna>0,],
            GRanges(Rle(paste0("chr",chr)), IRanges(start,end))); length(gr2)
m[,in.erna.up:=addgr(gr2)]
## in downregulated erna
gr2 <- with(erna[FDR.erna<0.01 & logFC.erna<0,],
            GRanges(Rle(paste0("chr",chr)), IRanges(start,end))); length(gr2)
m[,in.erna.down:=addgr(gr2)]


## categorise each
m[,int.cat:=ifelse(pmax(Total_CD4_Activated,Total_CD4_NonActivated)<5, "back",
            ifelse(sig.up,"gain",
            ifelse(sig.down,"loss","invar")))]
m[,int.cat:=factor(int.cat, levels=c("invar","gain","loss"))]

m[,enh.cat:=ifelse(!in.erna,"no",
            ifelse(in.erna.up,"up",
            ifelse(in.erna.down,"down","invar")))]
m[,enh.cat:=factor(enh.cat, levels=c("no","invar","up","down"))]
with(m,table(enh.cat))

m[,bi.cat:=factor(paste0("int.",int.cat,"-enh.",enh.cat))]
m[,bi.cat:=relevel(bi.cat, "int.invar-enh.no")]

modelcmp(list("logFC.expr ~ int.cat",
              "logFC.expr ~ enh.cat",
              "logFC.expr ~ int.cat + enh.cat",
              "logFC.expr ~ bi.cat",
              "logFC.expr ~ int.cat * enh.cat"), m)

mod <- model("logFC.expr ~ bi.cat",m) %>%as.data.frame()
mod.ols <- ols(as.formula("logFC.expr ~ bi.cat"),x=TRUE,y=TRUE,m)
## NB bug in robcov, bi-cat with single observation becomes strongly significant
se <-mod.ols$var["bi.cat=int.gain-enh.down","bi.cat=int.gain-enh.down"]
mod["bi.cat=int.gain-enh.down","P"] <- 1

mod$x <- sub("bi.cat=","",rownames(mod))
mod <- rbind(mod[1,],mod)
mod$x[1] <- levels(m$bi.cat)[1]
mod[1,1:3] <- 0
## mod <- rbind(mod[1,],mod)
## mod$x[1] <- "int.up-enh.down"
## mod[1,1:3] <- NA
mod$n <- table(m$bi.cat)[mod$x] %>% as.vector()
mod$x <- factor(mod$x,levels=c("int.invar-enh.no", "int.invar-enh.invar",
                               "int.invar-enh.down", "int.invar-enh.up",
                               "int.loss-enh.no", "int.loss-enh.invar",
                               "int.loss-enh.down", "int.loss-enh.up",
                               "int.gain-enh.no", "int.gain-enh.invar",
                               "int.gain-enh.down","int.gain-enh.up"))
mod$enhancer <- sub(".*-enh.","",mod$x)
mod$interaction <- gsub("int.|-enh.*","",mod$x)
mod$y.enhancer <- 2
mod$y.interaction <- 1

#library(ggbio)

mod <- as.data.table(mod)
mod[,int2:=ifelse(n==1,"one",interaction)]

mod[,1:3] <- mod[,1:3] + 6
plots <- list(
p1=
ggplot(mod,aes(x=x,y=Coef,ymin=Lower.95,ymax=Upper.95)) + 
geom_hline(yintercept=6,col="grey") +
geom_pointrange(aes(col=int2)) +
geom_bar(aes(y=y.enhancer,fill=enhancer),col="grey20",stat="identity",data=mod) +
geom_bar(aes(y=y.interaction,fill=interaction),col="grey20",stat="identity",data=mod) +
geom_text(aes(y=y.enhancer-0.5,label=enhancer),size=3) +
geom_text(aes(y=y.interaction-0.5,label=interaction),size=3) +
scale_y_continuous("Effect: log fold change\nat protein-coding gene",breaks=c(-5.5,-4.5,-2,0,2)+6,labels=c("PIR","eRNA",-2,0,2)) +
scale_colour_manual(values=c("no"="white","invar"="grey40","down"="magenta","up"="DeepSkyBlue","loss"="magenta","gain"="DeepSkyBlue","one"="grey80")) +
scale_fill_manual(values=c("no"="white","invar"="grey","down"="magenta","up"="cyan","loss"="magenta","gain"="cyan"))
,
p2=
ggplot(mod,aes(x=x,y=n+1)) + geom_bar(stat="identity",col="grey20",fill="white") +
geom_text(aes(x=x,y=n+1,label=n),hjust=1.1,size=4,angle=90) +
scale_y_log10("Interaction\nCount",breaks=c(1,10,100,1000,10000)+1,labels=c(1,10,100,1000,10000))
)
plots <- lapply(plots, function(p) 
    p + background_grid() +
    theme(legend.position="none",
              axis.text.x=element_blank(),
#              axis.title.x=element_blank()
              axis.ticks.x=element_blank()
              )) #text(angle=45,hjust=1)))
mod[,1:3] <- mod[,1:3] - 6

library(cowplot)
plot_grid(plotlist=rev(plots),nrow=2,rel_heights=c(1,2))

#tracks(plots$p2,plots$p1,heights=c(2,4))


f <- file.path(CD4CHIC.OUT,"paper/figure-change-expr-prey-erna.tiff")
ggsave(file.path(CD4CHIC.OUT,"paper/figure-change-expr-prey-enra-screen.pdf"),height=8,width=10)
ggsave(file.path(CD4CHIC.OUT,"paper/figure-change-expr-prey-enra.pdf"),height=4.5,width=6)
ggsave(f,height=4.5,width=6)
#write.table(s2,file=file.path(CD4CHIC.OUT,"paper/figure-change-expr-number-prey.csv"))

system(paste("display",f))

## alternative plot - for talks
head(mod)
mod[,x:=as.numeric(factor(interaction,levels=c("loss","invar","gain")))]
mod[,y:=as.numeric(factor(enhancer,levels=c("down","no","invar","up")))]
mod[,sumn:=paste("n =",sum(n)),by=c("y")]
mod$xn <- 5
ggplot(mod, aes(xmin=x,
                xmax=x+1,
                ymin=y,
                ymax=y+1,
                fill=Coef-6)) +
geom_rect(col="grey") +
geom_text(aes(x=x+0.5,y=y+0.5,label=ifelse(Coef>6,"⬆","⬇")),
          data=subset(mod,Coef!=6 & P<0.02),size=5) +
geom_text(aes(x=xn,y=y+0.5,label=as.character(sumn)),adj=1) +
scale_fill_gradient2("Fold change in expression",low="red",high="blue") +
scale_x_continuous("Change in PIRs",breaks=1:3+0.5,labels=c("loss","no change","gain"),limits=c(1,5)) +
scale_y_continuous("Enhancer activity",breaks=1:4+0.5,labels=c("down","not detected","no change","up")) +
theme(legend.position="bottom")
ggsave("~/prot-erna-alt.jpg",height=3,width=5)
system("display ~/prot-erna-alt.jpg")

## alternative plots - for paper
head(mod)
mod[,x:=as.numeric(factor(interaction,levels=c("loss","invar","gain")))]
mod[,y:=as.numeric(factor(enhancer,levels=c("down","no","invar","up")))]
mod[,sumn:=paste("n =",sum(n)),by=c("y")]
mod[,stars:=""]
mod[P<0.01,stars:="*"]
mod[P<0.001,stars:="**"]
mod[P<0.0001,stars:="***"]
mod[,label:=sprintf("%6.2f %s\n[%6.2f,%6.2f ]",(Coef),stars,(Lower.95),(Upper.95))]
mod[enhancer=="no" & interaction=="invar",label:="1.00"]
mod[enhancer=="down" & interaction=="gain",label:=sprintf("%6.2f",(Coef))]

mod$xn <- 5
ggplot(mod, aes(xmin=x,
                xmax=x+1,
                ymin=y,
                ymax=y+1,
                fill=Coef)) +
geom_rect(col="grey") +
#geom_text(aes(x=x+0.5,y=y+0.5,label=ifelse(Coef>6,"⬆","⬇")),
                                        #          data=subset(mod,Coef!=6 & P<0.02),size=5) +
geom_text(aes(x=x+0.5,y=y+0.5,label=label)) + 
#geom_text(aes(x=xn,y=y+0.5,label=as.character(sumn)),adj=1) +
scale_fill_gradient2("log2 fold\nchange\nin gene\nexpression",low="darkorchid3",high="darkolivegreen3") +
scale_x_continuous("Change in PIRs",breaks=1:3+0.5,labels=c("loss","no change","gain")) +
scale_y_continuous("Enhancer activity",breaks=1:4+0.5,labels=c("down","not detected","no change","up")) +
theme(legend.position="right")


ggsave(file.path(CD4CHIC.OUT,"paper/figure-prot-erna-alt.pdf"),height=4*1.2,width=6*1.2)
#system(paste("display", file.path(CD4CHIC.OUT,"paper/figure-prot-erna-alt.pdf"))

## IRF4
m[gene=="IRF4" & enh.cat=="up" & int.cat=="gain",]




## *** some plots of a few things

data.plotter <- function(v) {
    message(v)
N <- 5
mm <- copy(m)
  mm[,test.sig.up:= sig.up & tmp[[v]]>0]
  mm[,test.sig.down:= sig.down & tmp[[v]]>0]
mm[,c("nsig.up","nsig.down","vnsig.up","vnsig.down"):=
    list(sum(sig.up & !test.sig.up,na.rm=TRUE),
         sum(sig.down & !test.sig.down,na.rm=TRUE),
         sum(test.sig.up,na.rm=TRUE),
         sum(test.sig.down,na.rm=TRUE)),
   by=c("baitID","id")]
mm <- unique(mm[nsig.up<=N & nsig.down<=N,],by=c("baitID","id"))
mm[,c("fnsig.up","fnsig.down","fvnsig.up","fvnsig.down"):=
    list(factor(nsig.up),
         factor(nsig.down),
         factor(vnsig.up),
         factor(vnsig.down))]
sf <- model(logFC.expr ~ fnsig.up + fnsig.down + fvnsig.up + fvnsig.down -1, mm, FALSE)
sc <- model(logFC.expr ~ nsig.up + nsig.down + vnsig.up + vnsig.down, mm, FALSE)
#print( sc )

s2 <- data.frame(Estimate=coefficients(sf),
                 se=vcov(sf) %>% diag() %>% sqrt())
s2$Estimate[-1] <- s2$Estimate[-1] + s2$Estimate[1]
s2$x <- rownames(s2) %>% as.character()
s2$group <- ifelse(grepl("down",s2$x),"loss","gain")
s2$group <- with(s2, ifelse(grepl("v",x), paste(group, "(+)"), paste(group, "(-)")))
s2$x %<>% sub("as.factor(up.ac27)","",.,fixed=TRUE)
s2$x %<>% sub("as.factor(up)","",.,fixed=TRUE)
s2$x %<>% sub("as.factor(down)","-",.,fixed=TRUE)
s2$x %<>% gsub("fnsig.up=","",.,fixed=TRUE)
s2$x %<>% gsub("fnsig.down=","",.,fixed=TRUE)
s2$x %<>% gsub("fvnsig.up=","",.,fixed=TRUE)
s2$x %<>% gsub("fvnsig.down=","",.,fixed=TRUE)
s2$x[1] <- 0
s2$x <- as.numeric(s2$x)
s2 <- s2[c(1,1,1,1:nrow(s2)),]
s2$group[1:4] <- c("gain (-)","gain (+)","loss (-)","loss (+)")
df <- with(s2, data.frame(nsig.up=ifelse(group=="gain (-)",x,0),
                          nsig.down=ifelse(group=="gain (-)",0,x),
                          vnsig.up=ifelse(group=="gain (+)",x,0),
                          vnsig.down=ifelse(group=="gain (+)",0,x)
                          ))
p <- predict(sc, newdata=df,se=TRUE)
s2$Pred <- p$linear.predictors
s2$lci <- p$linear.predictors - 1.96 * p$se.fit
s2$uci <- p$linear.predictors + 1.96 * p$se.fit
s2$subcat=grepl("+",s2$group,fixed=TRUE)
s2$x <- with(s2,
             ifelse(grepl("gain",s2$group),
                    s2$x, # + ifelse(subcat, 0.05,-0.05),
                    -1 * s2$x )) # - ifelse(subcat, 0.1,-0.1)))
s2$mark <- v
    return(s2)
}

 df <- lapply(marks.test,data.plotter) %>% do.call("rbind",.)
   
ggplot(df,aes(x=x,y=Estimate,ymax=Estimate+1.96*se,ymin=Estimate-1.96*se,group=group,col=subcat,fill=subcat)) +
geom_smooth(method="lm",aes(y=Pred,ymin=lci,ymax=uci),data=subset(s2,group=="gain (+)"),stat="identity",size=0.7) +
geom_smooth(method="lm",aes(y=Pred,ymin=lci,ymax=uci),data=subset(s2,group=="loss (+)"),stat="identity",size=0.7) +
geom_smooth(method="lm",aes(y=Pred,ymin=lci,ymax=uci),data=subset(s2,group=="gain (-)"),stat="identity",size=0.7) +
geom_smooth(method="lm",aes(y=Pred,ymin=lci,ymax=uci),data=subset(s2,group=="loss (-)"),stat="identity",size=0.7) +
  geom_hline(yintercept=0,linetype="dashed",size=0.5) +
geom_pointrange(size=0.8) +
scale_x_continuous(breaks=c(-5:5),limits=c(-5.5,5.5)) +
facet_wrap(~mark) +
labs(x="Number of prey lost (negative) or gained (positive)",
     y="log2 fold change in gene expression") +
background_grid(size.major=0.3,colour.major="grey80") + ggtitle(v) +
ylim(-5,5) 

ggsave(file.path(CD4CHIC.OUT,"paper/figure-change-expr-number-prey-chromatin.pdf"),height=8,width=8)
ggsave(file.path(CD4CHIC.OUT,"paper/figure-change-expr-number-prey-chromatin.tiff"),height=8,width=8)
write.table(df,file=file.path(CD4CHIC.OUT,"paper/figure-change-expr-number-prey-chromatin.csv"))

######################################################################

## * Expr, diffchic and chromatin

## ** What chromatin states predict chic.up, chic.down?

v2 <- grep("_Act$|_NAct$|_up$|_down$|^TE$|^SE$",names(nonb2b),value=TRUE) %>% setdiff(.,c("H3K9me3_Act","H3K9me3_NAct"))
up <- down <- vector("list",length(v2))
for(v in v2) {
    f <- as.formula(paste("chic.up ~ 1",v,sep="+"))
    up[[v]] <- lmodel(f,nonb2b)
    f <- as.formula(paste("chic.down ~ 1",v,sep="+"))
    down[[v]] <- lmodel(f,nonb2b)    
}
do.call("rbind",up)
do.call("rbind",down)

##' best predictor of gaining prey is H3K27ac_up
##' best predictors of losing prey are H3K4me3_up, H3K4me1_down and H3K27me3_NAct

rp.up <- rpart(paste("as.factor(chic.up) ~ ",covars) %>% as.formula(),
                      data=nonb2b,method="class")
rp.down <- rpart(paste("as.factor(chic.down) ~ ",covars) %>% as.formula(),
                      data=nonb2b,method="class")

library(glmnet)
library(lars)
library(covTest)

v3 <- grep("prom",v2,invert=TRUE,value=TRUE)
NFOLDS <- 10

LASSO <- list()

myfit <- function(x,y) {
    x[is.na(x)] <- 0
    a <- lars(x,y)
    b <- covTest(a,x,y)
    n.target <- which(b$results[,"P-value"]>0.1) %>% min()    
    gfit <- glmnet(x,y,family="binomial")
    cf <- coef(gfit) %>% as.matrix()
    n.sel <- apply(cf!=0,2,sum)
    cf[,which(n.sel==n.target) %>% min()]
}

tmp <- unique(nonb2b[logFC.expr>0,],by=c("oeID","chic.up","chic.down"))
x <- tmp[,v3,with=FALSE] %>% as.matrix()
yup <- tmp$chic.up
ydown <- tmp$chic.down
LASSO$chicup.eup <- myfit(x,yup)
LASSO$chicdown.eup <- myfit(x,ydown)

tmp <- unique(nonb2b[logFC.expr<0,],by=c("oeID","chic.up","chic.down"))
x <- tmp[,v3,with=FALSE] %>% as.matrix()
yup <- tmp$chic.up
ydown <- tmp$chic.down
LASSO$chicup.edown <- myfit(x,yup)
LASSO$chicdown.edown <- myfit(x,ydown)

tmp <- unique(nonb2b,by=c("oeID","chic.up","chic.down"))
x <- tmp[,v3,with=FALSE] %>% as.matrix()
yup <- tmp$chic.up
ydown <- tmp$chic.down
LASSO$chicup <- myfit(x,yup)
LASSO$chicdown <- myfit(x,ydown)

X <- do.call("cbind",LASSO)
#X <- lapply(LASSO,coef) %>% do.call("cbind",.) 
colnames(X) <- names(LASSO)
## print(X[-1,grep("1se",colnames(X))],col.names=TRUE)
## print(X[-1,grep("min",colnames(X))],col.names=TRUE)

library(pheatmap)
XX <- X[-c(1,grep("prom",rownames(X))),
#        intersect(grep("1se",colnames(X)), grep("eup|edown",colnames(X)))] %>% as.matrix
#        grep("1se",colnames(X))] %>% as.matrix
        ] %>% as.matrix()
nm1 <- ifelse(grepl("chicup",colnames(XX)),"CHi-C Gain","CHi-C Loss")
nm2 <- ifelse(grepl("eup",colnames(XX)),"Expr up",
              ifelse(grepl("edown",colnames(XX)),"Expr down",""))
colnames(XX) <- paste(nm1,nm2,sep="\n")
v <- apply(XX,1,var)
XX[XX>2] <- 2
XX[XX< -2] <- -2
library(RColorBrewer)
psimmap <- function(X,lim=NULL,...) {
    cols <- rev(colorRampPalette(brewer.pal(n = 7, name ="RdYlBu"))(99))
    cols[50] <- "#ffffff"
    if(is.null(lim))
        lim <- max(c(abs(X),0.001))
    breaks=c(seq(-lim,-0.001,length=50),seq(0.001,lim,length=50))
    use <- breaks>=min(X) & breaks<=max(X)
    pheatmap(X,...,
             color=cols[use],
             breaks=breaks[use])
}
psimmap(XX[v>0,],cluster_cols=FALSE,file=file.path(CD4CHIC.OUT,"paper/figure-oe-gain-loss-vs-chipseq.pdf"),height=8,width=6)
psimmap(XX[v>0,],cluster_cols=FALSE,file=file.path(CD4CHIC.OUT,"paper/figure-oe-gain-loss-vs-chipseq.tiff"),height=8,width=6)
write.table(XX[v>0,],file=file.path(CD4CHIC.OUT,"paper/figure-oe-gain-loss-vs-chipseq.csv"))

## ** do any of these alter effect of prey on expression?

tmp <- copy(nonb2b)

v4 <- setdiff(v3,c("H3K27me3_down","H3K27me3_NAct","enh_poised_NAct"))
MODELS <- vector("list",length(v4))
names(MODELS) <- v4
for(v in v4) {
    tmp[,paste0(v,c(".sig.up",".sig.down")) := list(sig.up & tmp[[v]],
                                                     sig.down & tmp[[v]])]
    changes <- c(".sig.up",".sig.down")
    if(v %in% c("enh_poised_Act","H3K27me3_up","H3K4me1_up"))
        changes <- ".sig.up"
    if(v %in% c("H3K27me3_down"))
        changes <- ".sig.down"
    MODELS[[v]] <- paste0("logFC.expr ~ sig.up + sig.down + ",v,changes)
}
modelcmp(unlist(MODELS),
         tmp)

## ** diff chic + chromhmm
v1 <- c("chic.up","chic.down")
v2 <- grep("^E",names(nonb2b),value=TRUE) %>% setdiff(.,c("E13.act","E13.non","E9.non","up.H3K9me3_1","down.H3K9me3_1","avgZ.H3K9me3_1"))
tmp <- nonb2b
terms <- character()

BICS <- vector("list",length(v2))
names(BICS) <- v2
for(v in v2) {
    message(v)
   tmp <- subset(nonb2b,complete.cases(sig.up,sig.down))
    mm <- median(tmp[[v]],na.rm=TRUE)
    tmp[,c("v.up","v.down") := list(chic.up * tmp[[v]]>mm,
                                    chic.down * tmp[[v]]>mm)]
     f0 <- c("logFC.expr ~ 1", v1) %>% paste(.,collapse="+") %>% as.formula()
f1 <- c("logFC.expr ~ 1", v1,c("v.up","v.down")) %>% paste(.,collapse="+") %>% as.formula()
BICS[[v]] <- modelcmp(list(f1,f0), tmp, do.print=FALSE) %>% print() 
}
BICS <- do.call("rbind",BICS)
BICS <- BICS[order(BICS[,1]),]
(BICS)[ BICS[,2] - BICS[,1] > 50, ]

for(va in grep("act",rownames(BICS)[ BICS[,2] - BICS[,1] > 50],value=TRUE)) {
    for(vn in grep("non",rownames(BICS)[ BICS[,2] - BICS[,1] > 50],value=TRUE)) {
        tmp <- subset(nonb2b,complete.cases(sig.up,sig.down))
        mma <- median(tmp[[va]],na.rm=TRUE)
        mmn <- median(tmp[[vn]],na.rm=TRUE)
        tmp[,c("v.up","v.down") := list(chic.up * (tmp[[va]] > mma & tmp[[vn]] > mmn),
                                    chic.down * (tmp[[va]] > mma & tmp[[vn]] > mmn))]
        f1 <- c("logFC.expr ~ 1", v1,c("v.up","v.down")) %>% paste(.,collapse="+") %>% as.formula()
        ret <- tryCatch(model(f1, tmp), error=function(e) NULL)
        if(!is.null(ret)) {
            message("\n",vn,"\t",va)
            print(ret)
        }
        
    }
}

tmp <- copy(nonb2b)
tmp[,c("chic.up","chic.down"):=list(sum(chic.up),sum(chic.down)),by="baitID"]
tmp <- unique(tmp,by="baitID")
ggplot(tmp[chic.down<=5,], aes(x=log(chic.up+1),y=logFC.expr,col=factor(chic.down))) + geom_point() + geom_smooth(method="lm")
ggplot(tmp[chic.up<=5,], aes(x=log(chic.down+1),y=logFC.expr,col=factor(chic.up),fill=factor(chic.up)))  + geom_smooth(method="lm")

modelcmp(models, subset(nonb2b,complete.cases(sig.up,sig.down,HL.328)))

ff <- paste("logFC.expr ~ ",
           paste(c(paste0(v2,".up"), paste0(v2,".down")),collapse=" + ")) %>% as.formula()
ffup <- paste("logFC.expr ~ chic.down +",
           paste(c(paste0(v2,".up")),collapse=" + ")) %>% as.formula()
ffdown <- paste("logFC.expr ~ chic.up +",
           paste(c(paste0(v2,".down")),collapse=" + ")) %>% as.formula()
model(ff, tmp[Total_CD4_NonActivated>5 | Total_CD4_Activated>5,])
model(ffup, tmp[Total_CD4_NonActivated>5 | Total_CD4_Activated>5,])
model(ffdown, tmp[Total_CD4_NonActivated>5 | Total_CD4_Activated>5,])

modelcmp(list(logFC.expr ~ chic.up + chic.down,
         ffup,
         ffdown,
         ff),
         tmp[Total_CD4_NonActivated>5 | Total_CD4_Activated>5,])

modelcmp(logFC.expr ~ chic.up + chic.down,
         ff,
         tmp[Total_CD4_NonActivated>5 | Total_CD4_Activated>5,])


outer(v1,v2,paste,sep=":") %>% as.vector()

ff <- paste("logFC.expr ~ ",
paste(c(v1,v2),collapse="+"),
"+",
(outer(v1,v2,paste,sep=":") %>% as.vector() %>% paste(.,collapse="+"))) %>%
as.formula()
ff

model(ff, nonb2b[Total_CD4_NonActivated>5 | Total_CD4_Activated>5,])



## check method
lm(logFC.expr ~ log(HL.328) + log(HL.M), data=expr, subset=ngenesperbait==1) %>% summary()
model(logFC.expr ~ log(HL.328) + log(HL.M), data)

## looks ok
model(logFC.expr ~ chic.up + chic.down, nonb2b)

model(logFC.expr ~ chic.up*ac27.any + chic.down*ac27.any, nonb2b)
model(logFC.expr ~ chic.up*K4me3.any + chic.down*K4me3.any, nonb2b)
model(logFC.expr ~ chic.up*K4me1.any + chic.down*K4me1.any, nonb2b)
model(logFC.expr ~ chic.up*ac27.any + chic.down*ac27.any+ chic.up*me27.any + chic.down*me27.any, nonb2b)



model(logFC.expr ~ sig.up + sig.down, nonb2b)
model(logFC.expr ~ sig.up + sig.down + HL.M + HL.328, nonb2b)
model(logFC.expr ~ sig.up + sig.down + HL.328, nonb2b)

nonb2b[,c("sig","dir"):=list(FDR.chic<0.1, logFC.chic>0)]
model(logFC.expr ~ sig*up, nonb2b)

ggplot(expr,aes(x=HL.328/60,y=logFC.expr)) + geom_point() + geom_smooth() + scale_x_log10() + xlab("Half life in activated cells (hours)") + ylab("log2 FC expr")
ggsave("hl-expr.pdf",height=5,width=6)
system("scp hl-expr.pdf chrisw@salthouse:~/Words/talks/labtalk-dec15")


model(logFC.expr ~ logFC.chic*enh.any + log(HL.328) + log(HL.M), nonb2b)
model(logFC.expr ~ logFC.chic*enh.any + log(HL.328) + log(HL.M), nonb2b)


## is logFC predicted by chrom states at oe at baseline? NO.

model(logFC.expr ~ ac27.any, nonb2b)
anys <- grep("any",names(nonb2b),value=TRUE)
ff <- paste("logFC.expr ~ ", paste(anys,collapse="+")) %>% as.formula()
model(ff, nonb2b[Total_CD4_NonActivated>5,])
model(ff, nonb2b[Total_CD4_Activated>5,])

## is logFC predicted by chrom states at differential interactions?

model(logFC.expr ~ ac27.any, nonb2b)
anys <- grep("any",names(nonb2b),value=TRUE)
ff <- paste("logFC.expr ~ logFC.chic * ", paste(anys,collapse="+")) %>% as.formula()
model(ff, nonb2b[Total_CD4_NonActivated>5 | Total_CD4_Activated>5,])

v1 <- c("chic.up","chic.down")
v2 <- anys
outer(v1,v2,paste,sep=":") %>% as.vector()

ff <- paste("logFC.expr ~ ",
paste(c(v1,v2),collapse="+"),
"+",
(outer(v1,v2,paste,sep=":") %>% as.vector() %>% paste(.,collapse="+"))) %>%
as.formula()
ff

model(ff, nonb2b[Total_CD4_NonActivated>5 | Total_CD4_Activated>5,])

## gene level:


thr <- 0.1
genes <- m[id!="." & !is.na(id) & !b2b & ngenesperbait==1, ] %>%
 group_by(., id) %>%
    summarise(.,
              act.ints=count5(Total_CD4_Activated),
              non.ints=count5(Total_CD4_NonActivated),
              int.up=ndiff5(Total_CD4_Activated, Total_CD4_NonActivated),
              int.down=ndiff5(Total_CD4_NonActivated, Total_CD4_Activated),
              ## n1.up=cond.ndiff5(Total_CD4_Activated, Total_CD4_NonActivated,nonClass==1),
              ## n1.down=cond.ndiff5(Total_CD4_NonActivated, Total_CD4_Activated,nonClass==1),
              ## n2.up=cond.ndiff5(Total_CD4_Activated, Total_CD4_NonActivated,nonClass==2),
              ## n2.down=cond.ndiff5(Total_CD4_NonActivated, Total_CD4_Activated,nonClass==2),
              ## n3.up=cond.ndiff5(Total_CD4_Activated, Total_CD4_NonActivated,nonClass==3),
              ## n3.down=cond.ndiff5(Total_CD4_NonActivated, Total_CD4_Activated,nonClass==3),
              ## n4.up=cond.ndiff5(Total_CD4_Activated, Total_CD4_NonActivated,nonClass==4),
              ## n4.down=cond.ndiff5(Total_CD4_NonActivated, Total_CD4_Activated,nonClass==4),
              ## a1.up=cond.ndiff5(Total_CD4_Activated, Total_CD4_NonActivated,actClass==1),
              ## a1.down=cond.ndiff5(Total_CD4_NonActivated, Total_CD4_Activated,actClass==1),
              ## a2.up=cond.ndiff5(Total_CD4_Activated, Total_CD4_NonActivated,actClass==2),
              ## a2.down=cond.ndiff5(Total_CD4_NonActivated, Total_CD4_Activated,actClass==2),
              ## a3.up=cond.ndiff5(Total_CD4_Activated, Total_CD4_NonActivated,actClass==3),
              ## a3.down=cond.ndiff5(Total_CD4_NonActivated, Total_CD4_Activated,actClass==3),
              ## a4.up=cond.ndiff5(Total_CD4_Activated, Total_CD4_NonActivated,actClass==4),
              ## a4.down=cond.ndiff5(Total_CD4_NonActivated, Total_CD4_Activated,actClass==4),
              sig.up=nasum(logFC.chic>0 & FDR.chic<thr),
              sig.down=nasum(logFC.chic<0 & FDR.chic<thr),
              bsig.up=nasum(logFC.chic>1 & FDR.chic<thr),
              bsig.down=nasum(logFC.chic< -1 & FDR.chic<thr))
dim(genes)
setkey(genes,id)

## add modules
## mods <- get.modules()
## setkey(mods,id)
## mods <- unique(mods)
## genes <- merge(genes,mods,all.x=TRUE)

## add rnaseq
genes <- merge(genes,expr[,.(id,logFC.expr,FDR.expr,HL.M,HL.328)],by="id")
ggplot(genes, aes(x=sig.up,y=logFC.expr)) + geom_point() + geom_smooth(method="lm") + facet_wrap(~sig.down)

library(pheatmap)
tmp <- genes
tmp[,id:=NULL]
as.matrix(tmp) %>% cor(., use="pair") %>% pheatmap()

## rnadecay
xx <- get.rnadecay()
setkey(xx,id)
genes <- merge(genes,xx,all.x=TRUE)

genes[,c("y.act","y.non"):=list(log2(expr.reads.act+1),log2(expr.reads.non+1))]
genes$adjHL.M <- scale(log(genes$HL.M),scale=FALSE,center=TRUE)
genes$adjHL.3 <- scale(log(genes$HL.3),scale=FALSE,center=TRUE)
genes$adjHL.328 <- scale(log(genes$HL.328),scale=FALSE,center=TRUE)

lm(y.act ~ y.non + int.up + int.down, data=genes) %>% summary()
lm(y.act ~ y.non + int.up + int.down + adjHL.M*adjHL.328, data=genes) %>% summary()
lm(logFC.expr ~ int.up + int.down + adjHL.M*adjHL.328, data=genes) %>% summary()
lm(logFC.expr ~ int.up + int.down + sig.up + sig.down + bsig.up + bsig.down, data=genes) %>% summary()
lm(logFC.expr ~ sig.up + sig.down + adjHL.M*adjHL.328, data=genes) %>% summary()
lm(logFC.expr ~ int.up + int.down + sig.up + sig.down + adjHL.M*adjHL.328, data=genes) %>% summary()
lm(logFC.expr ~ int.up + int.down + sig.up + sig.down + adjHL.M*adjHL.328, data=genes) %>% summary()
lm(logFC.expr ~ int.up + int.down + sig.up + sig.down + adjHL.M+adjHL.328, data=genes) %>% summary()
lm(logFC.expr ~ int.up + int.down + sig.up + sig.down + bsig.up + bsig.down, data=genes) %>% summary()
lm(logFC.expr ~ sig.up + sig.down, data=genes) %>% summary()
lm(logFC.expr ~ n1.up + n1.down + sig.up + sig.down, data=genes) %>% summary()
lm(logFC.expr ~ int.up + int.down + n1.up + n1.down + sig.up + sig.down, data=genes) %>% summary()

genes[,.(y.act,y.non,logFC.expr,int.up,int.down,sig.up,sig.down,adjHL.M,adjHL.328)] %>% as.matrix() %>% cor(.,use="pair")



lm(y.non ~ y.act + int.up + int.down + adjHL.M*adjHL.328, data=genes) %>% summary()
lm(y.act ~ y.non + sig.up + sig.down, data=genes) %>% summary()
lm(y.act ~ y.non + int.up + int.down + sig.up + sig.down, data=genes) %>% summary()
lm(y.act ~ y.non + int.up + int.down + sig.up + sig.down + bsig.up + bsig.down, data=genes) %>% summary()
lm(y.act ~ y.non + sig.up + sig.down + bsig.up + bsig.down, data=genes) %>% summary()
lm(logFC.expr ~ sig.up + sig.down + bsig.up + bsig.down -1, data=genes) %>% summary()
lm(logFC.expr ~ sig.up + sig.down -1, data=genes) %>% summary()
lm(logFC.expr ~ sig.up * sig.down -1, data=genes) %>% summary()

## check as FDR thr changes
thr <- 10^seq(-1,-6,by=-0.5)[-1]
thr <- seq(0.5,5,by=0.5)
RESULTS <- vector("list",length(thr))
for(i in seq_along(thr)) {
tmpgenes <- m[id!="." & !is.na(id) & !b2b & ngenesperbait==1, ] %>%
 group_by(., id) %>%
    summarise(.,
              act.ints=count5(Total_CD4_Activated),
              non.ints=count5(Total_CD4_NonActivated),
              int.up=ndiff5(Total_CD4_Activated, Total_CD4_NonActivated),
              int.down=ndiff5(Total_CD4_NonActivated, Total_CD4_Activated),
              sig.up=nasum(logFC.chic>0 & FDR.chic<0.1),
              sig.down=nasum(logFC.chic<0 & FDR.chic<0.1),
              ## bsig.up=nasum(logFC.chic>0 & FDR.chic<thr[[i]]),
              ## bsig.down=nasum(logFC.chic< 0 & FDR.chic<thr[[i]]))
              bsig.up=nasum(logFC.chic>thr[[i]] & FDR.chic<0.1),
              bsig.down=nasum(logFC.chic< -thr[[i]] & FDR.chic<0.1))
      setkey(rnaseq, id)
      setkey(tmpgenes,id)
      tmpgenes <- merge(tmpgenes,
               rnaseq, all.x=TRUE)
    modsum <- lm(logFC.expr ~ sig.up + sig.down + bsig.up + bsig.down -1, data=tmpgenes) %>% summary()
tmp <- cbind(modsum$coefficients,thr=thr[[i]]) %>% as.data.frame()
tmp$var=rownames(modsum$coefficients)
RESULTS[[i]] <- tmp
}

RESULTS2 <- do.call("rbind",RESULTS) %>% as.data.frame()
ggplot(RESULTS2,aes(x=log10(thr),y=Estimate)) + geom_point() + facet_wrap(~var)
ggplot(RESULTS2,aes(x=thr,y=Estimate)) + geom_point() + facet_wrap(~var)

lapply(RESULTS,"[","Estimate")


ggplot(subset(genes,sig.down==0 & sig.up<=10 & logFC.expr<5), aes(x=sig.up,y=logFC.expr)) + geom_point() + geom_smooth() + geom_smooth(method="lm",col="green")

ggplot(subset(genes,sig.up<=10 & logFC.expr<5), aes(x=sig.up,y=logFC.expr)) + geom_point() + geom_smooth() + geom_smooth(method="lm",col="green") + facet_wrap(~sig.down)

ggplot(subset(genes,sig.up==0 & sig.up<=10 & logFC.expr<5), aes(x=sig.down,y=logFC.expr)) + geom_point() + geom_smooth() + geom_smooth(method="lm",col="green")

ggplot(subset(genes,sig.up<=10 & logFC.expr<5), aes(x=sig.down,y=logFC.expr)) + geom_point() + geom_smooth() + geom_smooth(method="lm",col="green") + facet_wrap(~sig.up)



mod <- lm(logFC.expr ~ sig.up + sig.down - 1, data=genes)
summary(mod)

kk

## ggplot(subset(genes,!is.na(logFC.expr)),aes(x=sig.up,y=sig.down,col=logFC.expr)) + geom_jitter() + scale_colour_gradient2(low="red",mid="grey",high="green",midpoint=0)
## ggplot(subset(genes,!is.na(logFC.expr)),aes(x=sig.up,y=logFC.expr)) + geom_point() + facet_wrap(~sig.down) + geom_smooth(method="lm")
## ggplot(subset(genes,!is.na(logFC.expr)),aes(x=sig.down,y=logFC.expr)) + geom_point() + facet_wrap(~sig.up) + geom_smooth(method="lm")

## library(scatterplot3d)
## genes[, pcol:=ifelse(logFC.expr>0, ifelse(FDR.expr<0.1, "darkgreen", "green"),
##             ifelse(FDR.expr<0.1, "red", "pink")) ]
## sp <- with(subset(genes,!is.na(logFC.expr)), {
##    scatterplot3d(sig.up,   # x axis
##                  -sig.down,     # y axis
##                  logFC.expr,    # z axis
##                  type="p",
##                  pch=19,
##                  color=pcol,
##                  main="3-D Scatterplot Example 1",
##                  zlim=c(-6,12),
##                  grid=FALSE)
## })
## sp$plane3d(0,coefficients(mod)[1], -coefficients(mod)[2],col="lightblue")
## sp$plane3d(0,0,0,col="grey")

## library(lattice)
## with(subset(genes,!is.na(logFC.expr)), {
##     cloud(logFC.expr ~ sig.up * sig.down, type = c('p','h'), pch=16,
##           scales=list(arrows=FALSE),
##           screen = list(z = 50, x = -70, y = -10),
##           col=pcol)
## })

## ## =======================================================================
## ## scatter3D with fitted surface : the mtcars dataset
## ## implemented by Karline Soetaert
## ## =======================================================================
 
## require(plot3D)
 
## # predict on x-y grid, for surface
## nup <- max(genes$sig.up)
## ndown <- max(genes$sig.down)
## xp <- seq(0,nup,by=5)
## yp <- seq(0,ndown,by=5)
## xy <- expand.grid(sig.up = xp, 
##                  sig.down = yp)
 
## pred <- matrix (nrow = length(xp), ncol = length(yp), 
##   data = predict(mod, newdata = data.frame(xy), interval = "prediction"))
 
## # predicted z-values, fitted points for droplines to surface

## make.ci <- function(z) {
##     cbind(ifelse(z<0, z, 0),
##           ifelse(z<0, 0, z))
## }

## with(subset(genes,!is.na(logFC.expr)), {
## wh <- c(which.max(logFC.expr),which.min(logFC.expr),which.min(sig.up),which.max(sig.up),which.min(sig.down),which.max(sig.down))
## scatter3D(#z = logFC.expr[wh], x = sig.up[wh], y = sig.down[wh],
##           z = logFC.expr, x = sig.up, y = sig.down, 
##           pch = 16, cex = 1, 
##       theta = 20, phi = 20, ticktype = "detailed",
##       xlab = "nup", ylab = "ndown", zlab = "log2 FC", #clab = "mpg", 
##       ## surf = list(x = xp, y = yp, z = pred, 
##       ##                 facets = NA),
## #          CI=list(z=make.ci(logFC.expr[wh]),alen=0),
##    bty = "g",
##           ## surf = list(x = seq(0,nup), y = seq(0,ndown), z = matrix(0,nup+1,ndown+1), 
##           ##     facets = NA, fit = fitpoints),
##       colkey = list(length = 0.8, width = 0.4))
## scatter3D(z = logFC.expr[wh], x = sig.up[wh], y = sig.down[wh], 
##           pch = 16, 
##       theta = 20, phi = 20, ticktype = "detailed",
##       surf = list(x = xp, y = yp, z = pred, #matrix(0,length(xp),length(yp)),
##                       facets = NA),
##           add=TRUE)
## })


## scatter3D(pch = 16, cex = 1, 
##       theta = 20, phi = 20, ticktype = "detailed",
##       xlab = "wt", ylab = "disp", zlab = "log2 FC", #clab = "mpg", 
##       surf = list(x = seq(0,nup), y = seq(0,ndown), z = pred, 
##                       facets = NA),
## add=TRUE,
##           ## surf = list(x = seq(0,nup), y = seq(0,ndown), z = matrix(0,nup+1,ndown+1), 
##           ##     facets = NA, fit = fitpoints),
##       colkey = list(length = 0.8, width = 0.4))
## })

## with(subset(genes,!is.na(logFC.expr)), {
## scatter3D(z = logFC.expr, x = sig.up, y = sig.down, pch = 16, cex = 1, 
##       theta = 20, phi = 20, ticktype = "detailed", type="h")

##       xlab = "wt", ylab = "disp", zlab = "log2 FC", #clab = "mpg", 
##       ## surf = list(x = seq(0,nup), y = seq(0,ndown), z = pred, 
##           ##             facets = NA),
##           surf = list(x = seq(0,nup), y = seq(0,ndown), z = matrix(0,nup+1,ndown+1), 
##               facets = NA, fit = fitpoints),
##       colkey = list(length = 0.8, width = 0.4))
## })


c1 <- plotcor(x=log2(genes$non.ints+1), y=log2(genes$expr.reads.non+1), y.label="top")
c2 <- plotcor(x=log2(genes$act.ints+1), y=log2(genes$expr.reads.act+1), y.label="top")
c3 <- with(genes[!is.na(logFC.expr),], plotcor(x=int.fc, y=logFC.expr, xmax=4, y.label="bottom"))
## dropping 57/6978 observations above 4
## full range of data is: 0.0769230769230769 - 23

## p <- ggplot(genes, aes(x=logFC, y=log2(int.fc))) + geom_point() + geom_vline(xintercept=0,col="grey",size=0.5) + geom_hline(xintercept=0,col="grey",size=0.5) + geom_smooth(se=FALSE,method="lm") +
##     labs(x="expression fold change (log2 scale)", y="fold change in number of interacting fragments (log2 scale)") + colScale + theme(legend.position="none")
## p
## png(file.path(CD4CHIC.OUT, "comparisons", "fig-exprFC-nintFC.png"),res=300,units="in",height=6,width=6)
## print(p)
## dev.off()

## png(file.path(CD4CHIC.OUT, "comparisons", "fig-exprFC-nintFC-bymodule.png"),res=300,height=6,width=6,units="in")
## ggplot(genes, aes(x=logFC, y=log2(int.fc))) + geom_point(aes(col=module)) + geom_vline(xintercept=0,col="grey",size=0.5) + geom_hline(xintercept=0,col="grey",size=0.5) + geom_smooth(se=FALSE,method="lm") +
##     labs(x="expression fold change (log2 scale)", y="fold change in number of interacting fragments (log2 scale)") + facet_wrap(~module) + colScale + theme(legend.position="none")
## dev.off()


################################################################################


## * Does the strength of interaction change with gene expression upon activation for OEs present in both?

################################################################################

## collapse to one row per gene, counting interactions per gene
setkey(m,id)
genes.5 <- m[id!="." & !is.na(id) & !is.na(logFC.chic) & Total_CD4_Activated>=5 & Total_CD4_NonActivated>=5,
             c("n.stable","act.score5.50","non.score5.50",
               "sfc5.50","sfc5.25","sfc5.75","sfc5.0","sfc5.100",
               "fc5.50","fc5.25","fc5.75","fc5.0","fc5.100"
               ) :=
             list(sum(Total_CD4_Activated>=5),
                 median(Total_CD4_Activated), median(Total_CD4_NonActivated),
                  median(log(Total_CD4_Activated/Total_CD4_NonActivated),na.rm=TRUE),
                  quantile(log(Total_CD4_Activated/Total_CD4_NonActivated),0.25,na.rm=TRUE),
                  quantile(log(Total_CD4_Activated/Total_CD4_NonActivated),0.75,na.rm=TRUE),
                  min(log(Total_CD4_Activated/Total_CD4_NonActivated),na.rm=TRUE),
                  max(log(Total_CD4_Activated/Total_CD4_NonActivated),na.rm=TRUE),
                  median(logFC.chic,na.rm=TRUE),
                  quantile(logFC.chic,0.75,na.rm=TRUE),
                  quantile(logFC.chic,0.75,na.rm=TRUE),
                  min(logFC.chic,na.rm=TRUE),
                  max(logFC.chic,na.rm=TRUE)
                  ),
             by=id] %>% unique()
genes.a <- m[id!="." & !is.na(id) & !is.na(logFC.chic),
             c("n.act","n.non","act.score.50","non.score.50",
               "sfc.50","sfc.25","sfc.75","sfc.0","sfc.100",
               "fc.50","fc.25","fc.75","fc.0","fc.100"
               ) :=
             list(sum(Total_CD4_Activated>=5),
                 sum(Total_CD4_NonActivated>=5),
                 median(Total_CD4_Activated), median(Total_CD4_NonActivated),
                  median(log(Total_CD4_Activated/Total_CD4_NonActivated),na.rm=TRUE),
                  quantile(log(Total_CD4_Activated/Total_CD4_NonActivated),0.25,na.rm=TRUE),
                  quantile(log(Total_CD4_Activated/Total_CD4_NonActivated),0.75,na.rm=TRUE),
                  min(log(Total_CD4_Activated/Total_CD4_NonActivated),na.rm=TRUE),
                  max(log(Total_CD4_Activated/Total_CD4_NonActivated),na.rm=TRUE),
                  median(logFC.chic,na.rm=TRUE),
                  quantile(logFC.chic,0.75,na.rm=TRUE),
                  quantile(logFC.chic,0.75,na.rm=TRUE),
                  min(logFC.chic,na.rm=TRUE),
                  max(logFC.chic,na.rm=TRUE)
                  ),
             by=id] %>% unique()
setkey(genes.5,id)
setkey(genes.a,id)
genes <- merge(genes.5[,.(id,n.stable,act.score5.50,non.score5.50,
               sfc5.50,sfc5.25,sfc5.75,sfc5.0,sfc5.100,
               fc5.50,fc5.25,fc5.75,fc5.0,fc5.100)],
               genes.a[,.(id,n.act,n.non,act.score.50,non.score.50,
               sfc.50,sfc.25,sfc.75,sfc.0,sfc.100,
               fc.50,fc.25,fc.75,fc.0,fc.100)])

## add rnaseq
setkey(genes,id)
setkey(rnaseq,id)
genes <- merge(genes,
               rnaseq[,.(id, FDR.expr, logFC.expr)],
                   all.x=TRUE, by="id")
dim(genes)
genes$sup <- genes$sfc5.50>0
genes$up <- genes$fc5.50>0


mg <- melt(genes, c("id","n.stable","FDR.expr","logFC.expr","n.act","n.non","sup","up"))
mg <- subset(mg,!is.na(value) & is.finite(value))
mg$n.stable[is.na(mg$n.stable)] <- 0
mg$cut.stable <- cut(mg$n.stable, c(0, 3, 10, 200))

       
ggplot(subset(mg, grepl("sfc5",mg$variable) & n.stable>0),
       aes(y=value,x=logFC.expr,col=variable)) + geom_point(size=1,alpha=0.1) + geom_hline(yintercept=0) + geom_vline(xintercept=0) + geom_smooth(se=FALSE) + geom_rug(alpha=0.1) +
    facet_grid(sup~cut.stable)

ggplot(subset(mg, grepl("^fc5",mg$variable) & n.stable>0 & !is.na(up)),
       aes(y=value,x=logFC.expr,col=variable)) + geom_point(size=1,alpha=0.1) + geom_hline(yintercept=0) + geom_vline(xintercept=0) + geom_smooth(se=FALSE,method="lm") + geom_rug(alpha=0.1) +
    facet_grid(up~cut.stable) + ylim(-1,1)

mg$diff.cut <- cut(mg$n.act - mg$n.non, c(-150, -10, -2, -1,0,1,2, 10, 150))

ggplot(subset(mg, grepl("sfc\\.",mg$variable) & n.stable==0 &
              logFC.expr<5 & !is.na(value) & is.finite(value) & value>0),
       aes(y=value,x=logFC.expr,col=variable)) + geom_point(size=1,alpha=0.1) + geom_hline(yintercept=0) + geom_vline(xintercept=0) + geom_smooth(se=FALSE) +
   facet_grid(sup~diff.cut) + ylim(-2,2)

ggplot(subset(mg, grepl("^fc\\.",mg$variable) & !is.na(up)),
       aes(y=value,x=logFC.expr,col=variable)) + geom_point(size=1,alpha=0.1) + geom_hline(yintercept=0) + geom_vline(xintercept=0) + geom_smooth(se=FALSE) + geom_rug(alpha=0.1) +
    facet_wrap(~diff.cut)


       

## genes <- subset(m, id!="." & !is.na(id) & Total_CD4_Activated>=5 & Total_CD4_NonActivated>=5)
## genes <- group_by(genes, id) %>%
##     summarise(.,
##               act.int=mean(Total_CD4_Activated), non.int=mean(Total_CD4_NonActivated),
##               score.fc=median(Total_CD4_Activated/Total_CD4_NonActivated,na.rm=TRUE),
##               med.logFC.chic=median(logFC.chic,na.rm=TRUE))

## ## add modules
## ## SHOULD DO THIS BY ENSEMBL ID
## genes <- merge(genes,modules[,c("ensGeneId","geneSymbol","module")], by.x="id", by.y="ensGeneId", all.x=TRUE)
## dim(genes)
## dim(modules)
## #mod.ints <- split(gene.ints, gene.ints$module)

mycols <- unique(as.character(genes$module))
names(mycols) <- mycols
mycols[mycols=="NA"] <- "grey20"
colScale <- scale_colour_manual(name = "module",values = mycols, na.value="purple")


genes <- subset(genes, !is.na(FDR.expr))
## genes$int.fc <- (genes$act.int+1)/(genes$non.int+1)
## genes$int.delta <- genes$act.int-genes$non.int

c4 <- with(genes, plotcor(x=int.fc, y=logFC.expr,xmax=2))
c4$plot
with(genes, plotcor(x=med.score.fc, y=logFC.expr))
with(genes, plotcor(x=med.score.fc, y=logFC.expr,xmax=2))
with(genes, plotcor(x=log(int.fc), y=logFC.expr,xmax=2))
with(genes, plotcor(x=log(score.fc), y=logFC.expr,xmax=2))
with(genes[!is.na(med.logFC.chic),], plotcor(x=med.logFC.chic, y=logFC.expr,xmax=2))

p1<-with(genes[med.logFC.chic<0,], plotcor(x=exp(med.logFC.chic), y=logFC.expr,xmax=2))
p2<-with(genes[med.logFC.chic>0,], plotcor(x=exp(med.logFC.chic), y=logFC.expr,xmax=2))
library(gridExtra)
grid.arrange(p1,p2)
pdf("~/kk.pdf")
grid.arrange(p1,p2)
dev.off()


with(genes[int.fc<1 & !is.na(med.logFC.chic),], bootcor(x=log(int.fc),y=logFC.expr))
with(genes[int.fc>1 & !is.na(med.logFC.chic),], bootcor(x=log(int.fc),y=logFC.expr))

with(genes[score.fc<1 & !is.na(med.logFC.chic),], bootcor(x=log(score.fc),y=logFC.expr))
with(genes[score.fc>1 & !is.na(med.logFC.chic),], bootcor(x=log(score.fc),y=logFC.expr))

with(genes[med.logFC.chic<0 & !is.na(med.logFC.chic),], bootcor(x=med.logFC.chic,y=logFC.expr))
with(genes[med.logFC.chic>0 & !is.na(med.logFC.chic),], bootcor(x=med.logFC.chic,y=logFC.expr))


c1$plot <- c1$plot + labs(x="Log2 number of interactions",y="Log2 read depth")
c2$plot <- c2$plot + labs(x="Log2 number of interactions",y="Log2 read depth")
c3$plot <- c3$plot +  labs(x="Ratio of number of interacting regions",y="Log2 fold change in expression")
c4$plot <- c4$plot + labs(x="Ratio of mean interaction scores",y="Log2 fold change in expression")

library(cowplot)
png(file.path(CD4CHIC.OUT,"comparisons","fig-interactions-vs-expression.png"),res=300,height=8,width=10,units="in")
plot_grid(c1$plot,
          c2$plot,
          c3$plot,
          c4$plot,
          cols=2,
          labels=letters[1:4] %>% toupper())
dev.off()
   
png(file.path(CD4CHIC.OUT,"comparisons","fig-interactions-vs-expression.png-A"),res=300,height=8,width=10,units="in")
print(c1$plot)
dev.off()

ggsave(file.path(CD4CHIC.OUT,"comparisons","fig-interactions-vs-expression-B.png"),
       c2$plot,
       dpi=300,height=8,width=10)
ggsave(file.path(CD4CHIC.OUT,"comparisons","fig-interactions-vs-expression-C.png"),
       c3$plot,
       dpi=300,height=8,width=10)
ggsave(file.path(CD4CHIC.OUT,"comparisons","fig-interactions-vs-expression-D.png"),
       c4$plot,
       dpi=300,height=8,width=10)

lapply(list(c1,c2,c3,c4),"[[","stats")


################################################################################

## * does the strength asymetry persist if no thresholding on score?

head(diffchic)
gi2 <- function(threshold=5) {
    require(data.table)#merged_samples_12Apr2015_full_denorm_bait2baits_e75.tab
##    int <- fread(file.path(CD4CHIC.DATA,"merged_samples_12Apr2015_full_denorm_bait2baits_e75.tab")) # %>% as.data.frame()
    int <- fread(file.path(CD4CHIC.DATA,"merged_samples_12Apr2015_full.txt"))

    ## link baits to genes
    int.genes <- fread(file.path(CD4CHIC.DATA,"merged_samples_12Apr2015_full_denorm_bait2baits_e75.tab"),select=c(1:4,8))
#    nm <- scan(file.path(CD4CHIC.DATA,"merged_samples_12Apr2015_full_denorm_bait2baits.txt"), nlines=1, what="")
    setnames(int.genes, names(int.genes), c("id","gene","biotype","strand","baitID"))
    setkey(int.genes,baitID,id)
    int.genes <- unique( int.genes )
    setkey(int,baitID)
    setkey(int.genes,baitID)
    int <- merge(int,int.genes,allow.cartesian=TRUE)
    
    ## int[ oeID==120852 & baitID==120843, .(baitChr,baitStart,baitID,oeChr,oeStart,oeID,id,gene,Total_CD4_Activated,Total_CD4_NonActivated)]
    int <- int[,b2b:= oeID %in% int$baitID]
    if(!is.null(threshold)) {
        maxscore <- apply(int[,c(12:21,25:28),with=FALSE],1,max)
        int <- int[maxscore>=threshold, ]
    }
    int <- int[, c("baitChr","oeChr") := list(paste0("chr",baitChr), paste0("chr",oeChr)) ]
    int <- int[, c("baitLength","oeLength"):= list(abs(baitStart-baitEnd), abs(oeStart-oeEnd)) ]
    int <- int[,.(baitChr,baitStart,baitID,oeChr,oeStart,oeID,id,gene,Total_CD4_Activated,Total_CD4_NonActivated,baitLength,oeLength,b2b)]
    return(unique(int,by=NULL)) # sometimes one bait maps to multiple promoters for the same gene
}

ints <- get.interactions(threshold=5)

setkeyv(gene.ints, c("baitID","oeID"))
setkeyv(diffchic, c("baitID","oeID"))
m <- merge(gene.ints,
           diffchic[,.(baitID,oeID,chic.reads.non,chic.reads.act,logFC,FDR)], all.x=TRUE)
#           by.x=c("baitID","oeID"),by.y=c("bait","oe"))
setnames(m,c("logFC","FDR"), c("logFC.chic","FDR.chic"))
## summarise one row per gene
setkey(m,id)

thr <- 0
genes <- m[id!="." & !is.na(id) & !is.na(logFC.chic) & !b2b &
           (Total_CD4_Activated>=thr & Total_CD4_NonActivated>=thr) &
          FDR.chic<0.1,
             c("n.act","n.non","n.stable","act.score5.50","non.score5.50",
               "sfc5.50",
               "fc5.50"
               ) :=
             list(sum(Total_CD4_Activated>=5),
                  sum(Total_CD4_NonActivated>=5),
                  sum(Total_CD4_Activated>=5 & Total_CD4_NonActivated>=5),
                  median(Total_CD4_Activated),
                  median(Total_CD4_NonActivated),
                  median(log((Total_CD4_Activated+1)/(Total_CD4_NonActivated+1)),na.rm=TRUE),
                  median(logFC.chic,na.rm=TRUE)
                  ),
             by=id] %>% unique()
nsgenes <- m[id!="." & !is.na(id) & !is.na(logFC.chic) & !b2b &
           (Total_CD4_Activated>=thr & Total_CD4_NonActivated>=thr) &
          FDR.chic>0.1,
             c("n.act","n.non","n.stable","act.score5.50","non.score5.50",
               "sfc5.50",
               "fc5.50"
               ) :=
             list(sum(Total_CD4_Activated>=5),
                  sum(Total_CD4_NonActivated>=5),
                  sum(Total_CD4_Activated>=5 & Total_CD4_NonActivated>=5),
                  median(Total_CD4_Activated),
                  median(Total_CD4_NonActivated),
                  median(log((Total_CD4_Activated+1)/(Total_CD4_NonActivated+1)),na.rm=TRUE),
                  median(logFC.chic,na.rm=TRUE)
                  ),
             by=id] %>% unique()
genes <- genes[sfc5.50> -1.5,]
nsgenes <- nsgenes[sfc5.50> -1.5,]
setkey(genes,id)
setkey(nsgenes,id)
setkey(rnaseq,id)
genes <- merge(genes,
               rnaseq[,.(id, FDR.expr, logFC.expr)],
                   all.x=TRUE, by="id")
nsgenes <- merge(nsgenes,
               rnaseq[,.(id, FDR.expr, logFC.expr)],
                   all.x=TRUE, by="id")
plots <- list(ggplot(genes[FDR.expr<0.01,],aes(x=fc5.50,y=logFC.expr)) + geom_point() + geom_smooth() + ggtitle("FC"),
ggplot(genes[FDR.expr<0.01,],aes(x=sfc5.50,y=logFC.expr)) + geom_point() + geom_smooth() + ggtitle("score"),
ggplot(nsgenes[FDR.expr<0.01 & FDR.chic>0.1,],aes(x=fc5.50,y=logFC.expr)) + geom_point() + geom_smooth() + ggtitle("FC, nsig"),
ggplot(nsgenes[FDR.expr<0.01 & FDR.chic>0.1,],aes(x=sfc5.50,y=logFC.expr)) + geom_point() + geom_smooth() + ggtitle("score, nsig"))
library(gridExtra)
do.call(grid.arrange,plots)

pdf(file.path(CD4CHIC.OUT,"comparisons","asymmetry.pdf"),height=8,width=8)
do.call(grid.arrange,plots)
dev.off()

## only signif FC
setkey(m,id)
setkey(rnaseq,id)
m2 <- merge(m,rnaseq[,.(id, FDR.expr, logFC.expr)])
m2 <- m2[FDR.chic<0.05 & abs(logFC.chic)<5,]

ggplot(m2,aes(x=Total_CD4_NonActivated,y=Total_CD4_Activated,col=sign(logFC.chic))) + geom_point() + geom_smooth() + geom_abline()
ggplot(m2[abs(logFC.chic)<5,],aes(x=logFC.chic,y=logFC.expr)) + geom_point() + geom_smooth()
ggplot(m2[abs(logFC.chic)<5 & logFC.chic<0,],aes(x=logFC.chic,y=logFC.expr)) + geom_point() + geom_smooth()
ggplot(m2[abs(logFC.chic)<5 & logFC.chic>0,],aes(x=logFC.chic,y=logFC.expr)) + geom_point() + geom_smooth()

library(boot)
bcor <- function(x,y,by,R=200) {
    f <- function(di,i){    
        d2.idx <- di[i,]
        d2 <- d[unlist(indices[d2.idx]),]
        cor.test(d2$x, d2$y)$statistic %>% as.numeric()
    }
    d <- cbind(x=x,y=y,by=by) %>% as.data.frame()
    d <- d[!is.na(x) & !is.na(y),]
    indices <- split(1:nrow(d),d$by)
    d.idx <- data.frame(i=1:length(indices))
    bb <- boot(d.idx, f, R=R) #, parallel="multicore",ncpus=5)
    s <- sd(bb$t)
    c(rho=bb$t0, se.rho=s, p=pnorm(abs(bb$t0),sd=s,lower.tail=FALSE)*2 )
}

with(m2,bcor(logFC.chic,logFC.expr,id))
with(m2[logFC.chic<0,],bcor(logFC.chic,logFC.expr,id))
with(m2[logFC.chic>0,],bcor(logFC.chic,logFC.expr,id))


library(lme4)
m2 <- m2[!is.na(logFC.expr) & !is.na(logFC.chic) & !b2b,]
m2$group <- as.factor(m2$id)
m2$sign <- m2$logFC.chic>0

## promiscuous baits
prom <- copy(m2)
setkey(prom,baitID)
prom[,ngene:=length(unique(id)),by=baitID]
prom <- unique(prom[ngene==1,])
m2 <- m2[baitID %in% prom$baitID,]

setkey(m2,id)
with(unique(m2[logFC.chic<0,]),mean(logFC.expr)) # -0.1105277
with(unique(m2[logFC.chic>0,]),mean(logFC.expr)) # 0.8939688

m2[,allup:=all(sign(logFC.chic)>0),by=id]
m2[,alldown:=all(sign(logFC.chic)<0),by=id]
with(m2,table(allup,alldown))
m2[,chic.dir:=ifelse(allup,"up", ifelse(alldown,"down","inconsistent"))]

with(unique(m2[chic.dir=="down",]), t.test(logFC.expr))
with(unique(m2[chic.dir=="up",]), t.test(logFC.expr))

with(unique(m2[chic.dir=="down" & Total_CD4_Activated>5 & Total_CD4_NonActivated>5,]), t.test(logFC.expr))
with(unique(m2[chic.dir=="up" & Total_CD4_Activated>5 & Total_CD4_NonActivated>5,]), t.test(logFC.expr))

bmean <- function(x,by,R=200) {
    f <- function(di,i){    
        d2.idx <- di[i,]
        d2 <- d[unlist(indices[d2.idx]),]
        mean(as.numeric(d2$x)) %>% as.numeric()
    }
    d <- cbind(x=x,by=by) %>% as.data.frame()
    d <- d[!is.na(x),]
    indices <- split(1:nrow(d),d$by)
    d.idx <- data.frame(i=1:length(indices))
    bb <- boot(d.idx, f, R=R) #, parallel="multicore",ncpus=5)
    s <- sd(bb$t)
    c(mean=bb$t0, se.mean=s, p=pnorm(abs(bb$t0),sd=s,lower.tail=FALSE)*2 )
}
with(m2,bmean(logFC.expr,id))
with(m2[logFC.chic<0,],bmean(logFC.expr,id))
with(m2[logFC.chic>0,],bmean(logFC.expr,id))

ggplot(m2[chic.dir!="inconsistent",],
       aes(x=logFC.chic,y=logFC.expr,col=chic.dir)) + geom_point() +
geom_smooth(method="lm")

summary(lmer( logFC.expr ~ logFC.chic + (1 | group), data=m2))
summary(lmer( logFC.expr ~ logFC.chic + (1 | group), data=m2[logFC.chic>0,]))
summary(lmer( logFC.expr ~ logFC.chic + (1 | group), data=m2[logFC.chic<0,]))

summary(lmer( logFC.expr ~ logFC.chic * sign + (sign | group), data=m2))

ggplot(m2,aes(x=logFC.chic,y=logFC.expr)) + geom_point() + geom_smooth(method="lm")
ggplot(m2[logFC.chic<0,],aes(x=logFC.chic,y=logFC.expr)) + geom_point() + geom_smooth(method="lm")
ggplot(m2[logFC.chic>0,],aes(x=logFC.chic,y=logFC.expr)) + geom_point() + geom_smooth(method="lm")


################################################################################

## use GenomicRanges to reduce
library(GenomicRanges)
library(parallel); options(mc.cores=20)

split.reduce <- function(x) {
  x.gr <- with(x,
               GRanges(seqnames=Rle(oeChr),
                       ranges=IRanges(start=oeStart,
                         width=oeLength)))
  mcols(x.gr) <- x[,c("baitID","id")]
  length.1 <- length(x.gr)
  
  x.gr <- split(x.gr,x.gr$id)
  x.gr <- mclapply(x.gr,reduce)
  x.gr <- GRangesList(x.gr)

  ## add ids as meta data
  G <- DataFrame(id=rep(names(x.gr),times=sapply(x.gr,length)))
  
  x.gr <- unlist(x.gr)
  mcols(x.gr) <- G
  length.2 <- length(x.gr)
  
  message("Reduction: ",length.1, " -> ", length.2)
  return(x.gr)
}

## next two lines are s.l.o.w....
f <- "act-non-gr.RData"
if(file.exists(f)) {
  load(f)
} else {
  act.gr <- split.reduce(subset(gene.ints, Total_CD4_Activated>5))
  non.gr <- split.reduce(subset(gene.ints, Total_CD4_NonActivated>5))
  colnames(mcols(act.gr)) <- "id"
  colnames(mcols(non.gr)) <- "id"
  save(act.gr, non.gr, file=f)
}


## NB little reduction!
## 369398 -> 337675 for act
## 373117 -> 339117 for non

## HERE

## collapse to one row per gene, counting interactions per gene
genes.act <- table(act.gr$id)
genes.non <- table(non.gr$id)
genes <- data.frame(row.names=unique(c(names(genes.act)))) #,names(genes.non))),act=NA,non=NA)
genes$id <- rownames(genes)
genes$act <- genes.act[rownames(genes)] %>% as.numeric()
genes$non <- genes.non[rownames(genes)] %>% as.numeric()
genes[is.na(genes)] <- 0
genes <- subset(genes,id!="." & !is.na(id))
dim(genes)

## add modules
## SHOULD DO THIS BY ENSEMBL ID
## genes <- merge(genes,modules[,c("geneSymbol","module")], by.x="gene", by.y="geneSymbol", all.x=TRUE)
## dim(genes)
## dim(modules)
## #mod.ints <- split(gene.ints, gene.ints$module)

## add rnaseq
genes <- merge(genes,
                   subset(rnaseq, id %in% gene.ints$id)[,c("id","act.total","non.total","logFC","adj.P.Val")],
                   all.x=TRUE, by="id")
dim(genes)

mycols <- unique(as.character(genes$module))
names(mycols) <- mycols
mycols[mycols=="NA"] <- "grey20"
colScale <- scale_colour_manual(name = "module",values = mycols, na.value="purple")

genes <- subset(genes, !is.na(adj.P.Val))
genes$int.fc <- (genes$act+1)/(genes$non+1)
with(genes, cor.test(log(int.fc), logFC,use="pair"))
with(genes, cor.test(act, act.total,use="pair"))
## data:  act and act.total
## t = 2.5593, df = 6968, p-value = 0.01051
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  0.007173003 0.054083232
## sample estimates:
##        cor 
## 0.03064499 

with(genes, cor.test(non, non.total,use="pair"))
## data:  non and non.total
## t = 2.0736, df = 6968, p-value = 0.03815
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  0.001357375 0.048282733
## sample estimates:
##        cor 
## 0.02483373 

png(file.path(paper.dir, "fig-exprFC-nintFC.png"),res=300,height=6,width=6)
ggplot(genes, aes(x=logFC, y=log2(int.fc))) + geom_point() + geom_smooth(se=FALSE,method="lm") +
    labs(x="expression fold change (log2 scale)", y="fold change in number of interacting fragments (log2 scale)") + colScale + theme(legend.position="none")
dev.off()
png(file.path(paper.dir, "fig-exprFC-nintFC-bymodule.png"),res=300,height=6,width=6)
ggplot(genes, aes(x=logFC, y=log2(int.fc))) + geom_point(aes(col=module)) + geom_smooth(se=FALSE,method="lm") +
    labs(x="expression fold change (log2 scale)", y="fold change in number of interacting fragments (log2 scale)") + facet_wrap(~module) + colScale + theme(legend.position="none")
dev.off()


################################################################################

## ** Does length of interacting fragments change with gene expression upon activation for OEs present in both?


## collapse to one row per gene, counting interaction fragment length per gene
summ <- function(x,y) {
  mean(x[y>5])
}
genes <- subset(gene.ints, id!="." & !is.na(id)) %>%
  group_by(., id) %>%
    summarise(., act=summ(oeLength,Total_CD4_Activated), non=summ(oeLength,Total_CD4_NonActivated))

## add modules
## SHOULD DO THIS BY ENSEMBL ID
genes <- merge(genes,modules[,c("ensGeneId","geneSymbol","module")], by.x="id", by.y="ensGeneId", all.x=TRUE)
dim(genes)
dim(modules)
#mod.ints <- split(gene.ints, gene.ints$module)

## add rnaseq

genes <- merge(genes,
                   subset(rnaseq, id %in% gene.ints$id)[,c("id","non.total","act.total","logFC","adj.P.Val")],
                   all.x=TRUE, by="id")
dim(genes)

mycols <- unique(as.character(genes$module))
names(mycols) <- mycols
mycols[mycols=="NA"] <- "grey20"
colScale <- scale_colour_manual(name = "module",values = mycols, na.value="purple")


genes <- subset(genes, !is.na(adj.P.Val))
genes$int.fc <- (genes$act+1)/(genes$non+1)
with(genes, cor.test(log(int.fc), logFC,use="pair"))

## 	Pearson's product-moment correlation

## data:  log(int.fc) and logFC
## t = 0.7207, df = 6441, p-value = 0.4711
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  -0.01544218  0.03339087
## sample estimates:
##         cor 
## 0.008979696 

with(genes, cor.test(act,act.total))
with(genes, cor.test(non,non.total))


png(file.path(CD4CHIC.OUTPUT,"paper", "fig-exprFC-intStrength.png"),res=300,height=6,width=6)
ggplot(genes, aes(x=logFC, y=log2(int.fc))) + geom_point() + geom_smooth(se=FALSE,method="lm") +
    labs(x="expression fold change (log2 scale)", y="fold change in average interaction score (log2 scale)") + colScale + theme(legend.position="none")
dev.off()
png(file.path(paper.dir, "fig-exprFC-intStrength-bymodule.png"),res=300,height=6,width=6)
ggplot(genes, aes(x=logFC, y=log2(int.fc))) + geom_point(aes(col=module)) + geom_smooth(se=FALSE,method="lm") +
    labs(x="expression fold change (log2 scale)", y="fold change in average interaction score (log2 scale)") + facet_wrap(~module) + colScale + theme(legend.position="none")
dev.off()

################################################################################


