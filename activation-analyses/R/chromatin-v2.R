library(ggplot2)

## * Chromatin peaks at baits/oe

source(file.path(CD4CHIC.ROOT,"activation-analyses/R/common.R"))

ints <- get.interactions()
setkey(ints,baitID,oeID)
ints[,ngenesperbait:=length(unique(id)),by=baitID]
tmp <- subset(ints,baitID==223)[,ngenesperbait] %>% unique()
if(tmp!=3)
    stop("ngenesperbait should == 3 for baitID 223") 

diffchic <- get.diffchic()
setnames(diffchic,c("non","act"),c("chic.reads.non","chic.reads.act"))
setkey(diffchic,baitID,oeID)
ints <- merge(ints,diffchic)
setnames(ints,c("logFC","FDR"),c("logFC.chic","FDR.chic"))
ints <- ints[biotype=="protein_coding",]
ints[,ngenesperbait:=length(unique(id)),by="baitID"]
ints <- ints[ngenesperbait==1,]    

######################################################################
## * data prep
    chromatin <- get.chromatin()
(load(file.path(CD4CHIC.OUT,"chipseq","hind-TE-SE.RData")))
tese <- as.data.frame(hind) %>% as.data.table()
tese$name %<>% as.integer()
setkey(tese,name)
setnames(tese,"seqnames","Chr")
chromatin <- merge(chromatin, tese[,.(name,TE,SE,bait,Chr,start,width)])

## nob2b <- ints[!(oeID %in% baitID) & !(baitID %in% oeID),]
## xx <- chromatin[name %in% c(nob2b$baitID,nob2b$oeID),]
## everything is either a prey or a bait, not both

cats <- list(
back.bait = ints[b2b==FALSE & pmax(Total_CD4_Activated,Total_CD4_NonActivated)<3,]$baitID,
back.oe = ints[b2b==FALSE & pmax(Total_CD4_Activated,Total_CD4_NonActivated)<3,]$oeID,
act.bait = ints[b2b==FALSE & Total_CD4_Activated>5 & Total_CD4_NonActivated<5 & FDR.chic<0.1,]$baitID,
act.oe = ints[b2b==FALSE & Total_CD4_Activated>5 & Total_CD4_NonActivated<5 & FDR.chic<0.1,]$oeID,
non.bait = ints[b2b==FALSE & Total_CD4_NonActivated>5 & Total_CD4_Activated<5 & FDR.chic<0.1,]$baitID,
non.oe = ints[b2b==FALSE & Total_CD4_NonActivated>5 & Total_CD4_Activated<5 & FDR.chic<0.1,]$oeID,
cd4.bait = ints[b2b==FALSE & Total_CD4_NonActivated>5 & Total_CD4_Activated>5 & FDR.chic>0.1,]$baitID,
cd4.oe = ints[b2b==FALSE & Total_CD4_NonActivated>5 & Total_CD4_Activated>5 & FDR.chic>0.1,]$oeID,
back.b2b = ints[b2b==TRUE & pmax(Total_CD4_Activated,Total_CD4_NonActivated)<3,]$baitID,
act.b2b = ints[b2b==TRUE & Total_CD4_Activated>5 & Total_CD4_NonActivated<5 & FDR.chic<0.1,]$baitID,
non.b2b = ints[b2b==TRUE & Total_CD4_NonActivated>5 & Total_CD4_Activated<5 & FDR.chic<0.1,]$baitID,
cd4.b2b = ints[b2b==TRUE & Total_CD4_NonActivated>5 & Total_CD4_Activated>5 & FDR.chic>0.1,]$baitID)

xx <- chromatin
xx[,cat:="--"]
for(nm in names(cats))
    xx[ name %in% cats[[nm]], cat:=nm]

## xx[ name %in% setdiff(c(back.bait,back.oe), c(act.bait,non.bait,act.oe,non.oe)), cell:="back"]
## xx[ name %in% setdiff(c(act.bait,act.oe), c(non.bait,non.oe)), cell:="act"]
## xx[ name %in% setdiff(c(non.bait,non.oe),c(act.bait,act.oe)), cell:="non"]
## xx[ name %in% c(intersect(act.bait,non.bait),intersect(act.oe,non.oe)), cell:="cd4"]
## xx[ name %in% setdiff(c(back.bait,act.bait,non.bait),c(back.oe,act.oe,non.oe)), class:="bait"]
## xx[ name %in% setdiff(c(back.oe,act.oe,non.oe),c(back.bait,act.bait,non.bait)), class:="oe"]
## xx[ name %in% c(back.b2b,act.b2b,non.b2b), class:="b2b"]
## xx <- xx[class!="--" & cell!="--",]
## xx[,cat:=paste(class,cell,sep="/")]

## add fragment length for stratification in testing
## * Check for confounding by fragment width etc (can skip once done)

## ** baited regions are longer than oe regions are longer than background regions

xx$cat2 <- "--"
xx[bait=="bait" & grepl("act|cd4|non",cat),cat2:="interacting bait"]
xx[bait=="oe" & grepl("act|cd4|non",cat),cat2:="interacting PIR"]
xx[bait=="bait" & grepl("--|back",cat),cat2:="non-interacting bait"]
xx[bait=="oe" & grepl("--|back",cat),cat2:="non-interacting, not bait"]
xx[,cat2:=factor(cat2, levels=rev(c("non-interacting, not bait","interacting PIR",
                                    "non-interacting bait","interacting bait")))]

## * Plot width of hind III by interaction category

library(ggplot2)
library(cowplot)
ggplot(xx,aes(x=width/1000)) +
geom_histogram(binwidth=0.1,col="grey") + #aes(fill=baitprey),col="grey",) +
#  geom_density(aes(y=..count..)) +
facet_grid(cat2 ~ .,scales="free") +
scale_x_log10("HindIII Fragment width (kb)",
              breaks=c(0.1,1,10,1000,10000),
              labels=c("0.1","1","10","1,000","10,000")) +
ylab("Fragments") +
background_grid() + 
theme(legend.position="bottom")

ggsave(file.path(CD4CHIC.OUT,"paper/hind-width-by-cat.pdf"),height=8,width=8)
ggsave(file.path(CD4CHIC.OUT,"paper/hind-width-by-cat.tiff"),height=8,width=8)
write.table(xx[,.(name,Chr,start,width,cat)],
            file=file.path(CD4CHIC.OUT,"paper/hind-width-by-cat.csv"))


################################################################################

## downsample cat="--"
table(xx$cat)
set.seed(7654321)
drop <- sample(which(xx$cat=="--"),sum(xx$cat=="--") - 50000)
xx <- xx[-drop,]
table(xx$cat)

## library(devtools)
## install_github("chr1swallace/genomic.autocorr")

library(genomic.autocorr)
options(mc.cores=20)

quantile(xx$width,seq(0,1,0.1))
quantile(xx$width,seq(0,1,0.25))
S <- stratum<-c(0,1000,2000,5000,10000,1e7)

const <- c("name","cat","cat2","Chr","start","width","bait")
    xxm <- melt(xx, const)
    setnames(xxm,c("variable","value"),c("mark","score"))    
    xxm$strat <- cut(xxm$width,stratum)
    xxm[,mark:=sub("_NAct$|_Act$","",mark)]
    xxm[,score:=ifelse(is.na(score) | is.nan(score), 0, score)]
    xxm[,dscore:=as.integer(score>0)]    

glm(dscore ~ cat + strat, data=xxm[mark=="H3K27ac",], family="binomial") %>% summary()


## * There is autocorrelation

summ <- acf.summary(xx,variables=setdiff(names(xx),const),order.by="name")
ggplot(melt(summ,c("lag","variable"),variable.name="ac.type"),
       aes(x=lag,y=value,col=variable)) + geom_path() + facet_grid(ac.type~.)

##' looks ok by lag=20

## * Plot/test diffchrom

xx.bak <- xx

xxd <- xx[,c("name","cat","Chr","start","width",grep("_up$|_down$",names(xx),value=TRUE)),with=FALSE]
xxa <- xx[,c("name","cat","Chr","start","width",grep("_Act$|TE|SE",names(xx),value=TRUE)),with=FALSE]
xxn <- xx[,c("name","cat","Chr","start","width",grep("_NAct$|TE|SE",names(xx),value=TRUE)),with=FALSE]
xxd[,cell:="diff"]
xxa[,cell:="Act"]
xxn[,cell:="NonAct"]
setnames(xxa,sub("_Act","",names(xxa)))
setnames(xxn,sub("_NAct","",names(xxn)))
xx <- rbindlist(list(xxa,xxn))

## * Test

library(cowplot)
plotter <- function(rd) {
    rd <- rd[grep("cat",x),]
    rd[,int:=sub(".*\\.","",sub("oe","PIR",x))]
    rd[,cell:=factor(sub("\\..*.","",sub("cat","",x)),levels=c("back","cd4","non","act"))]
    rd[,cat:=factor(sub("cat","",x),
                    levels=c("back.b2b","cd4.b2b","non.b2b","act.b2b",
                             "back.bait","cd4.bait","non.bait","act.bait",
                             "back.oe","cd4.oe","non.oe","act.oe"))]
    levels(rd$cat) <- c("B2B back","B2B CD4", "B2B CD4(n)","B2B CD4(a)",
                        "bait-oe-cmp.R back","bait CD4", "bait CD4(n)","bait CD4(a)",
                        "PIR back","PIR CD4", "PIR CD4(n)","PIR CD4(a)")
    ggplot(rd,
           aes(x=cell,y=beta,ymin=beta.025,ymax=beta.975)) +
    geom_hline(aes(yintercept=beta),lty="dotted",col="grey",data=rd[cell=="back",]) +
    geom_pointrange() +
    background_grid() +
    panel_border() +
    facet_grid(y ~ int, scales="free") +
    labs(x="Interaction category",y="Log OR of mark overlap")
}

xx[,width2:=width^2]
xx[,width3:=width^3]

options(mc.cores=25)
r <- block.glm(f.lhs=c("H3K27ac","H3K4me3","H3K4me1","TE","SE"),
          f.rhs="width + width2 + width3 + cell + cat",
          data=xx,
          order.by="name",
          strat.by="Chr",
          block.size=40,
          B=50,
          family="binomial")
   
          
xxd[,width2:=width^2]
xxd[,width3:=width^3]

options(mc.cores=25)
rd <- block.glm(f.lhs=c("H3K27ac_up","H3K4me3_up","H3K4me1_up","H3K27ac_down","H3K4me3_down","H3K4me1_down"),
          f.rhs="width + width2 + width3 + cat",
          data=xxd,
          order.by="name",
          strat.by="Chr",
          block.size=40,
          B=50,
          family="binomial")


plotter(r)
ggsave(file=file.path(CD4CHIC.OUT,"paper/figure-chrom-baitprey.pdf"),width=9,height=12)

plotter(rd)
ggsave(file=file.path(CD4CHIC.OUT,"paper/figure-diffchrom-baitprey.pdf"),width=9,height=12)

p1 <- plotter(r)
p2 <- plotter(rd)
plot_grid(p1,p2,nrow=1,labels=c("a","b"))
ggsave(file=file.path(CD4CHIC.OUT,"paper/figure-chrom-diffchrom-baitprey.pdf"),width=15,height=10)
    
