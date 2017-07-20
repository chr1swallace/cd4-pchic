library(data.table)
library(magrittr)

source(file.path(CD4CHIC.ROOT,"activation-analyses/R/common.R"))
DO <- c()

## * bait-oe

header <- function(...)
    message("\n\t---\t",...,"\t---\n")

## interactions
header("Interactions")
ints <- get.interactions()
setkey(ints,baitID,oeID)
ints[,ngenesperbait:=length(unique(id)),by=baitID]
tmp <- subset(ints,baitID==223)[,ngenesperbait] %>% unique()
if(tmp!=3)
    stop("ngenesperbait should == 3 for baitID 223") 

## * pchic diff analysis
header("Differential Interactions")
diffchic <- get.diffchic()
setnames(diffchic,c("non","act"),c("chic.reads.non","chic.reads.act"))
setkey(diffchic,baitID,oeID)
ints <- merge(ints,diffchic)
setnames(ints,c("logFC","FDR"),c("logFC.chic","FDR.chic"))

## BAIT LINKED

## modules
header("modules")
spower <- "combined"
modules <- get.modules()

## RNA seq
header("RNAseq")
rnaseq <- get.rnaseq() # %>% subset(., !is.na(adj.P.Val))
subset(rnaseq, id %in% head(rnaseq$id[ duplicated(rnaseq$id) ]))
setnames(rnaseq,c("non.total","act.total"),c("expr.reads.non","expr.reads.act"))
setnames(rnaseq,c("logFC","adj.P.Val"), c("logFC.expr","FDR.expr"))

b2g <- get.b2gene()
setkey(b2g,id)
setkey(rnaseq,id)
setkey(modules,id)
expr <- merge(rnaseq,b2g[,.(id,gene,baitID)])
expr <- merge(expr,modules[,.(id,module)])
setkey(expr,baitID)
expr[,ngenesperbait:=length(unique(id)),by=baitID]

header("decay")
decay <- get.rnadecay()
setkey(expr,gene)
expr <- merge(expr,decay[,.(geneSymbol,HL.M,HL.328)],all.x=TRUE, by.x="gene",by.y="geneSymbol")


## FRAG LINKED

## SE/TE
header("TE/SE")
(load(file.path(CD4CHIC.OUT,"chipseq","hind-TE-SE.RData")))
tese <- as.data.frame(hind) %>% as.data.table()
tese$name %<>% as.integer()
setkey(tese,name)

## chromatin
header("Peaks and differential chipseq")

if("chromatin" %in% DO) {
    chromatin <- get.chromatin()
    chromatin <- merge(chromatin, tese[,.(name,TE,SE)])
}


## chromhmm
if("chromhmm" %in% DO) {
header("chromhmm")
f <- function(n) {
        (load(file.path(CD4CHIC.OUT,"chipseq",paste0("hind-",n,".RData"))))
        hind.dt <- mcols(hind) %>% as.data.table()
        hind.dt[,bait:=NULL]
        hind.dt$name %<>% as.integer()
        setkey(hind.dt,name)
        return(hind.dt)
    }
    hind.act <- f("Act_15")
    nm <- names(hind.act)[-1]
    setnames(hind.act,nm,paste0(nm,".act"))
    hind.non <- f("NAct_15")
    setnames(hind.non,nm,paste0(nm,".non"))
    chromhmm <- merge(hind.act, hind.non)
    chromatin <- merge(chromatin,chromhmm)
}
