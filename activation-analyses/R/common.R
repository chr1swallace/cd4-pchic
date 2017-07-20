library(magrittr)

## global variables

## where all output is written
plot.dir <- out.dir <- file.path(CD4CHIC.OUT,"modules")
paper.dir <- file.path(CD4CHIC.OUT,"paper")
## file prefix
file.prefix<-'t-cell-act'
##variable that defines whether a gene is expressed - genes that are below this value in all samples are excluded
expression.threshold<-6
## network.type - whether we want a biological pathway where modules contain things that are both correlated and anti correlated or just correlated.
network.type <- 'signed' ## argument is that this better abstracts underlying biological pathway

## expression matrix and associated data
file.expression <- file.path(out.dir,"expression-data.RData")

make.file <- function(nm,dir,suffix) {
    spower <- "combined"
    f <- file.path(dir,
                   paste0(file.prefix,"_",network.type,"_power_",spower,"_",nm,suffix))
    message("file: ",f)
    if(any(file.exists(f)))
    warning(f," will be overwritten")
    return(f)
  
}
plot.file <- function(nm,suffix=".pdf") make.file(nm,dir=plot.dir,suffix=suffix)
out.file <- function(nm,suffix=".RData") make.file(nm,dir=out.dir,suffix=suffix)

## handy data loaders
pmean <- function(...) {
  x <- do.call("cbind",list(...))
  rowSums(x)/ncol(x)
}
get.rnaseq <- function() {
    rnaseq <- fread(file.path(CD4CHIC.DATA,"rna-with-regrna-diff-v2.csv"))
    #rnaseq<-fread(file.path(CD4CHIC.DATA,"rna-kallisto-gene-diff.csv"),header=TRUE,stringsAsFactors=FALSE,sep="\t")
    rnaseq <- rnaseq[!is.na(adj.P.Val),]
    setnames(rnaseq, c("1_non","1_act","2_non","2_act"), c("unstm.1","activ.1","unstm.2","activ.2"))
    #names(rnaseq),sub("gene","id",make.names(names(rnaseq))))

    ## restrict to protein coding
    ## if(!exists("b2g"))
    ##     b2g <- get.bait2gene()
    ## pg <- b2g[biotype=="protein_coding",]$id %>% unique()
    ## rnaseq <- rnaseq[ id %in% pg,]
    rnaseq[, c("act.total","non.total") := list(pmean(activ.1,activ.2), pmean(unstm.1,unstm.2)) ]
    return(rnaseq)
}

get.modules <- function() {
   geneSummary.file<-file.path(CD4CHIC.OUT,"modules/t-cell-act_signed_power_combined_geneSummary.RData")
    (load(geneSummary.file))
    if("colors" %in% colnames(geneInfo)) {
        geneInfo$module <- geneInfo$colors
    } else {
        geneInfo$module <- geneInfo$all_12
    }
    return(geneInfo[ , -grep("^all_|^treated_|^colors", colnames(geneInfo)), with=FALSE])
}
get.baitints <- function(baits) {
## files <- file.path("..", paste(args$id,".gz",sep=""))
myread <- function(bait,actnon) {
  f <- file.path(CD4CHIC.DATA, paste0(actnon,"-bybait"), paste0(bait,".gz"))
  if(!file.exists(f))
    stop("file not found: ",f)
  xx <- read.table(f,as.is=TRUE,sep="\t",header=FALSE)
  names(xx) <- c("numPairs",	"nTrans",	"otherEndID",	"otherEndChr",	"baitID",	"baitChr",	"s_j",	"otherEndLen",	"distSign",	"isBait2bait",	"N.1",	"N.2",	"N.3",	"N",	"refBinMean",	"NNb",	"distbin",	"s_i",	"NNboe",	"transBaitLen",	"tlb",	"tblb",	"Tmean",	"Bmean",	"log.p",	"log.w",	"log.q",	"score",	"newScore")
  xx$actnon <- actnon
  return(xx)
}

non <- lapply(baits,myread,"non")
act <- lapply(baits,myread,"act")
data <- do.call("rbind",c(non,act))
rm(non,act)

data$cscore <- cut(data$newScore,c(0,3,5,100),include.lowest=TRUE)
data$isBait2bait <- ifelse(data$isBait2bait,"bait2bait","bait2prey")
return(data)
}

get.b2gene <- function() {
    require(data.table)
    int.genes <- fread(file.path(CD4CHIC.DATA,"merged_samples_12Apr2015_full_denorm_bait2baits_e75.tab"),select=c(1:8))
#    nm <- scan(file.path(CD4CHIC.DATA,"merged_samples_12Apr2015_full_denorm_bait2baits.txt"), nlines=1, what="")
    setnames(int.genes, names(int.genes), c("id","gene","biotype","strand","baitChr","baitStart","baitEnd","baitID"))
    setkey(int.genes,baitID,id)
    int.genes <- unique( int.genes )
    return(int.genes)
}

get.interactions <- function(threshold=5) {
    require(data.table)#merged_samples_12Apr2015_full_denorm_bait2baits_e75.tab
##    int <- fread(file.path(CD4CHIC.DATA,"merged_samples_12Apr2015_full_denorm_bait2baits_e75.tab")) # %>% as.data.frame()
    int <- fread(file.path(CD4CHIC.DATA,"merged_samples_12Apr2015_full.txt"))
    int.genes <- get.b2gene()
    setkey(int,baitID)
    setkey(int.genes,baitID)
    int <- merge(int,int.genes[,.(id,gene,biotype,strand,baitID)],allow.cartesian=TRUE)
    ## ## subset to protein coding only
    ## int <- int[biotype=="protein_coding",]
    ## int[ oeID==120852 & baitID==120843, .(baitChr,baitStart,baitID,oeChr,oeStart,oeID,id,gene,Total_CD4_Activated,Total_CD4_NonActivated)]
    int <- int[,b2b:= oeID %in% int$baitID]
    int[,Background:=apply(int[,.(Erythroblasts,Megakaryocytes)],1,max)]
    if(!is.null(threshold))
        int <- int[pmax(Background,Total_CD4_Activated,Total_CD4_NonActivated)>=threshold, ]
    int <- int[, c("baitChr","oeChr") := list(paste0("chr",baitChr), paste0("chr",oeChr)) ]
    int <- int[, c("baitLength","oeLength"):= list(abs(baitStart-baitEnd), abs(oeStart-oeEnd)) ]
    int <- int[,.(baitChr,baitStart,baitID,oeChr,oeStart,oeID,id,gene,Total_CD4_Activated,Total_CD4_NonActivated,Background,biotype,baitLength,oeLength,b2b)]
    return(unique(int,by=NULL)) # sometimes one bait maps to multiple promoters for the same gene
}

get.diffchic <- function() {
    require(data.table)
    tmp <- fread(file.path(CD4CHIC.DATA, "pchic-diff.csv"))
    setnames(tmp,colnames(tmp),make.names(colnames(tmp)))
    tmp <- tmp[, non:=P1.non+P2.non+P3.non]
    tmp <- tmp[, act:=P1.act+P2.act+P3.act]
    tmp <- tmp[, total:=non+act]
    ss <- strsplit(tmp[,oe.bait],":")
    tmp[, oeID:=as.integer(sub(":.*","",oe.bait))]
    tmp[, baitID:=as.integer(sub(".*:","",oe.bait))]
    tmp
}

get.rnadecay <- function() {
    xx <- fread(file.path(CD4CHIC.DATA,"Bohjanen-T-cell.csv"),skip=15,header=TRUE,
                select=c(1:4,11,13:14,17,21,23:24,27,31,33:34,36))
    names(xx) <- c("PROBE_ID", "call.M", "call.3", "call.328", "HL.M", "P.M.3",
                   "P.M.328", "FC.3.M", "HL.3", "P.3.M", "P.3.328",
                   "FC.328.M", "HL.328", "P.328.M", "P.328.3",
                   "geneSymbol")
    mods <- get.modules()
    setkey(mods,geneSymbol)
    imedian <- function(x) { median(x) %>% as.integer() }
    xx[, c("MHL.3","MHL.M","MHL.328"):=
       list(imedian(HL.3),imedian(HL.M),imedian(HL.328)), by=geneSymbol]
    setkey(xx,geneSymbol)
    xx <- unique(xx)
    
    merge(xx, unique(mods[,.(geneSymbol,id)]))
}

get.hind <- function() {
    require(rtracklayer)
    f <- file.path(CD4CHIC.DATA, "Digest_Human_HindIII.bed")
    hind <- import(f)
    seqlevels(hind) <- paste0("chr",seqlevels(hind))
    hind <- keepStandardChromosomes(hind)
#    int0 <- get.interactions(threshold=0)
    b2g <- get.b2gene()
    baits <- unique(b2g$baitID)
    
    hind$bait <- ifelse(hind$name %in% baits, "bait", "oe")
    ## hind$cell <- ifelse(hind$name %in% subset(int,Total_CD4_Activated>5 & Total_CD4_NonActivated>5,select="oeID",drop=TRUE), "both",
    ##              ifelse(hind$name %in% subset(int,Total_CD4_Activated>5,select="oeID",drop=TRUE),
    ##                     "act","non"))
    return(hind)
}

make.colScale <- function(mycols) {
require(ggplot2)
mycols <- as.character(mycols)
names(mycols) <- mycols
mycols[mycols=="NA"] <- "grey20"
mycols["yellow"] <- "goldenrod"
mycols["brown"] <- "SaddleBrown"
mycols["grey"] <- "grey40"
mycols["green"] <- "DarkGreen"
mycols["blue"] <- "DarkBlue"
scale_colour_manual(name = "module",values = mycols, na.value="purple",guide=FALSE)
}


Cochran.Armitage.test <- function(exposure, cc, stratum=rep(1,length(cc)),quick=FALSE) {
  N <- length(cc)
  if (is.factor(exposure))
    exposure <- as.numeric(exposure)
  cl <- match.call()
  narg <- length(cl) - 1
  arguments <- character(narg-1)
#   for (i in 1:narg)
#     arguments[i] <- as.character(cl[[i+1]])
  use <- complete.cases(cc,exposure,stratum)
  if (any(!use)) {
    cc <- cc[use]
    exposure <- exposure[use]
    stratum <- stratum[use]
  }
  dh <- table(stratum, cc)
  if (ncol(dh)!=2)
    stop("cc argument must have two levels")
  nt <- table(stratum)
  zm <- tapply(exposure, list(stratum, cc), mean)
  ut <- dh[,1]*dh[,2]*(zm[,2] - zm[,1])/nt
  if(quick)
  	return(sum(ut,na.rm=T))
  zv <- tapply(exposure, stratum, var)
  vt <-  dh[,1]*dh[,2]*zv/nt
  x2 <- sum(ut,na.rm=T)^2/sum(vt,na.rm=T)
  names(x2) <- "Chi-squared"
  df <- 1
  names(df) <- "df"
  res <- list(statistic=x2, parameter=df,
              p.value=pchisq(x2, 1, lower.tail=FALSE),
              method="Cochran-Armitage test with Mantel's extension",
              data.name=paste(cl,collapse=" "),
              score=ut,
              sum.score=sum(ut,na.rm=T),
              score.variance=vt,
              sum.variance=sum(vt,na.rm=T))
  class(res) <- "htest"
  res
}
##' Sanitizer function for xtable that converts scientific notation to latex style
##'
##' @param x character or numeric vector 
##' @return sanitized character
##' @export
##' @author Chris Wallace
##' @examples
##' p<-10^(-seq(1,10))
##' format.pval(p) # scientific notation
##' mysan(p) # latex style
##'
##' library(xtable)
##' df <- data.frame(numeric=p,scientific=format.pval(p),p.san=mysan(p))
##' print(xtable(df)) # doesn't work
##' print(xtable(df), sanitize.text.function=mysan)  # looks nice!
mysan <- function(x) {
  if(is.numeric(x))
    x <- base::format.pval(x,digits=3)
  sub("(<? ?[0-9\\.]+)e-([0-9]+)","$\\1\\\\times10^{-\\2}$",x)
}

get.translength <- function() {
  trans <- fread(file.path(CD4CHIC.DATA,"transcript-length.csv"),header=TRUE)
  gene <- fread(file.path(CD4CHIC.DATA,"max-transcript-length.csv"),header=TRUE)
  setnames(gene,c("id","gene.length","name"))
  link <- fread(file.path(CD4CHIC.DATA,"gene-transcript-link.csv"),header=TRUE)
  trans <- merge(trans,link,by="transcript")
  trans[,max.trans.length:=max(bp),by="gene"]
  gene <- merge(gene, unique(trans[,.(gene,max.trans.length)],by="gene"), by.x="id",by.y="gene")
  gene[,.(id,max.trans.length)]
}

library(GenomicRanges)
add.hind <- function(erna,start="start",end="end",id="id") {
    gr <- GRanges(seqnames=Rle(paste0("chr",erna[["chr"]])),
                       ranges=IRanges(start=erna[[start]],end=erna[[end]]))
    mcols(gr) <- erna[,id,with=FALSE]
    gr <- mergeByOverlaps(gr, hind)
    gr <- as.data.table(as.data.frame(gr[,c(id,"name","bait")]))
    gr[,nfrag:=.N,by="id"]
    merge(erna,gr,by=c("id"))
}
## get.erna <-  function() {
##  erna <- fread(file.path(CD4CHIC.DATA,"rna-with-regrna-diff-v2.csv"))
##     setnames(erna,c("logFC","adj.P.Val","1_act","2_act","1_non","2_non"),
##              c("logFC.erna","FDR.erna","act_1","act_2","non_1","non_2"))
##  message("note column name changes:")
##  message("logFC -> logFC.erna")
##  message("adj.P.Val -> FDR.erna")
##  erna[,expr:= (as.numeric(non.cpm1>=0.4) + as.numeric(non.cpm2>=0.4) +
##                as.numeric(act.cpm1>=0.4) + as.numeric(act.cpm2>=0.4)) > 1 ]
##  erna[,non.expr:=pmin(non.cpm1,non.cpm2)>=0.4]
##  erna[,act.expr:=pmin(act.cpm1,act.cpm2)>=0.4]
##  erna <- erna[type %in% c("regulatory","lincRNA","protein_coding","pseudogene"), ]
##     intergenic.ids <- scan(file.path(CD4CHIC.DATA, "distant-erna.csv"),what="")
##  erna[,intergenic:=id %in% intergenic.ids]
## erna[id %in% intergenic.ids,type:="intergenic.reg"]
## message("expressed in non-activated")
## with(erna,table(non.expr,type)) # 3897, want 3897 - agree
## message("expressed overall")
## with(erna,table(non.expr | act.expr,type)) # 5898, want 6133
## with(erna,table(expr,type)) # 6147, want 6133

##  ## add hindIII ids
##   if(!exists("hind") || !("GRanges" %in% class(hind)))
##     hind <- get.hind()

##   erna <- add.hind(erna)
##  erna[,name:=as.integer(name)]
##  return(erna)
## }

## get.erna <- function() {
##   library(GenomicRanges)
##   erna.type <- fread(file.path(CD4CHIC.DATA,"rna-with-regrna-diff-v3.csv"))
##   erna <- fread(file.path(CD4CHIC.DATA,"rna-with-regrna-diff.csv"))
##   m <- merge(erna,erna.type[,.(id,type)],by="id",all.x=TRUE)
##   erna <- erna[type=="regulatory" & !is.na(adj.P.Val),]

##   ## add hindIII ids
##   if(!exists("hind"))
##     hind <- get.hind()

##   gr <- with(erna,
##              GRanges(seqnames=Rle(paste0("chr",chr)),
##                      ranges=IRanges(start=start,end=end)))
##   mcols(gr) <- erna[,.(id)]
##   gr <- mergeByOverlaps(gr, hind)
##   gr <- as.data.table(as.data.frame(gr[,c("id","name","bait")]))
##   gr[,nfrag:=.N,by="id"]

##   merge(erna,gr,by=c("id"))
## }

get.erna <-  function() {
    erna <- fread(file.path(CD4CHIC.DATA,"rna-with-regrna-diff-v2.csv"))
    setnames(erna,c("logFC","adj.P.Val","1_act","2_act","1_non","2_non"),
             c("logFC.erna","FDR.erna","act_1","act_2","non_1","non_2"))
    message("note column name changes:")
    message("logFC -> logFC.erna")
    message("adj.P.Val -> FDR.erna")
    erna[,non.cpm1:=1e6 * non_1/sum(non_1)]
    erna[,non.cpm2:=1e6 * non_2/sum(non_2)]
    erna[,act.cpm1:=1e6 * act_1/sum(act_1)]
    erna[,act.cpm2:=1e6 * act_2/sum(act_2)]
    erna[,expr:= (as.numeric(non.cpm1>=0.4) + as.numeric(non.cpm2>=0.4) +
                  as.numeric(act.cpm1>=0.4) + as.numeric(act.cpm2>=0.4)) > 1 ]
    erna[,non.expr:=pmin(non.cpm1,non.cpm2)>=0.4]
    erna[,act.expr:=pmin(act.cpm1,act.cpm2)>=0.4]
    erna <- erna[type %in% c("regulatory","lincRNA","protein_coding","pseudogene"), ]
    intergenic.ids <- scan(file.path(CD4CHIC.DATA, "distant-erna.csv"),what="")
    erna[,intergenic:=id %in% intergenic.ids]
    erna[id %in% intergenic.ids,type:="intergenic.reg"]
    message("expressed in non-activated")
    with(erna,table(non.expr,type)) # 3897, want 3897 - agree
    message("expressed overall")
    with(erna,table(non.expr | act.expr,type)) # 5898, want 6133
    with(erna,table(expr,type)) # 6147, want 6133
    
    ## add hindIII ids
    if(!exists("hind") || !("GRanges" %in% class(hind)))
        hind <- get.hind()
    
    erna <- add.hind(erna)
    erna[,name:=as.integer(name)]
    return(erna)
}


## misc functions
bootcor <- function(x,y,R=200) {
    f <- function(d,i){    
        d2 <- d[i,]
        cor(d2$x, d2$y)
    }
    d <- cbind(x=x,y=y) %>% as.data.frame()
    bb <- boot(d, f, R=R) #, parallel="multicore",ncpus=5)
    s <- sd(bb$t)
    c(rho=bb$t0, se.rho=s, p=pnorm(abs(bb$t0),sd=s,lower.tail=FALSE)*2 )
}
plotcor <- function(x,y,xmax=NULL,y.label="bottom") {
    stats <- bootcor(x,y)
    label <- paste("rho == ", round(stats[["rho"]],4))
    if(!is.null(xmax)) {
        n <- sum(x>xmax)
        message("dropping ",n,"/",length(x)," observations above ",xmax)
        message("full range of data is: ",min(x)," - ",max(x))
        y <- y[x<xmax]
        x <- x[x<xmax]
    }
    label.y <- if(y.label=="bottom") { min(y) } else { 0.95*max(y) }
    p <- qplot(x,y) + geom_point() + geom_smooth() + annotate("text",x=max(x),y=label.y,label=label,parse=TRUE,hjust=1,vjust=0)
    return(list(plot=p,stats=stats))
}


get.chromatin <- function() {
    require(GenomicRanges)
    require(S4Vectors)
    (load(file.path(CD4CHIC.OUT,"chipseq","hind-our-peaks.RData")))
    ## if(!exists("hind"))
    ##     hind <- get.hind()
    hind <- mcols(hind)[,-2] %>% as.data.frame() %>% as.data.table()
hind$name %<>% as.integer()
setkey(hind,name)
hind

## (load(file.path(CD4CHIC.OUT,"chipseq","hind-our-marks.RData")))
## marks <- as.data.table(as.data.frame(mcols(hind))[,-c(2,3)])
## marks[,name:=as.integer(name)]
## keep <- "name"
## names.orig <- names(marks)
## for(mm in c("H3K27ac","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K9me3")) {
##     for(cell in c("_Act","_NAct")) {
##         nm <- intersect(grep(mm,names.orig,value=TRUE),
##                         grep(paste0(cell,"_"),names.orig,value=TRUE))
##         nm2 <- paste0(mm,cell)
##         if(length(nm)>1) {
##             marks[,c(nm2):=list(rowMeans(marks[,nm,with=FALSE]))]
##         } else {
##             setnames(marks,nm,nm2)
##         }
##         keep <- c(keep,nm2)        
##     }
## }
## marks <- marks[,keep,with=FALSE]
## head(marks)
## merge(chromatin,marks)
}
