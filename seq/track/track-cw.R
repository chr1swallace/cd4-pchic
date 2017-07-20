library(Gviz)
library(GenomicInteractions)
setwd(file.path(CD4CHIC.ROOT,"seq/track"))
source("stranded.r")
options(ucscChromosomeNames=FALSE)
hindt <- function(chromosome, from, to) {
    data <- hind[seqnames == chromosome & end > from & start < to,]
    an <- AnnotationTrack(chromosome=data$seqnames, start=data$start, end=data$start, name="HindIII", stacking="dense", stackHeight=0.15)
    displayPars(an)$min.width <- 0.1
    an
    ## AnnotationTrack(chromosome=data$seqnames, start=data$start, end=data$start, name="HindIII", group="H",col="black",
    ##                 groupAnnotation="group", genome="ENSEMBL", stacking="dense", stackHeight=0.15)
}
hict <- function(rdat, name, chromosome, from, to, i.baitID, col) {
    data <- rdat[rdat$baitID %in% i.baitID,] # &
                                        #                 rdat$oeChr == chromosome & rdat$oeEnd >= from & rdat$oeStart <= to,]
    if(!nrow(data))
        return(NULL)
    bait <- with(data, GRanges(baitChr, IRanges(start=baitStart, width=baitLength)))
    oe <- with(data, GRanges(oeChr, IRanges(start=oeStart, width=oeLength)))
    interaction <- GenomicInteractions(bait, oe, counts=1)
    hic <- InteractionTrack(name=name, interaction, chromosome=chromosome)
    displayPars(hic)$anchor.height <- 0.2
    displayPars(hic)$col.anchors.line <- "darkgrey"
    displayPars(hic)$col.anchors.fill <- "lightgrey"
    displayPars(hic)$col.interactions <- col
    hic
}
## hict <- function(file, name, chromosome, from, to, col) {
##     rdat <- read.table(file, header=T)
##     data <- rdat[rdat$baitChr == chromosome & rdat$baitStart >= from & rdat$baitEnd <= to &
##                  rdat$oeChr == chromosome & rdat$oeStart >= from & rdat$oeEnd <= to,]
##     bait <- GRanges(data$baitChr, IRanges(start=data$baitStart, end=data$baitEnd))
##     oe <- GRanges(data$oeChr, IRanges(start=data$oeStart, end=data$oeEnd))
##     interaction <- GenomicInteractions(bait, oe, counts=1)
##     hic <- InteractionTrack(name=name, interaction, chromosome=chromosome)
##     displayPars(hic)$anchor.height <- 0.2
##     displayPars(hic)$col.anchors.line <- "lightgrey"
##     displayPars(hic)$col.anchors.fill <- "orange"
##     displayPars(hic)$col.interactions <- col
##     hic
## }
rnat <- function(file, name, ylim) {
    logt <- function(x) sign(x) * log2(abs(x)+0.1)
#    kk <- strandedBamImport("data_srt.2_non.bam")
    d <- DataTrack(range=file, importFunction=strandedBamImport, stream=T, baseline=0, col.baseline="grey",
              col=c("#e31a1c","#1f78b4"), groups=c("Forward", "Reverse"), type=c("histogram", "g"), col.histogram=NA,
              name=name, ylim=ylim, transformation=logt)
    displayPars(d)$col.grid <- adjustcolor("grey", alpha.f = 0.2)
    displayPars(d)$v <- 0
    d
}
## rnat <- function(file, name, ylim) {
##     logt <- function(x) sign(x) * log2(abs(x)+0.1)
##     DataTrack(range=file, importFunction=strandedBamImport, stream=T, baseline=0, col.baseline="grey",
##               col=c("#e31a1c","#1f78b4"), groups=c("Forward", "Reverse"), type="histogram", col.histogram=NA,
##               name=name, ylim=ylim, transformation=logt)
## }
chpt <- function(file, name, ylim, col) {
                                        #    AlignmentsTrack(file, isPaired=F, type=c("mountain"), name=name, ylim=c(0, ylim))
    d <-
    ## AlignmentsTrack(file, isPaired=F, type=c("coverage"), name=name, ylim=c(0, ylim),
    DataTrack(file, isPaired=F, type=c("polygon"), window=-1, windowSize=1000, baseline=0, name=name, ylim=c(0, ylim)
             ## ,transformation=function(x){x[x<10] <- NA; x}
              )
    displayPars(d)$col=col
    d
}
regt <- function(file, name, chromosome) {
    d <- read.table(file)
    colnames(d) <- c("chromosome", "start", "end", "state")
    d$start <- d$start+1
    s <- d[d$chromosome==chromosome,]
    txn <- NA
    enhancer <- "#cab2d6"
    promoter <- "#cab2d6"
    rep <- NA
    AnnotationTrack(chromosome=s$chromosome, start=s$start, end=s$end, feature=s$state,
                    name=name, genome="ENSEMBL", stacking="dense", stackHeight=0.15,
                    E1=txn, E2=txn, E3=txn, E4=enhancer, E5=enhancer,
                    E6=enhancer, E7=promoter, E8=promoter, E9=promoter, E10=promoter,
                    E11=enhancer, E12=rep, E13=rep, E14=rep, E15=rep,
                    col=NA)
}
rgnt <- function(file, name, chromosome) {
    d <- read.table(file)
    colnames(d) <- c("chromosome", "start", "end")
    d$start <- d$start+1
    s <- d[d$chromosome==chromosome,]
    enhancer <- "#cab2d6"
    AnnotationTrack(chromosome=s$chromosome, start=s$start, end=s$end,
                    name=name, genome="ENSEMBL", stacking="dense", stackHeight=0.15)
}
ant <- function(range) AnnotationTrack(range=range, genome="ENSEMBL", name="GWAS", stackHeight=0.3, col=NA)
plotregion <- function(name, chromosome, from, to, gwast, hicn, hica, rlim, clim) {
#    pdf(file.path("~/",paste(name, "pdf", sep=".")), paper="a4")
#    dev.off()
}
plotter <- function(nm,plus=NULL) {
    while(length(w <- which(search()=="params[[nm]]"))) 
        detach(pos=w[1])
    attach(params[[nm]])
    tracks <- TRACKS[[nm]]
    trt.l <- with(tracks, list(hd3, tracks$hicn, tracks$hica))
    trt.r <- with(tracks, list(gwast, genes, rnan, rnaa, ac27n, ac27a, 
                               reg))
    if(!is.null(plus) & length(plus)) {
        for(p in plus) {
            a <- TRACKS[[p]]$hicn
            b <- TRACKS[[p]]$hica
            trt.l <- c(trt.l, list(a,b))
        }
    }
    drop <- sapply(trt.l,is.null)
    if(any(drop))
        trt.l <- trt.l[!drop]
    sz.l <- c(rep(0.3,1),rep(0.5,length(trt.l)-1))
    sz.r <- c(rep(0.3,2),rep(0.5,length(trt.r)-3),0.3)
    trt <- c(trt.l, trt.r)
    sz <- c(0.3,sz.l,sz.r)
    ht <- HighlightTrack(trackList=trt,
                         chromosome=chromosome, range=range(tracks$gwast), col="grey", alpha=0.25)
    plotTracks(c(axis, ht),
               chromosome=chromosome, from=from, to=to,
               sizes=sz,
               background.title=NA, col.axis="darkgrey", fontcolor.title="black", col.axis="black")
    detach(params[[nm]])
}

## common params
setwd("/home/ar756/plots")
rlim=c(-12, 12)
clim=220

library(biomaRt)
e75.genemart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",  host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

library(data.table)
pmi <- list(CRO=fread("/home/oliver/JAVIERRE_ICHIP/out/pmi/CRO_IC.pmi"),
            T1D=fread("/home/oliver/JAVIERRE_GWAS/out/pmi/T1D.pmi"),
            RA=fread("/home/oliver/JAVIERRE_GWAS/out/pmi/RA_OKADA_IMB.pmi"),
            MS=fread("/home/oliver/JAVIERRE_ICHIP/out/pmi/MS_IC.pmi"),
            CEL=fread("/home/oliver/JAVIERRE_ICHIP/out/pmi/CEL_IC.pmi"),
            PBC=fread("/home/oliver/JAVIERRE_ICHIP/out/pmi/PBC_IC.pmi"),
            SLE=fread("/home/oliver/JAVIERRE_GWAS/out/pmi/SLE.pmi"),
            UC=fread("/home/oliver/JAVIERRE_GWAS/out/pmi/UC.pmi")
            )
#multiple sclerosis, primary billiary cirrhosis, RA,systemic lupus erythematosus , ulcerative colitis and Crohn’s disease

for(nm in names(pmi)) {
    pmi[[nm]][,group:="A"]
    pmi[[nm]][,disease:=nm]
    setnames(pmi[[nm]],c("rsid","start"),c("snp","pos"))
}

source(file.path(CD4CHIC.ROOT,"activation-analyses/R/common.R"))
b2g <- get.b2gene()


## make plot
## common objects
hind <- get.hind()
hind <- as.data.table(as.data.frame(hind))
hind[,seqnames:=sub("chr","",as.character(seqnames))]
hind[,start:=start-1]
hind[,end:=end-1]
axis <- GenomeAxisTrack(col="black", fontcolor="black")
    genes <- BiomartGeneRegionTrack(genome="hg19", name="Genes", transcriptAnnotation="symbol",
                                    mart=e75.genemart,
                                    collapseTranscripts="meta",
                                    stackHeight=0.2#, filters=list(with_ox_refseq_mrna=T)
                                    )
    rnan <- rnat("data_srt.2_non.bam", "RNA-seq (n)", rlim)
    rnaa <- rnat("data_srt.2_act.bam", "RNA-seq (a)", rlim)
    ac27n <- chpt("H3K27ac_PoolQTW_NAct.bam", "H3K27ac (n)", clim, "green4")
    ac27a <- chpt("H3K27ac_PoolQTW_Act.bam", "H3K27ac (a)", clim, "mediumpurple4")
ints <- get.interactions()
ints[,oeChr:=sub("chr","",oeChr)]
ints[,baitChr:=sub("chr","",baitChr)]
dat.hicn <- ints[Total_CD4_NonActivated>5,]
dat.hica <- ints[Total_CD4_Activated>5,]


## specific objects
## specific params
## > IRF8  irf8 (16:85928981-86053661) use the Crohn's disease IChip data to show disease association.
## >
## > AHR ahr (7:16400468-18016562) use RA (Okada) GWAS data to show disease association.
## >
## > CCR7 ccr7 (17:38433123-38804536) use T1D IChip data to show disease association.

params <- list(## il2ra=list(snps=read.table("il2ra-gwas.csv", header=T),
               ##            f.hicn="hic-non-il2ra-s.csv",
               ##            f.hica="hic-act-il2ra-s.csv",
               ##            name="il2ra",
               ##            chromosome="10",
               ##            from=6040000,
               ##            to=6160000),
## ptprc=list(disease=c("MS","CEL","T1D"),
##             chromosome="1",
##             from=196000000,
##             to=201000000),
trove2=list(disease=c("CEL","MS","T1D"),
            chromosome="1",
            from=192400000,
            to=193120000),
               ccr7=list(disease="T1D",
                         chromosome="17",
                         from=38433123,
                         to=38804536),
               irf8=list(disease=c("CRO","MS","PBC","RA","SLE","UC"),
                         chromosome="16",
                         from=85790000,
                         to=86100000),
               ahr=list( disease="RA",
                         chromosome="7",
                         from=16900000,
                        to=17416562),
  ## add TROVE2/T1D
  ## trove2=list(disease="T1D",
  ##             chromosome="1",
  ##             from=192400000,
  ##             to=193120000),
  ## add PTPRC/T1D
  ptprc=list(disease="T1D",
             chromosome="1",
             from=196000000,
             to=198800000),
  il10=list(disease=c("T1D", "CRO", "SLE","UC"),
            chromosome="1",
            from=206934000,
            to=207136000))
params$cox4i1 <- params$irf8
params$rara <- params$ccr7
params$il19 <- params$il20 <- params$il24 <- params$pigr <- params$il10
params$rgs1 <- params$trove2
## auto add snps, baits
for(nm in names(params)) {
    message(nm)
    params[[nm]]$name <- nm
    params[[nm]]$snps <- lapply(params[[nm]]$disease,function(d) {
        pmi[[d]][chr==params[[nm]]$chromosome &
                 pos >= params[[nm]]$from &
                 pos <= params[[nm]]$to &
                 ppi > 0.01,
                 .(disease,group,snp,chr,pos,pval,ppi)] }) %>% do.call("rbind",.)
    params[[nm]]$baitID <- with(params[[nm]],
                              unique(b2g[gene==toupper(name),]$baitID))
}

params[["irf8"]]$snps[,group:=disease]
params[["irf8"]]$disease <- "GWAS"

params[["trove2"]]$snps[,group:=disease]
params[["trove2"]]$disease <- "GWAS"

params[["il10"]]$snps[,group:=disease]
params[["il10"]]$disease <- "GWAS"

## generate plot components
TRACKS <- vector("list",length(params))
names(TRACKS) <- names(params)

nm <- "trove2"

for(nm in names(params)) {
    while(length(w <- which(search()=="params[[nm]]"))) 
        detach(pos=w[1])
    message(nm)
    attach(params[[nm]])
    tracks <- list()
    tracks$gwast <- AnnotationTrack(range=data.frame(chromosome=snps$chr, start=snps$pos, end=snps$pos), ## genome="ENSEMBL",
                                    stacking="squish",
                                    name=disease,
                             group=snps$group, groupAnnotation="group",
                             col.line=NA, stackHeight=0.15, showOverplotting=T, col=NA,
                             feature=snps$group, A="brown", B="blue", E="cyan")
    tracks$hd3 <- hindt(chromosome, from, to)
    tracks$hicn <- hict(dat.hicn, "PCHi-C (n)", chromosome, from, to, baitID, "green4")
    tracks$hica <- hict(dat.hica, "PCHi-C (a)", chromosome, from, to, baitID, "mediumpurple4")
    tracks$reg <- rgnt("chrom.bed", name="regRNA", chromosome)
    TRACKS[[nm]] <- tracks
    detach(params[[nm]])
}

plotter("il10")

plotter("ptprc")
plotter("trove2")

## for IL10, want variant view ie
TRACKS$il10v <- TRACKS$il10
params$il10v <- params$il10
attach(params$il10)
TRACKS$IL10v$hicn <- hict(dat.hicn, "PCHi-C (n)", chromosome, from, to,
                         i.baitID=unlist(sapply(params[c("il10","il19","il20","il24","pigr")], "[[", "baitID")),
                         "green4")
TRACKS$IL10v$hica <- hict(dat.hica, "PCHi-C (n)", chromosome, from, to,
                         i.baitID=unlist(sapply(params[c("il10","il19","il20","il24","pigr")], "[[", "baitID")),
                         "mediumpurple4")
detach(params$il10)
plotter("il10v")

a <- 1
w <- 11 * a
h <- 8 * a

doplot <- list("irf8"=list("irf8","cox4i1"),
               "ptprc"=list("ptprc"),
               trove2=list("trove2","rgs1"),
               ahr=list("ahr"),
               ccr7=list("ccr7","rara"),
               il10=list("il10","il19","il20","il24","faim3","pigr"),
               il10v=list("il10v"))

for(nm in names(doplot)) {
    ## f <-file.path(CD4CHIC.OUT,paste0("paper/",nm,".pdf")); message(f)
    ## pdf(f,height=h,width=w)
    ## do.call("plotter",doplot[[nm]])
    ## dev.off()

    ## svg(sub("pdf","svg",f),height=h,width=w)
    ## do.call("plotter",doplot[[nm]])
    ## dev.off()
    f <-file.path(CD4CHIC.OUT,paste0("paper/",nm,"-square-allgwas.pdf")); message(f)
    if(nm=="il10") {
        cairo_pdf(f,height=h*1.5,width=h)
    } else {
        cairo_pdf(f,height=h,width=h)
    }
    plotter(nm, plus=doplot[[nm]][-1])
    dev.off()

    ## svg(sub("pdf","svg",f),height=h,width=h)
    ## do.call("plotter",doplot[[nm]])
    ## dev.off()

    ## ppi <- 300
    ## png(sub("pdf","png",f), height=h*ppi, width=w*ppi, res=ppi)
    ## do.call("plotter",doplot[[nm]])
    ## dev.off()
    
}
