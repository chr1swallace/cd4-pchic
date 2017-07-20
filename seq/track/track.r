library(Gviz)
library(GenomicInteractions)
source("stranded.r")
options(ucscChromosomeNames=FALSE)

hict <- function(file, name, chromosome, from, to, col) {
    rdat <- read.table(file, header=T)
    data <- rdat[rdat$baitChr == chromosome & rdat$baitStart >= from & rdat$baitEnd <= to &
                 rdat$oeChr == chromosome & rdat$oeStart >= from & rdat$oeEnd <= to,]
    bait <- GRanges(data$baitChr, IRanges(start=data$baitStart, end=data$baitEnd))
    oe <- GRanges(data$oeChr, IRanges(start=data$oeStart, end=data$oeEnd))
    interaction <- GenomicInteractions(bait, oe, counts=1)
    hic <- InteractionTrack(name=name, interaction, chromosome=chromosome)
    displayPars(hic)$anchor.height <- 0.2
    displayPars(hic)$col.anchors.line <- "lightgrey"
    displayPars(hic)$col.anchors.fill <- "orange"
    displayPars(hic)$col.interactions <- col
    hic
}

rnat <- function(file, name, ylim) {
    logt <- function(x) sign(x) * log2(abs(x)+0.1)
    DataTrack(range=file, importFunction=strandedBamImport, stream=T, baseline=0, col.baseline="grey",
              col=c("#e31a1c","#1f78b4"), groups=c("Forward", "Reverse"), type="histogram", col.histogram=NA,
              name=name, ylim=ylim, transformation=logt)
}

chpt <- function(file, name, ylim) AlignmentsTrack(file, isPaired=F, type=c("coverage"), name=name, ylim=c(0, ylim))

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

plotregion <- function(name, chromosome, from, to, gwast, hicn, hica, rlim, clim) {
    pdf(paste(name, "pdf", sep="."), paper="a4")
    axis <- GenomeAxisTrack(col="black", fontcolor="black")
    genes <- BiomartGeneRegionTrack(genome="hg19", name="Genes", transcriptAnnotation="symbol",
                                    collapseTranscripts="longest",
                                    stackHeight=0.2, filters=list(with_ox_refseq_mrna=T))
    hicn <- hict(hicn, "Hi-C (n)", chromosome, from, to, "green4")
    hica <- hict(hica, "Hi-C (a)", chromosome, from, to, "mediumpurple4")
    rnan <- rnat("data_srt.2_non.bam", "RNA-seq (n)", rlim)
    rnaa <- rnat("data_srt.2_act.bam", "RNA-seq (a)", rlim)
    ac27n <- chpt("H3K27ac_PoolQTW_NAct.bam", "H3K27ac (n)", clim)
    ac27a <- chpt("H3K27ac_PoolQTW_Act.bam", "H3K27ac (a)", clim)
    reg <- rgnt("chrom.bed", name="regRNA", chromosome)
    seqt <- list(hicn, hica, rnan, rnaa, ac27n, ac27a, reg)
    trt <- c(gwast, genes, seqt)
    ht <- HighlightTrack(trackList=trt, chromosome=chromosome, range=range(gwast), col="grey", alpha=0.25)
    plotTracks(c(axis, ht),
               chromosome=chromosome, from=from, to=to, sizes=c(0.3, 0.3, 0.3, 0.3, 0.3, 0.5, 0.5, 0.5, 0.5, 0.5),
               background.title=NA, col.axis="darkgrey", fontcolor.title="black", col.axis="black")
    dev.off()
}

snps <- read.table("il2ra-gwas.csv", header=T)
gwast <- AnnotationTrack(range=data.frame(chromosome=snps$chr, start=snps$pos, end=snps$pos), genome="ENSEMBL", name="GWAS",
                         group=snps$group, groupAnnotation="group", col.line=NA, stackHeight=0.15, showOverplotting=T, col=NA,
                         feature=snps$group, A="brown", B="blue", E="cyan")
displayPars(gwast)$C <- "black"
displayPars(gwast)$D <- "magenta"
displayPars(gwast)$F <- "green"
plotregion(name="il2ra", chromosome="10", from=6040000, to=6160000, gwast=gwast, hicn="hic-non-il2ra-s.csv", hica="hic-act-il2ra-s.csv",
           rlim=c(-10, 10), clim=220)

pos <- c(206964952, 206944645, 206955041, 206943968, 206942413, 206939904)
gw <- AnnotationTrack(range=data.frame(chromosome="1", start=pos, end=pos), genome="ENSEMBL", name="GWAS", stackHeight=0.15, col=NA)
plotregion(name="il10", chromosome="1", from=206930000, to=207137000, gwast=gw, hicn="hic-non-il10-s.csv", hica="hic-act-il10-s.csv",
           rlim=c(-12, 12), clim=220)
