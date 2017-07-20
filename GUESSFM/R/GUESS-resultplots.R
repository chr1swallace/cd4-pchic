#!/usr/bin/env Rscript
## follows on from GUESS-qc.R, may need some code from that inserted

library(randomFunctions)
library(annotSnpStats)
library(magrittr)
library(RColorBrewer)
library(reshape)
library(plyr)
library(igraph)
library(ggplot2)
theme_set(theme_bw())
library(gridExtra)
setwd("~/Projects/cd4chic/GUESSFM")
source("~/DIRS.txt")

library(annotSnpStats)

## files
args <- getArgs(default=list(d=file.path(ROOT,"FM-GUESS/2q-204446380-204816382"))) #10p-6030243-6169685 1q-172650685-172940450
 ## 1: rs4646925.34852546.A.G     86            --- rs113313247:34897582:TTCTC:T
 ## 2: rs4646925.34852546.A.G     79            ---   rs150252190:34840426:CTG:C
 ## 3: rs4646925.34852546.A.G     82 imm_6_34960524       rs4646925:34852546:A:G

## functions and set some parameters
source("R/GUESS-common.R")
devtools::load_all("~/RP/GUESSFM")
(load(f.geno))

DISEASES <- intersect(#c("ICOELIAC","COELIAC","RA","T1D","GRAVES","JIA","MS", "RA.UK","RA.US"),
                      c("ICOELIAC","RA","T1D","GRAVES","MS","JIA"),
                      names(SM2))
SM2 <- SM2[DISEASES]
SP2 <- SP2[DISEASES]

## ## union of SNP groups
## s <- "rs10912265.172673583.C.T" #is problem snp
## ## issue occurs when adding i=3
## gSP2 <- SP2[[1]]
##   gSP2 <- GUESSFM::union(gSP2,SP2[[2]])
##    x=gSP2
##  y=SP2[[3]]
## ## union?
## as(x,"groups")
## as(y,"groups")
## (z <- as(GUESSFM::union(x,y),"groups"))
## str(z)
## length(z@tags)
## length(unique(z@tags))
## snpin(s,z)

gSP2 <- SP2[[1]]
if(length(SP2)>1)
  for(i in 2:length(SP2)) {
    gSP2 <- GUESSFM::union(gSP2,SP2[[i]])
    tg<-as(gSP2,"tags")    
    message(i,"\t",length(tg@.Data),"\t",length(unique(tg@.Data)))
 }
plot(gSP2)
save(gSP2,file=file.path(args$d,"snp-picker-grouped2.RData"))
groups2 <- as(gSP2,"groups")
if(args$d==file.path(ROOT,"FM-GUESS/6p-34507311-35327577")) {
    ## bad merge in GRAVES, fix
    groups2@.Data <- c(lapply(groups2@.Data,setdiff,"rs4646925.34852546.A.G"),list("rs4646925.34852546.A.G"))
    groups2@tags <- c(sub("rs4646925.34852546.A.G","rs9469923.34850411",groups2@tags),"rs4646925.34852546.A.G")
}

## results=SM2
## groups=groups2
## snps=DATA@snps
## position="position"
## tag.thr=0.8
## pp.thr=0.01
## method="complete"

summx <- guess.summ(SM2,groups=groups2,snps=DATA@snps,position="position")
summx <- scalepos(summx,position="position")
chr <- sub("[pq].*","",basename(args$d))
summx$chr <- paste0("chr",chr)
save(summx, file=file.path(args$d,"summx2.RData"))


manhattan <- do.call("rbind",MANHATTAN)
manhattan <- merge(manhattan,DATA@snps[,c("type"),drop=FALSE],by.x="snp",by.y=0,all.x=TRUE)
if(all(grepl("samples",colnames(samples(DATA))))) {
  s <- samples(DATA)
  colnames(s) <- sub("samples.","",colnames(s))
  samples(DATA) <- s
}
data.uk <- as(DATA[samples(DATA)$country=="UK" & samples(DATA)$phenotype=="CONTROL", ], "SnpMatrix")
man <- subset(manhattan, trait %in% DISEASES & !is.na(p))
man <- merge(man,summx[,c("snp","tag")],all.x=TRUE)
man <- man[c(which(is.na(man$tag)),which(!is.na(man$tag))),]

library(ggplot2)
library(gtable)
library(grid) # low-level grid functions are required
## gene locations
(load(file.path(args$d,"genes.RData")))
any.genes <- length(genes)>0
if(any.genes) {
    library(ggbio)
    p.genes <- ggplot(genes) + geom_alignment(aes(group=genename)) + scale_x_sequnit("Mb")
}
revaxis <- function(p) {
# extract gtable
g <- ggplot_gtable(ggplot_build(p))
    xl <- which(g$layout$name == "axis-l")
    xr <- which(g$layout$name == "axis-r")
    ax <- g$grobs[[xl]]
    ## ax <- g$grobs[[xl]]$children[[2]]
    ## ax$widths <- rev(ax$widths)
    ## ax$grobs <- rev(ax$grobs)
    #ax$grobs[[1]]$x <- ax$grobs[[1]]$x - unit(1, "npc") + unit(0.15, "cm")
    g$grobs[c(xl,xr)] <- g$grobs[c(xr,xl)]
#    g$grobs[xr]$childrenOrder <- rev(g$grobs[xr]$childrenOrder)
    g$grobs[xr] <- t(g$grobs[xr])
    g$widths[c(2,3,5,6)] <- g$widths[c(6,5,3,2)]
    ## gtable_show_layout(g)
    ## plot(g)
## ia <- which(g$layout$name == "axis-l")
## ax <- g$grobs[[ia]]$children[[2]]
## ax$widths <- rev(ax$widths)
## ax$grobs <- rev(ax$grobs)
## ax$grobs[[1]]$x <- ax$grobs[[1]]$x - unit(1, "npc") + unit(0.15, "cm")
## pp <- c(subset(g$layout, name == "panel", select = t:r))
## g <- gtable_add_cols(g, g$widths[g$layout[ia, ]$l], length(g$widths) - 1)
## g <-  gtable_add_grob(g, ax, pp$t, length(g$widths) - 1, pp$b)
## g$grobs[[ia]]$children[[2]] <- NULL
## ##############################
## ia <- which(g$layout$name == "ylab-l")
## ylab <- g$grobs[[ia]]
## #g <- gtable_add_cols(g, g$widths[g$layout[ia, ]$l], length(g$widths) - 1)
## g <-  gtable_add_grob(g, ylab, pp$t, length(g$widths) - 2, pp$b)
## g$grobs[[ia]]$label = ''
    g
}


plotter <- function(man,summx) {

    xl <- min(man$position,na.rm=TRUE)
    xm <- max(man$position,na.rm=TRUE)

    pvals <- ggplot(man,aes(x=position,y=-log10(p),col=tag)) + geom_point() + facet_grid(trait ~ .) + addlines(summx) + geom_hline(yintercept=8,linetype="dotted",col="grey") + theme(legend.position="none",strip.text.y = element_text(size = 9, angle = 0)) + ggtitle(REGION)
signals <- signal.plot(summx)
chr <- ggchr(summx)
snps <- ggsnp(summx)
lds <- ggld(data.uk, summx)

nd <- length(DISEASES)
if(any.genes) {
    G <- list(p.genes + ggtitle(paste(REGION,NAME)),pvals,chr,signals,snps)
    h <- c(2,rep(5/nd,nd),1,rep(3/nd,nd),1)
} else {
    G <- list(pvals,chr,signals,snps)
    h <- c(rep(5/nd,nd),1,rep(3/nd,nd),1)
}
if(nrow(summx)>1) {
    h <- c(h,list(2))
    G <- c(G,list(lds))
}
    G %<>% lapply(., function(p) p + xlim(xl,xm) +
                                 theme(axis.title.x = element_blank(), 
                                       axis.text.x = element_blank(), 
                                       axis.ticks.x= element_blank(), 
                                       plot.margin= unit(c(1, 1, -0.5, 0.5), "lines")))
if(any.genes) {
    G[-1] <- lapply(G[-1], ggplotGrob)
    G[[1]] <- revaxis(G[[1]])
} else {
    G <- lapply(G,ggplotGrob)
}
maxcol <- sapply(G,ncol) %>% max()
while(any(sapply(G,ncol)< maxcol)) {
    wh <- which(sapply(G,ncol) < maxcol)
    if(length(wh)) {
        for(w in wh) {
            G[[w]] <- gtable_add_cols(G[[w]], unit(0,"mm"))
        }
    }
}
g <- do.call("rbind",c(G,size="max")) # stack the plots
#g$widths <- G[[1]]$widths #unit.pmax(unlist(lapply(G, function(g) g$widths))) # use the largest widths
# center the legend vertically
g$layout[grepl("guide", g$layout$name),c("t","b")] <- c(1,nrow(g))
panels <- g$layout$t[grep("panel", g$layout$name)]
# set the relative panel heights 1/3 for the top two
g$heights[panels] <- unit(h, "null")
grid.newpage()
grid.draw(g)
}


f <- start.plot("summx2",height=10)
plotter(man,summx)
dev.off()


snps.use <- subset(summx,
                   trait %in% DISEASES &
                   pp>0.01,
                   select="snpnum",drop=TRUE) %>%
  unique() %>% setdiff(., NA)

if(length(snps.use)) {
    summx2 <- subset(summx,snpnum %in% snps.use & trait %in% DISEASES)
summx2 <- GUESSFM:::summ.setminmax(summx2)
summx2 <- scalepos(summx2,position="position")
man2 <- subset(manhattan, trait %in% DISEASES & !is.na(p))
man2 <- merge(man2,summx2[,c("snp","tag")],all.x=TRUE)
man2 <- man2[c(which(is.na(man2$tag)),which(!is.na(man2$tag))),]
f <- start.plot("summx2-001",height=10)
plotter(man2,summx2)
dev.off()
}

bm <- best.snps(SM2,0.03)
  hsnps <- lapply(bm[DISEASES],rownames) %>% unlist() %>% unique()
if(length(hsnps)) {
    summx3 <- subset(summx,snp %in% hsnps & trait %in% DISEASES)
summx3 <- GUESSFM:::summ.setminmax(summx3)
summx3 <- scalepos(summx3,position="position")
man2 <- subset(manhattan, trait %in% DISEASES & !is.na(p))
man2 <- merge(man2,summx3[,c("snp","tag")],all.x=TRUE)
man2 <- man2[c(which(is.na(man2$tag)),which(!is.na(man2$tag))),]
##plotter(man2,summx3)
f <- start.plot("summx2-best",height=10)
plotter(man2,summx3)
dev.off()
}
## doplot("manhattan-bestsnps",pvals + theme(legend.position="top",legend.direction="vertical"))
## doplot("signals-bestsnps",signals + theme(legend.position="top",legend.direction="vertical"))
## doplot("snps-bestsnps",snps)
## doplot("summx-bestsnps",p,height=10)

tmp <- summx2[,c("tag","snp","position","trait","pp","ppsum")]
w1 <- reshape::recast(tmp, tag + snp + position ~ trait,  measure.var="pp")
w2 <- reshape::recast(tmp, tag + snp + position ~ trait,  measure.var="ppsum")
colnames(w2)[-c(1:3)] <- paste0(colnames(w2)[-c(1:3)],".sum")
w <- merge(w1,w2)
max.sum <- apply(w[,grep(".sum$",colnames(w)),drop=FALSE],1,max)
w <- w[order(w$tag,w$position),]
write.table(w, file=file.path(plot.dir,"summw2.csv"))


## if(interactive())
##   system(paste("rsync -av ",plot.dir," chrisw@salthouse:~/Projects/FineMapping/output"))
