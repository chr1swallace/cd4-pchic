library(WGCNA)
library(data.table)
library(parallel); options(mc.cores=2)
library(ggplot2); theme_set(theme_bw())
library(magrittr)
library(gridExtra)
library(reshape)
setwd(CD4CHIC.ROOT)

options(stringsAsFactors = FALSE);
enableWGCNAThreads(10)
source("common.R")

(load(file.expression))

######################
### NETWORK ANALYSIS #
######################

## do some rudimentary cleaning
check.ok <- function(x) {
    gsg = goodSamplesGenes(x, verbose = 3)
    message(paste("All OK",gsg$allOK,sep="=="))
    invisible()
}
lapply(s.exp,check.ok)


## pick a soft threshold
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = mclapply(s.exp, pickSoftThreshold,
    powerVector = powers, verbose = 5,networkType=network.type)

## plot results
df <- lapply(sft, "[[", "fitIndices") %>% do.call("rbind",.)
df$samples <- ifelse(grepl("all",rownames(df)), "all","treated")
p1 <- ggplot(df, aes(x=Power,y=-sign(slope)*truncated.R.sq)) +
geom_path() +
geom_path(aes(y=-sign(slope)*SFT.R.sq),col="blue") +
facet_wrap(~samples) +
geom_hline(yintercept=0.9,col="red") + ggtitle("Scale-free topology fit index as a function of the soft-thresholding power")

p2 <- ggplot(df, aes(x=Power,y=-sign(slope)*truncated.R.sq)) +
geom_path() +
geom_path(aes(y=-sign(slope)*SFT.R.sq),col="blue") +
facet_wrap(~samples) + ggtitle("Mean connectivity as a function of the soft-thresholding power")
soft.threshold.plot <- file.path(plot.dir,
                                 paste0(file.prefix,'_',network.type,'_softThreshold.pdf'))
pdf(file=soft.threshold.plot,width=8,height=8)
grid.arrange(p1,p2,nrow=2)
dev.off()

## spower<-lapply(sft, function(x) {
##     xx <- x$powerEstimate
##     if(is.na(xx)) {
##         ## defaults from http://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/faq.html
##         return(switch(network.type, "signed"=12, "unsigned"=6))
##     }
##     return(xx)
## })
## spower <- unlist(spower) %>% max()
## spower <- 22
## message("choosing power: ",spower)

powers <- c(8,12,22) # manually chosen from graph

################################################################################

getmods <- function(x,spower) {
    blockwiseModules(x, maxBlockSize = 20000,
                     networkType=network.type,
                     power = spower, TOMType = network.type, minModuleSize = 30,
                     reassignThreshold = 0, mergeCutHeight = 0.25,
                     numericLabels = TRUE,
                     saveTOMs = FALSE,
                     verbose = 3)
}
plotter <- function(bwnet,...) {
    mergedColors = labels2colors(bwnet$colors)
    ## Plot the dendrogram and the module colors underneath
    plotDendroAndColors(bwnet$dendrograms[[1]], mergedColors[bwnet$blockGenes[[1]]],
                        "Module colors",
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05,...)
}
getMEs <- function(exp,net) {
    cols <- labels2colors(net$colors)
    MEs <- moduleEigengenes(exp, cols)$eigengenes
    pheno <- pheno.data[match(rownames(exp), pheno.data$sname), pheno.cols ]
    df <- cbind(pheno,MEs)
}
nGenes = lapply(s.exp,ncol)
nSamples = lapply(s.exp,nrow)
pheno.cols <- c('uniqueID','pair','sex','time','treatment')

powers <- 12 # final, manual choice

for(spower in powers) {
    message("running spower: ",spower)
    f <- out.file("bwnet")
    if(file.exists(f)) {
        (load(f))
    } else {
        library(parallel)
        bwnet <- mclapply(s.exp, getmods, spower=spower, mc.cores=2)
        save(bwnet,file=f)
    }
    colors = lapply(bwnet, function(x) labels2colors(x$colors))
    message("gene counts per module")
    table(all=colors$all,treated=colors$treated)

    ## plot the dendrogram	
    pdf(plot.file("dendrogram"),width=12,height=9)
    lapply(names(bwnet), function(nm)
           plotter(bwnet[[nm]],main=paste(network.type,nm,"samples")))
    dev.off()

    ## Module eigengenes
    MEs <- mapply("getMEs",s.exp,bwnet,SIMPLIFY=FALSE)
    MEs$all$samples <- "all"
    MEs$treated$samples <- "treated"

    ## Module patterns
    
    df <- lapply(MEs, melt, c(pheno.cols,"samples")) %>% do.call("rbind", .)
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

    p <- ggplot(df,aes(x=time,y=value,col=treatment)) +
    geom_path(mapping=aes(group=id.treat),size=0.1,alpha=0.2) +
    facet_wrap(samples~variable) +
    geom_path(data=dfs,aes(x=time,y=value,colour=treatment,group=treatment),size=1) +
    scale_colour_manual("Treated",values=c("grey20","darkblue")) + ggtitle(paste("Eigengenes in each individual and averaged (thick line) by module,",network.type,"power =",spower))
    p

    ## save using all samples
    p <- ggplot(subset(df,samples=="all"),aes(x=time,y=value,col=treatment)) +
    geom_path(mapping=aes(group=id.treat),size=0.1,alpha=0.2) +
    facet_wrap(~variable) +
    geom_path(data=subset(dfs,samples=="all"),aes(x=time,y=value,colour=treatment,group=treatment),size=1) +
    scale_colour_manual("Treated",values=c("grey20","darkblue")) + ggtitle(paste("Eigengenes in each individual and averaged (thick line) by module,",network.type,"power =",spower))
    p
    
    
    pdf(plot.file("modulePattern"),height=8,width=8)
    print(p)
    dev.off()
    
}
