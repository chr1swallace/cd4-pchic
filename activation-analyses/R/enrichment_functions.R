## pathways
extract.genesets <- function(str,myid,gmt,minsize=10) {
    gmt.sub<-subset(gmt,grepl(str,gmt$set))
    nm <- gmt.sub[,"set"]
    genes <- gmt.sub[,"genes"] %>%
        strsplit(.,split="&") %>% lapply(.,unique) %>% lapply(.,intersect,myid)
    use <- sapply(genes,length) > minsize
    structure( genes[ use ], names=nm[ use ])
}
cattest_genesetlist <- function(test,setlist,universe) {
    tmp <- lapply(setlist, function(set) cattest_geneset(test,set,universe))
    tmp <- tmp[!sapply(tmp,is.null)]
    Z <- lapply(tmp,"[[",2) %>% do.call("rbind",.)
    LOR <- lapply(tmp,"[[",1) %>% do.call("rbind",.)
    list(LOR=LOR,Z=Z)
}
cattest_geneset<-function(test,set,universe){
    ## message("test: ",length(test))
    ## message("set: ",length(set))
    ## message("universe: ",length(universe))
    ## m <- sum(universe %in% set)
    ## n <- length(universe) - m
    ## k <- sum(test %in% universe)
    ## q <- sum(universe %in% test & universe %in% set)
    test <- intersect(test,universe)
    tt <- table(in.test=universe %in% test, in.set=universe %in% set)
    if(length(tt)<4 || nrow(tt)!=2 || ncol(tt)!=2 || any(tt<5))
        return(c(NA,NA))
    if(sum(tt)>100) {
        p <- chisq.test(tt)$p.value
    } else {
        p <- fisher.test(tt)$p.value
    }
    M <- matrix(c(1,-1,-1,1),nrow=2)
    LOR <- sum(log(tt) * M)
    c(LOR,qnorm(p/2,lower.tail=FALSE)*sign(LOR))
}
valtest_genesetlist <- function(test.id, test.val, setlist,...) {
    use <- test.id %in% unlist(setlist)
    tmp <- lapply(setlist, function(set) valtest_geneset(test.id[use],test.val[use],set,...))
    names(tmp) <- names(setlist)
    Z <- lapply(tmp,"[[",2) %>% do.call("rbind",.)
    EST <- lapply(tmp,"[[",1) %>% do.call("rbind",.)
    list(EST=EST,Z=Z)
}
valtest_geneset<-function(test.id,test.val,set,stat="W"){
    ## message("test: ",length(test))
    ## message("set: ",length(set))
    ## message("universe: ",length(universe))
    ## m <- sum(universe %in% set)
    ## n <- length(universe) - m
    ## k <- sum(test %in% universe)
    ## q <- sum(universe %in% test & universe %in% set)
    inset <- which(test.id %in% set)
    if( ## (length(set) - length(inset) < 5) || 
        ( length(inset) < 5 ) )
        return(c(NA,NA))
    kk <- t.test(test.val[inset], test.val[-inset])
    w <- wilcox.test(test.val[inset], test.val[-inset])
    p <- if(stat=="W") { w$p.value } else { kk$p.value }
    c(diff(kk$estimate),
      qnorm(p/2,lower.tail=FALSE) * sign(kk$statistic))
}

get.pathways <- function(myid) {
    L <- list(
	gmt_file = file.path(CD4CHIC.DATA,"ob_stats/RESOURCES","ensg.msigdb.v5.0.symbols.gmt"),
	protein_coding_genes = file.path(CD4CHIC.DATA,'ob_stats/RESOURCES/ensg2genename_protein_coding.tab'))

    pc<-unique(read.delim(L[['protein_coding_genes']],sep=",",stringsAsFactors=FALSE)[[1]])
    pc <- intersect(pc,myid)
    gmt<-read.delim(L$gmt_file,sep="\t",stringsAsFactors=FALSE)
    names(gmt)<-c('set','url','genes')
    return(gmt)
}


