#!/usr/bin/env Rscript
library(magrittr)
TRAITS <- c("RA","JIA","MS","ICOELIAC","GRAVES","MS","T1D")
## files <- list.files("~/FM-GUESS",full=TRUE)
## for(f in files) {
##     f.geno <- file.path(ROOT,"FM-impute_qc",paste0("imputed-",basename(f),".RData"))
##     if(!file.exists(f.geno))
##         next
##     (load(f.geno))
##     file.manhattan <- file.path(f,"manhattan.RData")
##     if(!file.exists(file.manhattan)) # should exist by now if anything interesting
##         next

library(data.table)
files <- list.files("~/FM-GUESS",recursive=TRUE,pattern="manhattan.RData",full=TRUE)
MANN <- vector("list",length(files))
for(i in seq_along(files)) {
    (load(files[[i]]))
    nc <- lapply(MANHATTAN,ncol)
    nc[ sapply(nc,is.null) ] <- 0
    use <- nc==11
    if(any(use)) {
        MANN[[i]] <- do.call("rbind",MANHATTAN[use])
        MANN[[i]]$region <- basename(dirname(files[[i]]))
    }
}
MANN <- do.call("rbind",MANN)
MANN <- as.data.table(MANN)
MANN$snp <- as.character(MANN$snp)

## bad snps
files <- list.files(".",pattern="snps-to-drop")
#bad <- lapply(files, function(f) eval(as.symbol(load(f))))
bad <- lapply(files,scan,what="") %>% unlist()

MANN <- MANN[!(snp %in% bad),]

save(MANN, file="~/FM/manhattan.RData")
MANN[region=="19p-10396336-10628468",]
