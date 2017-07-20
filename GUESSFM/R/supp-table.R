#!/usr/bin/env Rscript
library(data.table)
agg <- function(patt) {
    files <- list.files("~/scratch/FM-GUESS",recursive=TRUE,pattern=paste0(patt,".RData"),full=TRUE)
    VEC <- vector("list",length(files))
    for(i in seq_along(files)) {
        f <- files[i]
        message(f)
        obj <- eval(as.symbol(load(f)))
        obj <- subset(obj, trait %in% c("ICOELIAC","RA","T1D","GRAVES"))
        obj <- as.data.table(obj)
        obj$region <- basename(dirname(f))
        VEC[[i]] <- obj
    }
    do.call("rbind",VEC)
}

## OUTPUT FILE:
summx <- agg("summx")
(load(file = "~/FM/manhattan.RData"))
(load(file="~/FM/nsnps-fixed2.RData"))
(load(file="~/FM/mppi-fixed2.RData"))

m <- merge(subset(MPPI,var!=1),
           MANN,
           by.x=c("var","trait","region"),by.y=c("snp","trait","region"),all.x=TRUE)
colnames(m) <- sub("var","snp",colnames(m))

with(subset(m, grepl("19p",MPPI$region)),tapply(mppi,trait,max))



dim(summx)
m <- merge(m,summx[,.(snp,tag,trait,pp,ppsum)],all.x=TRUE,by=c("snp","trait"))
m <- as.data.table(m)
m[,chr:=as.integer(sub("[pq]-.*","",region))]
m <- m[order(m$chr,m$position),]
m <- subset(m,trait %in% c("T1D","ICOELIAC","RA","GRAVES"))
m[trait=="ICOELIAC",trait:="CEL"]
m[trait=="GRAVES",trait:="ATD"]
write.table(m[,.(chr,position,snp,snp_id,trait,p,exp_freq_a1,info,A1,A2,region,mppi,tag,ppsum)],
            file="/scratch/cew54/cd4chic/figures/sdata-fine-mapping-results.csv",
            row.names=FALSE)
## how many regions have > 1 disease signal?
## see GUESS-roundup.R
