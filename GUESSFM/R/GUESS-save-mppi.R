
library(annotSnpStats)
library(devtools)
pdir="~/RP/GUESSFM"
load_all(pdir)
library(magrittr)
library(reshape)

files <- list.files("~/FM-GUESS",recursive=TRUE,pattern="snpmod-fixed2.RData",full=TRUE)
DATA <- vector("list",length(files))
for(i in seq_along(files)) {
    f <- files[i]
    message(f)
    (load(f)) # SM2
    mppi <- lapply(SM2, function(x) x@snps[,c("var","Marg_Prob_Incl")])
    for(nm in names(mppi))
        mppi[[nm]]$trait <- nm
    mppi <- do.call("rbind",mppi)
    mppi$region <- basename(dirname(f))
    DATA[[i]] <- mppi
    colnames(DATA[[i]]) %<>% sub("Marg_Prob_Incl","mppi",.)
}

DATA <- do.call("rbind",DATA)
MPPI <- DATA
save(MPPI,file="~/FM/mppi-fixed2.RData")

