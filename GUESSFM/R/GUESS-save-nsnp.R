
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
    nsnp <- pp.nsnp(SM2,expected=2)
    save(nsnp,file=sub("snpmod","nsnp",f))
    data <- do.call("rbind",nsnp@.Data) %>% as.data.frame()
    data$trait <- nsnp@traits
    data$region <- basename(dirname(f))
    DATA[[i]] <- melt(data, variable_name="nsnp")
    colnames(DATA[[i]]) %<>% sub("value","pp",.)
}

DATA <- do.call("rbind",DATA)
DATA$nsnp <- as.numeric(as.character(DATA$nsnp))

library(ggplot2)
DATA$region.trait <- paste(DATA$region,DATA$trait,sep=".")
ggplot(DATA,aes(x=nsnp,y=pp,group=region.trait,col=trait)) + geom_path()

save(DATA,file="~/FM/nsnps-fixed2.RData")
save(DATA,file="/scratch/cew54/FM/nsnps-fixed2.RData")
