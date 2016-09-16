library(abind)
library(pipeR)
read.imputeinfo <- function(f.info) {
  info <- mclapply(f.info,read.table,header=TRUE,as.is=TRUE)
  data <- lapply(info,"[",-c(1:3)) %>>% { abind(., along=0) }
  data.mean <- apply(data, c(2,3), mean)
  data.sd <- apply(data, c(2,3), sd)
  return(cbind(info[[1]][,1:3],data.mean))
}

read.geno <- function(f.geno,samples,alleles,maxafield=115,maxnfield=63) {
    message("reading ",f.geno)
    ## any alleles > 125 characters (large indels)
    n1 <- nchar(alleles$A1)
    n2 <- nchar(alleles$A2)
    nr <- nchar(alleles$rs_id)
    if(max(c(n1,n2))>=maxafield || max(nr)>maxnfield) {
        message("fixing long alleles")
        f.temp <- tempfile(tmpdir="./tmp")
        ff <- readLines(f.geno)
        wh <- which(pmax(n1,n2)>maxafield | nr>maxnfield)
        ss <- strsplit(ff[wh]," ")
        ss.orig <- ss
        ss <- lapply(ss,function(x) {
            x[[2]] <- substr(x[[2]],1,maxnfield)
            x[[4]] <- substr(x[[4]],1,maxafield)
            x[[5]] <- substr(x[[5]],1,maxafield)
            return(x)})
        ff[wh] <- sapply(ss,paste,collapse=" ")
        writeLines(ff,f.temp)
        X <- read.impute(f.temp,rownames=samples,nsnp=nrow(alleles))
        colnames(X)[wh] <- sapply(ss.orig,"[[",2)
        unlink(f.temp)
    } else {
        X <- read.impute(f.geno,rownames=samples,nsnp=nrow(alleles))
    }
    return(X)
}

## sample info
read.sinfo <- function(d.pq,chrarm) {
  f.allsamples <- paste0(d.pq,"/samples-",chrarm,".csv.gz")
  if(!file.exists(f.allsamples))
    return(NULL)
  sinfo <- read.table(f.allsamples,header=TRUE,row.names=1,sep="\t",as.is=TRUE)
  ra.map <- c("RA.RA1"="RA.ES",
              "RA.RA2"="RA.NL",
              "RA.RA3"="RA.SEE",
              "RA.RA4"="RA.SEU",
              "RA.RA5"="RA.UK",
              "RA.RA6"="RA.US")
  wh <- which(sinfo$file %in% names(ra.map))
  sinfo[wh,"file"] <- ra.map[sinfo[wh,"file"]]
  table(sinfo$file)
  sinfo$country <- ifelse(grepl("RA.",sinfo$file,fixed=TRUE),sub("RA.","",sinfo$file),"UK")
  sinfo$phenotype <- ifelse(sinfo$affected==1,"CONTROL",sub("\\..*","",sinfo$file))
  return(sinfo)
}
