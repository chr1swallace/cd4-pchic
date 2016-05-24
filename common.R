library(magrittr)

pmean <- function(...) {
  x <- do.call("cbind",list(...))
  rowSums(x)/ncol(x)
}
get.rnaseq <- function() {
    rnaseq <- fread(file.path(CD4CHIC.EXTDATA,"rnaseq.csv"))
    #rnaseq<-fread(file.path(CD4CHIC.DATA,"rna-kallisto-gene-diff.csv"),header=TRUE,stringsAsFactors=FALSE,sep="\t")
    rnaseq[, c("act.total","non.total") := list(pmean(act_1,act_2), pmean(non_1,non_2)) ]
    return(rnaseq)
}

make.colScale <- function(mycols) {
  require(ggplot2)
  mycols <- as.character(mycols)
  names(mycols) <- mycols
  mycols[mycols=="NA"] <- "grey20"
  mycols["yellow"] <- "goldenrod"
  mycols["brown"] <- "SaddleBrown"
  mycols["grey"] <- "grey40"
  mycols["green"] <- "DarkGreen"
  mycols["blue"] <- "DarkBlue"
  scale_colour_manual(name = "module",values = mycols, na.value="purple",guide=FALSE)
}

Cochran.Armitage.test <- function(exposure, cc, stratum=rep(1,length(cc)),quick=FALSE) {
  N <- length(cc)
  if (is.factor(exposure))
    exposure <- as.numeric(exposure)
  cl <- match.call()
  narg <- length(cl) - 1
  arguments <- character(narg-1)
#   for (i in 1:narg)
#     arguments[i] <- as.character(cl[[i+1]])
  use <- complete.cases(cc,exposure,stratum)
  if (any(!use)) {
    cc <- cc[use]
    exposure <- exposure[use]
    stratum <- stratum[use]
  }
  dh <- table(stratum, cc)
  if (ncol(dh)!=2)
    stop("cc argument must have two levels")
  nt <- table(stratum)
  zm <- tapply(exposure, list(stratum, cc), mean)
  ut <- dh[,1]*dh[,2]*(zm[,2] - zm[,1])/nt
  if(quick)
  	return(sum(ut,na.rm=T))
  zv <- tapply(exposure, stratum, var)
  vt <-  dh[,1]*dh[,2]*zv/nt
  x2 <- sum(ut,na.rm=T)^2/sum(vt,na.rm=T)
  names(x2) <- "Chi-squared"
  df <- 1
  names(df) <- "df"
  res <- list(statistic=x2, parameter=df,
              p.value=pchisq(x2, 1, lower.tail=FALSE),
              method="Cochran-Armitage test with Mantel's extension",
              data.name=paste(cl,collapse=" "),
              score=ut,
              sum.score=sum(ut,na.rm=T),
              score.variance=vt,
              sum.variance=sum(vt,na.rm=T))
  class(res) <- "htest"
  res
}
