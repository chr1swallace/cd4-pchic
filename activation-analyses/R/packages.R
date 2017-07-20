##' Recipe of which packages need to be installed from where for code in this directory to work

##' CRAN
install.packages(c("data.table","magrittr","devtools"))
install.packages(c("ggplot2","dplyr"))
# Error: package 'Rcpp' 0.11.6 was found, but >= 0.12.0 is required by 'dplyr'
install.packages("Rcpp")
install.packages("dplyr")
install.packages("cowplot")
install.packages("pheatmap")
# ERROR: dependencies 'impute', 'preprocessCore', 'GO.db' are not available for package 'WGCNA'
install.packages("rms")

##' BIOCONDUCTOR
source("http://bioconductor.org/biocLite.R")
biocLite("ggbio")
biocLite("GenomicRanges")
biocLite("preprocessCore")
biocLite(c("impute","GO.db"))

##' CRAN with BIOC dependencies
install.packages("WGCNA")
