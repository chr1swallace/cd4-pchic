### final code for combining frag level prioritisation.

library(data.table)
library(ggplot2)
library(GenomicRanges)

DATA.DIR<-'/home/ob219/scratch/DATA/'
expressed.t.cells.only<-FALSE

## GWAS

gwas<-fread(file.path(DATA.DIR,'JAVIERRE_GWAS/out/highscore_fragments_JUNE_13/Autoimmune_13_06_17_gs_0.01.csv'))
gwas[gwas$Disease=="RA_OKADA_IMB",]$Disease<-'RA'
gwas$analysis<-'GWAS_PMI'

## ok filter so that only include those genes that are expressed in T-Cells

if(expressed.t.cells.only)
  gwas<-subset(gwas,!is.na(Log_FC_Gene_Expression_On_Activation))

## ICHIP GUESSFM

ichip.gf<-fread(file.path(DATA.DIR,'/JAVIERRE_ICHIP/out/highscore_fragments_GUESSFM_JUNE_13/frags_13_06_2017_0.01.csv'))
ichip.gf$analysis<-'ICHIP_GF'

## merge ICHIP and PMI

## first remove those traits that are in GUESSFM from PMI
keep<-setdiff(unique(ichip.pmi$Disease),unique(ichip.gf$Disease))
ichip.pmi<-subset(ichip.pmi,Disease %in% keep)

## need to remove some columns so we can join.

keep.cols<-intersect(names(ichip.pmi),names(ichip.gf))

ichip.all<-rbind(ichip.pmi[,keep.cols,with=FALSE],ichip.gf[,keep.cols,with=FALSE])

if(expressed.t.cells.only)
 ichip.all<-subset(ichip.all,!is.na(Log_FC_Gene_Expression_On_Activation))

ic.d<-subset(ichip.all,!COGS_Category %in% c('coding','promoter','noncoding','overall'))
gwas.d<-subset(gwas,!COGS_Category %in% c('coding','promoter','noncoding','overall'))
gwas.d<-subset(gwas.d,!is.na(Distance_between_Bait_Max_MPPI_Fragment))

## Total unique genes 
length(unique(gwas$Ensembl_GeneID))
## 505
## Genes per trait
table(gwas$Disease)
## CD CEL  MS PBC  RA SLE T1D  UC 
## 95  32  87  56 121 100  53 164 

## Genes per trait per category
table(gwas$COGS_Category,gwas$Disease)

## Number of genes in known region
table(gwas$Trait_Susceptibility_Region_Gene,gwas$Disease)
# CEL CRO  MS PBC  RA SLE T1D  UC
# FALSE  23  47  55  44  47  86  36 134
# TRUE    9  46  32  12  40  10  16  27

## Number of genes already causal candidate
table(gwas$Trait_Causal_Candidate_Gene,gwas$Disease)
# CEL CRO  MS PBC  RA SLE T1D  UC
# FALSE  25  59  69  48  65  88  43 134
# TRUE    7  34  18   8  22   8   9  27

## Mean number of genes per trait
nrow(gwas)/length(unique(gwas$Disease))

## number of protein coding genes skipped

load(file.path(DATA.DIR,'JAVIERRE_GWAS/RDATA/javierre_interactions.RData'))
int<-subset(int,biotype=='protein_coding')
int<-int[,.(ensg,name,biotype,strand,baitChr,baitStart,baitEnd,baitID)]
setkey(int,ensg)
int<-unique(int)

fstart<-with(gwas.d,ifelse(sign(Distance_between_Bait_Max_MPPI_Fragment)==-1,Max_MPPI_Fragment_Start,Bait_Start))
fstop<-with(gwas.d,ifelse(sign(Distance_between_Bait_Max_MPPI_Fragment)==-1,Bait_End,Max_MPPI_Fragment_End))

gwas.d.bait.gr<-with(gwas.d,GRanges(seqnames=Rle(Bait_Chr),ranges=IRanges(start=fstart,end=fstop),baitID=Bait_ID,oeID=Max_MPPI_Fragment_ID,name=Ensembl_GeneID,Disease=Disease,uid=1:nrow(gwas.d)))
int.gr<-with(int,GRanges(seqnames=Rle(baitChr),ranges=IRanges(start=baitStart,end=baitEnd),baitID=baitID,name=name,ensg=ensg))

foo<-mergeByOverlaps(gwas.d.bait.gr,int.gr)
blah<-split(foo,foo$uid)

skip.gene.count.gwas<-sapply(blah,function(t){
  ## remove those genes that lie in the actual fragments
  all.baits<-unique(c(t$baitID,t$oeID))
  nrow(t)-length(which(t$baitID.1 %in% all.baits))
})
skip.gene.count.gwas<-data.table(id=as.numeric(names(skip.gene.count.gwas)),skip.gene.count=skip.gene.count.gwas)
setkey(skip.gene.count.gwas,id)
gwas.d$uid<-1:nrow(gwas.d)
setkey(gwas.d,uid)
merge(gwas.d,skip.gene.count.gwas,by.x='uid',by.y='id',all.x=TRUE)
gwas.d<-merge(gwas.d,skip.gene.count.gwas,by.x='uid',by.y='id',all.x=TRUE)

gwas.d[,list(in.region=sum(Trait_Susceptibility_Region_Gene)),by=Disease]

median(gwas.d$skip.gene.count,na.rm = TRUE)

### next for iCHIP

## Total unique genes 
length(unique(ichip.all$Ensembl_GeneID))
## 256
## Genes per trait
table(ichip.all$Disease)
# ATD CEL CRO  MS NAR PBC  PS  RA T1D  UC 
# 15  40 114  76   1  50  37  25  56  50 

## Genes per trait per category
table(ichip.all$COGS_Category,ichip.all$Disease)

table(ichip.all$Trait_Susceptibility_Region_Gene,ichip.all$Disease)
# ATD CEL CRO MS NAR PBC PSO RA T1D UC
# FALSE  11  28  70 45   1  32  21  9  24 36
# TRUE    4  12  42 31   0  17  14 16  32 14

## Number of genes already causal candidate
table(ichip.all$Trait_Causal_Candidate_Gene,ichip.all$Disease)
# CEL CRO  MS PBC  RA SLE T1D  UC
# FALSE  25  59  69  48  65  88  43 134
# TRUE    7  34  18   8  22   8   9  27

## Mean number of genes per trait
nrow(ichip.all)/length(unique(ichip.all$Disease))

#ggplot(ichip.all,aes(x=Disease)) + geom_bar()

## number of protein coding genes skipped

fstart<-with(ic.d,ifelse(sign(Distance_between_Bait_Max_MPPI_Fragment)==-1,Max_MPPI_Fragment_Start,Bait_Start))
fstop<-with(ic.d,ifelse(sign(Distance_between_Bait_Max_MPPI_Fragment)==-1,Bait_End,Max_MPPI_Fragment_End))
## note that if we have bait and max_mppi on different chr then this causes NULL
fstart[is.na(fstart)]<-1
fstop[is.na(fstop)]<-1

ic.d.bait.gr<-with(ic.d,GRanges(seqnames=Rle(Bait_Chr),ranges=IRanges(start=fstart,end=fstop),baitID=Bait_ID,oeID=Max_MPPI_Fragment_ID,name=Ensembl_GeneID,Disease=Disease,uid=1:nrow(ic.d)))
int.gr<-with(int,GRanges(seqnames=Rle(baitChr),ranges=IRanges(start=baitStart,end=baitEnd),baitID=baitID,name=name,ensg=ensg))

foo<-mergeByOverlaps(ic.d.bait.gr,int.gr)
blah<-split(foo,foo$uid)

skip.gene.count.ic<-sapply(blah,function(t){
  ## remove those genes that lie in the actual fragments
  all.baits<-unique(c(t$baitID,t$oeID))
  nrow(t)-length(which(t$baitID.1 %in% all.baits))
})
skip.gene.count.ic<-data.table(id=as.numeric(names(skip.gene.count.ic)),skip.gene.count=skip.gene.count.ic)
setkey(skip.gene.count.ic,id)
ic.d$uid<-1:nrow(ic.d)
setkey(ic.d,uid)
merge(ic.d,skip.gene.count.ic,by.x='uid',by.y='id',all.x=TRUE)
ic.d<-merge(ic.d,skip.gene.count.ic,by.x='uid',by.y='id',all.x=TRUE)
median(ic.d$skip.gene.count,na.rm = TRUE)


big.list<-rbind(ichip.all,gwas[,intersect(names(gwas),names(ichip.all)),with=FALSE])
big.list$COGS_Overall_MPPI<-pmin(big.list$COGS_Overall_MPPI,1)
write.csv(big.list,file="/scratch/ob219/transfer/big_list_final_20_6_17_0.01.csv",row.names = FALSE)
## Total Genes for each analysis per disease
