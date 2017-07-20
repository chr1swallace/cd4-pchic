## set things up using chrisw libraries

#.libPaths("/home/chrisw/R/lib")


library(randomFunctions)
library(devtools)
library(data.table)
library(GenomicRanges)
library(magrittr)
#load_all("/home/chrisw/RP/annotSnpStats")
#load_all("/home/chrisw/RP/GUESSFM")
CD4CHIC.OUT<-'/scratch/cew54/cd4chic/'
CD4CHIC.DATA<-CD4CHIC.OUT

## the universe of interactions including promoter etc.
#(load("/home/oliver/JAVIERRE_ICHIP/RDATA/GUESSFM_intervals_tnact.RData"))
## coding SNPs
(load("/home/ob219/scratch/DATA/JAVIERRE_ICHIP/RDATA/javierre_tnact_csnps.by.ld.RData"))
## unique fragments
load("/home/ob219/scratch/DATA/JAVIERRE_ICHIP/RDATA/GUESSFM_intervals_tnact.RData")
all<-unlist(GUESS.ichip.grl)
all$uid<-paste(all$id,all$ld.id,sep=':')
all.u<-all[!duplicated(all$uid),]
## next filter so we only include those interactions that overlap
## densely mapped region - the marginal (mPPi) for those that do not
## is zero.

## use Chris' designations from running GUESSFM over the diseases
## for which we have genotyping information.

cw.regions<-list.files("/scratch/cew54/FM-GUESS")
cw.regions<-cw.regions[cw.regions!='finished']
cw.chr<-sub("(^[0-9]+)[pq].*","\\1",cw.regions)
cw.start<-as.numeric(sapply(strsplit(cw.regions,"-"),'[[',2))
cw.end<-as.numeric(sapply(strsplit(cw.regions,"-"),'[[',3))
	
cw.dr<-GRanges(seqnames=Rle(cw.chr),ranges=IRanges(start=cw.start,end=cw.end),fname=cw.regions)

all.u.id<-subsetByOverlaps(all.u,cw.dr)$id

GUESS.ichip.grl.f<-lapply(GUESS.ichip.grl,function(x){
	x[x$id %in% all.u.id,]
})

logsum <- function(x) {
    my.max <- max(x) ##take out the maximum value in log form)
    my.res <- my.max + log(sum(exp(x - my.max )))
    return(my.res)
}

splitByOverlaps<-function(gr,s.gr){
	ol<-as.matrix(findOverlaps(gr,s.gr))
	gr$split<-'none'
	gr[ol[,1],]$split<-s.gr[ol[,2],]$fname
	tmp<-split(gr,gr$split)
	tmp[names(tmp) != 'none']
}

cs.gr<-subsetByOverlaps(cs.gr,cw.dr)
cs.grl<-splitByOverlaps(cs.gr,cw.dr)
int.grl<-splitByOverlaps(all.u,cw.dr)

MPPI.overlap <- function(SM2, in.gr, out.gr) {
    ## extract PP from SM2 over SNPs in in.gr regions, excluding any in out.gr regions
    if(!length(in.gr))
        return(NULL)
    use.gr <- setdiff(in.gr,out.gr)
    if(!length(use.gr))
        return(NULL)
    o <- findOverlaps(use.gr,snps)
    #isnps <- o@subjectHits    
    isnps <- subjectHits(o)

    snps.overlap <- snps[isnps,]$snp.name
    lapply(SM2, function(sm) {
        wh <- which(sapply(sm@model.snps, function(ss) any(ss %in% snps.overlap)))
        if(!length(wh))
            return(0)
        logsum(sm@models[wh,"lPP"]) %>% exp()
    }) %>% unlist()
}

g.f<-"/scratch/cew54/FM-GUESS"
library(parallel)
out<-mclapply(seq_along(int.grl),function(i){
	region<-names(int.grl)[i]
	print(paste("Processing",region))
	ints<-int.grl[[i]]
	fpath<-file.path(g.f,region)
	#f.geno <- paste0("/stats/chrisw/FM-impute_qc/imputed-",region,".RData")
	f.geno <- paste0("/scratch/cew54/FM-impute_qc/imputed-",region,".RData")		
	f.sm <- file.path(fpath, "snpmod-fixed2.RData")
	if(!file.exists(f.sm))
		return(NA) 
	message("\n",region)
	chr <- sub("[pq].*","",basename(region))
	(load(f.sm))
	(load(f.geno))
	snps <- with(DATA@snps, GRanges(seqnames=Rle(chr),
	  ranges=IRanges(start=position,width=pmax(nchar(A1),nchar(A2)))))
	mcols(snps) <- data.frame(snp.name=rownames(DATA@snps),stringsAsFactors=FALSE)
	ints<-subsetByOverlaps(ints,snps)
	if(!is.null(cs.grl[[region]])){
		ol<-as.matrix(findOverlaps(snps,cs.grl[[region]]))
		## remove location of all coding SNPs
		## these should not contribute to fragment marginal PPi
		if(nrow(ol)>0)
			snps<-snps[-ol[,1],]
	}
	m<-mergeByOverlaps(snps,ints)
	m<-split(m$snp.name,m$id)
	names(SM2)<-sub("[.]","_",names(SM2))
	foo<-lapply(SM2, function(sm) {
		## build a fast lookup for all the models for a disease
		p<-sm@model.snps
		names(p)<-paste0('model',1:length(p),'_')
		p<-unlist(p)
		rp<-split(as.numeric(sub("model([0-9]+)[_].*","\\1",names(p))),p)
		lapply(m,function(sl){
        	wh<-unique(unlist(rp[sl]))
        	if(!length(wh))
            	return(0)
            logsum(sm@models[wh,"lPP"]) %>% exp()
    	}) %>% unlist()
    }) %>% unlist()
    nf<-names(foo)
    res<-data.table(do.call("rbind",strsplit(nf,'\\.')),mppi=foo)
    setnames(res,c('disease','frag.id','mppi'))
    res$GUESSFM.region<-region
	res	
},mc.cores=24)
names(out)<-names(int.grl)

(load("/home/ob219/scratch/DATA/JAVIERRE_ICHIP/RDATA/javierre_tnact_interactions.RData"))


f<-rbindlist(out[!is.na(out)])
f$frag.id<-as.numeric(f$frag.id)
setkey(f,frag.id)
setkey(int,oeID)


### this merge is fine for interactions but for fragments that are based on promoters
### we will miss the correct things - we need to add !


## this is the bit where we create dummy lines for promoter fragments
names(all.u)<-NULL
proms<-as.data.frame(GUESS.ichip.grl$promoter_only)
ify<-with(proms,data.table(ensg=ensg,baitChr=seqnames,baitStart=start,baitEnd=end,
baitID=id,baitName='.',oeChr=seqnames,oeStart=start,oeEnd=end,oeID=id,oeName='.',
dist=1,Total_CD4_Activated=FALSE,Total_CD4_NonActivated=FALSE))

xtra<-int[,.(ensg,biotype,strand,name)]
setkey(xtra,ensg)
xtra<-unique(xtra)
setkey(ify,ensg)
ify<-xtra[ify]
inty<-rbind(int,ify)


results<-merge(f,inty,by.x='frag.id',by.y='oeID')
v.high.score<-subset(results,mppi>0.9)
ic.stub<-c('AA','AS','ATD','CEL','CRO','IBD','IGE','JIA','MS','NAR','PBC','PSC','PSO','RA','SJO','SLE','SSC','T1D','UC','VIT')
trait.genes<-lapply(seq_along(ic.stub),function(i){
  t<-fread(paste0('https://www.immunobase.org/downloads/regions-files-archives/latest_1.11/Hs_GRCh37-',ic.stub[i],'-assoc_genesTAB'))
  t$disease<-ic.stub[i]
  t
})


	




tg<-do.call("rbind",trait.genes)
setnames(tg,sub(" ","_",names(tg)))
cc.genes<-unique(subset(tg,Cand_Gene==1)$Ensembl_ID)
in.region.genes<-unique(subset(tg,In_Region==1)$Ensembl_ID)

pos.reg<-unique(tg$Region_Pos)
disease_region<-lapply(split(tg$Region_Pos,tg$disease),function(pos.reg){
	pchr<-sub("([^:]+):([0-9]+)\\.\\.([0-9]+)","\\1",unique(pos.reg))
	pstart<-sub("([^:]+):([0-9]+)\\.\\.([0-9]+)","\\2",unique(pos.reg))
	pend<-sub("([^:]+):([0-9]+)\\.\\.([0-9]+)","\\3",unique(pos.reg))
	imb.reg<-GRanges(seqname=Rle(pchr),ranges=IRanges(start=as.numeric(pstart),end=as.numeric(pend)))
	unique(subsetByOverlaps(all,imb.reg)$id)
})

disease_gene<-lapply(split(tg,tg$disease),function(t){
	unique(subset(t,Cand_Gene==1)$Ensembl_ID)
})

disease_region_gene<-lapply(split(tg,tg$disease),function(t){
	unique(subset(t,In_Region==1)$Ensembl_ID)
})

results<-subset(results,disease %in% c('ICOELIAC','JIA','RA','T1D','MS','GRAVES'))
results[results$disease=='ICOELIAC',]$disease<-'CEL'
results[results$disease=='GRAVES',]$disease<-'ATD'



g<-results


hgenes<-fread("/home/ob219/scratch/DATA/JAVIERRE_ICHIP/out/GUESSFM/gene_prioritisation_0.01.csv")
hgenes<-subset(hgenes,!disease %in% c('MS','JIA'))

### this needs to be matched by disease otherwise we include silly fragments in latter analysis

rf<-hgenes[hgenes$node!='coding',]
rf[rf$disease=='ICOELIAC',]$disease<-'CEL'
rf[rf$disease=='GRAVES',]$disease<-'ATD'
disease.by.gene<-split(rf$ensg,rf$disease)

#tf<-subset(g,g$ensg %in% hgenes[hgenes$node!='coding',]$ensg)
g.lu<-split(rf$ensg,rf$disease)
tf<-rbindlist(lapply(split(g,g$disease),function(a){
	d<-unique(a$disease)
	print(d)
	subset(a,ensg %in% g.lu[[d]])
}))


## next add in closest loci

boho<-do.call("rbind",lapply(split(tf,tf$disease),function(d){
	do.call("rbind",lapply(split(d,d$ensg),function(b){
		b[which.max(b$mppi),]
	}))
}))


### add region and snp annotations

mesh.lu<-list(AA='Alopecia Areata',
AS='Spondylitis, Ankylosing',
ATD='Thyroiditis, Autoimmune',
CEL='Celiac Disease',
CRO='Crohn Disease',
IBD='Inflammatory Bowel Disease',
IGE='',
JIA='Arthritis, Juvenile Rheumatoid',
MS='Multiple Sclerosis',
NAR='Narcolepsy',
PBC='Liver Cirrhosis, Biliary',
PSC='Cholangitis, Primary Sclerosing',
PSO='Psoriasis',
RA='Arthritis, Rheumatoid',
SJO="Sjogren's Syndrome",
SLE='Lupus Erythematosus, Systemic',
SSC='Scleroderma, Systemic',
T1D='Diabetes Mellitus, Type 1',
UC='Colitis, Ulcerative',
VIT='Vitiligo')

trait.snps<-lapply(seq_along(ic.stub),function(i){
  t<-fread(paste0('https://www.immunobase.org/downloads/regions-files-archives/latest_1.11/Hs_GRCh37-',ic.stub[i],'-assoc_variantsTAB'))
  mesh<-mesh.lu[[ic.stub[i]]]
  print(mesh)
  t<-t[t$Disease==mesh,]
  colnames(t)<-make.names(colnames(t))
  t$rchr<-sub("([^:]+):([0-9]+)\\.\\.([0-9]+)","\\1",t$Coords.)
  t$rstart<-sub("([^:]+):([0-9]+)\\.\\.([0-9]+)","\\2",t$Coords.)
  t$rend<-sub("([^:]+):([0-9]+)\\.\\.([0-9]+)","\\3",t$Coords.)
  t$disease<-ic.stub[i]
  t[order(t$P.Value,decreasing=FALSE),]
})
names(trait.snps)<-ic.stub

getClosestSNP<-function(snp.df,frag.chr,frag.start,frag.end){
	fm<-subset(snp.df,rchr==frag.chr)
	if(nrow(fm)==0)
		return(data.table(Region='NONE',Coords.='NONE',Rs.Id='NONE',Alleles='NONE',P.Value=NA,Odds.Ratio=NA,Position=NA,lead.snp.frag.dist=NA,lead.snp.in.frag=FALSE))
	frag.size<-frag.end-frag.start
	mid<-((frag.size)/2)+frag.start
	idx<-which.min(abs(fm$Position-mid))
	fm<-fm[idx,.(Region,Coords.,Rs.Id,Alleles,P.Value,Odds.Ratio,Position)]
	fm$lead.snp.frag.dist<-abs(mid-fm$Position)
	if(fm$lead.snp.frag.dist>5e6)
		return(data.table(Region='NONE',Coords.='NONE',Rs.Id='NONE',Alleles='NONE',P.Value=NA,Odds.Ratio=NA,Position=NA,lead.snp.frag.dist=NA,lead.snp.in.frag=FALSE))
	fm$lead.snp.in.frag<-fm$lead.snp.frag.dist<(frag.size/2)
	fm
		
}

reg<-rbindlist(lapply(1:nrow(boho),function(i){
	b<-boho[i,]
	getClosestSNP(trait.snps[[b$disease]],b$oeChr,b$oeStart,b$oeEnd)
}))

out<-cbind(boho,reg)

## remove MS and JIA
#out<-subset(out,!disease %in% c('MS','JIA'))

## add in the number of GUESSFM snps.
(load("/scratch/cew54/FM/nsnps-fixed2.RData"))
DATA<-as.data.table(DATA)
nexp<-DATA[,.(nexp=sum(nsnp*pp)),by=c("trait","region")]
## normal subsetting does not work not sure why
idx<-which(nexp$trait=='ICOELIAC')
nexp[idx,]$trait<-'CEL'
idx<-which(nexp$trait=='GRAVES')
nexp[idx,]$trait<-'ATD'
nexp<-subset(nexp,trait %in% g$disease)
nexp$jid<-paste(nexp$trait,nexp$region,sep=':')
## merge in 
out$jid<-paste(out$disease,out$GUESSFM.region,sep=':')
gt<-merge(out,nexp[,.(jid,nexp)],by.x='jid',by.y='jid',all.x=TRUE)

## chris has a bunch of helper functions

source("/home/cew54/Projects/cd4chic/activation-analyses/R/common.R")
if(!exists("hind") || !("GRanges" %in% class(hind)))
    hind <- get.hind()

erna<-get.erna()
#erna<-subset(erna,type=='regulatory')
erna<-subset(erna,type %in% c("regulatory","intergenic.reg") & expr)
merna<-erna[,.(name,logFC.erna,FDR.erna)]
merna<-merna[,.SD[which.min(FDR.erna),],by=name]
setnames(merna,c('name','erna.logFC','erna.adj.P.Val'))
merna$name<-as.integer(merna$name)
gt<-merge(gt,merna,by.x='frag.id',by.y='name',all.x=TRUE)

dchic<-get.diffchic()
dchic$uid<-paste(dchic$baitID,dchic$oeID,sep=':')
gt$uid<-paste(gt$baitID,gt$frag.id,sep=':')
dchic<-dchic[,.(uid,logFC,FDR)]
setnames(dchic,c('uid','chic.logFC','chic.FDR'))

gt<-merge(gt,dchic,by.x='uid',by.y='uid',all.x=TRUE)
gt$hid<-paste(gt$disease,gt$ensg,sep=':')
gtx<-gt[,names(gt)[!names(gt) %in% names(hgenes)],with=FALSE]
hgenes$hid<-paste(hgenes$disease,hgenes$ensg,sep=':')

final<-merge(hgenes,gtx,by.x='hid',by.y='hid',all.x=TRUE)

rnaseq<-get.rnaseq()

final<-merge(final,rnaseq[,.(id,logFC,adj.P.Val)],by.x='ensg',by.y='id',all.x=TRUE)

final<-subset(final,biotype=='protein_coding')

g<-lapply(split(final,final$disease),function(d){
	dis<-unique(d$disease)
	d$frag.overlaps.reg<-d$frag.id %in% disease_region[[dis]]
	d$disease.gene<-d$ensg %in% disease_gene[[dis]]
	d$disease.region.gene<-d$ensg %in% disease_region_gene[[dis]]
	d
})

final<-rbindlist(g)

## remove columns that were used to build table

final$hid<-NULL
final$jid<-NULL
final$uid<-NULL

scols<-fread("/scratch/ob219/DATA/JAVIERRE_ICHIP/support/supp_table_output_columns.csv")
setnames(final,scols$raw.colname,scols$pretty.colname)



#write.csv(final,file="/home/oliver/JAVIERRE_ICHIP/out/highscore_fragments_GUESSFM/frags_13_06_2016.csv",row.names=FALSE)
write.csv(final,file="/scratch/ob219/DATA/JAVIERRE_ICHIP/out/highscore_fragments/frags_0.01.csv",row.names=FALSE)

