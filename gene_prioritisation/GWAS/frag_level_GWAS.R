### Stuff for Roman

library(data.table)
library(GenomicRanges)

source("/home/ob219/git/CHIGP/R/common.R")
args<-getArgs(verbose=TRUE)
CD4CHIC.DATA<-'/scratch/cew54/cd4chic/'
CD4CHIC.OUT<-'/scratch/cew54/cd4chic-output/'



GS.THRESHOLD <- 0.01
DATA.DIR <- '/home/ob219/scratch/DATA/JAVIERRE_GWAS/'
disease<-'SLE'
#disease<-args[['disease']]
#odir<-args[['out_dir']]
odir<-file.path(DATA.DIR,'out/highscore_fragments_JUNE_13')
tissues<-c('Total_CD4_NonActivated','Total_CD4_Activated')
#tissues<-paste0(tissues,'_interactions_only_gene_score')

#pri<-fread(file.path(DATA.DIR,'out/tnact_hierarchical_geneScore',paste(disease,'pmi_prioritised','tab',sep='.')))
pri<-fread(file.path(DATA.DIR,'out/hierarchical_geneScore_tnact_0.01',paste(disease,'pmi_prioritised','tab',sep='.')))
#pri<-t[,c('ensg','disease','name','strand','biotype',tissues),with=FALSE]

## these are precomputed data objects
(load(file.path(DATA.DIR,'RDATA','javierre_interactions.RData')))
(load(file.path(DATA.DIR,'RDATA','javierre_frags.by.ld.RData')))
(load(file.path(DATA.DIR,'RDATA','javierre_csnps.by.ld.RData')))
cs.lu<-paste(seqnames(cs.gr),start(cs.gr),sep=':')
 

t<-fread(file.path(DATA.DIR,'out/pmi',paste(disease,'pmi',sep='.')))
t$disease<-disease
## remove coding SNPs
t<-subset(t,!paste(chr,end,sep=':') %in% cs.lu)
pmi<-with(t,GRanges(seqnames=Rle(chr),ranges=IRanges(start=start+1,end=end),ppi=ppi,disease=disease))



#idx<-unique(do.call('c',lapply(tissues,function(x){
#  which(pri[[x]]>GS.THRESHOLD)
#})))

#pri<-pri[idx,]
##table(pri$disease)



## load in pmi data and remove any coding SNPs as not relevant to this analysis.



library(parallel)
all.res<-lapply(unique(pri[node!='coding',]$ensg),function(g){
  print(paste('Processing',g))
  int.g<-subset(int,ensg==g)
  frags.g<-subset(frags.gr,ensg==g & type %in% c('interaction','promoter'))
  spri<-subset(pri,ensg==g)
  p<-pmi
  #abvt<-subset(pri,ensg==g & disease == disease )[,tissues,with=FALSE]>GS.THRESHOLD
  #st<-sub('_interactions_only_gene_score','',colnames(abvt)[which(abvt)])
  fb<-do.call('rbind',lapply(tissues,function(ti){
      fr<-frags.g[frags.g$id %in% int.g[int.g[[ti]],]$oeID,]
      fr<-c(fr,subset(frags.g,type=='promoter'))
      #do.call('rbind',lapply(pmi,function(p){
      d<-disease
      gsp<-mergeByOverlaps(fr,p)
      gsp<-with(gsp,data.table(ld.id=ld.id,id=id,ppi=ppi,type=type))
      frag.score<-gsp[,list(ld.frag.score=sum(ppi),type=unique(type)),by=c('id','ld.id')]
      frag.score<-frag.score[,list(frag.score=1-prod(1-ld.frag.score),type=unique(type)),by=id]
      frag.score$disease<-disease
      frag.score$gene<-g
      frag.score$tissue<-ti

      frag.score
  }))
  ### this should return the max frag out of all
  fb[,.SD[which.max(frag.score)],]
})

hiscore.frags<-do.call("rbind",all.res)
setkey(hiscore.frags,id)
hsf<-unique(hiscore.frags)
hsf<-hsf[,.(id,frag.score,disease)]
setnames(hsf,c('frag.id','ppi','disease'))

ify<-with(as.data.frame(subset(frags.gr,type=='promoter')),data.table(ensg=ensg,baitChr=seqnames,baitStart=start,baitEnd=end,
baitID=id,baitName='.',oeChr=seqnames,oeStart=start,oeEnd=end,oeID=id,oeName='.',
dist=1,Total_CD4_Activated=FALSE,Total_CD4_NonActivated=FALSE))

## format int

tint<-int[Total_CD4_Activated==TRUE | Total_CD4_NonActivated==TRUE,.(ensg,biotype,strand,name,baitChr,baitStart,baitEnd,baitID,baitName,oeChr,oeStart,oeEnd,oeID,oeName,dist,Total_CD4_Activated,Total_CD4_NonActivated)]

xtra<-tint[,.(ensg,biotype,strand,name)]
setkey(xtra,ensg)
xtra<-unique(xtra)
setkey(ify,ensg)
ify<-xtra[ify]
inty<-rbind(tint,ify)

results<-merge(hsf,inty,by.x='frag.id',by.y='oeID')

ic.stub<-c('AA','AS','ATD','CEL','CRO','IBD','IGE','JIA','MS','NAR','PBC','PSC','PSO','RA','SJO','SLE','SSC','T1D','UC','VIT')
#trait.genes<-lapply(seq_along(ic.stub),function(i){
#  t<-fread(paste0('http://www.immunobase.org/regions/htdocs/downloads/Hs_GRCh37-',ic.stub[i],'-assoc_genesTAB'))
#  t$disease<-ic.stub[i]
#  t
#})
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
	unique(subsetByOverlaps(frags.gr,imb.reg)$id)
})

disease_gene<-lapply(split(tg,tg$disease),function(t){
	unique(subset(t,Cand_Gene==1)$Ensembl_ID)
})

disease_region_gene<-lapply(split(tg,tg$disease),function(t){
	unique(subset(t,In_Region==1)$Ensembl_ID)
})

results[results$disease=='CD_IMB',]$disease<-'CRO'
results[results$disease=='RA',]$disease<-'RA_EYRE'
results[results$disease=='RA_OKADA_IMB',]$disease<-'RA'

## g<-lapply(split(results,results$disease),function(d){
## 	dis<-unique(d$disease)
## 	d$frag.overlaps.reg<-d$frag.id %in% disease_region[[dis]]
## 	d$disease.gene<-d$ensg %in% disease_gene[[dis]]
## 	d$disease.region.gene<-d$ensg %in% disease_region_gene[[dis]]
## 	d
## })
## 
## g<-rbindlist(g)
g<-results
hgenes<-pri
rf<-hgenes[hgenes$node!='coding',]
rf[rf$disease=='CD_IMB',]$disease<-'CRO'
rf[rf$disease=='RA',]$disease<-'RA_EYRE'
rf[rf$disease=='RA_OKADA_IMB',]$disease<-'RA'
hgenes[hgenes$disease=='CD_IMB',]$disease<-'CRO'
hgenes[hgenes$disease=='RA',]$disease<-'RA_EYRE'
hgenes[hgenes$disease=='RA_OKADA_IMB',]$disease<-'RA'
disease.by.gene<-split(rf$ensg,rf$disease)

g.lu<-split(rf$ensg,rf$disease)
tf<-rbindlist(lapply(split(g,g$disease),function(a){
	d<-unique(a$disease)
	print(d)
	subset(a,ensg %in% g.lu[[d]])
}))

boho<-do.call("rbind",lapply(split(tf,tf$disease),function(d){
	do.call("rbind",lapply(split(d,d$ensg),function(b){
		b[which.max(b$ppi),]
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
  #t<-fread(paste0('http://www.immunobase.org/regions/htdocs/downloads/Hs_GRCh37-',ic.stub[i],'-assoc_variantsTAB'))
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
	if(is.null(snp.df))
		return(data.table(Region='NONE',Coords.='NONE',Rs.Id='NONE',Alleles='NONE',P.Value=NA,Odds.Ratio=NA,Position=NA,lead.snp.frag.dist=NA,lead.snp.in.frag=FALSE))
		
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
gt<-out

#source("/home/oliver/GIT_REPOS/cd4chic/activation-analyses/R/common.R")

source("/home/cew54/Projects/cd4chic/activation-analyses/R/common.R")
if(!exists("hind") || !("GRanges" %in% class(hind)))
    hind <- get.hind()
#source("/home/cew54/Projects/cd4chic/activation-analyses/R/eRNA-common.R")
## override get.erna found in this file
#get.erna<-function() {
 ## this loads into dt
# (load("/scratch/wallace/old-new-erna.RData"))
# erna <- fread(file.path(CD4CHIC.DATA,"rna-with-regrna-diff-v2.csv"))
# setnames(erna,c("logFC","adj.P.Val","1_act","2_act","1_non","2_non"),
#            c("logFC.erna","FDR.erna","act_1","act_2","non_1","non_2"))
# ## we don't trust the FDR  from arcadio so add in the values from dt
# erna$FDR.erna<-dt$adj.P.Val
# erna[,non.cpm1:=1e6 * non_1/sum(non_1)]
# erna[,non.cpm2:=1e6 * non_2/sum(non_2)]
# erna[,act.cpm1:=1e6 * act_1/sum(act_1)]
# erna[,act.cpm2:=1e6 * act_2/sum(act_2)]
# erna[,expr:= (as.numeric(non.cpm1>=0.4) + as.numeric(non.cpm2>=0.4) +
#               as.numeric(act.cpm1>=0.4) + as.numeric(act.cpm2>=0.4)) > 1 ]
# erna[,non.expr:=pmin(non.cpm1,non.cpm2)>=0.4]
# erna[,act.expr:=pmin(act.cpm1,act.cpm2)>=0.4]
# #erna <- erna[type %in% c("regulatory","lincRNA","protein_coding","pseudogene"), ]
# erna <- erna[type %in% c("regulatory"), ]
# intergenic.ids <- scan(file.path(CD4CHIC.DATA, "distant-erna.csv"),what="")
# erna[,intergenic:=id %in% intergenic.ids]
#erna[id %in% intergenic.ids,type:="intergenic.reg"]
#message("expressed in non-activated")
#with(erna,table(non.expr,type)) # 3897, want 3897 - agree
#message("expressed overall")
#with(erna,table(non.expr | act.expr,type)) # 5898, want 6133
#with(erna,table(expr,type)) # 6147, want 6133
#
# ## add hindIII ids
#  if(!exists("hind") || !("GRanges" %in% class(hind)))
#    hind <- get.hind()
#
#  erna <- add.hind(erna)
# erna[,name:=as.integer(name)]
# return(erna)
#}


#erna<-get.erna()
#merna<-erna[,.(name,logFC,adj.P.Val)]
#merna<-merna[,.SD[which.min(adj.P.Val),],by=name]
#setnames(merna,c('name','erna.logFC','erna.adj.P.Val'))
#merna$name<-as.integer(merna$name)
#gt<-merge(gt,merna,by.x='frag.id',by.y='name',all.x=TRUE)

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


### based on the gene so should go last

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

#scols<-fread("/home/oliver/JAVIERRE_ICHIP/support/supp_table_output_columns.csv")
scols<-fread("/scratch/ob219/DATA/JAVIERRE_ICHIP/support/supp_table_output_columns.csv")
scols[scols$raw.colname=='mppi',]$pretty.colname<-'Fragment_PPI'
scols[scols$raw.colname=='mppi',]$raw.colname<-'ppi'
scols<-scols[scols$raw.colname %in% names(final),]
setcolorder(final,scols$raw.colname)
setnames(final,scols$pretty.colname)


write.table(final,file=file.path(odir,paste(disease,'13_06_17_gs_0.01.csv',sep='_')),sep=",",quote=FALSE,row.names=FALSE)  


 library(data.table)
 AI<-c('RA_OKADA_IMB','SLE')
 catres<-rbindlist(lapply(AI,function(a){
 	#fname<-paste0(odir,'/',a,'_13_06_17_gs_0.5.csv')
 	fname<-paste0(odir,'/',a,'_13_06_17_gs_0.01.csv')
 	fread(fname)  
 }))
 
 #write.table(catres,file=file.path(odir,paste('Autoimmune','13_06_17_gs_0.5.csv',sep='_')),sep=",",quote=FALSE,row.names=FALSE)  
 write.table(catres,file=file.path(odir,paste('Autoimmune','13_06_17_gs_0.01.csv',sep='_')),sep=",",quote=FALSE,row.names=FALSE)  
