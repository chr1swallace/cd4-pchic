## cheplev/martin and javierre overlap

## load in cheplev data and see what it's all about
## the data comes from supp table 1

library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(reshape2)
library(magrittr)
options(stringAsFactors=FALSE)

## alter this and alter chicago thresh to use
thresh<-5
setwd(CD4CHIC.EXTDATA)

source(file.path(CD4CHIC.ROOT,"activation-analyses/R/common.R"))
prom.ass<-fread('HindIII_baits_e75.bed')
setnames(prom.ass,c('chr','start','end','frag.id','gene_details'))

chepfile <- "http://www.nature.com/cr/journal/v22/n3/extref/cr201215x1.xlsx"
system(paste("wget",chepfile))
library(readxl)
chep <- basename(chepfile) %>%
read_excel(.,skip=1) %>% as.data.frame() %>% as.data.table()
chep$line<-1:nrow(chep)
basename(chepfile) %>% unlink()

chep[,R1.start:=as.integer(R1.start)]
chep[,R2.start:=as.integer(R2.start)]
chep[,R1.end:=as.integer(R1.end)]
chep[,R2.end:=as.integer(R2.end)]
chep[,N.12:=as.integer(N.12)]
chep[,N.1:=as.integer(N.1)]
chep[,N.2:=as.integer(N.2)]

chepr1<-with(chep,GRanges(seqnames=Rle(R1.chrom),ranges=IRanges(start=R1.start,end=R1.end),id=line))
chepr2<-with(chep,GRanges(seqnames=Rle(R2.chrom),ranges=IRanges(start=R2.start,end=R2.end),id=line))

## lift these over to build 37
f <- "http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz"
system(paste("wget",f))
system(paste("gunzip",basename(f)))
f <- sub(".gz","",f)
c<-basename(f) %>% import.chain() 
chepr1.37<-unlist(liftOver(chepr1,c))  
chepr2.37<-unlist(liftOver(chepr2,c))
basename(f) %>% unlink()

chep.37<-GRangesList(r1=chepr1.37,r2=chepr2.37)
seqlevels(chep.37)<-sub("^chr","",seqlevels(chep.37))


## mean frag size of 3.7kb

## how do hindIII fragments map to these fragments ?

h3.file<-'Digest_Human_HindIII.bed'
h3<-fread(h3.file,header=FALSE)
setnames(h3,c('chr','start','end','id'))
## get rid of y chromosome
h3<-subset(h3,chr!='Y')
h.gr<-with(h3,GRanges(seqnames=Rle(chr),ranges=IRanges(start=start,end=end),id=id))

## which of these has a promoter associated with it ?
## note the cheplev paper uses a 5KB window around the middle of a fragment to assign promoters.

mapChepToHindIII<-function(chep.gr,h.gr,geneDetails){
	#add in line number to hind3 frag
	tmp<-mergeByOverlaps(h.gr,chep.gr)
	tmp[['h.gr']]$line<-tmp$id.1
	tmp<-tmp[['h.gr']]
	tmp<-data.table(as.data.frame(tmp[!duplicated(paste(tmp$id,tmp$line,sep=":"))]))
	##add in gene details
	tmp<-merge(tmp,geneDetails,by.x="id",by.y="frag.id",all.x=TRUE)
	tmp<-tmp[,.(line,seqnames,start.x,end.x,id,gene_details)]
	setnames(tmp,c('line','chr','start','end','frag.id','gene_details'))
}

library(parallel)
options(mc.cores=2)
m<-mclapply(seq_along(chep.37),function(i){
	gr<-chep.37[[i]]
	mapChepToHindIII(gr,h.gr,prom.ass)
})

setnames(m[[1]],paste('r1',names(m[[1]]),sep="_"))
setnames(m[[2]],paste('r2',names(m[[2]]),sep="_"))

m.chep<-merge(m[[1]],m[[2]],by.x="r1_line",by.y="r2_line",allow.cartesian=TRUE)
m.chep$uid<-with(m.chep,paste(r1_frag.id,r2_frag.id,sep=":"))
m.chep$ruid<-with(m.chep,paste(r2_frag.id,r1_frag.id,sep=":"))

## load in jav data
    require(data.table)#merged_samples_12Apr2015_full_denorm_bait2baits_e75.tab
##    int <- fread(file.path(CD4CHIC.DATA,"merged_samples_12Apr2015_full_denorm_bait2baits_e75.tab")) # %>% as.data.frame()
    int <- fread(file.path(CD4CHIC.DATA,"merged_samples_12Apr2015_full.txt"))
    int.genes <- get.b2gene()
    setkey(int,baitID)
    setkey(int.genes,baitID)
    int <- merge(int,int.genes[,.(id,gene,biotype,strand,baitID)],allow.cartesian=TRUE)
    ## ## subset to protein coding only
    ## int <- int[biotype=="protein_coding",]
    ## int[ oeID==120852 & baitID==120843, .(baitChr,baitStart,baitID,oeChr,oeStart,oeID,id,gene,Total_CD4_Activated,Total_CD4_NonActivated)]
    int <- int[,b2b:= oeID %in% int$baitID]
    int[,Background:=apply(int[,.(Erythroblasts,Megakaryocytes)],1,max)]
threshold <- 5
        int <- int[pmax(Background,Total_CD4_Activated,Total_CD4_NonActivated)>=threshold, ]
    int <- int[, c("baitChr","oeChr") := list(paste0("chr",baitChr), paste0("chr",oeChr)) ]
    int <- int[, c("baitLength","oeLength"):= list(abs(baitStart-baitEnd), abs(oeStart-oeEnd)) ]
    int <- int[,.(baitChr,baitStart,baitID,oeChr,oeStart,oeID,id,gene,Total_CD4_Activated,Total_CD4_NonActivated,Erythroblasts,Megakaryocytes,biotype,baitLength,oeLength,b2b)]
int$uid<-with(int,paste(baitID,oeID,sep=":"))
int$ruid<-with(int,paste(oeID,baitID,sep=":"))

## chromatin data
## (load("/ipswich/data/chrisw/cd4chic-output/peaks/hind-our-marks.RData"))

## act <- pmax(hind$H3K4me3_PoolQTW_Act_2, hind$H3K4me3_PoolRG_Act_2) 
## non <- pmax(hind$H3K4me3_PoolQTW_NAct_2, hind$H3K4me3_PoolRG_NAct_2) 
## median(act,na.rm=TRUE)
## median(non,na.rm=TRUE)
## hind$enh.act <- act
## hind$enh.non <- non
## hind <- as.data.table(as.data.frame(hind))
## hind$name <- as.integer(hind$name)
## setkey(hind,name)
## setkey(int,oeID)
## int <- merge(int,hind[,.(name,enh.act,enh.non)],by.x="oeID",by.y="name")

int[,mega.nocd4:=ifelse(pmax(Total_CD4_Activated,Total_CD4_NonActivated)>5, 0, Megakaryocytes)]
int[,ery.nocd4:=ifelse(pmax(Total_CD4_Activated,Total_CD4_NonActivated)>5, 0, Erythroblasts)]
comp.jav<-list(
	## all.jav=header[16:32],
	## cd4pos=c('Naive_CD4','Total_CD4_MF','Total_CD4_Activated','Total_CD4_NonActivated'),
	## b=c('Total_B','Naive_B'),
	tact="Total_CD4_Activated",
	tnact="Total_CD4_NonActivated",
	erythroblasts="Erythroblasts",
	megakaryocytes="Megakaryocytes",
	ery.nocd4="ery.nocd4",
	mega.nocd4="mega.nocd4"
	)
	
############################ COMPARE CHEPELEV TO JAV ########################################################

compare_chep<-function(chep,pm,comp.pm,name1="all.chep",name2,strat.label="NONE"){
	if(length(comp.pm)>1){
		t2<-which(rowSums(pm[,comp.pm,with=FALSE]>=thresh)>0)
	}else{
		t2<-which(pm[[comp.pm]]>=thresh)
	}
	t1uid<-unique(chep$uid)
	t2uid<-unique(c(pm[t2,]$uid,pm[t2,]$ruid))
	perc<-round((sum(t1uid %in% t2uid)/length(t1uid))*100,digits=1)
	data.frame(n1=name1,n2=name2,overlap=perc,stratification=strat.label)
}

## we have no chance to find e2e
m.chep[,e2e:=is.na(r1_gene_details) & is.na(r2_gene_details)]
## also flag promoter to promoter?
m.chep[,p2p:=!is.na(r1_gene_details) & !is.na(r2_gene_details)]

options(mc.cores=length(comp.jav))
chep_vs_jav <- do.call("rbind",mclapply(seq_along(comp.jav),function(j){
		message(paste("Comparing all.chep with",names(comp.jav)[j]))
		compare_chep(m.chep[! m.chep$e2e,],int,comp.jav[[j]],'all.chep',names(comp.jav)[j])
	}))

## stratified analysis

s<-c(1e3,1e4,5e4,1e5,5e5,1e6,5e6,1e7,1e8,1e9)

int[, dist:= abs( (baitStart + baitLength/2) - (oeStart + oeLength/2) )]
int.strat<-split(1:nrow(int),cut(int$dist,s))
m.chep$dist<-with(m.chep,ceiling(abs((((r1_end-r1_start)/2) + r1_start) - (((r2_end-r2_start)/2) + r2_start))))
## for this analysis remove e2e
chep.f<-subset(m.chep,!e2e)

chep.strat<-split(1:nrow(chep.f),cut(chep.f$dist,s))

chep_vs_jav_strat<-do.call("rbind",mclapply(seq_along(chep.strat),function(z){
	intz<-int[int.strat[[z]],]
	chepz<-chep.f[chep.strat[[z]],]
	slabel<-names(chep.strat)[z]
	do.call("rbind",lapply(seq_along(comp.jav),function(j){
		message(paste("Comparing all.chep with",names(comp.jav)[j]),slabel)
		compare_chep(chepz,intz,comp.jav[[j]],'all.chep',names(comp.jav)[j],slabel)
	}))
}))


pretty.chep_vs_jav_strat<-melt(chep_vs_jav_strat,id=c('n1','n2','stratification'))
pretty.chep_vs_jav_strat<-recast(pretty.chep_vs_jav_strat,n1+n2~stratification)
pretty.chep_vs_jav_strat


## plot these

library(ggplot2)
library(magrittr)
df <- chep.f[,.(dist),drop=FALSE] %>% as.data.frame()
df$study <- "Chepelev"
df2 <- int[Total_CD4_NonActivated>=5,.(dist),drop=FALSE] %>% as.data.frame()
df2$study <- "Javierre"
df <- rbind(df,df2)
df$cdist <- cut(df$dist,s)
xlabs <- c("<10\nkb","10-50\nkb","50-100\nkb","100-500\nkb","500 - 1\n kb  mb","1-5\nmb","5-10\nmb","10-100\nmb","100 - 1\n mb  Gb",">1\nGb")

dfn <- as.data.table(df)
dfn <- dfn[,N:=length(dist),by=c("study","cdist")] %>% unique(.,by=c("study","cdist"))
dfn[,dist:=NULL]
#dfn[,N:=min(N),by="cdist"]

df2 <- melt(chep_vs_jav_strat,id=c('n1','n2','stratification'))
df2$study <- "Chepelev"
df2$n2 %<>% sub("tact","CD4-act",.)
df2$n2 %<>% sub("tnact","CD4-nonact",.)
df2$n1 <- df2$n2
df2 <- merge(df2,dfn,by.x=c("study","stratification"),by.y=c("study","cdist"))
df2 <- as.data.table(df2)

dff <- df2 #rbind(df2,df3,df4)

library(cowplot)

p3 <- ggplot(subset(dff,n1 %in% c("CD4-act","CD4-nonact","ery.nocd4","mega.nocd4") & study=="Chepelev"),
       aes(x=stratification,ymin=0,y=value,ymax=value/100,fill=n1)) +
    geom_bar(stat="identity",position="dodge",col="grey") +
theme(legend.position="bottom") +
labs(x="Interaction distance",y="Overlap (%)") +
scale_fill_brewer("Javierre cell",palette="Set2",labels=c("CD4 (a)", "CD4 (n)", "Erythroblasts, not CD4", "Megakaryocytes, not CD4")) +
background_grid() +
scale_x_discrete(labels=xlabs) +
scale_y_continuous(breaks=c(0,10,20,30,40,50,60))

p2 <- ggplot(subset(dff,n1 %in% c("CD4-act","CD4-nonact","ery.nocd4","mega.nocd4") & study=="Chepelev"),
       aes(x=stratification,y=value*N/100,fill=n1)) +
geom_bar(mapping=aes(y=N),stat="identity",fill="black",data=dff[study=="Chepelev" & !duplicated(paste(stratification,study)),]) +
    geom_bar(stat="identity",position="dodge",col="grey") +
theme(legend.position="none") +
labs(x="Interaction distance",y="Count") +
scale_fill_brewer("Javierre cell",palette="Set2",labels=c("CD4 (a)", "CD4 (n)", "Erythroblasts, not CD4", "Megakaryocytes, not CD4")) +
scale_x_discrete(labels=xlabs) +
scale_y_continuous(breaks=c(0,5000,1e+4),labels=c("0","5,000","10,000"))

p1 <- ggplot(dfn[study=="Javierre",],aes(x=cdist,y=N)) + geom_bar(stat="identity",fill="black") + labs(x="Interaction distance",y="Count")+ scale_x_discrete(labels=xlabs)  + scale_y_continuous(breaks=c(0,2e+5,4e+5,6e+5),labels=c("0","200,000","400,000","600,000"))

plot_grid(p1,p2,p3,nrow=3,labels=c("a","b","c"),align="v")

f <- file.path(CD4CHIC.OUT,"paper/chepelev-comp.pdf")
ggsave(f,height=10,width=10)
system(paste("display",f))



