library(data.table)
library(reshape2)
library(yaml)
library(data.tree)


DATA.DIR <- '/home/ob219/scratch/DATA/JAVIERRE_ICHIP/'

#all.thresh<-0.5
all.thresh<-0.01

args <- list(
  out_dir = file.path(DATA.DIR,'out/hierachical_geneScore/'),
  sets = file.path(DATA.DIR,'support/tnact_tree.yaml'),
  decompose = 1,
 # ppi.thresh = 0.5,
  ppi.thresh = 0.01,
  BF.thresh = 3,
  geneName = 'PRKCQ'
)

NON_ZERO<-10e-50 ## software does not cope well with zeros


## FUNCTIONS 

getTreeYAML<-function(yfile){
  osList <- yaml.load_file(yfile)
  as.Node(osList)
}

GetSetsFromTree<-function(osNode){
  sets<-lapply(Traverse(osNode,filterFun=isNotLeaf),function(x) c(x$Get("name",filterFun=isLeaf),use.names = FALSE))
  names(sets)<-names(osNode$Get('levelName', filterFun = isNotLeaf))
  sets
}

DATA.DIR <- '/home/ob219/scratch/DATA/JAVIERRE_ICHIP/'

fil<-list.files(file.path(DATA.DIR,'out/GUESSFM_JUNE_13'),pattern="*.RData",full.names = TRUE)

dat<-lapply(fil,function(f) get(load(f)))
names(dat)<-sub("MPPI-([^\\.]+)\\.RData","\\1",basename(fil))

lu.df<-data.frame(guessn=c('ICOELIAC','COELIAC','RA','T1D','GRAVES','JIA','MS','RA.UK','RA.US','PBC'),
                  fname=paste0(c('CEL','CEL','RA','T1D','ATD','JIA','MS','RA','RA','PBC'),'_IC.pmi.tab'))
lu.df<-split(lu.df$fname,lu.df$guessn)


get_disease_gene_score<-function(d,dname){
  ## not sure that this is required but just in case
  setkeyv(d,c('ensg','ld.id','type'))
  ua<-unique(d)
  t<-ua[,list(ld.score=sum(get(dname))),by=c('ensg','ld.id')]
  #t<-subset(t,ld.score != 0)
  t[,list(gs=1-prod(1-ld.score)),by=ensg]
}

all.results<-do.call("rbind",lapply(seq_along(dat),function(i){
  d<-dat[[i]]
  print(paste("Processing",names(dat)[i]))
  na<-names(lu.df)
  o<-do.call("rbind",lapply(seq_along(na),function(j){
    print(paste("Disease",na[j]))
    r<-get_disease_gene_score(d,na[j])
    r$disease<-na[j]
    r
  }))
  o$ctype<-names(dat)[i]
  o
}))

all.results[all.results$gs==0,]$gs<-1e-13

## split by disease

ds<-split(all.results,all.results$disease)

tmp.stub<-fread(file.path('/home/ob219/scratch/DATA/JAVIERRE_GWAS/out/geneScore/RA_OKADA_IMB.pmi.tab'))[,(2:6),with=FALSE]
setkey(tmp.stub,ensg)
tmp.stub<-unique(tmp.stub)
foo<-lapply(seq_along(ds),function(i){
  print(i)
  test<-ds[[i]]
  test<-melt(test,id.vars = c('ensg','ctype'), measure.vars=c('gs'),variable.name = 'gene_score')
  test<-data.table(dcast(test,ensg~ctype,fill=0))
  test$disease<-names(ds)[i]
  setkey(test,ensg)
  ## filter to remove genes with no score
  test<-test[rowSums(test[,c(2:8),with=FALSE])!=0,]
  test<-merge(tmp.stub,test)
  setnames(test,c('ensg','name','biotype','strand','baitChr','coding_gene_score','Total_CD4_Activated_prey_only_gene_score',
                  'Total_CD4_NonActivated_prey_only_gene_score','interaction_gene_score','noncoding_gene_score','overall_gene_score',
                  'promoter_gene_score','disease'))
  test
  
})
names(foo)<-names(ds)

setTree <- getTreeYAML(file.path(DATA.DIR,'support/tnact_tree.yaml'))

par_yaml<-'
name: overall
coding:
  name: coding
noncoding:
  promoter:
    name: promoter
  interaction:
    name: interaction
'

osList <- yaml.load(par_yaml)
n <- as.Node(osList)
chil <- setTree$children
n$noncoding$interaction$AddChildNode(chil[[1]])
n$noncoding$interaction$AddChildNode(chil[[2]])
#n$noncoding$interaction$AddChildNode(setTree$Lymphoid)
#n$noncoding$interaction$AddChildNode(setTree$Myeloid)

foo<-foo[!names(foo) %in% c('PBC','RA.UK','RA.US','COELIAC')]
all.pri<-lapply(names(foo),function(disease){
  print(paste("Processing",disease))

merged<-foo[[disease]]

merged.f<-subset(merged,overall_gene_score >= args[['ppi.thresh']])


setTree<-getTreeYAML(file.path(DATA.DIR,'support/tnact_tree.yaml'))

par_yaml<-'
name: overall
coding:
  name: coding
noncoding:
  promoter:
    name: promoter
  interaction:
    name: interaction
'

osList <- yaml.load(par_yaml)
n<-as.Node(osList)
chil<-setTree$children
n$noncoding$interaction$AddChildNode(chil[[1]])
n$noncoding$interaction$AddChildNode(chil[[2]])
#n$noncoding$interaction$AddChildNode(setTree$Lymphoid)
#n$noncoding$interaction$AddChildNode(setTree$Myeloid)



convertNodeToDTName<-function(node){
  if(node$name %in% c('overall','coding','noncoding','promoter','interaction'))
    return(paste0(node$name,'_gene_score'))
  if(node$isLeaf)
    return(paste0(node$name,'_prey_only_gene_score'))
  return(paste0('set.',node$name,'_prey_only_gene_score'))
}

avPPi<-function(node){
  return(rowMeans(merged.f[,node$Get('dtName',filterFun = isLeaf),with=FALSE],na.rm=TRUE))
}

maxPPi<-function(node){
  ## return max and adjust for numerical issues
  ## first get the leaf columns
  #tmp<-as.data.frame(merged.f[,node$Get('dtName',filterFun = isLeaf),with=FALSE])
  tmp<-as.data.frame(merged.f[,node$Get('dtName'),with=FALSE])
  indx <- max.col(tmp, ties.method='first')
  pmin(tmp[cbind(1:nrow(tmp), indx)],1)
  #tmp<-tmp[,list(maxPPi=max(tmp,na.rm=TRUE)),by=1:nrow(tmp)]
  #tmp<-merged.f[,list(maxPPi=max(.(get(node$Get('dtName',filterFun = isLeaf)),na.rm=TRUE))),by=1:nrow(merged.f)]
  #return(tmp[order(tmp$nrow),]$maxPPi)
  #return(pmin(merged.f[,max(get(node$Get('dtName',filterFun = isLeaf)),na.rm=TRUE),by=1:nrow(merged.f)]$V1,1))
}

ppi<-function(node){
  pmax(merged.f[[node$dtName]],NON_ZERO)
  #merged.f[[node$dtName]]
}

BF<-function(node){
  if(!node$isLeaf){
    ch<-node$children
    return(ch[[1]]$ppi/ch[[2]]$ppi)
  }
  return(rep(-1,length(node$avPPi)))
}

## find the best path by scoring as to whether avPPi improves

bscore<-function(node){
  #if(!node$isLeaf & !node$isRoot){
  if(!node$isRoot){
    parent.avppi<-node$parent$avPPi
    current.avppi<-node$avPPi
    #we effectively terminate branch where avppi stops increasing
    term.branch.idx<-which(parent.avppi>=current.avppi)
    current.avppi[term.branch.idx]<-0
    return(current.avppi)
  }
  return(rep(-1,length(node$avPPi)))
}

wBF<-function(node){
  if(!node$isRoot){
    par<-node$parent
    ch<-par$children
    sib.name<-names(ch)[names(ch) != node$name]
    sib<-ch[[sib.name]]
    wBF<-node$ppi/sib$ppi
    # for those below threshold set to zero
    wBF[wBF<args[['BF.thresh']]]<-0
    return(wBF)
  }
  return(args[['BF.thresh']])
}

cumBF<-function(node){
  if(node$isRoot)
    return(args[['BF.thresh']])
  ancwBF<-do.call("cbind",node$Get('wBF',traversal = "ancestor"))
  return(apply(ancwBF, 1, prod,na.rm=TRUE))
}




## this allows us to easilly map between nodes and the datatable containing values
n$Do(function(node) node$dtName=convertNodeToDTName(node))
n$Do(function(node) node$ppi=ppi(node) )
n$Do(function(node) node$avPPi=avPPi(node) )
n$Do(function(node) node$maxPPi=maxPPi(node) )
n$Do(function(node) node$BF=BF(node) )
n$Do(function(node) node$wBF=wBF(node) )
n$Do(function(node) node$score=cumBF(node))
n$Do(function(node) node$nBF=ifelse(node$BF>1,node$BF,1/node$BF) )



## only consider nodes where the ppi is above our threshold
res.dt <- do.call("rbind",lapply(Traverse(n),function(node) {
  dt <-
    data.table(
      ensg = merged.f$ensg, biotype=merged.f$biotype,name = merged.f$name,score = node$score,appi = node$avPPi,maxppi=pmin(node$maxPPi,1),node =
        node$name,wBF = node$wBF,nBF = node$nBF, isLeaf = node$isLeaf, overall_ppi=merged.f$overall_gene_score
    )
}))

## next we
#res.dt<-res.dt[,nBF:=ifelse(BF<1,1/BF,BF)]
h<-res.dt[,.SD[which.max(.SD$score),],by=ensg]
h<-h[order(h$score),]
h$disease<-disease
h
})
all.pri<-rbindlist(all.pri)
# priortised.results.file<-file.path(args[['out_dir']],paste0(file.stub,'_prioritised.tab'))
# write.table(h,file=priortised.results.file,sep="\t",row.names=FALSE,quote=FALSE)
# message(paste('written',priortised.results.file))

## add immunobase annotationsassoc_table

ic.stub<-c('AA','AS','ATD','CEL','CRO','IBD','JIA','MS','NAR','PBC','PSC','PSO','RA','SJO','SLE','SSC','T1D','UC','VIT')


trait.regions<-lapply(seq_along(ic.stub),function(i){
  t<-fread(paste0('https://www.immunobase.org/downloads/regions-files-archives/latest_1.11/Hs_GRCh37-',ic.stub[i],'-assoc_table'))
  t$disease<-ic.stub[i]
  colnames(t)<-make.names(colnames(t))
  t
})


#ic.stub<-c('AA','AS','ATD','CEL','CRO','IBD','IGE','JIA','MS','NAR','PBC','PSC','PSO','RA','SJO','SLE','SSC','T1D','UC','VIT')

trait.genes<-lapply(seq_along(ic.stub),function(i){
  print(paste("Getting",ic.stub[i]))
  t<-fread(paste0('https://www.immunobase.org/downloads/regions-files-archives/latest_1.11/Hs_GRCh37-',ic.stub[i],'-assoc_genesTAB'))
  t$disease<-ic.stub[i]
  t
})

tg<-do.call("rbind",trait.genes)
setnames(tg,sub(" ","_",names(tg)))
cc.genes<-unique(subset(tg,Cand_Gene==1)$Ensembl_ID)
in.region.genes<-unique(subset(tg,In_Region==1)$Ensembl_ID)

# pos.reg<-unique(tg$Region_Pos)
# disease_region<-lapply(split(tg$Region_Pos,tg$disease),function(pos.reg){
#   pchr<-sub("([^:]+):([0-9]+)\\.\\.([0-9]+)","\\1",unique(pos.reg))
#   pstart<-sub("([^:]+):([0-9]+)\\.\\.([0-9]+)","\\2",unique(pos.reg))
#   pend<-sub("([^:]+):([0-9]+)\\.\\.([0-9]+)","\\3",unique(pos.reg))
#   imb.reg<-GRanges(seqname=Rle(pchr),ranges=IRanges(start=as.numeric(pstart),end=as.numeric(pend)))
#   unique(subsetByOverlaps(all,imb.reg)$id)
# })

disease_gene<-lapply(split(tg,tg$disease),function(t){
  unique(subset(t,Cand_Gene==1)$Ensembl_ID)
})

disease_region_gene<-lapply(split(tg,tg$disease),function(t){
  unique(subset(t,In_Region==1)$Ensembl_ID)
})

## add closest region information using gene locations




all.pri<-subset(all.pri,disease %in% c('ICOELIAC','JIA','RA','T1D','MS','GRAVES'))
all.pri[all.pri$disease=='ICOELIAC',]$disease<-'CEL'
all.pri[all.pri$disease=='GRAVES',]$disease<-'ATD'

g<-lapply(split(all.pri,all.pri$disease),function(d){
  dis<-unique(d$disease)
  #d$frag.overlaps.reg<-d$frag.id %in% disease_region[[dis]]
  d$disease.gene<-d$ensg %in% disease_gene[[dis]]
  d$disease.region.gene<-d$ensg %in% disease_region_gene[[dis]]
  d
})

g<-rbindlist(g)



## to start with look at all genes identified in AI traits
# tg<-do.call("rbind",trait.genes)
# setnames(tg,sub(" ","_",names(tg)))
# cc.genes<-unique(subset(tg,Cand_Gene==1)$Ensembl_ID)
# in.region.genes<-unique(subset(tg,In_Region==1)$Ensembl_ID)
# all.pri$imb.candidate.causal<-all.pri$ensg %in% cc.genes
# all.pri$imb.region<-all.pri$ensg %in% in.region.genes
write.table(g,file="/home/ob219/scratch/DATA/JAVIERRE_ICHIP/out/GUESSFM_JUNE_13/gene_prioritisation_0.01.csv",quote=FALSE,sep=",",row.names = FALSE)

exit()
### PLOT STUFF

merged<-rbindlist(foo)
merged.f <- subset(merged,overall_gene_score >= args[['ppi.thresh']] & biotype=='protein_coding')

convertNodeToDTName <- function(node) {
  if (node$name %in% c('overall','coding','noncoding','promoter','interaction'))
    return(paste0(node$name,'_gene_score'))
  if (node$isLeaf)
    return(paste0(node$name,'_prey_only_gene_score'))
  return(paste0('set.',node$name,'_prey_only_gene_score'))
}

avPPi<-function(node){
  return(rowMeans(merged.f[,node$Get('dtName',filterFun = isLeaf),with=FALSE],na.rm=TRUE))
}

ppi<-function(node){
  pmax(merged.f[[node$dtName]],NON_ZERO)
  #merged.f[[node$dtName]]
}

BF<-function(node){
  if(!node$isLeaf){
    ch<-node$children
    return(ch[[1]]$ppi/ch[[2]]$ppi)
  }
  return(rep(-1,length(node$avPPi)))
}

## find the best path by scoring as to whether avPPi improves

bscore <- function(node) {
  #if(!node$isLeaf & !node$isRoot){
  if (!node$isRoot) {
    parent.avppi <- node$parent$avPPi
    current.avppi <- node$avPPi
    #we effectively terminate branch where avppi stops increasing
    term.branch.idx <- which(parent.avppi >= current.avppi)
    current.avppi[term.branch.idx] <- 0
    return(current.avppi)
  }
  return(rep(-1,length(node$avPPi)))
}

wBF<-function(node){
  if(!node$isRoot){
    par<-node$parent
    ch<-par$children
    sib.name<-names(ch)[names(ch) != node$name]
    sib<-ch[[sib.name]]
    wBF<-node$ppi/sib$ppi
    # for those below threshold set to zero
    wBF[wBF<BF.thresh]<-0
    return(wBF)
  }
  return(BF.thresh)
}

cumBF<-function(node){
  if(node$isRoot)
    return(BF.thresh)
  ancwBF<-do.call("cbind",node$Get('wBF',traversal = "ancestor"))
  return(apply(ancwBF, 1, prod,na.rm=TRUE))
}


BF.thresh<-3
## this allows us to easilly map between nodes and the datatable containing values
n$Do(function(node)
  node$dtName = convertNodeToDTName(node))
n$Do(function(node)
  node$ppi = ppi(node))
n$Do(function(node)
  node$avPPi = avPPi(node))
n$Do(function(node)
  node$BF = BF(node))
n$Do(function(node)
  node$wBF = wBF(node))
n$Do(function(node)
  node$score = cumBF(node))
n$Do(function(node)
  node$nBF = ifelse(node$BF > 1,node$BF,1 / node$BF))



## only consider nodes where the ppi is above our threshold
res.dt <- do.call("rbind",lapply(Traverse(n),function(node) {
  dt <-
    data.table(
      ensg = merged.f$ensg, name = merged.f$name,score = node$score,appi = node$avPPi,node =
        node$name,wBF = node$wBF,nBF = node$nBF, isLeaf = node$isLeaf
    )
}))

## next we
#res.dt<-res.dt[,nBF:=ifelse(BF<1,1/BF,BF)]
h <- res.dt[,.SD[which.max(.SD$score),],by = c('ensg')]
h <- h[order(h$score),]




#lapply(foo,function(x) if(is.data.table(x)) subset(x,name=="BACH2") )

## plotting fun ;)

## here we implement for a given gene a decision tree and then spit it out

ppi.sg <- function(n,geneName) {
  return(pmax(subset(merged.f,name == geneName)[[n$dtName]],NON_ZERO))
}

BF.sg <- function(n,geneName) {
  if (!n$isLeaf)
    return(oBF <- subset(res.dt,name == geneName &
                           node == n$name)$nBF)
  return(0)
}

getBranchPath <- function(node,geneName) {
  node.name <- subset(h,name == geneName)$node
  path <-
    node$Get(
      function(x)
        return(x$pathString),traversal = "post-order",filterFun = function(y)
          y$name == node.name
    )
  print(path)
  return(unlist(strsplit(path,split = '/')))
}

selected <- function(n,node.list) {
  return(n$name %in% node.list)
}

above.thresh <- function(n,thresh) {
  if (is.na(n$ppi))
    return(FALSE)
  return(n$ppi >= thresh)
}

gn <- Clone(n)
selected.branch.path <- getBranchPath(gn,args[['geneName']])
gn$Do(function(node)
  node$ppi = ppi.sg(node,args[['geneName']]))
gn$Do(function(node)
  node$BF = BF.sg(node,args[['geneName']]))
gn$Do(function(node)
  node$selected = selected(node,selected.branch.path))
gn$Do(function(node)
  node$above.thresh = above.thresh(node,all.thresh))

library(ape)
gn$Revert()
gnp <- as.phylo(gn)

Nodelabel <- function(node) {
  if (node$isLeaf)
    return(paste0(
      node$name,': PPi ', format(
        signif(node$ppi,digits = 2), scientific = ifelse(node$ppi < 0.001,TRUE,FALSE), big.mark = "'"
      )
    ))
  if (node$above.thresh)
    return(paste0(node$name,'\n BF ', ifelse(
      is.na(node$BF),'NA',format(
        signif(node$BF,digits = 2), scientific = ifelse(node$BF > 100,TRUE,FALSE), big.mark = "'"
      )
    )))
  return(paste0(node$name,'\n'))
}

getFormat <- function(node) {
  if (node$selected)
    return(list(col = 'red',lw = 2,tcol = 'red'))
  if (node$above.thresh)
    return(list(col = 'black',lw = 1,tcol = 'black'))
  return(list(col = 'grey',lw = 1,tcol = 'white'))
}
disease<-'AI'
par(mar = c(1,1,1,1))
plot(
  gnp, show.tip.label = TRUE, type = "cladogram",edge.color = "grey",edge.width =
    0.5,plot = FALSE,main = paste(toupper(disease),args[['geneName']],sep =
                                    ':')
)

for (node in Traverse(gn)) {
  f <- getFormat(node)
  print(f)
  if (!node$isRoot)
    edges(
      GetPhyloNr(node$parent, "node"), GetPhyloNr(node, "node"), lwd = f[['lw']],col = f[['col']]
    )
  if (!node$isLeaf) {
    nodelabels(
      Nodelabel(node), GetPhyloNr(node, "node"), frame = 'none', adj = c(0.5, 0.2),col = f[['col']],cex =
        0.8
    )
  }
  if (node$isLeaf)
    tiplabels(
      Nodelabel(node), GetPhyloNr(node, "node"), frame = "none", adj = c(0, 0),col =
        f[['col']],cex = 0.8
    )
  
  
}

## another thing one might do is plot the pathweight for a complete disease for all selected genes

## get unique nodes

wc <- h[,list(ncount = nrow(.SD)),by = "node"]
wn <- Clone(n)
wnp <- as.phylo(wn)
assignWeights <- function(n,dt) {
  print(n$name)
  w <- subset(dt,node == n$name)$ncount
  if (length(w) == 0)
    w <- 0
  return(w)
}

wn$Do(function(node)
  node$weight = assignWeights(node,wc))

Nodelabel2 <- function(node) {
  cw <- node$cWeight
  w <- paste0('(',node$weight,')')
  if (!node$isLeaf)
    return(paste0(node$name,' ',ifelse(is.na(cw),0,cw),w))
  return(paste0(node$name,' ',ifelse(is.na(cw),0,cw)))
}

cumWeight <- function(node) {
  w <- sum(sum(node$Get('weight')))
  if (w == 0)
    w <- return(NA)
  #return(log(w))
  return(w)
}

lineWidthFactor <- function(node) {
  w <- sum(sum(node$Get('weight')))
  if (w == 0)
    return(0)
  if (w == 1)
    w <- 2
  ceiling(log(w)) * 1.5
}

wn$Do(function(node)
  node$cWeight = cumWeight(node))
wn$Do(function(node)
  node$lw = lineWidthFactor(node))



#tw<-pmin(wn$Get(cumWeight), mean(sapply(wn$noncoding$interaction$children,cumWeight)))
tw <- wn$Get(lineWidthFactor)
rbPal <- colorRampPalette(c('green','red'))
tw <- split(rbPal(10)[as.numeric(cut(tw,breaks = 10))],names(tw))


par(mar = c(1,1,1,1))
plot(
  wnp, show.tip.label = TRUE, type = "cladogram",edge.color = "grey",edge.width =
    0.5,plot = FALSE,main = toupper(disease)
)
lw.scale.factor <- 20 / nrow(h)
cex.scale.factor <- 2 / nrow(h)
for (node in Traverse(wn)) {
  pw <- sum(sum(node$Get('weight')))
  if (!node$isRoot)
    #edges(GetPhyloNr(node$parent, "node"), GetPhyloNr(node, "node"), lwd=pw*lw.scale.factor,col = tw[[node$name]])
    if (node$lw == 0) {
      edges(
        GetPhyloNr(node$parent, "node"), GetPhyloNr(node, "node"), lwd = 0.5,col = 'grey'
      )
    }else{
      edges(
        GetPhyloNr(node$parent, "node"), GetPhyloNr(node, "node"), lwd = node$lw,col = tw[[node$name]]
      )
    }
  if (!node$isLeaf & node$lw != 0)
    nodelabels(
      Nodelabel2(node), GetPhyloNr(node, "node"), frame = 'none', adj = c(-0.3, 0.1),col = "black",cex =
        0.7
    )
  if (node$isLeaf) {
    if (node$lw != 0) {
      tiplabels(
        Nodelabel2(node), GetPhyloNr(node, "node"), frame = "none", adj = c(0, 0),col = tw[[node$name]],cex =
          0.8
      )
    }else{
      tiplabels(
        Nodelabel2(node), GetPhyloNr(node, "node"), frame = "none", adj = c(0, 0),col = 'grey',cex =
          0.8
      )
    }
  }
}

