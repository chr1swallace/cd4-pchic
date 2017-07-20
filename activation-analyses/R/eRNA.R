#b# NB run first three stages of activation-analyses first.
## This script analyses ernas.

quick <- TRUE

## * basic summary plots

erna <- get.erna()
urna <- unique(erna,by="id")
table(expr=urna$expr,type=urna$type,P=!is.na(urna$FDR.erna))
table(urna$expr,!is.na(urna$FDR.erna))

## ernap <- fread(file.path(CD4CHIC.DATA,"rna-with-regrna-diff-v2.csv"))
## erna[,AveExpr:=NULL]
## erna[,t:=NULL]
## erna[,P.Value:=NULL]
## erna <- merge(erna,ernap[,.(id,logFC,AveExpr,t,P.Value,adj.P.Val)],by="id",all.x=TRUE)
## setnames(erna,c("logFC","adj.P.Val"),c("logFC.erna","FDR,erna"))
## table(erna$type)

if(!quick) {
##' volcano plots
ggplot(unique(erna,by="id"),
       aes(x=logFC.erna,y=-log10(FDR.erna))) + geom_point() + facet_wrap(~type) 

##' distance distributions
ggplot(unique(erna,by="id"), aes(x=(end-start))) +
  geom_density() +
  geom_histogram(aes(y=..density..), binwidth=1000) +
  facet_grid(type ~ .,scales="free") + xlim(0,50000)
}

cutter <- function(x,n=4) {
  q <- quantile(x,seq(0,1,length=n+1),na.rm=TRUE)
  cut(x,breaks=q,include.lowest=TRUE)
}

mods <- get.modules()

## get oe frags for each module
f <- function(mod) {
    ids <- unique(mods[module==mod,]$id)
    oes <- ints[pmax(Total_CD4_NonActivated,Total_CD4_Activated)>5 &
                id %in% ids, ]$oeID %>% unique()
    e.mod <- erna[name %in% oes,]
    e.mod$module <- mod
    return(unique(e.mod,by="id"))
}

library(parallel)
data <- mclapply(unique(mods$module), f, mc.cores=8) %>% do.call("rbind",.)
if(!quick) {
ggplot(data[type %in% c("regulatory","intergenic.reg"),],
       aes(x=logFC.erna,y=-log10(FDR.erna),col=module)) +
geom_point() + facet_grid(type~module) 
}

## * Add gene, hind annotations
b2g[,baitID:=as.integer(baitID)]

hind <- get.hind()

h <- as.data.table(as.data.frame(hind))
h[,mid:=(start+end)/2]
h[,name:=as.integer(name)]
h <- h[order(seqnames,start),]

## TODO: take account of overlaps when defining nearest - should be
## nearest bait with non-overlapping feature
feat <- copy(erna)
feat[,c("minfrag","maxfrag"):=list(min(name),max(name)), by="id"]
feat[,mid:=pmean(start,end)]
## feat <- unique(feat[,.(id,chr,type,start,end,minfrag,maxfrag)],
##                 by="id")
setkey(feat,id)
overlap <- function(id1,id2) {
    o.length <- (feat[id1,start] < feat[id2,end]) &  (feat[id1,end] > feat[id1,start])
    o.frag <- (feat[id1,minfrag] <= feat[id2,maxfrag]) & (feat[id1,maxfrag] >= feat[id1,minfrag])
    o.length | o.frag
}

## for each regulatory eRNA, add fold change at + distance to all nearby protein coding genes
## NB a reg feature can map to multiple hindIII, but a protein coding gene has generally a single bait
feat.reg <- feat[type %in% c(## "regulatory",
"intergenic.reg") & expr,]#!is.na(FDR.erna),]
#feat.reg <- feat[type=="intergenic.reg" & !is.na(FDR.erna),]
feat.prot <- unique(feat[type=="protein_coding" & expr,],#!is.na(FDR.erna),],
                    by="id")

nrow(feat.reg)
nrow(feat.prot)
feat.reg <- split(feat.reg,feat.reg$chr)
sapply(feat.reg,nrow)
feat.prot <- split(feat.prot,feat.prot$chr)
sapply(feat.prot,nrow)

## * Correlation as a function of distance for different classes

## all pairs of features by chromosome
DATA <- structure(vector("list",length(feat.reg)), names=names(feat.reg))
for(v in unique(names(feat.reg))) {
    myf <- function(i, j, col) 
        abs(feat.reg[[v]][i][[col]] - feat.prot[[v]][j][[col]])
    O <- outer(1:nrow(feat.reg[[v]]), 1:nrow(feat.prot[[v]]), myf, col="mid")
    myreg <- function(col)
        matrix(feat.reg[[v]][[col]], nrow=nrow(O), ncol=ncol(O)) %>% as.vector()
    myprot <- function(col)
        matrix(feat.prot[[v]][[col]], ncol=nrow(O), nrow=ncol(O)) %>% t() %>% as.vector()
    DATA[[v]] <- data.table(dist=as.numeric(O),
                            REG=myreg("logFC.erna"),
                            PROT=myprot("logFC.erna"),
                            REG.fdr=myreg("FDR.erna"),
                            PROT.fdr=myprot("FDR.erna"),
                            REG.id=myreg("id"),
                            PROT.id=myprot("id"),
                            REG.start=myreg("start"),
                            PROT.start=myprot("start"),
                            REG.end=myreg("end"),
                            PROT.end=myprot("end"),
                            oeID=myreg("name"),
                            chr=v)
}
DATA <- do.call("rbind",DATA)

b2g <- get.b2gene()
table(b2g$id %in% DATA$PROT.id)
DATA <- merge(DATA, b2g[,.(id,baitID,gene)], by.x="PROT.id", by.y="id", allow.cartesian=TRUE)

## add max interaction scores
DATA <- merge(DATA,
              ints[,.(baitID,oeID,Total_CD4_Activated,Total_CD4_NonActivated,Background)],
              by.x=c("baitID","oeID"), by.y=c("baitID","oeID"), all.x=TRUE)
setnames(DATA, c("Total_CD4_Activated","Total_CD4_NonActivated","Background"), c("act","non","back"))
myf <- function(x) {
    x <- max(x, na.rm=TRUE)
    if(is.na(x))
        return(0)
    return(x)
}
with(DATA, tapply(dist, pmax(act,non,back)>5, summary))
DATA <- DATA[dist<1e+7,]
setkeyv(DATA, c("baitID","REG.id"))
DATA[, c("act","non","back") := list(max(ifelse(is.na(act), 0, act), na.rm=TRUE),
                                     max(ifelse(is.na(non), 0, non), na.rm=TRUE),
                                     max(ifelse(is.na(back), 0, back), na.rm=TRUE)),
     by=c("baitID","REG.id")]
DATA <- unique(DATA, by=c("baitID","REG.id"))

## group dist
DATA[,cutdist:=droplevels(cut(dist,breaks=50,include.lowest=TRUE))]
table(DATA$cutdist, DATA$non>5)
with(DATA, table(cutdist, pmax(non,act,back)>5))

if(!quick)
    ggplot(DATA[act>5,], aes(x=REG,y=PROT)) + geom_point() + geom_hline(yintercept=0,col="grey") + geom_vline(xintercept=0,col="grey")

## * Correlation between bait and prey
myf <- NULL
mycor <- function(x,y,use=rep(TRUE,length(x))) {
    if(sum(use)>=3)
        return( cor(x[use], y[use], use="pair") )
    return(as.numeric(NA))
}
myf <- function(x) {
    ss <- strsplit(gsub("\\[|\\(|\\]","",x, ","), ",")
    ss <- lapply(ss, as.numeric)
    sapply(ss, mean)
}

Dcor <- copy(DATA)
Dcor <- Dcor[PROT.fdr<0.1,]
Dcor[, all := mycor(PROT,REG), by="cutdist"]    
Dcor[, act := mycor(PROT,REG,act>5), by="cutdist"]    
Dcor[, non := mycor(PROT,REG,non>5), by="cutdist"]    
Dcor[, back := mycor(PROT,REG, back>5), by="cutdist"]    
Dcor <- unique(Dcor,by="cutdist")
Dcor[,mid:=myf(cutdist)]
Dcor[is.na(mid), mid:=1.3e+6]
Dcor <- melt(Dcor, c("mid","cutdist"))
Dcor <- Dcor[order(Dcor$mid),]
if(!quick)
    ggplot(Dcor,aes(x=mid,y=value,col=variable)) + geom_point() + geom_smooth() + xlim(0,2.0E+06)

## * Categorise reg/prot pairs

DATA[,int.cat:=paste0(ifelse(act>5 & non<5,"act",
                      ifelse(act<5 & non>5, "non",
                      ifelse(act>5 & non>5, "cd4",
                      ifelse(back>5 & act<3 & non<3, "back", "--")))))]
table(DATA$int.cat)
with(unique(DATA,by=c("baitID","oeID")),table(int.cat))

DATA[,enh.cat:=ifelse(REG.fdr>0.01,"invar",
                      ifelse(REG>0,"up","down"))]
with(unique(DATA,by="REG.id"),table(enh.cat))

D <- DATA[REG.fdr<0.01 & PROT.fdr<0.01, ]
#D <- DATA
D[,cat:=int.cat]
table(D$cat)
D <- DATA[int.cat!="--",]
save(D, file=file.path(CD4CHIC.DATA,"prot-reg-table-gb.RData"))

D <- DATA#[REG.fdr<0.05, ]
D[,cat:=int.cat]
#D[, .( agree=sum(sign(REG)==sign(PROT))/.N ),
#  by="cat"]


source(file.path(CD4CHIC.ROOT, "activation-analyses/R/robust-models.R"))

D[,y:=sign(REG)==sign(PROT)]
D[,mb:=log10(dist/1e+6)]
D <- D[!is.na(y),]
## linear


D$cat <- factor(D$cat, levels=c("--","back","non","act","cd4"))
mod <- model(PROT ~ mb*REG*cat, x=as.data.frame(D),do.print=FALSE)
pmod <- modelprint(mod)  %>% as.data.frame()
pmod$int <- sub(".*cat=","",rownames(pmod)) %>% factor(., levels=c("back","non","act","cd4"))
pmod

## logistic
D[,mb2:=mb^2]
D[,mb3:=mb^3]

P <- vector("list",4)
names(P) <- setdiff(D$cat,"--")
for(l in names(P)) {
    message(l)
    D[,mb2:=ifelse(cat==l,mb,0)]
    P[[l]] <- lmodel(y ~ mb + mb2, x=as.data.frame(D[cat %in% c("--",l),]))["mb2","P"]
    m0 <- lm(y ~ mb,data=as.data.frame(D[cat %in% c("--",l),]))
    m1 <- lm(y ~ mb*cat,data=as.data.frame(D[cat %in% c("--",l),]))
    print(anova(m0,m1))
}
P


mod <- lmodel(y ~ mb*cat, x=as.data.frame(D),do.print=FALSE)
pmod <- modelprint(mod)  %>% as.data.frame()
pmod$int <- sub("cat=","",rownames(pmod)) %>% factor(., levels=c("back","non","act","cd4"))
pmod

modelcmp(list(y ~ mb*cat, y ~ mb), x=D)

pdata <- expand.grid(mb=seq(-2,0.5,0.1), cat=unique(D$cat)) %>% as.data.table()
#pdata[,mb2:=mb^2]
#pdata[,mb3:=mb^3]
mm <- D[,.(min=min(mb),max=max(mb)),by="cat"]
mn <- mm$min
mx <- mm$max
names(mn) <- mm$cat
names(mx) <- mm$cat
pdata <- pdata[mb >= mn[cat] & mb <= mx[cat], ]
pr <- predict(mod, newdata=pdata, se.fit=TRUE)
pdata$ypred <- pr$linear.predictors
pdata$ylo <- pdata$ypred - 1.96 * pr$se.fit
pdata$yhi <- pdata$ypred + 1.96 * pr$se.fit
pdata <- split(as.data.table(pdata), pdata$cat)
for(nm in rev(names(pdata))) {
    pdata[[nm]] <- rbind(pdata[[nm]],pdata[["--"]])
    pdata[[nm]][,cat2:=nm]
}
pdata <- do.call("rbind",pdata)
pdata <- pdata[cat2!="--",]
pdata[,cat2:=factor(cat2, levels=c("back","non","act","cd4"))]
levels(pdata$cat2) <- c("control cells","CD4 lost","CD4 gained", "CD4 invariant")
pdata[,cat:=ifelse(cat=="--","no","yes")]

dP <- t(as.data.frame(P))
colnames(dP) <- "P"
dP <- as.data.frame(dP)
dP$cat <- rownames(dP)
dP <- as.data.table(dP)
dP$P <- format.pval(dP$P)
dP$P <- sub("9e-05","0.00009",dP$P); dP$P
dP$P <- sub("0e.*","P<10^{-16}",dP$P)
dP$P <- ifelse(grepl("<",dP$P),dP$P, paste0("'P='*'",dP$P,"'"))
dP[,cat2:=factor(cat, levels=c("back","non","act","cd4"))]
levels(dP$cat2) <- c("control cells","CD4 lost","CD4 gained", "CD4 invariant")

cols <- c("grey40",scales:::hue_pal()(1))
ggplot(pdata, aes(x=exp(mb),y=ypred)) +
background_grid() +
geom_ribbon(aes(ymin=ylo,ymax=yhi,fill=cat),alpha=0.1) +
geom_path(aes(col=cat,fill=cat)) +
facet_wrap(~cat2) +
geom_text(aes(label=P),x=1.4,y=1.25,data=dP,parse=TRUE,size=4) +
scale_x_continuous("Distance between eRNA and gene promoter (Mb)") +
scale_y_continuous("log odds of agreement in fold change direction") +
scale_colour_manual("PCHi-C interaction detected",labels=c("no","yes"),values=cols) +
scale_fill_manual("PCHi-C interaction detected",labels=c("no","yes"),values=cols) +
theme(legend.position="top")

f <- file.path(CD4CHIC.OUT,"paper/figure-erna-prot-agreement.tiff")
ggsave(f,height=4.5,width=6,units="in")
ggsave(sub(".tiff",".pdf",f),height=4.5,width=6,units="in")
system(paste("display",f))
write.table(pmod,file=file.path(CD4CHIC.OUT,"paper/table-erna-prot-agreement.csv"))


## * basic summary stats

head(erna)
head(ints)
erna[,act:=name %in% ints[Total_CD4_Activated>5,][["oeID"]]]
erna[,non:=name %in% ints[Total_CD4_NonActivated>5,][["oeID"]]]
erna <- erna[pmax(erna.reads.act,erna.reads.non)>1,]
cerna <- erna
cerna[,c("act","non"):=list(any(act),any(non)),by="id"]
cerna <- unique(cerna,by="id")
with(cerna,ftable(diff=!is.na(logFC.erna),act,non,type=type))


for(l in unique(erna$type)) {
    tt <- with(erna[type==l,], table(int=act|non, diff=!is.na(logFC.erna)))
    print(tt)
    message(l,"\t", format.pval(chisq.test(tt)$p.value)) 
}

with(erna[type=="intergenic.reg" & FDR.erna<0.01,],
     table(sign(logFC.erna)))
table(DATA$cat)
D <- DATA[cat!="--",]

cat(sort(unique(D$gene)),sep="\n") # AHR, BACH2, BATF, 




