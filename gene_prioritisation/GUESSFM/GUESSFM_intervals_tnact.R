library(GenomicRanges)
library(data.table)
library(rtracklayer)

data.dir<-file.path(GRPATH,'cd4chic','ichip_paper')

(load(file.path(data.dir,'DATA','RDATA','javierre_tnact_frags.by.ld.RData')))
(load(file.path(data.dir,'DATA','RDATA','javierre_tnact_csnps.by.ld.RData')))
(load(file.path(data.dir,'DATA','RDATA','javierre_tnact_interactions.RData')))

#dense.regions.gr<-import.bed(file.path(data.dir,'support','dense.ic.regions.shuffled.bed'))

## only include genes that have at least one fragment overlapping a dense region

## cSNPs
csnp.ichip.gr<-subset(cs.gr,ensg %in% int$ensg)

## noncoding (anything in promoters and or interaction)

nc.frags.gr<-subset(frags.gr,type=='promoter' || type=='interaction')

nc.frags.gr$category<-'non_coding'

## promoters

promoter.frags.gr<-subset(frags.gr,type=='promoter')

promoter.frags.gr$category<-'promoter_only'

## any interaction

interaction.frags.gr<-subset(frags.gr,type=='interaction')

interaction.frags.gr$category<-'interaction_only'

## act interaction

interaction.act.frags.gr<-subset(frags.gr,type=='interaction' & id %in% subset(int,int$Total_CD4_Activated)$oeID)

interaction.act.frags.gr$category<-'act_interaction_only'

## nact interaction

interaction.nact.frags.gr<-subset(frags.gr,type=='interaction' & id %in% subset(int,int$Total_CD4_NonActivated)$oeID)

interaction.nact.frags.gr$category<-'non_act_interaction_only'

## overall

overall.frags.gr<-frags.gr
overall.frags.gr$category<-'overall'

GUESS.ichip.grl<-GRangesList(
  non_coding=nc.frags.gr,
  promoter_only=promoter.frags.gr,
  interaction_only=interaction.frags.gr,
  interaction_act_only=interaction.act.frags.gr,
  interaction_nact_only=interaction.nact.frags.gr,
  overall=overall.frags.gr
  )

save(GUESS.ichip.grl,file=file.path(data.dir,'DATA','RDATA','GUESSFM_intervals_tnact.RData'))
