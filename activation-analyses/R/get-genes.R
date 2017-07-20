library(biomaRt)
library(GenomicRanges)

mart <- useMart('ENSEMBL_MART_ENSEMBL',host="grch37.ensembl.org")
ensembl.gene <- useDataset("hsapiens_gene_ensembl",mart=mart)
gene.details<-getBM(
#        filters= c("chromosome_name","start","end"),
        attributes= c('ensembl_gene_id','chromosome_name','exon_chrom_start','exon_chrom_end','ensembl_transcript_id','external_gene_name','strand','gene_biotype'),
#        values= ss,
        mart= ensembl.gene)
head(gene.details)

transcript.gr<-with(gene.details,
                    GRanges(seqnames=Rle(chromosome_name),
                            ranges=IRanges(start=exon_chrom_start,end=exon_chrom_end),
                            tid=ensembl_transcript_id,
                            strand=strand,
                            genename=external_gene_name,
                            biotype=gene_biotype))
genes<-transcript.gr[transcript.gr$biotype=="protein_coding",]

