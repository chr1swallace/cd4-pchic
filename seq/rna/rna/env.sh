P=/dunwich/scratch/ar756/.local

export PATH=$P/parallel/20150322/bin:$PATH
export PATH=$P/fastqc/0.11.3/bin:$PATH
export PATH=$P/cutadapt/1.8.1/bin:$PATH
export PATH=$P/star/2.4.1a/bin:$PATH
export PATH=$P/samtools/1.2/bin:$PATH
export PATH=$P/htseq/0.6.1p1/bin:$PATH
export PATH=$P/subread/1.4.6-p1/bin:$PATH

export PYTHONPATH=$P/cutadapt/1.8.1/lib/python2.7/site-packages:$PYTHONPATH
export PYTHONPATH=$P/htseq/0.6.1p1/lib/python2.7/site-packages:$PYTHONPATH
export R_LIBS=$P/r

SRC=/chiswick/data/store/sequencing/RNA-seq/P140369/140707_SN847_0422_BC4H6RACXX
GEN=$P/ensembl/75/star
GTF=$P/ensembl/75/Homo_sapiens.GRCh37.75.gtf
PICARD=$P/picard/1.130/picard.jar

MIN=201
MAX=204
