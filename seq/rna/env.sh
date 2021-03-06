source path.sh

export PATH=$P/bedtools/2.25.0/bin:$PATH
export PATH=$P/emacs/24.5/bin:$PATH
export PATH=$P/fastqc/0.11.5/bin:$PATH
export PATH=$P/parallel/20160422/bin:$PATH
export PICARD=$P/picard/2.0.1/picard.jar
export PATH=$P/python/2.7.11/bin:$PATH
export PYTHONPATH=$P/python/2.7.11/lib/python2.7/site-packages:$PYTHONPATH
export PATH=$P/r/3.2.4/bin:$PATH
export PATH=$P/ruby/2.3.1/bin:$PATH
export PATH=$P/samtools/1.3.1/bin:$PATH
export PATH=$P/star/2.5.1b/bin:$PATH

ENS=ftp://ftp.ensembl.org/pub/grch37/release-84
GTF=$ENS/gtf/homo_sapiens/Homo_sapiens.GRCh37.82.chr_patch_hapl_scaff.gtf.gz
DNA=$ENS/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz

export PYTHONPATH=$P/cutadapt/1.8.1/lib/python2.7/site-packages:$PYTHONPATH
export PYTHONPATH=$P/htseq/0.6.1p1/lib/python2.7/site-packages:$PYTHONPATH
export R_LIBS=$P/r

SRC=/chiswick/data/store/sequencing/RNA-seq/P140369/140707_SN847_0422_BC4H6RACXX
GEN=$P/ensembl/75/star
GTF=$P/ensembl/75/Homo_sapiens.GRCh37.75.gtf
PICARD=$P/picard/1.130/picard.jar

MIN=201
MAX=204
