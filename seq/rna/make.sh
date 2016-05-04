source env.sh

IDX1='GATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
IDX2='ATCTCGTATGCCGTCTTCTGCTTG'
IDX="${IDX1}NNNNNNNN${IDX2}"
PRM='AAGCAGTGGTATCAACGCAGAGTAC'
STR='--runThreadN 24 --genomeDir index'
SAM='-@ 6'

ruby gend.rb > design.csv
cat design.csv | sed '1d' | cut -f 1 > ids.csv

curl -s $DNA | zcat > dna.fa
curl -s $GTF | zcat > ens.gtf
rm -Rf index && mkdir index
STAR $STR --runMode genomeGenerate --genomeFastaFiles dna.fa --sjdbGTFfile ens.gtf

cat design.csv | sed '1d' | cut -f 1,4 | parallel --col-sep $'\t' 'ln -s {2} data_raw.{1}_1.fastq.gz'
cat design.csv | sed '1d' | cut -f 1,5 | parallel --col-sep $'\t' 'ln -s {2} data_raw.{1}_2.fastq.gz'
cat ids.csv | parallel 'fastqc data_raw.{}_?.fastq.gz'

cat ids.csv | parallel "cutadapt -b $IDX -b $IDX1 -b $IDX2 -o data_tr1.{}_1.fastq.gz data_raw.{}_1.fastq.gz > data_tr1.{}_1.log"
cat ids.csv | parallel "cutadapt -b $PRM -O 10 --discard -o data_tr2.{}_2.fastq.gz -p data_tr2.{}_1.fastq.gz data_raw.{}_2.fastq.gz data_tr1.{}_1.fastq.gz > data_tr2.{}.log"
cat ids.csv | parallel 'fastqc data_tr2.{}_?.fastq.gz'

cat ids.csv | parallel -P 1 "STAR $STR --outStd SAM --readFilesIn <(zcat data_tr2.{}_1.fastq.gz) <(zcat data_tr2.{}_2.fastq.gz) | samtools view $SAM -S - -b -o data_aln.{}.bam"
cat ids.csv | parallel 'fastqc data_aln.{}.bam'

cat ids.csv | parallel "samtools view $SAM -q 10 -f 3 -F 256 -b data_aln.{}.bam > data_flt.{}.bam"
cat ids.csv | parallel 'fastqc data_flt.{}.bam'

cat <(cat chrom.bed | ruby genr.rb) <(cat ens.gtf | grep -v '^#') | bedtools sort > cd4.gtf
cat ids.csv | parallel "samtools view data_flt.{}.bam | htseq-count -s reverse -m intersection-strict - ens.gtf > data_gne.{}.csv 2> data_gne.{}.log"
cat ids.csv | parallel "samtools view data_flt.{}.bam | htseq-count -t gene -s reverse -m intersection-strict - cd4.gtf > data_ncd.{}.csv 2> data_ncd.{}.log"
cat ids.csv | parallel 'cat data_gne.{}.csv <(cat data_ncd.{}.csv | grep CD4R) | grep -v '^__'  > data_exp.{}.csv'
