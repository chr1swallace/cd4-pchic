source env.sh

IDX1='GATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
IDX2='ATCTCGTATGCCGTCTTCTGCTTG'
IDX="${IDX1}NNNNNNNN${IDX2}"
PRM='AAGCAGTGGTATCAACGCAGAGTAC'

cat <(echo $'id\tgene\tsource\tbiotype') <(cat $GTF | grep -v '^#' | awk '$3 == "gene"' | cut -f 9 | tr -d ';' | tr -d '"' | tr ' ' $'\t' | cut -f 2,4,6,8 | sort -k 1,1) > data_ref.ens

seq $MIN $MAX | parallel "ln -s $SRC/FASTQ/*_{}_1.fastq.gz data_raw.{}_1.fastq.gz"
seq $MIN $MAX | parallel "ln -s $SRC/FASTQ/*_{}_2.fastq.gz data_raw.{}_2.fastq.gz"
seq $MIN $MAX | parallel 'fastqc -t 2 data_raw.{}_?.fastq.gz'

seq $MIN $MAX | parallel "cutadapt -b $IDX -b $IDX1 -b $IDX2 -o data_tr1.{}_1.fastq.gz data_raw.{}_1.fastq.gz > data_tr1.{}_1.log"
seq $MIN $MAX | parallel "cutadapt -b $PRM -O 10 --discard -o data_tr2.{}_2.fastq.gz -p data_tr2.{}_1.fastq.gz data_raw.{}_2.fastq.gz data_tr1.{}_1.fastq.gz > data_tr2.{}.log"
seq $MIN $MAX | parallel 'fastqc -t 2 data_tr2.{}_?.fastq.gz'

STR="--runThreadN 24 --genomeDir $GEN --outStd SAM --outSAMtype BAM Unsorted"
seq $MIN $MAX | parallel -P 1 "STAR $STR --outFileNamePrefix data_aln.{}. --readFilesIn <(zcat data_tr2.{}_1.fastq.gz) <(zcat data_tr2.{}_2.fastq.gz)"
seq $MIN $MAX | parallel 'fastqc data_aln.{}.*.bam'
cat <(echo 'input,unique,multiple,many,mismatch,short,other')\
    <(seq $MIN $MAX | parallel -k "cat data_aln.{}.Log.final.out | grep 'reads' | egrep 'input|\%' | cut -f 2 | tr -d '%' | paste -d, -s") > data_aln.csv
Rscript star.r

seq $MIN $MAX | parallel "featureCounts -a $GTF -Q 40 --primary -p -B -C -o data_ftc.{}.csv data_aln.{}.*.bam &> data_ftc.{}.log"
cat <(cat data_ftc.${MIN}.csv.summary | sed '1d' | cut -f 1 | sed 's/.*_//g' | tr '[:upper:]' '[:lower:]' | paste -d, -s)\
    <(seq $MIN $MAX | parallel -k "cat data_ftc.{}.csv.summary | sed '1d' | cut -f 2 | paste -d, -s") > data_ftc.csv

seq $MIN $MAX | parallel 'samtools view -@6 -q 10 -f 3 -F 256 -b data_aln.{}.*.bam > data_flt.{}.bam'
seq $MIN $MAX | parallel 'fastqc data_flt.{}.bam'
seq $MIN $MAX | parallel "java -jar $PICARD CollectMultipleMetrics I=data_flt.{}.bam O=data_flt.{}"

seq $MIN $MAX | parallel "samtools view data_flt.{}.bam | htseq-count -s reverse -m intersection-strict - $GTF > data_gne.{}.csv 2> data_gne.{}.log"
Rscript gdf.r 2> data_gdf.log
