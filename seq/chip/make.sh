source env.sh

cat $S/*.csv | tr -d '"' | tr ',' $'\t' | grep -v '^Pool' | cut -f 1,2,4 | sort -u | ruby gend.rb > design.csv

curl -s $GEN | zcat > ensembl.fa
bwa index -p ensembl ensembl.fa
cat ensembl.fa | grep '^>' | sed 's/^>//g' | tr ' :' $'\t' | cut -f 1,8 > ensembl.sizes

cat design.csv | sed '1d' | cut -f 1,6- | parallel -P 1 --col-sep $'\t' 'bwa mem -t 20 ensembl <(zcat {2} {3}) | samtools sort -@8 - {1}'
cat design.csv | sed '1d' | cut -f 1 | parallel 'samtools index {}.bam; samtools stats {}.bam > {}.txt'

ruby gent.rb
ls *.thor | parallel -P 1 'rgt-THOR -n {.} -m {}'

ruby genm.rb | parallel -P 5 --tmpdir . --colsep $' ' 'macs2 callpeak {} --bdg --trackline --SPMR'
