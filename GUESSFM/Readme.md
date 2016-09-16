# R packages

R scripts depend on some packages not yet on CRAN, but at https://github.com/chr1swallace - 
```
randomFunctions
annotSnpStats
GUESSFM
```


# GET 1000 GENOMES GENOTYPES INTO RDATA

get data

```
mkdir ./1000GP_Phase3_vcf
cd ./1000GP_Phase3_vcf
for chr in `seq 1 22`; do
  wget -nc ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
  wget -nc ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
done
wget -nc ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples.20130502.ALL.ped
```

reindex

```
cd ./1000GP_Phase3_vcf
for f in ALL.chr1*.vcf.gz ALL.chr[3456789]*.vcf.gz ALL.chr2[012]*.vcf.gz; do
  echo $f
  uf=`basename $f .gz`
  t=`tempfile -d ./tmp`
cat > $t << EOF
#!/bin/bash
gunzip $f
rm $f
/stats/oliver/bin/bgzip $uf 
/stats/oliver/bin/tabix -p vcf $f
EOF
echo $t
bash $t
done
```

## RUN SHAPEIT AND IMPUTE

first run runs shapeit
```
Ruby/local-shapeit-impute.rb -n -b
```

second time, if shapeit output exists, runs IMPUTE2
```
Ruby/local-shapeit-impute.rb
```

## Remove monomorphic SNPs from impute output
This makes reading files into R faster, easier on memory.

Will only run if output file does not exist.

`Ruby/dropmono.rb`


## QC imputed data

```
## to start over:
## rm ./FM-impute_qc/*
	
		Ruby/find-missing-impute-qc.rb | wc # want to be 188
	
		for f in `Ruby/find-missing-impute-qc.rb` ; do
			bf=`basename $f .gen.gz.gz`
			if=./FM-impute/${bf}
			Rscript R/join-qc-imputed.R --args file=${if} 
		done
	
		ls ./FM-impute_qc/imputed*RData | wc # 188 - want 188
```

## QC part 2 - LD with 1000 genomes

```
  Ruby/find-missing-bydir.rb ./1000Genomes-by-ichip-region/

  for f in `Ruby/find-missing-bydir.rb ./1000Genomes-by-ichip-region/`; do
      bf=`basename $f`
          Rscript R/qc-check-snps-by-1000genomes-ld-process-each-region.R \
          --args ifile=./FM-impute_qc/${bf}
  done

  Ruby/find-missing-bydir.rb ./1000Genomes-by-ichip-region/

  ## summarise
  Rscript R/qc-check-snps-by-1000genomes-ld-combined.R > log/qc-check-snps-by-1000genomes-ld.log 2>&1
```



# Stepwise results, for comparison
```
    for f in ./FM-impute_qc/imputed-*.RData; do
        region=`basename $f .RData | sed 's/imputed-//'`
        echo stepwise $region
        of=./FM-stepwise/$region
        [[ -e $of ]] || \
            ./R/stepwise.R --args file=$f > log/stepwise-${region}.log 2>&1
    done
```

# Run GUESSFM
```
for f in ./FM-impute_qc/imputed-*.RData; do
        region=`basename $f .RData | sed 's/imputed-//'`
        echo $region
        of=./FM-GUESS/$region
##        [[ -e $of ]] || \
            nohup ./R/prepguess.R --args file=$f > log/prepguess-${region}.log 2>&1
        njobs=`ps -afe | grep GUESS | wc -l`
        while [ "$njobs" -gt 10 ]; do sleep 120; njobs=`ps -afe | grep GUESS | wc -l`; done
    done
```

Checks that everything that ran was ok:

```
  ## 188 input regions
  ls ./FM-GUESS | wc -l

  ## 116 regions have at least one output file 
  find ./FM-GUESS -name out_50000_sweeps_output_best_visited_models.txt | \
  sed 's/[^\/]*$//'| \
  awk -F"/" '{print $5}'| sort | uniq -c|wc -l

  ## 237 output files
  find ./FM-GUESS -name out_50000_sweeps_output_best_visited_models.txt | wc -l

  ## failures with small output files: 0
  find ./FM-GUESS -name out_50000_sweeps_output_best_visited_models.txt -size -2b | \
  sed 's/\/out.*//' |wc -l

```


** QC checks on results

```
## delete any qc output files from an earlier run
      # find ./FM-GUESS/ -name qc2-ppnsnp.csv | xargs rm
      # find ./FM-GUESS/ -name qc2-snps.csv | xargs rm

      ## run QC
      for d in ./FM-GUESS/*; do
          bd=`basename $d`
	./R/GUESS-qc.R --args d=$d 
      done

      ## highlight any qc issues
      find ./FM-GUESS/ -name qc-ppnsnp.csv | xargs awk '$4=="TRUE"'
      find ./FM-GUESS/ -name qc-snps.csv | xargs wc

      ## examine manually, eg
      find ./FM-GUESS/ -name qc3-ppnsnp.csv | \
          xargs awk '$4=="TRUE" {
             gsub(/\"/,"",$5); 
           print;
           system("display ./FM-GUESS/" $5 "/qc3-plots.pdf");
          }' 

      ## generate a list of snps to remove and delete bad runs
      Rscript R/check-guess-qc.R
  find ./FM-GUESS/ -name qc2-snps.csv | xargs awk '{print $5}' | uniq
#+END_SRC

## Rerun QC failures

NB, this part needs to be run repeatedly until no unexplainable QC
failures.
The full set of SNPs excluded through this are to be found in the
files:

```
  ls snps-to-drop-*.txt
```

# Collect positions of genes to help plotting results
```
   for d in ./FM-GUESS/*; do
      echo $d
      region=`basename $d`
      [[ -e $d/genes.RData ]] || \
          Rscript R/get_genes.R --args d=$d > /dev/null 2>&1
  done
```

# Process GUESSFM results

```
    for d in ./FM-GUESS/*; do
      echo $d
      region=`basename $d`
     [[ -e $d/snpmod.RData ]] || [[ -e $d/skip ]] || \
          R/GUESS-expand-models.R --args d=$d > log/GUESS-expand-models-$region.Rout 2>&1
    done

    for d in ./FM-GUESS/*; do
      echo $d
      region=`basename $d`
      echo $region
      d=./FM-GUESS/$region
R/GUESS-resultplots.R --args d=$d 
R/GUESS-group-overlap.R --args d=$d 
  done
```

Wait for the above to complete!

```
for d in ./FM-GUESS/*; do
      echo $d
      region=`basename $d`
      echo $region
      d=./FM-GUESS/$region
         R/GUESS-haplotypes.R --args d=$d > log/GUESS-haplotypes-$region.Rout 2>&1
  done
```

## General plots
```
# rm ./FM/output/*/summx-* ./FM/output/*/manhattan-* 
  for d in `Ruby/find-expanded-guess.rb`; do
    echo $d
    region=`basename $d`
    [[ -e ./FM/output/$region/summx-bestsnps.pdf ]] || \
 R/GUESS-resultplots.R --args d=$d > log/GUESS-group-overlap-$region.Rout 2>&1
  done
```

## Plot networks of SNP groups
```
  for d in  `Ruby/find-expanded-guess.rb`; do
    echo $d
    region=`basename $d`
    [[ -e ./FM/output/$region/group-overlaps.pdf ]] || \
 R/GUESS-group-overlap.R --args d=$d > log/GUESS-group-overlap-$region.Rout 2>&1
  done
```

## Haplotype analyses

```

  for d in `Ruby/find-expanded-guess.rb`; do
    echo $d
    region=`basename $d`
    [[ -e ./FM/output/$region/haplotypes.pdf ]] || \
 R/GUESS-haplotypes.R --args d=$d > log/GUESS-haplotypes-$region.Rout 2>&1
  done
```




## Collect and collate summaries of analysis results

```
## Manhattan data
Rscript R/GUESSFM-save-manhattan.R
## MPPI
Rscript R/GUESS-save-mppi.R
## nsnp posteriors
Rscript R/GUESS-save-nsnp.R
```

## Compare to PMI

```
Rscript R/compare-GUESSFM-scv.R
```

## Gene ids
```
: mkdir ./FM
: mysql  --user=genome --host=genome-mysql.cse.ucsc.edu -A -D hg19 -e "select distinct G.gene,N.value from ensGtp as G, ensemblToGeneName as N where    G.transcript=N.name" > ./FM/ens2symbol.txt
```
