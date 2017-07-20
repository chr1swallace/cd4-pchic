cat <(echo $'baitID\toeID\tactivated\tnonactivated') <(cat raw-hic.csv | grep -v '^#' | sed '1d' | cut -f 4,9,23,24) > hic.csv
cat <(echo $'oeID\tgene') <(cat raw-hic.csv | grep -v '^#' | cut -f 9,10 | awk '$2 != "."' | sed '1d' | ruby melt.rb) > oeg.csv
cat raw-rna.csv | cut -f 2,4,6,10- > rna.csv
