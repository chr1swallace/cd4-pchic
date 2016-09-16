#!/usr/bin/ruby

jobs=[]
BASEDIR = "./FM-impute_qc/"

## read regions
fregions="regions.csv"
lines = File.readlines(fregions)
## drop first line, corresponds to headings
lines.shift()

def outfile(obase)
  BASEDIR + obase + ".RData"
end

## existing outfiles
outfiles = Dir.glob(BASEDIR + 'imputed-*.RData')
wantedfiles = []

lines.each do |line|
    ss = line.split(/\t/)
   # command = "\n## " + line
    chr=ss[1]
    chrarm=ss[8]
    rstart=ss[2]
    rend=ss[3]
    obase = "imputed-#{chrarm}" + "-" + rstart + "-" + rend

    wantedfiles.push(outfile(obase))
end

wantedfiles = wantedfiles - outfiles

wantedfiles=wantedfiles.map {|s| s.gsub(/_qc/, '')}
wantedfiles=wantedfiles.map {|s| s.gsub(/.RData/, '.gen.gz.gz')}

puts wantedfiles
