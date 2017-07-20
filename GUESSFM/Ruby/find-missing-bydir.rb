#!/usr/bin/ruby

jobs=[]
# BASEDIR = "/chiswick/data/ncooper/imputation/COMBINED/TMP/"
BASEDIR = ARGV[0]
## ensure ends in '/'
if BASEDIR[-1,1] != "/"
  BASEDIR = BASEDIR + "/"
end

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

puts wantedfiles


# require '/home/cew54/bin/Qsub.rb'

# q = Qsub.new("impute.sh", :account=>"TODD-SL2", :time=>"24:00:00")
# jobs.each do |j|
#   q.add(j)
# end
# q.close()
