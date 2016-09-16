#!/usr/bin/ruby
require 'optparse'

INPUTDIR = "./FM-impute/"
COMMONDIR = "./1000GP_Phase3/"
SHAPEITDIR = "./FM-phased/"

OPTIONS = {}
OptionParser.new do |opts|
  opts.banner = "Usage: local-shapeit-impute.rb [OPTIONS]"

  opts.on("-v", "--[no-]verbose", "Run verbosely") do |v|
    OPTIONS[:verbose] = v
  end
  opts.on("-q", "--[no-]queue", "Use qCom.sh in output") do |v|
    OPTIONS[:queue] = v
  end
  opts.on("-n", "--nohup", "Prepend with nohup") do |v|
    OPTIONS[:nohup] = v
  end
  opts.on("-b", "--[no-]background", "Append with &") do |v|
    OPTIONS[:background] = v
  end
end.parse!

def file_in(chrarm)
  INPUTDIR + "impute-#{chrarm}.gen.gz"
end
def file_phased(chrarm,rstart,rend)
  INPUTDIR + "phased-#{chrarm}-#{rstart}-#{rend}"
end
def file_imputed(chrarm,rstart,rend)
  INPUTDIR + "imputed-#{chrarm}-#{rstart}-#{rend}.gen.gz"
end
def file_exclude(chrarm,rstart,rend)
  INPUTDIR + "exsamples-#{chrarm}-#{rstart}-#{rend}"
end
def command_shapeit(ss) 
  chr=ss[1]
  chrarm=ss[8]
  rstart=ss[2]
  rend=ss[3]
  wstart=Integer(rstart) - 50000
  wend=Integer(rend) + 50000
  ifile=file_in(chrarm)
  sfile=ifile.sub(".gen.gz",".sample")
  ofile=file_phased(chrarm,rstart,rend)
  efile=file_exclude(chrarm,rstart,rend)
  estr = if File.exist?(efile)
           "--exclude-ind #{efile} "
         else
           ""
         end
  map = COMMONDIR + "genetic_map_chr" + chr + "_combined_b37.txt"
  threads=12
  threads=1 if OPTIONS[:queue]
  comm="/home/chrisw/local/bin/shapeit --thread #{threads} --window 0.5 --input-gen #{ifile} #{sfile} -M #{map} -O #{ofile} --input-from #{wstart} --input-to #{wend} --output-from #{rstart} --output-to #{rend} #{estr} --output-log #{ofile}.log"
  wrap_command(comm)
end
def command_impute(ss) 
  chr=ss[1]
  chrarm=ss[8]
  rstart=ss[2]
  rend=ss[3]
  pfile=file_phased(chrarm,rstart,rend)
  ofile=file_imputed(chrarm,rstart,rend)
  map = COMMONDIR + "genetic_map_chr" + chr + "_combined_b37.txt"
  haps= COMMONDIR + "1000GP_Phase3_chr" + chr + ".hap.gz"
  legend = COMMONDIR + "1000GP_Phase3_chr" + chr + ".legend.gz"
  comm="/home/chrisw/local/bin/impute2 -filt_rules_l 'EUR<0.001' -use_prephased_g -m #{map} -h #{haps} -l #{legend} -known_haps_g #{pfile}.haps -int #{rstart} #{rend} -Ne 20000 -o #{ofile} -o_gz"
  wrap_command(comm)
end

def wrap_command(comm) 
  comm = "qCom.sh " + comm   if OPTIONS[:queue]
  comm = "nohup " + comm   if OPTIONS[:nohup]
  comm = comm + " &"   if OPTIONS[:background]
  comm
end    
def check_shapeit_complete(file)
  logfile = file.sub(".haps",".log")
  File.open(logfile) do |f|
    f.each_line do |line|
      return file if line.include?('Running time')
    end
  end
  return nil
end
# def logfile(obase,i)
#   obase + "-#{i}.log"
# end

## read regions
fregions="regions.csv"
lines = File.readlines(fregions)
## drop first line, corresponds to headings
lines.shift()

## existing input files
infiles = Dir.glob(INPUTDIR + 'impute-*.gen.gz')
puts "INFILES" if OPTIONS[:verbose]
puts infiles if OPTIONS[:verbose]

## existing phased files
pfiles = Dir.glob(INPUTDIR + 'phased-*.haps')
pfiles = pfiles.map { |file| check_shapeit_complete(file) }
pfiles = pfiles.compact
puts "PHASEDFILES" if OPTIONS[:verbose]
puts pfiles if OPTIONS[:verbose]

## existing imputed files
ifiles = Dir.glob(INPUTDIR + 'imputed-*.gen.gz.gz')
puts "IMPUTEDFILES" if OPTIONS[:verbose]
puts ifiles if OPTIONS[:verbose]

scommands = []
icommands = []
shaped = []

#line = lines[40]
lines.each do |line|
  ss = line.split(/\t/)
   # command = "\n## " + line 
    chr=ss[1]
    chrarm=ss[8]
    rstart=ss[2]
    rend=ss[3]
  
  ## does the input file exist?
  ifile = file_in(chrarm)
  next unless (infiles.include? ifile) 

  ## does the shapeit file exist?
  pfile = file_phased(chrarm,rstart,rend) + ".haps"
  if ( pfiles.include? pfile ) 
    ## run impute
    next if (ifiles.include?(file_imputed(chrarm,rstart,rend)+".gz") )
    icommands.push( command_impute(ss) )
  else
    ## run shapeit
#    shaped.push(chrarm)
    scommands.push( command_shapeit(ss) )
  end
  
end

puts scommands
puts
puts icommands

