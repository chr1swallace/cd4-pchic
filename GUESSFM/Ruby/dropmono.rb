#!/usr/bin/ruby

# Search a directory for impute output files, and step through them,
# deleting SNPs which have expected MAF==0, and creating new output
# and info files:
#
# { root/*.gen.gz.gz,     root/*.gen.gz_info } ->
# { root/*.nomono_gen.gz, root/*.nomono_info }

require 'zlib'

load('Ruby/DIRS.rb')

# Given impute genotype file, corresponding info file
def ifile(f)
  f.sub('.gen.gz.gz','.gen.gz_info')
end

# Given impute genotype file, corresponding nomono info file
def oifile(f)
  f.sub('.gen.gz.gz','.nomono_info')
end

# Given impute genotype file, corresponding nomono genotype file
def odfile(f)
  f.sub('.gen.gz.gz','.nomono_gen.gz')
end

files=Dir.glob(IMPUTEDIR + '/imputed-*.gen.gz.gz')

files.each do |f|

  next if File.exist?(oifile(f))

  puts "dropping monomorph SNPs from #{f} -> " + oifile(f)
  
  ## infiles
  fi=File.open(ifile(f))
  fd=open(f)

  ## outfiles
  ofi=File.open(oifile(f), "w")
  ofd=File.open(odfile(f), "w")

  ## copy header of fi to ofi
  info=fi.readline
  ofi.write(info)

  ## read
  gz = Zlib::GzipReader.new(fd)
  ogz = Zlib::GzipWriter.new(ofd)
  gz.each_line do |line|
    info=fi.readline
    ss=info.split(/\s/)
    if ss[3] != "0.000" then
      ogz.write(line)
      ofi.write(info)
    end
  end


  ## close files
  ogz.close
  gz.close
  fi.close
  ofi.close

end
