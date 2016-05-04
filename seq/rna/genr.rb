i = ARGF.each_line.lazy
f = i.map { |l| l.chomp.split("\t") }.each_with_index
d = f.map do |l,i|
  e = "CD4R" + i.to_s.rjust(11, "0")
  [l[0], 'regulatory', 'gene', *l[1..2], '.', '.', '.', "gene_id \"#{e}\"; gene_name \"#{e}\"; gene_biotype \"regulatory\""]
end
d.each do |e|
  puts [*e[0..5], '+', *e[7..-1]].join("\t")
  puts [*e[0..5], '-', *e[7..-1]].join("\t")
end
