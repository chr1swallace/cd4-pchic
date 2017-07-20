ARGF.lines.each do |e|
  l = e.split("\t")
  b = l.first
  l.last.split(';').each { |f| puts [b,f].join("\t") }
end
