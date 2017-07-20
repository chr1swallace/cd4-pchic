h = %w(id type sample condition batch lane1 lane2)
d = ARGF.each_line.map { |l| l.chomp.split("\t") }.map { |e| e[0..-2] + e[-1].split('_') - %w(TC) }
f = d.map do |e| 
  b = e[0].match(/\-8/) ? '1' : '2'
  [[e[5], e[2], e[3], b].join('_'), e[5], e[2], e[3], b, *Dir["#{ENV['S']}/#{e[0]}.#{e[1]}.*.fq.gz"].sort]
end
puts ([h] + f.sort).map { |l| l.join("\t") }.join("\n")
