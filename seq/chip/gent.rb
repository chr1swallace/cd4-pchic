def listf(type, sample, cond, batch)
  [type].product(sample).product([cond]).product(batch).map { |e| e.join('_') + '.bam' } 
end

def build(type, sample, condition, batch)
  { :rep1 => listf(type, sample, condition[0], batch),
    :rep2 => listf(type, sample, condition[1], batch),
    :genome => 'ensembl.fa',
    :chrom_sizes => 'ensembl.sizes',
    :inputs1 => listf('Input', sample, condition[0], [1] * batch.size),
    :inputs2 => listf('Input', sample, condition[1], [1] * batch.size) }
end

def pexp(type, sample, condition, batch)
  e = build(type, sample, condition, batch)
  f = [:rep1, :rep2, :inputs1, :inputs2].map { |k| e[k] }.reduce(&:+)
  if f.all? { |i| File.exists?(i) }
    n = [type, batch.join].join('_') + '.thor'
    File.open(n, 'w') { |f| f.puts(e.to_a.flatten.map { |e| e.is_a?(Symbol) ? '#' + e.to_s : e }.join("\n")) }
  end
end

d = File.read('design.csv').each_line.map { |l| l.chomp.split("\t") }
f = d.slice(1..-1).map { |e| d[0].map(&:to_sym).zip(e).to_h }
t = f.map { |e| e[:type] }.uniq - ["Input"]
s = f.map { |e| e[:sample] }.uniq
c = f.map { |e| e[:condition] }.uniq
t.product([[1], [2], [1,2]]).map { |t,b| pexp(t, s, c, b) }
