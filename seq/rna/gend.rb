H = %w(id sample condition pair1 pair2)
N = %w(202 201 204 203)
S = %w(1 1 2 2)
C = %w(non act non act)

d = N.zip(S).zip(C).map(&:flatten).map { |e| [[e[1], e[2]].join('_'), e[1], e[2], "#{ENV['SRC']}/WTCHG_133843_#{e[0]}_1.fastq.gz", "#{ENV['SRC']}/WTCHG_133843_#{e[0]}_2.fastq.gz"] }
puts ([H] + d).map { |e| e.join("\t") }.join("\n")
