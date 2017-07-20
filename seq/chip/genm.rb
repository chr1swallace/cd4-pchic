require 'set'

BROAD = %w(H3K27me3 H3K36me3 H3K9me3 H3K4me1).to_set

d = File.read('design.csv').each_line.map { |l| l.chomp.split("\t") }
f = d.slice(1..-1).map { |e| d[0].map(&:to_sym).zip(e).to_h }
p = f.reject { |e| e[:type] == 'Input' }.map { |e| [e[:id] + '_macs', e[:id] + '.bam', 
                                                    ['Input', e[:sample], e[:condition], '1'].join('_') + '.bam', 
                                                   BROAD.member?(e[:type]) ? '--broad' : nil] }
puts p.map { |e| ['-n', '-t', '-c', nil].zip(e).map { |e| e.join(' ').gsub(/^ /, '') }.reject { |e| e.empty? }.join(' ') }
