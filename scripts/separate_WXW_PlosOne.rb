#! /bin/env ruby

pgdd_file = ARGV[0]
gff_file = ARGV[1]

gff_info = Hash.new{|h,k|h[k]=Hash.new}
pgdd_info = Hash.new


##################################################################3
File.open(gff_file, 'r').each_line do |line|
  line.chomp!
  chr, feature, attributes = line.split("\t").values_at(0,2,8)
  next if feature != 'gene'
  gene_name = nil
  if attributes =~ /Name=([^;]+)/
    gene_name = $1
    gff_info[gene_name]['chr'] = chr
  end
end


chr = ''
block_name = 0
File.open(pgdd_file,'r').each_line do |line|
  line.chomp!
  at_gene, brapa_gene = line.split(',').values_at(3,4)
  if gff_info[brapa_gene]['chr'] != chr
    block_name += 1
    chr = gff_info[brapa_gene]['chr']
  end
  puts [block_name,[nil]*2,at_gene,brapa_gene,[nil]*2].join(',')
end


