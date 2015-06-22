#! /bin/env ruby

PGDD_file = ARGV[0]
gff_file = ARGV[1]

gff_info = Hash.new
pgdd_info = Hash.new


##################################################################3
File.open(PGDD_file,'r').each_line do |line|
  line.chomp!
  gene_name = line.split(',')[4]
  pgdd_info[gene_name] = ''
end

File.open(gff_file, 'r').each_line do |line|
  line.chomp!
  feature, attributes = line.split("\t").values_at(2,8)
  next if feature != 'gene'
  gene_name = nil
  if attributes =~ /Name=([^;]+)/
    gene_name = $1
    gff_info[gene_name] = ''
  end
end

puts pgdd_info.size
pgdd_info.delete_if{|h,k|! gff_info.include?(h)}
puts pgdd_info.size

