#! /bin/env ruby

require 'getoptlong'

##########################################################################
input = nil
output = nil
attribute_prefix = 'ID'
feature = nil
output_lines = Array.new

opts = GetoptLong.new(
  ['-i', '--in', '--input', GetoptLong::REQUIRED_ARGUMENT],
  ['-o', '--out', '--output', GetoptLong::REQUIRED_ARGUMENT],
  ['--intput-format', GetoptLong::REQUIRED_ARGUMENT],
  ['--output-format', GetoptLong::REQUIRED_ARGUMENT],
  ['--feature', GetoptLong::REQUIRED_ARGUMENT],
  ['--attribute', '--attribute_prefix', GetoptLong::REQUIRED_ARGUMENT],
  ['--chr_prefix', GetoptLong::REQUIRED_ARGUMENT],
)

opts.each do |opt,value|
  case opt
    when '-i', '--in', '--input'
      input = value
    when '-o', '--out', '--output'
      output = value
    when '--attribute', '--attribute_prefix'
      attribute_prefix = value
    when '--feature'
      feature = value
    end
end

raise "input has not been given! Exiting ......" if input.nil?

##########################################################################
File.open(input, 'r').each_line do |line|
  line.chomp!
  next if line =~ /^#/
  line_array = line.split("\t")
  strand = ""
  feature = "" if feature.nil?

  chr, gene_name, start, stop = line_array[0,4]
  if line_array.size == 5
    strand = line_array[4]
    if strand == "1" or strand == "-1"
      strand = strand=="1" ? "+" : "-" 
    end
  end

  attribute = attribute_prefix + "=" + gene_name
  output_lines << [chr, '', feature, start, stop, '', strand, '', attribute].join("\t")
end


if output
  out_fh = File.open(output, 'w')
  output_lines.each do |line| out_fh.puts line end
  out_fh.close
else
  output_lines.each{|line|puts line}
end

