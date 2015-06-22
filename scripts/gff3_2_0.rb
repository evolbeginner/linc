#! /bin/env ruby

# Gff0 is defined as a format as follows:
# chr start stop gene_name [strand]
# Gff0 is used in some programmes as the specified formats for input files, e.g. MCScan.

require 'getoptlong'

########################################################################
## ----------- trash -------------- ##
class GFF
  attr_accessor :input, :input_format, :output_format

  def initialize()

  end

  def read_input(input, input_format=nil)
    raise "input not given" if input.nil?
    if input_format.nil?
      detect_format_result = detect_format(input)
      if not detect_format_result
        self.input_format = detect_format_result
      else
        raise "gff input #{input} format error!"
      end
    end
  end

  def detect_format(input)
    input_format = nil
    line_array = Array.new
    File.open(input).each_line do |line|
      next if line =~ /^\#/
      line_array = line.split("\t")
      break
    end
    if line_array.size == 4 or line_arrya.size == 5
      input_format = 'gff0'
    else
      start, stop = line_array.values_at(3,4)
      if start =~ /^[0-9]+$/ and stop =~ /^[0-9]+$/
        input_format = 'gff3'
      else 
        input_format = false
      end
    end
    return(input_format)
  end

  def read_input()
    gff_lines = Array.new
    File.open(input, 'r').each_line do |line|
      line.chomp!
      next if line =~ /^#/
      attrs = Array.new
      line_eles = Array.new
      eles = Hash.new

      case input_format
        when 'gff0'
          attrs = %w[chr name start stop strand]
          line_eles = line.split("\t")
        when 'gff3'
          attrs = %w[chr feature start stop strand attributes]
          line_eles = line.split("\t").values_at(0,2,3,4,6,8)
      end

      line.split("\t").each_with_index do |ele, index|
        eles[attrs[index]] = ele
      end

      gff_lines.push eles
    end
    return gff_lines
  end


  def convert_format()
    case output_format
      when 'gff0'
        
      when 'gff3'
        required_eles = %w[chr start stop strand attributes]
    end
  end

end
## ----------------------------------------------------- ##


########################################################################
input = nil
input_format = nil
output_format = nil
output = nil
attribute_types = Array.new
chr_prefix = nil
feature_hash=Hash.new

output_lines = Array.new

opts = GetoptLong.new(
  ['-i', '--in', '--input', GetoptLong::REQUIRED_ARGUMENT],
  ['-o', '--out', '--output', GetoptLong::REQUIRED_ARGUMENT],
  ['--intput-format', GetoptLong::REQUIRED_ARGUMENT],
  ['--output-format', GetoptLong::REQUIRED_ARGUMENT],
  ['--feature', GetoptLong::REQUIRED_ARGUMENT],
  ['--attribute', '--attributes', '--attribute_types', GetoptLong::REQUIRED_ARGUMENT],
  ['--chr_prefix', GetoptLong::REQUIRED_ARGUMENT],
)

opts.each do |opt,value|
  case opt
    when '-i', '--in', '--input'
      input = value
    when '-o', '--out', '--output'
      output = value
    when '--input-format'
      input_format = value
    when '--output-format'
      output_format = value
    when '--feature'
      feature_hash[value] = ''
    when '--attribute', '--attributes', '--attribute_types'
      value.split(',').each do
        attribute_types << value
      end
    when '--chr_prefix'
      chr_prefix = value
  end
end

raise "input has not been given! Exiting ......" if input.nil?

attribute_types = %w[Name ID] if attribute_types.empty?

########################################################################
File.open(input, 'r').each_line do |line|
  line.chomp!
  chr,feature,start,stop,strand,attributes = line.split("\t").values_at(0,2,3,4,6,8)

  if ! feature_hash.empty?
    next if not feature_hash.include? feature
  end

  attribute_types.each do |attribute_type|
    if attributes=~/#{attribute_type}=([^;]+)/
      attr = $1
      break
    end
  end
  chr = chr_prefix + chr if ! chr_prefix.nil?
  output_lines << [chr,attr,start,stop,strand].join("\t")
end


if output
  out_fh = File.open(output, 'w')
  output_lines.each do |line| out_fh.puts line end
  out_fh.close
else
  output_lines.each{|line|puts line}
end


