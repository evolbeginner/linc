#! /bin/env ruby

require 'getoptlong'

################################################################
mcscan_file = nil

opts = GetoptLong.new(
  ['-i', '--in', GetoptLong::REQUIRED_ARGUMENT],
)

opts.each do |opt,value|
  case opt
    when '-i', '--in'
      mcscan_file = value
  end
end

raise "input file has not been given by '-i|--in'!" if mcscan_file.nil?

################################################################
fh = File.open(mcscan_file, 'r')
fh.each_line do |line|
  line.chomp!
  next if line =~ /^#/;
  # 0-  0:	AT1G02130	470156	  1e-49
  block_full_name, paralog1, paralog2 = line.split("\t")[0,3]
  if block_full_name =~ /(\d+)/
    block = $1
  else
    STDERR.puts "Warning: block name error at #{block_full_name}"
    next
  end
  puts [block, "", "", paralog1, paralog2, "", ""].join(",")
end
fh.close

