#! /usr/bin/ruby
#0) l suource
#1) b source

load ENV["AGILE"] + "/scripts/conf.rb"

if ARGV[0].to_s == "help" || ARGV[0] == nil || ARGV[0] == "h"
	system("head -3 " + $0 );
	exit;
end

l = ARGV[0].to_f
b = ARGV[1].to_f

datautils = DataUtils.new

min = 9999
dirmin = ""
Dir["0*.*"].each do | dirs |
	lc = dirs.split(".")[0].to_f
	bc = dirs.split(".")[1].to_f
	puts lc.to_f
	d = datautils.distance(lc, bc, l, b)
	if d.to_f < min.to_f
		min = d
		dirmin = dirs
	end
end
puts format("%.2f", min) + " " + dirmin.to_s

