#! /usr/bin/ruby
#0) input file resfinal
#1) output file
#2) ring XXXX.YYYY

load ENV["AGILE"] + "/scripts/conf.rb"

if ARGV[0].to_s == "help" || ARGV[0].to_s == "h" || ARGV[0] == nil
	system("head -4 " + $0 );
	exit;
end

fout = File.new(ARGV[1], "w")

File.open(ARGV[0]).each_line do | line |
	
	l = line.split(" ")
	lout = l[1] + " " + l[2] + " " + l[4] + " " + l[5] + " " + l[15] + " " + l[16] + " " + l[21] + " " + l[22] + " " + l[29] + " " + l[30] + " " + l[7] + " " + l[8] + " " + l[10] + " " + l[11] + " " + l[12] + " " + l[13] + " " + l[27] + " " + ARGV[2] + "\n"
	fout.write(lout)
	
end
