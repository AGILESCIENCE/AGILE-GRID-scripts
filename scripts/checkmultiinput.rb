#! /usr/bin/ruby
#0) file name input
#1) file name output (optional): if specified, check the fixflag and rewrite the file
#2) fixflag of analysis

load ENV["AGILE"] + "/scripts/conf.rb"

filename = ARGV[0]

datautils = DataUtils.new

if ARGV[0].to_s == "help" || ARGV[0] == nil || ARGV[0] == "h"
	system("head -3 " + $0 );
	exit;
end

if ARGV[1] != nil
	fout = File.new(ARGV[1], "w")
end

if ARGV[2] != nil
	ff = ARGV[2]
else
	ff = 1
end

index1 = 0

File.open(ARGV[0]).each do | line1 |
	index1 = index1 + 1
	noinc = false
	lll1 = line1.split(" ")
	l1 = lll1[1]
	b1 = lll1[2]
	name = lll1[6]
	fixflag = lll1[4]
	dist = "0.0"
	if lll1.size() == 8
		dist = "0.0"
	end
	fixflag = ff.to_s
	if name[0] == 49 || name[0] == "_"
		fixflag = "0"
	end
	index2 = 0
	File.open(ARGV[0]).each do | line2 |
		index2 = index2 + 1
		if line1 != line2 && index2 > index1
				
				lll2 = line2.split(" ")
 			 	
			 	l2 = lll2[1]
			 	b2 = lll2[2]
			 	
				d = datautils.distance(l1, b1, l2, b2)
				
				if(d.to_f < 1)
					puts line1.chomp + " - " +  line2
					noinc = true
				end
				if lll1[6] == lll2[6]
					puts "Found duplicated" + line1.chomp + " - " +  line2
					noinc = true
				end
				
				
		end
	end
	
	if noinc == false
		fout.write(lll1[0] + " " + lll1[1] + " " + lll1[2] + " " + lll1[3] + " " + fixflag + " " + lll1[5] + " " + lll1[6] + " " + dist + "\n")
	end
end




