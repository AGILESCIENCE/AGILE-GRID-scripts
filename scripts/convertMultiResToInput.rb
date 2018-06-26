#! /usr/bin/ruby
#0) input file name
#1) output file name
#2) analysis type (optional, default 1)
#3) minsqrtTS to cut (optional, default 0)
#4) radious search for alike (optional, default 0)
#5) max off-axis (optional, if appliable ( column 9 is present), default 90)
#6) fixflag=0 for sources that starts with _

load ENV["AGILE"] + "/scripts/conf.rb"

if ARGV[0].to_s == "help" || ARGV[0] == nil || ARGV[0] == "h"
	system("head -9 " + $0 );
	exit;
end

fileinp1 = ARGV[0]
fileout1 = ARGV[1]
if ARGV[2] != nil
	ant = ARGV[2]
else
	ant = 1
end
if ARGV[3] != nil
	mins = ARGV[3]
else
	mins = 0
end
if ARGV[4] != nil
	radiuossearch = ARGV[4]
else
	radiuossearch = 0
end
if ARGV[5] != nil
	maxoff = ARGV[5]
else
	maxoff = 90
end
if ARGV[6] != nil
	flag = ARGV[6]
else
	flag = 0
end

alikeutils = AlikeUtils.new

alikeutils.convertMultiResToMulti(fileinp1, fileout1, ant, 2.0, mins, maxoff, radiuossearch);

if flag.to_i != 0
        fout = File.new(ARGV[1] + ".tmp", "w")
        File.open(ARGV[1]).each do | line1 |
                lll1 = line1.split(" ")
                name = lll1[6]
                fixflag = lll1[4]
                dist = "0.0"
                if lll1.size() == 8
                        dist = lll1[7]
                end
                if name[0] == 49 || name[0] == "_"
                        fixflag = "0"
                end
                fout.write(lll1[0] + " " + lll1[1] + " " + lll1[2] + " " + lll1[3] + " " + fixflag + " " + lll1[5] + " " + lll1[6] + " " + dist + "\n")
        end
        fout.close()
        system("mv " + ARGV[1] + ".tmp " + ARGV[1])
        system("cat " + ARGV[1]);
end




