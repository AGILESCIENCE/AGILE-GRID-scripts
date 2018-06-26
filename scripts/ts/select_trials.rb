#0 input file name
#1 l center
#2 b center
#3 radius
#4 output filename

load "~/grid_scripts2/conf.rb"
load "~/grid_scripts2/MultiOutput.rb"

if ARGV[0].to_s == "help" || ARGV[0] == nil || ARGV[0] == "h"
	system("head -6 " + $0 );
	exit;
end


d = DataUtils.new
l1 = ARGV[1].to_f
b1 = ARGV[2].to_f
radius = ARGV[3].to_f
fout = File.new(ARGV[4], "w")

File.open(ARGV[0]).each_line do | line  |
	ll = line.split(" ")
	l = ll[1].to_f
	b = ll[2].to_f
	dist  = d.distance(l, b, l1, b1)
	
	if dist.to_f < radius.to_f
		fout.write(line)
	end
	
end
fout.close()