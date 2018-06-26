#! /usr/bin/ruby
#0) catalog (~/grid_scripts2/catalogs/integral_lmxb.multi) - .multi format
#1) l center
#2) b center
#3) outfile
#4) 0) d0 distance from (l,b)
#5) 0) fixflag for source with dist <= d0
#6) 1) d1 distance from (l,b)
#7) 1) fixflag for source with d0 < dist <= d1
#8) 2) d2 distance from (l,b)
#9) 2) fixflag for source with d1 < dist <= d2
#10) radius of search (last column of the .multi) - 0 do not specify
#11) min flux e.g. 25e-08
#12) min radius e.g. 1 -> do not select source with a distance less than minradius

load ENV["AGILE"] + "/scripts/DataUtils.rb"

if ARGV[0].to_s == "help" || ARGV[0] == nil || ARGV[0] == "h"
		system("head -14 " + $0 );
		exit;
	end

catal = ARGV[0];
l = ARGV[1]
b = ARGV[2]
out = ARGV[3];

d0 = ARGV[4]
fixflag0 = ARGV[5]
d1 = ARGV[6]
fixflag1 = ARGV[7]
d2 = ARGV[8]
fixflag2 = ARGV[9]
radius = ARGV[10]
minflux = ARGV[11]
minradius = ARGV[12]


datautils = DataUtils.new
input = catal;
outfile = File.new(out, "w");

File.open(input).each_line do |x|
	
	a = x.split(" ");
	
	if a[0].to_f < minflux.to_f
		next
	end

	d = datautils.distance(a[1], a[2], l, b);
	
	if minradius.to_f > 0
			if d.to_f != 0
				if d.to_f <= minradius.to_f
					next
				end
			end
	end
	
	fixflag=-1;
	
	if d.to_f <= d0.to_f
		fixflag = fixflag0;
	end
	if d.to_f > d0.to_f and d.to_f <= d1.to_f
		fixflag = fixflag1;
	end
	if d.to_f > d1.to_f and d.to_f <= d2.to_f
		fixflag = fixflag2;
	end
	
	outline = "";
	outline = a[0].to_s;
	outline += " ";
	outline += a[1].to_s;
	outline += " ";
	outline += a[2].to_s;
	outline += " "
	outline += a[3].to_s;
	outline += " "
	outline += fixflag.to_s; 
	outline += " ";
	outline += a[5].to_s;
	outline += " ";
	if fixflag.to_i == 0
		outline += "_";
	end
	outline += a[6].to_s;
	outline += " ";
	if fixflag.to_i >= 3
		outline += radius.to_f.to_s;
	else
		outline += "0.0"
	end
	
	for iii in 8...a.size
		outline += " ";
		outline += a[iii].to_s;
	end
		
	outline += "\n";
		
	if(fixflag.to_i != -1)
		outfile.write(outline.to_s);
		puts outline;
	end
	

end
outfile.close();
