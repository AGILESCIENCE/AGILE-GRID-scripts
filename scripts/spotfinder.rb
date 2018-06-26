#! /usr/bin/ruby
#0) outputfilename
#1) input map (int o cts) 
#-----
#2) spot finder: use cts (0) or int (1) 
#3) spot finder: smooth (e.g. 2)
#4) spot finder: number of source to be found (e.g. 10)
#5) spot finder: remove spots too near (radious), ( e.g. 0.7)
#6) spot finder: type of algorithm (optional, default 1) 0 original, 1 new with baricenter calculation, 2 new without baricenter calculation e.g. (1)
#7) min exposure (e.g. 50)
#-----
#8) mapsize (e.g. 120)
#9) exp map 

# spotfinder 0 2 10 0.7 1 50

load ENV["AGILE"] + "/scripts/conf.rb"
datautils = DataUtils.new
fits = Fits.new

if ARGV[0].to_s == "help" || ARGV[0] == nil || ARGV[0] == "h"
	system("head -16 " + $0 );
	exit;
end

outputfilename = ARGV[0]
map = ARGV[1]

fits.readFitsHeader(map.to_s);

binsize= fits.binsize


useint = ARGV[2]
smooth = ARGV[3]
numberofsource = ARGV[4]
radiousremove = ARGV[5]
algspotfiner = ARGV[6]
minExposure = ARGV[7]

mapsize = ARGV[8]

exp = ARGV[9]

#1) esecuzione di AG_spotfinder e generazione delle liste
strmt = "";
mapci = map;  
if useint.to_i == 1
	strmt = "I";
else
	strmt = "C";
end

	#listsource = prefix.to_s + "." + strmt.to_s + "sm" + format("%02d", smooth.to_i) + "ns" + format("%02d", numberofsource.to_i) + "al" + format("%02d", algspotfiner.to_i) + "rr" +  radiousremove.to_s + ".spotmulti";
	
	listsource = outputfilename;
	
	numberofsource2 = (numberofsource.to_i / 2).to_i;
	makebinshift = 0;
	ranal = 10;
	rextract = mapsize.to_i / 2.0 - ranal.to_i;
	
	cmd = PATH + "bin/AG_spotfinder5 " + mapci.to_s + " " + binsize.to_s + " " + smooth.to_s + " " + numberofsource2.to_s + " " + listsource.to_s + " " + algspotfiner.to_s + " " + radiousremove.to_s + " 1 " + makebinshift.to_s + " " + rextract.to_s + " " + exp.to_s + " " + minExposure.to_s;
	puts cmd; system(cmd);
	cmd = PATH + "bin/AG_spotfinder5 " + mapci.to_s + " " + binsize.to_s + " " + smooth.to_s + " " + numberofsource.to_s + " " + listsource.to_s + " " + algspotfiner.to_s + " " + radiousremove.to_s + " 2 " + makebinshift.to_s + " " + rextract.to_s + " " + exp.to_s + " " + minExposure.to_s;
	puts cmd; system(cmd);
	cmd = PATH + "bin/AG_spotfinder5 " + mapci.to_s + " " + binsize.to_s + " " + smooth.to_s + " " + numberofsource2.to_s + " " + listsource.to_s + " " + algspotfiner.to_s + " " + radiousremove.to_s + " 3 " + makebinshift.to_s + " " + rextract.to_s + " " + exp.to_s + " " + minExposure.to_s;
	
	puts cmd; system(cmd);

	
