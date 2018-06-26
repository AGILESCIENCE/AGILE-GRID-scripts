#! /usr/bin/ruby
#0) prefix input 1 of files to be merged
#1) prefix output
#2) filter
#2) emin (optional, default 100)
#3) emax (optional, default 50000)
#4) gal coeff
#5) iso coeff
#6) skytype : 0 SKY000-1 + SKY000-5, 1 gc_allsky maps + SKY000-5, 2 SKY000-5, 3 SKY001 (old galcenter, binsize 0.1, full sky), 4 SKY002 (new galcenter, binsize 0.1, full sky)

load ENV["AGILE"] + "/scripts/conf.rb"

datautils = DataUtils.new

if ARGV[0].to_s == "help" || ARGV[0] == nil || ARGV[0] == "h"
	system("head -10 " + $0 );
	exit;
end


prefix1 = ARGV[0];
prefix2 = ARGV[1];
filter = ARGV[2];

if ARGV[3] != nil
	emin1 = ARGV[3];
else
	emin1 = 100
end

if ARGV[4] != nil
	emax1 = ARGV[4];
else
	emax1 = 50000
end

if ARGV[5] != nil
        galcoeff = ARGV[5];
else
        galcoeff = -1
end

if ARGV[6] != nil
        isocoeff = ARGV[6];
else
        isocoeff = -1
end

if ARGV[7] != nil
        skytype = ARGV[7];
else
        skytype = 4
end

a = Dir[prefix1.to_s + "*.cts.gz"].sort
index = 0;
output = "temp" + index.to_s + ".cts.gz"
cmd = "cp " + a[0].to_s + " " + output.to_s
puts cmd
system(cmd);
output2=""
a.each do | file |
	if index == 0
		index = index + 1
		next
	end
	output = "temp" + (index-1).to_s + ".cts.gz"
	output2 = "temp" + (index).to_s + ".cts.gz"
	cmd = "fimgmerge " + file.to_s + " " + output.to_s + " " + output2.to_s + " 0 0 "
	puts cmd;
	system(cmd);
	index = index + 1
end
cmd = "mv " + output2.to_s + " " + prefix2.to_s + ".cts.gz"
system(cmd)

a = Dir[prefix1.to_s + "*.exp.gz"].sort
index = 0;
output = "temp" + index.to_s + ".exp.gz"
cmd = "cp " + a[0].to_s + " " + output.to_s
puts cmd
system(cmd);
output2=""
a.each do | file |
	if index == 0
		index = index + 1
		next
	end
	output = "temp" + (index-1).to_s + ".exp.gz"
	output2 = "temp" + (index).to_s + ".exp.gz"
	cmd = "fimgmerge " + file.to_s + " " + output.to_s + " " + output2.to_s + " 0 0 "
	puts cmd;
	system(cmd);
	index = index + 1
end
cmd = "mv " + output2.to_s + " " + prefix2.to_s + ".exp.gz"
system(cmd)

cmd = "rm temp*"
system(cmd);

datautils.getSkyMatrix(filter, emin1, emax1, skytype)
skymap =  datautils.skymatrix;
skymapL =  datautils.skymatrixL;
skymapH = datautils.skymatrixH;
#skymap =  format("%01d_%01d", emin1, emax1) + ".0.1.conv.sky ";

cmd = "cp " + PATH + "share/AG_gasmapgen5.par . "
system(cmd)
cmd = "export PFILES=.:$PFILES; " + PATH + "bin/AG_gasmapgen5 " + prefix2.to_s + ".exp.gz " + prefix2.to_s + ".gas.gz "  + skymapL.to_s + " " + skymapH.to_s;
puts cmd
system(cmd)

cmd = "cp " + PATH + "share/AG_intmapgen5.par . "
system(cmd)
cmd = "export PFILES=.:$PFILES; " + PATH + "bin/AG_intmapgen5 " + prefix2.to_s + ".exp.gz" + " " + prefix2.to_s + ".int.gz" + " " + prefix2.to_s + ".cts.gz";
puts cmd
system(cmd)

ff=File.new(prefix2.to_s + ".maplist4", "w")
ff.write(prefix2.to_s + ".cts.gz " + prefix2.to_s + ".exp.gz " + prefix2.to_s + ".gas.gz 30 " + galcoeff.to_s + " " + isocoeff.to_s)
ff.close() 
