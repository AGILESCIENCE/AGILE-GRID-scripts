#! /usr/bin/ruby
#script compatible with BUILD20
#0) filter DIR (es: FT3ab_2_I0007, FM3.119_2c, F4_2c_I0010)
#1) output file name prefix
#2) exp map
#3) skytype: 0 standard hires diffuse maps, 1 gc_allsky maps, 2 lowres
#optional
#emin
#emax

load ENV["AGILE"] + "/scripts/conf.rb"

datautils = DataUtils.new
parameters = Parameters.new
fits = Fits.new

if ARGV[0].to_s == "help" || ARGV[0] == nil || ARGV[0] == "h"
	system("head -10 " + $0 );
	exit;
end

filter = ARGV[0];
gas = ARGV[1];

exp = ARGV[2];

skytype = ARGV[3];

fits.readFitsHeader(exp);

datautils.extractFilterDir(filter)
filterdir = datautils.filterdir
filterall = filterdir;

filterbase2 = filter.split("_")[0] + "_" + filter.split("_")[1];

parameters.processInput(4, ARGV, filter)

emin1 = fits.minenergy;
emax1 = fits.maxenergy;
datautils.getSkyMatrix(filter, emin1, emax1, skytype)
skymap =  datautils.skymatrix;
skymapH = datautils.skymatrixH;
skymapL = datautils.skymatrixL;
puts "Sky map H: " + skymapH.to_s;
puts "Sky map L: " + skymapL.to_s;

cmd = "cp " + PATH + "share/AG_gasmapgen5.par . "
datautils.execute("", cmd);
cmd = "export PFILES=.:$PFILES; " + PATH + "bin/AG_gasmapgen5 " + exp.to_s + " " + gas.to_s + " " + skymapL.to_s + " " + skymapH.to_s;
datautils.execute("", cmd);
cmd = "rm AG_gasmapgen5.par"
datautils.execute("", cmd);
		
				


