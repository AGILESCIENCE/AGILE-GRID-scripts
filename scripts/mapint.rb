#! /usr/bin/ruby
#script compatible with BUILD20
#0) filter DIR (es: FT3ab_2_I0007, FM3.119_2c, F4_2c_I0010)
#1) output file name 
#2) prefix input map

load ENV["AGILE"] + "/scripts/conf.rb"

datautils = DataUtils.new
parameters = Parameters.new

if ARGV[0].to_s == "help" || ARGV[0] == nil || ARGV[0] == "h"
	system("head -9 " + $0 );
	exit;
end

filter = ARGV[0];
int  = ARGV[1];

prefix  = ARGV[2];



datautils.extractFilterDir(filter)
filterdir = datautils.filterdir
filterall = filterdir;

filterbase2 = filter.split("_")[0] + "_" + filter.split("_")[1];




parameters.processInput(3, ARGV, filter)

exp = prefix.to_s + ".exp.gz"
cts = prefix.to_s + ".cts.gz"

cmd = PATH + "bin/AG_intmapgen5 " + exp.to_s + " " + int.to_s + " " + cts.to_s;
datautils.execute(prefix, cmd);

		
				


