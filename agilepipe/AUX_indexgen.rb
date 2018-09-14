#! /usr/bin/ruby

# Copyright (c) 2016, AGILE team
# Authors: Andrea Bulgarelli <bulgarelli@iasfbo.inaf.it>,
#          Andrea Zoli <zoli@iasfbo.inaf.it>
#
# Any information contained in this software is property of the AGILE TEAM
# and is strictly private and confidential. All rights reserved.

ROOT_DATA = ENV["PATH_DATA"]
ROOT_DATABASE = "DATA_" + ENV["ARCHIVE"] + "/"

#AUX_PATH="/AGILE_DATA/AUX/"
AUX_PATH= ROOT_DATA + ROOT_DATABASE + "AUX/";

#EARTH_NAME="/AGILE_DATA/INDEX/EARTH.index"
EARTH_NAME= ROOT_DATA + ROOT_DATABASE + "INDEX/tmp_EARTH.index"

#SAS_NAME="/AGILE_DATA/INDEX/SAS.index"
SAS_NAME= ROOT_DATA + ROOT_DATABASE + "INDEX/tmp_SAS.index"

#ACSMAN_NAME="/AGILE_DATA/INDEX/ACSMAN.index"
ACSMAN_NAME= ROOT_DATA + ROOT_DATABASE + "INDEX/tmp_ACSMAN.index"

# EARTH searching
Dir.chdir(AUX_PATH)
dir0 = Dir[AUX_PATH+"*.EARTH"]
findex = File.new(EARTH_NAME,"w");
dir0.sort.each do | nomefile |
	File.open(nomefile).each_line do | fileline |
		a = fileline.split(" ");
		if a[0] == "<Process_TR>"
			starttime =  a[1] + " " + a[2];
			endtime = a[3] + " " + a[4];
			bb = nomefile + " " + starttime + " " + endtime;
			puts "ADDED"
			puts bb;
			findex.write(bb);
			findex.write("\n");
		end
	end
end
findex.close();

# SAS searching
Dir.chdir(AUX_PATH)
dir0 = Dir[AUX_PATH+"*.SAS"]
findex = File.new(SAS_NAME,"w");
dir0.sort.each do | nomefile |
	File.open(nomefile).each_line do | fileline |
		a = fileline.split(" ");
		if a[0] == "<Process_TR>"
			starttime =  a[1] + " " + a[2];
			endtime = a[3] + " " + a[4];
			bb = nomefile + " " + starttime + " " + endtime;
			puts "ADDED"
			puts bb;
			findex.write(bb);
			findex.write("\n");
		end
	end
end
findex.close();

# ACSMAN searching
Dir.chdir(AUX_PATH)
dir0 = Dir[AUX_PATH+"*.ACSMAN"]
findex = File.new(ACSMAN_NAME,"w");
dir0.sort.each do | nomefile |
	File.open(nomefile).each_line do | fileline |
		a = fileline.split(" ");
		if a[0] == "<Process_TR>"
			starttime =  a[1] + " " + a[2];
			endtime = a[3] + " " + a[4];
			bb = nomefile + " " + starttime + " " + endtime;
			puts "ADDED"
			puts bb;
			findex.write(bb);
			findex.write("\n");
		end
	end
end
findex.close();

system("sync");
