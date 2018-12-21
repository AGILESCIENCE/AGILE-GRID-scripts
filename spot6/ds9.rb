#! /usr/bin/ruby
#0) name input
#1) name output
#2) scale low (-1 dont scale) If scale high = -1 and this != -1 => scale squared
#3) scale high (-1 dont scale)
#4) smooth radius (default 3)
#5) cmap (default B)
#6) zoom (0 to fit, 4, 8 and so on)
#7) type of file output (png, jpg, default jpg)
#8) dimension of window (default 1100x1100)
#9) region file name to load
#10) additional region file
#11) additional region file 2
#12) additional command

load ENV["AGILE"] + "/scripts/conf.rb"
load ENV["AGILEPIPE"] + "/env.rb"

if ARGV[0].to_s == "help" || ARGV[0] == nil || ARGV[0] == "h"
	system("head -14 " + $0 );
	exit;
end

names = ARGV[0];
nameo = ARGV[1];

scale1 = -1
scale2 = -1

if( ARGV[2] != nil)
	scale1 = ARGV[2].to_f;
end
if( ARGV[3] != nil)
	scale2 = ARGV[3].to_f;
end
if( ARGV[4] != nil)
	smooth_radius = ARGV[4].to_i;
else
	smooth_radius = 3;
end

if( ARGV[5] != nil)
	cmap = ARGV[5];
 	#puts cmap
 	if cmap.include?("inv")
 		cmap = cmap.split("inv")[0] + " -cmap invert yes "
 	end
else
	cmap = "B";
end

if( ARGV[6] != nil)
        zzz = ARGV[6];
else
        zzz = "to fit";
end

if zzz.to_i == 0
	zzz = "to fit"
end

if( ARGV[7] != nil)
	typeout = ARGV[7];
else
	typeout = "jpg";
end

if( ARGV[8] != nil)
	wdim = ARGV[8];
else
	wdim = "1100x1100";
end

if( ARGV[9] != nil)
	region = ARGV[9];
else
	region = "";
end

if( ARGV[10] != nil)
        addregion = ARGV[10];
else
        addregion = "";
end

if( ARGV[11] != nil)
        addregion2 = ARGV[11];
else
        addregion2 = "";
end

if( ARGV[12] != nil)
        addcom = ARGV[12];
else
        addcom = "";
end

puts "region 1 " + region.to_s
puts "region 2 " + addregion.to_s
puts "region 3 " + addregion2.to_s

command = "";

if scale2.to_f != -1
	command += " -scale limits ";
	command += scale1.to_s;
	command += " ";
	command += scale2.to_s
else
	if scale1.to_f == 1
                command += " -scale linear "
        end
	if scale1.to_f == 2
		command += " -scale squared "
	end
end



#900x900
command += " -geometry " + wdim.to_s + " -zoom " + zzz.to_s + "  -smooth radius " + smooth_radius.to_s + " -smooth yes -cmap " + cmap.to_s
#command += " -catalog load ~/Projects_script/catalogs/3eg.cat "
command += " -wcs sky galactic -wcs skyformat degrees -grid yes -single "

# command += " -pan to 98 131 physical"

if region != ""
	if region.include?(".reg") == true
        command += " -region " + region.to_s + " ";
    end
    if region.include?(".con") == true
        command += " -contour load " + region.to_s + " wcs galactic red 4 yes"
    end
end

if addregion.to_s != ""
	if addregion.include?(".reg") == true
        command += " -region " + addregion.to_s + " ";
    end
    if addregion.include?(".con") == true
        command += " -contour load " + addregion.to_s + " wcs galactic yellow 2 yes"
    end
end

if addregion2.to_s != ""
	if addregion2.include?(".reg") == true
        command += " -region " + addregion2.to_s + " ";
    end
    if addregion2.include?(".con") == true
        command += " -contour load " + addregion2.to_s + " wcs galactic yellow 2 yes"
    end
end


if typeout.to_s == "png"
	command += " -saveimage png ";
end

if typeout.to_s == "jpg"
	command += " -saveimage jpeg ";
end

cmd = ""

output = %x[xvfb-run -h]
lines = output.split("\n")

xvfb_option = "-a"
lines.each do |line|
	#puts "line-"+line
	if(line.start_with?('-d'))
		#puts "preso -d"
		xvfb_option = "-d"
	end
end


cmd = "xvfb-run "+xvfb_option+" -s \"-screen 0 2000x2000x24\""
#cmd = "DISPLAY=:3 "


cmd += " ds9 -fits ";
cmd += names;
cmd += command;
cmd += nameo
if typeout.to_s == "png"
	cmd += ".png";
end
if typeout.to_s == "jpg"
	cmd += ".jpg";
end

if addcom.to_s != ""
	cmd += " "
	cmd += addcom;
	cmd += " "
end

cmd += " -exit";

puts cmd;
system(cmd);
#File.new(abspath=PATH_RES + "/commands/ds9_" + (rand * 100000000).to_i.to_s + ".genimg", "w")
#File.write(Dir.pwd)
#File.write(cmd)
#File.close()
