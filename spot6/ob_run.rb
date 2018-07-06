#! /usr/bin/ruby
#1) .conffile
#2) .ob
# where .ob contains
#tstart tstop obname galcoeff (-1) isocoeff (-1) lpoiting (-1) bpointing (-1)

load ENV["AGILE"] + "/scripts/sor/sorpaths.rb"
load ENV["AGILE"] + "/scripts/conf.rb"

if ARGV[0].to_s == "help" || ARGV[0].to_s == "h" || ARGV[0] == nil
	system("head -6 " + $0 );
	exit;
end

conffile = ARGV[0]
obfile = ARGV[1]


abspath=PATH_RES

File.open(obfile).each_line do | line |
	if line[0] == 35 
		next
	end
	ob = line.split(" ")
	tstart = ob[0]
	tstop = ob[1]
	obname = ob[2]
	galcoeff = ob[3]
	isocoeff = ob[4]
	lpointing = ob[5]
	bpointing = ob[6]
	
	newconffilename = "/tmp/ob" + (rand * 100000000).to_i.to_s + ".conf"
	
	fo = File.new(newconffilename, "w")
	
	index = 0
	File.open(conffile).each_line do | line |
			out = line
			
			if index.to_i == 2
				out = tstart.to_i.to_s + "\n"
			end
			if index.to_i == 3
				out = tstop.to_i.to_s + "\n"
			end
			if index.to_i == 8
				out = galcoeff  + "\n"
			end
			if index.to_i == 9
				out = isocoeff  + "\n"
			end
			if index.to_i == 10
				mapparam = line.chomp
				mapparam = mapparam + " lpointing="+lpointing + " bpointing=" + bpointing  + "\n"
				out = mapparam
			end
			

			if index.to_i == 25
				out = line.chomp + "_" + tstart.to_i.to_s + "_" +  tstop.to_i.to_s + "_" + obname + "\n"
				#abspath += out.chomp + "/"
			end
			fo.write(out);
			index = index + 1
		
		
	end
	fo.close()
	
	system("mv " + newconffilename + " " + abspath + "/commands/");
	
end