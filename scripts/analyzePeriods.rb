#! /usr/bin/ruby
# Prende una lista di periodi da analizzare (param 1) e li analizza in base ad una lista passata come param 2
#00) filter
#01) list of | tstart tstop name gal iso l_pointing b_pointing 
# 202152000.000040 203016000.000040 OB7500	-1	1 -1 -1
#02) type of analysis: generate only maps (1) / multi5 (2) / map + multi5 (3)  
#03) l of the map center
#04) b of the map center
#05) list of sources (.multi) for multi5
#06) cluster, 1 activate cluster, 0 do not activate
#07) fixisogal, 0 do not fix, 1 fix (at the value of the .ob)
#08) binsize
#09) mapsize
#10) optional parameters for map.rb " "
#11) optional parameters for multi5.rb

load  ENV["AGILE"] + "/scripts/conf.rb"


if ARGV[0].to_s == "help" || ARGV[0].to_s == "h" || ARGV[0] == nil
	system("head -16 " + $0 );
	exit;
end
alikeutils = AlikeUtils.new;
datautils = DataUtils.new;
parameters = Parameters.new;

system("hostname")

filter = ARGV[0];
inputfileob = ARGV[1];
analysistype = ARGV[2]
l = ARGV[3];
b = ARGV[4];
sourcelist = ARGV[5]
cluster = ARGV[6]
fixisogal = ARGV[7]
binsize= ARGV[8]
mapsize = ARGV[9]
optmap = ARGV[10]
optmulti = ARGV[11]

index = 0;
File.open(inputfileob).each_line do |line|
	index = index + 1;
end

puts index
startt = Array.new(index)
endt = Array.new(index)
vp = Array.new(index)
gal = Array.new(index)
iso = Array.new(index)
ll1 = Array.new(index) # l pointing
bb1 = Array.new(index) # b pointing

index = 0
File.open(inputfileob).each_line do |line|
	puts line
	ll = line.split(" ");
	startt[index] = ll[0];
	endt[index] = ll[1];
	vp[index] = ll[2];
	gal[index] = ll[3];
	iso[index] = ll[4];
	ll1[index] = ll[5];
	bb1[index] = ll[6];
	index = index + 1;
end
puts "size " + startt.size().to_s


for i in 0..startt.size()-1
	puts "index " + i.to_s + " " + endt[i].to_s
	
	if binsize.to_f >= 0.1
		prefix = vp[i].to_s +  "_" + filter.to_s + "_b0" + format("%0d", binsize.to_f * 100);
	else
		prefix = vp[i].to_s +  "_" + filter.to_s + "_b00" + format("%0d", binsize.to_f * 100);
	end
	puts prefix
	
	datautils.executeprefixcluster(prefix, cluster)

	puts "FIXISOGAL " + fixisogal.to_s
	if fixisogal.to_i == 1
		galtmp = gal[i].to_s;
		isotmp = iso[i].to_s;
	else
		galtmp = "-1";
		isotmp = "-1";
	end
	puts gal[i].to_s + " " + iso[i].to_s
	puts galtmp.to_s + " " + isotmp.to_s
	
	name = vp[i].to_s + "_" + filter.to_s + "_b0" + format("%0d", binsize.to_f * 100);
	cts = name.to_s + ".cts.gz"


	if analysistype.to_i >= 1 
		outputres = name + ".res"
		
		if cluster.to_i >= 1
			system("cat run.ll > " + prefix.to_s + ".run")
		end

		if analysistype.to_i ==1 || analysistype.to_i == 3
			cmd = "ruby " + ENV["AGILE"] + "/scripts/map.rb " + filter.to_s + " " + name + " " + startt[i].to_s + " " + endt[i].to_s + " " + l.to_s + " " + b.to_s + " timetype=TT " 
			if ll1[i].to_s != "-1" and bb1[i].to_s != "-1"
				cmd += " lpointing=" + ll1[i].to_s + " bpointing=" + bb1[i].to_s + " "
			end
			cmd += " binsize=" + binsize.to_s
			cmd += " mapsize=" + mapsize.to_s + " "
			cmd += optmap.to_s
			datautils.executecluster(prefix, cmd.to_s, cluster)
		end
		if analysistype.to_i ==2 || analysistype.to_i == 3
			cmd = "ruby " + ENV["AGILE"] + "/scripts/multi5.rb " + filter.to_s + " " + name + ".maplist4 " + sourcelist.to_s + " " + name.to_s + ".res "
			if fixisogal.to_i == 1
				cmd += " galcoeff=" + galtmp.to_s + " isocoeff=" + isotmp.to_s + " "
			end
			cmd += " binsize=" + binsize.to_s
			cmd += " mapsize=" + mapsize.to_s + " "
			cmd += " prefix=" + prefix.to_s + " " + optmulti.to_s
			datautils.executecluster(prefix, cmd.to_s, cluster)
		end
	end
	
	if cluster.to_i >= 1
		system("llsubmit " + prefix.to_s + ".run")
	end
end
