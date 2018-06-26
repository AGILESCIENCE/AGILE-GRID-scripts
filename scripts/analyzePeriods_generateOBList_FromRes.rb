#! /usr/bin/ruby
#0) input res file name
#1) output file name

load ENV["AGILE"] + "/scripts/conf.rb"

if ARGV[0].to_s == "help" || ARGV[0].to_s == "h" || ARGV[0] == nil
	system("head -7 " + $0 );
	exit;
end

res = ARGV[0]
out = ARGV[1]

fout = File.open(out, "w")
index = 0
File.open(res).each_line do | line |
	ll = line.split(" ");
	tstart = ll[15]
	tend = ll[16]
	lpointing = ll[17]
	bpointing = ll[18]
	gal = ll[8]
	iso = ll[9]
	name = ll[11]
		
	cm =  tstart.to_s + "\t" + tend.to_s + "\t" + name.to_s + "\t" + gal.to_s + "\t" +  iso.to_s + "\t" + lpointing.to_s + "\t" + bpointing.to_s + "\n"	
	fout.write(cm)
	puts cm
end
fout.close()