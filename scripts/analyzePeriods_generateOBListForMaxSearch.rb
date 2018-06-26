#! /usr/bin/ruby
#0) tcenter
#1) tbinsize (1 day = 86400)
#2) N
#3) output file name

#optional
#1) tmin
#2) tmax
#3) obfile: .ob file with gal e iso

#due modi di ricerca
#1) da tcenter	-> [tcenter, tcenter + tbinsize * k]
#				-> [tcenter - tbinsize * k, tcenter]
#2) da tcenter  -> [tcenter - tbinsize * k, tcenter + tbinsize * k]
# k = 1..N e fino a quando non supero tmin e tmax come seconda condizione)

load ENV["AGILE"] + "/scripts/conf.rb"

datautils = DataUtils.new
agilefov = AgileFOV.new

if ARGV[0].to_s == "help" || ARGV[0].to_s == "h" || ARGV[0] == nil
	system("head -16 " + $0 );
	exit;
end

tcenter = ARGV[0]
tbinsize = ARGV[1]
n = ARGV[2]
out = ARGV[3]

tmin = -999
tmax = -999
obfile = nil

for i in 4...ARGV.size
	if ARGV[i] == nil
		break;
	else
		keyw = ARGV[i].split("=")[0];
		value = ARGV[i].split("=")[1];
		puts keyw.to_s + " " + value.to_s
		case keyw
			when "tmin"
				tmin = value;
			when "tmax"
				tmax = value;
			when "obfile"
				obfile = value;
		end
	end
end



index = 0

fout = File.open(out, "w")

#1
for i in 1..n.to_i
	

	time = tcenter
	t1 = time
	t2 = time.to_f + tbinsize.to_f * i.to_f
	
	index = index + 1

	lobs = agilefov.longitudeFromPeriod2(t1, t2);
 	bobs = agilefov.latitudeFromPeriod2(t1, t2);

	gal = -1
	iso = -1
	if obfile != nil
		File.open(obfile).each_line do | line |
			t = line.split(" ")
			if t1.to_f >= t[0].to_f and t1.to_f <= t[1].to_f
				gal = t[3]
				iso = t[4]
			end
		end
	end
	
	cm =  t1.to_s + "\t" + t2.to_s + "\tOB" + format("%05d", index) + "\t" + gal.to_s + "\t" + iso.to_s + "\t" + lobs.to_s + "\t" + bobs.to_s + "\n"
	
	fout.write(cm)
	
end

#2
for i in 1..n.to_i
	

	time = tcenter
	t1 = time.to_f - tbinsize.to_f * i.to_f
	t2 = time
	
	
	index = index + 1

	lobs = agilefov.longitudeFromPeriod2(t1, t2);
 	bobs = agilefov.latitudeFromPeriod2(t1, t2);

	gal = -1
	iso = -1
	if obfile != nil
		File.open(obfile).each_line do | line |
			t = line.split(" ")
			if t1.to_f >= t[0].to_f and t1.to_f <= t[1].to_f
				gal = t[3]
				iso = t[4]
			end
		end
	end
	
	cm =  t1.to_s + "\t" + t2.to_s + "\tOB" + format("%04d", index) + "\t" + gal.to_s + "\t" + iso.to_s + "\t" + lobs.to_s + "\t" + bobs.to_s + "\n"
	
	fout.write(cm)
	
end

#3
for i in 1..n.to_i
	

	time = tcenter
	t1 = time.to_f - tbinsize.to_f * i.to_f
	t2 = time.to_f + tbinsize.to_f * i.to_f
	
	index = index + 1

	lobs = agilefov.longitudeFromPeriod2(t1, t2);
 	bobs = agilefov.latitudeFromPeriod2(t1, t2);

	gal = -1
	iso = -1
	if obfile != nil
		File.open(obfile).each_line do | line |
			t = line.split(" ")
			if t1.to_f >= t[0].to_f and t1.to_f <= t[1].to_f
				gal = t[3]
				iso = t[4]
			end
		end
	end
	
	cm =  t1.to_s + "\t" + t2.to_s + "\tOB" + format("%04d", index) + "\t" + gal.to_s + "\t" + iso.to_s + "\t" + lobs.to_s + "\t" + bobs.to_s + "\n"
	
	fout.write(cm)
	
end

fout.close()