#! /usr/bin/ruby
#0) tstep 
#1) tbinsize (1 day = 86400)
#2) .ob file with gal e iso (generated with readOB.rb)
#3) output file name
#4) prefix (OP, OS, OB)
#5) endattend: end the period of analysis at tend of the ob file (1), otherwise close the period of analysis using tbinsize 

load ENV["AGILE"] + "/scripts/conf.rb"

datautils = DataUtils.new
agilefov = AgileFOV.new

if ARGV[0].to_s == "help" || ARGV[0].to_s == "h" || ARGV[0] == nil
	system("head -7 " + $0 );
	exit;
end

tstep = ARGV[0]
tbinsize = ARGV[1]
galisofile = ARGV[2]
out = ARGV[3]
prefix=ARGV[4]
endattend=ARGV[5]
 
fout = File.open(out, "w")
index = 0
File.open(galisofile).each_line do | line |
	ll = line.split(" ");
	tstart = ll[0]
	tend = ll[1]
	lo = ll[5]
	bo = ll[6]

	time = tstart.to_f - tstep.to_f
	
	t2 = 0
	t1 = 0
	
	while t1.to_f <= tend.to_f
		
	
		time = time.to_f + tstep.to_f
		t1 = time
		t2 = time.to_f + tbinsize.to_f
		index = index + 1
	
		lobs = agilefov.longitudeFromPeriod2(t1, t2);
		if lobs == -999
			lobs = lo
		end
		bobs = agilefov.latitudeFromPeriod2(t1, t2);
		if bobs == -999
			bobs = bo
		end
		gal = -1;
		iso = -1;
		File.open(galisofile).each_line do | line |
			t = line.split(" ")
			if t1 >= t[0].to_f and t1 <= t[1].to_f
				gal = t[3]
				iso = t[4]
			end
			
		end
		if t2.to_f > tend.to_f
			if endattend.to_i == 1
				t2 = tend
			end
		end

		if gal.to_f == 0
			gal=-1
		end
		if iso.to_f == 0
			iso=-1
		end
		cm =  t1.to_s + "\t" + t2.to_s + "\t" + prefix.to_s + format("%05d", index) + "\t" + gal.to_s + "\t" + iso.to_s + "\t" + lobs.to_s + "\t" + bobs.to_s + "\n"
		#cm = cm + "\n"		
		if t1.to_f < tend.to_f
			fout.write(cm)
			puts cm
		end
		
	end
	
	

	
end
fout.close()
