#0) num of bin
#1) tbinsize (1 day = 86400)
#2) tstart
#3) tstop
#4) output file name
#5) optional: ob file input
#6) optional: list of time to exclude

load ENV["AGILE"] + "/scripts/conf.rb"

datautils = DataUtils.new
agilefov = AgileFOV.new

if ARGV[0].to_s == "help" || ARGV[0].to_s == "h" || ARGV[0] == nil
	system("head -7 " + $0 );
	exit;
end
numbin = ARGV[0]

tbinsize = ARGV[1]
tstep =  tbinsize.to_f * (numbin.to_f )

tstartall = ARGV[2]
tend = ARGV[3]
out = ARGV[4]


if ARGV[5] != nil
	galisofile = ARGV[5];
else
	galisofile = ""
end

if ARGV[6] != nil
	timeexclude = ARGV[6].to_s;
else
	timeexclude = "";
end


for i in 1..numbin.to_i

	fout = File.open(out.to_s + i.to_s, "w")
	
	tstart = tstartall.to_f + (tbinsize.to_f * (i.to_i - 1) ).to_f

	time = tstart.to_f - tstep.to_f
	
	t2 = 0
	t1 = 0
	
	while t1.to_f <= tend.to_f
		
		time = time.to_f + tstep.to_f
		t1 = time
		t2 = time.to_f + tbinsize.to_f

		
		
		if t1.to_f < tend.to_f
			outok = false
			if galisofile.to_s != ""
				File.open(galisofile).each_line do | line |
					ll = line.split(" ");
					tA = ll[0]
					tB = ll[1]
					if t1.to_f >= tA.to_f and t1.to_f <= tB.to_f
						outok = true
					end
				end
			else
				outok = true
			end
			
			if timeexclude.to_s != ""
				File.open(timeexclude).each_line do | line |
					ll = line.split(" ");
					tA = ll[0]
					tB = ll[1]
					
					if t1.to_f >= tA.to_f and t1.to_f <= tB.to_f
						t1 = tB
					end
					if t2.to_f >= tA.to_f and t2.to_f <= tB.to_f
						t2 = tA 
					end
					if tB <= tA
						outok = false
					end
				end
				
			end
		
			if outok  then
				cm =  t1.to_s + "\t" + t2.to_s + "\n"
				fout.write(cm)
				puts cm
			end
		end
		
	end
	
	fout.close()
	
end