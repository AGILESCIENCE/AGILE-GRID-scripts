#! /usr/bin/ruby
#0) l
#1) b
#2) fov (use and additional 4-5 degree)
#3) add or remove time (+ - secs)
#4) gal (use -999 to keep it free)
#5) iso (use -999 to keep it free)


load ENV["AGILE"] + "/scripts/conf.rb"

if ARGV[0].to_s == "help" || ARGV[0] == nil
	system("head -6 " + $0 );
	exit;
end

l = ARGV[0]
b = ARGV[1]
fov = ARGV[2]
gal = -999
iso = -999

if ARGV[3] != nil
	addtime = ARGV[3].to_f;
else
	addtime = 0;
end

if ARGV[4] != nil
        gal = ARGV[4].to_f;
else
        gal = -999;
end
if ARGV[5] != nil
        iso = ARGV[5].to_f;
else
        iso = -999;
end
datautils = DataUtils.new
agilefov = AgileFOV.new

File.open(ENV["AGILE"]+"/scripts/AGILEOBPOINT.list").each_line do | line |
	af = line.split(" ");
	lc = af[1];
	bc = af[2];
	obname = af[0];
	tstart = af[3];
	tstop = af[4];
	dist = datautils.distance(lc, bc, l, b)
	if dist.to_f <= fov.to_f
		t0 = datautils.time_utc_to_tt(tstart)
		t1 = datautils.time_utc_to_tt(tstop)
		#puts line.chomp.to_s + " " + t0.to_s + " " + t1.to_s	
		puts (t0.to_f + addtime.to_f).to_s + "\t" + (t1.to_f).to_s + "\t" + obname.to_s + "\t"+gal.to_s+"\t"+iso.to_s+"\t" + lc.to_s + "\t" + bc.to_s
		
# 		t0mjd = datautils.time_utc_to_mjd(tstart)
# 		t1mjd = datautils.time_utc_to_mjd(tstop)
# 		puts obname.to_s + " & " + lc.to_s + " & " + bc.to_s + " & " + t0mjd.to_s + "-" + t1mjd.to_s + " & " + format("%.2f", dist) + "\\\\"
	end
end
