#! /usr/bin/ruby
#0) pattern (OB*.res)
#1) output file name
#2) source name
#3) min sqrt(TS) (optional, default 0)
#quelli di seguito sono parametri opzionali, quindi da impostare tipo parametro=valore
#6) radius: max distance from starting source (using l_peak, b_peak) (optional, default -999)
#7) maxR: remove source if outside contour level (optional, default -999 none) - only if calculated - if >= 0 enable this option and use this number as systematic error
#9) t0: for calc phase (in MJD)
#10) period: for calc phase
#11) minexp: remove source if exposure < minexp, default -1

#generate a LC list reading the files with the pattern contained in prefix

load ENV["AGILE"] + "scripts/conf.rb"
load ENV["AGILE"] + "scripts/MultiOutput.rb"

datautils = DataUtils.new
agilefov = AgileFOV.new

if ARGV[0].to_s == "help" || ARGV[0].to_s == "h" || ARGV[0] == nil
	system("head -14 " + $0 );
	exit;
end

output = ARGV[1];
sourcename = ARGV[2]

a = Dir[ARGV[0] + "_" + sourcename + ".src"]
a.sort!
ndim = a.size();

minsqrtTS = ARGV[3]


radius = -999;
maxR = -999;
t0 = 0;
period = 0;
minexp = -1;

clean=0
for i in 4...ARGV.size
	if ARGV[i] == nil
		break;
	else
		keyw = ARGV[i].split("=")[0];
		value = ARGV[i].split("=")[1];
		puts keyw.to_s + " " + value.to_s
		case keyw
			when "radius"
				radius = value;
			when "maxR"
				maxR = value;
			when "t0"
				t0 = value;
			when "period"
				period = value;
			when "minexp"
				minexp = value;
		end
	end

end

out1file = output.to_s + ".lc"
outlcfile = File.new(out1file, "w");
out6file = output.to_s + ".resfinalfull"
outrfffile = File.new(out6file, "w");
out3file = output.to_s + ".ob"
outobfile = File.new(out3file, "w");
out5file = output.to_s + ".tlist"
outtlistfile = File.new(out5file, "w");
outafile = output.to_s + ".res"
outresfile = File.new(outafile, "w");

startmjd = -1;
endmjd = -1;
index = 0;
a.each do | xx |
	puts xx
			
	obname = xx.split("_")[0];

	mo = MultiOutput.new;
	mo.setCalcPhase(t0, period);
	mo.readDataSingleSource(xx);
	
	#identificazione dell'off-axis
	lobs = agilefov.longitudeFromPeriod2(mo.timestart_tt, mo.timestop_tt);
	bobs = agilefov.latitudeFromPeriod2(mo.timestart_tt, mo.timestop_tt);
	
	#distance from center of fov
	distfov = -1
	if(lobs != -999 and bobs != -999)
		distfov = datautils.distance(lobs, bobs, mo.l_peak, mo.b_peak)
	end
	
	#calcolo fase
	phase = mo.orbitalphase;
	
	#dati ellisse
	ellipseres = mo.fullellipseline
	
	#remove source if exposure < minexp
	if minexp.to_f >= 0
		puts mo.exposure.to_f 
		if mo.exposure.to_f < minexp.to_f
			next
		end
	end
	
	
	#remove source if outside contour level (optional, default -999 none) - if >= 0 enable this option and use this number as systematic error
	criteriaellipse = false
	if maxR.to_f >= 0
		if mo.r.to_f > 0
			if (mo.r.to_f + maxR.to_f) > mo.distellipse.to_f
				criteriaellipse = true
			end
		end
	end
	
	#max distance from starting source (using l_peak, b_peak) 
	criteriadist = false
	if radius.to_f >= 0
		if mo.dist.to_f > radius.to_f
				criteriadist = true
		end
	end
	
	#se la sorgente è troppo distante ma sta dentro l'ellisse -> do nothing
	if criteriadist == true and (maxR.to_f > 0 and criteriaellipse == false)
		
	else
		if criteriaellipse == true or criteriadist == true	
			next
		end
	end
	
	
	if not (mo.sqrtTS.to_f >= minsqrtTS.to_f && mo.flux_ul.to_f != 0.0)
		next
	end 
	
	puts mo.sqrtTS;
	
	#scrittura resfinalfull
	outrfffile.write(mo.multiOutputLineFull4(format("%05d ", index),  obname, format("%.2f ", distfov)) + "\n")

	#scrittura res
	outresfile.write(mo.multiOutputLineShort4(format("%05d ", index),  obname, format("%.2f ", distfov)) + "\n")

	#.lc
	outlc = ""
	mints1 = 0.0
	mints2 = 2.0
	outlc = mo.flux_ul.to_s + "\t0\t\t1\t"
	if mo.sqrtTS.to_f < mints1.to_f || mo.sqrtTS.to_s == "nan"
		outlc = "0\t0\t1\t"
	end
	if mo.sqrtTS.to_f > mints1.to_f and mo.sqrtTS.to_f < mints2.to_f and 
		outlc = mo.flux_ul.to_s + "\t0\t\t1\t"
	end
	if mo.sqrtTS.to_f >= mints2.to_f
		outlc = mo.flux.to_s + "\t" + mo.flux_error.to_s + "\t0\t";
	end

	outlc = outlc + format("%.2f", mo.timestart_mjd) + "\t" + format("%.2f", mo.timestop_mjd.to_f-mo.timestart_mjd.to_f) + "\t" + format("%.2f",distfov) + "\t" + obname.to_s + "\t" + format("%.2f",mo.sqrtTS) + "\t" + mo.exposure.to_s + "\t" +  mo.galcoeff.to_s + "\t" + mo.isocoeff.to_s + "\t" + format("%3.2f", mo.dist.to_f) + "\t" + "( " + ellipseres.chomp.to_s + " ) \t" + mo.sicalc.to_s + "\t" + mo.sicalc_error.to_s + "\n";
 
	outlcfile.write(outlc);
	
	#scirttura .tlist
	outtlistfile.write(format("%.6f", mo.timestart_tt.to_f) + "\t" + format("%.6f", mo.timestop_tt.to_f)+ "\n");
	
	#scrittura .ob
	outlcob = "";
	outlcob = format("%.6f", mo.timestart_tt.to_f) + "\t" + format("%.6f", mo.timestop_tt.to_f) + "\t" + obname.to_s + "\t" + mo.galcoeff.to_s + "\t" + mo.isocoeff.to_s + "\t" + lobs.to_s + "\t" + bobs.to_s + "\t" + mo.exposure.to_s + "\n";
	outobfile.write(outlcob);
	
	index = index + 1;
	
end
outobfile.close();
outlcfile.close();
outtlistfile.close();
outrfffile.close();

