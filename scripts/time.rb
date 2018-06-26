#! /usr/bin/ruby
#0) input type (0 = TT, 1 = UTC, 2=MJD, 3= CONTACT)
#1) time
#2) time2 (optional)
#3) title

load ENV["AGILE"] + "/scripts/conf.rb"

if ARGV[0].to_s == "help" || ARGV[0] == nil || ARGV[0] == "h"
	system("head -6 " + $0 );
	exit;
end

type = ARGV[0]
time = ARGV[1]
if ARGV[2] != nil
	time2 = ARGV[2].to_s;
else
	time2 = 0.0;
end
if ARGV[3] != nil
	title  = ARGV[3].to_s;
else
	title  = "";
end

datautils = DataUtils.new


if type.to_i == 0
	tt1_A = time
	utc1_A = datautils.time_tt_to_utc(time)
	mjd1_A = datautils.time_tt_to_mjd(time)
end
if type.to_i == 1
	tt1_A = datautils.time_utc_to_tt(time)
	utc1_A = time
	mjd1_A = datautils.time_utc_to_mjd(time)
end
if type.to_i == 2
	tt1_A = datautils.time_mjd_to_tt(time)
	utc1_A = datautils.time_mjd_to_utc(time)
	mjd1_A = time
end
if type.to_i == 3
	index_name_cor = "/AGILE_PROC3/DATA/INDEX/3901.cor.index"
	datautils.extractTimeMinMaxForContact(index_name_cor, time);
    time1_A = datautils.tmin;
    time1_B = datautils.tmax;
	tt1_A = time1_A
	tt1_B = time1_B
	utc1_A = datautils.time_tt_to_utc(time1_A)
    mjd1_A = datautils.time_utc_to_mjd(utc1_A)
	utc1_B = datautils.time_tt_to_utc(time1_B)
    mjd1_B = datautils.time_utc_to_mjd(utc1_B)
end

if type.to_i < 3
	datautils.extractContactLocalIndex(ENV['AGILE'] + "/scripts/EVT.index", tt1_A) 
	orbit = datautils.ec.to_i
else
	orbit = time
end

if time2 != 0

	datautils2 = DataUtils.new
	if type.to_i == 0
		tt2_A = time2
		utc2_A = datautils2.time_tt_to_utc(tt2_A)
		mjd2_A = datautils2.time_utc_to_mjd(utc2_A)
	end
	if type.to_i == 1
		tt2_A = datautils2.time_utc_to_tt(time2)
		utc2_A = time2
		mjd2_A = datautils2.time_utc_to_mjd(time2)
	end
	if type.to_i == 2
		tt2_A = datautils2.time_mjd_to_tt(time2)
		utc2_A = datautils2.time_mjd_to_utc(time2)
		mjd2_A = time2
	end
	if type.to_i == 3
		
		index_name_cor = "/AGILE_PROC3/DATA/INDEX/3901.cor.index"
		datautils2.extractTimeMinMaxForContact(index_name_cor, time2);
        time2_A = datautils2.tmin;
        time2_B = datautils2.tmax;
		tt2_A = time2_A
		tt2_B = time2_B
		utc2_B = datautils2.time_tt_to_utc(tt2_B)
        mjd2_B= datautils2.time_utc_to_mjd(utc2_B)
	end
		
	if type.to_i < 3
		datautils2.extractContactLocalIndex("EVT.index", tt2_A) 
		orbit2 = datautils2.ec.to_i
	else
		orbit2 = time2
	end
	
	puts tt1_A.to_s + " " + tt2_B.to_s + " " + utc1_A.to_s + " " + utc2_B.to_s + " " + format("%.2f", mjd1_A).to_s + " " + format("%.2f", mjd2_B).to_s + " " + (orbit.to_i).to_s + " " + (orbit2.to_i).to_s
	
else

	puts "TT " + tt1_A.to_s
	puts "UTC " + utc1_A.to_s
	puts "MJD " + format("%.7f", mjd1_A).to_s

	puts tt1_A.to_s + " " + utc1_A.to_s + " " + format("%.7f", mjd1_A).to_s + " " + orbit.to_s

end


