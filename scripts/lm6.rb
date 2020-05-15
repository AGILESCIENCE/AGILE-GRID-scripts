#! /usr/bin/ruby
#script for BUILD25
#00) filter DIR (es: FM3.119_2_I0023, FM3.119_ASDCe_I0023, FM3.119_ASDCSTDf_I0023, FM3.119_ASDCSTDk_I0023)
#01) output file name prefix
#02) l GRB
#03) b GRB
#04) T0

#optional

#05) t1s
#06) t2s
#07) t1b
#08) t1bshift
#09) t2b
#10) t2bshift
#11) fovradmax: fov rad max, to be used also with fovbinnumber, default 60
#12) albedorad: default 80
#13) radius: of search, default 10
#14) timelist: a file with a list of tstart/stop
#15) timetype: CONTACT, MJD, UTC, TT, default TT

#scan

#16) timeslot
#17) timeslotstart
#18) timeslotstop

load ENV["AGILE"] + "/scripts/conf.rb"

datautils = DataUtils.new
parameters = Parameters.new

if ARGV[0].to_s == "help" || ARGV[0] == nil || ARGV[0] == "h"
	system("head -28 " + $0 );
	exit;
end

filter = ARGV[0];
output = ARGV[1];
l = ARGV[2];
b = ARGV[3];
t0 = ARGV[4];

datautils.extractFilterDir(filter)
filterdir = datautils.filterdir

filterbase2 = filter.split("_")[0] + "_" + filter.split("_")[1];

parameters.processInput(5, ARGV, filter)

emin1 = parameters.emin;
emax1 = parameters.emax;
indexlog = datautils.logindex(filterbase2)
indexfilter = datautils.evtindex(filterbase2)
puts "indexfilter: " + indexfilter.to_s


index_name_cor = BASEDIR_ARCHIVE.to_s + "/DATA/INDEX/3901.cor.index"
puts "index name cor: " + index_name_cor;
tstart = 0

if(parameters.timetype == "CONTACT")
	#estrazione dei tempi min e max dal corfileindex

	datautils.extractTimeMinMaxForContact(index_name_cor, t0);
	tstart = datautils.tmin;

end
if(parameters.timetype == "TT")
	tstart = t0;

end
if(parameters.timetype == "MJD")
	tstart = datautils.time_mjd_to_tt(t0);
end

if(parameters.timetype == "UTC")
	tstart = datautils.time_utc_to_tt(t0)
end

puts "T0 (TT): " + tstart.to_s;


if tstart.to_f == 0
	puts "Error in TMIN, exit"
	exit(1)
end


#phasecode

parameters.setPhaseCode(tstart)

lonpole = 180.0

#selezione della sar matrix
filterbase = filter.split("_")[0]
datautils.getResponseMatrix(filter);
sarmatrix = datautils.sarmatrix;
edpmatrix = datautils.edpmatrix;

	# execute lm
	listfile=output + ".lm"
	sarmatrixfull = PATHMODEL + sarmatrix
	edpmatrixfull = " None "
	if parameters.useEDPmatrixforEXP.to_i == 1
		edpmatrixfull =  PATHMODEL + edpmatrix
	end

	# create 3  map for AG_lm6
	puts "emax"+parameters.emax.to_s
	#S
	tmin = tstart.to_f - parameters.t1s.to_f
	tmax = tstart.to_f + parameters.t2s.to_f
	cmd = "map.rb "+filter+" " + output + "_S " + tmin.to_s + " " + tmax.to_s + " "+l.to_s+" "+b.to_s+" emin="+parameters.emin.to_s+" emax="+parameters.emax.to_s+" filtercode="+parameters.filtercode.to_s+" timestep=1 mapsize=25 binsize=0.5 fovradmax="+parameters.fovradmax.to_s+" albedorad="+parameters.albedorad.to_s
	datautils.execute("", cmd);

	#B1
	tmin = tstart.to_f - parameters.t1s.to_f - parameters.shiftt1b.to_f - parameters.t1b.to_f
	tmax = tstart.to_f - parameters.t1s.to_f - parameters.shiftt1b.to_f
	cmd = "map.rb "+filter+" " + output + "_B1 " + tmin.to_s + " " + tmax.to_s + " "+l.to_s+" "+b.to_s+" emin="+parameters.emin.to_s+" emax="+parameters.emax.to_s+" filtercode="+parameters.filtercode.to_s+" timestep=1 mapsize=25 binsize=0.5 fovradmax="+parameters.fovradmax.to_s+" albedorad="+parameters.albedorad.to_s
	datautils.execute("", cmd);

	#B2
	if t2b != 0
		tmin = tstart.to_f + parameters.t2s.to_f + parameters.shiftt2b.to_f
		tmax = tstart.to_f + parameters.t2s.to_f + parameters.shiftt2b.to_f + parameters.t2b.to_f
		cmd = "map.rb "+filter+" " + output + "_B2 " + tmin.to_s + " " + tmax.to_s + " "+l.to_s+" "+b.to_s+" emin="+parameters.emin.to_s+" emax="+parameters.emax.to_s+" filtercode="+parameters.filtercode.to_s+" timestep=1 mapsize=25 binsize=0.5 fovradmax="+parameters.fovradmax.to_s+" albedorad="+parameters.albedorad.to_s
		datautils.execute("", cmd);
		t2bmapc = output+"_B2.cts.gz"
		t2bmape = output+"_B2.exp.gz"
	else
		t2bmapc = "None"
		t2bmape = "None"
	end

	cmd = "cp " + PATH + "share/AG_lm6.par . "
	datautils.execute("", cmd);
	cmd = "export PFILES=.:$PFILES; "+PATH+"bin/AG_lm6 "+output+".lm "+output+"_S.cts.gz "+output+"_S.exp.gz "+output+"_B1.cts.gz "+output+"_B1.exp.gz "+ t2bmapc + " " + t2bmape + " no " +l.to_s+" "+b.to_s+" 10 no no no 0 15 10"
	puts cmd
	datautils.execute("", cmd);
