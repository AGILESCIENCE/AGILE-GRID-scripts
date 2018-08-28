#! /usr/bin/ruby
#0) config file name


#TODO
#/scratch/prod/agile-B23/scripts/polygonfilter/polygonfilter_2.py -l 1 -b 2 -g MLE0000display.con.galactic TSMAPMLE0001.scanlist TSMAPMLE0001_filtered.scanlist
require 'rubygems'
require 'mysql'

load ENV["AGILE"] + "/scripts/conf.rb"
load ENV["AGILE"] + "/AGILEPIPE/env.rb"
load ENV["AGILE"] + "/AGILEPIPE/spot6/Conf.rb"

def getcoeff(path, coefftype, timeend, timetype)
	if File.exists?(PATH_RES + "/" + path + "/analysisres.ob")
		finalcoeff = "-1"
		datautils = DataUtils.new
		tstop = 0
		if(timetype == "TT")
			tstop = timeend
		end
		if(timetype == "MJD")
			tstop = datautils.time_mjd_to_tt(timeend);
		end
		if(timetype == "UTC")
			tstop = datautils.time_utc_to_tt(timeend)
		end

		tdistance = -1

		File.open(PATH_RES + "/" + path + "/analysisres.ob").each_line do | line |
			lsplit = line.split(" ")
			if tstop.to_f >= lsplit[0].to_f and tstop.to_f <= lsplit[1].to_f
				tdist2 = lsplit[1].to_f - tstop.to_f
				found = false
				if tdistance == -1
					tdistance = tdist2
					found = true
				end
				if tdistance > tdist2
					tdistance = tdist2
					found = true
				end
				if found == true
					if coefftype == "gal"
						finalcoeff = lsplit[3]
					end
					if coefftype == "iso"
						finalcoeff = lsplit[4]
					end
				end
			end
		end
		if finalcoeff.split(",").size == 1
			finalcoeff = format("%.3f", finalcoeff.to_s)
		end
		return finalcoeff
	else
		return "-1"
	end
end

def extractcat(hypothesisgen, l, b, outfile)
	h = hypothesisgen.split(" ")
	cmd = "extract_catalog.rb " + PATH_RES + "/catalogs/" + h[1].to_s + " " + l.to_s + " " + b.to_s + " " + outfile + " " + h[2].to_s + " " + h[3].to_s + " " + h[4].to_s + " " + h[5].to_s + " " + h[6].to_s + " " + h[7].to_s + " " + h[8].to_s
	if h.size >= 10
		cmd = cmd + " " + h[9].to_s
	end
	if h.size >= 11
		cmd = cmd + " " + h[10].to_s
	end
	puts cmd
	system cmd

end

def spotfinder(hypothesisgen, outfile, eb, mapsize)
	h = hypothesisgen.split(" ")

	prefix = "MAP"

	if eb.to_i != 0
		prefix = "MAP"
		index = 1
		if h.size() == 8
			index = h[7].to_i
		end
		indexfile = 1
		File.open("MAP.maplist4").each_line do | line |
			if index.to_i == indexfile.to_i
				prefix = line.split(" ")[0].split(".cts.gz")[0]
			end
			indexfile = indexfile.to_i + 1
		end
	end

	map = ""
	if h[1].to_i == 0
		map = prefix + ".cts.gz"
	else
		map = prefix + ".int.gz"
	end

	cmd = "spotfinder.rb " + outfile + " " + map + " " + h[1].to_s + " " + h[2].to_s + " " + h[3].to_s + " " + h[4].to_s + " " + h[5].to_s + " " + h[6].to_s + " " + mapsize.to_s + " " + prefix + ".exp.gz"
	puts cmd
	system cmd
end

#2017 N. Parmiggiani, GPL
def cat2db(hypothesisgen, filename)

end

#2017 N. Parmiggiani, GPL
def spot6db(hypothesisgen, t_start, t_stop, filename,source_tec_name)

	h = hypothesisgen.split(" ")
	l_search = h[1]
	b_search = h[2]
	r_search = h[3]
	lc = h[4]
	fixflag = h[5]
	minexp = h[6]
	maxexpratio = h[7]
	minsqrtts = h[8]
	minflux = h[9]
	server = h[10]
	usn = h[11]
	pwd = h[12]
	dbname = h[13]

	puts "database info"+server+" "+usn+" "+pwd+" "+dbname

	con = Mysql.new server, usn, pwd, dbname

  	file = File.open(filename,"w")

	#get sqrtts of source
	#query = "select sqrtts from source_info where label = '"+source_tec_name+"'"
	#source_sqrtts=nil

	#rs = con.query(query)
	#rs.each_hash do |row|
	#	source_sqrtts = row['sqrtts']
	#end
	#puts "source_sqrtts = "+source_sqrtts


	query = "select d.sqrtts,d.counts_err,d.counts,d.id_detection,d.spectral_index,d.flux,d.l_peak,d.b_peak,d.label, d.exposure, d.exp_ratio,d.ella,d.ellb from detection d join run r on(r.id_run = d.id_run) join analysis an on (an.id_analysis = r.id_analysis)  where an.name = '"+lc+"' and t_start_tt = "+t_start+" and t_stop_tt = "+t_stop
  puts "spot6cat2-"+query
	rs = con.query(query)
	count = 0
	deltaT = t_stop.to_f - t_start.to_f

	rs.each_hash do |row|

			id_detection = row['id_detection']
			l_peak = row['l_peak']
			b_peak = row['b_peak']
			s_label = row['label']
			exposure = row['exposure']
			exp_ratio = row['exp_ratio']
			flux = row['flux']
			ella = row['ella']
			ellb = row['ellb']
			spectral_index= row['spectral_index']
			counts = row['counts']
			counts_err = row['counts_err']
			sqrtts=row['sqrtts']

			#exposure filters
			if (exposure.to_f / deltaT.to_f < minexp.to_f)
				next
			end
			if (exp_ratio.to_f > maxexpratio.to_f)
				next
			end

			#data quality filters
			if (ella.to_f > 10 || ellb.to_f >10)
				next
			end

			cts = 0;
			cts_err = 0;
			if(counts.to_f == 0 || counts_err.to_f == 0)

					cts_err = 1;
					cts =1;

			else
					cts_err = counts_err;
					cts = counts;
			end

			if(sqrtts.to_f > (5 * (cts.to_f/cts_err.to_f + Math.log(cts.to_f/cts_err.to_f))))
				next
			end

			#mints and minflux filters
			if(sqrtts.to_f <= minsqrtts.to_f)
				next
			end

			if(flux.to_f <= minflux.to_f)
				next
			end

			#distance filters
			datautils = DataUtils.new

			d = datautils.distance(l_peak.to_f, b_peak.to_f, l_search.to_f, b_search.to_f)

			if (d <= r_search.to_f) then # la prendo
					#puts "id"+id_detection
					count = count+1
					file.puts(flux+" "+l_peak+" "+b_peak+" "+spectral_index+" " + fixflag + " 2.0 "+s_label+"_"+count.to_s+" 0.0")
			end
	end

	file.close
end


datautils = DataUtils.new
alikeutils = AlikeUtils.new
parameters = Parameters.new

filenameconf = ARGV[0];



mleindex = 0;
ml = Dir["MLE????"].sort
puts "index: " + ml.size().to_s
if ml.size() > 0
	mleindex = ml[ml.size()-1].split("MLE")[1].to_i;
	mleindex = mleindex.to_i + 1
else
	cmd = "cp MLE0000.conf tmp.conf"
	puts cmd
	#system(cmd)
	cmd = "rm MAP* MLE*"
	puts cmd
	#system(cmd)
	cmd = "mv tmp.conf MLE0000.conf"
	puts cmd
	#system(cmd)
end

mle = "MLE" + format("%04d", mleindex)

filenameconfext = filenameconf
filenameconf = mle + ".conf"

if File.exists?(filenameconf) == false
	cmd = "cp " + filenameconfext  + " " + mle + ".conf "
	puts cmd
	system(cmd)
end

#estrazione lista sorgenti
fndisplayreg = mle + "display"

fnhyp0 = mle+"hypothesis0.multi"
fnhyp = mle+"hypothesis.multi"

conffile = Conf.new

conffile.process(filenameconf, fnhyp0, fndisplayreg)
analysis_name = conffile.analysis_name
filter = conffile.filter
tstart = conffile.tstart
tstop = conffile.tstop
timetype = conffile.timetype
l = conffile.l
b = conffile.b
proj = conffile.proj
galcoeff = conffile.galcoeff
isocoeff = conffile.isocoeff
mapparam = conffile.mapparam
hypothesisgen1 = conffile.hypothesisgen1
hypothesisgen2 = conffile.hypothesisgen2
radmerger = conffile.radmerger
multiparam = conffile.multiparam
tsmapparam = conffile.tsmapparam
iddisp = conffile.iddisp
dir_run_output = conffile.dir_run_output
mail = conffile.mail
run_name = conffile.run_name
ds91 = conffile.ds91
ds92 = conffile.ds92
ds93 = conffile.ds93
ds94 = conffile.ds94
regfile = conffile.regfile
detGIF = conffile.detGIF
comments = conffile.comments
reg = conffile.reg
binsize = conffile.binsize
queue = conffile.queue
dir_analysis_result = conffile.dir_analysis_result
analysis_result_minSqrtTS = conffile.analysis_result_minSqrtTS
analysis_result_sourcename = conffile.analysis_result_sourcename

l = l.to_f + 0.0001
l = l.to_s

if proj.to_s == "AIT"
	l = 0.to_s;
	b = 0.to_s;
	mapparam = mapparam.to_s + " proj=AIT mapsize=360 ";
end

mapsize = parameters.mapsize

if mapparam.split("mapsize").size() > 1
	mapsize  = mapparam.split("mapsize")[1].split("=")[1].split(" ")[0]
end

eb = parameters.energybin

if mapparam.split("eb").size() > 1
	eb  = mapparam.split("eb")[1].split("=")[1].split(" ")[0]
end

puts mleindex

if(mleindex.to_i == 0)
	cmd = "map.rb " + filter.to_s + " MAP " + tstart.to_s + " " + tstop.to_s + " " + l.to_s + " " + b.to_s + " timetype=" + timetype.to_s + " " + mapparam.to_s;
	puts cmd
	system(cmd)
end

#MAP.maplist4 hypothesis.multi MLE$10 galcoeff=$8 isocoeff=$9 $11

#TODO hypothesis1.multi
op = hypothesisgen1.split(" ")[0]

fnh1 = mle + "hypothesis1.multi"
system("touch " + fnh1)
if op != "nop"
	if op == "cat"
		extractcat(hypothesisgen1, l, b, fnh1);
	end
	if op == "cat2db"
		cat2db(hypothesisgen1, fnh1)
	end
	if op == "spotfinder"
		spotfinder(hypothesisgen1, fnh1, eb, mapsize);
	end
	if op == "spot6db"
		#spot6db l b r tstart tstop
		spot6db(hypothesisgen1, tstart, tstop, fnh1,analysis_result_sourcename)
	end

end

#TODO hypothesis2.multi
op = hypothesisgen2.split(" ")[0]

fnh2 = mle + "hypothesis2.multi"
system("touch " + fnh2)
if op != "nop"
	if op == "cat"
		extractcat(hypothesisgen2, l, b, fnh2);
	end
	if op == "cat2db"
		cat2db(hypothesisgen1, fnh2)
	end
	if op == "spotfinder"
		spotfinder(hypothesisgen2, fnh2, eb, mapsize);
	end
	if op == "spot6db"
		#spot6db l b r tstart tstop
		spot6db(hypothesisgen2, tstart, tstop, fnh2,analysis_result_sourcename)
	end
end


alikeutils.appendMulti(fnh2, fnh1, mle+"hypothesisM1.multi", radmerger );
alikeutils.appendMulti(mle+"hypothesisM1.multi", fnhyp0, mle+"hypothesisM0.multi", radmerger );

#aggiungi la selezione della parte variabile e rimuovi le sorgenti di CAT2 che hanno parte variabile


#copy input files
cmd = "mv " + mle + "hypothesisM0.multi " + fnhyp
puts cmd
system(cmd)

if detGIF != "" and detGIF != "tbd" and detGIF != "nop"
	if proj.to_s == "ARC"
		if Dir["*GIFMAP.cts.gz"].size() == 0
			deltatime = 0;
			tstartgif = 0
			if timetype.to_s == "MJD"
				deltatime = 7
				tstartgif = tstart.to_f - deltatime.to_f
			end
			if timetype.to_s == "TT"
				deltatime = 7 * 86400
				tstartgif = tstart.to_f - deltatime.to_f
			end
			if timetype.to_s == "CONTACT"
				deltatime = 7 * 14
				tstartgif = tstart.to_f - deltatime.to_f
			end
			if timetype.to_s == "UTC"
				tstarttmp = datautils.time_utc_to_tt(tstart)
				deltatime = 7 * 86400
				tstarttmp2 = tstarttmp.to_f - deltatime.to_f
				tstartgif = datautils.time_tt_to_utc(tstarttmp2)
			end
			cmd = "map.rb " + filter.to_s + " GIFMAP " + tstartgif.to_s  + " " + tstart.to_s + " " + l.to_s + " " + b.to_s + " timetype=" + timetype.to_s + " " + mapparam.to_s;
			puts cmd
			system(cmd)
		end
		cmd = "multi.rb " + filter + " GIFMAP.maplist4 " + fnhyp + " GIF" + mle
		if galcoeff != "-1"
			cmd = cmd + " galcoeff=" + galcoeff
		end
		if isocoeff != "-1"
			cmd = cmd + " isocoeff=" + isocoeff
		end
		cmd = cmd + " " + multiparam
		puts cmd
		system(cmd)
		giffot = "GIF" + mle + "_" + detGIF
		mo = MultiOutput6.new
		mo.readDataSingleSource(giffot);
		if galcoeff == "-1"
			galcoeff = mo.galcoeff
		end
		if isocoeff == "-1"
			isocoeff = mo.isocoeff
		end
	end
end

if not (multiparam.to_s == "nop" || proj.to_s == "AIT")

	if galcoeff.split("/").size > 1
		galcoefffinal = getcoeff(galcoeff, "gal", tstop, timetype)
	else
		galcoefffinal = galcoeff
	end

	if isocoeff.split("/").size > 1
		isocoefffinal = getcoeff(isocoeff, "iso", tstop, timetype)
	else
		isocoefffinal = isocoeff
	end

	cmd = "multi.rb " + filter + " MAP.maplist4 " + fnhyp + " " + mle + " galcoeff=" + galcoefffinal + " isocoeff=" + isocoefffinal + " " + multiparam
	puts cmd
	system(cmd)

	#cmd = "convertMultiResToReg.rb " + mle + " white 0.1"
	#puts cmd
	#system(cmd)
end

if not (tsmapparam.to_s == "nop")
# || proj.to_s == "AIT"
	cmd = "iterative5.rb " + filter + " MAP.maplist4 outfile=TSMAP_" + mle;
	puts cmd
	system (cmd)

	if proj.to_s == "ARC"
		conffile.plotjpgmap_arc("TSMAP_" + mle + ".lst_00.TS.fits.gz")
		conffile.plotjpgmap_arc("TSMAP_" + mle + ".lst_00.flux.fits.gz")
		conffile.plotjpgmap_arc("TSMAP_" + mle + ".lst_00.fluxul.fits.gz")
		conffile.plotjpgmap_arc("TSMAP_" + mle + ".lst_00.ISO.fits.gz")
		conffile.plotjpgmap_arc("TSMAP_" + mle + ".lst_00.GAL.fits.gz")
		conffile.plotjpgmap_arc("TSMAP_" + mle + ".lst_00.GAS.fits.gz")
	end

	if proj.to_s == "AIT"
		conffile.plotjpgmap_ait("TSMAP_" + mle + ".lst_00.TS.fits.gz")
		conffile.plotjpgmap_ait("TSMAP_" + mle + ".lst_00.flux.fits.gz")
		conffile.plotjpgmap_ait("TSMAP_" + mle + ".lst_00.fluxul.fits.gz")
		conffile.plotjpgmap_ait("TSMAP_" + mle + ".lst_00.ISO.fits.gz")
		conffile.plotjpgmap_ait("TSMAP_" + mle + ".lst_00.GAL.fits.gz")
		conffile.plotjpgmap_ait("TSMAP_" + mle + ".lst_00.GAS.fits.gz")
	end
end



#send mail
if mail.to_s != "" or mail != nil
	puts mail
	cmd = "mail -s \"end RUN\" " + mail.to_s + " < " +  mle
	puts cmd
	if mail.to_s != "nop"
        system(cmd)
    end
	#cat T1.sh | mutt -a T.res.html  bulgarelli@iasfbo.inaf.it -s 'res'
end

#TODO send notification


#generate .jpg

if proj.to_s == "ARC" and File.exists?(mle + ".reg") and File.exists?(mle + ".multi.reg")
	conffile.detsmooth()

	conffile.plotjpgcts1(mle, conffile.smooth)

	conffile.plotjpgint(mle, conffile.smooth)

	conffile.plotjpgexp(mle)

	conffile.plotjpgcts2(mle, conffile.smooth)

	conffile.plotjpgcts2(mle + ".step0", conffile.smooth)

	conffile.plotjpgcts2(mle + ".step1", conffile.smooth)


end

if analysis_name == "ap"

end

if analysis_name == "spot6"
	begin
		rttype = run_name.split("_")[3]
		pathalerts = PATH_RES + "/alerts/" + rttype + "_" + tstart.to_i.to_s + "_" + tstop.to_i.to_s;
		#copy .conf
		system("mkdir -p " + pathalerts);
		cmd = "cp MLE0000.conf " + pathalerts + "/" + run_name + "_MLE0000.conf"
		puts cmd
		system cmd
		cmd = "cp MLE0000.ll " + pathalerts + "/" + run_name + "_MLE0000.ll"
		puts cmd
		system cmd
		cmd = "cp MLE0000.multi " + pathalerts + "/" + run_name + "_MLE0000.multi"
		puts cmd
		system cmd
		cmd = "cp MLE0000 " + pathalerts + "/" + run_name + "_MLE0000.res"
		puts cmd
		system cmd
		cmd = "cp MLE0000.html " + pathalerts + "/" + run_name + "_MLE0000.html"
		puts cmd
		system cmd
		warningthrmin = 4
		alertthrmin_gal = 4.1
		alertthrmin_egal = 5
		Dir["MLE0000_*.source"].each do | file |
			#rttype = file.split("_")[3]
			mo = MultiOutput6.new
			mo.readDataSingleSource(file)
			pref = "_="
			if mo.sqrtTS.to_f > 3
				#create a dir with the time
				if mo.sqrtTS.to_f < warningthrmin.to_f
					pref = "_-"
				end
				if mo.b_peak.to_f > 5 or mo.b_peak.to_f < -5
					if mo.sqrtTS.to_f > alertthrmin_egal.to_f
						pref = "_+"
					end
				end
				if mo.b_peak.to_f <= 5 and mo.b_peak.to_f >= -5
					if mo.sqrtTS.to_f > alertthrmin_gal.to_f
						pref = "_+"
					end
				end

				puts "prefix: " + pref + " " + mo.b_peak.to_s + " " + mo.sqrtTS.to_s

				#system("mkdir -p " + pathalerts);

				if Dir[pathalerts + "/*.source"].size() == 0
					system("cp " + file.to_s + " " + pathalerts + "/" + pref + run_name + "_" + file);
				else
					snear = false

					puts "copy results"
					Dir[pathalerts + "/*.source"].each do | fsource |

						mo2 = MultiOutput6.new
						mo2.readDataSingleSource(fsource)
						if datautils.distance(mo2.l_peak, mo2.b_peak, mo.l_peak, mo.b_peak).to_f < 1
							snear = true
							if  mo.sqrtTS.to_f > mo2.sqrtTS.to_f
								#copy .source in the dir, appending the name of this dir
								system("cp " + file.to_s + " " + pathalerts + "/" + pref + run_name + "_" + file);
								system("rm " + fsource);
								break
							end
						end
					end
					if snear == false
						system("cp " + file.to_s + " " + pathalerts + "/" + pref + run_name + "_" + file);
					end
				end
			end
		end

	rescue
		puts "error SPOT6 results"
	end
end




if proj.to_s == "AIT"
	smooth = 7
	if binsize.to_f == 1
		smooth = 3
	end
	if binsize.to_f == 0.5
		smooth = 7
	end
    if binsize.to_f == 0.3
    	smooth = 10
    end
    if binsize.to_f == 0.2
    	smooth = 12
    end
    if binsize.to_f == 0.1
    	smooth = 15
    end
	#cmd = "export DISPLAY=localhost:3.0; ds9.rb MAP.cts.gz " + mle  + ".ctsall 2 -1 " + smooth.to_s + " B 2 jpg 1500x1000 " +  regcat.to_s;
	if File.exists?("MAP.cts.gz") and ds91 != "none"
		if ds91 == "default"
			cmd = ENV["AGILE"] + "/AGILEPIPE/spot6/ds9.rb MAP.cts.gz " + mle  + "_MAP.ctsall 1 -1 " + smooth.to_s + " B 2 png 1400x1000 ";
		else
			cmd = ENV["AGILE"] + "/AGILEPIPE/spot6/ds9.rb MAP.cts.gz " + mle  + "_MAP.ctsall " + ds91.to_s + " png 1400x1000 ";
		end
		puts "reg: " + reg + " fndisplayreg: " + fndisplayreg
		if reg == "yes" or reg == "reg" or reg == "con"
			cmd += " "
			cmd += fndisplayreg
			cmd += "."
			cmd += reg
		end
		puts cmd
		system(cmd)
	end
	if File.exists?("MAP.int.gz") and ds92 != "none"
		if ds92 == "default"
			cmd = ENV["AGILE"] + "/AGILEPIPE/spot6/ds9.rb MAP.int.gz " + mle  + "_MAP.intall 0 0.0010 " + smooth.to_s + " B 2 png 1400x1000 ";
		else
			cmd = ENV["AGILE"] + "/AGILEPIPE/spot6/ds9.rb MAP.int.gz " + mle  + "_MAP.intall " + ds92.to_s + " png 1400x1000 ";
		end
		puts "reg: " + reg + " fndisplayreg: " + fndisplayreg
		if reg == "yes" or reg == "reg" or reg == "con"
			cmd += " "
			cmd += fndisplayreg
			cmd += "."
			cmd += reg
		end
        puts cmd
        system(cmd)
	end
	if File.exists?("MAP.exp.gz") and ds93 != "none"
		if ds93 == "default"
			cmd = ENV["AGILE"] + "/AGILEPIPE/spot6/ds9.rb MAP.exp.gz " + mle  + "_MAP.expall 1 -1 1 B 2 png 1400x1000 ";
		else
			cmd = ENV["AGILE"] + "/AGILEPIPE/spot6/ds9.rb MAP.exp.gz " + mle  + "_MAP.expall " + ds93.to_s +  " png 1400x1000 ";
		end
		puts "reg: " + reg + " fndisplayreg: " + fndisplayreg
		if reg == "yes" or reg == "reg" or reg == "con"
			cmd += " "
			cmd += fndisplayreg
			cmd += "."
			cmd += reg
		end
		puts cmd
		system(cmd)
	end
end

if analysis_name != "single" and analysis_name != "spot6"
	conffile.copyresults(mle)
end

#generate auxiliary file orbit
#mjd_mjd
#tstart utc
#tstop utc
tstarttt = 0
tstartutc = 0
tstartmjd = 0
tstoptt = 0
tstoputc = 0
tstoptmjd = 0
if timetype == "TT"
	tstarttt = tstart
	tstoptt = tstop
	tstartutc = datautils.time_tt_to_utc(tstart)
	tstoputc = datautils.time_tt_to_utc(tstop)
	tstartmjd = datautils.time_tt_to_mjd(tstart)
	tstopmjd = datautils.time_tt_to_mjd(tstop)
end
if timetype == "UTC"
	tstarttt = datautils.time_utc_to_tt(tstart)
	tstoptt = datautils.time_utc_to_tt(tstop)
	tstartutc = tstart
	tstoputc = tstop
	tstartmjd = datautils.time_utc_to_mjd(tstartutc)
	tstopmjd = datautils.time_utc_to_mjd(tstoputc)
end
if timetype == "MJD"
	tstarttt = datautils.time_mjd_to_tt(tstart)
	tstoptt = datautils.time_mjd_to_tt(tstop)
	tstartutc = datautils.time_mjd_to_utc(tstart)
	tstoputc = datautils.time_mjd_to_utc(tstop)
	tstartmjd = tstart
	tstopmjd = tstop
end

forbit = File.new("orbit", "w")
forbit.write(format("%.2f", tstartmjd) + "_" + format("%.2f", tstopmjd) + "\n")
forbit.write(tstartutc.to_s + "\n")
forbit.write(tstoputc.to_s + "\n")
forbit.write(tstarttt.to_s + "\n")
forbit.write(tstoptt.to_s + "\n")
forbit.close()
