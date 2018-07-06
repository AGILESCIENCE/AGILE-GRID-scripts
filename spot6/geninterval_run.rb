#! /usr/bin/ruby
# 0) .conffile for LC
# 1) tstart (in TT)
# 2) tend (in TT)
# 3) tinterval (in seconds)
# 4) tshift (in second)
# 5) l source
# 6) b source
# 7) source_name
# 8) analysisname
# 9) source_flux
#10) source_si
#optionals:
#11) ap_conffile (optional) conffile for AP -> only for pointing_mode=0 or pointing_mode=1
#12) pointing_mode (optional):
#0 -> consider a linear division of the time between [tstart, tend] with tinterval and tshift. if tend.to_f > 182692800 and commandline.pointing_mode == 1 then commandline.pointing_mode == 0
#1 -> if = 1 use the full OB.
#2 -> if = 2 use the bounadries of the .lc file: PATH_RES + "/" + source_name + "/GI/analysisres.lc"
#13) dbserver
#14) dbusn
#15) dbpwd
#16) dbname (to select detections in the same period)
#17) immediatesubmitrun (0 do not submit, 1 submit)
#18) gianalysisname 
#19) galcoeff (optional) -> set this galcoeff into the .conf (fixed)


load ENV["AGILE"] + "/AGILEPIPE/env.rb"
load ENV["AGILE"] + "/AGILEPIPE/spot6/Conf.rb"
load ENV["AGILE"] + "/scripts/conf.rb"

class CommandLineParam

	def initialize
		@ap_conffile = nil
		@pointing_mode = 0
		@galcoeff = nil
		@lcfile = nil
		@dbserver = nil
		@dbusn = nil
		@dbpwd = nil
		@dbname = nil
		@immediatesubmitrun = 0
		@gianalysisname = "GI"
		@addcurvedspectra = "0 0.0 0.0"
	end

	attr_accessor :ap_conffile

	attr_accessor :pointing_mode

	attr_accessor :galcoeff

	attr_accessor :lcfile

	attr_accessor :dbserver

	attr_accessor :dbusn

	attr_accessor :dbpwd

	attr_accessor :dbname

	attr_accessor :immediatesubmitrun
	
	attr_accessor :gianalysisname
	
	attr_accessor :addcurvedspectra

	def processLine(argv)
		keyw = argv.split("=")[0];
		value = argv.split("=")[1];
		puts keyw.to_s + " " + value.to_s
		case keyw

		when "ap_conffile"
			@ap_conffile = value;
		when "pointing_mode"
			@pointing_mode = value;
		when "galcoeff"
			@galcoeff = value;
		when "lcfile"
			@lcfile = value;
		when "dbserver"
			@dbserver = value;
		when "dbusn"
			@dbusn = value;
		when "dbpwd"
			@dbpwd = value;
		when "dbname"
			@dbname = value;
		when "immediatesubmitrun"
			@immediatesubmitrun = value;
		when "gianalysisname"
			@gianalysisname = value;
		when "addcurvedspectra"
			@addcurvedspectra = value;
		else
			puts "Keyword " + argv.to_s + " error."
			#exit;
		end
	end

	def processInput(startindex, s)
		for i in startindex...s.size
			if s[i] == nil
				break;
			else
				processLine(s[i]);
			end
		end
	end
end


def makeAP(source_name, ap_conffile, tstart, tend, lsource, bsource, tshift, abspath, ob, analysisname,commandline)
	newconffilename = "/tmp/ob" + (rand * 1000000000).to_i.to_s + ".conf"

	system("ruby -p -e \"gsub(/SOURCENAME/, '"+source_name+"')\" " + ap_conffile + " > " + newconffilename + ".bak")
	newconffilename2 = "/tmp/ob" + (rand * 1000000000).to_i.to_s + ".conf"
	if analysisname != nil
		system("ruby -p -e \"gsub(/ANALYSISNAME/, '"+analysisname+"')\" " + newconffilename + ".bak" + " > " + newconffilename2 + ".bak")
		system("mv " + newconffilename2 + ".bak " + newconffilename + ".bak")
	end

	newconffilename2 = "/tmp/ob" + (rand * 1000000000).to_i.to_s + ".conf"
	system("ruby -p -e \"gsub(/SOURCEL/, '"+lsource.to_s+"')\" " + newconffilename + ".bak" + " > " + newconffilename2 + ".bak")
	system("mv " + newconffilename2 + ".bak " + newconffilename + ".bak")

	newconffilename2 = "/tmp/ob" + (rand * 1000000000).to_i.to_s + ".conf"
	system("ruby -p -e \"gsub(/SOURCEB/, '"+bsource+"')\" " + newconffilename + ".bak" + " > " + newconffilename2 + ".bak")
	system("mv " + newconffilename2 + ".bak " + newconffilename + ".bak")


	if ob != nil
		apob = "AP_" + ob.to_s
		newconffilename2 = "/tmp/ob" + (rand * 1000000000).to_i.to_s + ".conf"
		system("ruby -p -e \"gsub(/AP/, '"+apob+"')\" " + newconffilename + ".bak > " + newconffilename2 + ".bak")

		system("mv " + newconffilename2 + ".bak " + newconffilename + ".bak")

	end
	puts newconffilename
	fo = File.new(newconffilename, "w")

	index = 0
	File.open(newconffilename + ".bak").each_line do | line |
		out = line

		if index.to_i == 2
			out = tstart.to_f.to_s  + "\n"
		end
		if index.to_i == 3
			out = tend.to_f.to_s + "\n"
		end
		if index.to_i == 4
			out = "TT\n"
		end
		#if index.to_i == 5
		#	out = lsource.to_s + "\n"
		#end
		#if index.to_i == 6
		#	out = bsource.to_s + "\n"
		#end

		if index.to_i == 10
			out = out.chomp + " timeslot=" + tshift.to_s + "\n"
		end

		fo.write(out);
		index = index + 1
	end
	fo.close()

	system("mv " + newconffilename + " " + abspath + "/commands/");
	if commandline.immediatesubmitrun.to_i == 1
		cmd = ENV["AGILE"] + "/AGILEPIPE/spot6/immediatesubmitrun.rb " + abspath + "/commands/" + newconffilename.split("/").last
		puts cmd
		system cmd
	end
	system("rm " + newconffilename + ".bak")
end

if ARGV[0].to_s == "help" || ARGV[0].to_s == "h" || ARGV[0] == nil
	system("head -27 " + $0 );
	exit;
end

datautils = DataUtils.new

conffile = ARGV[0]
tstart = ARGV[1]
tend = ARGV[2]
tinterval = ARGV[3]
tshift = ARGV[4]
lsource = ARGV[5]
bsource = ARGV[6]
source_name = ARGV[7]
analysisname = ARGV[8]
source_flux = ARGV[9]
source_si = ARGV[10]

commandline =  CommandLineParam.new
commandline.processInput(11, ARGV)

if tend.to_f > 182692800 and commandline.pointing_mode == 1
	commandline.pointing_mode = 0
end

abspath=PATH_RES

if commandline.pointing_mode.to_i == 0
	tcurrent = tstart.to_f + tinterval.to_f
end

if commandline.pointing_mode.to_i == 1
	tcurrent = tstart.to_f
end

if commandline.pointing_mode.to_i == 2
	tcurrent = tstart.to_f
end

indexrun = 0
while tcurrent.to_f <= tend.to_f

	lpointing = -1
	bpointing = -1

	if commandline.pointing_mode.to_i == 0
		tstarttmp = tcurrent.to_f - tinterval.to_f
		tendtmp = tcurrent.to_f
	end

	if commandline.pointing_mode.to_i == 1
		puts "############ pointing mode 1 -> use full OB" + tendtmp.to_s
		tstartob = -1
		tendob = -1
		tstarttmp = -1
		tendtmp = -1
		File.open(ENV["AGILE"] + "/scripts/AGILEOBPOINT.list").each_line do | line |
			#puts line
			ll = line.split(" ")
			if ll[6].to_f < tcurrent.to_f
				next
			end
			if tcurrent.to_f >= ll[5].to_f  and tcurrent.to_f <= ll[6].to_f
				tstartob = ll[5]
				tendob = ll[6]
				lpointing = ll[1]
				bpointing = ll[2]
				if datautils.distance(lsource, bsource, lpointing, bpointing) < 60
					tcurrent = tendob.to_f + 1
					tstarttmp = tstartob
					tendtmp = tendob
					if commandline.ap_conffile != nil

						makeAP(source_name, commandline.ap_conffile, tstartob, tendob, lsource, bsource, tshift, abspath, ll[0].to_s, analysisname,commandline)

					end

					break
				else
					tcurrent = tendob.to_f + 1
				end
			end
		end
	end

	if commandline.pointing_mode.to_i == 2
		tstartob = -1
		tendob = -1
		tstarttmp = -1
		tendtmp = -1
		#get lpointing and bpointing
		File.open(ENV["AGILE"] + "/scripts/AGILEOBPOINT.list").each_line do | line |
			ll = line.split(" ")
			if ll[6].to_f < tcurrent.to_f
				next
			end
			if tcurrent.to_f >= ll[5].to_f  and tcurrent.to_f <= ll[6].to_f
				tstartob = ll[5]
				tendob = ll[6]
				lpointing = ll[1]
				bpointing = ll[2]
				
				break;
			end
		end
		tfound = false
		puts "pointing: " + lpointing.to_s + " " + bpointing.to_s

		lcfile = PATH_RES + "/" + source_name + "/"+commandline.gianalysisname+"/analysisres.lc"

		File.open(lcfile).each_line do | line |
			ll = line.split(" ")
			tstartob = ll[13]
			tendob = ll[14]

			if tendob.to_f < tcurrent.to_f
				next
			end
			if tcurrent.to_f >= tstartob.to_f  and tcurrent.to_f <= tendob.to_f
				#found ob
				#if (tcurrent.to_f - tstartob.to_f) < tshift.to_f
				#	tcurrent = tstartob
				#end
				tfound = true
				tstarttmp = tcurrent.to_f - tinterval.to_f
				tendtmp = tcurrent.to_f
				if tstarttmp.to_f < tstartob.to_f
					tstarttmp = tstartob
					tendtmp = tstarttmp + tinterval
					if commandline.ap_conffile != nil

						makeAP(source_name, commandline.ap_conffile, tstartob, tendob, lsource, bsource, tshift, abspath, ll[6].to_s, analysisname,commandline)

					end
				end
				if tendtmp.to_f > tendob.to_f
					if commandline.pointing_mode.to_i == 2
						#analizza i dati solo se l'ob contiene tutti il periodo temporale (per il pointing)
						tendtmp = -1
					end
				end

				#tcurrent = tendtmp
			end

		end
		puts tfound
		if tfound == false

			tcurrent = tcurrent.to_f + tshift.to_f
			puts "not found ob " + tcurrent.to_s
			next
		end
		#exit with tcurrent, tstarttmp and tendtmp
	end


	puts tcurrent

	if tstarttmp.to_i == tendtmp.to_i
		tcurrent = tcurrent.to_f + tshift.to_f
		puts "discard conf code 1"
		next
	end

	if tendtmp.to_i == -1
		puts "discard conf code 2"
		break
	end

	newconffilename = "/tmp/ob" + (rand * 1000000000).to_i.to_s + ".conf"
	puts newconffilename
	system("ruby -p -e \"gsub(/SOURCENAME/, '"+source_name+"')\" " + conffile + " > " + newconffilename + ".bak")
	
	newconffilename2 = "/tmp/ob" + (rand * 1000000000).to_i.to_s + ".conf"
	system("ruby -p -e \"gsub(/ANALYSISNAME/, '"+analysisname+"')\" " + newconffilename + ".bak" + " > " + newconffilename2 + ".bak")
	system("mv " + newconffilename2 + ".bak " + newconffilename + ".bak")
	
	newconffilename2 = "/tmp/ob" + (rand * 1000000000).to_i.to_s + ".conf"
	system("ruby -p -e \"gsub(/ANALYSISGI/, '"+commandline.gianalysisname+"')\" " + newconffilename + ".bak" + " > " + newconffilename2 + ".bak")
	system("mv " + newconffilename2 + ".bak " + newconffilename + ".bak")

	newconffilename2 = "/tmp/ob" + (rand * 1000000000).to_i.to_s + ".conf"
	system("ruby -p -e \"gsub(/SOURCEL/, '"+lsource.to_s+"')\" " + newconffilename + ".bak" + " > " + newconffilename2 + ".bak")
	system("mv " + newconffilename2 + ".bak " + newconffilename + ".bak")

	newconffilename2 = "/tmp/ob" + (rand * 1000000000).to_i.to_s + ".conf"
	system("ruby -p -e \"gsub(/SOURCEB/, '"+bsource+"')\" " + newconffilename + ".bak" + " > " + newconffilename2 + ".bak")
	system("mv " + newconffilename2 + ".bak " + newconffilename + ".bak")

	newconffilename2 = "/tmp/ob" + (rand * 1000000000).to_i.to_s + ".conf"
	system("ruby -p -e \"gsub(/SOURCEFLUX/, '"+source_flux+"')\" " + newconffilename + ".bak" + " > " + newconffilename2 + ".bak")
	system("mv " + newconffilename2 + ".bak " + newconffilename + ".bak")

	newconffilename2 = "/tmp/ob" + (rand * 1000000000).to_i.to_s + ".conf"
	system("ruby -p -e \"gsub(/SOURCESI/, '"+source_si+"')\" " + newconffilename + ".bak" + " > " + newconffilename2 + ".bak")
	system("mv " + newconffilename2 + ".bak " + newconffilename + ".bak")
	
	newconffilename2 = "/tmp/ob" + (rand * 1000000000).to_i.to_s + ".conf"
	system("ruby -p -e \"gsub(/SOURCECS/, '"+commandline.addcurvedspectra+"')\" " + newconffilename + ".bak" + " > " + newconffilename2 + ".bak")
	system("mv " + newconffilename2 + ".bak " + newconffilename + ".bak")

	fo = File.new(newconffilename, "w")

	index = 0

	File.open(newconffilename + ".bak").each_line do | line |
		out = line
		if index.to_i == 2
			out = tstarttmp.to_f.to_s  + "\n"
		end
		if index.to_i == 3
			out = tendtmp.to_f.to_s + "\n"
		end
		if index.to_i == 4
			out = "TT\n"
		end
		#if index.to_i == 5
		#	out = lsource.to_s + "\n"
		#end
		#if index.to_i == 6
		#	out = bsource.to_s + "\n"
		#end
		if index.to_i == 8
			if commandline.galcoeff != nil
				out = commandline.galcoeff.to_s + "\n"
			end
		end
		if index.to_i == 10
			out = out.chomp + " lpointing=" + lpointing.to_s + " bpointing=" + bpointing.to_s + "\n"
		end
		if (index.to_i == 11 || index.to_i == 12)
			if line.strip == "spot6db"
				# spot6db l b r DB-LCname fixflag_of_analysis minexp maxexpratio minsqrtts minflux server usn pwd dbname
				out = "spot6db " + lsource.to_s + " " + bsource.to_s + " 10 " + analysisname + " 1 50 10  4 20e-08 " + commandline.dbserver + " " + commandline.dbusn + " " + commandline.dbpwd + " " + commandline.dbname + "\n";
			end
		end

		if index.to_i == 25
			out = line.chomp + "_" + format("%.2f", tstarttmp.to_f.to_s) + "_" + format("%.2f", tendtmp.to_f.to_s) + "_" + format("%06d", indexrun.to_s) + "\n"
			#abspath += out.chomp + "/"
		end

		#if index.to_i == 29
		#	lll = out.split(" ")
		#	if lll.size > 1
		#		#change the only source to be analysed
		#		out = lll[0].to_s + " " + lsource.to_s + " " + bsource.to_s + " " + lll[3].to_s + " " + lll[4].to_s + " " + lll[5].to_s + " " + lll[6].to_s
		#		if lll.size == 8
		#			out = out.chomp + " " + lll[7].to_s
		#		end
		#		out = out + "\n"
		#	end
		#end

		fo.write(out);
		index = index + 1

	end
	fo.close()

	system("mv " + newconffilename + " " + abspath + "/commands/");
	if commandline.immediatesubmitrun.to_i == 1
		cmd = ENV["AGILE"] + "/AGILEPIPE/spot6/immediatesubmitrun.rb " + abspath + "/commands/" + newconffilename.split("/").last
		puts cmd
		system cmd
	end
	system("rm " + newconffilename + ".bak")

	if commandline.pointing_mode.to_i == 0
		tcurrent = tcurrent.to_f + tshift.to_f
	end
	if commandline.pointing_mode.to_i == 2
		tcurrent = tcurrent.to_f + tshift.to_f
	end
	indexrun = indexrun.to_i + 1

end

if commandline.ap_conffile != nil
	obname = nil
	if commandline.pointing_mode.to_i == 0
		if tend.to_f > 182692800
			obname = "OB09000"
			makeAP(source_name, commandline.ap_conffile, tstart, tend, lsource, bsource, tshift, abspath, obname, analysisname,commandline)
		end
		#makeAP(source_name, commandline.ap_conffile, tstart, tend, lsource, bsource, tshift, abspath, nil)
	end

end
