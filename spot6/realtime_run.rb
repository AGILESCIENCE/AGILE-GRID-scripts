#! /usr/bin/ruby

load ENV["AGILE"] + "/AGILEPIPE/env.rb"
load ENV["AGILE"] + "/AGILEPIPE/spot6/Conf.rb"
load ENV["AGILE"] + "/scripts/conf.rb"

def sendmail(subject, bodyfile)
	agilepipeconf = ENV["AGILE"] + "/AGILEPIPE/agilepipe.conf"
	recipients = []
	lines = File.open(agilepipeconf).to_a
	sendit = false
	for line in lines do
		if line.include?("enable_email =")
			sendmail = line.split("=")[1].gsub!(/\s+/,'')
			if sendmail.downcase.include?("true")
				sendit = true
			end
		end
		if line.include?("from =")
			from = line.split("=")[1].gsub!(/\s+/,'')
		end
		if line.include?("to =")
			recipients << line.split("=")[1].gsub!(/\s+/,'').split(",")
			to = line.split("=")[1].gsub!(/\s+/,'').split(",")[0]
		end
		if line.include?("cc =")
			recipients << line.split("=")[1].gsub!(/\s+/,'').split(",")
		end
		if line.include?("bcc_spot6 =")
			recipients << line.split("=")[1].gsub!(/\s+/,'').split(",")
		end
		if line.include?("user =")
			user = line.split("=")[1].gsub!(/\s+/,'')
		end
		if line.include?("password =")
			pass = line.split("=")[1].gsub!(/\s+/,'')
		end
		if line.include?("smtp =")
			smtp = line.split("=")[1].gsub!(/\s+/,'')
		end
	end
	if sendit == true
		puts "Sending email '" + subject + "' to " + to
		cmd = "mail -v -s \"" + subject + "\" -S smtp-use-starttls -S ssl-verify=ignore -S smtp-auth=login -S smtp=" + smtp +
			  " -S from=" + from + " -S smtp-auth-user=" + from + " -S smtp-auth-password=" + pass + " -S nss-config-dir=~/.certs" +
			  " -S ssl-verify=ignore " + recipients.join(" ") + " < " + bodyfile
		system(cmd)
	end
end

def runait(lastcontacttime, day, hours_shift)

	tstart = 0
	tstop = 0
	abspath=PATH_RES

	#read last processed time
	lastprocessing2 = 0
	if File.exists?(abspath + "/commands/lastprocessing_aitoff_rt"+format("%02i", day)+"")
		File.open(abspath + "/commands/lastprocessing_aitoff_rt"+format("%02i", day)+"", "r").each_line do | line |
			lastprocessing2 = line.to_i
		end
	else
		fout = File.new(abspath + "/commands/lastprocessing_aitoff_rt"+format("%02i", day)+"", "w")
		fout.write(lastcontacttime.to_s + "\n")
		fout.close();
	end

	if lastprocessing2.to_i == 0
		tstart = lastcontacttime.to_i - 86400 * day.to_i
		tstop = lastcontacttime.to_i
		lastprocessing2 = tstop
	end

	if REALTIMEMODE.to_i == 1
		tstop = lastprocessing2 + 3600 * hours_shift.to_i
		tstart = tstop.to_i - 86400 * day.to_i
		while (lastprocessing2.to_i + 3600 * hours_shift.to_i).to_i <= lastcontacttime.to_i
			tstop = lastprocessing2 + 3600 * hours_shift.to_i
			tstart = tstop.to_i - 86400 * day.to_i
			lastprocessing2 = tstop
		end
	else
		tstop = lastprocessing2 + 3600 * hours_shift.to_i
		tstart = tstop.to_i - 86400 * day.to_i
		lastprocessing2 = tstop
	end

	puts "TIME: " + lastcontacttime.to_s + " - " + lastprocessing2.to_s + " - " + tstart.to_s + " " + tstop.to_s

	if tstop.to_i <= lastcontacttime.to_i
		#change and copy the card
		confoutfilename = "/tmp/aitoff_rt"+format("%02i", day)+"_"+rand(10000000).to_i.to_s+".conf"
		fo = File.new(confoutfilename, "w")
		index = 0
		File.open(ENV["AGILE"] + "/AGILEPIPE/spot6/cards/ait_rt/aitoff_rt.conf").each_line do | line |
			out = line
			if index.to_i == 1
				out1 = line.split(",")[0]

				out = out1.split("_")[0] + "_" + ARCHIVE.to_s + "_" + out1.split("_")[2].chomp
				out += ","
				if out1.size >= 2
					out += line.split(",")[1].chomp
					out += format("%02i", day)
				else
					out += "none"
				end
				out += ",ALLSKY"
				out += "\n"
			end
			if index.to_i == 2
				out = tstart.to_i.to_s + "\n"
			end
			if index.to_i == 3
				out = tstop.to_i.to_s + "\n"
			end
			if index.to_i == 23
				outdir = line.split(",")[0]
				outqueue = SPOT6QUEUE.to_s
				outmodule = MODULELOAD.to_s
				out = outdir.to_s + "," + outqueue.to_s + "," + outmodule + "\n"
			end
			if index.to_i == 25
				out = line.chomp + format("%02i", day) + "_" + tstart.to_i.to_s + "_" + tstop.to_i.to_s + "\n"
				#abspath += out.chomp + "/"
			end
			fo.write(out);
			index = index + 1
		end
		fo.close()
		lastprocessing2 = tstop

		puts "New run AIT: " + "cp " + confoutfilename.to_s + " " + abspath + "/commands/";
		system("mv " + confoutfilename.to_s + " " + abspath + "/commands/");

		runspot6(lastcontacttime, day, hours_shift, tstart, tstop);

		f = File.new(abspath + "/commands/lastprocessing_aitoff_rt"+format("%02i", day), "w")
		f.write(lastprocessing2);
		f.close();

	end
end

def runspot6(lastcontacttime, day, hours_shift, tstart, tstop)

	abspath=PATH_RES

	#change and copy the card
	indexfile = 0;
	Dir[ENV["AGILE"] + "/AGILEPIPE/spot6/cards/"+SPOT6CARDDIR.to_s+"/*.conf"].sort.each do | file |
		indexring = 0;
		File.open(ENV["AGILE"] + "/AGILEPIPE/spot6/cards/"+SPOT6CARDDIR.to_s+"/ringsB.coord").each_line do | coords |
			if coords.split(" ").size < 6
				next
			end
			outfileconf = "/tmp/spot6_" + format("%02i", day) + "_" + format("%02d", indexfile) + "_" + coords.split(" ")[0].chomp + "_"+rand(10000000).to_i.to_s+".conf";

			fo = File.new(outfileconf, "w")
			index = 0
			File.open(file).each_line do | line |
				out = line
				if index.to_i == 1
					out1 = line.split(",")[0]

					out = out1.split("_")[0] + "_" + ARCHIVE.to_s + "_" + out1.split("_")[2].chomp
					out += ","
					if out1.size >= 2
						out += line.split(",")[1].chomp
					else
						out += "none"
					end
					out += ","
					out += coords.split(" ")[0].chomp
					out += "\n"
				end
				if index.to_i == 2
					out = tstart.to_i.to_s + "\n"
				end
				if index.to_i == 3
					out = tstop.to_i.to_s + "\n"
				end
				if index.to_i == 4
					out = "TT\n"
				end
				if index.to_i == 5
					out = coords.split(" ")[1].to_f.to_s + "\n"
				end
				if index.to_i == 6
					out = coords.split(" ")[2].to_f.to_s + "\n"
				end
				if index.to_i == 7
					out = coords.split(" ")[5].to_s + "\n"
				end
				if index.to_i == 10
					out = line.chomp + " mapsize=" + coords.split(" ")[3].chomp + " binsize=" + coords.split(" ")[4].chomp + "\n"
				end
				if index.to_i == 23
					outdir = line.split(",")[0]
					outqueue = SPOT6QUEUE.to_s
					outmodule = MODULELOAD.to_s
					out = outdir.to_s + "," + outqueue.to_s + "," + outmodule + "\n"
				end
				if index.to_i == 25
					out = line.chomp + format("%02i", day) + "_" + tstart.to_i.to_s + "_" + tstop.to_i.to_s + "_" + format("RING%02d", indexring) + "\n"
				end
				fo.write(out);
				index = index + 1
			end
			fo.close()
			puts "New run SPOT6: " + "cp " + outfileconf + " " + abspath + "/commands/";
			system("mv " + outfileconf + " " + abspath + "/commands/");
			indexring = indexring.to_i + 1
		end
		indexfile = indexfile.to_i + 1
	end
end

def existsFile(filename)
	if filename == nil
		return ""
	end
	if File.exists?(filename)
		return filename
	else
		return ""
	end
end

def genaitoffspot6(rttype)
	if SPOT6GEN_HTMLANDPNG.to_i == 0
		return
	end
	abspath=PATH_RES + "/aitoff_rt/"
	apppath=PATH_RES + "/app/"
	b02=Dir[abspath + "*" + rttype + "*/orbit"].sort()
	last02 = b02[b02.size() - 1].split("orbit")[0]
	pathaitoff = last02 + "/MAP.cts.gz";
	pathaitoffint = last02 + "/MAP.int.gz";
	if File.exists?(pathaitoff)
		#build path
		aitname = last02.split("/")[last02.split("/").size - 1]
		aitsplitname = aitname.split("_");
		pathalerts = PATH_RES + "/alerts/" + aitsplitname[1] + "_" + aitsplitname[2] + "_" + aitsplitname[3] + "/"
		#find the AITOFF - questo va rimosso da qui e messo in un task a parte che fa la scansione della dir alert e prende l'ultimo

		begin
			mout = MultiOutputList.new
			mout.readSourcesInDir(pathalerts, "spot6", "SPOT6", 4);
		rescue
			puts "gen .reg and .html error"
		end

		cmd = ENV["AGILE"] + "/AGILEPIPE/spot6/ds9.rb " + pathaitoff + " " + pathalerts + "/" + rttype + "spot6.ctsall 1 -1 7 B 2 png 1400x1000 " + existsFile(pathalerts + "/spot6.reg");
		puts cmd
		system(cmd)
		system("cp " + pathalerts + "/" + rttype + "spot6.ctsall.png " + apppath + "lastaitspot6_"+rttype+".png")

		system("cp " + pathalerts + "/spot6.html " + apppath + "lastaitspot6_"+rttype+".html")
		system("cp " + pathalerts + "/spot6.sel.html " + apppath + "lastaitspot6_"+rttype+".sel.html")
		system("cp " + b02[b02.size() - 1] + " " + apppath);

		if File.exists?(pathaitoffint)
			cmd = ENV["AGILE"] + "/AGILEPIPE/spot6/ds9.rb " + pathaitoffint + " " + pathalerts + "/" + rttype + "spot6.intall 0 0.0010 7 B 2 png 1400x1000 " + existsFile(pathalerts + "/spot6.reg");
			puts cmd
			system(cmd)
			system("cp " + pathalerts + "/" + rttype + "spot6.intall.png " + apppath + "lastaitspot6_"+rttype+".int.png")
		end

		puts("send alerts + " + pathalerts + "/_+*")
		begin
			puts Dir[pathalerts + "/_+*"]
			Dir[pathalerts + "/_+*"].each do | file |
				nfile = file.sub("_+", "++")
				system("mv " + file + " " + nfile);

				basefile = nfile.split("_MLE0000")[0].split("++")[1]
				basedirres="http://agilepipe.iasfbo.inaf.it/analysis/spot6/" + basefile + "/"

				mout = MultiOutput6.new
				mout.readDataSingleSource(nfile)
				mout.assoccat(",");

				subject = "S6" + basefile.split("_")[3] + " " + format("%.2f", mout.sqrtTS) + " (" + format("%.2f", mout.l_peak) + "," + format("%.2f", mout.b_peak) + "," + format("%.2E", mout.exposure) + ") " + format("[%.2f-%.2f]", mout.timestart_mjd, mout.timestop_mjd) + " " + format("%.1E", mout.flux) + "+/-" + format("%.1E", mout.flux_error) + " " + mout.assoc.to_s

				system("cat " + nfile + " > alert")
				system("echo \" \" >> alert")
				puts "cat /ANALYSIS3/spot6/" + basefile + "/MLE0000 >> alert"
				system("cat /ANALYSIS3/spot6/" + basefile + "/MLE0000 >> alert")
				falert = File.open("alert", "a")
				falert.write("\n")

				stringoutput = "Dir: " + basedirres + "\n"
				falert.write(stringoutput);

				stringoutput = "Img S1: " + basedirres + "MLE0000.step1_MAP.cts2.png \n"
				falert.write(stringoutput);
				stringoutput = "Img S2: " + basedirres + "MLE0000_MAP.cts2.png \n"
				falert.write(stringoutput);
				stringoutput = "Img S2: " + basedirres + "MLE0000_MAP.ctsall.png \n"
				falert.write(stringoutput);
				stringoutput = "Img S2: " + basedirres + "MLE0000_MAP.intall.png \n"
				falert.write(stringoutput);

				falert.write("\n")
				falert.write(mout.assoc.to_s)
				falert.write("\n")

				stringoutput = "\nSIMBAD: "
		stringoutput += "http://simbad.u-strasbg.fr/simbad/sim-coo?CooDefinedFrames=none&CooEpoch=2000&Coord=";
				stringoutput += format("%.2f", mout.l_peak);
				if mout.b_peak.to_f > 0
						stringoutput += "%2b";
				end
				stringoutput += format("%.2f", mout.b_peak);
				stringoutput += "&submit=submit%20query&Radius.unit=arcmin&CooEqui=2000&CooFrame=Gal&Radius=60"
				falert.write(stringoutput);

				stringoutput = "\nNED: ";
				stringoutput += "http://nedwww.ipac.caltech.edu/cgi-bin/nph-objsearch?in_csys=Galactic&in_equinox=J2000.0&lon=" + format("%.2f", mout.l_peak) + "&lat=" + format("%.2f", mout.b_peak) + "&radius=60&search_type=Near+Position+Search&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=Distance+to+search+center&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp=YES&z_constraint=Unconstrained&z_value1=&z_value2=&z_unit=z&ot_include=ANY&nmp_op=ANY"
				stringoutput += "\n";
				falert.write(stringoutput);

				falert.write("\nAGILE SPOT6 alert system for gamma-ray transients developed by Andrea Bulgarelli (INAF)\n")

				falert.close()

				datautils = DataUtils.new
				#distance from Vela e Geminga
				if datautils.distance(mout.l_peak, mout.b_peak, 263.5, -2.8) > 1 and datautils.distance(mout.l_peak, mout.b_peak, 195.13, 4.26) > 1
					sendmail(subject, "alert")
				end
			end
		rescue
			puts "mail problem"
		end
	end
end

begin
	#b=1

	#while b == 1

		#crea le dir necessarie se non esistono

		#alerts
		alertspath=PATH_RES + "/alerts/"
		if File.exists?(alertspath) == false
			system("mkdir " + alertspath);
			system("chmod -R g+w " + alertspath);
		end

		#spot6
		spot6path=PATH_RES + "/spot6/"
		if File.exists?(spot6path) == false
			system("mkdir " + spot6path);
			system("chmod -R g+w " + spot6path);
		end

		#aitoff_rt
		abspath=PATH_RES + "/aitoff_rt/"
		if File.exists?(abspath) == false
			system("mkdir " + abspath);
			system("chmod -R g+w " + abspath);
		end

		#app
		apppath=PATH_RES + "/app/"
		if File.exists?(apppath) == false
			system("mkdir " + apppath);
			system("chmod -R g+w " + apppath);
		end

		#cmd = "sort --key=3 " + PATH_DATA + "/FM3.119_" + ARCHIVE + "/INDEX/EVT.index | tail -1 > " + PATH_RES + "/commands/lastorbit "
		#puts cmd
		#system(cmd)
		lastcontacttime = 0

		#cmd = "sort --key=3 " + PATH_DATA + "/FM3.119_" + ARCHIVE + "/INDEX/EVT.index | head -1 > " + PATH_RES + "/commands/firstorbit "
		#puts cmd
		#system(cmd)

		File.open(PATH_RES + "/commands/lastorbit", "r").each_line do | line |
			if line.split(" ").size >= 3
				lastcontacttime = line.split(" ")[2].to_i
			else
				lastcontacttime = line.to_i
			end
		end
		puts "lastcontacttime " + lastcontacttime.to_s

		#run02
		runait(lastcontacttime, 2, 1);

		#run01
		runait(lastcontacttime, 1, 1);

		#run04
		runait(lastcontacttime, 4, 8);

		#run07
		runait(lastcontacttime, 7, 12);

		#copy aitoff
		begin
			b02=Dir[abspath + "*RT02*/orbit"].sort()
			if b02.size() > 1
				last02 = b02[b02.size() - 1].split("orbit")[0]
				system("cp " + last02 + "/MLE000*.ctsall.png " + apppath + "lastait2.png")
				system("cp " + last02 + "/orbit " + apppath + "lastait2.orb")
				#system("cp " + last02 + "/MLE000*.ctsall.png " + apppath + "public.png")
				#system("cp " + last02 + "/MLE000*.ctsall.png " + apppath + "public.jpg")
				system("cp " + last02 + "/orbit " + apppath + "public.orb")
			end

			b01=Dir[abspath + "*RT01*/orbit"].sort()
			if b01.size() > 1
				last01 = b01[b01.size() - 1].split("orbit")[0]
				system("cp " + last01 + "/MLE000*.ctsall.png " + apppath + "lastait1.png")
				system("cp " + last01 + "/orbit " + apppath + "lastait1.orb")
			end

			b04=Dir[abspath + "*RT04*/orbit"].sort()
			if b04.size() > 1
				last04 = b04[b04.size() - 1].split("orbit")[0]
				system("cp " + last04 + "/MLE000*.ctsall.png " + apppath + "lastait4.png")
				system("cp " + last04 + "/orbit " + apppath + "lastait4.orb")
			end

			b07=Dir[abspath + "*RT07*/orbit"].sort()
			if b07.size() > 1
				last07 = b07[b07.size() - 1].split("orbit")[0]
				system("cp " + last07 + "/MLE000*.ctsall.png " + apppath + "lastait7.png")
				system("cp " + last07 + "/orbit " + apppath + "lastait7.orb")
			end
		rescue
			puts "warning: no file found"
		end

		begin
			genaitoffspot6("RT02")
			genaitoffspot6("RT01")
			genaitoffspot6("RT04")
			genaitoffspot6("RT07")
		rescue
			puts "error or warning in file system "
		end

		#copy images to agile.iasfbo.inaf.it
		#system("scp " + apppath + "* marlin:/var/www/html/AGILEApp/RT/")

	#end
end
