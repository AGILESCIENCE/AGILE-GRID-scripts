
class DataUtils
		def initialize
			@lines = [];
		end

		def extractlc(file, prefix, l, b, t0, t1)
			cts2 = prefix.to_s  + ".cts.gz"
			exp2 = prefix.to_s  + ".exp.gz"
			gas2 = prefix.to_s  + ".gas.gz"
			int2 = prefix.to_s  + ".int.gz"
			the2 = prefix.to_s  + ".the.gz"
			
			cts = pixel_value(cts2, l, b);
			exp = pixel_value(exp2, l, b);
			gas = pixel_value(gas2, l, b);
			int = pixel_value(int2, l, b);
			the = pixel_value(the2, l, b);
			t0m = time_tt_to_mjd(t0)
			t1m = time_tt_to_mjd(t1)
			out = t0.to_s + " " + t1.to_s + " " + t0m.to_s + " " + t1m.to_s + " " + cts.to_s + " " + exp.to_s + " " + gas.to_s + " " + int.to_s  + " " + the.to_s + "\n"
			file.write(out)
		end
		
		def mergeMaps(prefix1, prefixoutput, pattern)
			a = Dir[prefix1.to_s + "*" + pattern.to_s].sort
			index = 0;
			output = "temp" + index.to_s + pattern.to_s
			cmd = "cp " + a[0].to_s + " " + output.to_s
			system(cmd);
			output2=""
			a.each do | file |
				if index == 0
					index = index + 1
					next
				end
				output = "temp" + (index-1).to_s + pattern.to_s
				output2 = "temp" + (index).to_s + pattern.to_s
				cmd = "fimgmerge " + file.to_s + " " + output.to_s + " " + output2.to_s + " 0 0 "
				puts cmd;
				system(cmd);
				index = index + 1
			end
			cmd = "mv " + output2.to_s + " " + prefixoutput.to_s + pattern.to_s
			system(cmd)
			cmd = "rm temp*"
			system(cmd);
		end

		def execute(prefixlog, cmd)
			writeLog(prefixlog, cmd)
			puts cmd
			system(cmd)
		end
		
		def executeprefixcluster(prefix, cluster)
			if cluster.to_i == 1
				writeCommands(prefix.to_s + ".run", " ", "w")
				#writeCommands(prefix.to_s + ".clwait", "source /usr/local_agile/profile_agile_x64", "a")
				system("chmod u+x " + prefix.to_s + ".run")
			end
			if cluster.to_i == 2
				writeCommands(prefix.to_s + ".run", "#!/bin/bash", "w")
				system("chmod u+x " + prefix.to_s + ".run")
			end
			
		end
		
		def executecluster(prefixlog, cmd, cluster)
			if cluster.to_i >= 1
				writeCommands(prefixlog.to_s + ".run", cmd, "a")
			else
				writeLog(prefixlog, cmd)
			end
			puts "This is the command:"
			puts cmd
			if cluster.to_i == 0
				valret = system(cmd)
				if valret == false
					puts "EXECUTION ERROR"
					exit(1)
				end
			end
		end
		
		def writeCommands(prefix, line, mode)
			logname = prefix.to_s
			a = File.open(logname, mode)
			a.write(line.to_s + "\n")
			a.close();
		end

		def writeLog(prefix, line)
			logname = prefix.to_s + ".command.log"
			a = File.open(logname, "a")
			a.write(line.to_s + "\n")
			a.close();
		end

		def pixel_value(filename, l, b) 
			ss = rand(10000000)
			cmd2 = "/tmp/pv." + ss.to_s + ".out"
			cmd = PATH + "/bin/AG_checkMapValue " + filename.to_s + " " + l.to_s + " " + b.to_s + " > " + cmd2.to_s;
			system(cmd);
			puts cmd
			File.open(cmd2).each_line do | line |
				@pixel_value = line.to_f
			end
			cmd = "rm " + cmd2.to_s;
			system(cmd)
			@pixel_value 
		end

		def extractTimeMinMaxForContactFromFLGIndex(index_name_flg, contact) 
			begin
				#estrazione dei tempi min e max dal corfileindex
				@tmin = 0.0;
				@tmax = 0.0;
				corname = "";
				contact = contact.to_i;
				contact = format("%06d", contact);
				File.open(index_name_flg).each_line do |x|
					#puts x.size();
					if x.size() != 1
						x1 = x.split(" ");
						
						corname = x1[0];
						x2 = x1[0].split("/");
				
						cname=x2[x2.size()-1];
						x3=cname.split("_");
						x4=x3[0].split("PKP")[1];
						
						if x4 == contact
				
							if @tmin == 0.0
								@tmin = x1[1].to_f;
							end
							@tmax = x1[2].to_f;
							break;	
						end	
					end
				end
# 				puts "MIN TIME: " + @tmin.to_s + " MAX TIME: " + @tmax.to_s;
			rescue SystemCallError
				puts "DataUtils::extractTimeMinMaxForContact error";
			end

		end

		def extractTimeMinMaxForContact(index_name_cor, contact) 
			begin
				#estrazione dei tempi min e max dal corfileindex
				@tmin = 0.0;
				@tmax = 0.0;
				corname = "";
				contact = contact.to_i;
				contact = format("%06d", contact);
				File.open(index_name_cor).each_line do |x|
					#puts x.size();
					if x.size() != 1
						x1 = x.split(" ");
						
						corname = x1[0];
						x2 = x1[0].split("/");
				
						cname=x2[x2.size()-2];
						
						if cname == contact
				
							if @tmin == 0.0
								@tmin = x1[1].to_f;
							end
							@tmax = x1[2].to_f;
							break;	
						end	
					end
				end
# 				puts "MIN TIME: " + @tmin.to_s + " MAX TIME: " + @tmax.to_s;
			rescue SystemCallError
				puts "DataUtils::extractTimeMinMaxForContact error";
			end

		end
		def lines
			@lines
		end
		def extractContact(time)
			@extractContact = ""
			@ec = ""
			if @lines.size() <= 1
				File.open(BASEDIR_ARCHIVE.to_s + "/DATA/INDEX/LOG.log.index").each_line do |x1|
					@lines += [x1];
				end
				@lines.sort!
			end
			@lines.each do |x|
				c = x.split(" ");
				file = c[0];
				tmin = c[1];
				tmax = c[2];
				if (tmin.to_f <= time.to_f && time.to_f <= tmax.to_f) || time.to_f < tmax.to_f
					f1 = file.split("/");
					f2 = f1[f1.size()-1];
					f3 = f2.split(".")[0];
					@extractContact = f3.to_s
					@ec = f3
# 					puts @ec
#  					puts time.to_s + " " + tmax.to_s
					break;
# 					puts "EC " + @extractContact.to_s
				end
			end
		end
		
		def extractContactLocalIndex(index, time)
			@extractContactLocalIndex = ""
			@ec = ""
			if @lines.size() <= 1
				File.open(index).each_line do |x1|
					@lines += [x1];
				end
				@lines.sort!
			end
			@lines.each do |x|
				c = x.split(" ");
				file = c[0];
				tmin = c[1];
				tmax = c[2];
				
				#if (tmin.to_f <= time.to_f && time.to_f <= tmax.to_f) || time.to_f < tmax.to_f
				if (time.to_f >= tmin.to_f && time.to_f <= tmax.to_f) 
					
					@extractContactLocalIndex = x.split(" ")[0].split("/")[4].split(".")[0]
					@ec = @extractContactLocalIndex
					#puts x.split(" ")[0].split("/")[4].split(".")[0]
					
# 					puts @ec
#  					puts time.to_s + " " + tmax.to_s
					break;
# 					puts "EC " + @extractContact.to_s
				end
			end
		end

		def time_mjd_to_utc(mjd)
			time_mjd_to_tt(mjd)
			time_tt_to_utc(@time_mjd_to_tt)
			@time_mjd_to_utc = @time_tt_to_utc
		end
		
		def time_utc_to_mjd(utc)
			time_utc_to_tt(utc)
			time_tt_to_mjd(@time_utc_to_tt)
			@time_utc_to_mjd = @time_tt_to_mjd
		end

		def time_mjd_to_tt(mjd)
			@time_mjd_to_tt = (mjd.to_f - 53005.0) *  86400.0;
		end

		def time_tt_to_mjd(time)
			@time_tt_to_mjd = (time.to_f / 86400.0)+53005.0
		end
		
		def julian? (jd, sg)
			case sg
			when Numeric
				jd < sg
			else
				not sg
			end
		end
		
		# Convert a fractional day +fr+ to [hours, minutes, seconds,
		# fraction_of_a_second]
		def day_fraction_to_time(fr)
			ss,  fr = fr.divmod(Rational(1, 86400)) # 4p
			h,   ss = ss.divmod(3600)
			min, s  = ss.divmod(60)
			return h, min, s, fr
		end
		
		# Convert a Julian Day Number to a Civil Date.  +jd+ is
		# the Julian Day Number. +sg+ specifies the Day of
		# Calendar Reform.
		#
		# Returns the corresponding [year, month, day_of_month]
		# as a three-element array.
		def jd_to_civil(jd, sg=Date::GREGORIAN)
			if julian?(jd, sg)
				a = jd
			else
				x = ((jd - 1867216.25) / 36524.25).floor
				a = jd + 1 + x - (x / 4.0).floor
			end
			b = a + 1524
			c = ((b - 122.1) / 365.25).floor
			d = (365.25 * c).floor
			e = ((b - d) / 30.6001).floor
			dom = b - d - (30.6001 * e).floor
			if e <= 13
				m = e - 1
				y = c - 4716
			else
				m = e - 13
				y = c - 4715
			end
			return y, m, dom
		end

		def time_tt_to_utc(time)
			#TODO
			@time_tt_to_utc = -1;
			#return
			
			utc_offset = 2400000.5.to_f
			mjdref = 53005.0.to_f
			sod = 86400.0.to_f
			sec_offset = 43200.0.to_f
			
			a=jd_to_civil(utc_offset + mjdref + ( (time.to_f+sec_offset)/86400.0 ))
			b=day_fraction_to_time(a[2].modulo(1).to_f)
			utc = [a[0] , a[1] , a[2].to_i.to_f , b[0] , b[1] , b[2] , b[3]]
			utc_s = a[0].to_s+"-"+format("%02d", a[1].to_i)+"-"+format("%02d", a[2].to_i)+"T"+format("%02d", b[0])+":"+format("%02d", b[1])+":"+format("%02d", b[2]);
			#+":"+b[3].to_s
			#puts utc
			#puts "INPUT TT TIME "+ARGV[0].to_s+"      correspond to following UTC time "+utc_s

			
			@time_tt_to_utc = utc_s
			@time_tt_to_utc
		end
		
		def time_to_day_fraction(h, min, s)
			if Integer === h && Integer === min && Integer === s
			  Rational(h * 3600 + min * 60 + s, 86400) # 4p
			else
			  (h * 3600 + min * 60 + s).to_r/86400 # 4p
			end
   		end

		def time_utc_to_tt(time)
		
			@time_utc_to_tt = -1;
			
			if time == nil or  time == ""
				return
			end
			#puts time
			#ss = rand(10000000)
			ls1 = time.split("T")[0];
			ls2 = time.split("T")[1];
			year=ls1.split("-")[0];
			month=ls1.split("-")[1];
			day=ls1.split("-")[2];
			hour=ls2.split(":")[0];
			min=ls2.split(":")[1];
			sec=ls2.split(":")[2];
			
			utc_offset = 2400000.5.to_f
			mjdref = 53005.0.to_f
			sod = 86400.0.to_f
			sec_offset = 43200.0.to_f
			fod = time_to_day_fraction( hour.to_i , min.to_i , sec.to_i ).to_f
			#mjd = Date.civil_to_jd( year.to_i , month.to_i , day.to_f ).to_f + fod.to_f
			jd = DateTime.strptime(time, '%Y-%m-%dT%H:%M:%S').jd + fod.to_f
			#puts jd
			tt = (jd - utc_offset - mjdref) * sod
			tt -= sec_offset

			@time_utc_to_tt = tt
			@time_utc_to_tt
		end

		def ec
			@ec
		end

		def tmin
			@tmin
		end

		def tmax
			@tmax
		end

		def fitskeyword
			@fitskeyword
		end

		

		def extractALIKEFIXED_FLUX(input)
			index = 0;
			File.open(input).each_line do | line |
				if index == 3
					@flux = line.split(" ")[8].to_f
				end
				index = index.to_i + 1;
			end
		end

		def extractALIKEFIXED_LINE(input)
			index = 0;
			File.open(input).each_line do | line |
				if index == 3
					@alikesingle_line = line
				end
				index = index.to_i + 1;
			end
		end


		def extractALIKEFIXED_SQRTS(input)
			index = 0;
			File.open(input).each_line do | line |
				if index == 3
					@sqrts = line.split(" ")[7].to_f
				end
				index = index.to_i + 1;
			end
		end

		def flux
			@flux
		end

		def sqrts
			@sqrts
		end

		def alikesingle_line
			@alikesingle_line
		end

		def distance(ll, bl, lf, bf)
			begin
				#verify parameters
				if ll.to_f < 0 or ll.to_f > 360 or lf.to_f < 0 or lf.to_f > 360
					@distance = -2
					
				else
					if bl.to_f < -90 or bl.to_f > 90 or bf.to_f < -90 or bf.to_f > 90
						@distance = -2
					
					else

						d1 = bl.to_f - bf.to_f
						d2 = ll.to_f - lf.to_f;
				
						bl1 = Math::PI / 2.0 - (bl.to_f * Math::PI  / 180.0)
						bf1 = Math::PI / 2.0 - (bf.to_f * Math::PI  / 180.0)
						m4 = Math.cos(bl1) * Math.cos(bf1)  + Math.sin(bl1) * Math.sin(bf1) * Math.cos(d2.to_f * Math::PI  / 180.0);			
						if m4.to_f > 1 
							m4 = 1
						end
						#puts "DEBUG " + m4.to_s
						d4 = Math.acos(m4.to_f) *  180.0 / Math::PI;
						@distance = d4;
					end
				end
				
			rescue SystemCallError
				puts "ERROR IN ACOS"
				d3 = Math.sqrt(d1.to_f*d1.to_f + d2.to_f * d2.to_f);
				@distance = d3
			end
		end

		def extractFilterDir(filter)
			l0 = filter.split("_")[0];
			l1 = filter.split("_")[1];
			l2 = filter.split("_")[2];
			if l0 != nil
				@filterdir = l0;
			end
			if l1 != nil
				@filterdir = l0.to_s + "_"  + l1.to_s;
			end
		end

		def filterdir
			@filterdir
		end

		#OK
		def getResponseMatrix(filter)
			filterbase = filter.split("_")[0];
			version = filter.split("_")[2];
			if version == nil
				version = TYPE_MATRIX
			end
			@sarmatrix = "AG_GRID_G0017_S0001_I0001_NEW3.sar"
			@edpmatrix = "AG_GRID_G0017_S0001_I0001_NEW3.edp"
			@psdmatrix = "AG_GRID_G0017_S0001_I0001_NEW3.psd"
			@expcorr = "expcorr.fits.gz"
			#selezione della matrix
			
			if filterbase == "FT3ab"  
				@sarmatrix = "AG_GRID_G0017_SFT3abG_" + version.to_s + ".sar.gz"
				@edpmatrix = "AG_GRID_G0017_SFT3abG_" + version.to_s + ".edp.gz"
				@psdmatrix = "AG_GRID_G0017_SFT3abG_" + version.to_s + ".psd.gz"
			end
			if filterbase == "F4"
				@sarmatrix = "AG_GRID_G0017_SF4G_" + version.to_s + ".sar.gz"
				@edpmatrix = "AG_GRID_G0017_SF4G_" + version.to_s + ".edp.gz"
				@psdmatrix = "AG_GRID_G0017_SF4G_" + version.to_s + ".psd.gz"
			end
			if filterbase == "FT3"  
				@sarmatrix = "AG_GRID_G0017_S000FT3G_I0003.sar.gz"
				@edpmatrix = "AG_GRID_G0017_S0001_I0002.edp.gz"
				@psdmatrix = "AG_GRID_G0017_S0001_I0001_NEW3.psd"
			end
			if filterbase == "FM3.119"
				@sarmatrix = "AG_GRID_G0017_SFMG_" + version.to_s + ".sar.gz"
				@edpmatrix = "AG_GRID_G0017_SFMG_" + version.to_s + ".edp.gz"
				@psdmatrix = "AG_GRID_G0017_SFMG_" + version.to_s + ".psd.gz"
			end
		end

		def sarmatrix
			@sarmatrix
		end
		def edpmatrix
			@edpmatrix
		end
		def psdmatrix
			@psdmatrix
		end
		def expcorr
			@expcorr
		end

		#OK
		def getResponseMatrixString(filter)
			getResponseMatrix(filter);
			basicpath = PATHMODEL
			finalstring = " " + basicpath.to_s + @sarmatrix.to_s + " " + basicpath.to_s + @edpmatrix.to_s + "  " + basicpath.to_s + @psdmatrix.to_s + " ";

			@getResponseMatrixString = finalstring
		end

		def getSkyMatrix(filter, emin, emax, skytype)
			ext = ".conv."
			if Dir[PATHMODEL+"/*.disp.*"].size.to_i > 0
				ext = ".disp.conv."
			end
			emin = emin.to_i;
			emax = emax.to_i
			filterbase = filter.split("_")[0];
			version = filter.split("_")[2];
			if version == "I0025"
				skytype = 4
			end
			if version == nil
				version = TYPE_MATRIX
			end
			binsize = filter.split("_")[3];
			if binsize == nil
				binsize = "0.1"
			end
			@skymatrix = PATHMODEL + "" + format("%01d_%01d", emin, emax) + ".0.1.conv.sky ";
			@skymatrixL = @skymatrix
			@skymatrixH = @skymatrix
			binsizeL = 0.5
			if skytype.to_i == 0
				binsizeH = 0.1
			end
			if skytype.to_i == 1
				binsizeH = "gc_allsky"
			end
			if skytype.to_i == 2
				binsizeL = 0.5
				binsizeH = 0.5
			end
			if skytype.to_i == 3
				binsizeL = "SKY001"
				binsizeH = "SKY001"
			end
			if skytype.to_i == 4
				binsizeL = "SKY002"
				binsizeH = "SKY002"
			end
			
			if filterbase == "FT3ab"  
					@skymatrix = PATHMODEL + format("%01d_%01d", emin, emax) + "." + binsize.to_s + ".SFT3abG_" + version.to_s + ext.to_s + "sky.gz "
					@skymatrixL = PATHMODEL + format("%01d_%01d", emin, emax) + "." + binsizeL.to_s + ".SFT3abG_" + version.to_s + ext.to_s + "sky.gz "
					@skymatrixH = PATHMODEL + format("%01d_%01d", emin, emax) + "." + binsizeH.to_s + ".SFT3abG_" + version.to_s + ext.to_s + "sky.gz "
			end
			if filterbase == "F4"
					@skymatrix = PATHMODEL + format("%01d_%01d", emin, emax) + "." + binsize.to_s + ".SF4G_" + version.to_s + ext.to_s + "sky.gz "
			end
			if filterbase == "FT3"  
					@skymatrix = PATHMODEL + format("%01d_%01d", emin, emax) + "." + binsize.to_s + ".SFT3G_" + version.to_s + ext.to_s + "sky.gz "
			end
			if filterbase == "FM3.119"
					@skymatrix = PATHMODEL + format("%01d_%01d", emin, emax) + "." + binsize.to_s + ".SFMG_" + version.to_s + ext.to_s + "sky.gz "
					@skymatrixL = PATHMODEL + format("%01d_%01d", emin, emax) + "." + binsizeL.to_s + ".SFMG_" + version.to_s + ext.to_s + "sky.gz "
					@skymatrixH = PATHMODEL + format("%01d_%01d", emin, emax) + "." + binsizeH.to_s + ".SFMG_" + version.to_s + ext.to_s + "sky.gz "
			end
				
			
			
		end
		
		def skymatrix
			@skymatrix
		end
		
		def skymatrixL
			@skymatrixL
		end
		
		def skymatrixH
			@skymatrixH
		end

		def logindex(filterall)
			extractFilterDir(filterall)
			
			fbase = filterall.split("_");
			fbaseext = ""
			#puts fbase.size();
			for i in 2..fbase.size()
				fbaseext += "_"; fbaseext += fbase[i-1].to_s
				puts fbaseext
			end
			@logindex = ""
			if fbaseext.size() == 2
				@logindex =  BASEDIR_ARCHIVE.to_s + "/DATA" + fbaseext.to_s + "/INDEX/LOG.log.index";
			end
			if fbaseext.size() == 3
				@logindex =  BASEDIR_ARCHIVE.to_s + "/DATA" + fbaseext.to_s + "/INDEX/LOG.log_" + @filterdir.to_s + ".index";
			end
			if fbaseext.size() > 3
				@logindex =  BASEDIR_ARCHIVE.to_s + "/DATA" + fbaseext.to_s + "/INDEX/LOG.log.index";
			end
			puts @logindex
			@logindex
		end

		def evtindex(filterall)
			extractFilterDir(filterall)
			fbase = filterall.split("_");
			fbaseext = ""
			#puts fbase.size();
			for i in 2..fbase.size()
				fbaseext += "_"; fbaseext += fbase[i-1].to_s
				puts fbaseext
			end
			@evtindex = ""
			if fbaseext.size() == 2
				@evtindex =  BASEDIR_ARCHIVE.to_s + @filterdir.to_s  + "/INDEX/EVT.index";
			end
			if fbaseext.size() == 3
				@evtindex =  BASEDIR_ARCHIVE.to_s + @filterdir.to_s  + "/INDEX/EVT_" + @filterdir.to_s + ".index";
			end
			if fbaseext.size() > 3
				@evtindex =  BASEDIR_ARCHIVE.to_s + @filterdir.to_s  + "/INDEX/EVT.index";
			end
			@evtindex
		end
end
