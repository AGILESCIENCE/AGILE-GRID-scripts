


class AgileFOV
	#14 contatti per ogni giornata
	# N contatti = N day * 14
	# end contact = start contact + N contact -1
	# end contact = start contact + N day * 14 - 1
	public
		def initialize() 
			@pathbase = BASEDIR_ARCHIVE + "/DATA/INDEX/LB.index";
		end

		def distanceFromPointing(lp, bp, mindistance, maxdistance)
			begin
				datautils = DataUtils.new
				File.open(@pathbase).each_line do |x|
					cc = x.split(" ")[0];
					lc = x.split(" ")[1];
					lce = x.split(" ")[2];
					if lce == "nan"
						lce = -1
					end
					bc = x.split(" ")[3];
					bce = x.split(" ")[4];
					if bce == "nan"
						bce = -1
					end
					dist = datautils.distance(lc, bc, lp, bp).to_f;
					if (	dist <= maxdistance.to_f &&
						dist >= mindistance.to_f)
						puts cc.to_s + "\t" + format("%.2f", lc) + "\t(" + format("%.2f", lce) +  ")\t" + format("%.2f", bc) + "\t(" + format("%.2f", bce) + ")\t" + format("%.2f", dist).to_s;
					end
				end
			rescue SystemCallError
				puts "distanceFromPointing error"
			end
		end

		def getRow(c)
			begin
				File.open(@pathbase).each_line do |x|
					cc = x.split(" ")[0];
# 					puts cc.to_s + " " + c.to_s
					if cc.to_i == c.to_i
						@line = x;
						break;
					end
				end
			rescue SystemCallError
				@line = "";
			end
		end

		def longitudeFromPeriod(tstart, tend)
			datautils = DataUtils.new
			ss = rand(10000000)
			fileout = "/tmp/coordinates." + ss.to_s + ".out"

			cmd = ENV["AGILE"]+"/scripts/lb_log   " + BASEDIR_ARCHIVE + "/DATA/INDEX/LOG.log.index " + tstart.to_s + " " + tend.to_s + " > " + fileout.to_s
			puts cmd
			system(cmd);

			lon = 0;
			File.open(fileout).each_line do | line |
				if line.split("=").size() == 2
					lon = line.split("=")[1].strip;
					break;
				end
			end

			cmd = "rm " + fileout.to_s; system(cmd)

			@longitudeFromPeriod = lon;
		end

		def latitudeFromPeriod(tstart, tend)
			datautils = DataUtils.new
			ss = rand(10000000)
			fileout = "/tmp/coordinates." + ss.to_s + ".out"

			cmd = ENV["AGILE"]+"/scripts/lb_log   /AGILE_DATA/INDEX/LOG.log.index " + tstart.to_s + " " + tend.to_s + " > " + fileout.to_s
			system(cmd);

			lat = 0;
			File.open(fileout).each_line do | line |
				if line.split("=").size() == 2
					lat = line.split("=")[1].strip;
				end
			end

			cmd = "rm " + fileout.to_s; system(cmd)

			@latitudeFromPeriod = lat;
		end

		def longitudeFromPeriod2(tstart, tend)
			begin
				timef = 0
				File.open(ENV["AGILE"]+"/scripts/AGILEPOINTING").each_line do | line |
					time = line.split(" ");
					@longitudeFromPeriod2 = -999
					if tstart.to_f >= time[0].to_f && tstart.to_f <= time[1].to_f
						timef = time[2].to_f
						break;
					end
				end
				@longitudeFromPeriod2 = timef.to_f
			rescue SystemCallError
				@longitudeFromPeriod2  = -999
			end
		end

		def latitudeFromPeriod2(tstart, tend)
			begin
				timef = 0
				File.open(ENV["AGILE"]+"/scripts/AGILEPOINTING").each_line do | line |
					time = line.split(" ");
					@latitudeFromPeriod2 = -999
					if tstart.to_f >= time[0].to_f && tstart.to_f <= time[1].to_f
						timef = time[3].to_f
						break;
					end
				end
				@latitudeFromPeriod2 = timef
			rescue SystemCallError
				@latitudeFromPeriod2  = -999
			end
		end

		def longitudeFromTime(time)
			datautils = DataUtils.new
			puts "time " + time.to_s
			datautils.extractContact(time)
			ec = datautils.ec
# 			puts ec
			@longitudeFromTime = longitude(ec);
			@longitudeFromTime
# 			puts @longitudeFromTime
		end

		def latitudeFromTime(time)
			datautils = DataUtils.new
			datautils.extractContact(time)
			ec = datautils.ec
			@latitudeFromTime = latitude(ec);
			@latitudeFromTime
		end

		#si ottiene il centro del puntamento (galactic longitude l)
		def longitude(contact)
			begin
				
				datautils = DataUtils.new
				if contact.to_i > 100000
					datautils.extractContact(contact)
					c = datautils.ec;
				else
					c = contact;
				end
				@longitude = 0;
				getRow(c);
				if @line == nil
					@longitude = -999;
					return
				end	
				l = @line.split(" ")[1];
				@longitude = l.to_f;
				@longitude ;
			rescue SystemCallError
				@longitude = -999;
			end
		end


		#si ottiene il centro del puntamento (galactic latitude b)
		def latitude(contact)
			begin
				datautils = DataUtils.new
				if contact.to_i > 100000
					datautils.extractContact(contact)
					c = datautils.ec;
				else
					c = contact;
				end
				@latitude = 0;
				getRow(c);
				if @line == nil
					@latitude = -999;
					return
				end	
				b = @line.split(" ")[3];
				@latitude = b.to_f;
				@latitude;
			rescue SystemCallError
				@latitude = -999;
			end
		end

		def pathbase=(pb)
			@pathbase=pb;
			@pathbase;
		end

		def line
			@line
		end
	
end
