
class Fits

		def readFitsHeader(filename)
			fexta = filename.split(".")
			fext = fexta[fexta.size() -1]
			out = filename
			if fext == "gz"
				value = ""; 20.times{value  << (65 + rand(25)).chr}
				out = "/tmp/" + value
				cmd = "gunzip -c " + filename + " > " + out
				#puts cmd
				system(cmd)
				
			else
				out = filename
			end
			#puts filename
			
			s = File.open(out, 'r') { |io| io.read };
			str = s[0..5000].to_str;
			str2 = str.split(" ")
			@header = Hash.new
			for i in 0..str2.size()
				if str2[i] != nil and str2[i].include?("=")
					@header[str2[i].split("=")[0]] = str2[i+1].gsub("'", " ").strip
				end
				if str2[i] == "="
					@header[str2[i-1]] = str2[i+1].gsub("'", " ").strip
				end
			end
			#puts @header
			if(out != filename)
				system("rm " + out)
			end
			
			#AGILE FITS
			@binsize = @header["CDELT2"].to_f;
			if @binsize.to_f < 0
				@binsize = -@binsize
			end
			

			naxis1 = @header["NAXIS1"].to_f;
			@mapsize = ((naxis1.to_i) * binsize.to_f).to_i
			
			@lcenter = @header["CRVAL1"].to_f;
			@bcenter = @header["CRVAL2"].to_f;
			
			@minenergy = @header["MINENG"].to_f;
			@maxenergy = @header["MAXENG"].to_f;
			
			@utc_start = @header["DATE-OBS"];
			
			@utc_end = @header["DATE-END"];
			
			@tt_start = @header["TSTART"];
			
			@tt_end = @header["TSTOP"];
			
		end	
		
		def header
			@header
		end	
		
		def binsize
			@binsize
		end	
		
		def mapsize
			@mapsize
		end
		
		def lcenter
			@lcenter
		end	
		
		def bcenter
			@bcenter
		end	
		
		def utc_start
			@utc_start
		end	
		
		def utc_end
			@utc_end
		end	
		
		def tt_start
			@tt_start
		end	
		
		def tt_end
			@tt_end
		end	
		
		def minenergy
			@minenergy
		end	
		
		def maxenergy
			@maxenergy
		end	
end