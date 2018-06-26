

class DataConversion
	def convertMultiToList(inputfilename, addfile)
		a = File.open(inputfilename)
		binsize = 0
		index = 0;
		lastline = "";
		c = File.new(inputfilename.to_s + ".list", "w");
		a.each_line do | l |
			if l.split("!").size() >= 2
				next
			end
			if l.split("'").size() >= 2
				next
			end
			if index.to_i == 0
				outline = "         L         B  Src cnts       Err  sqrt(TS)      Flux       Err\n"
				c.write(outline)
				outline = "		 deg       deg                              cm^-2 s^-1\n"
				c.write(outline)
			end
			if index.to_i >= 2
				ee = l.split(" ");
	
				a=-1; #semiassi .con
				b=-1;
				r = -1;
				ulv = 0;
				#legge i semiassi maggiore e minore dai file di output della multi2
				if File.exists?(inputfilename.to_s + "_" + ee[0].to_s + ".con") == true
					indexll = 0;
					File.open(inputfilename.to_s + "_" + ee[0].to_s).each_line do | ll2 |
						if ll2.split("!").size() >= 2
							next
						end
						if indexll.to_i == 3
							r = ll2.split(" ")[2]
							a = ll2.split(" ")[3]
							b = ll2.split(" ")[4]
						end
						if indexll.to_i == 5
							ulv = ll2.split(" ")[4]
						end
						indexll = indexll + 1;
					end
				else
					indexll = 0;
					File.open(inputfilename.to_s + "_" + ee[0].to_s).each_line do | ll2 |
						if ll2.split("!").size() >= 2
							next
						end
						if indexll.to_i == 4
							ulv = ll2.split(" ")[4]
						end
						indexll = indexll + 1;
					end
				end
				outline = "\t" + ee[2].to_s + "\t " + (ee[3].to_f + binsize.to_f).to_s + "\t " + ee[4].to_s + "\t " + ee[5].to_s + "\t " + ee[1].to_s + "\t " + ee[6].to_s + "\t " + ee[7].to_s;
	# 			puts outline
				c.write(outline);
				if addfile != ""
					if inputfilename.split("res")[1] != nil
						c.write("\t" + addfile.to_s + inputfilename.split("res")[1] + "_" + ee[0].to_s + " " + r.to_s + " " + a.to_s + " " + b.to_s + " " + ulv.to_s);
					else
						c.write("\t" + addfile.to_s  + "_00_" + ee[0].to_s + " " + r.to_s + " " + a.to_s + " " + b.to_s + " " + ulv.to_s);
					end
				end
				c.write("\n");
			end
	
			index = index + 1;
		end
		c.close();
	end
end