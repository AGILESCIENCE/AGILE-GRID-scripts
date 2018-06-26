load ENV["AGILE"] + "/scripts/conf.rb"

prefix=ARGV[0]

fout = File.open(ARGV[1], "w")

util = DataUtils.new

Dir["[0-9]*"].each do | dir |
        begin
        	if File.exists?(dir + "/" + prefix + "/" + prefix + ".resfinalfull")	
				File.open(dir + "/" + prefix + "/" + prefix + ".resfinalfull").each_line do | line |
					#puts line
					ll = line.split(" ")
					lpeak = ll[4].to_f
					bpeak = ll[5].to_f
					l = ll[7].to_f
					b = ll[8].to_f
					ring = dir
					lring = dir.split("_")[0].to_f
					bring = dir.split("_")[1].to_f
					dpeakring = util.distance(lpeak, bpeak, lring, bring);
					dellring = util.distance(l, b, lring, bring);
					dellpeak = util.distance(lpeak, bpeak, l, b);
					#puts dellpeak
					fout.write(line.chomp + " " + dir + " " + format("%.6f %.6f DIST %.4f %.4f %.4f\n",lring.to_f, bring.to_f, dpeakring.to_f, dellring.to_f, dellpeak.to_f))
				end
		else
			puts "File not found of dir: " + dir
		end

        #rescue
                #puts "File not found of dir: " + dir
        end
end
fout.close()
