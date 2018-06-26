#! /usr/bin/ruby
#0) start file of source list
#1) maplist
#2) filter
#3) dir output
#4) fixflag of analysis of the main source
#5) distanceToFixFlag0: see next parameters
#6) fixflagneighbour: if the source of the list is < distanceToFixFlag0 put its fixflag=fixflagneighbour
#7) additional commands to multi5.rb (optional)
#8) fixisogalstep0 = 0 none, 1 apply
#9) update results for each step (0, 1 - default 1)
#10) galcoeff (apply this galcoeff for source with | b | > 7
#11) maxdistancefromcenter: default 360. Use fixflag=1 if the source position is maxdistancefromcenter of the map

#Questo script si usa per fare la scansione su una lista di sorgenti (source list) fissandole tutte tranne una. Alla fine si raccoglie il risultato finale nelle directory che viene creata (dir output)
#NB: tutte le sorgenti sono messe con fixflag=0 di default prima di inziare l'analisi
#Le sorgenti che hanno il nome che inizia con _ non sono analizzate, ma lasciate nel fondo
#Per le sorgenti che hanno il nome che inizia con # si mette fixflag=3

load ENV["AGILE"] + "/scripts/conf.rb"

def savesourcelist(listfile, sources2)
	fo = File.new(listfile, "w")
	#scrivi il file delle sorgenti
	sources2.each { |s2|
		out = s2.output
		s2.print
		fo.write(out + "\n")
	}
	fo.close()
end

class Source < 
	Struct.new(:flux, :l, :b, :si, :fixflag, :minsqrtts, :name, :rmax)
	 	
	def print
 		puts flux.to_s + " " + l.to_s + " " + b.to_s + " " + si.to_s + " " + fixflag.to_s + " " + minsqrtts.to_s + " " + name.to_s + " " + rmax.to_s
	end

	def output
		@output = flux.to_s + " " + l.to_s + " " + b.to_s + " " + si.to_s + " " + fixflag.to_s + " " + minsqrtts.to_s + " " + name.to_s + " " + rmax.to_s  
	end
end

begin 
	if ARGV[0].to_s == "help" || ARGV[0] == nil || ARGV[0] == "h"
		system("head -19 " + $0 );
		exit;
	end

	datautils = DataUtils.new

	sourcelist = ARGV[0]
	maplist = ARGV[1]
	puts "maplist="+maplist
	filter = ARGV[2]
	puts "filter=" + filter
	diroutput = ARGV[3]
	multilist = diroutput
	fixflaganalysis=ARGV[4]

	distanceToFixFlag0 = ARGV[5]
	fixflagneighbour = ARGV[6]

	additionalcmd=nil
	if ARGV[7] != nil
		additionalcmd=ARGV[7]
	end

	fixisogalstep0 = 0
	if ARGV[8] != nil
		fixisogalstep0  = ARGV[8]
	end

	updateres = 1
	if ARGV[9] != nil
		updateres  = ARGV[9]
	end

	galcoeff = "-1"
	if ARGV[10] != nil
		galcoeff  = ARGV[10]
	end
	
	maxdistancefromcenter = 360.0
	if ARGV[11] != nil
		maxdistancefromcenter  = ARGV[11]
	end
	
	cts=""
	File.open(maplist).each_line do | line |
		cts = line.split(" ")[0]
	end
	fits = Fits.new
	fits.readFitsHeader(cts);

	outlog = diroutput + ".log"

	cmd = "mkdir " + diroutput
	datautils.execute(outlog, cmd)

	resfilename = diroutput.to_s + ".res"

	sources = Array.new

	a = File.open(sourcelist, "r")

	a.each_line do | line |
		 words = line.split(" ")
		 s = Source.new
		 s.flux = words[0].to_f
		 s.l = words[1].to_f
		 s.b = words[2].to_f
		 s.si = words[3].to_f
		 s.fixflag = words[4].to_i
		 s.minsqrtts = words[5].to_f
		 s.name = words[6].to_s
		 s.rmax = words[7].to_f
		 sources.push(s)
		#   s.print
	end

	sources2 = sources.sort_by { |a| [ -a.flux ] }

	sources2.each { |s|
		s.print
		s.fixflag=0
		#s.rmax = 2.0
	}

	index = 0
	ffinal = diroutput.to_s + ".resfinal"
	ffinalfull = diroutput.to_s + ".resfinalfull"
	ffmulti = diroutput.to_s + ".res.multi"

	fout1 = File.new(ffinal, "w")
	fout2 = File.new(ffinalfull, "w")
	fout3 = File.new(ffmulti, "w")

	sources2.each { |s|

		namesource = s.name
		#salta quelle che iniziano con _
		if namesource[0] == "_" or namesource[0] == 95
			next
		end
		puts namesource
	
		if distanceToFixFlag0.to_f > 0
			sources2.each { |s1|
				d = datautils.distance(s1.l, s1.b, s.l, s.b)
				if d.to_f > distanceToFixFlag0.to_f
					s1.fixflag=0
				else
					s1.fixflag=fixflagneighbour
				end
			}
		end
	
		#metti la sorgente da analizzare con fixflag passato da input
		s.fixflag = fixflaganalysis
		#se le sorgenti iniziano con # mett fixflag=3
		if (namesource[0] == "#" or namesource[0] == 35) and s.fixflag.to_i == 7
			s.fixflag = 3;
		end
	
		#prepare additional commands for multi5.rb
		addcmd = ""
		if additionalcmd != nil
			addcmd = additionalcmd
		end
	
		if fixisogalstep0.to_i == 1
			addcmd = addcmd + " fixisogalstep0=" + namesource.to_s + " "
		end
	
		if s.b >= 6 or s.b <= -6
			addcmd += " galcoeff=";
			addcmd += galcoeff.to_s;
			addcmd += " "
		end

		#----------
		resfilename = diroutput.to_s + "_" + format("%03d", index) + ".res"
		listfile = multilist.to_s + "_" + format("%03d", index) + ".multi"
		savesourcelist(listfile, sources2)
		cmd = "ruby " + ENV["AGILE"] + "/scripts/multi5.rb " + filter + " " + maplist.to_s + " " + listfile.to_s + " " + resfilename.to_s + " " + addcmd.to_s
		datautils.execute(outlog, cmd)
		#--------
	
		#modifico ora i dati della sorgente che e' stata calcolata
		sout = MultiOutput.new
		sout.readDataSingleSource(resfilename.to_s + "_" + namesource.to_s);
	
		#se la significativita' con fixflag=7 e' troppo bassa, metti fixflag=3 e ripeti l'analisi
		#e crea i file che finiscono con .ff3.res
		postfix=""
		if sout.sqrtTS.to_f < 4.0 and s.fixflag.to_i == 7
	
			s.fixflag = 3;
		
			resfilename = diroutput.to_s + "_" + format("%03d", index) + ".ff3.res"
			listfile = multilist.to_s + "_" + format("%03d", index) + ".ff3.multi"
			postfix += "_ff3"
			savesourcelist(listfile, sources2)
			cmd = "ruby " + ENV["AGILE"] + "/scripts/multi5.rb " + filter + " " + maplist.to_s + " " + listfile.to_s + " " + resfilename.to_s + " " + addcmd.to_s
			datautils.execute(outlog, cmd)
		
			sout.readDataSingleSource(resfilename.to_s + "_" + namesource.to_s);
	
		end	
		
		#verifica la distanza della sorgente calcolata
		d1 =  datautils.distance(sout.l_peak, sout.b_peak, fits.lcenter, fits.bcenter)
		if d1.to_f > maxdistancefromcenter.to_f
			s.fixflag = 1;
		
			resfilename = diroutput.to_s + "_" + format("%03d", index) + ".ff1.res"
			listfile = multilist.to_s + "_" + format("%03d", index) + ".ff1.multi"
			postfix += "_ff1"
			savesourcelist(listfile, sources2)
			cmd = "ruby " + ENV["AGILE"] + "/scripts/multi5.rb " + filter + " " + maplist.to_s + " " + listfile.to_s + " " + resfilename.to_s + " " + addcmd.to_s
			datautils.execute(outlog, cmd)
		
			sout.readDataSingleSource(resfilename.to_s + "_" + namesource.to_s);
		end
	
		d = datautils.distance(sout.l_peak, sout.b_peak, s.l, s.b)
	
		#eseguita la alike, copio i risulati
		cmd = "cp " + resfilename.to_s + "_" + namesource.to_s + "* " + diroutput
		datautils.execute(outlog, cmd)
		cmd = "cp " + resfilename.to_s + " " + diroutput + "/" + namesource.to_s + ".res"
		datautils.execute(outlog, cmd)

		fout1.write(sout.multiOutputLine + "\n")
		fout2.write(sout.multiOutputLineFull3(diroutput + postfix) + " " + d.to_s + "\n")

		#aggiorna i valori
		if updateres.to_i == 1
			s.flux = sout.flux
			s.l = sout.l_peak
			s.b = sout.b_peak
			s.si = sout.sicalc
		end
		#alla fine rimettilo a 0
		sources2.each { |s|
			s.print
			s.fixflag=0
		}

		fout3.write(sout.flux.to_s + " " + sout.l_peak.to_s + " " + sout.b_peak.to_s + " " + sout.sicalc.to_s + " " + " 0 2 " + namesource.to_s + "\n")
		
		index = index + 1
	}
	puts index
	fout1.close()
	fout2.close()
	fout3.close()

	cmd = "cp " + ffinal.to_s + " " + diroutput.to_s + "/" + diroutput.to_s + ".res"
	datautils.execute(outlog, cmd)

	cmd = "cp " + ffinalfull.to_s + " " + diroutput.to_s
	datautils.execute(outlog, cmd)

	cmd = "cp " + ffmulti.to_s + " " + diroutput.to_s
	datautils.execute(outlog, cmd)

end
