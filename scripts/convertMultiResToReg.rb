#! /usr/bin/ruby
#0) input file name
#1) color (optiona, red)
#2) systematic error
#3) min sqrt(TS) 

load ENV["AGILE"] + "/scripts/conf.rb"

if ARGV[0].to_s == "help" || ARGV[0] == nil || ARGV[0] == "h"
	system("head -4 " + $0 );
	exit;
end

fileinp1 = ARGV[0]
a = File.open(fileinp1)

if ARGV[1] != nil
	color = ARGV[1].to_s;
else
	color = "white";
end
if ARGV[2] != nil
	systematic = ARGV[2].to_s;
else
	systematic = 0.0;
end

if ARGV[3] != nil
	minsqrtTS = ARGV[3].to_s;
else
	minsqrtTS = 0.0;
end

	multioutput = MultiOutput.new
	
	index = 0;
	lastline = "";
	c = File.new(ARGV[0] + ".reg", "w");
	c1 = File.new(ARGV[0] + ".tex", "w");
	#c2 = File.new(ARGV[0] + ".multi", "w");
	c3 = File.new(ARGV[0] + ".resconv", "w");
	modefile = 3
	a.each_line do | l |
		if l.split("!").size() >= 2
			next
		end
		if l.split("'").size() >= 2
			next
		end
		fs = l.split("Galactic")
		if fs.size == 2
			modefile = 4
			next
		end
		fs = l.split("Isotropic")
		if fs.size == 2
			modefile = 4
			next
		end
		fs = l.split(" ")
		if fs.size == 4
			next
		end
		if index.to_i == 0
			outline = "global color=green font=\"helvetica 9 normal\" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n"
			c.write(outline)
			outline = "galactic\n"
			c.write(outline)
			c1.write("\\begin {table}[!htb]\n")
			c1.write("\\caption {\\em{}}\n")
			c1.write("\\label{table_}\n")
			c1.write("\\renewcommand{\\arraystretch}{1.2}\n");
			c1.write("\\begin{tabular}{@{}llllllll} \\hline\n");
			c1.write("\\textbf{ Name } & \\textbf{ l } & \\textbf{ b } & Counts & \\textbf{ $\\sqrt(TS)$ } & \\textbf{ Flux 10$^{-6}$ }  & \\textbf{ r } & \\textbf{ S.Index } \\\\\n")
		end
		if index.to_i >= 0
			puts l
			ee = l.split(" ");
			if ee.size() <= 2
				next
			end
			
			#new
			multioutput.readDataSingleSource2(fileinp1, ee[0].to_s)
			dist = multioutput.sicalc
			disterr = multioutput.sicalc_error
			#dist = ""	
			#disterr = ""	
			#if ee[8] != nil
			#	dist = ee[8] #distanza o spectral index in multi4
			#end
			#if ee[9] != nil
			#	disterr = ee[9]
			#end
			if ee[1] != "nan" && ee[6].to_s != "inf" && ee[6] != "nan"
				if ee[4].to_f > 0 && ee[1].to_f > minsqrtTS.to_f
					#recupera il raggio
					nameout1 = fileinp1.to_s + "_" + ee[0].to_s
					nameout2 = nameout1 + ".con"
					rrr = "500\""
					lll = multioutput.l_peak
					bbb = multioutput.b_peak
					marker = "";
					if File.exists?(nameout2) == true
						#indexll1 = 0
						#iidr = 9
						#File.open(nameout1).each_line do | ll1 |
						#	if modefile.to_i == 4
						#		iidr = 10
						#	end
						#	if indexll1.to_i == iidr.to_i
						#		ll11 = ll1.split(" ");
						#		rrr = ll11[2].to_s
								#puts rrr
						#		break;
						#	else
						#		indexll1 = indexll1 + 1
						#	end
						#end
						#if rrr.to_s == "nan"
						#	rrr = "500\""
						#end
						rrr = multioutput.r.to_f + systematic.to_f
						lll = multioutput.l
						bbb = multioutput.b
					else
						rrr = "500\""
						marker = "*"
					end
					
					#fidl = File.new("idlcom", "w")
					#fidl.write("getsource_name, " + lll.to_s + ", " + bbb.to_s)
					#fidl.close()
					#system("$IDL_DIR/bin/idl < idlcom")
					#ssname = ""
					#File.open("source_name.prt", "r").each_line do |line2|
					#	ssname = line2
					#end
					#ssname = "2AGLJ" + ssname.chomp
					ssname = ee[0].to_s
					#outline = "ellipse(" + ee[2].to_s + "," + ee[3].to_s + "," + rrr.to_s + "," + rrr.to_s + ",0) #color=" + color.to_s + " width=1 text={" + ee[0].to_s + " " + format("%.2f", ee[1])  + marker.to_s + " (" + format("%.2f", ee[2]) + "," + format("%.2f", ee[3]) + "," + dist.to_s + ")}";
					outline = "ellipse(" + lll.to_s + "," + bbb.to_s + "," + rrr.to_s + "," + rrr.to_s + ",0) #color=" + color.to_s + " width=1 text={" + ssname + " " + format("%.2f", ee[1])  + marker.to_s + " (" + format("%.2f", lll) + "," + format("%.2f", bbb) + "," + dist.to_s + ")}";
					c.write(outline);
					c.write("\n");
					mmou = 0
					if multioutput.r.to_f == -1
						mmou = -1
					else
						mmou = rrr
					end
					c1.write(ssname.to_s + " & " + format("%.2f", ee[2]) + " & " + format("%.2f", ee[3]) + " & " + format("%.2f", ee[4]) + " $\\pm$ " + format("%.2f", ee[5]) + " & " + format("%.2f", ee[1]) + " & " + format("%.2f",ee[6].to_f * 10e7) + " $\\pm$ " + format("%.2f",ee[7].to_f * 10e7) +  " & " + format("%.2f",mmou) +  " & " + format("%.2f", ee[8]) + " $\\pm$ " + format("%.2f",ee[9]) + "\\\\" + "\n");
					#c2.write( ee[6].to_s + " " + lll.to_s + " " + bbb.to_s + " 2.1 0 2 " + ssname.to_s + "\n") 
					#puts  ee[0].to_s + " & " + ee[2].to_s + " & " + ee[3].to_s + " & " + ee[4].to_s + " $\\pm$ " + ee[5].to_s + " & " + ee[1].to_s + " & " + ee[6].to_s + " $\\pm$ " + ee[7].to_s +  " & " + dist.to_s + " $\\pm$ " + disterr.to_s + "\\\\";
					
					eeindex = 0
					outress = ""
					ee.each do | e |
						
						
						if eeindex.to_i == 0
							outress = ssname
							eeindex = eeindex + 1
							next
						end
						if eeindex.to_i != 12 || eeindex.to_i != 13 || eeindex.to_i != 14
							outress = outress.to_s + " " + e.to_s
						else
							outress = outress.to_s + " " + (e.to_f+systematic.to_f).to_s
						end
						eeindex = eeindex + 1
					end
					c3.write(outress + "\n")
					
				end
			end
		end

		index = index + 1;
	end
	c.close();
	c1.write(  "\\hline \\end{tabular} \\end{table}\n" )
	c1.close()
	#c2.close()
	c3.close()

