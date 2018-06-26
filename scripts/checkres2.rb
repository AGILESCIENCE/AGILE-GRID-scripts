#! /usr/bin/ruby
#0) file name input
#1) file name output (optional): if specified, generate a multi file
#2) mints
#3) color
#4) 2AGL name generation (1 yes, 0 no)
#5) minimum exposure (default 0, e.g. 5.4e8)
#6) systematic error (default 0.2)
#7) list of sources to write to output

load ENV["AGILE"] + "/scripts/conf.rb"

filename = ARGV[0]

datautils = DataUtils.new

if ARGV[0].to_s == "help" || ARGV[0] == nil || ARGV[0] == "h"
	system("head -10 " + $0 );
	exit;
end

if ARGV[5] != nil
	expmin = ARGV[5]
else
	expmin = 0
end

if ARGV[6] != nil
	systematicerror = ARGV[6]
else
	systematicerror = 0.2
end


File.open(ARGV[0]).each do | line1 |

	File.open(ARGV[0]).each do | line2 |
		if line1 != line2
				lll1 = line1.split(" ")
				lll2 = line2.split(" ")
 			 	l1 = lll1[2]
			 	b1 = lll1[3]
			 	l2 = lll2[2]
			 	b2 = lll2[3]
				ring1_l = lll1[17].split(".")[0].to_f
				ring1_b = lll1[17].split(".")[1].to_f
				ring2_l = lll2[17].split(".")[0].to_f
				ring2_b = lll2[17].split(".")[1].to_f

				d = datautils.distance(l1, b1, l2, b2)
				d1 = datautils.distance(l1, b1, ring1_l, ring1_b)
				d2 = datautils.distance(l2, b2, ring2_l, ring2_b)

				if(d.to_f < 1)
					puts "\n" + line1.chomp + " " + d1.to_s + " - " +  line2.chomp + " " + d2.to_s + "\n"
				end
				if lll1[0] == lll2[0]
					puts "\nFound duplicated: \n" + line1.chomp + " " + d1.to_s + "\n" +  line2.chomp + " " + d2.to_s +  "\n"
				end
		end
	end
end


if ARGV[1] != nil
	if ARGV[2] != nil
		mints = ARGV[2]
	else
		mints = 2
	end
	if ARGV[3] != nil
		color = ARGV[3]
	else
		color = "green"
	end
	fout = File.new(ARGV[1] + ".multi", "w")
	fout2 = File.new(ARGV[1] + ".reg", "w");
	fout3 = File.new(ARGV[1] + ".ell", "w");
	fout4 = File.new(ARGV[1] + ".res2", "w");
	File.open(ARGV[0]).each do | line1 |
		lll1 = line1.split(" ")
		name = lll1[0]
		
		ts = lll1[1]
		flux = lll1[6]
		fixflag = "0"
		dist = "0.0"
		l1 = lll1[2]
		b1 = lll1[3]
		si1 = lll1[8]
		exp = lll1[16]
		if ts.to_f >= mints.to_f && exp.to_f > expmin.to_f
			if lll1[10].to_i != -1
				lll = lll1[10]
				bbb = lll1[11]
			else
				lll = l1
				bbb = b1
			end
			if(ARGV[4].to_i == 1)
				fidl = File.new("idlcom", "w")

                        	fidl.write("getsource_name, " + lll.to_s + ", " + bbb.to_s)

                        	fidl.close()

                        	system("idl < idlcom > /dev/null 2> /dev/null")

                        	ssname = ""

                        	File.open("source_name.prt", "r").each_line do |line2|

                            ssname = line2

                        	end

                        	ssname = "2AGLJ" + ssname.chomp

				name = ssname;
				#puts name
			end
			wr = false
			if ARGV[7] != nil
				File.open(ARGV[7]).each_line do | line2 |
					if name.chomp == line2.chomp
						wr = true
					end
				end
			else
				wr = true
			end
			if wr == false
				next
			end
			fout.write(flux + " " + l1 + " " + b1 + " " + si1.to_s + " 0 2 " + name + " 0.0\n")
			out2 = name
			for i in 1..lll1.size()
				out2 = out2 + " " + lll1[i].to_s
			end
			fout4.write(out2 + "\n")
			
			
			
			if lll1[10].to_i != -1
				fout2.write("galactic;ellipse(" + lll1[10].to_s + "," + lll1[11].to_s + "," + (lll1[13].to_f+systematicerror.to_f).to_s + "," + (lll1[14].to_f+systematicerror.to_f).to_s + ", " + (- lll1[15].to_f).to_s + ") # color=" + color + " text={" + name + "}")
				fout3.write(name + " " + lll1[10].to_s + " " + lll1[11].to_s + " " + (lll1[13].to_f+systematicerror.to_f).to_s + " " + (lll1[14].to_f+systematicerror.to_f).to_s + " " + (- lll1[15].to_f).to_s)
				
			else
				fout2.write("galactic;ellipse(" + l1.to_s + "," + b1.to_s + ",1800\",1800\",0.0) # color=" + color + " text={" + name + "}")
				fout3.write(name + " " + l1.to_s + " " + b1.to_s + " 0.6 0.6 0.0 ")
			end
			fout2.write("\n")
			fout3.write("\n")
		end
	end
	fout.close()
	fout2.close()
	fout3.close()
	fout4.close()
end

