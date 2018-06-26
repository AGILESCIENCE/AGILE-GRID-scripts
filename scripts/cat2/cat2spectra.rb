#! /usr/bin/ruby
sourcename = ARGV[0]
spectratype = ARGV[1] #pl plec plsec lp
energyrange = ARGV[2] #00100-10000 00100-50000
analysisname = ARGV[3] #EDP1-EB01-FB01
irf = ARGV[4]
inttype = ARGV[5] #1-8
minradius = ARGV[6]

fixsi = nil
if ARGV[7] != nil
	fixsi = ARGV[7].to_s;
end

fan = "R" + inttype.to_s  + "_C" + format("%02d", minradius.to_f*10) + "-" + ARGV[1] + "-" + ARGV[2] + "-" + ARGV[3] + "-" + ARGV[4]

maplist4name = fan + "_FM3.119_ASDCe_"+irf+"_B01_"+energyrange+".maplist4"

f1 = File.new(maplist4name, "w")
gcf = ""
if irf == "H0025"
	f1.write("EMIN00100_EMAX00300_FM3.119_ASDCe_H0025_B01.cts.gz EMIN00100_EMAX00300_FM3.119_ASDCe_H0025_B01.exp.gz EMIN00100_EMAX00300_FM3.119_ASDCe_H0025_B01.gas.gz 25 -1 -1\n")
	f1.write("EMIN00300_EMAX01000_FM3.119_ASDCe_H0025_B01.cts.gz EMIN00300_EMAX01000_FM3.119_ASDCe_H0025_B01.exp.gz EMIN00300_EMAX01000_FM3.119_ASDCe_H0025_B01.gas.gz 25 -1 -1\n")
	f1.write("EMIN01000_EMAX03000_FM3.119_ASDCe_H0025_B01.cts.gz EMIN01000_EMAX03000_FM3.119_ASDCe_H0025_B01.exp.gz EMIN01000_EMAX03000_FM3.119_ASDCe_H0025_B01.gas.gz 25 -1 -1\n")
	f1.write("EMIN03000_EMAX10000_FM3.119_ASDCe_H0025_B01.cts.gz EMIN03000_EMAX10000_FM3.119_ASDCe_H0025_B01.exp.gz EMIN03000_EMAX10000_FM3.119_ASDCe_H0025_B01.gas.gz 25 -1 -1\n")
	gcf = "0.7,0.7,0.7,0.7"
	if energyrange.split("-")[1].to_i == 50000
		f1.write("EMIN10000_EMAX50000_FM3.119_ASDCe_H0025_B01.cts.gz EMIN10000_EMAX50000_FM3.119_ASDCe_H0025_B01.exp.gz EMIN10000_EMAX50000_FM3.119_ASDCe_H0025_B01.gas.gz 25 -1 -1\n")
		gcf = "0.7,0.7,0.7,0.7,0.7"
	end
	
end

if irf == "I0025"
	f1.write("EMIN00100_EMAX00200_FM3.119_ASDCe_I0025_B01.cts.gz EMIN00100_EMAX00200_FM3.119_ASDCe_I0025_B01.exp.gz EMIN00100_EMAX00200_FM3.119_ASDCe_I0025_B01.gas.gz 25 -1 -1\n")
	f1.write("EMIN00200_EMAX00400_FM3.119_ASDCe_I0025_B01.cts.gz EMIN00200_EMAX00400_FM3.119_ASDCe_I0025_B01.exp.gz EMIN00200_EMAX00400_FM3.119_ASDCe_I0025_B01.gas.gz 25 -1 -1\n")
	f1.write("EMIN00400_EMAX01000_FM3.119_ASDCe_I0025_B01.cts.gz EMIN00400_EMAX01000_FM3.119_ASDCe_I0025_B01.exp.gz EMIN00400_EMAX01000_FM3.119_ASDCe_I0025_B01.gas.gz 25 -1 -1\n")
	f1.write("EMIN01000_EMAX03000_FM3.119_ASDCe_I0025_B01.cts.gz EMIN01000_EMAX03000_FM3.119_ASDCe_I0025_B01.exp.gz EMIN01000_EMAX03000_FM3.119_ASDCe_I0025_B01.gas.gz 25 -1 -1\n")
	f1.write("EMIN03000_EMAX10000_FM3.119_ASDCe_I0025_B01.cts.gz EMIN03000_EMAX10000_FM3.119_ASDCe_I0025_B01.exp.gz EMIN03000_EMAX10000_FM3.119_ASDCe_I0025_B01.gas.gz 25 -1 -1\n")
	gcf = "0.7,0.7,0.7,0.7"
	if energyrange.split("-")[1].to_i == 50000
		f1.write("EMIN10000_EMAX50000_FM3.119_ASDCe_I0025_B01.cts.gz EMIN10000_EMAX50000_FM3.119_ASDCe_I0025_B01.exp.gz EMIN10000_EMAX50000_FM3.119_ASDCe_I0025_B01.gas.gz 25 -1 -1\n")
		gcf = "0.7,0.7,0.7,0.7,0.7"
	end
end

f1.close()

#detcatline
catline = " "
galcoeff = "-1"
#192 314
File.open("/ANALYSIS3/catalogs/cat2.multi").each do | line |
	ll = line.split(" ")
	if ll[6] == sourcename
		catline = ll[0] + " " + ll[1] + " " + ll[2] + " " + ll[3]
		# " 2.1 "
		#+ ll[3]
		
		if ll[2].to_f < -10 or ll[2].to_f > 10
			galcoeff = gcf
		end
		
		if spectratype == "pl"
			fixflag = "4" #4
			endline = "0 0.0 0.0"
		end
		if spectratype == "plec"
			fixflag = "12"
			endline = "1 2000.0 0.0"
		end
		if spectratype == "plsec"
			fixflag = "28"
			endline = "2 2000.0 1.0"
		end
		if spectratype == "lp"
			fixflag = "28"
			endline = "3 2000.0 1.0"
		end
		
		if fixsi != nil
			catline = ll[0] + " " + ll[1] + " " + ll[2]  + " " + fixsi.to_s + " 1 " + ll[5] + " " + ll[6] + " " + ll[7] + " 0 0.0 0.0"
		else
			catline = catline + " " + fixflag.to_s + " " + ll[5] + " " + ll[6] + " " + ll[7]
			catline += " "
			catline += endline
		end
		
		if ll.size > 11
			for i in 11..ll.size-1
				catline += " "
				catline += ll[i].to_s
			end
		end
		
		
		
		
		break
	end
end



system("rm INT_"+fan+"*")
cmd = "multi6.rb FM3.119_ASDC2_"+irf+" " + maplist4name +" none INT_"+fan+" addcat=\""+ catline +"\" catminradius=" + minradius.to_s + " catminflux=25e-08 fluxcorrection=1 scanmaplist=" + sourcename + "," + fan + " minimizertype=Minuit minimizeralg=Migrad minimizerdefstrategy=2 fluxcorrection=1 galmode2=3 isomode2=3 galmode2fit=0 isomode2fit=0 integratortype=" + inttype.to_s + " galcoeff=" + galcoeff.to_s
#galmode2=3 isomode2=3
puts cmd
system cmd
#fluxcorrection=1 isomode2=1 isomode2fit=2  # galmode2=1 isomode2=1 galmode2fit=1 isomode2fit=2

