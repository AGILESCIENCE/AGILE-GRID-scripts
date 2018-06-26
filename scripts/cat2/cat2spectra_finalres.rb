#! /usr/bin/ruby
#0) cat sourcename
#1) spectral shape: pl plec plsec lp
#2) energyrange where to perform optimisation parameters: 00100-10000 00100-50000 00100-01000
#3) analysisname, e.g. EDP1-EB01-FB01
#4) IRF name, e.g. H0025
#5) integrator type: 1..8
#6) selection from cat multi: minradius around cat sourcename
#7) prefix (to be added to analysis name) e.g. FINAL
#8) add fix flag: 1 (only flux free) or 3 (flux and position free)
#9) fix spectral index (optional, set the value if >=0 or do not fix if < 0, e.g. -1)
#10) fix galactic coeff (optional, set the value if >=0 or keep free if < 0, e.g. -1)
#11) fix galiso step2 (optional, fix value if >=0 or keep free if < 0, e.g. -1)
#12) selection from cat multi: mincatflux (optional or e.g 25e-08)

#Get all sources for ff=1
# R1-G10 2AGL0161 localhost cat2 0 pl FINAL 1 2.1 -1 1 0e-08
#Standard set
# R1-G10 2AGL0161 localhost cat2 0 pl FINAL 1
# R1-G10 2AGL0161 localhost cat2 0 pl FINAL 3
# R1-G10 2AGL0161 localhost cat2 0 pl FINALgcf 3 -1 -1
# R1-G10 2AGL0161 localhost cat2 0 pl FINALgif 3 -1 -1 -1

sourcename = ARGV[0]
spectratype = ARGV[1] #pl plec plsec lp
energyrange = ARGV[2] #00100-10000 00100-50000 00100-01000
analysisname = ARGV[3] #EDP1-EB01-FB01
irf = ARGV[4]
inttype = ARGV[5] #1-8
minradius = ARGV[6]
prefix = ARGV[7]
addff = ARGV[8] #1 (only flux free) or 3 (flux and position free)

fixsi = nil
if ARGV[9] != nil
	fixsi = ARGV[9].to_s;
end

fixgal = 1;
fixgalcoeff = "0.7"
if ARGV[10] != nil
	fixgalcoeff = ARGV[10].to_s;
	if fixgalcoeff.to_f >= 0
		fixgal = 1
	else
		fixgal = 0
	end
end
puts fixgalcoeff

fixgalisostep2 = 1;
if ARGV[11] != nil
	fixgalisostep2 = ARGV[11].to_s
end

mincatflux = nil
if ARGV[12] != nil
	mincatflux = ARGV[12].to_s;
end

load ENV["AGILE"] + "/scripts/conf.rb"

alikeutils = AlikeUtils.new

##############################
# DETCATLINE
##############################

catline = " "
galcoeff = "-1"
galcoefffull = "-1"
galcoeffhe = "-1"
fixflag = 1
coordb = 0
#192 314
File.open("/ANALYSIS3/catalogs/cat2.multi").each do | line |
	ll = line.split(" ")
	if ll[6] == sourcename
		#catline = ll[0] + " " + ll[1] + " " + ll[2] + " 2.1 " # + ll[3]
		catline = ll[0] + " " + ll[1] + " " + ll[2] + " " + ll[3]
		coordb = ll[2].to_f
		
		if spectratype == "pl"
			fixflag = 4 #4
			if fixsi != nil and fixsi.to_i > 0
				fixflag = 0
			end
			endline = "0 0.0 0.0"
		end
		if spectratype == "plec"
			fixflag = 12
			endline = "1 2000.0 0.0"
		end
		if spectratype == "plsec"
			fixflag = 28
			endline = "2 2000.0 1.0"
		end
		if spectratype == "lp"
			fixflag = 28
			endline = "3 2000.0 1.0"
		end
		
		fixflag = fixflag.to_i +  addff.to_i
		
		if fixsi != nil and fixsi.to_i > 0
			catline = ll[0] + " " + ll[1] + " " + ll[2]  + " " + fixsi.to_s + " " + fixflag.to_s + " " + ll[5] + " " + ll[6] + " " + ll[7] + " 0 0.0 0.0"
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

##############################
# MAPLIST
##############################

suffix = "R" + inttype.to_s  + "_C" + format("%02d", minradius.to_f*10) + "-" + ARGV[1] + "-" + ARGV[2] + "-" + ARGV[3] + "-" + ARGV[4]

if fixflag.to_i == 1
	fan = prefix + "_FF2_" +  suffix
else
	fan = prefix + "_FF2"+fixflag.to_s+"_" +  suffix
end


maplist4name = fan + "_FM3.119_ASDCe_"+irf+"_B01_"+energyrange+".maplist4"

maplist4namefull = fan + "_FM3.119_ASDCe_"+irf+"_B01_"+energyrange+".full.maplist4"

maplist4namehe = fan + "_FM3.119_ASDCe_"+irf+"_B01_"+energyrange+".he.maplist4"

f1 = File.new(maplist4name, "w")
gcf = ""
if irf == "H0025"
	f1.write("EMIN00100_EMAX00300_FM3.119_ASDCe_H0025_B01.cts.gz EMIN00100_EMAX00300_FM3.119_ASDCe_H0025_B01.exp.gz EMIN00100_EMAX00300_FM3.119_ASDCe_H0025_B01.gas.gz 25 -1 -1\n")
	f1.write("EMIN00300_EMAX01000_FM3.119_ASDCe_H0025_B01.cts.gz EMIN00300_EMAX01000_FM3.119_ASDCe_H0025_B01.exp.gz EMIN00300_EMAX01000_FM3.119_ASDCe_H0025_B01.gas.gz 25 -1 -1\n")
	if energyrange.split("-")[1].to_i > 1000
		f1.write("EMIN01000_EMAX03000_FM3.119_ASDCe_H0025_B01.cts.gz EMIN01000_EMAX03000_FM3.119_ASDCe_H0025_B01.exp.gz EMIN01000_EMAX03000_FM3.119_ASDCe_H0025_B01.gas.gz 25 -1 -1\n")
		f1.write("EMIN03000_EMAX10000_FM3.119_ASDCe_H0025_B01.cts.gz EMIN03000_EMAX10000_FM3.119_ASDCe_H0025_B01.exp.gz EMIN03000_EMAX10000_FM3.119_ASDCe_H0025_B01.gas.gz 25 -1 -1\n")
		gcf = fixgalcoeff.to_s + "," + fixgalcoeff.to_s + "," + fixgalcoeff.to_s + "," + fixgalcoeff.to_s
		if energyrange.split("-")[1].to_i == 50000
			f1.write("EMIN10000_EMAX50000_FM3.119_ASDCe_H0025_B01.cts.gz EMIN10000_EMAX50000_FM3.119_ASDCe_H0025_B01.exp.gz EMIN10000_EMAX50000_FM3.119_ASDCe_H0025_B01.gas.gz 25 -1 -1\n")
			gcf = fixgalcoeff.to_s + "," + fixgalcoeff.to_s + "," + fixgalcoeff.to_s + "," + fixgalcoeff.to_s + "," + fixgalcoeff.to_s
		end
	else
		gcf = fixgalcoeff.to_s + "," + fixgalcoeff.to_s
	end
end
f1.close()

f2 = File.new(maplist4namefull, "w")
gcffull = ""
if irf == "H0025"
	f2.write("EMIN00030_EMAX00050_FM3.119_ASDCe_H0025_B01.cts.gz EMIN00030_EMAX00050_FM3.119_ASDCe_H0025_B01.exp.gz EMIN00030_EMAX00050_FM3.119_ASDCe_H0025_B01.gas.gz 25 -1 -1\n")
	f2.write("EMIN00050_EMAX00100_FM3.119_ASDCe_H0025_B01.cts.gz EMIN00050_EMAX00100_FM3.119_ASDCe_H0025_B01.exp.gz EMIN00050_EMAX00100_FM3.119_ASDCe_H0025_B01.gas.gz 25 -1 -1\n")
	f2.write("EMIN00100_EMAX00300_FM3.119_ASDCe_H0025_B01.cts.gz EMIN00100_EMAX00300_FM3.119_ASDCe_H0025_B01.exp.gz EMIN00100_EMAX00300_FM3.119_ASDCe_H0025_B01.gas.gz 25 -1 -1\n")
	f2.write("EMIN00300_EMAX01000_FM3.119_ASDCe_H0025_B01.cts.gz EMIN00300_EMAX01000_FM3.119_ASDCe_H0025_B01.exp.gz EMIN00300_EMAX01000_FM3.119_ASDCe_H0025_B01.gas.gz 25 -1 -1\n")
	f2.write("EMIN01000_EMAX03000_FM3.119_ASDCe_H0025_B01.cts.gz EMIN01000_EMAX03000_FM3.119_ASDCe_H0025_B01.exp.gz EMIN01000_EMAX03000_FM3.119_ASDCe_H0025_B01.gas.gz 25 -1 -1\n")
	f2.write("EMIN03000_EMAX10000_FM3.119_ASDCe_H0025_B01.cts.gz EMIN03000_EMAX10000_FM3.119_ASDCe_H0025_B01.exp.gz EMIN03000_EMAX10000_FM3.119_ASDCe_H0025_B01.gas.gz 25 -1 -1\n")
	f2.write("EMIN10000_EMAX50000_FM3.119_ASDCe_H0025_B01.cts.gz EMIN10000_EMAX50000_FM3.119_ASDCe_H0025_B01.exp.gz EMIN10000_EMAX50000_FM3.119_ASDCe_H0025_B01.gas.gz 25 -1 -1\n")
	gcffull = fixgalcoeff.to_s + "," + fixgalcoeff.to_s + "," + fixgalcoeff.to_s + "," + fixgalcoeff.to_s + "," + fixgalcoeff.to_s + "," + fixgalcoeff.to_s + "," + fixgalcoeff.to_s
end
f2.close()

f3 = File.new(maplist4namehe, "w")
gcfhe = ""
if irf == "H0025"
	f3.write("EMIN00100_EMAX00300_FM3.119_ASDCe_H0025_B01.cts.gz EMIN00100_EMAX00300_FM3.119_ASDCe_H0025_B01.exp.gz EMIN00100_EMAX00300_FM3.119_ASDCe_H0025_B01.gas.gz 25 -1 -1\n")
	f3.write("EMIN00300_EMAX01000_FM3.119_ASDCe_H0025_B01.cts.gz EMIN00300_EMAX01000_FM3.119_ASDCe_H0025_B01.exp.gz EMIN00300_EMAX01000_FM3.119_ASDCe_H0025_B01.gas.gz 25 -1 -1\n")
	f3.write("EMIN01000_EMAX03000_FM3.119_ASDCe_H0025_B01.cts.gz EMIN01000_EMAX03000_FM3.119_ASDCe_H0025_B01.exp.gz EMIN01000_EMAX03000_FM3.119_ASDCe_H0025_B01.gas.gz 25 -1 -1\n")
	f3.write("EMIN03000_EMAX10000_FM3.119_ASDCe_H0025_B01.cts.gz EMIN03000_EMAX10000_FM3.119_ASDCe_H0025_B01.exp.gz EMIN03000_EMAX10000_FM3.119_ASDCe_H0025_B01.gas.gz 25 -1 -1\n")
	f3.write("EMIN10000_EMAX50000_FM3.119_ASDCe_H0025_B01.cts.gz EMIN10000_EMAX50000_FM3.119_ASDCe_H0025_B01.exp.gz EMIN10000_EMAX50000_FM3.119_ASDCe_H0025_B01.gas.gz 25 -1 -1\n")
	gcfhe = fixgalcoeff.to_s + "," + fixgalcoeff.to_s + "," + fixgalcoeff.to_s + "," + fixgalcoeff.to_s + "," + fixgalcoeff.to_s
end
f3.close()

if fixgal.to_i >= 0
	if coordb.to_f < -10 or coordb.to_f > 10
		galcoeff = gcf
		galcoefffull = gcffull
		galcoeffhe = gcfhe
	end
end

puts "##############################################"
puts "### STEP 1 - determination of the parameters"
puts "##############################################"

if addff.to_i == 1
	minf = "25e-08"
else
	minf = "0e-08"
end
if mincatflux != nil
	minf = mincatflux.to_s
end

if fixgalisostep2.to_i >= 0
	cmdadd = " galmode2=3 isomode2=3 "
else
	cmdadd = " galmode2=0 isomode2=0 "
end

system("rm INT_"+fan+"*")
cmd = "multi6.rb FM3.119_ASDC2_"+irf+" " + maplist4name +" none INT_"+fan+" addcat=\""+ catline +"\" catminradius=" + minradius.to_s + " catminflux="+minf.to_s+" fluxcorrection=1 emin_sources=100 emax_sources=50000 edpcorrection=0.75 minimizertype=Minuit minimizeralg=Migrad minimizerdefstrategy=2 scanmaplist=" + sourcename + "," + fan + " fluxcorrection=1 " + cmdadd.to_s + " integratortype=" + inttype.to_s + " galcoeff=" + galcoeff.to_s
#galmode2=3 isomode2=3
puts cmd
system cmd
#fluxcorrection=1 isomode2=1 isomode2fit=2  # galmode2=1 isomode2=1 galmode2fit=1 isomode2fit=2

puts "####################################################"
puts "### multi preparation
puts "####################################################"


#reprocess with fan + ".multi" as input
newlistsource = fan + ".multi"
newlistsource2 = fan + ".ff1.multi"
alikeutils.rewriteMultiInputWithSingleSourcenewFixFlag(newlistsource, newlistsource2, sourcename, "1");
puts "############ NEW MULTI: " + newlistsource2

puts "####################################################"
puts "### STEP 2 - analysis of energy bins 00030-50000"
puts "####################################################"


suffix = "R" + inttype.to_s  + "_C" + format("%02d", minradius.to_f*10) + "-" + ARGV[1] + "-00030-50000-" + ARGV[3] + "-" + ARGV[4]
fan1 = prefix + "_FF1_" + suffix

cmd = "multi6.rb FM3.119_ASDC2_"+irf+" " + maplist4namefull +" " + newlistsource2 + " INT_"+fan1+" fluxcorrection=1 emin_sources=100 emax_sources=50000 edpcorrection=0.75 scanmaplist=" + sourcename + "," + fan1 + " minimizertype=Minuit minimizeralg=Migrad minimizerdefstrategy=2 fluxcorrection=1 " + cmdadd.to_s + " galmode2fit=0 isomode2fit=0 integratortype=" + inttype.to_s + " galcoeff=" + galcoefffull.to_s
puts cmd
system cmd

puts "####################################################"
puts "### STEP 2 - analysis of energy bins 00100-50000"
puts "####################################################"


suffix = "R" + inttype.to_s  + "_C" + format("%02d", minradius.to_f*10) + "-" + ARGV[1] + "-00100-50000-" + ARGV[3] + "-" + ARGV[4]
fan1 = prefix + "_FF1_" + suffix

cmd = "multi6.rb FM3.119_ASDC2_"+irf+" " + maplist4namehe +" " + newlistsource2 + " INT_"+fan1+" fluxcorrection=1 emin_sources=100 emax_sources=50000 edpcorrection=0.75 scanmaplist=" + sourcename + "," + fan1 + " minimizertype=Minuit minimizeralg=Migrad minimizerdefstrategy=2 fluxcorrection=1 " + cmdadd.to_s + " galmode2fit=0 isomode2fit=0 integratortype=" + inttype.to_s + " galcoeff=" + galcoeffhe.to_s
puts cmd
system cmd

