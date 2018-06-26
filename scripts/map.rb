#! /usr/bin/ruby
#script for BUILD24
#00) filter DIR (es: FM3.119_2_I0023, FM3.119_ASDCe_I0023, FM3.119_ASDCSTDf_I0023, FM3.119_ASDCSTDk_I0023)
#01) output file name prefix
#02) time start: contact or TT or MJD or UTC (default TT, use timetype)
#03) time end: contact or TT or MJD or UTC (default TT, use timetype)
#04) l (map center)
#05) b (map center)
# Optional
### COMMON
#06) timetype: CONTACT, MJD, UTC, TT, default TT
#07) emin: energy min in MeV, if energybin=0, default 100
#08) emax: energy max in MeV, if energybin=0, default 500000
#09) fovradmin: fov rad min, to be used also with fovbinnumber, default 0
#10) fovradmax: fov rad max, to be used also with fovbinnumber, default 60
#11) mapsize: map size (diameter of the map in degree), default 40
#12) binsize: bin size, default 0.3
#13) proj: projection ARC or AIT, default ARC
#14) albedorad: default 80
#15) (SEL) fovbinnumber: number of bins between fovradmin and fovradmax. Dim = (fovradmax-fovradmin)/fovbinnumber, default 1
#16) (SEL) energybin or eb: default 0, 
#	=1 activate [00030-00050] [00050-00100], [00100-00400], [00400-01000], [01000-03000], [03000, 10000], [10000, 50000])
#	=2 activate [00030-00050] [00050-00100], [00100-00200], [00200-00400], [00400-01000], [01000-03000], [03000, 500000])
#	=3 activate [00100-00200], [00200-00400], [00400-01000], [01000-03000])
#	=4 activate [00100-00200], [00200-00400], [00400-01000])
#	=5 activate [00050-00100], [00100-00200], [00200-00400], [00400-01000], [01000-03000])
#	=6 activate [00100-00200], [00200-00400], [00400-01000], [01000-03000], [03000, 500000])
#	=7 activate [00050-00100], [00100-00200], [00200-00400], [00400-01000], [01000-03000], [03000, 500000])
#	=8 activate [00050-00100], [00100-00400], [00400-01000], [01000-03000], [03000-50000])
#	=9 activate [00050-00100], [00100-00400], [00400-01000], [01000-03000], [03000, 100000], [10000, 50000])
#	=10 activate [00100-00400], [00400-01000], [01000-03000], [03000, 100000], [10000, 50000])
#	=11 activate [00400-01000], [01000-03000], [03000, 100000], [10000, 50000])
#	=12 activate [00400-01000], [01000-03000], [03000, 100000])
#   =13 activate [00100-00400], [00400-01000], [01000-03000]
#   =14 activate [00100-00400], [00400-01000]
#   =15 activete [00030-00050], [00050-00100], [00100-50000]
#   =16 activate [00050-00100], [00100-00[3|4]00], [00[3|4]00-01000], [01000-03000], [03000, 100000], [10000, 50000])
#   =17 activate [00050-00100], [00100-00200], [00200-00400], [00400-01000], [01000-03000], [03000, 10000], [10000, 50000])
#17) phasecode: optional, default 2. If -1 => automatic determination ==>  if (time end  > 182692800.0 (MJD 55119.5, UTC 2009-10-15T12:00:00, fine pointing) && phasecode == -1) then phasecode = 6 (SPIN) else phasecode = 18 (POIN)
#18) timelist: a file with a list of tstart/stop
#19) timebinsize: optional, default 999999999
#20) makelc: optional, default 0
#21) lpointing: optinal, default -1
#22) bpointing: optional, default -1
#23) maplistgen: filename of a file with  mapspec.fovradmin >> mapspec.fovradmax >> mapspec.emin >> mapspec.emax >> mapspec.index
### EXPOSURE
#23) (SEL) useEDPmatrixforEXP: use the EDP matrix to generate expmap, default 0. WARNING: 1 does not work
#24) expstep: step size of exp map gen, it depends by binsize (e.g. 0.3->3, 0.25->4, 0.1->10)
#25) spectralindex: spectral index of exp map, default 2.1
#26) timestep: LOG file step size, default 160
### GAS MAP
#27) (SEL) skytype: 0 SKY000-1 + SKY000-5, 1 gc_allsky maps + SKY000-5, 2 SKY000-5, 3 SKY001 (old galcenter, binsize 0.1, full sky), 4 SKY002 (new galcenter, binsize 0.1, full sky)
#28) skymapL: sky map low resolution
#29) skymapH: sky map high resolution
#30) dq: data quality, default 0. dq = 1 -> albedorad=80,fovradmax=60. dq = 2 -> albedorad=80,fovradmax=50. dq = 3 -> albedorad=90,fovradmax=60. dq = 4 -> albedorad=90,fovradmax=50. dq=0 use standard albedorad and fovradmax
#31) filtercode = 0 G+L+S, filtercode=5 only G
#32) execap, AP: exec aperture photometry (default 0) 
#33) timeslot, AP: timeslot for aperture photometry
#34) ranal, AP: radius of analysis to extract event, for aperture photometry

#Lo script crea le mappe mancanti, e se ne crea almeno uno aggiunge la corrispondente riga nel .maplitsX. Attenzione quindi alle duplicazioni

#Note sul phasecode:
#2 -> per lo spinning, esclude la SAA con metodo conteggi AC
#6 -> per lo spinngin, esclude la SAA in base ad intensità campo magnetico (TPZ)
#18 -> per il pointing, esclude la SAA e il recovery
# Normalmente usate il phasecode = 6 nei dati in spinning. Questo phasecode esclude i fotoni presenti nella SAA ridefinita con i conteggi dell'AC. Se invece vuoi usare la vecchia definizione della SAA (in base all'intensità del campo magnetico così come definito da TPZ) devi usare il phasecode = 6.

load ENV["AGILE"] + "/scripts/conf.rb"

datautils = DataUtils.new
parameters = Parameters.new

if ARGV[0].to_s == "help" || ARGV[0] == nil || ARGV[0] == "h"
	system("head -68 " + $0 );
	exit;
end

filter = ARGV[0];
name = ARGV[1];
contact0 = ARGV[2];
contact1 = ARGV[3];
l = ARGV[4].to_f + 0.000001;
b = ARGV[5];

datautils.extractFilterDir(filter)
filterdir = datautils.filterdir

filterbase2 = filter.split("_")[0] + "_" + filter.split("_")[1];

parameters.processInput(6, ARGV, filter)

emin1 = parameters.emin;
emax1 = parameters.emax;
datautils.getSkyMatrix(filter, emin1, emax1, parameters.skytype)
skymap =  datautils.skymatrix;
skymapL = parameters.skymapL;
skymapH = parameters.skymapH;
indexlog = datautils.logindex(filterbase2)
indexfilter = datautils.evtindex(filterbase2)
puts "indexfilter: " + indexfilter.to_s
puts "Sky map: " + skymap.to_s;
eminarr = Array.new
emaxarr = Array.new
puts "parameters.energybin " + parameters.energybin.to_s
if parameters.energybin.to_i == 0
	energybinnumber = 1;
end
puts "parameters.eboundaryIF " + parameters.eboundaryIF.to_s

if parameters.energybin.to_i == 1
	eminarr = [30, 50, 100, parameters.eboundaryIF, 1000, 3000, 10000]
	emaxarr = [50, 100, parameters.eboundaryIF, 1000, 3000, 10000, 50000]
	energybinnumber = 7
end

if parameters.energybin.to_i == 2
	eminarr = [30,  50, 100, 200,  parameters.eboundaryIF, 1000,  3000, 10000]
	emaxarr = [50, 100, 200, parameters.eboundaryIF, 1000, 3000, 10000, 50000]
	energybinnumber = 8
end

if parameters.energybin.to_i == 3
	eminarr = [100, 200,  parameters.eboundaryIF, 1000]
	emaxarr = [200, parameters.eboundaryIF, 1000, 3000]
	energybinnumber = 4
end

if parameters.energybin.to_i == 4
	eminarr = [100, 200,  parameters.eboundaryIF]
	emaxarr = [200, parameters.eboundaryIF, 1000]
	energybinnumber = 3
end

if parameters.energybin.to_i == 5
	eminarr = [50,  100, 200,  parameters.eboundaryIF, 1000,  3000]
	emaxarr = [100, 200, parameters.eboundaryIF, 1000, 3000, 10000]
	energybinnumber = 6
end

if parameters.energybin.to_i == 6
	eminarr = [100, 200,  parameters.eboundaryIF, 1000,  3000, 10000]
	emaxarr = [200, parameters.eboundaryIF, 1000, 3000, 10000, 50000]
	energybinnumber = 6
end

if parameters.energybin.to_i == 7
	eminarr = [50,  100, 200,  parameters.eboundaryIF, 1000,  3000, 10000]
	emaxarr = [100, 200, parameters.eboundaryIF, 1000, 3000, 10000, 50000]
	energybinnumber = 7
end

if parameters.energybin.to_i == 8
        eminarr = [50,  100, parameters.eboundaryIF,  1000, 3000]
        emaxarr = [100, parameters.eboundaryIF, 1000, 3000, 50000]
        energybinnumber = 5
end

if parameters.energybin.to_i == 9
	eminarr = [50,  100, parameters.eboundaryIF,  1000, 3000,  10000]
	emaxarr = [100, parameters.eboundaryIF, 1000, 3000, 10000, 50000]
	energybinnumber = 6
end

if parameters.energybin.to_i == 10
	eminarr = [100,  parameters.eboundaryIF, 1000, 3000, 10000]
	emaxarr = [parameters.eboundaryIF, 1000, 3000, 10000, 50000]
	energybinnumber = 5
end

if parameters.energybin.to_i == 11
	eminarr = [parameters.eboundaryIF, 1000, 3000, 10000]
	emaxarr = [1000, 3000, 10000, 50000]
	energybinnumber = 4
end

if parameters.energybin.to_i == 12
	eminarr = [parameters.eboundaryIF, 1000, 3000]
	emaxarr = [1000, 3000, 10000]
	energybinnumber = 3
end

if parameters.energybin.to_i == 13
	eminarr = [100,  parameters.eboundaryIF, 1000]
	emaxarr = [parameters.eboundaryIF, 1000, 3000]
	energybinnumber = 3
end

if parameters.energybin.to_i == 14
	eminarr = [100,  parameters.eboundaryIF]
	emaxarr = [parameters.eboundaryIF, 1000]
	energybinnumber = 2
end

if parameters.energybin.to_i == 15
	eminarr = [30, 50, 100]
	emaxarr = [50, 100, 50000]
	energybinnumber = 3
end

if parameters.energybin.to_i == 16
	eminarr = [50, 100,  parameters.eboundaryIF, 1000, 3000, 10000]
	emaxarr = [100, parameters.eboundaryIF, 1000, 3000, 10000, 50000]
	energybinnumber = 6
end

if parameters.energybin.to_i == 17
	eminarr = [ 50, 100, 200,  400, 1000,  3000, 10000]
	emaxarr = [100, 200, 400, 1000, 3000, 10000, 50000]
	energybinnumber = 7
end

index_name_cor = BASEDIR_ARCHIVE.to_s + "/DATA/INDEX/3901.cor.index"
puts "index name cor: " + index_name_cor;
tstart = 0
tstop = 0	
if(parameters.timetype == "CONTACT")
	#estrazione dei tempi min e max dal corfileindex

	datautils.extractTimeMinMaxForContact(index_name_cor, contact0);
	tstart = datautils.tmin;

	datautils.extractTimeMinMaxForContact(index_name_cor, contact1);
	tstop = datautils.tmax;
end
if(parameters.timetype == "TT")
	tstart = contact0;
	tstop = contact1
end
if(parameters.timetype == "MJD")
	tstart = datautils.time_mjd_to_tt(contact0);
	tstop = datautils.time_mjd_to_tt(contact1);
end

if(parameters.timetype == "UTC")
	tstart = datautils.time_utc_to_tt(contact0)
	tstop = datautils.time_utc_to_tt(contact1)
end

puts "TMIN: " + tstart.to_s;
puts "TMAX: " + tstop.to_s;

if tstart.to_f == 0
	puts "Error in TMIN, exit"
	exit(1)
end
if tstop.to_f == 0
	puts "Error in TMAX, exit"
	exit(1)
end

#phasecode

parameters.setPhaseCode(tstop)

#selezione della sar matrix
filterbase = filter.split("_")[0]
datautils.getResponseMatrix(filter);
sarmatrix = datautils.sarmatrix;
edpmatrix = datautils.edpmatrix;

nameori = name

flc = nil
if parameters.makelc != nil
	flc = File.new(parameters.makelc, "w")
end

time = tstart
index = 0;
puts time
puts tstop
while time.to_f < tstop.to_f
	t0 = time.to_f
	t1 = time.to_f + parameters.timebinsize.to_f
	if t1.to_f > tstop.to_f
		t1 = tstop
	end
	puts index.to_s + " " + t0.to_s + " " + t1.to_s;

	time = t1;
	index = index + 1;

	if index.to_i == 1 && t1 == tstop
		name = nameori
	else
		name = format("%04d", index) + "_" + nameori.to_s;
	end
	cts = name + ".cts.gz";
	exp = name + ".exp.gz";
	gas = name + ".gas.gz";
	int = name + ".int.gz";

	conffile4 = name + ".maplist4";
	
	prefix = name


	if parameters.execap.to_i == 0
		fconf4 = File.new(conffile4, "w")
	end

	#if File.exists?(conffile4) == false
	#	fconf4 = File.new(conffile4, "w")
	#else
	#	fconf4 = File.new(conffile4, "a")
	#end
	
	puts "energybinnumber=" + energybinnumber.to_s

	createdexpmap = false
	for stepi in 1..parameters.fovbinnumber.to_i

		for stepe in 1..energybinnumber.to_i
			emin = parameters.emin;
			emax = parameters.emax; 
			prefixi = ""
			if energybinnumber != 1
				emin = eminarr[stepe-1]
				emax = emaxarr[stepe-1]
				prefixi = format("EMIN%05d", emin.to_i) + "_EMAX" + format("%05d", emax.to_i) + "_"
			end
			datautils.getSkyMatrix(filter, emin, emax, parameters.skytype)
			skymap =  datautils.skymatrix;
			#if skymapL == ""
				skymapL =  datautils.skymatrixL;
			#end
			#if skymapH == ""
				skymapH =  datautils.skymatrixH;
			#end
			puts "skymap: " + skymap.to_s
			puts "skymap LOW res: " + skymapL.to_s
			puts "skymap HIGH res: " + skymapH.to_s
		
			if parameters.fovbinnumber.to_i == 1
				fovmin = parameters.fovradmin
				fovmax = parameters.fovradmax
				bincenter = 30
				if parameters.lpointing.to_i != -1 and parameters.lpointing.to_i != -999
					bincenter = datautils.distance(l, b, parameters.lpointing, parameters.bpointing);
				end
			else
				bincenter = ( ((parameters.fovradmax.to_f-parameters.fovradmin.to_f) / parameters.fovbinnumber.to_f) * stepi.to_f ) - ((parameters.fovradmax.to_f-parameters.fovradmin.to_f) / parameters.fovbinnumber.to_f) / 2.0
				fovmin = bincenter.to_f - ((parameters.fovradmax.to_f-parameters.fovradmin.to_f) / parameters.fovbinnumber.to_f) / 2.0
				fovmax = bincenter.to_f + ((parameters.fovradmax.to_f-parameters.fovradmin.to_f) / parameters.fovbinnumber.to_f) / 2.0
				prefixi = prefixi.to_s + format("%02d", stepi.to_i) + "_"
			end
			puts "Map generation for " + stepi.to_s + " fovmin " + fovmin.to_s + " and fovmax " + fovmax.to_s + " with center " + bincenter.to_s

			cts2 = prefixi.to_s  + cts.to_s
			exp2 = prefixi.to_s  + exp.to_s
			gas2 = prefixi.to_s  + gas.to_s
			int2 = prefixi.to_s  + int.to_s
			
			createdmap = false
		
			#BUILD22
			lonpole = 180.0
			
			if parameters.execap.to_i == 0

				# generate maps
				if File.exists?(cts2) == false
					cmd = "cp " + PATH + "share/AG_ctsmapgen5.par . "
					datautils.execute(prefix, cmd);
				
					cmd = "export PFILES=.:$PFILES; " + PATH + "bin/AG_ctsmapgen5 " + cts2.to_s  + " " + indexfilter.to_s + " " + parameters.timelist.to_s  + "  " + parameters.mapsize.to_s + " " + parameters.binsize.to_s + " "  + l.to_s + " " + b.to_s + " " + lonpole.to_s + " " + " " + parameters.albedorad.to_s + " " + parameters.phasecode.to_s + " " + parameters.filtercode.to_s + " "  + parameters.proj.to_s + " "+ t0.to_s + " " + t1.to_s + " " + emin.to_s + " " + emax.to_s + " " + fovmin.to_s + " " + fovmax.to_s;
					datautils.execute(prefix, cmd);
					createdmap = true
					cmd = "rm ./AG_ctsmapgen5.par"
					#datautils.execute(prefix, cmd);
				end
				sarmatrixfull = PATHMODEL + sarmatrix
				edpmatrixfull = " None "
				if parameters.useEDPmatrixforEXP.to_i == 1
					edpmatrixfull =  PATHMODEL + edpmatrix
				end
				if File.exists?(exp2) == false
			
					if createdexpmap == false
						cmd = "cp " + PATH + "share/AG_expmapgen5.par . "
						datautils.execute(prefix, cmd);
				
						#campionamento ogni 0.1 sec del file di LOG
						#maplist = mapstream >> mapspec.fovradmin >> mapspec.fovradmax >> mapspec.emin >> mapspec.emax >> mapspec.index;
						cmd = "export PFILES=.:$PFILES; " + PATH + "bin/AG_expmapgen5 " + exp2.to_s + " " + indexlog.to_s  +  " " + sarmatrixfull.to_s +  " " + edpmatrixfull.to_s + " " + parameters.maplistgen.to_s + " " + " " + parameters.timelist.to_s + " " + parameters.mapsize.to_s + " " + parameters.binsize.to_s  + " " + l.to_s + " " + b.to_s + " " + lonpole.to_s + " " + parameters.albedorad.to_s + " 0.5 360.0 5.0 " + parameters.phasecode.to_s + " " +  parameters.proj.to_s + " " + parameters.expstep.to_s + " " + parameters.timestep.to_s +  " " + parameters.spectralindex.to_s + " " + t0.to_s + " " + t1.to_s + " " + emin.to_s + " " + emax.to_s  + " " + fovmin.to_s + " " + fovmax.to_s;
						datautils.execute(prefix, cmd);
						createdmap = true
						if parameters.maplistgen != "None" 
							createdexpmap = true
						end
					
						cmd = "rm ./AG_expmapgen5.par"
						#datautils.execute(prefix, cmd);
					end
				end
		
				if File.exists?(gas2) == false
					cmd = "cp " + PATH + "share/AG_gasmapgen5.par . "
					datautils.execute(prefix, cmd);
					cmd = "export PFILES=.:$PFILES; " + PATH + "bin/AG_gasmapgen5 " + exp2.to_s + " " + gas2.to_s + " " + skymapL.to_s + " " + skymapH.to_s;
					datautils.execute(prefix, cmd);
					createdmap = true
					cmd = "rm ./AG_gasmapgen5.par"
					#datautils.execute(prefix, cmd);
				end
		
				if File.exists?(int2) == false
					cmd = "cp " + PATH + "share/AG_intmapgen5.par . "
					datautils.execute(prefix, cmd);
					cmd = "export PFILES=.:$PFILES; " + PATH + "bin/AG_intmapgen5 " + exp2.to_s + " " + int2.to_s + " " + cts2.to_s;
					datautils.execute(prefix, cmd);
					cmd = "rm ./AG_intmapgen5.par"
					#datautils.execute(prefix, cmd);
				end
				
				if createdmap == true
					fconf4.write(cts2.to_s + " " + exp2.to_s + " " + gas2.to_s + " " + bincenter.to_s + " -1 -1 \n")
				end
		
				if parameters.makelc != nil
					datautils.extractlc(flc, prefixi.to_s + prefix.to_s, l, b, t0, t1)
					system("rm " + prefixi.to_s + prefix.to_s + "*")
				end
			else

				# execute ap
				listfile=prefixi.to_s + name + ".ap"
				sarmatrixfull = PATHMODEL + sarmatrix
				edpmatrixfull = " None "
				if parameters.useEDPmatrixforEXP.to_i == 1
					edpmatrixfull =  PATHMODEL + edpmatrix
				end
				cmd = "cp " + PATH + "share/AG_ap5.par . "
				datautils.execute(prefix, cmd);
				cmd = "export PFILES=.:$PFILES; "+PATH+"bin/AG_ap5 "+listfile+" "+indexlog.to_s+" "+indexfilter.to_s+" "+sarmatrixfull.to_s+" "+edpmatrixfull.to_s+" "+
					  parameters.timelist.to_s+" "+parameters.ranal+" "+l.to_s+" "+b.to_s+" "+lonpole.to_s+" "+" "+parameters.albedorad.to_s+" 0.5 360.0 5.0 "+
					  parameters.phasecode.to_s+" "+parameters.timestep.to_s+" "+parameters.spectralindex.to_s+" "+t0.to_s+" "+t1.to_s+" "+emin.to_s+" "+emax.to_s+" "+
					  fovmin.to_s+" "+fovmax.to_s+" "+parameters.filtercode.to_s+" "+parameters.timeslot.to_s
				datautils.execute(prefix, cmd);
			end
		
			puts "###"
			puts stepi-1
		end
	end

	if parameters.execap.to_i == 0
		fconf4.close()
	end
end

if parameters.makelc != nil
	flc.close()
end

exit(0)

