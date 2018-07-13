#! /usr/bin/ruby
#script for BUILD24
################ Mandatory and ordered parameters:
#00) filter (default FM3.119_2_I0025)
#01) maplist - Note that the extension .maplist4 is mandatory. See below for more details.
#02) listsource - extension .multi - file name of the list (in alike multi format - example 73.0e-08 80.27227 0.73223 2.1 1 2 3EGJ2033+4118). See below for more details. Or none if the .multi must be generated with the addcat procedure (otherwise will be appended)
#03) outfile - output file name
################ Optional parameters
#04) prefix (if specify a prefix, use only 1 map, if -999 use all the maps in the directory) - disabled with parameter maplist
#05) offaxis - off axix pointing (default 30) - set into .maplist4
#06) ranal   - radius of analysis (default 10)
#07) galcoeff     - gal coefficient (default -1) - set into .maplist4. If -2, use the same coeff of the input maplist4
#08) isocoeff     - iso coefficient (default -1) - set into .maplist4. If -2, use the same coeff of the input maplist4
#09) galmode     - gal mode, default 1 - See below for more details.
#10) isomode     - iso mode, default 1 - See below for more details.
#11) ulcl    - upper limit confidence level (default 2),  espressed as sqrt(TS)
#12) loccl   - source location contour confidence level (default 95 (%)confidence level) Vales: 99, 95, 68, 50
#13) flag    - a flag of the analysis (that is written in the final file)
#15) fixisogalstep0 - default 0 = do not calculate, otherwise specify the name of the source to calculate the gal and iso (if not specified by galcoeff or isocoeff parameters, otherwise use these values also in this step): calculate gal and iso setting fixflag=1 for the source under analysis
#16) findermultimode - default 0 = do not use, or the name of the source to be found. Analysis in 2 steps:
#	(1) ulcl=0, loccl=0, fixflag=3 to perform the first search
#	(2) use standard ulcl and loccl but with the new position of the source found in step (1)
#17) doublestep - default none = do not perform double step, otherwise doublestep=minsqrttsthr,secondstepfixflag,secondstepmaxradius (e.g. 3,3,0) and perform analysis in two steps (spot6 mode), where secondstepmaxradius is the last column of the .multi
#	(1) fixflag=1 for all the sources of the list
#	(2) generate a new list selecting the sources of the first list with sqrt(TS) > minsqrttsthr. The new list has fixflag = secondstepfixflag
#18) emin_sources, default 100: energy min of the input .multi
#19) emax_sources, default 50000: energy min of the input .multi
#20) listsourceextended
#21) checksourceposition - default 0 = do not use, or specify the name of the source. if fixflag > 1, check the position of the source. If the position of the source is too far, set fixflag=1 and recalculate the result. The parameters are the following:
#- name of the source
#- max distance wrt the initial position. If the ellipse is present, the radius is used instead of the user parameter
#22) galmode2
#23) galmode2fit
#24) isomode2
#25) isomode2fit
#26) minimizertype. Default Minuit
#27) minimizeralg. Default Migrad
#28) minimizerdefstrategy. Default 2 for Minuit
#29) mindefaulttolerance. Defaul 0.01
#xx) integratortype (1 gauss 2 gaussht 3 gausslefevre 4 gausslefevreht)
#30) edpcorrection, default 0 (no), 1 yes. EDP cleaning correction NOT IMPLEMENTED YET
#31) fluxcorrection, defaul 0 (no), 1 yes. Flux calculation correction for spectral shape in output - 2 input and output
#32) scanmaplist - default 0. Calculate one TS for each map of the maplist4 provided as input, or group by some set of maps (e.g. for fovbinnumer > 1) -> specify the name of the source and the prefix e.g. VELA,pl -> one MLE using pl law for each map of the maplist4 for VELA source. e.g. VELA,pl,5 -> one MLE each 5 maps of the maplist4
#33) (CAT) addcat. Specify the string of the source to be analysed". e.g. addcat="2.0e-07 34.7  -0.5  2.5 12 2 W44 0.0 1 2000.0 0.0". Remove sources with the same name, to avoid duplicate
#34) (CAT) catpath, the path of the cat file list (.multi). Default is /ANALYSIS3/catalogs/cat2.multi
#35) (CAT) catminflux, the min flux to be selected from the cat list
#36) (CAT) catminradius, the min radius to be selected from the cat list
#37) contourpoints, Number of points to determine the contour (0-400)
#38) testmode 0, 1 (gal-1sigma), 2 (gal+1sigma), 3 (iso-1sigma), 4 (iso+1sigma)

# MAPLIST
#Each line contains a set of maps:
#	<countsMap> <exposureMap> <gasMap> <offAxisAngle> <galCoeff> <isoCoeff>
#	where:
#	offAxisAngle is in degrees;
#	galcoeff and isocoeff are the coefficients for the galactic and isotropic diffuse components. 
#	If positive they will be considered fixed (but see galmode and isomode below).
#	The file names are separated by a space, so their name should not contain one.

# SOURCE LIST
#Each source is described by a line containing space separated values, in the following order:
# <flux> <l> <b> <spectral index> <fixFlag> <minSqrt(TS)> <name> <location limitation> <funtype> <par2> <par3> <par2 limit min> <par2 limit max> <par3 limit min> <par3 limit max>
# The first 4 values, flux in cm^-2 s^-1, galactic longitude and latitude in degrees, and spectral index of each source, represent the initial estimates of the values for that source (a positive number). According to the fixflag some or all of those values will be optimized by being allowed to vary.
#The flux estimates are relevant in the fitting process, as the sources are considered one by one starting with the one with the brightest initial flux value, regardless of the order they are given in the source file.
#The fixflag is a bit mask, each bit indicating whether the corresponding value is to be allowed to vary:
# fixFlag = 1 indicates the flux,
# fixFlag = 2 the position
# fixFlag = 4 the spectral index
# fixFlag = 8 par2 variable (and flux variable)
# fixFlag = 16 par3 variable (and flux variable)
# fixFlag = 32 force position to be variable only in Loop2
#The user may combine these values, but the flux will always be allowed to vary if at least one of the other values are.
#Examples:
#fixFlag = 0: everything is fixed. This is for known sources which must be included in order to search for other nearby sources.
#fixFlag = 1: flux variable, position fixed
#fixFlag = 2: only the position is variable, but AG_multi will let the flux vary too, so this is equivalent to 3.
#fixFlag = 3: flux and position variable, index fixed
#fixFlag = 4: index variable (and flux variable)
#fixFlag = 5: flux and index variable, position fixed
#fixFlag = 7: flux, position and index variable and also
#fixFlag = 28: index, par2 and par3 variable (and flux variable)
#fixFlag = 30: position index, par2 and par3 variable (and flux variable)

#minSqrt(TS) is the minimum acceptable value for the square root of TS: if the optimized significance of a source lies below this value, the source is considered undetected and will be ignored (set to flux = 0) when considering the other sources.
#The name of the source cannot contain any space. It is used to refer to this source in all the output files, and as part of the name of the source-specific output file.

#typefun
#0) "PL", "x^(-[index])"
#1) "PLExpCutoff", "x^(-[index]) * e^(- x / [par2])"
#2) "PLSuperExpCutoff", "x^(-[index]) * e^(- pow(x / [par2], [par3]))"
#3) "LogParabola", "( x / [par2] ) ^ ( -( [index] + [par3] * log ( x / [par2] ) ) )"

#GALMODE e ISOMODE
#   galmode and isomode are an integer value saying how the corresponding coefficients galcoeff or 
#   isoCoeff found in all the lines of the maplist are to be used:
#		0: all the coefficients are fixed.
#		1: all the coefficients are fixed if positive, variable if negative (the absolte value is the 
#		initial value). This is the default behaviour.
#		2: all the coefficients are variable, regardless of their sign.
#		3: all the coefficients are proportionally variable, that is the relative weight of their absolute value is kept.

#GALMODE2
#0 none
#1 set gal0 for L0 and gal1 for L1
#2 set gal0 for L0 and L1
#3 set gal1 for L0 and L1
#4 set gal1 - gal1err for L0 and L1
#5 set gal1 + gal1err for L0 and L1

#GALMODE2FIT
#	0 do not fit
#	1 pol0 fit
#	2 powerlaw fit

#ISOMODE2
#0 none
#1 set iso0 for L0 and gal1 for L1
#2 set iso0 for L0 and L1
#3 set iso1 for L0 and L1
#4 set iso1 - iso1err for L0 and L1
#5 set iso1 + iso1err for L0 and L1

#ISOMODE2FIT
#	0 do not fit
#	1 pol0 fit
#	2 powerlaw fit

#Extract list from CAT2
#extract_catalog.rb /ANALYSIS3/catalogs/cat2.multi 195.09 4.28 list.multi 0.1 1 5 0 10 0 0 25e-06

#Examples
#multi6.rb FM3.119_ASDCe_H0025 FM3.119_ASDCe_H0025_B01_00100-50000.maplist4 none RES1 addcat="2.0e-07 0.647894 61.9818  2.1 3 2 NS 0.0 0 0.0 0.0" catminflux=0e-08

#########################################Minimizers
#*** * Minuit (library libMinuit). Old version of Minuit, based on the TMinuit class. The list of possible algorithms are:
#Migrad (default one)
#Simplex
#Minimize (it is a combination of Migrad and Simplex)
#MigradImproved
#Scan
#Seek
#*** * Minuit2 (library libMinuit2). New C++ version of Minuit. The list of possible algorithm is :
#Migrad (default)
#Simplex
#Minimize
#Scan
#*** *Fumili . This is the same algorithm of TFumili, but implemented in the Minuit2 library.
#*** *GSLMultiMin (library libMathMore). Minimizer based on the Multidimensional Minimization routines of the Gnu Scientific Library (GSL). The list of available algorithms is
#BFGS2 (default) : second version of the vector Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm;
#BFGS : old version of the vector Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm;
#ConjugateFR : Fletcher-Reeves conjugate gradient algorithm;
#ConjugatePR : Polak-Ribiere conjugate gradient algorithm;
#SteepestDescent: steepest descent algorithm;
#*** * GSLMultiFit (library libMathMore). Minimizer based on the Non-Linear Least-Square routines of GSL. This minimizer can be used only for least-square fits.
#*** * GSLSimAn (library libMathMore). Minimizer based on simulated annealing.
#*** * Genetic (library libGenetic). Genetic minimizer based on an algorithm implemented in the TMVA package.

#Each minimizer can be configured using the ROOT::Math::MinimizerOptions class. The list of possible option that can be set are:
#* Minimizer type (MinimizerOptions::SetMinimizerType(const char *)) .
#* Minimizer algorithm (MinimizerOptions::SetMinimizerAlgorithm(const char *)).
#* Print Level (MinimizerOptions::SetPrintLevel(int )) to set the verbose printing level (default is 0).
#* Tolerance (MinimizerOptions::SetTolerance(double )) tolerance used to control the iterations.
#* Maximum number of function calls (MinimizerOptions::SetMaxFunctionCalls(int )).
#* Maximum number of iterations (MinimizerOptions::SetMaxIterations(int )). Note that this is not used by Minuit
#FCN Upper value for Error Definition (MinimizerOptions::SetMaxIterations(int )). Value in the minimization function used to compute the parameter errors. The default is to get the uncertainties at the 68% CL is a value of 1 for a chi-squared function minimization and 0.5 for a log-likelihood function.
#* Strategy (MinimizerOptions::SetStrategy(int )), minimization strategy used. For each minimization strategy Minuit uses different configuration parameters (e.g. different requirements in computing derivatives, computing full Hessian (strategy = 2) or an approximate version. The default is a value of 1. In this case the full Hessian matrix is computed only after the minimization.
#* Precision (MinimizerOptions::SetTolerance(double )). Precision value in the evaluation of the minimization function. Default is numerical double precision.

load ENV["AGILE"] + "/scripts/conf.rb"
load ENV["AGILE"] + "/scripts/MultiOutput6.rb"
datautils = DataUtils.new
fits = Fits.new

if ARGV[0].to_s == "help" || ARGV[0].to_s == "h" || ARGV[0] == nil
	system("head -128 " + $0 );
	exit;
end

filter = ARGV[0]
inputmaplist = ARGV[1];
listsource = ARGV[2];
baseoutfile = ARGV[3];
baseoutfile2 = ARGV[3];

p = Parameters.new
p.processInput(4, ARGV, filter)
alikeutils = AlikeUtils.new
multioutput = MultiOutput6.new

prefix = p.prefix
cts = ""
if prefix != -1 then
	cts = prefix.to_s + ".cts.gz";
	exp = prefix.to_s + ".exp.gz";
	gas = prefix.to_s + ".gas.gz";
end



#check energy range of input maps
emin_sin = 50000;
emax_sin = 0;
File.open(inputmaplist).each_line do | line |
	cts = line.split(" ")[0];
	fits.readFitsHeader(cts);
	if(emin_sin.to_f > fits.minenergy.to_f)
		emin_sin = fits.minenergy.to_f
	end
	if(emax_sin.to_f < fits.maxenergy.to_f)
		emax_sin = fits.maxenergy.to_f;
	end
	exp = line.split(" ")[1];
	gas = line.split(" ")[2];
end
puts "MAPS: emin " + emin_sin.to_s + " emax " + emax_sin.to_s

if p.addcat != ""
	if listsource == "none"
		listsource = baseoutfile + ".tmp.multi"
		system("touch " + listsource)
		puts listsource
	end
	
	listsourcetmp1 = listsource + ".tmp1"
	fout1 = File.new(listsourcetmp1, "w")
	fout1.write(p.addcat + "\n")
	fout1.close()
	
	ll = p.addcat.split(" ")
	listsourcetmp2 = listsource + ".tmp2"
	cmd = "extract_catalog.rb " + p.catpath + " " + ll[1].to_s + " " + ll[2].to_s + " " + listsourcetmp2 + " 0.1 1 5 0 10 0 0 " + p.catminflux.to_s + " " +  p.catminradius.to_s
	puts cmd
	system cmd
	
	listsourcetmp3 = listsource + ".tmp3"
	listsourcetmp4 = listsource + ".tmp4"
	alikeutils.appendMulti(listsourcetmp2, listsourcetmp1, listsourcetmp3, 0.1)
	alikeutils.appendMulti(listsource, listsourcetmp3, listsourcetmp4, 0.1)
	listsource = listsourcetmp4
end


#if different from default energy range of input source list, generate a new sourcelist with the energy range of the input maps
if(emin_sin.to_f != p.emin_sources.to_f or emax_sin.to_f != p.emax_sources.to_f)
	puts "change the flux"
	listsourceold = listsource
	listsource = format("en_%05d_%05d_", emin_sin, emax_sin) + listsource 
	fouts = File.new(listsource, "w")
	File.open(listsourceold).each_line do | line |
        l = line.split(" ")
        gamma = l[3].to_f
        fl = l[0].to_f
		
		c = p.emin_sources.to_f
		d = p.emax_sources.to_f
		a = emin_sin.to_f
		b = emax_sin.to_f
        
        p1 = fl * (gamma-1) / ( c ** (1-gamma) - d ** (1-gamma) )
        
        f1 = (p1 / (gamma-1)) * ( a ** (1-gamma) - b ** (1-gamma) )
        
        fouts.write(f1.to_s)
        for ii in 1..l.size()
        	fouts.write(" " + l[ii].to_s)
        end
        fouts.write("\n")
	end
	fouts.close()
end

stepi=1
prefixi = ""
ulcl = p.ulcl
loccl = p.loccl

if p.findermultimode != nil
	stepi=2
	prefixi = ".step1"
	ulcl = 0
	loccl = 0
end

doublestep_thr = nil
doublestep_fixflag = nil
if p.doublestep != nil
	stepi=2
	prefixi = ".step1"
	begin
		doublestep_thr = "3.0";
		doublestep_fixflag = "1";
		doublestep_maxradius = "0.0";
		doublestep_thr = p.doublestep.split(",")[0]
		doublestep_fixflag = p.doublestep.split(",")[1]
		doublestep_maxradius = p.doublestep.split(",")[2]
	rescue
		puts "error in doublestep parameters"
	end
end

lastoutfile = ""
outfile2 = baseoutfile2

cmd = "cp " + PATH + "share/AG_multi.par . "
datautils.execute(outfile2, cmd);
cmd = "cp " + PATH + "share/AG_multiext.par . "
datautils.execute(outfile2, cmd);

for i in 1..stepi
	#outfile
	outfile = baseoutfile
	

	#selezione delle calibration matrix
	filterbase = filter.split("_")[0];
	datautils.getResponseMatrix(filter);
	matrixconf = datautils.getResponseMatrixString(filter);
	
	#per prima cosa, se richiesto, si cerca il valore di gal e iso
	inputfilemaps = inputmaplist
	if p.fixisogalstep0 != nil
		outfile22 = outfile.to_s + ".step0"
		inputfilemaps22 = outfile.to_s + ".step0" + ".maplist4"
		alikeutils.rewriteMaplist(inputfilemaps, inputfilemaps22, p.galcoeff, p.isocoeff)
		
		#aggiorna il .multi mettendo a fixflag=1 solo la sorgente specificata in p.fixisogalstep0
		newlistsource = outfile.to_s + ".step0" + ".multi"
		alikeutils.rewriteMultiInputWithSingleSourceToAnalyze(listsource, newlistsource, p.fixisogalstep0, "1");
		
		if p.listsourceextended == "" 
			cmd = "export PFILES=.:$PFILES; " + PATH + "bin/AG_multi " + inputfilemaps22.to_s + " " + matrixconf.to_s + " "  + p.ranal.to_s + " " + p.galmode.to_s + " " + p.isomode.to_s +  " " + newlistsource.to_s + " " + outfile22.to_s + " " + ulcl.to_s + " " + loccl.to_s + " " + p.galmode2.to_s + " " + p.galmode2fit.to_s + " " + p.isomode2.to_s + " " + p.isomode2fit.to_s + " " + p.edpcorrection.to_s + " " + p.fluxcorrection.to_s + " " + p.minimizertype.to_s +  " " + p.minimizeralg.to_s + " " + p.minimizerdefstrategy.to_s + " " + p.mindefaulttolerance.to_s + " " + p.integratortype + " " + p.expratioevaluation.to_s + " " + p.minThreshold.to_s + " " +  p.maxThreshold.to_s + " " + p.squareSize.to_s + " " + p.contourpoints.to_s;
			datautils.execute(outfile2, cmd)	
		else
			cmd = "export PFILES=.:$PFILES; " + PATH + "bin/AG_multiext " + inputfilemaps22.to_s + " " + matrixconf.to_s + " "  + p.ranal.to_s + " " + p.galmode.to_s + " " + p.isomode.to_s +  " " + newlistsource.to_s + " " + p.listsourceextended + " " + outfile22.to_s + " " + ulcl.to_s + " " + loccl.to_s + " " + p.edpcorrection.to_s + " " + p.fluxcorrection.to_s + " " + p.minimizertype.to_s +  " " + p.minimizeralg.to_s + " " + p.minimizerdefstrategy.to_s + " " + p.mindefaulttolerance.to_s + " " + p.integratortype;
			datautils.execute(outfile2, cmd)
		end
		
		#step0b: prendi il valore di gal e iso calcolati e genera il maplist4 per le analisi successive
		multioutput.readDataSingleSource2(outfile22, p.fixisogalstep0);
		maplist = outfile22 + ".tmpmaplist4"
		
		alikeutils.rewriteMaplist(inputfilemaps, maplist, multioutput.galcoeff, multioutput.isocoeff)
		
	end
	
	if p.fixisogalstep0 == nil
		#copy the .maplist4 in .maplist
		maplist = outfile.to_s + ".maplist4"
		if p.galcoeff.to_f >= -1 and p.isocoeff.to_f >= -1
			alikeutils.rewriteMaplist(inputmaplist, maplist, p.galcoeff, p.isocoeff);
		else	
			system("cp " + inputmaplist + " " + maplist)
		end
	end
	
	#si esce dai due step precedenti con il maplist corretto

	if p.findermultimode != nil && i.to_i == 2
		prefixi = ""
		ulcl = p.ulcl
		loccl = p.loccl
	end
	
	if p.doublestep != nil && i.to_i == 2
		prefixi = ""
	end

	if maplist != nil
		inputfilemaps = outfile.to_s + prefixi.to_s + ".maplist4"
		#if p.findermultimode != nil
		if maplist != inputfilemaps
			if File.exists?(inputfilemaps) == false
				datautils.execute(outfile2, "cp " + maplist.to_s + " " + inputfilemaps.to_s)
			end
		end
		maplist = inputfilemaps
	end
	
	##....

	#list source
	newlistsource = outfile.to_s + prefixi.to_s + ".multi" 
	#if p.findermultimode != nil && i.to_i == 1
	#if File.exists?(newlistsource) == false #WARNING, BE CAREFUL
			datautils.execute(outfile2, "cp " + listsource.to_s + " " + newlistsource.to_s)
		#end
	#end
	if p.findermultimode != nil && i.to_i == 2
		#rewrite newlistsource
	
		multioutput.readDataSingleSource2(lastoutfile, p.findermultimode.split(",")[0])
	
		alikeutils.rewriteMultiInputWithNewCoordinatesSource(listsource, newlistsource, p.findermultimode.split(",")[0], multioutput.l_peak, multioutput.b_peak);
	
	end

	if p.doublestep != nil && i.to_i == 1
		#rewrite newlistsource with fixflag == 1 for all the sources except source names starting with _
		alikeutils.rewriteMultiListWithFixflag(listsource, newlistsource, 1)
	end
	
	if p.doublestep != nil && i.to_i == 2
		#rewrite newlistsource with fixflag = doublestep_fixflag, sqrt(TS) > doublestep_thr except for sources with source name starting with _
		cmd = "convertMultiResToInput.rb " + lastoutfile.to_s + " " + newlistsource.to_s + " " + doublestep_fixflag.to_s + " " + doublestep_thr.to_s + " " + doublestep_maxradius.to_s + " 90 1"
		puts cmd
		datautils.execute(outfile2, cmd)
	end

	newoutfile = outfile.to_s + prefixi.to_s
	
	if p.listsourceextended == ""
	
		cmd = "export PFILES=.:$PFILES; " + PATH + "bin/AG_multi " + inputfilemaps.to_s + " " + matrixconf.to_s + " "  + p.ranal.to_s + " " + p.galmode.to_s + " " + p.isomode.to_s +  " " + newlistsource.to_s + "  " + newoutfile + " " + ulcl.to_s + " " + loccl.to_s + " " + p.galmode2.to_s + " " + p.galmode2fit.to_s + " " + p.isomode2.to_s + " " + p.isomode2fit.to_s + " " + p.edpcorrection.to_s + " " + p.fluxcorrection.to_s + " " + p.minimizertype.to_s +  " " + p.minimizeralg.to_s + " " + p.minimizerdefstrategy.to_s + " " + p.mindefaulttolerance.to_s + " "  + p.integratortype + " " + p.expratioevaluation.to_s + " " + p.minThreshold.to_s + " " +  p.maxThreshold.to_s + " " + p.squareSize.to_s + " " + p.contourpoints.to_s;
		datautils.execute(outfile2, cmd)
	
	else
		
		cmd = "export PFILES=.:$PFILES; " + PATH + "bin/AG_multiext " + inputfilemaps.to_s + " " + matrixconf.to_s + " "  + p.ranal.to_s + " " + p.galmode.to_s + " " + p.isomode.to_s +  " " + newlistsource.to_s + " " + p.listsourceextended.to_s + " " + newoutfile + " " + ulcl.to_s + " " + loccl.to_s + " " + p.edpcorrection.to_s + " " + p.fluxcorrection.to_s + " " + p.minimizertype.to_s +  " " + p.minimizeralg.to_s + " " + p.minimizerdefstrategy.to_s + " " + p.mindefaulttolerance.to_s + " " + p.integratortype;
		datautils.execute(outfile2, cmd)
		
	end
	
	
	
	if p.checksourceposition != nil
		checksource_name = p.checksourceposition.split(",")[0]
		checksource_maxR = p.checksourceposition.split(",")[1]
		multioutput.readDataSingleSource2(newoutfile, checksource_name);
		#usa r se e' presente l'ellisse
		if multioutput.r.to_f > 0
			checksource_maxR = multioutput.r.to_f
		end
		
		#check if ff>1 and check if too far
		if multioutput.fix.to_i > 1 and multioutput.dist.to_f > checksource_maxR.to_f
			#repeat analysis
			newlistsource2 = "far." + newlistsource.to_s
			alikeutils.rewriteMultiInputWithSingleSourcenewFixFlag(newlistsource, newlistsource2, checksource_name, "1");
			
			#clean results
			puts "rm " + newoutfile + ".res* *.source"
            system("rm " + newoutfile + ".res*")
            system("rm " + newoutfile + "*.source*")
			
			if p.listsourceextended == ""
	
				cmd = "export PFILES=.:$PFILES; " + PATH + "bin/AG_multi " + inputfilemaps.to_s + " " + matrixconf.to_s + " "  + p.ranal.to_s + " " + p.galmode.to_s + " " + p.isomode.to_s +  " " + newlistsource2.to_s + "  " + newoutfile + " " + ulcl.to_s + " " + loccl.to_s + " " + p.galmode2.to_s + " " + p.galmode2fit.to_s + " " + p.isomode2.to_s + " " + p.isomode2fit.to_s + " " + p.edpcorrection.to_s + " " + p.fluxcorrection.to_s + " " + p.minimizertype.to_s +  " " + p.minimizeralg.to_s + " " + p.minimizerdefstrategy.to_s + " " + p.mindefaulttolerance.to_s + " " + p.integratortype.to_s + " " + p.expratioevaluation.to_s + " " + p.minThreshold.to_s + " " +  p.maxThreshold.to_s + " " + p.squareSize.to_s + " " + p.contourpoints.to_s;
				datautils.execute(outfile2, cmd)
				
			else
		
				cmd = "export PFILES=.:$PFILES; " + PATH + "bin/AG_multiext " + inputfilemaps.to_s + " " + matrixconf.to_s + " "  + p.ranal.to_s + " " + p.galmode.to_s + " " + p.isomode.to_s +  " " + newlistsource2.to_s + " " + p.listsourceextended.to_s + " " + newoutfile + " " + ulcl.to_s + " " + loccl.to_s + " " + p.edpcorrection.to_s + " " + p.fluxcorrection.to_s + " " + p.minimizertype.to_s +  " " + p.minimizeralg.to_s + " " + p.minimizerdefstrategy.to_s + " " + p.mindefaulttolerance.to_s + " " + p.integratortype.to_s;;
				datautils.execute(outfile2, cmd)
		
			end
		end
	end

	#cmd = "ruby ~/grid_scripts3/convertMultiResToReg.rb " + outfile.to_s + " white";
	#datautils.execute(outfile2, cmd)
	
	mout = MultiOutputList.new
	mout.readSources(newoutfile, newlistsource, p.flag);

	cmd = "ruby " + ENV["AGILE"] + "/scripts/convertMultiInputToReg.rb " + newlistsource.to_s + " green";
	datautils.execute(outfile2, cmd)

	
	fits.readFitsHeader(cts.to_s);
	cmd = "echo \"" + fits.utc_start + "\" >> " + newoutfile.to_s + "";
	datautils.execute(outfile2, cmd)

	cmd = "echo \"" + fits.utc_end + "\" >> " + newoutfile.to_s + "";
	datautils.execute(outfile2, cmd)
	
	cmd = "ruby " + ENV["AGILE"] + "/scripts/convertMultiResToInput.rb " + newoutfile.to_s + "  " + newoutfile.to_s + ".resfull.multi";
	datautils.execute(outfile2, cmd)
	
	cmd = "ruby " + ENV["AGILE"] + "/scripts/extractres2.rb " + newoutfile.to_s + ".resfull " + newoutfile.to_s + ".res2 none";
	datautils.execute(outfile2, cmd)

	#cmd = "ruby ~/grid_scripts3/convertMultiResToRegData.rb " + outfile.to_s;
	#datautils.execute(outfile2, cmd)

	system("cat " + newoutfile.to_s + "")

	lastoutfile = newoutfile
	
	if p.scanmaplist != 0
		
		sourcename = p.scanmaplist.split(",")[0]
		prefixscan = p.scanmaplist.split(",")[1]
		nfovbins = 1
		if p.scanmaplist.split(",").size == 3
			nfovbins = p.scanmaplist.split(",")[2].to_i
		end
		mouthe = MultiOutput6.new
		mouthe.readDataSingleSource(newoutfile + "_"+sourcename+".source")
		puts mouthe.sicalc
		puts mouthe.galcoeff.to_s + " +/- " + mouthe.galcoeff_err.to_s
		puts mouthe.isocoeff.to_s  + " +/- " + mouthe.isocoeff_err.to_s
		
		newoutfile2 = prefixscan + "_" + newoutfile;
		system("cp " + lastoutfile + " " + newoutfile2 + ".originalres")
		system("cp " + lastoutfile + ".multi " + newoutfile2 + ".multi.originalres")
		system("cp " + newoutfile + "_"+sourcename+".source " + newoutfile2 + "_"+sourcename+".source.originalres")
		fheso = File.new(newoutfile2 + ".multi", "w")
		File.open(newoutfile + ".multi").each_line do | line |
			ll = line.split(" ")
			if ll[6] == sourcename
				fheso.write(ll[0].to_s + " " + ll[1].to_s + " " + ll[2].to_s + " " + mouthe.sicalc + " 1 " + ll[5].to_s + " " + ll[6].to_s + " " + ll[7].to_s + " " + mouthe.typefun + " " + mouthe.par2 + " " + mouthe.par3 + "\n")
			else
				fheso.write(line)
			end
		end
		fheso.close()
		system("cp " + newoutfile2 + ".multi " + prefixscan + ".multi")
		
		cmd = $0 + " " +  filter + " " + newoutfile + ".maplist4 " + newoutfile2 + ".multi " + newoutfile2 + " galcoeff=" + mouthe.galcoeff + " isocoeff=" + mouthe.isocoeff + " fluxcorrection=" + p.fluxcorrection.to_s + " edpcorrection=" + p.edpcorrection.to_s + " emin_sources=" + emin_sin.to_s + " emax_sources=" + emax_sin.to_s + " contourpoints=" + p.contourpoints.to_s
		
		puts cmd
		system cmd
		
		
		indexmapl = 0;
		fres = File.new(prefixscan + ".spe", "w")
		puts "open " + newoutfile + ".maplist4"
		
		lines = File.readlines(newoutfile + ".maplist4")
		
		#File.open(newoutfile + ".maplist4").each_line do | line |
		
		iline = 0
		while iline < lines.size
			galcoeffout = ""
			isocoeffout = ""
			
			line = lines[iline]
			lname = prefixscan + "_" + line.split(" ")[0].split(".cts.gz")[0]
			puts "Processing ... " + line + " in file " + lname + "_s"
			fml = File.new(lname + "_s.maplist4", "w")
			for iii in 0...nfovbins
				fml.write(lines[iline + iii])
				if galcoeffout != ""
					galcoeffout += ","
				end
				gcfd1 = mouthe.galcoeff.split(",")[iline + iii];
				gcfd1_err = mouthe.galcoeff_err.split(",")[iline + iii];
				if p.testmode.to_i == 1
					gcfd1 = (gcfd1.to_f - gcfd1_err.to_f).to_s
				end
				if p.testmode.to_i == 2
					gcfd1 = (gcfd1.to_f + gcfd1_err.to_f).to_s
				end
				galcoeffout += gcfd1.to_s
				
				if isocoeffout != ""
					isocoeffout += ","
				end
				icfd1 = mouthe.isocoeff.split(",")[iline + iii];
				icfd1_err = mouthe.isocoeff_err.split(",")[iline + iii];
				if p.testmode.to_i == 3
					icfd1 = (icfd1.to_f - icfd1_err.to_f).to_s
				end
				if p.testmode.to_i == 4
					icfd1 = (icfd1.to_f + icfd1_err.to_f).to_s
				end
				isocoeffout += icfd1.to_s
			end
			fml.close()
			puts "isocoeffout: "
			puts isocoeffout
			iline = iline + nfovbins
			
			#build galcoeff and isocoeff
			
			
			#lname = prefixscan + "_" + line.split(" ")[0].split(".cts.gz")[0]
			
			#cmd = $0 + " " +  filter + " " + lname + "_s.maplist4 " + newoutfile2 + ".multi " + lname + " galcoeff=" + mouthe.galcoeff.split(",")[indexmapl] + " isocoeff=" + mouthe.isocoeff.split(",")[indexmapl] + " fluxcorrection=" + p.fluxcorrection.to_s + " edpcorrection=" + p.edpcorrection.to_s + " emin_sources=" + emin_sin.to_s + " emax_sources=" + emax_sin.to_s
			cmd = $0 + " " +  filter + " " + lname + "_s.maplist4 " + newoutfile2 + ".multi " + lname + " galcoeff=" + galcoeffout + " isocoeff=" + isocoeffout + " fluxcorrection=" + p.fluxcorrection.to_s + " edpcorrection=" + p.edpcorrection.to_s + " emin_sources=" + emin_sin.to_s + " emax_sources=" + emax_sin.to_s
			puts cmd
			system cmd
			
			mouthe2 = MultiOutput6.new
			mouthe2.readDataSingleSource(lname + "_"+sourcename+".source")
			#(0)sqrtts (1)flux[ph/cm2/s] (2)flux_err (3)erg[erg/cm2/s] (4)Erg_err (5)Emin (6)Emax (7)E_log_center (8)exp (9)flux_ul (10)spectral_index (11)spectral_index_err (12)id_detection or expcorfactor
			alpha = -mouthe2.sicalc.to_f;
			e_min = mouthe2.energyrange.split("..")[0].to_f
			e_max = mouthe2.energyrange.split("..")[1].to_f
			ec = 0
			if alpha != -2
				ec = ( 1 / (alpha+2) ) *  (( e_max ** (alpha+2))  -  (e_min ** (alpha+2)) )  /  ( ( 1 / (alpha+1) ) * ( e_max ** (alpha+1))  -  e_min ** (alpha+1));
			else
				ec = Math.log(e_max/e_min) / ( ( ( e_max ** (alpha+1))  -  (e_min ** (alpha+1))) / (alpha + 1 ) );
			end
			
			ec = format("%.3f", e_min - ec)
			
			fres.write(mouthe2.sqrtTS.to_s + " " + mouthe2.flux.to_s + " " + mouthe2.flux_error.to_s + " " + " " + mouthe2.erglog.to_s + " " + mouthe2.erglog_error.to_s + " " + e_min.to_s + " " + e_max.to_s + " " + ec.to_s + " " + mouthe2.expspectracorfactor.to_s + " " + mouthe2.flux_ul.to_s + " " + mouthe2.sicalc.to_s + " " + mouthe2.sicalc_error.to_s + " " + mouthe2.likelihood1.to_s + " " + mouthe2.erglogul.to_s + "\n")
			indexmapl = indexmapl.to_i + 1
		end
		fres.close()
		
		puts "-----------------------------------------------"
		system("cat " + newoutfile2 + ".multi.originalres > " + prefixscan + ".sum")
		puts "-----------------------------------------------"
		system("cat " + newoutfile2 + ".originalres >> " + prefixscan + ".sum")
		puts "-----------------------------------------------"
		system("cat " + newoutfile2 + "_"+sourcename+".source.originalres >> " + prefixscan + ".sum")
		puts "-----------------------------------------------"
		system("cat " + prefixscan + ".spe >> " + prefixscan + ".sum")
		puts "-----------------------------------------------"
		system("cat " + prefixscan + ".sum")

	end

end

cmd = "rm ./AG_multi.par"
#datautils.execute(outfile2, cmd);

exit(0)
