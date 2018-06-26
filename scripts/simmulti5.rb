#! /usr/bin/ruby
#script for BUILD24
#0) filter/archive/calibMatrix (FM3.119_ASDCe_I0023)
#1) maplist sim
#2) list of sources (in .multi format) for simulation
#3) maplist ana
#4) list of sources (in .multi format) for analysis
#5) outfile
#6) calmatrix: format SKYXXX.SFILTER_CALMATRIX

#optional
#6) opmode	 Integer	 Operation Mode
#7) blocks	 Integer	 Block
#8) nruns	 Integer	 Number of runs
#9) seed	Integer	 Seed

#10) ranal
#11) galmode
#12) isomode
#13) ulcl
#14) loccl
#15) outputtype (defalt 0)


#Parametri di input del task
#Name	Type	Description
#opmode	 Integer	 Operation Mode
#block	 Integer	 Block
#nruns	 Integer	 Number of runs
#seed	Integer	 Seed
#sarfile	 String	 SAR file name
#edpfile	 String	 EDP file name
#psdfile	String	 PSD file name
#maplistsim	String	Simulation map list
#srclistsim	String	Simulation source list
#outfile	 String	 Output file name
#maplistanalysis	String	Analysis map list
#srclistanalysis	String	Analysis source list
#ranal	Real	Radius of analysis region
#galmode	Integer	 Diffuse emission mode
#isomode	Integer	 Isotropic emission mode
#ulcl	Real	 Upper limit confidence level
#loccl	 Real	 Source location contour confidence level

#opmode
#This program option is a bit mask with the following meaning:
#1 concise - Generate a single log file rather that a file for each source and each analysis
#2  skipanalysis -  Do the simulation only, skip any analysis
#4 doubleanalysis - Do the analysis in two steps, the second using the sources obtained from the first
#8 savemaps - Save to disk all the maps evaluated at run time and used for the analysis
#Please note the following:
#If the bit 1 is off the bits 0 and 2 are ignored.
#If the bit 1 is on and the bit 3 is off no output is generated
#If bit 1 is set, command line options 10 to 17 are optional and, if present, are ignored

#parameters.emin, parameters.emax, parameters.skytype

#Istruzioni: devono essere presenti la exp e la gas map e la maplist4. VERSIONE: loccl to 95

#AG_multisim performs a Monte Carlo simulation generating count maps to be analyzed with the maximum likelihood method.
#The counts maps can be analyzed on the fly and saved to disk.
#The simulation is based on a list of exposure and gas maps and a set of sources, and for each exposure map a counts map is generated.
#These maps can be either catenatd and analyzed as a set of separate maps, or grouped, summed up and analyzed as single maps.
#The map groups are obtained considering the first N maps, then again other N maps starting from the second, and so on.
#In either way, the counts map obtained are analyzed the same way AG_multi would do.

#The command line options includes:
#4 parameters to define how the simulation is to be performed
#3 file names that are calibration files used by both the simulation and the analysis
#3 parameters giving the input file names and part of the output file names.
#7 optional parameters contain the additional information required by the analysis, if this step is required

#block
#The block parameter defines how to group the maps. There are two cases:
#block is zero
#In this case all the maps are catenated together and analyzed the same way AG_multi normally does.
#maplistanalysis must have the same number of entries as maplistsim

#block greater than zero
#In this case the analysis is performed using as count map the sum of block maps, and as the exposure map the sum of block exposure maps.
#If there are N maps, then the maps are summed up N-block times, starting from the first, from the 2nd, and so on.
#In this case the maplistanalysis may contain only one line, and both its counts and exposure maps would be ignored and replaced with the ones #evaluated as explained above.

#nrun
#The simulation is repeated nrun times.

#seed
#This is the unsigned number used to start the random generator. You can give either zero, meaning the random generator is based on the machine #time, or any number you like.
#If you give zero AG_multisim prints to screen the number actually used, so that you can repeat the same simulation later on.
#To perform a different analysis on the same simulation you need to use the same seed but also the same calibration files, the same maps and the same sources, that are the parameters 5 to 9.

#The two MapLists
#The MapList is a text file whose format is described in the AG_multi User Manual.
#In the context of AG_multisim we have two map lists, one for the simulation and one for the analysis.
#The counts maps in the simulation MapList are ignored, so they can be set to 'None'.
#Everything else in the simulation MapList is used to generated the counts maps, one for each line of the MapList.
#If block is zero the counts maps will replace those of the analysis MapList, if any.
#If block is greater than zero the analysis MapList will be only one line long (other lines would be ignored), in which the counts and the exposure maps would be replaced by those generated by the simulation.

#Example 1. Both MapLists may loook like:
#None map1.exp.fits map1.gas.fits 30 -0.8 -7.2
#None map2.exp.fits map2.gas.fits 30 -0.8 -7.2
#The minus sign of the cooefficients would be ignored for the simulation, but marks them as variable during the analysis.
#Two counts map are generated and inserted in the map list for the analysis.
#In a case like this example the same MapList file name can be provided both as simulation and analysis MapList.
#If block>0, the second line of the MapLists would be ignored for the analysis, as well as the exposure file name.

#Example 1. If block is not zero, the two maps may look like:
#simulation MapList:
#None map1.exp.fits map1.gas.fits 30 0.82 7.2
#None map2.exp.fits map2.gas.fits 40 0.78 6.8
#analysis MapList
#None  None map.gas.fits 35 -0.7 -7.0
#Two counts map are generated by the simulation, and their sum replaces the corresponding place holder in the analysis MapList (the first 'None'.
#The sum if the corresponding exposure maps will replace the second 'None'.
#In this case it is necessary to provide two differrent MapLists because the angles and the diffuse coefficients differ from one another.
#The diffuse coefficients of the simulation MapList may or may not have the minus sign, but those of the analysis MapList must have it if we want #those coefficients to be variable.

#Output file names
#If the counts map are saved to disk, their name is given by:
#<iteration><mapProg><outfile>.cts.gz
#where:
#<iteration> is a 10 digits number between 1 and nruns
#<mapProg> is a 3 digits number between 1 and the number of exposure maps given
#<outfile> is the outfile command line option

#If block>0 <mapProg> is between 1 and mapCount-block+1. 001 is the sum of the first block counts map, 002 is from the second to block+1 and so on. In theis case also the exposure map is saved, and its name is the same with .exp.gz suffix.

#If opmode&concise==concise two files are generated at each iteration. Their name is
#<iteration><mapProg><outfile>
#and
#<iteration><mapProg><outfile>_<source_name>
#where <source_name> is the name of each variable source, and <mapProg> is present only if block>0.

#If opmode&concise==0 only one output file is generated for the entire simulation. The first entry of each line is the same as <iteration>.
#In this case, if block>0, there will be (mapCount-block+1)*varSourcesCount lines beginning with the same number (NOTE: This has probably to be changed)

load ENV["AGILE"] + "/scripts/conf.rb"
datautils = DataUtils.new

if ARGV[0].to_s == "help" || ARGV[0].to_s == "h" || ARGV[0] == nil
	system("head -125 " + $0 );
	exit;
end

filter = ARGV[0];
maplistsim = ARGV[1]
listsourcesim = ARGV[2];
maplistana = ARGV[3]
listsourceana = ARGV[4];
outfile = ARGV[5]
calmatrix = ARGV[6]
p = Parameters.new
p.processInput(7, ARGV, filter)

#selezione della matrix
filterbase = filter.split("_")[0];
datautils.getResponseMatrix(filter);
sarmatrix = datautils.sarmatrix
edpmatrix = datautils.edpmatrix
psdmatrix = datautils.psdmatrix
matrixconf = datautils.getResponseMatrixString(filter);



# outfile2 = prefix.to_s + "_" + listsourcesim.to_s + "_iso" + iso.to_s
logfile = outfile.to_s + ".log"

cmd = "cp " + PATH + "share/AG_multisim5.par . "
datautils.execute("", cmd);
						
cmd = "export PFILES=.:$PFILES; " + PATH + "bin/AG_multisim5 " + p.opmode.to_s + " " + p.blocks.to_s + " " + p.nruns.to_s + " " + p.seed.to_s + " "  + matrixconf.to_s + " " + maplistsim.to_s + " " + listsourcesim.to_s + " " + outfile.to_s + " " + maplistana.to_s + " " + listsourceana.to_s + " " + p.ranal.to_s + " " + p.galmode.to_s + " " + p.isomode.to_s + " " + p.ulcl.to_s + " " + p.loccl.to_s + " " + calmatrix.to_s + " " + ENV["AGILE"];
datautils.execute(logfile, cmd);

