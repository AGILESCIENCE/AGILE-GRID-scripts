class Parameters
	public
		def initialize()
			@addcat = ""
			@catminflux = "25e-08"
			@catminradius = "0"
			@integratortype = "1"
			@catpath = "/ANALYSIS3/catalogs/cat2.multi"
			@checksourceposition = nil
			@execap = 0
			@radius = 10
			@t1s = 0
			@t2s = 10
			@t1b = 600
			@shiftt1b = 0
			@t2b = 600
			@shiftt2b = 0
			@filtercode = 5
			@timeslot = 0
			@timeslotstart = 0
			@timeslotend = 0
			@filter = "FM3.119_ASDCe_" + TYPE_MATRIX
			@listsourceextended = ""
			@prefix=-1
			@timetype="TT"
			@useEDPmatrixforEXP = 0
			@flag="F0"
			@findermultimode = nil
			@enablespotfinder = 0;
			@gammaextractbin = 3600
			@disablespotfinder = 1;
			@makealikesingle = 0;
			@binsize = 0.3;
			@fixisogal = 0;
			@skytype = 4;
			@eboundaryIF = 400;
			@list = "none";
			@mapsize = 40;
			@skymapH = "";
			@skymapL = "";
			@spectralindex = 2.1;
			@emin = 100;
			@emax = 50000;
			@fovradmin = 0;
			@fovradmax = 60;
			@albedorad = 80;
			@minTS = SQRTTS_CUTTONEXTSTEP.to_f * SQRTTS_CUTTONEXTSTEP.to_f;
			@SFsmooth = 3;
			@SFnumberofsource = 5;
			@SFuseint = 0;
			@SFradiousremove = 0.7;
			@SFalgspotfinder = 1;
			@SFantype = 3;
			@SFrextract = -1;
			@ITbinstepiterative = 4;
			@ITscanit = 15;
			@ITradit = 0;
			@ranal = 10;
			@expstep = -1;
			@timestep = 160;
			@fovbinnumber = 1;
			@proj = "ARC";
			@phasecode = -1;
			@timebinsize = 99999999999999;
			@makelc = nil;
			@ulcl = 2.0;
			@loccl = 95;
			@token="localhost,-1"
			@cluster=0
			@onlymaps=0
			@outputtype=0
			@opmode=10
			@blocks=0
			@nruns=1
			@seed=0
			@dq = 0
			
			@fixisogalstep0 = nil;
			@doublestep = nil;
			@galcoeff = "-1";
			@isocoeff = "-1";
			@galmode = 1;
			@isomode = 1;
			@multitype = 4;
			@lpointing = -1;
			@bpointing = -1;
			@outputres = nil
			
			#only for multi
			@offaxis = 30;
			@maplist = nil;
			@timelist = "None";
			@maplistgen = "None";
			@outfile = nil
			@outfile2 = nil
			@energybin = 0;
			@emin_sources = 100
			@emax_sources = 50000
			@galmode2 = 0
			@galmode2fit = 0
			@isomode2 = 0
			@isomode2fit = 0
			@minimizertype = "Minuit"
			@minimizeralg = "Migrad"
			@minimizerdefstrategy = 2
			@mindefaulttolerance = 0.01
			@edpcorrection = 1.0
			@fluxcorrection = 1
			@scanmaplist = 0
		end
		
		def setPhaseCode(tstop)
			if @phasecode == -1
				if tstop.to_f >= 182692800.0 
					@phasecode = 6 #SPIN
				else
					@phasecode = 18 #POIN
				end
			end
		end
		
		def checksourceposition
			@checksourceposition
		end
		
		def integratortype
			@integratortype
		end
		
		def catminflux
			@catminflux
		end
		
		def catminradius
			@catminradius
		end
		
		def addcat
			@addcat
		end
		
		def catpath
			@catpath
		end
		
		def skytype
            @skytype
        end
        
        def eboundaryIF
        	@eboundaryIF
        end
		
		def scanmaplist
			@scanmaplist
		end
        
        def dq
            @dq
        end
        
        def emin_sources
            @emin_sources
        end
        
        def emax_sources
            @emax_sources
        end
        
        def useEDPmatrixforEXP
        	@useEDPmatrixforEXP
        end
        
        def timetype
            @timetype
        end
		
		def prefix
			@prefix
		end
		
		def flag
			@flag
		end
		
		def listsourceextended
			@listsourceextended
		end
		
		def filter
			@filter
		end

		def blocks
			@blocks
		end

		def opmode
			@opmode
		end
		
		def nruns
			@nruns
		end
		
		def seed
			@seed
		end
		
		def execap
			@execap
		end
		
		def radius
			@radius
		end
		
		def t1s
			@t1s
		end
		
		def t2s
			@t2s
		end
		
		def t1b
			@t1b
		end
		
		def t2b
			@t2b
		end
		
		def shiftt1b
			@shiftt1b
		end
		
		def shiftt2b
			@shiftt2b
		end
		
		def filtercode
			@filtercode
		end

		def timeslot
			@timeslot
		end
		
		def timeslotstart
			@timeslotstart
		end
		
		def timeslotend
			@timeslotend
		end

		def outputtype
			@outputtype
		end
		
		
		def findermultimode
			@findermultimode
		end
		
		def maplist
			@maplist
		end
		
		def maplistgen
			@maplistgen
		end
		
		def timelist
			@timelist
		end
		
		def onlymaps
			@onlymaps
		end
		
		def skymapH
			@skymapH
		end
		
		def skymapL
			@skymapL
		end
		
		def setonlymaps(val)
			@onlymaps=val
		end
		
		def cluster
			@cluster
		end
		
		def outfile
			@outfile
		end
		
		def outfile2
			@outfile2
		end
		
		def offaxis
			@offaxis
		end
		
		def ulcl
			@ulcl
		end
		
		def token
			@token
		end
		
		def loccl
			@loccl
		end
		
		def fixisogalstep0
			@fixisogalstep0 
		end
		
		def doublestep
			@doublestep
		end
		
		def galcoeff
			@galcoeff
		end
		
		def galmode
			@galmode
		end
		
		def isomode
			@isomode
		end
		
		def proj
			@proj
		end
		
		def phasecode
			@phasecode
		end
		
		def timebinsize
			@timebinsize
		end
		
		def makelc
			@makelc
		end

		def isocoeff
			@isocoeff
		end

		def multitype
			@multitype
		end

		def lpointing
			@lpointing
		end

		def bpointing
			@bpointing
		end

		def outputres
			@outputres
		end
		
		def enablespotfinder
			@enablespotfinder
		end
		
		def gammaextractbin
			@gammaextractbin
		end
		
		def disablespotfinder
			@disablespotfinder
		end
		
		def makealikesingle
			@makealikesingle
		end
		
		def binsize
			@binsize
		end
		
		def fixisogal
			@fixisogal
		end
		
		def list
			@list
		end
		
		def mapsize
			@mapsize
		end
		
		def spectralindex
			@spectralindex
		end
		
		def emin
			@emin
		end
		
		def emax
			@emax
		end
		
		def fovradmin
			@fovradmin
		end
		
		def fovradmax
			@fovradmax
		end
		
		def albedorad
			@albedorad
		end
		
		def minTS
			@minTS
		end
		
		def SFsmooth
			@SFsmooth
		end
		
		def SFnumberofsource
			@SFnumberofsource
		end
		
		def SFrextract
			@SFrextract
		end
		
		def SFuseint
			@SFuseint
		end
		
		def SFradiousremove
			@SFradiousremove
		end
		
		def SFalgspotfinder
			@SFalgspotfinder
		end
		
		def SFantype
			@SFantype
		end
		
		def ITbinstepiterative
			@ITbinstepiterative
		end
		
		def ITscanit
			@ITscanit
		end
		
		def ITradit
			@ITradit
		end
		
		def ranal
			@ranal
		end
		
		def expstep
			@expstep
		end
		
		def timestep
			@timestep
		end
		
		def fovbinnumber
			@fovbinnumber
		end
		
		def energybin
			@energybin
		end
		
		def galmode2
			@galmode2
		end
		
		def galmode2fit
			@galmode2fit
		end
		
		def isomode2
			@isomode2
		end
		
		def isomode2fit
			@isomode2fit
		end
		
		def minimizertype
			@minimizertype
		end
		
		def minimizeralg
			@minimizeralg
		end
		
		def minimizerdefstrategy
			@minimizerdefstrategy
		end
		
		def mindefaulttolerance
			@mindefaulttolerance
		end
		
		def edpcorrection
			@edpcorrection
		end
		
		def fluxcorrection
			@fluxcorrection
		end
		
		def processInput(startindex, s, filter)
			for i in startindex...s.size
				if s[i] == nil
					break;
				else
					processLine(s[i]);
				end
			end
			
			if filter != nil
				@filter = filter
			end
			
			initparam();
		end
		
		def buildCommandLine()
				a = "offaxis=" + @offaxis.to_s + " "
				if @outfile != nil
					a = a + "outfile=" + @outfile.to_s + " "
				end
				if @maplist != nil
					a = a + "maplist=" + @maplist.to_s + " "
				end
				if @maplistgen != nil
					a = a + "maplistgen=" + @maplistgen.to_s + " "
				end
				if @timelist != nil
					a = a + "timelist=" + @timelist.to_s + " "
				end
				if @checksourceposition != nil
					a = a + "checksourceposition=" + @checksourceposition + " "
				end
				a = a + "enablespotfinder=" + @enablespotfinder.to_s + " "
				a = a + "prefix=" + @prefix.to_s + " "
				a = a + "useEDPmatrixforEXP=" + @useEDPmatrixforEXP + " "
				a = a + "timetype=" + @timetype.to_s + " "
				a = a + "flag=" + @flag.to_s + " "
				a = a + "filter=" + @filter.to_s + " "
				a = a + "dq=" + @dq.to_s + " "
				a = a + "execap=" + @execap.to_s + " "
				a = a + "radius=" + @radius.to_s + " "
				a = a + "t1s=" + @t1s.to_s + " "
				a = a + "t2s=" + @t2s.to_s + " "
				a = a + "t1b=" + @t1b.to_s + " "
				a = a + "t2b=" + @t2b.to_s + " "
				a = a + "shiftt1b=" + @shiftt1b.to_s + " "
				a = a + "shiftt2b=" + @shiftt2b.to_s + " "
				a = a + "filtercode=" + @filtercode.to_s + " "
				a = a + "timeslot=" + @timeslot.to_s + " "
				a = a + "timeslotstart=" + @timeslotstart.to_s + " "
				a = a + "timeslotend=" + @timeslotend.to_s + " "
				a = a + "gammaextractbin=" + @gammaextractbin.to_s + " "
				a = a + "disablespotfinder=" + @disablespotfinder.to_s + " "
				a = a + "makealikesingle=" + @makealikesingle.to_s + " "
				a = a + "binsize=" + @binsize.to_s + " "
				a = a + "skymapH=" + @skymapH.to_s + " "
				a = a + "skymapL=" + @skymapL.to_s + " "
				a = a + "skytype=" + @skytype.to_s + " "
				a = a + "token=" + @token.to_s + " "
				a = a + "ulcl=" + @ulcl.to_s + " "
				a = a + "loccl=" + @loccl.to_s + " "
				a = a + "phasecode=" + @phasecode.to_s + " "
				a = a + "galmode=" + @galmode.to_s + " "
				a = a + "isomode=" + @isomode.to_s + " "
				a = a + "timebinsize=" + @timebinsize.to_s + " "
				a = a + "listsourceextended=" + @listsourceextended.to_s + " "
				if @makelc != nil
					a = a + "makelc=" + @makelc.to_s +  " "
				end
				a = a + "proj=" + @proj.to_s + " "
				a = a + "fixisogal=" + @fixisogal.to_s + " "
				a = a + "list=" + @list.to_s + " "
				a = a + "mapsize=" + @mapsize.to_s + " "
				a = a + "spectralindex=" + @spectralindex.to_s + " "
				a = a + "emin=" + @emin.to_s + " "
				a = a + "emax=" + @emax.to_s + " ";
				a = a + "fovradmin=" + @fovradmin.to_s + " "
				a = a + "fovradmax=" + @fovradmax.to_s + " "
				a = a + "albedorad=" + @albedorad.to_s + " "
				a = a + "minTS=" + @minTS.to_s + " "
				a = a + "SFsmooth=" + @SFsmooth.to_s + " "
				a = a + "SFnumberofsource=" + @SFnumberofsource.to_s + " "
				a = a + "SFuseint=" + @SFuseint.to_s + " "
				a = a + "SFrextract= " + @SFrextract.to_s + " " 
				a = a + "SFradiousremove=" + @SFradiousremove.to_s + " "
				a = a + "SFalgspotfinder=" + @SFalgspotfinder.to_s + " "
				a = a + "SFantype=" + @SFantype.to_s + " "
				a = a + "ITbinstepiterative=" + @ITbinstepiterative.to_s + " "
				a = a + "ITscanit=" + @ITscanit.to_s + " " 
				a = a + "ITradit=" + @ITradit.to_s + " "
				a = a + "expstep=" + @expstep.to_s + " "
				a = a + "timestep=" + @timestep.to_s + " "
				a = a + "ranal=" + @ranal.to_s + " "
				a = a + "fovbinnumber=" + @fovbinnumber.to_s + " "
				a = a + "galcoeff=" + @galcoeff.to_s + " "
				a = a + "fixisogalstep0=" + @fixisogalstep0.to_s + " "
				a = a + "doublestep=" + @doublestep.to_s + " "
				a = a + "isocoeff=" + @isocoeff.to_s + " "
				a = a + "multitype=" + @multitype.to_s + " "
				a = a + "lpointing=" + @lpointing.to_s + " "
				a = a + "bpointing=" + @bpointing.to_s + " "
				a = a + "energybin=" + @energybin.to_s + " "
				a = a + "cluster=" + @cluster.to_s + " "
				a = a + "onlymaps=" + @onlymaps.to_s + " "
				a = a + "outputtype=" + @outputtype.to_s + " "
				a = a + "opmode=" + @opmode.to_s + " "
				a = a + "blocks=" + @blocks.to_s + " "
				a = a + "nruns=" + @nruns.to_s + " "
				a = a + "seed=" + @seed.to_s + " "
				a = a + "emin_sources=" + @emin_sources.to_s + " "
				a = a + "emax_sources=" + @emax_sources.to_s + " "
				a = a + "findermultimode=" + @findermultimode.to_s + " "
				if @outputres != nil
					a = a + "outputres=" + @outputres.to_s + " "
				end
				@buildCommandLine = a
		end
		

		def processLine(argv)
			keyw = argv.split("=")[0];
			value = argv.split("=")[1];
			puts keyw.to_s + " " + value.to_s
			case keyw
				when "catminflux"
					@catminflux = value
				when "catminradius"
					@catminradius = value
				when "integratortype"
					@integratortype = value;
				when "addcat"
					@addcat = value;
				when "catpath"
					@catpath = value;
				when "scanmaplist"
					@scanmaplist = value;
				when "checksourceposition"
					@checksourceposition = value;
				when "outputtype"
					@outputtype = value;
				when "listsourceextended"
					@listsourceextended = value;
				when "opmode"
					@opmode = value;
				when "nruns"
					@nruns = value;
				when "dq"
					@dq = value;
				when "execap"
					@execap = value;
				when "radius"
					@radius = value;
				when "t1s"
					@t1s = value;
				when "t2s"
					@t2s = value;
				when "t1b"
					@t1b = value;
				when "t2b"
					@t2b = value;
				when "shiftt1b"
					@shiftt1b = value;
				when "shiftt2b"
					@shiftt2b = value;
				when "filtercode"
					@filtercode = value;
				when "timeslot"
					@timeslot = value;
				when "timeslotstart"
					@timeslotstart = value;
				when "timeslotend"
					@timeslotend = value;
				when "filter"
					@filter = value;
				when "useEDPmatrixforEXP"
					@useEDPmatrixforEXP = value;
				when "seed"
					@seed = value;
				when "timetype"
					@timetype = value;
				when "prefix"
					@prefix = value;
				when "flag"
					@flag = value;
				when "blocks"
					@blocks = value;
				when "findermultimode"
					@findermultimode = value;
				when "outfile"
					@outfile = value;
					@outfile2=value
				when "offaxis"
					@offaxis = value;
				when "skymapH"
					@skymapH = value;
				when "skymapL"
					@skymapL = value;
				when "skytype"
					@skytype = value;
				when "maplist"
					@maplist = value;
				when "maplistgen"
					@maplistgen = value;
				when "timelist"
					@timelist = value;
				when "cluster"
					@cluster = value;
				when "onlymaps"
					@onlymaps = value;
				when "gammaextractbin"
					@gammaextractbin = value;
				when "enablespotfinder"
					@enablespotfinder = value;
					if @enablespotfinder.to_i == 1
						@disablespotfinder = 0;
					else
						@disablespotfinder = 1;
					end
				when "disablespotfinder"
					@disablespotfinder = value;
					if @disablespotfinder.to_i == 1
						@enablespotfinder = 0;
					else
						@enablespotfinder = 1;
					end
				when "makealikesingle"
					@makealikesingle = value;
				when "binsize"
					@binsize = value;
				when "token"
					@token = value;
				when "ulcl"
					@ulcl = value;
				when "loccl"
					@loccl = value;
				when "phasecode"
					@phasecode = value;
				when "galmode"
					@galmode = value;
				when "isomode"
					@isomode = value;
				when "timebinsize"
					@timebinsize = value;
				when "makelc"
					@makelc = value;
				when "proj"
					@proj = value;
				when "fixisogal"
					@fixisogal = value;
				when "list"
					@list = value;
				when "mapsize"
					@mapsize = value;
				when "spectralindex"
					@spectralindex = value;
				when "emin"
					@emin = value;
				when "emax"
					@emax = value;
				when "emin_sources"
					@emin_sources = value;
				when "emax_sources"
					@emax_sources = value;
				when "fovradmin"
					@fovradmin = value;
				when "fovradmax"
					@fovradmax = value;
				when "albedorad"
					@albedorad = value;
				when "minTS"
					@minTS = value;
				when "SFsmooth"
					@SFsmooth = value;
				when "SFnumberofsource"
					@SFnumberofsource = value;
				when "SFuseint"
					@SFuseint = value;
				when "SFradiousremove"
					@SFradiousremove = value;
				when "SFrextract"
					@SFrextract = value;
				when "SFalgspotfinder"
					@SFalgspotfinder = value;
				when "SFantype"
					@SFantype = value;
				when "ITbinstepiterative"
					@ITbinstepiterative = value;
				when "ITscanit"
					@ITscanit = value;
				when "ITradit"
					@ITradit = value;
				when "expstep"
					@expstep = value;
				when "timestep"
					@timestep = value;
				when "ranal"
					@ranal = value;
				when "fovbinnumber"
					@fovbinnumber = value;
				when "galcoeff"
					@galcoeff = value;
				when "fixisogalstep0"
					@fixisogalstep0 = value;
				when "doublestep"
					@doublestep = value;
				when "isocoeff"
					@isocoeff = value;
				when "multitype"
					@multitype = value;
				when "lpointing"
					@lpointing = value;
				when "bpointing"
					@bpointing = value;
				when "outputres"
					@outputres = value;
				when "energybin"
					@energybin = value;
				when "eb"
					@energybin = value;
				when "galmode2"
					@galmode2 = value;
				when "galmode2fit"
					@galmode2fit = value;
				when "isomode2"
					@isomode2 = value;
				when "isomode2fit"
					@isomode2fit = value;
				when "minimizertype"
					@minimizertype = value;
				when "minimizeralg"
					@minimizeralg = value;
				when "minimizerdefstrategy"
					@minimizerdefstrategy = value;
				when "mindefaulttolerance"
					@mindefaulttolerance = value;
				when "edpcorrection"
					@edpcorrection = value;
				when "fluxcorrection"
					@fluxcorrection = value;
				else
					puts "Keyword " + argv.to_s + " error."
					#exit;
			end
			
			
		end
		
		def initparam()
			if @emax.to_f > 50000
				puts "Error in the energy range: the maximum energy should be 50000"
				exit(1)
			end
			
			if @ITradit.to_f == 0
				@ITradit = @mapsize.to_i / 2.0 - @ranal.to_i
			end
			
			if @dq.to_i == 1
				@albedorad = 80
				@fovradmax = 60
			end
			if @dq.to_i == 2
				@albedorad = 80
				@fovradmax = 50
			end
			if @dq.to_i == 3
				@albedorad = 90
				@fovradmax = 60
			end
			if @dq.to_i == 4
				@albedorad = 90
				@fovradmax = 50
			end
			
			if @expstep.to_f == -1
				@expstep = 4;
				if @binsize.to_f == 0.05
					@expstep = 20;
				end
				if @binsize.to_f == 0.1
					@expstep = 10;
				end
				if @binsize.to_f == 0.2
					@expstep = 5;
				end
				if @binsize.to_f == 0.25
					@expstep = 4;
				end
				if @binsize.to_f == 0.3
					@expstep = 3;
				end
				if @binsize.to_f == 0.5
					@expstep = 2;
				end
				if @binsize.to_f == 1
					@expstep = 1;
				end
			end
			
			if @loccl.to_f == 95
				@loccl = 5.99147
			end
			if @loccl.to_f == 99
				@loccl = 9.21034
			end
			if @loccl.to_f == 50
				@loccl = 1.38629
			end
			if @loccl.to_f == 68
				@loccl = 2.29575
			end

			
			#check only available combination of IRF matrices and SKY maps
			fconf = @filter.split("_");
			#skytype: 0 SKY000-1 + SKY000-5, 1 gc_allsky maps + SKY000-5, 2 SKY000-5, 3 SKY001 (old galcenter, binsize 0.1, full sky), 4 SKY002 (new galcenter, binsize 0.1, full sky)
			
			#I0023: skytype=0, skytype=1, skytype=2, skytype=4
			if @skytype.to_i == 3
				puts "Error: skytype=3 not available. Set skytype=4"
				@skytype = 4
			end
			
			#I0025: skyteype=4
			if fconf[2] == "I0025" 
				if @skytype.to_i != 4
					puts "Error: only skytype=4 with IRF=I0025 is available. Set skytype=4"
					@skytype = 4
				end
				@eboundaryIF = 400
			end
			
			if fconf[2] == "H0025" 
				if @skytype.to_i != 4
					puts "Error: only skytype=4 with IRF=H0025 is available. Set skytype=4"
					@skytype = 4
				end
				@eboundaryIF = 300
			end
			
		end
end

