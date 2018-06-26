load ENV["AGILE"] + "/scripts/DataUtils.rb"

class MultiOutput

	def setCalcPhase(t0, period)
		@calcorbitalphase_t0 = t0;
		@calcorbitalphase_period = period;
	end

	def readDataSingleSource2(res, sourcename)
		puts res.to_s + "_" + sourcename  + ".source"
		readDataSingleSource(res.to_s + "_" + sourcename);
	end
	
	#nameout = nome del file che contiene i dati	
	def readDataSingleSource(nameout)
			if nameout.include?(".source") == false
				nameout = nameout + ".source"
			end
			
			datautils = DataUtils.new
			#puts "nameout: " +  nameout;
			@l = -1;
			@b = -1;
			@r = -1;
			@ell_a = -1;
			@ell_b = -1;
			@ell_phi = -1;
			
			@sicalc = 0
			@sicalc_error = 0

			#read upper limit
			index2 = 0;
			@r = -1;
			@galcoeff = "-1"
			@galcoeff_err = "-1"
			@galcoeffzero = "-1"
			@galcoeffzero_err = "-1"	
			@isocoeff = "-1"
			@isocoeff_err = "-1"
			@isocoeffzero = "-1"
			@isocoeffzero_err = "-1"	
			@fitdata = "-1, -1, -1, -1, -1, -1, -1"
			@mlestep_1_ext1_res = -1
			@mlestep_2_step1_res = -1
			@mlestep_3_ext2_res = -1
			@mlestep_4_step2_res = -1
			@mlestep_5_contour_res = -1
			@mlestep_6_index_res = -1
			@mlestep_7_ul_res = -1
			@mlestep_1_ext1_cts = -1
			@mlestep_2_step1_cts = -1
			@mlestep_3_ext2_cts = -1
			@mlestep_4_step2_cts = -1
			@mlestep_5_contour_cts = -1
			@mlestep_6_index_cts = -1
			@mlestep_7_ul_cts = -1
			@mlestep_res = ""
			@mlestep_cts = ""
			indexstart = 0
			File.open(nameout).each_line do | line |
				index2 = index2 + 1;
				if line[0] == "!" or line[0] == 33
					indexstart = index2
					next
				end
				
				lll = line.split(" ")
				
				if index2.to_i == indexstart + 1
					@label =lll[0];
					@fix =lll[1];
					@si_start = lll[2];
					@ulconflevel = lll[3];
					@srcconflevel = lll[4];
					@startL = lll[5];
					@startB = lll[6];
					@startFlux = lll[7];
					@lmin = lll[9];
					@lmax = lll[11];
					@bmin = lll[14];
					@bmax = lll[16];
				end
				if index2.to_i == indexstart + 2
					@sqrtTS =lll[0];
				end
				if index2.to_i == indexstart + 3
					@l_peak = lll[0];
					@b_peak = lll[1];
					@dist = lll[2];
				end
				if index2.to_i == indexstart + 4
					@l = lll[0]
					@b = lll[1]
					@distellipse = lll[2];
					@r = lll[3];
					@ell_a = lll[4];
					@ell_b = lll[5];
					@ell_phi = lll[6];
					@fullellipseline = format("%.2f %.2f %.2f %.2f %.2f %.2f %.2f", @l, @b, @distellipse, @r, @ell_a, @ell_b, @ell_phi)
				end
				if index2.to_i == indexstart + 5
					@counts = lll[0]
					@counts_error = lll[1]
					@counts_error_p = lll[2]
					@counts_error_m = lll[3]
					@counts_ul = lll[4];
				end
				if index2.to_i == indexstart + 6
					@flux = lll[0]
					@flux_error = lll[1]
					@flux_error_p = lll[2]
					@flux_error_m = lll[3]
					@flux_ul = lll[4];
					@exposure = lll[5]
				end
				if index2.to_i == indexstart + 7
					@sicalc = lll[0]
					@sicalc_error = lll[1]
				end
				
				if index2.to_i == indexstart + 8
					@fit_cts = lll[0]
					@fit_fcn0 = lll[1]
					@fit_fcn1 = lll[2]
					@fit_edm0 = lll[3]
					@fit_edm1 = lll[4]
					@fit_iter0 = lll[5]
					@fit_iter1 = lll[6]
				end
				
				if index2.to_i == indexstart + 9
					@galcoeff = lll[0]
					@galcoeff_err = lll[1]
				end
				
				if index2.to_i == indexstart + 10
					@galcoeffzero = lll[0]
					@galcoeffzero_err = lll[1]
				end
				
				if index2.to_i == indexstart + 11
					@isocoeff = lll[0]
					@isocoeff_err = lll[1]
				end
				
				if index2.to_i == indexstart + 12
					@isocoeffzero = lll[0]
					@isocoeffzero_err = lll[1]
				end
				
				if index2.to_i == indexstart + 13
					@tstart = lll[0]
					@tstop = lll[1]
					
					@timestart_utc = @tstart
					@timestop_utc = @tstop
					if indexstart.to_i >= 16
						@timestart_tt = lll[2]
						@timestop_tt = lll[3]
						@timestart_mjd = lll[4]
						@timestop_mjd = lll[5]
					else
						@timestart_tt = datautils.time_utc_to_tt(@tstart);
						@timestop_tt = datautils.time_utc_to_tt(@tstop);
						@timestart_mjd = datautils.time_tt_to_mjd(@timestart_tt);
						@timestop_mjd = datautils.time_tt_to_mjd(@timestop_tt);
					end
					
					#calcolo fase orbitale
					@orbitalphase = -1;
					if(@calcorbitalphase_period.to_f != 0)
						timemjd = @timestart_mjd.to_f + (@timestop_mjd.to_f-@timestart_mjd.to_f)
						@orbitalphase = (timemjd.to_f - @calcorbitalphase_t0.to_f) / @calcorbitalphase_period.to_f;
						@orbitalphase = @orbitalphase.to_f - @orbitalphase.to_i;
					end
					
				end
				
				if index2.to_i == indexstart + 14
					@energyrange = lll[0]
					@fovrange = lll[1]
					@albedo = lll[2]
					@binsize  = lll[3] 
					@expstep  = lll[4]  
					@phasecode = lll[5]  
				end
				
				if index2.to_i == indexstart + 15
					#[-1 step skipped, 0 ok, 1 errors]
					@mlestep_res = line
					@mlestep_1_ext1_res = lll[0].to_i
					@mlestep_2_step1_res = lll[1].to_i
					@mlestep_3_ext2_res = lll[2].to_i
					@mlestep_4_step2_res = lll[3].to_i
					@mlestep_5_contour_res = lll[4].to_i
					@mlestep_6_index_res = lll[5].to_i
					@mlestep_7_ul_res = lll[6].to_i
				end
				
				if index2.to_i == indexstart + 16
					@mlestep_cts = line
					@mlestep_1_ext1_cts = lll[0].to_i
					@mlestep_2_step1_cts = lll[1].to_i
					@mlestep_3_ext2_cts = lll[2].to_i
					@mlestep_4_step2_cts = lll[3].to_i
					@mlestep_5_contour_cts = lll[4].to_i
					@mlestep_6_index_cts = lll[5].to_i
					@mlestep_7_ul_cts = lll[6].to_i
				end
			end

	end

	def multiOutputLine()
		@multiOutputLine = @label.to_s  + " " + @sqrtTS.to_s + " " + @l_peak.to_s + " " + @b_peak.to_s + " " + @counts.to_s + " " + @counts_error.to_s + " " + @flux.to_s + " " + @flux_error.to_s + " " + @sicalc.to_s + " " + @sicalc_error.to_s
	end
	
	def multiOutputLineFull()
		@multiOutputLineFull = @label.to_s  + " " + @sqrtTS.to_s + " " + @l_peak.to_s + " " + @b_peak.to_s + " " + @counts.to_s + " " + @counts_error.to_s + " " + @flux.to_s + " " + @flux_error.to_s + " " + @sicalc.to_s + " " + @sicalc_error.to_s + " " + @l.to_s + " " + @b.to_s + " " + @r.to_s + " " + @ell_a.to_s +  " " + @ell_b.to_s + " " + @ell_phi.to_s
	end
	
	def multiOutputLineFull2()
		@multiOutputLineFull2 = @label.to_s  + " " + @sqrtTS.to_s + " " + @l_peak.to_s + " " + @b_peak.to_s + " " + @dist.to_s + " " + @counts.to_s + " " + @counts_error.to_s + " " + @counts_ul.to_s + " " + @flux.to_s + " " + @flux_error.to_s + " " + @flux_ul.to_s + " " + @sicalc.to_s + " " + @sicalc_error.to_s + " " + @l.to_s + " " + @b.to_s + " " + @r.to_s + " " + @ell_a.to_s +  " " + @ell_b.to_s + " " + @ell_phi.to_s + " " + @galcoeff + " " + @galcoeff_err + " " + @isocoeff + " " + @isocoeff_err + " " + @galcoeffzero + " " + @galcoeffzero_err + " " + @isocoeffzero + " " + @isocoeffzero_err
		
	end
	
	def multiOutputLineFull3(flag)
		@multiOutputLineFull3 = flag + " " + @label.to_s + " " + format("%.2f", @sqrtTS) + " POS " + @l_peak.to_s + " " + @b_peak.to_s + " " + @dist.to_s + " " + @l.to_s + " " + @b.to_s + " " + @distellipse.to_s + " " + @r.to_s + " " + @ell_a.to_s +  " " + @ell_b.to_s + " " + @ell_phi.to_s + " CTS " + @counts.to_s + " " + @counts_error.to_s + " " + @counts_error_p.to_s + " " + @counts_error_m.to_s + " " + @counts_ul.to_s + " FL " + @flux.to_s + " " + @flux_error.to_s + " " + @flux_error_p.to_s + " " + @flux_error_m.to_s + " " + @flux_ul.to_s + " EXP " + @exposure.to_s + " SI " + @sicalc.to_s + " " + @sicalc_error.to_s + " TI " + @timestart_utc.to_s + " " + @timestop_utc.to_s + " " + @timestart_mjd.to_s + " " + @timestop_mjd.to_s + " " + @timestart_tt.to_s + " " + @timestop_tt.to_s + " GI " + @galcoeffzero + " " + @galcoeffzero_err + " " + @galcoeff + " " + @galcoeff_err + " "  + @isocoeffzero + " " + @isocoeffzero_err + " "  + @isocoeff + " " + @isocoeff_err + " FIT " + @fit_cts + " "  + @fit_fcn0 + " " + @fit_fcn1 + " " + @fit_edm0 + " " + @fit_edm1 + " " + @fit_iter0 + " " + @fit_iter1  + " ANA " + @fix.to_s + " " + @si_start.to_s + " " + @ulconflevel.to_s + " " + @srcconflevel.to_s + " " + @startL.to_s + " " + @startB.to_s + " " + @startFlux.to_s + " [ " + @lmin.to_s + " , " + @lmax.to_s + " ] [ " + @bmin.to_s + " , " + @bmax.to_s + " ] " + " " + @energyrange + " " + @fovrange + " " + @albedo + " " + @binsize + " " + @expstep + " " + @phasecode + " ORBPH " + format("%.3f", @orbitalphase);
		
	end
	
	def multiOutputLineFull4(flag)
		multiOutputLineFull3(flag);
		@multiOutputLineFull4 = @multiOutputLineFull3.chomp
		@multiOutputLineFull4 = @multiOutputLineFull4 + " " + @mlestep_res.chomp + " " + @mlestep_cts.chomp
		
	end
	
	def isnil(value)
		if value == nil
			@isnil = " "
		else
			@isnil = value.to_s
		end
	end
	
	def assoccat(endof)
		datautils = DataUtils.new
		@assoc = ""
		dircat = Dir[ENV["AGILE"] + "/catalogs/*.cat"]; 
		dircat.sort.each do | cat |
			
			File.open(cat).each_line do |x|
	
				a = x.split(" ");
				if a[2] == nil
					next
				end
				d = datautils.distance(a[0], a[1], @l_peak, @b_peak);
				
				if @r.to_f > 0 
					r = @r.to_f
				else 
					r = 1
				end
				if d.to_f < r.to_f	
					#puts a[0] + " " +  a[1] + " " +  @l_peak.to_s + " " + @b_peak.to_s + " " + a[2] + " " + d.to_s + " " + r.to_s				
					@assoc = @assoc + a[2] + endof
					#puts @assoc
				end
			end
		end
	end
	
	def assoc
		@assoc
	end
	
	def multiOutputLineFull3HTML(name, flag)
		
		assoccat("<br>")
		
		@multiOutputLineFull3HTML = "<tr><td>" + @label.to_s + "</td><td>" + name.to_s + "</td><td>" + format("%.2f", @sqrtTS) + "</td><td>" + format("%.2E", @flux) + " +/- " + format("%.2E", @flux_error) + "</td><td>" + format("%.2E", @flux_ul) + "</td><td>"  + @l_peak.to_s + "</td><td>" + @b_peak.to_s + "</td><td>" + format("%.2f", @dist) + "</td><td>" + isnil(@l).to_s + "</td><td>" + isnil(@b).to_s + "</td><td>" + format("%.2f", isnil(@distellipse)) + "</td><td>" + format("%.2f", isnil(@r)) + "</td><td>" + format("%.2f", isnil(@ell_a)) +  "</td><td>" + format("%.2f", isnil(@ell_b)) + "</td><td>" + format("%.2f", isnil(@ell_phi)) + "</td><td>" + format("%.2f", @counts) + " +/- " + format("%.2f", @counts_error) + "</td><td>" + format("%.2f", @counts_ul) + "</td><td>" + @assoc.to_s + "</td><td>" + format("%.3E", @exposure) + "</td><td>" + @sicalc.to_s + " +/- " + @sicalc_error.to_s + "</td><td>" + @timestart_utc.to_s + "</td><td>" + @timestop_utc.to_s + "</td><td>" + format("%.2f", @timestart_mjd) + "</td><td>" + format("%.2f", @timestop_mjd) + "</td><td>" + @timestart_tt.to_s + "</td><td>" + @timestop_tt.to_s + "</td><td>" + @galcoeffzero.to_s + "</td><td>" + @galcoeffzero_err.to_s + "</td><td>" + @galcoeff.to_s + "</td><td>" + @galcoeff_err.to_s + "</td><td>"  + @isocoeffzero.to_s + "</td><td>" + @isocoeffzero_err.to_s + "</td><td>"  + @isocoeff.to_s + "</td><td>" + @isocoeff_err.to_s + "</td><td>" + @fit_cts.to_s + "</td><td>"  + @fit_fcn0.to_s + "</td><td>" + @fit_fcn1.to_s + "</td><td>" + @fit_edm0.to_s + "</td><td>" + @fit_edm1.to_s + "</td><td>" + @fit_iter0.to_s + "</td><td>" + @fit_iter1.to_s  + "</td><td>" + @fix.to_s + "</td><td>" + @si_start.to_s + "</td><td>" + @ulconflevel.to_s + "</td><td>" + @srcconflevel.to_s + "</td><td>" + @startL.to_s + "</td><td>" + @startB.to_s + "</td><td>" + @startFlux.to_s + "</td><td>[ " + @lmin.to_s + " , " + @lmax.to_s + " ] [ " + @bmin.to_s + " , " + @bmax.to_s + " ] " + "</td><td>" + @energyrange.to_s + "</td><td>" + @fovrange.to_s + "</td><td>" + @albedo.to_s + "</td><td>" + @binsize.to_s + "</td><td>" + @expstep.to_s + "</td><td>" + @phasecode.to_s + "</td><td>" + @mlestep_res.to_s + "</td><td>" + @mlestep_cts.to_s + "</td></tr>";
		
	end
	
	def multiOutputLineFull3HTMLheader(flag)
		
		@multiOutputLineFull3HTMLheader = "<table border=1><tr><th>Name</th><th>dir</th><th>sqrt(TS)</th><th>Flux</th><th>Flux UL</th><th>l peak</th><th>b peak</th><th>dist</th><th>l</th><th>b</th><th>dist ell</th><th>R</th><th>a</th><th>b</th><th>phi</th><th>Counts</th><th>Counts UL</th><th>Assoc</th><th>Exp</th><th>SI</th><th>tstart UTC</th><th>tstop UTC</th><th>tstart MJD</th><th>tstop MJD</th><th>tstart TT</th><th>tstop TT</th><th>galcoeffzero</th><th>galcoeffzero_err</th><th>galcoeff</th><th>galcoeff_err</th><th>isocoeffzero</th><th>isocoeffzero_err</th><th>isocoeff</th><th>isocoeff_err</th><th>fit_cts</th><th>fit_fcn0</th><th>fit_fcn1</th><th>fit_edm0</th><th>fit_edm1</th><th>fit_iter0</th><th>fit_iter1</th><th>fix.to_s</th><th>si_start</th><th>ulconflevel</th><th>srcconflevel</th><th>startL</th><th>startB</th><th>startFlux</th><th>posrange</th><th>energyrange</th><th>fovrange</th><th>albedo</th><th>binsize</th><th>expstep</th><th>phasecode</th><th>MLE phase res</th><th>MLE phase cts</th></tr>";
		
	end
	
	#def multiOutputLineFull4(flag, ring, dist)
	#	multiOutputLineFull3(flag)
	#	@multiOutputLineFull4 = @multiOutputLineFull3 + " RING " + ring.to_s + " " + dist.to_s;
	#end
	
	def multiOutputLineShort3(flag)
		@multiOutputLineShort3 = flag + " " + @label.to_s + " " + format("%.2f", @sqrtTS) + " POS " + @l_peak.to_s + " " + @b_peak.to_s + " " + @dist.to_s + " " + @l.to_s + " " + @b.to_s + " " + @distellipse.to_s + " " + @r.to_s + " " + @ell_a.to_s +  " " + @ell_b.to_s + " " + @ell_phi.to_s + " FL " + @flux.to_s + " " + @flux_error.to_s + " " + @flux_ul.to_s + " EXP " + @exposure.to_s + " SI " + @sicalc.to_s + " " + @sicalc_error.to_s + " TI " + @timestart_mjd.to_s + " " + @timestop_mjd.to_s +  " CTS " + @counts.to_s + " " + @counts_error.to_s + " GI " + @galcoeff + " " + @galcoeff_err + " "  + @isocoeff + " " + @isocoeff_err + " FIT " + @fit_cts + " "  + @fit_fcn0 + " " + @fit_fcn1 + " " + @fit_edm0 + " " + @fit_edm1 + " " + @fit_iter0 + " " + @fit_iter1  + " ORBPH " + format("%.3f", @orbitalphase);
	end
	
	def regline(thr)
		l = @l
		b = @b
		ell_a = @ell_a
		ell_b = @ell_b
		ell_phi = @ell_phi
		check = ""
		if @l.to_f == -1
			l = @l_peak
			b = @b_peak
			ell_a = 1
			ell_b = 1
			ell_phi = 0
			check = "*"
		end
		
		if @sqrtTS.to_f > thr.to_f
			@regline = "galactic\nellipse("+l.to_s+","+b.to_s+","+ell_a.to_s+","+ell_b.to_s+",-"+ell_phi.to_f.abs.to_s+") #color=white width=2 text={"+@label.to_s+" "+format("%.2f",@sqrtTS.to_f)+""+check+" ("+format("%.2f,%.2f,%.2e", l, b, @flux)+")}\n"
		else
			@regline = "galactic\nellipse("+l.to_s+","+b.to_s+","+ell_a.to_s+","+ell_b.to_s+",-"+ell_phi.to_f.abs.to_s+") #color=white width=1\n"
		end
	end
	
	def multiOutputLineShort4(flag, ring, dist)
		multiOutputLineShort3(flag)
		@multiOutputLineShort4 = @multiOutputLineShort3 + " RING " + ring.to_s + " " + dist.to_s;
	end
	
	def mlestep_res
		@mlestep_res
	end
	
	def mlestep_1_ext1_res
		@mlestep_1_ext1_res
	end
	
	def mlestep_2_step1_res
		@mlestep_2_step1_res
	end
	
	def mlestep_3_ext2_res
		@mlestep_3_ext2_res
	end
	
	def mlestep_4_step2_res
		@mlestep_4_step2_res
	end
	
	def mlestep_5_contour_res
		@mlestep_5_contour_res
	end
	
	def mlestep_6_index_res
		@mlestep_6_index_res
	end
	
	def mlestep_7_ul_res
		@mlestep_7_ul_res
	end
	
	def mlestep_cts
		@mlestep_cts
	end
	
	def mlestep_1_ext1_cts
		@mlestep_1_ext1_cts
	end
	
	def mlestep_2_step1_cts
		@mlestep_2_step1_cts
	end
	
	def mlestep_3_ext2_cts
		@mlestep_3_ext2_cts
	end
	
	def mlestep_4_step2_cts
		@mlestep_4_step2_cts
	end
	
	def mlestep_5_contour_cts
		@mlestep_5_contour_cts
	end
	
	def mlestep_6_index_cts
		@mlestep_6_index_cts
	end
	
	def mlestep_7_ul_cts
		@mlestep_7_ul_cts
	end
	
	def calcorbitalphase_period
		@calcorbitalphase_period
	end
	
	def calcorbitalphase_t0
		@calcorbitalphase_t0
	end
	
	def orbitalphase
		@orbitalphase
	end
	
	def	phasecode
		@phasecode
	end
	
	def	expstep
		@expstep
	end
	
	def	binsize
		@binsize
	end
	
	def	albedo
		@albedo
	end
	
	def	fovrange
		@fovrange
	end
	
	def	energyrange
		@energyrange
	end
	
	def	lmax
		@lmax
	end
	
	def	lmin
		@lmin
	end
	
	def	bmax
		@bmax
	end
	
	def	bmin
		@bmin
	end
	
	def	tstart
		@tstart
	end
	
	def tstop
		@tstop
	end
	
	def	timestart_utc
		@timestart_utc
	end
	
	def	timestop_utc
		@timestop_utc
	end
	
	def	timestart_tt
		@timestart_tt
	end
	
	def	timestop_tt
		@timestop_tt
	end
	
	def	timestart_mjd
		@timestart_mjd
	end
	
	def	timestop_mjd
		@timestop_mjd
	end
	
	def	startL
		@startL
	end
	
	def startB
		@startB
	end
	
	def startFlux
		@startFlux
	end
	
	def dist
		@dist
	end
	
	def galcoeff 
		@galcoeff
	end
	
	def distellipse
		@distellipse
	end
	
	def fit_cts
		@fit_cts
	end
	
	def fit_fcn0
		@fit_fcn0
	end

	def fit_fcn1
		@fit_fcn1
	end

	def fit_edm0
		@fit_edm0
	end

	def fit_edm1
		@fit_edm1
	end
	
	def fit_iter0
		@fit_iter0
	end

	def fit_iter1
		@fit_iter1
	end
	
	def isocoeff
		@isocoeff
	end
	
	def galcoeff_err
		@galcoeff_err
	end
	
	def isocoeff_err
		@isocoeff_err
	end
	
	def galcoeffzero
		@galcoeffzero
	end
	
	def isocoeffzero
		@isocoeffzero
	end
	
	def galcoeffzero_err
		@galcoeffzero_err
	end
	
	def isocoeffzero_err
		@isocoeffzero_err
	end
	

	def label
		@label
	end
	
	def sources
		@sources
	end

	
	def sicalc
		@sicalc
	end
	
	def sicalc_error
		@sicalc_error
	end

	def fix
		@fix
	end

	def si_start
		@si_start
	end
	def ulconflevel
		@ulconflevel
	end
	def srcconflevel
		@srcconflevel
	end
	def sqrtTS
		if @sqrtTS == "-nan"
			@sqrtTS = 0;
		end
		@sqrtTS
	end
	def l_peak
		@l_peak
	end
	def b_peak
		@b_peak
	end
	def l
		@l
	end
	def b
		@b
	end
	def r
		@r
	end
	def ell_a
		@ell_a
	end
	def ell_b
		@ell_b
	end
	def ell_phi
		@ell_phi
	end
	def counts
		@counts
	end
	def counts_error
		@counts_error
	end
	def counts_error_p
		@counts_error_p
	end
	def counts_error_m
		@counts_error_m
	end
	def counts_ul
		@counts_ul
	end
	def flux
		@flux
	end
	def flux_error
		@flux_error
	end
	def flux_error_p
		@flux_error_p
	end
	def flux_error_m
		@flux_error_m
	end
	def flux_ul
		@flux_ul
	end
	def exposure
		@exposure
	end
	def fullellipseline
		@fullellipseline
	end
end


class MultiOutputList
	def readSources(resname, multilist, flag)
		f = File.new(resname + ".resfull", "w")
		f1 = File.new(resname + ".resfullsel", "w")
		freg = File.new(resname + ".reg", "w")
		File.open(multilist).each_line do | line |
			multioutput = MultiOutput.new()
			name = line.split(" ")[6];
			multioutput.readDataSingleSource2(resname, name)
			f.write(multioutput.multiOutputLineFull3(flag) + "\n");
			if multioutput.fix.to_i >= 1
				f1.write(multioutput.multiOutputLineFull3(flag) + "\n");
			end
			freg.write(multioutput.regline(4));
		end
		f.close();
		f1.close();
		freg.close();
	end
	
	def readSourcesInDir(dir, resname, flag, thrsqrtts)
		f = File.new(dir + "/" + resname + ".resfull", "w")
		f1 = File.new(dir + "/" + resname + ".resfullsel", "w")
		freg = File.new(dir + "/" + resname + ".reg", "w")
		fhtml = File.new(dir + "/" + resname + ".html", "w")
		fhtmlsel = File.new(dir + "/" + resname + ".sel.html", "w")
		fob = File.new(dir + "/" + resname + ".ob", "w")
		flc = File.new(dir + "/" + resname + ".lc", "w")
		multioutput = MultiOutput.new()
		fhtml.write(multioutput.multiOutputLineFull3HTMLheader(flag))
		fhtmlsel.write(multioutput.multiOutputLineFull3HTMLheader(flag))
		#puts Dir[dir + "/*.source"]	
		Dir[dir + "/*.source"].sort.each do | name |
			#puts name
			multioutput = MultiOutput.new()
			multioutput.readDataSingleSource(name)
			#puts name
			f.write(multioutput.multiOutputLineFull4(flag) + "\n"); #3 old, 4 new
			if multioutput.fix.to_i >= 1 and multioutput.sqrtTS.to_f > thrsqrtts.to_f
				f1.write(multioutput.multiOutputLineFull4(flag) + "\n"); #3 old, 4 new
			end
			
		    #puts multioutput.multiOutputLineFull3(flag)
			freg.write(multioutput.regline(thrsqrtts.to_f));
			fhtml.write(multioutput.multiOutputLineFull3HTML(name, flag))
			if multioutput.sqrtTS.to_f > 4.0
				fhtmlsel.write(multioutput.multiOutputLineFull3HTML(name, flag))
			end
			
			#write on and lc file
			runname = name.split("_MLE")[0]
			rn = runname.split("/")
			runname = rn[rn.size - 1]
			distpoint = " -1 -1 "
			distpointcalc = -1
			if multioutput.timestop_tt.to_f < 182692800
				File.open(ENV["AGILE"] + "/scripts/AGILEOBPOINT.list").each_line do | linetime |
					ll = linetime.split(" ")
					if multioutput.timestop_tt.to_f >= ll[5].to_f and multioutput.timestop_tt.to_f <= ll[6].to_f
						distpoint = ll[1].to_s + " " + ll[2].to_s + " "
						du = DataUtils.new
						distpointcalc = du.distance(multioutput.l_peak, multioutput.b_peak, ll[1].to_f, ll[2].to_f)
						break
					end
				end
			end
			#TODO: aggiungere qui la necessità di determinare lpointing e bpointing da un file di input
			fob.write(multioutput.timestart_tt.to_s + " " + multioutput.timestop_tt.to_s + " " + runname + " " + multioutput.galcoeff.to_s + " " + multioutput.isocoeff.to_s + " " + distpoint.to_s + " \n");
			
			#write lc file
			flux = 0
			fluxerror = "0"
			fluxtype = "0"
			if multioutput.sqrtTS.to_f < 3
				flux = multioutput.flux_ul
				fluxerror = "0"
				fluxtype = "0"
			else
				flux = multioutput.flux
				fluxerror = multioutput.flux_error
				fluxtype = "1"
			end
			mjdsize = multioutput.timestop_mjd.to_f - multioutput.timestart_mjd.to_f
			mjdcenter = multioutput.timestart_mjd.to_f + mjdsize.to_f / 2.0
			#calculate distance between lpointing, bpointing to l, b. NOW is -1
			flc.write(format("%.2E",flux.to_s) + " " + format("%.2E", fluxerror.to_s)  + " " + fluxtype.to_s  + " " + format("%.7f", mjdcenter.to_s)  + " " + format("%.4f",mjdsize.to_s) + " " + format("%.2f", distpointcalc.to_s) + " " + runname + " " +  format("%.2f", multioutput.sqrtTS.to_s) + " " + format("%.2E", multioutput.exposure.to_s) + " " + multioutput.galcoeff.chomp.to_s + " " + multioutput.isocoeff.chomp.to_s + " " + format("%.7f", multioutput.timestart_mjd.to_s) + " " + format("%.7f", multioutput.timestop_mjd.to_s) + " " + format("%.2f", multioutput.timestart_tt.to_s) + " " + format("%.2f", multioutput.timestop_tt.to_s) + " " + format("%.3f", multioutput.dist.to_s) + " ( " + multioutput.fullellipseline  + " ) " + multioutput.fix + "\n");
		end
		f.close();
		f1.close();
		freg.close();
		fhtml.close();
		fhtmlsel.close();
		fob.close();
		flc.close();
	end
	
	def sources
		@sources
	end
	
	
end
