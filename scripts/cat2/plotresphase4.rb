load ENV["AGILE"] + "/scripts/conf.rb"

datautils = DataUtils.new
mo = MultiOutput.new

mo.readDataSingleSource(ARGV[0]);
l1 = mo.l_peak;
b1 = mo.b_peak;

min = 9999
dirmin = ""

	
Dir["/scratch/bulgarelli/CAT2/CATALOG128_ASDCe_B01_I0023_FOVBINUMBER1_ENERGYBIN0/0*.*"].each do | dirs |
		dd = dirs.split("/")[5].to_s
		lc = dd.split(".")[0].to_f
		bc = dd.split(".")[1].to_f
		
		d = datautils.distance(lc, bc, l1, b1)
		#		puts l1.to_s + " " + b1.to_s + " " + lc.to_s + " " + bc.to_s + " " + d.to_s
		if d.to_f < min.to_f
			min = d
			dirmin = dirs
		end
end

regcat2 = ARGV[0].to_s + ".reg"
if File.exists?(regcat2) == false
	fout = File.new(regcat2, "w")
	fout.write(mo.regline(0) + "\n")
	fout.close()
end

#nohup Xvfb :1 -screen 0 2024x2024x16 &
c="export DISPLAY=localhost:1.0; "
cmd = "export DISPLAY=localhost:1.0; ds9 " + dirmin.to_s + "/FM3.119_ASDCe_I0023_B01.cts.gz  -zoom 8 -pan to " + l1.to_s + " " + b1.to_s + " wcs galactic -smooth radius 3 -smooth yes -cmap B -scale squared  -region " + regcat2 + " -region " + ENV["AGILE"] + "/catalogs/3FGL/3FGL_gll_psc_v14_ell.reg -region " +  ENV["AGILE"] +"catalogs/3FGL/3FGL_gll_psc_v14_assoc.reg -region " + ENV["AGILE"] + "/catalogs/cat2_res2_3_dash.reg -region " + ENV["AGILE"] + "/catalogs/phase3sel.reg -grid yes -geometry 1024x1024 -saveimage png " + ARGV[0] + "_fermiassoc.png -exit"

system(cmd)

puts cmd
