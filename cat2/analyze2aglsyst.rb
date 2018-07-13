load ENV['AGILE'] + "/scripts/DataUtils.rb"
load ENV['AGILE'] + "/scripts/MultiOutput6.rb"
datautils = DataUtils.new
mout = MultiOutput6.new


foutall = File.new("systFB.txt", "a")
fout100 = File.new("syst100.txt", "a")
fout300 = File.new("syst300.txt", "a")
fout1000 = File.new("syst1000.txt", "a")
fout3000 = File.new("syst3000.txt", "a")

File.open("/ANALYSIS3/catalogs/cat2.multi").each do | line |

lll=line.split(" ")

l = lll[1]
b = lll[2]
name = lll[6]
puts line

mind = 99999
dirmind = ""
Dir["[0-9]*"].each do | dir |
	l1 = dir.split("_")[0]
	b1 = dir.split("_")[1]
	d = datautils.distance(l.to_f, b.to_f, l1.to_f, b1.to_f)	
	if d < mind
		mind = d
		dirmind = dir
	end
end
puts dirmind
Dir.chdir(dirmind)


system(ENV["AGILE"] + "/scripts/cat2/cat2spectra_systematic.rb " + name + " pl 00100-10000 EDP1-EB01-FB01 H0025 3 0 SYSTEMATIC 1 -1 -1 -1 0e-08")



mout.readDataSingleSource("SYS0_SYSTEMATIC_FF2_R3_C00-pl-00100-10000-EDP1-EB01-FB01-H0025_" + name + ".source")
lout = "0 " + mout.label.to_s + " " + mout.l_peak.to_s + " " + mout.b_peak.to_s + " " + mout.sqrtTS.to_s + " " + mout.flux.to_s + " " + mout.flux_error.to_s + " " + mout.flux_ul + " " + mout.galcoeff.to_s + " " + mout.galcoeff_err.to_s + " " + mout.isocoeff.to_s + " " + mout.isocoeff_err.to_s + "\n";
foutall.write(lout)
mout.readDataSingleSource("SYS1_SYSTEMATIC_FF2_R3_C00-pl-00100-10000-EDP1-EB01-FB01-H0025_" + name + ".source")
lout = "1 " + mout.label.to_s + " " + mout.l_peak.to_s + " " + mout.b_peak.to_s + " " + mout.sqrtTS.to_s + " " + mout.flux.to_s + " " + mout.flux_error.to_s + " " + mout.flux_ul + " " + mout.galcoeff.to_s + " " + mout.galcoeff_err.to_s + " " + mout.isocoeff.to_s + " " + mout.isocoeff_err.to_s + "\n";;
foutall.write(lout)
mout.readDataSingleSource("SYS2_SYSTEMATIC_FF2_R3_C00-pl-00100-10000-EDP1-EB01-FB01-H0025_" + name + ".source")
lout = "2 " + mout.label.to_s + " " + mout.l_peak.to_s + " " + mout.b_peak.to_s + " " + mout.sqrtTS.to_s + " " + mout.flux.to_s + " " + mout.flux_error.to_s + " " + mout.flux_ul + " " + mout.galcoeff.to_s + " " + mout.galcoeff_err.to_s + " " + mout.isocoeff.to_s + " " + mout.isocoeff_err.to_s + "\n";;
foutall.write(lout)


mout.readDataSingleSource("SYS0_SYSTEMATIC_FF2_R3_C00-pl-00100-10000-EDP1-EB01-FB01-H0025_EMIN00100_EMAX00300_FM3.119_ASDCe_H0025_B01_" + name + ".source")
lout = "0 " + mout.label.to_s + " " + mout.l_peak.to_s + " " + mout.b_peak.to_s + " " + mout.sqrtTS.to_s + " " + mout.flux.to_s + " " + mout.flux_error.to_s + " " + mout.flux_ul + " " + mout.galcoeff.to_s + " " + mout.galcoeff_err.to_s + " " + mout.isocoeff.to_s + " " + mout.isocoeff_err.to_s + "\n";
fout100.write(lout)
mout.readDataSingleSource("SYS1_SYSTEMATIC_FF2_R3_C00-pl-00100-10000-EDP1-EB01-FB01-H0025_EMIN00100_EMAX00300_FM3.119_ASDCe_H0025_B01_" + name + ".source")
lout = "1 " + mout.label.to_s + " " + mout.l_peak.to_s + " " + mout.b_peak.to_s + " " + mout.sqrtTS.to_s + " " + mout.flux.to_s + " " + mout.flux_error.to_s + " " + mout.flux_ul + " " + mout.galcoeff.to_s + " " + mout.galcoeff_err.to_s + " " + mout.isocoeff.to_s + " " + mout.isocoeff_err.to_s + "\n";
fout100.write(lout)
mout.readDataSingleSource("SYS2_SYSTEMATIC_FF2_R3_C00-pl-00100-10000-EDP1-EB01-FB01-H0025_EMIN00100_EMAX00300_FM3.119_ASDCe_H0025_B01_" + name + ".source")
lout = "2 " + mout.label.to_s + " " + mout.l_peak.to_s + " " + mout.b_peak.to_s + " " + mout.sqrtTS.to_s + " " + mout.flux.to_s + " " + mout.flux_error.to_s + " " + mout.flux_ul + " " + mout.galcoeff.to_s + " " + mout.galcoeff_err.to_s + " " + mout.isocoeff.to_s + " " + mout.isocoeff_err.to_s + "\n";
fout100.write(lout)

mout.readDataSingleSource("SYS0_SYSTEMATIC_FF2_R3_C00-pl-00100-10000-EDP1-EB01-FB01-H0025_EMIN00300_EMAX01000_FM3.119_ASDCe_H0025_B01_" + name + ".source")
lout = "0 " + mout.label.to_s + " " + mout.l_peak.to_s + " " + mout.b_peak.to_s + " " + mout.sqrtTS.to_s + " " + mout.flux.to_s + " " + mout.flux_error.to_s + " " + mout.flux_ul + " " + mout.galcoeff.to_s + " " + mout.galcoeff_err.to_s + " " + mout.isocoeff.to_s + " " + mout.isocoeff_err.to_s + "\n";
fout300.write(lout)
mout.readDataSingleSource("SYS1_SYSTEMATIC_FF2_R3_C00-pl-00100-10000-EDP1-EB01-FB01-H0025_EMIN00300_EMAX01000_FM3.119_ASDCe_H0025_B01_" + name + ".source")
lout = "1 " + mout.label.to_s + " " + mout.l_peak.to_s + " " + mout.b_peak.to_s + " " + mout.sqrtTS.to_s + " " + mout.flux.to_s + " " + mout.flux_error.to_s + " " + mout.flux_ul + " " + mout.galcoeff.to_s + " " + mout.galcoeff_err.to_s + " " + mout.isocoeff.to_s + " " + mout.isocoeff_err.to_s + "\n";
fout300.write(lout)
mout.readDataSingleSource("SYS2_SYSTEMATIC_FF2_R3_C00-pl-00100-10000-EDP1-EB01-FB01-H0025_EMIN00300_EMAX01000_FM3.119_ASDCe_H0025_B01_" + name + ".source")
lout = "2 " + mout.label.to_s + " " + mout.l_peak.to_s + " " + mout.b_peak.to_s + " " + mout.sqrtTS.to_s + " " + mout.flux.to_s + " " + mout.flux_error.to_s + " " + mout.flux_ul + " " + mout.galcoeff.to_s + " " + mout.galcoeff_err.to_s + " " + mout.isocoeff.to_s + " " + mout.isocoeff_err.to_s + "\n";
fout300.write(lout)

mout.readDataSingleSource("SYS0_SYSTEMATIC_FF2_R3_C00-pl-00100-10000-EDP1-EB01-FB01-H0025_EMIN01000_EMAX03000_FM3.119_ASDCe_H0025_B01_" + name + ".source")
lout = "0 " + mout.label.to_s + " " + mout.l_peak.to_s + " " + mout.b_peak.to_s + " " + mout.sqrtTS.to_s + " " + mout.flux.to_s + " " + mout.flux_error.to_s + " " + mout.flux_ul + " " + mout.galcoeff.to_s + " " + mout.galcoeff_err.to_s + " " + mout.isocoeff.to_s + " " + mout.isocoeff_err.to_s + "\n";
fout1000.write(lout)
mout.readDataSingleSource("SYS1_SYSTEMATIC_FF2_R3_C00-pl-00100-10000-EDP1-EB01-FB01-H0025_EMIN01000_EMAX03000_FM3.119_ASDCe_H0025_B01_" + name + ".source")
lout = "1 " + mout.label.to_s + " " + mout.l_peak.to_s + " " + mout.b_peak.to_s + " " + mout.sqrtTS.to_s + " " + mout.flux.to_s + " " + mout.flux_error.to_s + " " + mout.flux_ul + " " + mout.galcoeff.to_s + " " + mout.galcoeff_err.to_s + " " + mout.isocoeff.to_s + " " + mout.isocoeff_err.to_s + "\n";;
fout1000.write(lout)
mout.readDataSingleSource("SYS2_SYSTEMATIC_FF2_R3_C00-pl-00100-10000-EDP1-EB01-FB01-H0025_EMIN01000_EMAX03000_FM3.119_ASDCe_H0025_B01_" + name + ".source")
lout = "2 " + mout.label.to_s + " " + mout.l_peak.to_s + " " + mout.b_peak.to_s + " " + mout.sqrtTS.to_s + " " + mout.flux.to_s + " " + mout.flux_error.to_s + " " + mout.flux_ul + " " + mout.galcoeff.to_s + " " + mout.galcoeff_err.to_s + " " + mout.isocoeff.to_s + " " + mout.isocoeff_err.to_s + "\n";;
fout1000.write(lout)

mout.readDataSingleSource("SYS0_SYSTEMATIC_FF2_R3_C00-pl-00100-10000-EDP1-EB01-FB01-H0025_EMIN03000_EMAX10000_FM3.119_ASDCe_H0025_B01_" + name + ".source")
lout = "0 " + mout.label.to_s + " " + mout.l_peak.to_s + " " + mout.b_peak.to_s + " " + mout.sqrtTS.to_s + " " + mout.flux.to_s + " " + mout.flux_error.to_s + " " + mout.flux_ul + " " + mout.galcoeff.to_s + " " + mout.galcoeff_err.to_s + " " + mout.isocoeff.to_s + " " + mout.isocoeff_err.to_s + "\n";;
fout3000.write(lout)
mout.readDataSingleSource("SYS1_SYSTEMATIC_FF2_R3_C00-pl-00100-10000-EDP1-EB01-FB01-H0025_EMIN03000_EMAX10000_FM3.119_ASDCe_H0025_B01_" + name + ".source")
lout = "1 " + mout.label.to_s + " " + mout.l_peak.to_s + " " + mout.b_peak.to_s + " " + mout.sqrtTS.to_s + " " + mout.flux.to_s + " " + mout.flux_error.to_s + " " + mout.flux_ul + " " + mout.galcoeff.to_s + " " + mout.galcoeff_err.to_s + " " + mout.isocoeff.to_s + " " + mout.isocoeff_err.to_s + "\n";
fout3000.write(lout)
mout.readDataSingleSource("SYS2_SYSTEMATIC_FF2_R3_C00-pl-00100-10000-EDP1-EB01-FB01-H0025_EMIN03000_EMAX10000_FM3.119_ASDCe_H0025_B01_" + name + ".source")
lout = "2 " + mout.label.to_s + " " + mout.l_peak.to_s + " " + mout.b_peak.to_s + " " + mout.sqrtTS.to_s + " " + mout.flux.to_s + " " + mout.flux_error.to_s + " " + mout.flux_ul + " " + mout.galcoeff.to_s + " " + mout.galcoeff_err.to_s + " " + mout.isocoeff.to_s + " " + mout.isocoeff_err.to_s + "\n";
fout3000.write(lout)


Dir.chdir("..")

end

