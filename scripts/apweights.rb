#!/usr/bin/env ruby

apfile = ARGV[0]

fd = File.open(apfile, "r+")

tstart = []
tstop = []
e = []
c = []
n = 0
while (line = fd.gets)
	f1, f2, f3, f4 = line.split(" ")
	tstart << f1.to_f
	tstop << f2.to_f
	e << f3.to_f
	c << f4.to_f
	if f3.to_f != 0.0
		n = n+1
	end
end
fd.close


# https://docs.google.com/document/d/1FLSHiYUU1Lcl2JFkwAMKMU0GG9GbX4Ou0jUADRYIPr8/edit?usp=sharing
# Section 2
#
# [1] Robin Corbet, LAT Light Curve Analysis: Aperture Photometry and Periodicity Searches, slides at COSPAR CBW 2010
# [2] Gehrels, 1986, ApJ, 303, 336
# [3] Kraft, Burrows, & Nousek, 1991, ApJ, 374, 344

#def_fluxw
# Variance weighting corrected for low count statistic - Weighted Power Spectrum	 [1] [2]
sum1 = 0.0
sum2 = 0.0
s = []
flux = []
def_fluxw = []
e.zip(c).each do |e_i, n_i|
	if (e_i != 0.0)
		flux_i = n_i / e_i
		s_ip = 0.5 + Math.sqrt(n_i + 0.25)
		s_im = -0.5 + Math.sqrt(n_i + 0.25)
		s_irms = Math.sqrt(s_ip*s_ip + s_im*s_im) / 2.0
		s_i = s_irms / e_i
		sum1 += flux_i / (s_i*s_i)
		sum2 += 1.0 / (s_i*s_i)
		flux << flux_i
		s << s_i
	else
		flux << 0.0
		s << 0.0
	end
end
fluxw_mean = sum1 / sum2
flux.zip(s).each do |flux_i, s_i|
	if (s_i != 0.0)
		#def_fluxw << (flux_i - fluxw_mean) / s_i
		def_fluxw << (flux_i - fluxw_mean) / (s_i*s_i)
	else
		def_fluxw << 0.0
	end
end

#AA
sum1 = 0.0
sum2 = 0.0
sum3 = 0.0
flux = []
aa_fluxw = []
e.zip(c).each do |e_i, n_i|
	if (e_i != 0.0)
		flux_i = n_i / e_i
		
		#method 1
		s_ip = 0.5 + Math.sqrt(n_i + 0.25)
		s_im = -0.5 + Math.sqrt(n_i + 0.25)
		s_irms = Math.sqrt(s_ip*s_ip + s_im*s_im) / 2.0
		s_i = s_irms / e_i
		
		
		sum1 += flux_i / (s_i*s_i)
		sum2 += 1.0 / (s_i*s_i)
		sum3 += flux_i
		flux << flux_i
	else
		flux << 0.0
	end
end
flux_mean = sum3 / n
fluxw_mean = sum1 / sum2
flux.zip(e).each do |flux_i, e_i|
	if (e_i != 0.0)
		fluxpred_i = e_i * flux_mean
	 	s_i = Math.sqrt(fluxpred_i) / e_i
		#aa_fluxw << (flux_i - fluxw_mean) / s_i
		aa_fluxw << (flux_i - fluxw_mean) / (s_i*s_i)
	else
		aa_fluxw << 0.0
	end
end



#AB
sum1 = 0.0
sum2 = 0.0
sum3 = 0.0
sum4 = 0.0
sum5 = 0.0
flux = []
ab1_fluxw = []
ab2_fluxw = []
ab3_fluxw = []
e.zip(c).each do |e_i, n_i|
	if (e_i != 0.0)
		flux_i = n_i / e_i
		s_ip = 0.5 + Math.sqrt(n_i + 0.25)
		s_im = -0.5 + Math.sqrt(n_i + 0.25)
		s_irms = Math.sqrt(s_ip*s_ip + s_im*s_im) / 2.0
		s_i = s_irms / e_i
		sum1 += flux_i / (s_i*s_i)
		sum2 += 1.0 / (s_i*s_i)
		sum3 += flux_i
		sum4 += n_i
		sum5 += e_i
		flux << flux_i
		s << s_i
	else
		flux << 0.0
		s << 0.0
	end
end
fluxw_mean1 = sum1 / sum2
fluxw_mean2 = sum3 / n
fluxw_mean3 = sum4 / sum5
flux.zip(e).each do |flux_i, e_i|
	if (e_i != 0.0)
		fluxpred_i = e_i * fluxw_mean1
	 	s_i = Math.sqrt(fluxpred_i) / e_i
		#ab1_fluxw << (flux_i - fluxw_mean1) / s_i
		ab1_fluxw << (flux_i - fluxw_mean1) / (s_i*s_i)
		
		fluxpred_i = e_i * fluxw_mean2
	 	s_i = Math.sqrt(fluxpred_i) / e_i
		#ab2_fluxw << (flux_i - fluxw_mean2) / s_i
		ab2_fluxw << (flux_i - fluxw_mean2) / (s_i*s_i)
		
		#red text method
		fluxpred_i = e_i * fluxw_mean3
	 	s_i = Math.sqrt(fluxpred_i) / e_i
		#ab3_fluxw << (flux_i - fluxw_mean3) / s_i
		ab3_fluxw << (flux_i - fluxw_mean3) / (s_i*s_i)
	else
		ab1_fluxw << 0.0
		ab2_fluxw << 0.0
		ab3_fluxw << 0.0
	end
end

#Write
str = ""
precision = 2
tstart.zip(tstop, e, c, flux, s, def_fluxw, aa_fluxw, ab1_fluxw, ab2_fluxw, ab3_fluxw).each do |tstart_i, tstop_i, e_i, c_i, flux_i, s_i, def_fluxw_i, aa_fluxw_i, ab1_fluxw_i, ab2_fluxw_i, ab3_fluxw_i|
	str = str + tstart_i.to_s + " " + tstop_i.to_s + " " + e_i.to_s + " " + c_i.to_i.to_s + " " +
		  '%.5e' % flux_i.to_s + " " +
		  '%.5e' % s_i.to_s + " " +
		  '%.5f' % def_fluxw_i + " " +
		  '%.5f' % aa_fluxw_i + " " +
		  '%.5f' % ab1_fluxw_i + " " + '%.5f' % ab2_fluxw_i + " " + '%.5f' % ab3_fluxw_i + "\n"
end
#print str

File.open(apfile + ".ap2", 'w') { |file| file.write(str) }
