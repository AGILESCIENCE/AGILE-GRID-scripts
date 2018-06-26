#0) start file name (e.g. "lst_00")
#1) l
#2) b
#3) distuser per iclNO
#4) input file format (1 = log multims4/5 B22, 5 = log multisim4 gtb3)

#Input  NUMIT L B TS FLUX SPECTRAL_INDEX FIXFLAG MINTS GAL ISO R ...
#Output file format is 1 but a mistake on gal and iso parameters: NUMIT L B TS FLUX SPECTRAL_INDEX FIXFLAG MINTS GAL ISO R


load "~/grid_scripts2/conf.rb"
load "~/grid_scripts2/MultiOutput.rb"

if ARGV[0].to_s == "help" || ARGV[0] == nil || ARGV[0] == "h"
	system("head -8 " + $0 );
	exit;
end

lock = true

#2AGLJ2033+4055 79.9031 0.558037
l = ARGV[1]
b = ARGV[2]
distuser = ARGV[3]
inputtype = ARGV[4]

m = MultiOutput.new
d = DataUtils.new

foutrC = File.new("lst_00.iclCALC", "w")
foutrC2 = File.new("lst_00.iclCALCMIX", "w")
foutrR2 = File.new("lst_00.iclCALC_without_r", "w")
foutrC1 = File.new("lst_00.iclC1", "w")
foutrC3 = File.new("lst_00.iclC2", "w")
foutrCC = File.new("lst_00.iclC", "w")
f1 = File.new("lst_00.iclNO", "w")
fout = File.new("lst_00.iclNO.outside", "w")

#SOURCE:L:B:TS:FLUX:SPECTRALINDEX:FIXFLAG:MINTS:GAL:ISO:R

File.open(ARGV[0]).each_line do | line  |

	#puts line

	ll = line.split(" ")
	
	l1 = ll[1]
	b1 = ll[2]
	dist  = d.distance(l, b, l1, b1)
	ts = ll[3]
	
	if inputtype.to_i == 1
		gal = ll[8]
		if gal == nil
			gal = -1
		end
		iso = ll[9]
		if iso == nil
			iso = -1
		end
		
		r = ll[10]
		if r == nil 
			r = -1
		end
	end
	
	if inputtype.to_i == 5
		gal = nil
		if gal == nil
			gal = -1
		end
		iso = nil
		if iso == nil
			iso = -1
		end
		
		r = ll[8]
		if r == nil 
			r = -1
		end
		gal = r; #questo riporta il formato file al 5
	end	
	
	
	#iclNO #######################################
	#per scrivere lst_00.iclNO ci metto tutte le sorgenti che cascano a dist distanza rispetto al posizionamento passato come argomento.
	#se stanno piu' lontano le metto in lst_00.iclNO.outside e le metto anche in lst_00.iclNO con TS = 0
	#qui non tengo conto del contour level calcolato e del TS
	if dist.to_f < distuser.to_f
		f1.write(ll[0]+ " " + ll[1] + " " + ll[2] + " " + ll[3] + " " + ll[4] + " " + ll[5] + " " + ll[6] + " " + ll[7] + " " + gal.to_s + " " + iso.to_s + " " + r.to_s + "\n")
	else
		f1.write(ll[0]+ " " + ll[1] + " " + ll[2] + " 0 " + ll[4] + " " + ll[5] + " " + ll[6] + " " + ll[7] + " " + gal.to_s + " " + iso.to_s + " " + r.to_s + "\n")
		fout.write(ll[0]+ " " + ll[1] + " " + ll[2] + " " + ll[3] + " " + ll[4] + " " + ll[5] + " " + ll[6] + " " + ll[7] + " " + gal.to_s + " " + iso.to_s + " " +  r.to_s + "\n")
	end
	
	#C1 ##########################################
	# unico file, TS = 0 per le sorgenti con TS>9 senza contour - prima si applica il criterio su R e poi su TS
	#se è stato calcolato il contour C
	#La differenza rispetto al caso B e' che tutte le sorgenti sono scritte in lst_00.iclC1. Quelle con TS>9 e senza contour hanno TS=0

	if r.to_f > 0
		if dist.to_f < r.to_f
			foutrC1.write(ll[0]+ " " + ll[1] + " " + ll[2] + " " + ll[3] + " " + ll[4] + " " + ll[5] + " " + ll[6] + " " + ll[7] + " " + gal.to_s + " " + iso.to_s + " " + r.to_s + "\n")
		else
			foutrC1.write(ll[0]+ " " + ll[1] + " " + ll[2] + " 0 " + ll[4] + " " + ll[5] + " " + ll[6] + " " + ll[7] + " " + gal.to_s + " " + iso.to_s + " " + r.to_s + "\n")
		end
	else
		if r.to_f < 0 and ts.to_f <= 9
			foutrC1.write(ll[0]+ " " + ll[1] + " " + ll[2] + " " + ll[3] + " " + ll[4] + " " + ll[5] + " " + ll[6] + " " + ll[7] + " " + gal.to_s + " " + iso.to_s + " " + r.to_s + "\n")
		else
			foutrC1.write(ll[0]+ " " + ll[1] + " " + ll[2] + " 0 " + ll[4] + " " + ll[5] + " " + ll[6] + " " + ll[7] + " " + gal.to_s + " " + iso.to_s + " " + r.to_s + "\n")
		end
	end


	#C ##########################################
	#come C1, ma tengo conto di 0.1 gradi del sistematico
	if r.to_f > 0
		if dist.to_f < (r.to_f + 0.1)
			foutrCC.write(ll[0]+ " " + ll[1] + " " + ll[2] + " " + ll[3] + " " + ll[4] + " " + ll[5] + " " + ll[6] + " " + ll[7] + " " + gal.to_s + " " + iso.to_s + " " + r.to_s + "\n")
		else
			foutrCC.write(ll[0]+ " " + ll[1] + " " + ll[2] + " 0 " + ll[4] + " " + ll[5] + " " + ll[6] + " " + ll[7] + " " + gal.to_s + " " + iso.to_s + " " + r.to_s + "\n")
		end
	else
		if r.to_f < 0 and ts.to_f <= 9
			foutrCC.write(ll[0]+ " " + ll[1] + " " + ll[2] + " " + ll[3] + " " + ll[4] + " " + ll[5] + " " + ll[6] + " " + ll[7] + " " + gal.to_s + " " + iso.to_s + " " + r.to_s + "\n")
		else
			foutrCC.write(ll[0]+ " " + ll[1] + " " + ll[2] + " 0 " + ll[4] + " " + ll[5] + " " + ll[6] + " " + ll[7] + " " + gal.to_s + " " + iso.to_s + " " + r.to_s + "\n")
		end
	end
	
	#CALCMIX ################################################
	#stimo il raggio di errore e uso questo invece di quello calcolato dalla multi se la multi non l'ha calcolato
	if r.to_f == -1
		if ts.to_f < 29
			#r = -0.1 * ts.to_f + 2.5
			r = -0.09 * ts.to_f + 3
		else
			r = 0.5
		end
	end
	r = r.to_f
	#0to95: r = r.to_f * 0.7
	#puts r.to_s + " " + dist.to_s

	if r.to_f >= dist.to_f  
		foutrC2.write(ll[0]+ " " + ll[1] + " " + ll[2] + " " + ll[3] + " " + ll[4] + " " + ll[5] + " " + ll[6] + " " + ll[7] + " " + gal.to_s + " " + iso.to_s + " " + r.to_s + "\n")
	else
		foutrC2.write(ll[0]+ " " + ll[1] + " " + ll[2] + " 0 " + ll[4] + " " + ll[5] + " " + ll[6] + " " + ll[7] + " " + gal.to_s + " " + iso.to_s + " " + r.to_s + "\n")
	end
	
	r=ll[10]
	if r.to_f == -1 and ts.to_f > 3
		foutrR2.write(ll[0]+ " " + ll[1] + " " + ll[2] + " " + ll[3] + " " + ll[4] + " " + ll[5] + " " + ll[6] + " " + ll[7] + " " + gal.to_s + " " + iso.to_s + " " + r.to_s + "\n")
	end
	
	#CALC ####################################################
	#stimo il raggio di errore senza mai utilizzare quello della multi
	r=0;
	if ts.to_f < 29
		#r = -0.1 * ts.to_f + 2.5
		r = -0.09 * ts.to_f + 3
	else
		r = 0.5
	end
	if r.to_f >= dist.to_f  
		foutrC.write(ll[0]+ " " + ll[1] + " " + ll[2] + " " + ll[3] + " " + ll[4] + " " + ll[5] + " " + ll[6] + " " + ll[7] + " " + gal.to_s + " " + iso.to_s + " " + r.to_s + "\n")
	else
		foutrC.write(ll[0]+ " " + ll[1] + " " + ll[2] + " 0 " + ll[4] + " " + ll[5] + " " + ll[6] + " " + ll[7] + " " + gal.to_s + " " + iso.to_s + " " + r.to_s + "\n")
	end
	
	
	#C2 ##########################################
	# unico file, TS = 0 per le sorgenti con TS>6 senza contour - prima si applica il criterio su R e poi su TS
	#se è stato calcolato il contour C
	#La differenza rispetto al caso B e' che tutte le sorgenti sono scritte in lst_00.iclC1. Quelle con TS>9 e senza contour hanno TS=0

	if r.to_f > 0 #c'è il contour
		if dist.to_f < r.to_f #la sorgente sta dentro il contour
			foutrC3.write(ll[0]+ " " + ll[1] + " " + ll[2] + " " + ll[3] + " " + ll[4] + " " + ll[5] + " " + ll[6] + " " + ll[7] + " " + gal.to_s + " " + iso.to_s + " " + r.to_s + "\n")
		else #altrimento sta fuori
			foutrC3.write(ll[0]+ " " + ll[1] + " " + ll[2] + " 0 " + ll[4] + " " + ll[5] + " " + ll[6] + " " + ll[7] + " " + gal.to_s + " " + iso.to_s + " " + r.to_s + "\n")
		end
	else #non c'è il contour
		if ts.to_f <= 6
			foutrC3.write(ll[0]+ " " + ll[1] + " " + ll[2] + " " + ll[3] + " " + ll[4] + " " + ll[5] + " " + ll[6] + " " + ll[7] + " " + gal.to_s + " " + iso.to_s + " " + r.to_s + "\n")
		else
			foutrC3.write(ll[0]+ " " + ll[1] + " " + ll[2] + " 0 " + ll[4] + " " + ll[5] + " " + ll[6] + " " + ll[7] + " " + gal.to_s + " " + iso.to_s + " " + r.to_s + "\n")
		end
	end
	
	
end

foutrC.close
foutrC2.close
foutrR2.close
foutrC1.close
foutrC3.close
foutrCC.close
f1.close
fout.close
