#0) source name
#1) l
#2) b
#3) dist
#questo script e' usato per unire i risultati delle simulazioni quando si usa AG_multi4
#Il caso è quello galattico in cui per unire i risultati serve sia il .res sia il file della singola sorgente 

load "~/grid_scripts2/conf.rb"
load "~/grid_scripts2/MultiOutput.rb"

if ARGV[0].to_s == "help" || ARGV[0] == nil || ARGV[0] == "h"
	system("head -7 " + $0 );
	exit;
end

lock = true

sname = ARGV[0]
#"2AGLJ2033+4055"
l = ARGV[1]
b = ARGV[2]
dist = ARGV[3]

m = MultiOutput.new
d = DataUtils.new

f0 = File.new("lst_00", "w")
f1 = File.new("lst_00.iclNO", "w")
fout = File.new("lst_00.iclNO.outside", "w")
foutr = File.new("lst_00.iclA", "w")
foutr2 = File.new("lst_00.iclA_without_r", "w")
foutrB = File.new("lst_00.iclB", "w")
foutr2B = File.new("lst_00.iclB_without_r", "w")
foutrC = File.new("lst_00.iclC", "w")

Dir["*res"].each do | file |
	nameout = file.to_s + "_" + sname.to_s;
	if !File.exists?(nameout)
		next
	end
	m.readDataSingleSource(file, sname)
	l1 = m.l_peak
	b1 = m.b_peak
	ddd = d.distance(l, b, l1, b1)
	ts = m.sqrtTS.to_f * m.sqrtTS.to_f
	if m.r.to_s == "nan"
		ad = " -1 ";
	else
		ad = " " + m.r.to_s + " "
	end
	
	f0.write(file.split(".")[0].to_i.to_s + " " + m.l_peak.to_s + " " + m.b_peak.to_s + " " + ts.to_s + " " + m.flux.to_s + " 2.1 1 0  " + m.galcoeff.to_s + " " + m.isocoeff.to_s + ad + "\n")

	#1) (0)
	#per scrivere lst_00 ci metto tutte le sorgenti che cascano a dist distanza rispetto al posizionamento passato come argomento.
	#se stanno piu' lontano le metto in lst_00.outside e le metto anche in lst_00 con TS = 0
	#qui non tengo conto del contour level calcolato e del TS
	if ddd.to_f < dist.to_f  
		f1.write(file.split(".")[0].to_i.to_s + " " + m.l_peak.to_s + " " + m.b_peak.to_s + " " + ts.to_s + " " + m.flux.to_s + " 2.1 1 0  " + m.galcoeff.to_s + " " + m.isocoeff.to_s + ad + "\n")
	else
		fout.write(m.l_peak.to_s + " " + m.b_peak.to_s + " " + ts.to_s + " " + m.flux.to_s + " " + m.galcoeff.to_s + " " + m.isocoeff.to_s + ad + "\n")
		f1.write(file.split(".")[0].to_i.to_s + " " + m.l_peak.to_s + " " + m.b_peak.to_s + " 0 " + m.flux.to_s + " 2.1 1 0  " + m.galcoeff.to_s + " " + m.isocoeff.to_s + ad + "\n")
	end
	
	#gli altri files (_r) li genero solo nel caso in cui TS <9 || esista un R calcolato
	
	#2) (A) divido le sorgenti in due categorie separate - Prima si applica il criterio sul TS e poi su R
	#scrivo lst_00.insidewith_r e lst_00.without_r (quelle senza contour con TS > 9)
	#se TS < 9 oppure TS > 9 ed esiste R
	if ts.to_f <= 9 ||  (ts.to_f >= 9 && m.r.to_f > 0)
		ddd = d.distance(l, b, m.l.to_f, m.b.to_f)
		#se sta dentro il contour level, scrivolo in lst_00.insidewith_r
		if ddd.to_f < (m.r.to_f + 0.1).to_f
			foutr.write(file.split(".")[0].to_i.to_s + " " + m.l_peak.to_s + " " + m.b_peak.to_s + " " + ts.to_s + " " + m.flux.to_s + " 2.1 1 0  " + m.galcoeff.to_s + " " + m.isocoeff.to_s + ad + "\n")
		else
			#put TS = 0, because no association with original source
			#altrimenti scrivilo sempre in lst_00.insidewith_r ma con TS = 0
			foutr.write(file.split(".")[0].to_i.to_s + " " + m.l_peak.to_s + " " + m.b_peak.to_s + " 0 " + m.flux.to_s + " 2.1 1 0  " + m.galcoeff.to_s + " " + m.isocoeff.to_s + ad + "\n")
		end
	else
		#altrimenti, se TS > 9 e non è stato calcolato R, scrivi la sorgente solo in lst_00.without_r
		foutr2.write(file.split(".")[0].to_i.to_s + " " + m.l_peak.to_s + " " + m.b_peak.to_s + " " + ts.to_s + " " + m.flux.to_s + " 2.1 1 0  " + m.galcoeff.to_s + " " + m.isocoeff.to_s + ad + "\n")
	end
	
	#3) (B) divido le sorgenti in due categorie separate - prima si applica il criterio su R e poi su TS
	#se è stato calcolato il contour B oppure TS < 9, scrivi la sorgente in lst_00.insidewith_rB, Se TS >9 e outisde contour level, metti TS=0
	#Se TS > 9 e non e' stato calcolato il contour, scrivi la sorgente in lst_00.without_rB
	#La differenza rispetto al caso 2 e' che qui prima verifico se c'è il contour, se esiste lavoro su questo, altrimenti entra in campo il criterio sul TS
	if  m.r.to_f > 0
		ddd = d.distance(l, b, m.l.to_f, m.b.to_f)
		#puts ts.to_s + " " + m.l.to_s + " " + m.b.to_s + " " + m.r.to_s
		if ddd.to_f < (m.r.to_f + 0.1).to_f
			#associazione ok
			foutrB.write(file.split(".")[0].to_i.to_s + " " + m.l_peak.to_s + " " + m.b_peak.to_s + " " + ts.to_s + " " + m.flux.to_s + " 2.1 1 0  " + m.galcoeff.to_s + " " + m.isocoeff.to_s + ad + "\n")
		else
			#put TS = 0, because no association with original source
			foutrB.write(file.split(".")[0].to_i.to_s + " " + m.l_peak.to_s + " " + m.b_peak.to_s + " 0 " + m.flux.to_s + " 2.1 1 0  " + m.galcoeff.to_s + " " + m.isocoeff.to_s + ad + "\n")
		end
	else
		#non è stato calcolato il contour ma TS<9
		if m.r.to_f <= 0 &&  ts.to_f <= 9 
			foutrB.write(file.split(".")[0].to_i.to_s + " " + m.l_peak.to_s + " " + m.b_peak.to_s + " " + ts.to_s + " " + m.flux.to_s + " 2.1 1 0  " + m.galcoeff.to_s + " " + m.isocoeff.to_s + ad + "\n")
		else
			#nonè stato calcolato il contour, ma TS>9, scarta per indeterminazione del risultato
			foutr2B.write(file.split(".")[0].to_i.to_s + " " + m.l_peak.to_s + " " + m.b_peak.to_s + " " + ts.to_s + " " + m.flux.to_s + " 2.1 1 0  " + m.galcoeff.to_s + " " + m.isocoeff.to_s + ad + "\n")
		end
	end
	
	#4) (C) unico file, TS = 0 per le sorgenti con TS>9 senza contour - prima si applica il criterio su R e poi su TS
	#se è stato calcolato il contour C
	#La differenza rispetto al caso B e' che tutte le sorgenti sono scritte in lst_00.insidewith_rC. Quelle con TS>9 e senza contour hanno TS=0
	if  m.r.to_f > 0
		ddd = d.distance(l, b, m.l.to_f, m.b.to_f)
		#puts ts.to_s + " " + m.l.to_s + " " + m.b.to_s + " " + m.r.to_s
		if ddd.to_f < (m.r.to_f + 0.1).to_f
			#associazione ok
			foutrC.write(file.split(".")[0].to_i.to_s + " " + m.l_peak.to_s + " " + m.b_peak.to_s + " " + ts.to_s + " " + m.flux.to_s + " 2.1 1 0  " + m.galcoeff.to_s + " " + m.isocoeff.to_s + ad + "\n")
		else
			#put TS = 0, because no association with original source
			foutrC.write(file.split(".")[0].to_i.to_s + " " + m.l_peak.to_s + " " + m.b_peak.to_s + " 0 " + m.flux.to_s + " 2.1 1 0  " + m.galcoeff.to_s + " " + m.isocoeff.to_s + ad + "\n")
		end
	else
		#non è stato calcolato il contour ma TS<9
		if m.r.to_f <= 0 &&  ts.to_f <= 9 
			foutrC.write(file.split(".")[0].to_i.to_s + " " + m.l_peak.to_s + " " + m.b_peak.to_s + " " + ts.to_s + " " + m.flux.to_s + " 2.1 1 0  " + m.galcoeff.to_s + " " + m.isocoeff.to_s + ad + "\n")
		else
			#nonè stato calcolato il contour, ma TS>9, TS = 0
			foutrC.write(file.split(".")[0].to_i.to_s + " " + m.l_peak.to_s + " " + m.b_peak.to_s + " 0 " + m.flux.to_s + " 2.1 1 0  " + m.galcoeff.to_s + " " + m.isocoeff.to_s + ad + "\n")
		end
	end
	
	if  (ts.to_f < 1 && m.r.to_f >= 0)
		puts file.split(".")[0].to_i.to_s + " " + m.l_peak.to_s + " " + m.b_peak.to_s + " " + ts.to_s + " " + m.flux.to_s + " 2.1 1 0  " + m.galcoeff.to_s + " " + m.isocoeff.to_s + ad + "\n"
	end
end
f0.close 
f1.close
fout.close
foutr.close
foutr2.close
foutrB.close
foutrC.close
foutr2B.close