#0) source name
#1) l
#2) b
#3) dist
#questo script e' usato per unire i risultati delle simulazioni quando si usa AG_multi4
#Il caso è quello galattico e per unire i risultati basta il .res e non anche il file singolo della sorgente 

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

f1 = File.new("lst_00", "w")
fout = File.new("lst_00.outside", "w")

Dir["*.res"].each do | file |
	m.readSingleSourceRes(file, sname)
	l1 = m.l
	b1 = m.b
	ddd = d.distance(l, b, l1, b1)
	ts = m.sqrtTS.to_f * m.sqrtTS.to_f
	if ddd.to_f < dist.to_f 
		f1.write(file.split(".")[0].to_i.to_s + " " + m.l.to_s + " " + m.b.to_s + " " + ts.to_s + " " + m.flux.to_s + " 2.1 1 0  " + m.galcoeff.to_s + " " + m.isocoeff.to_s + "\n")
	else
		fout.write(m.l.to_s + " " + m.b.to_s + " " + ts.to_s + " " + m.flux.to_s + " " + m.galcoeff.to_s + " " + m.isocoeff.to_s + "\n")
		f1.write(file.split(".")[0].to_i.to_s + " " + m.l.to_s + " " + m.b.to_s + " 0 " + m.flux.to_s + " 2.1 1 0  " + m.galcoeff.to_s + " " + m.isocoeff.to_s + "\n")
	end
end

f1.close
fout.close
