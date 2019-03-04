t0=54260.974
t0agile = 108516153.600 #54260.97400000000
p=8.4474*86400

puts ARGV[0]

File.open(ARGV[0]).each_line do | line |
	ll = line.split(" ")
	t1 = ll[0].to_f
	t2 = ll[1].to_f
	delta = t2-t1
	f1 = (t1-t0agile)/p
	f1 = f1 - f1.to_i
	f2 = (t2-t0agile)/p
        f2 = f2 - f2.to_i
	puts format("%.4f",t1) + " " + format("%.4f", t2) + " " + format("%.3f",delta) + " " + format("%.5f",f1) + " " + format("%.5f",f2)
end
