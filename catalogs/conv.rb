File.open(ARGV[0]).each_line do | line |
	l = line.split(" ")
	r = l[12]
	if r.to_f == -1
		r = "500\""
	else
		r = r.to_s
	end
	puts("ellipse(" +l[2] + "," + l[3] + "," + r + "," + r + ",0) #color=white width=1 text={" + l[0] + " " + l[1] + rs + "}")	
end
