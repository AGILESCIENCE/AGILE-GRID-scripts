#! /usr/bin/ruby
#0) input file name
#1) color (optiona, red)


if ARGV[0].to_s == "help" || ARGV[0].to_s == "h" || ARGV[0] == nil
	system("head -3 " + $0 );
	exit;
end

fileinp1 = ARGV[0]
a = File.open(fileinp1)

if ARGV[1] != nil
	color = ARGV[1].to_s;
else
	color = "red";
end


	index = 0;
	lastline = "";
	c = File.new(ARGV[0] + ".reg", "w");
	a.each_line do | l |
		if l.split("!").size() >= 2
			next
		end
		if index.to_i == 0
			outline = "global color=green font=\"helvetica 9 normal\" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n"
			c.write(outline)
			outline = "galactic\n"
			c.write(outline)
			puts "\\begin {table}[!htb]"
			puts "\\caption {\\em{}}"
			puts "\\label{table_}"
			puts "\\renewcommand{\\arraystretch}{1.2}";
			puts "\\begin{tabular}{@{}llll} \\hline";
			puts "\\textbf{ l } & \\textbf{ b } & \\textbf{ Flux } & \\textbf{ Name } \\\\"
		end
# 		if index.to_i >= 2
			ee = l.split(" ");
			if ee.size() <= 2
				next
			end
			dist = ""		
			if ee[8] != nil
				dist = ee[8]
			end
			if ee[1] != "nan"
				outline = "ellipse(" + ee[1].to_s + "," + ee[2].to_s + ",0.5,0.5,0) #color=" + color.to_s + " width=1 text={" + ee[6].to_s + " - " + ee[4].to_s + " - " + ee[0].to_s + "}";
				c.write(outline);
				c.write("\n");
				puts  ee[1].to_s + " & " + ee[2].to_s + " & " + ee[0].to_s + " & " + ee[6].to_s + "\\\\";
			end
# 		end

		index = index + 1;
	end
	c.close();
	puts "\\hline \\end{tabular} \\end{table}"



