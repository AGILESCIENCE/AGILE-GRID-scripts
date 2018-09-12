
load ENV["AGILE"] + "/AGILEPIPE/env.rb"

abspath=PATH_RES + "/commands/"

Dir[abspath + "*.genimg"].each do | file |
	index = 0
	File.open(file).each_line do | line |
		if index.to_i == 0
			Dir.cd(line)
		end
		if index.to_i == 1
			cmd = "xvfb-run -d -s \"-screen 0 2000x2000x24\""
			cmd += line
			puts cmd
			system(cmd)
		end
		index = index + 1
	end
	
end

