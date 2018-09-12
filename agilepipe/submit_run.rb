#! /usr/bin/ruby

# Copyright (c) 2016, AGILE team
# Authors: Andrea Bulgarelli <bulgarelli@iasfbo.inaf.it>,
#          Andrea Zoli <zoli@iasfbo.inaf.it>
#
# Any information contained in this software is property of the AGILE TEAM
# and is strictly private and confidential. All rights reserved.

load ENV["AGILE"] + "/AGILEPIPE/env.rb"
load ENV["AGILE"] + "/AGILEPIPE/spot6/Conf.rb"

$root = ENV['PATH_RES']

def clusteranalysis(filenameconf)
	agilepipe_path = ENV["AGILE"] + "/AGILEPIPE/"

	analysistype = nil
	instrument = nil
	triggerid = nil
	queue = nil
	alert_name = nil
	path_3908= nil
	path_3916 = nil
	orbit_number = nil

	f = File.open filenameconf
	analysistype = f.gets.chomp

	if analysistype == "GRB" or analysistype == "GW" or analysistype == "visCheckAux" or analysistype.include? "MCAL" or analysistype.include? "Ratemeters"
		instrument = f.gets.chomp
		triggerid = f.gets.chomp
		seqnum = f.gets.chomp
		queue = f.gets.chomp
		alert_name = f.gets.chomp
		run_name = f.gets.chomp
		basedir = $root + "/" + instrument + "/" + alert_name + "/"
		mle = run_name
		runname = run_name
		load_build_command  = ENV["MODULELOAD"]

	elsif analysistype=="CALPIPE"

		puts "calpipe"

		path_3908 = f.gets.chomp
		path_3916 = f.gets.chomp
		orbit_number = f.gets.chomp
		mle = orbit_number
		puts path_3908+" "+path_3916+" "+orbit_number
		runname = orbit_number
		basedir = "/ANALYSIS3/AGILE-MCAL/"+orbit_number
		load_build_command  = "agile-mcal"

	elsif analysistype=="GRID-PIPE-SHORT"

		load_build_command = "agile-preB25_iteration2-r5"

		puts "grid pipe short"

		basedir = f.gets.chomp
		puts basedir

		runname = f.gets.chomp
		puts runname

		mle = runname

	else
		index = 0;

		conffile = Conf.new

		conffile.process(filenameconf, nil, nil)

		basedir = $root + "/" + conffile.dir_run_output + "/" + conffile.run_name

		if basedir == $root + "//"
			puts "error in config file name"
			exit
		end

		mleindex = 0;
		if File.exists?(basedir)
			ml = Dir[basedir+"/MLE????"].sort
			puts basedir + " with index " + ml.size().to_s
			if ml.size() > 0
				mleindex = ml[ml.size()-1].split("MLE")[1].to_i;
				mleindex = mleindex.to_i + 1
			else
				cmd = "rm -rf " + basedir
				puts cmd
				system(cmd)
			end
		end

		mle = "MLE" + format("%04d", mleindex)
		puts mle

		queue = conffile.queue
		runname = conffile.run_name
		load_build_command = conffile.load_build_command
	end

	#after analysis type detection

	cmd = "mkdir -p " + basedir;
	puts cmd
	system(cmd);

	# create a .con file from the contour available in the .conf (if not exists)
	contourfile = basedir + "/" + triggerid.to_s + "_" + seqnum.to_s + ".con"
	cmd = "if [ ! -f " + contourfile + " ] ; then tail -n +$((1 + $(grep \"\\-\\-\\-\\-\" -n " + filenameconf + " | tail -n1 | cut -d':' -f 1) )) " + filenameconf + " > " + contourfile + "; fi"
	puts cmd
	system(cmd)

	#copy the .conf
	cmd = "cp " + filenameconf + " " + basedir + "/" + mle + ".conf"
	puts cmd
	puts "-"+mle+"-"
	system(cmd)

	newcmd = basedir + "/" + mle + ".ll"

	cmd = "cp " + agilepipe_path + "/template.ll " + newcmd
	puts cmd
	system(cmd)

	f = File.open(newcmd, "a")
	if queue != nil
		f.write("\#\@ class    = " + queue + "\n")
		#for slurm
		f.write("\#SBATCH --partition=" + queue + "\n")
	else
		f.write("\#\@ class    = large\n")
		#for slurm
		f.write("\#SBATCH --partition=large\n")
	end
	f.write("\#\@ job_name = sor4_" + runname + "\n")

	#enable/disable send mail when the task is finished
	#f.write("\#\@ notify_user = " + conffile.mail + "\n")

	f.write("\#\@ queue\n")



	#for slurm
	if analysistype.include? "CALPIPE"
		f.write("\#SBATCH --job-name=CALPIPE_"+runname+"\n")
	end
	if analysistype.include? "GRID-PIPE-SHORT"
		f.write("\#SBATCH --job-name=GRID-SHORT-"+runname+"\n")
	end

	f.write("\#SBATCH --share\n")


	f.write("date\n")
	f.write(". /usr/share/Modules/init/bash\n")
	f.write("module load " + load_build_command + "\n")

	output = %x[xvfb-run -h]
	lines = output.split("\n")

	xvfb_option = "-a"
	lines.each do |line|
		#puts "line-"+line
		if(line.start_with?('-d'))
			#puts "preso -d"
			xvfb_option = "-d"
		end
	end
	puts "xvfb-option"+xvfb_option

	if analysistype == "GRB" or analysistype == "GW"
		f.write("xvfb-run "+xvfb_option+" -s \"-screen 0 2000x2000x24\" "+agilepipe_path+"/grb_gw/grb_gw.py " + mle + ".conf\n")

	elsif analysistype == "visCheckAux"
		f.write("module load idl8.2sp2\n")
		f.write(agilepipe_path + "/visCheck/burst_position.sh " + mle + ".conf\n")

	elsif analysistype.include? "MCAL"
		f.write(agilepipe_path + "/mcal_pipe.sh " + mle + ".conf\n")

	elsif analysistype.include? "Ratemeters"
		f.write("module load heasoft-6.17\n")
		f.write("heainit\n")
		f.write(agilepipe_path + "/ratemeters/ratemeters_pipe.sh " + mle + ".conf\n")

	elsif analysistype.include? "CALPIPE"
		puts "include"
		f.write("xvfb-run "+xvfb_option+" -s \"-screen 0 2000x2000x24\" python /opt/prod/AGILE-MCAL-PIPE/pipe/MCAL-PIPE.py "+orbit_number+" "+path_3908+" "+path_3916+" "+basedir+" 7.0 4.0\n")
	elsif analysistype.include? "GRID-PIPE-SHORT"
		puts "include"
		f.write("ruby /opt/prod/AGILE-GRID-PIPE-SHORT/pipe/grid_short_analysis.rb "+basedir+" "+mle+".conf")

	else
		f.write(agilepipe_path + "/spot6/analysis.rb " + mle + ".conf" + "\n")
	end

	f.close()

	puts basedir

	Dir.chdir(basedir)
	cmd = "cd " + basedir + "; " + EXECCOM + " " + newcmd;
	puts cmd
	system(cmd)

end

while 1 do
	a = Dir[PATH_RES + "/commands/*.conf"]
	if a.size() > 0
		a.sort.each do | line |
			line = line.gsub(/["\\]/,'') # remove double slashes
			puts line
			cmd = "mv " + line + " /tmp/";
			puts cmd
			system(cmd);
			file = line.split("/")[line.split("/").size-1]
			clusteranalysis("/tmp/" + file)
			job_list =  %x[squeue --format=%j]
			while(job_list.split("\n").length >= 200)
				sleep(1)
				job_list = %x[squeue --format=%j]
			end
			ll = line.split("/")
			system("rm /tmp/" + ll[ll.size-1]);
		end
	end
	sleep(1)
end
