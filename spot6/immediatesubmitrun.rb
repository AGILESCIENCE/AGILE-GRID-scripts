#! /usr/bin/ruby

# Copyright (c) 2016, AGILE team
# Authors: Andrea Bulgarelli <bulgarelli@iasfbo.inaf.it>,
#          Andrea Zoli <zoli@iasfbo.inaf.it>
#
# Any information contained in this software is property of the AGILE TEAM
# and is strictly private and confidential. All rights reserved.

load ENV["AGILEPIPE"] + "/env.rb"
load ENV["AGILEPIPE"] + "/spot6/Conf.rb"

$root = ENV['PATH_RES']

def clusteranalysis(filenameconf)
	agilepipe_path = ENV["AGILEPIPE"] + "/"

	analysistype = nil
	instrument = nil
	triggerid = nil
	queue = nil
	alert_name = nil
	job_name = nil

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
	else
		index = 0;

		conffile = Conf.new

		conffile.process(filenameconf, nil, nil)

		job_name = conffile.dir_run_output

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
	puts "job_name"+job_name
	f.write("\#SBATCH --job-name="+ job_name + "\n")
	#enable/disable send mail when the task is finished
	#f.write("\#\@ notify_user = " + conffile.mail + "\n")

	f.write("\#\@ queue\n")



	#for slurm
	f.write("\#SBATCH --share\n")


	f.write("date\n")
	f.write(". /usr/share/modules/init/bash\n")
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
		f.write(agilepipe_path + "/visCheck/burst_position.sh " + mle + ".conf\n")
	elsif analysistype.include? "MCAL"
		f.write(agilepipe_path + "/mcal_pipe.sh " + mle + ".conf\n")
	elsif analysistype.include? "Ratemeters"
		f.write("module load heasoft-6.17\n")
		f.write("heainit\n")
		f.write(agilepipe_path + "/ratemeters/ratemeters_pipe.sh " + mle + ".conf\n")
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

	puts "file "+ARGV[0]
	clusteranalysis(ARGV[0])
	puts "delete "+ARGV[0]
	system("rm " + ARGV[0]);
