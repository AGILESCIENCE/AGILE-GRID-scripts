#! /usr/bin/ruby
#See map.rb parameters +
#00) timelist: a file with a list of tstart stop in TT
#01) filter DIR (es: FM3.119_2_I0023, FM3.119_ASDCe_I0023, FM3.119_ASDCSTDf_I0023, FM3.119_ASDCSTDk_I0023)
#02) output file name prefix
#03) l (map center)
#04) b (map center)
#05) emin for merge
#06) emax for merge
#07) a string with the remaining parameters of the map.rb

require 'fileutils'

if ARGV[0].to_s == "help" || ARGV[0] == nil || ARGV[0] == "h"
	system("head -11 " + $0 );
	exit;
end

timelist = ARGV[0]
filter = ARGV[1];
name = ARGV[2];
l = ARGV[3];
b = ARGV[4];
emin = ARGV[5];
emax = ARGV[6];
commands = ARGV[7]

tmpdir = timelist + rand(100000)
Dir.mkdir(tmpdir)
FileUtils.cp(timelist, tmpdir)
Dir.chdir(tmpdir)


i=0
File.open(timelist).each do | line |
	ll = line.split(" ")
	t1 = ll[0]
	t2 = ll[1]
	cmd = "map.rb "+filter+" "+name+format("%05d",i)+" "+ t1.to_s + " " + t2.to_s + " "+ l.to_s + " " + b.to_s + " " + commands
	puts cmd
	system cmd
	i = i + 1
end

cmd = "mapmerge.rb " + name + " " + name + " " + filter + " " + emin + " " + emax + " -1 -1 4"
puts cmd
system cmd

