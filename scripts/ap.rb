#!/usr/bin/env ruby

# NOTE: passing ARGV + execap + timeslot to the map.rb
ARGV << 'execap=1 timeslot=100'
load ENV["AGILE"] + "/scripts/map.rb"

apfile = ARGV[1] + ".ap"
ARGV.clear
ARGV[0] = apfile
load ENV['AGILE'] + '/scripts/weight.rb'
