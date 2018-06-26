#! /usr/bin/ruby
#0) l1
#1) b1
#2) l2
#3) b2

load ENV["AGILE"] + "/scripts/conf.rb"

if ARGV[0].to_s == "help" || ARGV[0].to_s == "h" || ARGV[0] == nil
	system("head -5 " + $0 );
	exit;
end

l1 = ARGV[0]
b1 = ARGV[1]
l2 = ARGV[2]
b2 = ARGV[3]

datautils = DataUtils.new

d = datautils.distance(l1, b1, l2, b2)

puts d
