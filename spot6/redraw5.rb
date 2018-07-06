#! /usr/bin/ruby
#0) config file name

load ENV["AGILE"] + "/scripts/conf.rb"
load ENV["AGILE"] + "/AGILEPIPE/env.rb"
load ENV["AGILE"] + "/AGILEPIPE/spot6/Conf.rb"

datautils = DataUtils.new
alikeutils = AlikeUtils.new
parameters = Parameters.new

filenameconf = ARGV[0];

mle = filenameconf.split(".conf")[0]

filenameconfext = filenameconf

#estrazione lista sorgenti
fndisplayreg = mle + "display"

fnhyp0 = mle+"hypothesis0.multi"
fnhyp = mle+"hypothesis.multi"

conf = Conf.new

conf.process(filenameconf, nil, nil);

conf.detsmooth()
	
conf.plotjpgcts1(mle, conf.smooth)

conf.plotjpgint(mle, conf.smooth)

conf.plotjpgexp(mle)

conf.plotjpgcts2(mle, conf.smooth)

conf.plotjpgcts2(mle + ".step0", conf.smooth)

conf.plotjpgcts2(mle + ".step1", conf.smooth)

conf.copyresults(mle)



