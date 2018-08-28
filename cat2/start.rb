

a = Dir.pwd.split("/")[6].split(".")
system("cp ../S3F1/OUTS3F1.res.multi list000.multi")
system("sh start.sh")
system("sh extractcat.sh " + a[0].to_f.to_s + " " + a[1].to_f.to_s + " >> list000.multi")
system("checkmultiinput.rb list000.multi list001.multi 1")
system("multi5.rb FM3.119_ASDCe_I0023_B05.maplist4 list001.multi RES001 galcoeff=0.7")
system("convertMultiResToInput.rb RES001 list002.multi 1 2.9 0.0 90 1")
system("makecatalog_phase3.rb list002.multi FM3.119_ASDCe_I0023_B05.maplist4 FM3.119_ASDCe_I0023 DIR003_00 3 0 0 \"flag=SCAN3_00\" 1")
system("cd DIR003_00; extractres2.rb DIR003_00.resfinalfull DIR003_00.res2 " + a[0] + "." + a[1].to_s + "; less -S DIR003_00.resfinalfull ")


