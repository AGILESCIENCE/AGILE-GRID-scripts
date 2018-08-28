a = Dir.pwd.split("/")[6].split(".")
system("rm -rf DIR003_00")
system("makecatalog_phase3.rb list002.multi FM3.119_ASDCe_I0023_B05.maplist4 FM3.119_ASDCe_I0023 DIR003_00 3 0 0 \"flag=SCAN3_00\" 1")
system("cd DIR003_00; extractres2.rb DIR003_00.resfinalfull DIR003_00.res2 " + a[0] + "." + a[1].to_s + ";  less -S DIR003_00.resfinalfull ")

