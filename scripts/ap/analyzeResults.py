import os
import sys
if __name__ == '__main__':
	filename = sys.argv[1]
	period = 8.4474
	fractionofday = 2/24.0
	with open(filename, "r") as ins:
		for line in ins:
			if(line != ""):
				val = line.split(" ")
				an1 = str(val[0])
				an2 = str(val[1])
				sig = float(val[2])
				time = float(val[5])

				if time >= period - fractionofday and time <= period + fractionofday:
					print(filename + " - " + line)
					fff = filename.split(".ap")[0]+".ap"
					vmt = an1.split("mu")[0]
					if vmt != "LS":
						os.system("python drawVM.py "+fff+"."+vmt+".resgf " + an2)
					else:
						os.system("python drawLS.py "+fff+".ap2 " + an2)

			#	if time >= period*2 - fractionofday and time <= period*2 + fractionofday:
                         #               print(filename + " - " + line)

			#	if time >= period/2 - fractionofday and time <= period/2 + fractionofday:
                         #               print(filename + " - " + line)
