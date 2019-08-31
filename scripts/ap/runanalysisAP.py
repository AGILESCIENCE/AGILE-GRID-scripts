import glob
import os

if __name__ == '__main__':
	regex="RES_*.ap"
	log_name_list =  sorted(glob.glob(regex), key=os.path.getmtime)
	for filename in log_name_list:	
		print(filename)
		ranal = filename.split("_")[2].split("R")[1]
		os.system("python analyzeAP.py " + filename + " " + ranal)
