import glob
import os

if __name__ == '__main__':
	regex="RES1_*.ap3"
	log_name_list =  sorted(glob.glob(regex), key=os.path.getmtime)
	for filename in log_name_list:	
		print(filename)
		ranal = filename.split("_")[2].split("R")[1]
		os.system("python analyzeLoadAP3.py " + filename + " " + ranal)
