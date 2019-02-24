import glob
import os

if __name__ == '__main__':
	regex="RES_*.ap"
	log_name_list =  sorted(glob.glob(regex), key=os.path.getmtime)
	for filename in log_name_list:	
		print(filename)
		os.system("python analyzeAP.py " + filename)
