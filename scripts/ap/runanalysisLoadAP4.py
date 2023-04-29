import glob
import os

def main():
	regex="RES1_*.ap4"
	log_name_list =  sorted(glob.glob(regex), key=os.path.getmtime)
	for filename in log_name_list:	
		print(filename)
		ranal = filename.split("_")[2].split("R")[1]
		os.system("python analyzeLoadAP4.py " + filename)

if __name__ == '__main__':
	main()