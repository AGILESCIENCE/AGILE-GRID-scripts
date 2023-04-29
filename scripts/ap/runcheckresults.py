import glob
import os

def main():
	regex="*.apres"
	log_name_list =  sorted(glob.glob(regex), key=os.path.getmtime)
	for filename in log_name_list:	
		print(filename)
		os.system("python analyzeResults.py " + filename)

if __name__ == '__main__':
	main()
