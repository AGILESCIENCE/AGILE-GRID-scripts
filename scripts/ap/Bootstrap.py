from random import shuffle
from time import time
import string, os, sys
import numpy as np
import sys
import math

class Bootstrap:

	def __init__(self):
		return

	def bootstrap(self, apfilename, apfilename2):
		rows_1_2_3 = []
		rows_4 = []

		with open(apfilename) as f:
				content = f.readlines()
		
		# Removing whitespace characters like `\n` at the end of each line
		lines = [x.strip() for x in content]
			
		for line in lines:
				columns = line.split(" ")
				rows_1_2_3.append(columns[0]+" "+columns[1]+" "+columns[2])
				rows_4.append(columns[3])
			
		shuffle(rows_4)
			
		with open(apfilename2, 'a') as f:
				for i in range(len(rows_4)):
				f.write(rows_1_2_3[i]+" "+rows_4[i]+"\n")

	def generateBootstrap(self, apfilename, nrun):
		for x in range(int(nrun)):
			apfilename2 = apfilename + "_"  + str(time())+".ap"
			self.bootstrap(apfilename, apfilename2)	
