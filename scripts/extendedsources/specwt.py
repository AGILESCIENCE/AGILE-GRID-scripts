#!/usr/bin/env python

def specwt(ei, ei1, emin, emax, index):
	if index < 0 :
		index = (-index)
	weight = (ei ** (1.0-index) - ei1 ** (1.0-index)) / (emin ** (1.0-index) - emax ** (1.0-index))
	return weight
if __name__ == "__main__":
	import sys
	print specwt(float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]))
