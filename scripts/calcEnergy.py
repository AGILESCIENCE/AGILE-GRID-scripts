import math
import sys

def ergLog(alpha, flux, flux_err, e_min, e_max):
	alpha = -alpha
	e_log =  (math.log10( e_max)-math.log10( e_min))/2 + math.log10( e_min)
	e_log = math.pow(10, e_log)
	erg  = (  flux / ( e_max -  e_min) ) *  e_log *  e_log * 1.6 / 1000000.
	erg_err = (  flux_err / ( e_max -  e_min) ) *  e_log *  e_log * 1.6 / 1000000.
	
	print('Log:')
	print(erg)
	print(erg_err)

def ergGeom(alpha, flux, flux_err, e_min, e_max):
	alpha = -alpha
	e_geom = (( e_max- e_min)/2 +  e_min)
	erg_two  = (  flux / ( e_max -  e_min) ) *  e_geom *  e_geom * 1.6 / 1000000.
	erg_two_err = (  flux_err / ( e_max -  e_min) ) *  e_geom *  e_geom * 1.6 / 1000000.
	
	print('Geom:')
	print(erg_two)
	print(erg_two_err)
	
def ergGiuliani(alpha, flux, flux_err, e_min, e_max):
	Ec = 0.0;
	alpha = -alpha

	if( alpha != -2):
		Ec = ( 1 / ( alpha+2) ) *  ( math.pow( e_max,( alpha+2))  -  math.pow( e_min,( alpha+2)) ) / ( ( 1 / ( alpha+1) ) * ( math.pow( e_max,( alpha+1))  -  math.pow( e_min,( alpha+1)) ) )
	else:
		Ec = log( e_max/ e_min) / ( ( math.pow( e_max,( alpha+1)) - math.pow( e_min,( alpha+1))) / ( alpha + 1 ) )

	erg_three  = ( flux / ( e_max -  e_min) ) *  Ec *  Ec * 1.6 / 1000000. ;
	erg_three_err = ( flux_err / ( e_max - e_min) ) *  Ec * Ec * 1.6 / 1000000. ;

	print('Giuliani:')
	print(erg_three)
	print(erg_three_err)

def ergChen(alpha, flux, flux_err, e_min, e_max):
	
	erg_four = 0
	erg_four_err = 0
	if(alpha == 2):
		erg_four =  flux *  e_min * math.log( e_max / e_min) / (1 - ( e_min /  e_max)) * 1.6 / 1000000.	
		erg_four_err =  flux_err *  e_min * math.log( e_max/ e_min) / (1 - ( e_min /  e_max)) * 1.6 / 1000000.
	else:
		erg_four =  flux *  e_min * ((1- alpha)/(2- alpha)) * ((math.pow(( e_max / e_min),(2- alpha))-1) / (pow(( e_max/ e_min), (1- alpha)) - 1)) * 1.6 / 1000000.
		erg_four_err =  flux_err *  e_min * ((1- alpha)/(2- alpha)) * ((math.pow(( e_max/ e_min),(2- alpha))-1) / (math.pow(( e_max/ e_min),(1- alpha))-1)) * 1.6 / 1000000.

	print('Chen:')
	print(erg_four)
	print(erg_four_err)

if __name__ == '__main__':

    # read the XML for the specific observation
    alpha = float(sys.argv[1])
    flux = float(sys.argv[2])
    flux_err = float(sys.argv[3])
    e_min = float(sys.argv[4])
    e_max = float(sys.argv[5])
    
    ergLog(alpha, flux, flux_err, e_min, e_max)
    ergGeom(alpha, flux, flux_err, e_min, e_max)
    ergGiuliani(alpha, flux, flux_err, e_min, e_max)
    ergChen(alpha, flux, flux_err, e_min, e_max)