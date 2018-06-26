#!/usr/bin/env python
import sys, argparse, struct
import numpy
from astropy import units as u
from astropy.coordinates import SkyCoord

# this is an hammer function normalized
def aitoff(l, b):
  """
    Parameters:
     - `l`, `b` - float, array
                   The longitude and latitude [deg]
    Returns:
      Aitoff-projected coordinates (`x`, `y`)

    .. note:: This function was ported from the IDL Astronomy User's Library.
    
    :IDL - Documentation:
  
    pro aitoff,l,b,x,y
    +
     NAME:
           AITOFF
     PURPOSE:
           Convert longitude, latitude to X,Y using an AITOFF projection.
     EXPLANATION:
           This procedure can be used to create an all-sky map in Galactic 
           coordinates with an equal-area Aitoff projection.  Output map 
           coordinates are zero longitude centered.
    
     CALLING SEQUENCE:
           AITOFF, L, B, X, Y 
    
     INPUTS:
           L - longitude - scalar or vector, in degrees
           B - latitude - same number of elements as L, in degrees
    
     OUTPUTS:
           X - X coordinate, same number of elements as L.   X is normalized to
                   be between -180 and 180
           Y - Y coordinate, same number of elements as L.  Y is normalized to
                   be between -90 and 90.
    
     NOTES:
           See AIPS memo No. 46, page 4, for details of the algorithm.  This
           version of AITOFF assumes the projection is centered at b=0 degrees.
    
     REVISION HISTORY:
           Written  W.B. Landsman  STX          December 1989
           Modified for Unix:
                   J. Bloch        LANL SST-9      5/16/91 1.1
           Converted to IDL V5.0   W. Landsman   September 1997
  """
  wasFloat = False
  if isinstance(l, float):
    l = numpy.array([l]); b = numpy.array([b])
    wasFloat = True
  
  sa = l.copy()
  x180 = numpy.where(sa > 180.0)[0]
  if len(x180) > 0: sa[x180] -= 360.
  alpha2 = sa/(2*(180.0/numpy.pi))
  delta = b/(180.0/numpy.pi)   
  r2 = numpy.sqrt(2.)    
  f = 2*r2/numpy.pi   
  cdec = numpy.cos(delta)    
  denom = numpy.sqrt(1. + cdec*numpy.cos(alpha2))
  x = cdec*numpy.sin(alpha2)*2.*r2/denom
  y = numpy.sin(delta)*r2/denom
  x = x*(180.0/numpy.pi)
  y = y*(180.0/numpy.pi)
  
  if wasFloat:
    return float(x), float(y)
  else:
    return x, y

def inverseAitoff(x, y):
  """
    Carry out an inverse Aitoff projection.
    
    Parameters:
     - `x` - float or array, A value between -180. and +180. (see convention in *aitoff* function).
     - `y` - float or array, A value between -90. and +90. (see convention in *aitoff* function).
     
    This function reverts to aitoff projection made by the function *aitoff*. The result is \
    either two floats or arrays (depending on whether float or array was used as input) representing \
    longitude and latitude. Both are given in degrees with -180 < longitude < +180 and \
    -90 < latitude < 90.
    
    If arrays are used for input, the function returns an array for the longitude, for the latitude, and \
    an index array containing those array indices for which the reprojection could be carried out.
  """
  wasFloat = False
  if isinstance(x, float):
    x = numpy.array([x]); y = numpy.array([y])
    wasFloat = True
  
  # First, rescale x and y
  x = x/180.0 * 2.0 * numpy.sqrt(2.0)
  y = y/90.0 * numpy.sqrt(2.0)
  
  zsqr = 1.0 - (x/4.0)**2 - (y/2.0)**2
  # Check whether x,y coordinates are within the ellipse of invertible values.
  indi = numpy.where((x**2/8. + y**2/2. - 1.0) <= 0.0)[0]
  if len(indi) == 0:
    return 0, 0, 0
  
  z = numpy.sqrt(zsqr[indi])
  l = 2.0 * numpy.arctan( z*x[indi]/(2.0 * (2.0*z**2 - 1.0)) )
  b = numpy.arcsin(z * y[indi])
  
  l = l*180.0 / numpy.pi
  b = b*180.0 / numpy.pi
  
  if wasFloat:
    return float(l), float(b)
  else:
    return l, b, indi

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("polygon", help="An input file containing the polygon as a list of vertices (default radec, see -g)")
    parser.add_argument("points", help="An input file containing a list of points to test (galactic coords only)")
    parser.add_argument("output", help="An output file containing only the points inside the polygon")
    parser.add_argument("-g", "-galactic", help="The input contour unit is galactic coordinates", action='store_true')
    parser.add_argument("-l", "-lcolumn", help="l column number inside the list of points file (default 0)", type=int, default=0)
    parser.add_argument("-b", "-bcolumn", help="b column number inside the list of points file (default 1)", type=int, default=1)
    args = parser.parse_args()

    print "Parsing contour file.."
    fc = numpy.loadtxt(args.polygon)
    c_gal_l = fc[:,0].astype(numpy.float)
    c_gal_b = fc[:,1].astype(numpy.float)
    if not args.g:
        print "Converting contour from radec to galactic.."
        for i in xrange(c_gal_l.size):
            tmp_radec = SkyCoord(ra=c_gal_l[i]*u.degree, dec=c_gal_b[i]*u.degree, frame='fk5')
            tmp_gal = tmp_radec.transform_to('galactic')
            c_gal_l[i] = tmp_gal.l.value
            c_gal_b[i] = tmp_gal.b.value
        numpy.savetxt(args.polygon+'.galactic', numpy.c_[c_gal_l, c_gal_b], delimiter=' ', fmt="%f")
    c_x, c_y = aitoff(c_gal_l, c_gal_b)
    nv = c_gal_l.shape[0]
#    numpy.savetxt(args.polygon+'.ait', numpy.c_[c_x, c_y], delimiter=' ', fmt="%f")

    print "Parsing points file.."
    lines=[]
    with open(args.points, 'r') as fl:
        for line in fl.readlines():
            lines.append(line.rstrip())

    fp = numpy.genfromtxt(args.points, delimiter=' ', dtype=None)
    p_gal_l=numpy.array([])
    p_gal_b=numpy.array([])
    for i in xrange(len(fp)):
        p_gal_l = numpy.append(p_gal_l, fp[i][args.l].astype(numpy.float))
        p_gal_b = numpy.append(p_gal_b, fp[i][args.b].astype(numpy.float))
    p_x, p_y = aitoff(p_gal_l, p_gal_b)

#    numpy.savetxt(args.points+'.ait', numpy.c_[p_x, p_y], delimiter=' ', fmt="%f")

    fo = open(args.output, 'w')

    print "Filtering.."
    for p in xrange(p_x.shape[0]):
        testx = p_x[p];
        testy = p_y[p];
        c = False
        j = nv-1
        for i in xrange(nv):
            if ( ((c_y[i]>testy) != (c_y[j]>testy)) and \
                 (testx < (c_x[j]-c_x[i]) * (testy-c_y[i]) / (c_y[j]-c_y[i]) + c_x[i]) ):
                c = not c
            j = i
#        print (p_x[p], p_y[p]),
        if c:
            fo.write(lines[p])
            fo.write('\n')
#            print "inside!"
#        else:
#            print "outside!"

    fo.close()
    print "Done."

if __name__ == "__main__":
   main(sys.argv[1:])
