#!/bin/bash

# requires heasoft ftools !

trap exit ERR

caldir="$AGILE/model/scientific_analysis/data"
newxbin=0.1

# Read the arguments
E_BADARGS=65
if [ $# -lt 2 ]
then
    echo "Usage: `basename $0` inputfile spectral_index"
    exit $E_BADARGS
fi

# Define output filenames
infile=$1
base=${infile%.*}
spectral_index=$2
convertedfile=${base}.converted.sky
rebinfile=${base}.rebinned.sky
padfile=${base}.pad.sky

convemin=( ${convemin:-00100 00400 01000 03000 10000} )
convemax=( ${convemax:-00400 01000 03000 10000 50000} )
convindex=( ${convindex:-${spectral_index} ${spectral_index} ${spectral_index} ${spectral_index} ${spectral_index}} )
echo "convemin  = ${convemin[*]}"
echo "convemax  = ${convemax[*]}"
echo "convindex = ${convindex[*]}"
echo ""
echo "convolving.."

# Rebin input file
AG_converttoSkyMap5 $infile \!$convertedfile
ftkeypar ${convertedfile} CDELT1
xbin=`pget ftkeypar rvalue`
multiple=`echo "(${xbin}>=0)*(${newxbin}/${xbin})+(${xbin}<0)*(${newxbin}/(-(${xbin})))" | bc`
#echo "$newxbin / $xbin = $multiple"
fimgbin $convertedfile+1 $rebinfile $multiple copyall=yes clobber=yes

# Pad file with blank pixels
AG_padMap5 $rebinfile \!$padfile
fparkey GLON $padfile CTYPE1
fparkey GLAT $padfile CTYPE2

# Save PFILES and get parameter file directory
PFILES_OLD=${PFILES}
PFILES_DIR=`echo ${PFILES} | awk -F \; '{print $1}'`

# Convolve template file for different energy bins
for i in ${!convemin[*]} ; do
    convbase[$i]=E${convemin[i]}_${convemax[i]}_${base}
    convfile[$i]=${convbase[$i]}.conv.sky
    PFILES=${convbase[$i]}_pfiles
    mkdir ${PFILES}
    cp ${AGILE}/share/AG_diff_conv5.par ${PFILES}
    cp ${HEADAS}/syspfiles/fparkey.par ${PFILES}
    outname=${convbase[$i]}.conv.out
    fparkey ${convemin[i]} $padfile E_MIN add=yes insert=CRVAL2
    fparkey ${convemax[i]} $padfile E_MAX add=yes insert=E_MIN
    fparkey ${convindex[i]} $padfile INDEX add=yes insert=E_MAX
    ( AG_diff_conv5 diffusefile=$padfile sarfile=${caldir}/AG_GRID_G0017_SFMG_I0025.sar.gz psdfile=${caldir}/AG_GRID_G0017_SFMG_I0025.psd.gz edpfile=${caldir}/AG_GRID_G0017_SFMG_I0025.edp.gz outfile=\!${convfile[$i]} > $outname 2>&1 ; fparkey 0 $padfile -E_MIN ; fparkey 0 $padfile -E_MAX ; fparkey 0 $padfile -INDEX ; rm ${PFILES}/*.par ; rmdir ${PFILES} )
done
PFILES=${PFILES_OLD}

