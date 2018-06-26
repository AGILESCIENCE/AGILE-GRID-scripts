#!/bin/bash

# requires heasoft ftools !

trap exit ERR

# Define default directories
caldir="$AGILE/model/scientific_analysis/data"

# Read the arguments
E_BADARGS=65
if [ $# -lt 2 ]
then
    echo "Usage: `basename $0` inputfile expmaplist"
    exit $E_BADARGS
fi

# Define output filenames
infile=$1
base=${infile%.*}
expmaplist=$2

# Define file names of convolved maps
convemin=( ${convemin:-00100 00400 01000 03000 10000} )
convemax=( ${convemax:-00400 01000 03000 10000 50000} )
echo "convemin  = ${convemin[*]}"
echo "convemax  = ${convemax[*]}"
echo ""
echo "adding.."
for i in ${!convemin[*]} ; do
    convbase[$i]=E${convemin[$i]}_${convemax[$i]}_${base}
    convfile[$i]=${convbase[$i]}.conv.sky
done

# Read exposure map file
expfiles=()
indices=()
while read -a line ; do
    expfiles=("${expfiles[@]}" "${line[0]}")
    indices=("${indices[@]}" "${line[1]}")
done < "${expmaplist}"

# For each exposure map, add convolved maps and generate template file
for i in ${!expfiles[*]} ; do
    echo ${expfiles[i]}
    ftkeypar ${expfiles[$i]} MINENG
    addemin=`pget ftkeypar value`
    addemin=`echo ${addemin} | sed 's/[.].*//'`
    ftkeypar ${expfiles[$i]} MAXENG
    addemax=`pget ftkeypar value`
    addemax=`echo ${addemax} | sed 's/[.].*//'`
    echo $addemin $addemax
    templatebase=$(basename ${expfiles[i]} .gz)
    templatebase=$(basename ${expfiles[i]} .exp)
    convfilelist=${templatebase}.convlist.in
    echo "${indices[$i]} ${#convemin[*]}" > ${convfilelist}
    for j in ${!convemin[*]} ; do
        weight=`$AGILE/scripts/extendedsources/specwt.py ${convemin[$j]} ${convemax[$j]} ${addemin} ${addemax} ${indices[$i]}`
        fextract ${convfile[$j]}+0 \!${convfile[$j]}.specwt
        fcarith ${convfile[$j]}+1 ${weight} \!${convfile[$j]}.specwttemp MUL
        fappend ${convfile[$j]}.specwttemp+0 ${convfile[$j]}.specwt
        rm ${convfile[$j]}.specwttemp
        echo "${convfile[$j]}.specwt" >> ${convfilelist}
    done
    dispfile=${templatebase}.disp.conv.sky
    AG_add_diff5 diffusefilelist=${convfilelist} edpfile=${caldir}/AG_GRID_G0017_SFMG_I0025.edp.gz sarfile=${caldir}/AG_GRID_G0017_SFMG_I0025.sar.gz outfile=\!${dispfile} emin=${addemin} emax=${addemax} > ${dispfile}.out 2>&1
    ftstat ${dispfile} chatter=0
    dispsum=`pget ftstat sum`
    fextract ${dispfile}+0 \!norm_${dispfile}
    fcarith ${dispfile}+1 ${dispsum} \!temp_${dispfile} DIV
    fappend temp_${dispfile}+0 norm_${dispfile}
    rm temp_${dispfile}
    # was outfile=\!
    AG_gasmapgen5 expfile=${expfiles[i]} outfile=${templatebase}.template.gz diffusefile=$AGILE/scripts/extendedsources/diffuse_null.fits hiresdiffusefile=norm_${dispfile} > ${templatebase}.template.out 2>&1
done
