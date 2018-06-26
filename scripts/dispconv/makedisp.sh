#!/bin/bash

trap exit ERR

usage() {
    cat << EOF
Usage: $(basename $0) [ options ... ]
This script search for a file .in in the current directory and generate
a .disp.conv.sky for each .conv.sky defined in the .in. To use different
energy ranges see the -r option.

  -h                print help
  -d PATH           model data path
  -n EMIN_STRING    additional user defined emin energies
  -x EMAX_STRING    additional user defined emax energies

  Default options: -d \$AGILE/model/scientific_analysis/data
EOF
    exit $1
}

# Check AGILE environment
if [ -z "$AGILE" ] ; then
    echo "\$AGILE environment variable is not set, quitting.."
    exit 1
fi

# Parse options
modelpath="$AGILE/model/scientific_analysis/data"
while getopts "hm:f:d:s:n:x:" opt; do
    case "${opt}" in
        h) usage 0;;
        d) modelpath="$OPTARG";;
        n) user_eminv=($OPTARG);;
        x) user_emaxv=($OPTARG);;
        *) usage 1;;
    esac
done
shift $(($OPTIND-1))

# Check user defined ranges
if [ ${#user_eminv[*]} != ${#user_emaxv[*]} ] ; then
    echo "Error on range definition, the two array emin and emax should have the same size."
    echo "user emin = \"${user_eminv[*]}\""
    echo "user emax = \"${user_emaxv[*]}\""
    exit
fi

# Check user defined ranges min < max
if [ ${#user_eminv[*]} != ${#user_emaxv[*]} ] ; then
    echo "Error on range definition, the two array emin and emax should have the same size."
    echo "user emin = \"${user_eminv[*]}\""
    echo "user emax = \"${user_emaxv[*]}\""
    exit
fi

# Check .in naming convention
in=$(find . -maxdepth 1 -name "AG_add_diff\.SKY*\.S*_*\.in" | head -n1)
if [ -z $in ] ; then
    in=$(find . -maxdepth 1 -name "AG_add_diff\.SKY*\.S*_*\.in" | head -n1)
fi
if [ -z $in ] ; then
    echo "Cannot find a valid .in file in the current directory."
    echo "Expected an input file 'AG_add_diff.<sky>.S<filter>_<matrix>.in', quitting.."
    exit
fi
in=$(basename $in)

# Get informations from .in file naming
tmp=${in#*.}
sky=${tmp%.*_*.*}
tmp=${in#*.*.S}
filter=${tmp%_*.*}
tmp=${in#*.*.S*_}
matrix=${tmp%.in}

# Check if model files exists
echo "Checking model files.."
modelnaming="AG_GRID_G0017_S${filter}_${matrix}"
sarfile="$modelpath/$modelnaming.sar.gz"
psdfile="$modelpath/$modelnaming.psd.gz"
edpfile="$modelpath/$modelnaming.edp.gz"
if [ ! -f $sarfile ] ; then
    echo "Cannot find the model file $sarfile}.. wrong model path?"
    exit 1
fi
if [ ! -f $psdfile ] ; then
    echo "Cannot find the model file ${psdfile}.. wrong model path?"
    exit 1
fi
if [ ! -f $edpfile ] ; then
    echo "Cannot find the model file ${edpfile}.. wrong model path?"
    exit 1
fi

# Check if sky maps exists
echo "Checking $in file.."
eminv=()
emaxv=()
for skymap in $(tail -n +2 $in); do
    mapname=$(basename $skymap)
    if [[ ! $mapname =~ [0-9]*_[0-9]*\.${sky}\.S${filter}_${matrix}\.conv\.sky\.gz ]] ; then
        echo "Error with sky map $skymap, bad naming."
        echo "Expected a conv map like 'EMIN_EMAX.${sky}.S${filter}_${matrix}.conv.sky.gz'."
        exit
    fi
    emin=${mapname%_*_*}
    tmp=${mapname#*_}
    emax=${tmp%%.*}
    eminv=("${eminv[@]}" "${emin}")
    emaxv=("${emaxv[@]}" "${emax}")
done

# Check matching between user defined ranges and map ranges
for i in ${!user_eminv[*]} ; do
    if [[ ${user_eminv[i]} -ge ${user_emaxv[i]} ]] ; then
        echo "Error with user defined energy ranges. The value of emin ${user_eminv[i]} is greater or equal to emax ${user_emaxv[i]}."
        exit
    fi
    found=false
    for j in ${!eminv[*]} ; do
        [ ${user_eminv[i]} == ${eminv[j]} ] && found=true && break
    done
    if [ ! "$found" = true ] ; then
        echo "Error with user defined energy ranges. The emin value \"${user_eminv[i]}\" is invalid."
        echo "The available emin values from the maps are: \"${eminv[*]}\""
        exit
    fi
    found=false
    for j in ${!emaxv[*]} ; do
        [ ${user_emaxv[i]} == ${emaxv[j]} ] && found=true && break
    done
    if [ ! "$found" = true ] ; then
        echo "Error with user defined energy ranges. The emax value \"${user_emaxv[i]}\" is invalid."
        echo "The available emax values from the maps are: \"${emaxv[*]}\""
        exit
    fi
done

echo "Energy ranges:"
echo "emin = ${eminv[*]}"
echo "emax = ${emaxv[*]}"

echo "User defined energy ranges:"
echo "emin = ${user_eminv[*]}"
echo "emax = ${user_emaxv[*]}"

for i in ${!eminv[*]} ; do
    naming="${eminv[i]}_${emaxv[i]}.${sky}.S${filter}_${matrix}.disp.conv.sky.gz"
    AG_add_diff5 diffusefilelist=$in edpfile=${edpfile} sarfile=${sarfile} outfile=\!$naming emin=${eminv[i]} emax=${emaxv[i]}
done

for i in ${!user_eminv[*]} ; do
    naming="${user_eminv[i]}_${user_emaxv[i]}.${sky}.S${filter}_${matrix}.disp.conv.sky.gz"
    AG_add_diff5 diffusefilelist=$in edpfile=${edpfile} sarfile=${sarfile} outfile=\!$naming emin=${user_eminv[i]} emax=${user_emaxv[i]}
done

echo "Done."
