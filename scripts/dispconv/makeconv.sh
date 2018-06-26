#!/bin/bash

trap exit ERR

usage() {
    cat << EOF
Usage: $(basename $0) [ options ... ] <sky map> ...

  -h            print help
  -m MATRIX     matrix to use
  -f FILTER     filter to use
  -d PATH       model data path
  -i INDEX      the spectral index

  Default options: -m I0025 -f FMG -d \$AGILE/model/scientific_analysis/data -i 2.0
EOF
    exit $1
}

# Check AGILE environment
if [ -z "$AGILE" ] ; then
    echo "\$AGILE environment variable is not set, quitting.."
    exit 1
fi

# Parse options
filter="FMG"
matrix="I0025"
modelpath="$AGILE/model/scientific_analysis/data"
index="2.0"
while getopts "hm:f:d:i:" opt; do
    case "${opt}" in
        h) usage 0;;
        m) matrix="$OPTARG";;
        f) filter="$OPTARG";;
        d) modelpath="$OPTARG";;
        i) index="$OPTARG";;
        *) usage 1;;
    esac
done
shift $(($OPTIND-1))

# Check if model files exists
echo "Checking model files.."
modelnaming="/AG_GRID_G0017_S${filter}_${matrix}"
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

# Check sky maps naming convention
echo "Checking sky maps.."
for skymap in $@; do
    if [[ ! $(basename $skymap) =~ \.[0-9]\.[0-9].sky.gz && ! $(basename $skymap) =~ \.SKY[0-9][0-9][0-9].sky.gz ]]
    then
        echo "Error with sky map $skymap: doesn't match the the old and new naming conventions"
        echo "Expected a name <emin>_<emax>.<sky>.sky.gz"
        echo "The <sky> name convention is like 'SKY002' for new maps, '0.1' or '0.5' for old ones"
        exit
    fi
done

tmp=${skymap#*.}
sky=${tmp%.*.*}
naming="${sky}.S${filter}_${matrix}"
in="AG_add_diff.${naming}.in"
if [ ! -f ${in} ] ; then
    echo "Creating new .in file .."
    echo "$index 0" > $in
fi

for skymap in $@; do
    name=$(basename $skymap .${sky}.sky.gz)
    outf="${name}.${naming}.conv.sky.gz"

    [ -f $outf ] && echo "Removing existing ${outf}.." &&  rm $outf

    echo "Convolving $outf .."
    AG_diff_conv5 diffusefile=$skymap sarfile=$sarfile psdfile=$psdfile edpfile=$edpfile outfile=$outf

    [ -z $(grep $outf $in) ] && echo "Adding ${outf} to ${in}.." && echo "$(readlink -m $outf)" >> $in
done

# Update the header of the .in file
let n=$(wc -l $in | cut -d' ' -f1)-1
sed -i "1 s/.*/$index $n/" $in

echo "Done."
