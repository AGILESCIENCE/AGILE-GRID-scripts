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
  -u EMAIL      user email to send jobs completition report.

  Default options: -m I0025 -f FMG -d \$AGILE/model/scientific_analysis/data -i 2.0
EOF
    exit $1
}

# Load the agile module
modulename=agile-B23

module load $modulename

# Parse options
filter="FMG"
matrix="I0025"
modelpath="$AGILE/model/scientific_analysis/data"
index="2.0"
user=
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

for skymap in $@; do
    fullpath=$(readlink -m $skymap)
    filename=$(basename $skymap)

    jobfile="#!/bin/bash
# @ shell =   /bin/bash
# @ job_name    = makeconv_$filename
# @ job_type    = serial
# @ class       = large
# @ resources   = ConsumableCpus(1)
# @ wall_clock_limit = 240:00:00
# @ environment = COPY_ALL
# @ error       = conv_$filename.err
# @ output      = conv_$filename.out
# @ notify_user = $user
# @ queue

date
module load $modulename
echo \"Input file: $skymap\"
cd $PWD
$(dirname $0)/makeconv.sh -f $filter -m $matrix -d \"$modelpath\" $fullpath
date
"
	echo "$jobfile" | llsubmit -

done
