#!/bin/bash
# @ shell = /bin/bash
# @ job_name = extended_analysis
# @ job_type = serial
# @ environment = COPY_ALL
# @ class = extended
# @ wall_clock_limit = 240:00:00
# @ resources = ConsumableCpus(1) 
# @ error = job1.$(jobid).err
# @ output = job1.$(jobid).out
# @ notify_user = bulgarelli@iasfbo.inaf.it
# @ queue

echo "Environment: ${SOURCE} ${L} ${B} ${R}"

if [[ -z "${SOURCE}" ]] ; then
    echo "Error, missing \$SOURCE."
    exit
fi
if [[ -z "${L}" ]] ; then
    echo "Error, missing \$L."
    exit
fi
if [[ -z "${B}" ]] ; then
    echo "Error, missing \$B."
    exit
fi
if [[ -z "${R}" ]] ; then
    echo "Error, missing \$R."
    exit
fi

analysisname="${SOURCE}-${R}"
archive="/data01/CAT2/CATALOG3072_ASDCe_B01_H0025_SKY002_FOVBINUMBER1_ENERGYBIN0/cat2edp5"
logfile="${archive}/run_extended.log"

if [[ -n ${LOADL_STEP_INITDIR} ]] ; then
    date >> $logfile
    echo "Source ${analysisname} submitted, errfile: $LOADL_STEP_INITDIR/$LOADL_STEP_ERR" >> $logfile
fi

module load agile-B25-r5
module load python2.7-sci
module load heasoft-6.17
#heainit
source /opt/sciencesoft/heasoft-6.17/x86_64-unknown-linux-gnu-libc2.12/headas-init.sh
export DISPLAY=:0
export HEADASNOQUERY=
export HEADASPROMPT=/dev/null
export PFILES=".;/opt/prod/agile-B25-r5/share:/opt/sciencesoft/heasoft-6.17/x86_64-unknown-linux-gnu-libc2.12/syspfiles"
echo $PFILES

date
ring=$(ruby distcenter.rb ${L} ${B} | cut -d ' ' -f1)
ringpath="${archive}/${ring}"
dist=$(ruby distcenter.rb ${L} ${B} | cut -d ' ' -f3)
echo "Ring: ${ringpath}" >> $logfile
echo "Ring: ${ring}" >> $logfile
echo "Distance: ${dist}" >> $logfile
echo "Source ${analysisname} Ring path: ${ringpath}" >> $logfile
echo "Source ${analysisname} Distance: ${dist}" >> $logfile

cd ${ringpath}
mkdir -p ${analysisname}
cd ${analysisname}

cp $AGILE/share/AG_circle.par .
AG_circle ../FM3.119_ASDCe_H0025_B01.exp.gz !${analysisname}.circle.gz ${L} ${B} ${R} >> commands.log 2>&1

echo "convolve........" >> commands.log 2>&1

$AGILE/scripts/extendedsources/convolve_template.sh ${analysisname}.circle.gz 2.1 >> commands.log 2>&1

echo "explist...." >> commands.log 2>&1

cp ${archive}/template_subdir.explist ${analysisname}.explist

### RUN1
### touch ${analysisname}.multi
### RUN2
touch ${analysisname}.multi

#extract_catalog.rb /ANALYSIS3/catalogs/cat2R.multi ${L} ${B} ${analysisname}.multi 0.1 1 5 0 7 0
#echo "0.0 ${L} ${B} 2.1 0 2 DUMMYFORIG 0.0" >> ${analysisname}.multi


echo "add_templates....." >> commands.log 2>&1
echo "${AGILE}/scripts/extendedsources/add_templates.sh ${analysisname}.circle.gz ${analysisname}.explist" >> commands.log 2>&1



$AGILE/scripts/extendedsources/add_templates.sh ${analysisname}.circle.gz ${analysisname}.explist >> commands.log 2>&1


echo "${analysisname} -2.0E-07 2.1 FM3.119_ASDCe_H0025_B01.exp.gz.template.gz " > ${analysisname}.multiext

sed "s/FM3/..\/FM3/g" ../FM3.119_ASDCe_H0025_B01.maplist4 > FM3.119_ASDCe_H0025_B01_100-50000.maplist4

### RUN1
### multi5.rb FM3.119_ASDCe_I0025 FM3.119_ASDCe_I0025_B01_100-50000.maplist4 ${analysisname}.multi ${analysisname}.res listsourceextended=${analysisname}.multiext >> commands.log 2>&1
### RUN2

echo "multi6.rb FM3.119_ASDCe_H0025 FM3.119_ASDCe_H0025_B01_100-50000.maplist4 ${analysisname}.multi ${analysisname}.res listsourceextended=${analysisname}.multiext " >> commands.log 2>&1

multi6.rb FM3.119_ASDCe_H0025 FM3.119_ASDCe_H0025_B01_100-50000.maplist4 ${analysisname}.multi ${analysisname}.res listsourceextended=${analysisname}.multiext  >> commands.log 2>&1

echo "Source ${analysisname} run completed." >> $logfile
