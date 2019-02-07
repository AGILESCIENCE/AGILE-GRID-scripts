source /opt/module/heasoft-6.25
source /opt/module/py_manual_27
source /opt/module/agile-B25-r5

cp $AGILE/share/AG_circle.par .
AG_circle gamma_cygni_point.exp.gz gamma_cygni_point.circle.gz 78.18 2.11 0.1

./convolve_template_single.sh gamma_cygni_point.circle.gz 2.1 00100 00400 I0025

echo "gamma_cygni_point.exp.gz -2.1" > gamma_cygni_point.explist

touch gamma_cygni_point.multi

echo "add template"
cp $AGILE/share/AG_add_diff.par .
cp $AGILE/share/AG_gasmapgen.par .
./add_templates_single.sh gamma_cygni_point.circle.gz gamma_cygni_point.explist 00100 00400 I0025

echo "gammacygni -2.0E-07 2.1 gamma_cygni_point.exp.gz.template.gz " > gamma_cygni_point.multiext

source /opt/module/agile-B25-r5
multi.rb FM3.119_ASDCe_I0025 gamma_cygni_point.maplist4 gamma_cygni_point.multi gamma_cygni_point.res listsourceextended=gamma_cygni_point.multiext
