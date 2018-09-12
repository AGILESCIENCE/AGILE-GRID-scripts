echo "install scripts"
test -d $AGILE/scripts || mkdir -p $AGILE/scripts
cp -rf scripts/* $AGILE/scripts

echo "install catalogs"
test -d $AGILE/catalogs || mkdir -p $AGILE/catalogs
cp -rf catalogs/* $AGILE/catalogs

echo "install spot6"
test -d $AGILE/AGILEPIPE/spot6 || mkdir -p $AGILE/AGILEPIPE/spot6
cp -rf spot6/* $AGILE/AGILEPIPE/spot6
cp -rf spot6/env.rb $AGILE/AGILEPIPE
cp -rf spot6/template.ll $AGILE/AGILEPIPE

echo "install agilepipe"
test -d /opt/prod/AGILEPIPE || mkdir -p /opt/prod/AGILEPIPE
cp -rf agilepipe/* /opt/prod/AGILEPIPE
