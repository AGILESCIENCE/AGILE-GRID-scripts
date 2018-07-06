echo "install scripts"
test -d $AGILE/scripts || mkdir -p $AGILE/scripts
cp -rf scripts/* $AGILE/scripts

echo "install catalogs"
test -d $AGILE/catalogs || mkdir -p $AGILE/catalogs
cp -rf catalogs/* $AGILE/catalogs

echo "install spot6"
test -d $AGILE/AGILEPIPE/spot6 || mkdir -p $AGILE/AGILEPIPE/spot6
cp -rf spot6/* $AGILE/AGILEPIPE/spot6
