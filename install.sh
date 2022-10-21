echo "install scripts"
test -d $AGILE/scripts || mkdir -p $AGILE/scripts
cp -rf scripts/* $AGILE/scripts

echo "install catalogs"
test -d $AGILE/catalogs || mkdir -p $AGILE/catalogs
cp -rf catalogs/* $AGILE/catalogs

