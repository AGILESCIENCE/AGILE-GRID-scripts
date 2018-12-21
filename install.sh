if [ -z "$AGILE" ] || [ -z $(env | grep "AGILE=") ] ; then
    echo "AGILE environment variable not set. Abort."
    exit
fi

echo "install scripts"
test -d $AGILE/scripts || mkdir -p $AGILE/scripts
cp -rf scripts/* $AGILE/scripts

echo "install catalogs"
test -d $AGILE/catalogs || mkdir -p $AGILE/catalogs
cp -rf catalogs/* $AGILE/catalogs



if [ -z "$AGILEPIPE" ] || [ -z $(env | grep "AGILEPIPE=") ] ; then
    echo "AGILEPIPE environment variable not set. Abort."
    exit
fi

echo "install spot6"
test -d $AGILEPIPE/spot6 || mkdir -p $AGILEPIPE/spot6
cp -rf spot6/* $AGILEPIPE/spot6
cp -rf spot6/env.rb $AGILEPIPE
cp -rf spot6/template.ll $AGILEPIPE

