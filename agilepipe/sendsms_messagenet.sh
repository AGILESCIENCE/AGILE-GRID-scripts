#!/bin/bash

if [ $# -ne 2 ] ; then
    echo "Error, missing arguments."
    echo "Usage: $0 \"Message to send\" \"italian telephone number\""
    exit 1
fi

id="andrea.bulgarelli@libero.it"
pass="gurumeu"

curl -X POST -F "redirect=https://www.messagenet.com/sms/invia/" -F "userid=$id" -F "Submit=Entra" -F "password=$pass" --cookie-jar /tmp/messagenetcookie https://www.messagenet.com/sms/invia/ &> /dev/null

curl -X POST --cookie /tmp/messagenetcookie -F "nazione_n=39" -F "numero_n=$2" -F "mittente=web" -F "messaggio=$1" -F"invia=Send" https://www.messagenet.com/sms/invia/ &> /dev/null

rm /tmp/messagenetcookie
