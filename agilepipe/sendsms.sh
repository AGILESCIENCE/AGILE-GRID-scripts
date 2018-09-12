#!/bin/bash

if [ $# -ne 2 ] ; then
    echo "Error, missing arguments."
    echo "Usage: $0 \"Message to send\" \"num1,num2,num3\""
    exit 1
fi

ssh smsextalert@giacal1.giano.iasfbo "echo $1 | /usr/local/bin/sendsms.pl $2"
