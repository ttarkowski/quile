#!/bin/bash

if test -f "$1"; then

    cat $1 | fold > tmp.txt
    mv tmp.txt $1

    lines=$(cat $1 | wc -l)
    max=50

    if [ ${lines} -gt ${max} ]; then
	tail -n ${max} $1 > tmp.txt
	(echo "[Output to this point was skipped.]" ; cat tmp.txt) > $1
	rm tmp.txt
    fi

fi

exit
