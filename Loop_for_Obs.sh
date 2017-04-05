#!/bin/bash

fieldname=$1
thefile=../Ana_Cadence/fieldIDs_minion_1016_${fieldname}.txt
list=`cat ${thefile} | grep fieldIds | cut -d ' ' -f2`

for file in $list
do
echo ${file}
python Loop_z.py --fieldname ${fieldname} --fieldid ${file} --runtype Observation --dbFile ../Fake_Rolling/Rolling.db 
done