#!/bin/bash

fieldname=$1
fieldid=$2

for num in 0 1 2 3 4 5 6 7 8 9 
do
python Loop_z.py --nevts 1000 --season ${num} --fieldname ${fieldname} --fieldid ${fieldid}

done