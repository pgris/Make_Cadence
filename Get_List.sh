#!/bin/bash

fieldname=$1
fieldid=$2
season=$3
nevts=$4
sntype=$5

output=List_${sntype}_${fieldname}_${fieldid}_${nevts}_season_${season}_minion.dat

thedir=Sim_minion_1016

ls ${thedir}/SuperNova_${sntype}_${fieldname}_${fieldid}_*_${nevts}_season_${season}_*.pkl | cut -d '/' -f2 >& ${output}

