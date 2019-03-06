#!/bin/bash

# ./RnPoCalculate.sh $1 $2 $3 $4 $5
# $1: Run MakeClass
# $2: Run RnPoVsCell.C
# $3: Run RnPoVsTime.C
# $4: Run RnPoCellVsTime.C
# $5: Run RnPoColVsTime.C and RnPoRowVsTime.C


promptPSDStdDev=3.5
delayPSDStdDev=3.5
promptEnStdDev=2.7
delayEnStdDev=2.0
dzStdDev=3.5

zLow=-1000
zHigh=1000

dtCut=0.0
timeBin=23.5

boolESmear=false

#==============================================
# Make Class
echo ========== Making TAc Class ==========

if [ $1 -eq 1 ]
then

root -l -b <<EOF
.L Calculate/MakeAcTreeClass.C+
MakeAcTreeClass()
.q
EOF

mv RNPO* Calculate/
fi

#==============================================
# Run RnPoVsCell
echo ========= Running RnPoVsCell =========

if [ $2 -eq 1 ]
then
root -l -b <<EOF 
.L Calculate/RnPoVsCell_Modified.C+
RnPoVsCell($zLow, $zHigh, $dtCut)
.q
EOF

fi

#==============================================
# Run RnPoVsTime
echo ========= Running RnPoVsTime =========

if [ $3 -eq 1 ]
then

root -l -b <<EOF 
.L Calculate/RnPoVsTime_Modified.C+
RnPoVsTime($zLow, $zHigh, $timeBin, $dtCut)
.q
EOF

fi


