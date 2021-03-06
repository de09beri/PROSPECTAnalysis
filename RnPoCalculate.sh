#!/bin/bash

# ./RnPoCalculate.sh $1 $2 $3 $4 $5
# $1: Run MakeClass
# $2: Run RnPoVsCell.C
# $3: Run RnPoVsTime.C
# $4: Run RnPoCellVsTime.C
# $5: Run RnPoColVsTime.C and RnPoRowVsTime.C


promptPSDStdDev=4.0
delayPSDStdDev=4.0
promptEnStdDev=4.0
delayEnStdDev=4.0
dzStdDev=4.0

zLow=-1000
zHigh=1000

dtCut=0.5
timeBin=47.5
#timeBin=1000000000

boolESmear=false

output=Top

#==============================================
# Make Class
echo ========== Making TAc Class ==========

if [ $1 -eq 1 ]
then

root -l -b <<EOF
.L Calculate/MakeAcTreeClass.C
MakeAcTreeClass()
.q
EOF

mv RNPO* Calculate/
fi

#==============================================
# Make list of cut parameters per cell
echo ========= Running CutsPerCell =========

if [ $2 -eq 1 ]
then

root -l -b <<EOF 
.L Calculate/CutsPerCell.C
CutsPerCell($dtCut)
.q
EOF

fi

#==============================================
# Run RnPoVsCell
echo ========= Running RnPoVsCell =========

if [ $3 -eq 1 ]
then
root -l -b <<EOF
string outputName = "$output"

.L Calculate/RnPoVsCell.C
RnPoVsCell($promptPSDStdDev, $delayPSDStdDev, $promptEnStdDev, $delayEnStdDev, $dzStdDev, $zLow, $zHigh, $dtCut, $boolESmear, outputName)
.q
EOF

fi

#==============================================
# Run RnPoVsTime
echo ========= Running RnPoVsCell_Calculations =========

if [ $4 -eq 1 ]
then

root -l -b <<EOF 
string outputName = "$output"

.L Calculate/RnPoVsCell_Calculations.C
RnPoVsCell_Calc(outputName)
.q
EOF

fi

#==============================================
# Make list of cut parameters versus time
echo ========= Running CutsVsTime =========

if [ $5 -eq 1 ]
then

root -l -b <<EOF 
.L Calculate/CutsVsTime.C
CutsVsTime($dtCut, $timeBin)
.q
EOF

fi

#==============================================
# Run RnPoVsTime
echo ========= Running RnPoVsTime =========

if [ $6 -eq 1 ]
then

root -l -b <<EOF 
string outputName = "$output"
.L Calculate/RnPoVsTime.C
RnPoVsTime($promptPSDStdDev, $delayPSDStdDev, $promptEnStdDev, $delayEnStdDev, $dzStdDev, $zLow, $zHigh, $timeBin, $dtCut, $boolESmear, outputName)
.q
EOF

fi

#==============================================
# Run RnPoVsTime
echo ========= Running RnPoVsTime_Calculations =========

if [ $7 -eq 1 ]
then

root -l -b <<EOF 
string outputName = "$output"
.L Calculate/RnPoVsTime_Calculations.C
RnPoVsTime_Calc(outputName)
.q
EOF

fi
