#!/bin/bash

# ./RnPoCalculate.sh $1 $2 $3 $4 $5
# $1: Run MakeClass
# $2: Run RnPoVsCell.C
# $3: Run RnPoVsTime.C
# $4: Run RnPoCellVsTime.C
# $5: Run RnPoColVsTime.C and RnPoRowVsTime.C


promptPSDStdDev=3.5
delayPSDStdDev=3.5
#promptEnStdDev=2.7
#delayEnStdDev=2.0
promptEnStdDev=3.5
delayEnStdDev=3.5
dzStdDev=3.5

zLow=-444
zHigh=444

dtCut=0.5
#timeBin=47.5
timeBin=1000000
timeBinCell=335.5

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
# Make list of cut parameters per cell
echo ========= Running CutsPerCell =========

if [ $2 -eq 1 ]
then

root -l -b <<EOF 
.L Calculate/CutsPerCell.C+
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
.L Calculate/RnPoVsCell.C+
RnPoVsCell($promptPSDStdDev, $delayPSDStdDev, $promptEnStdDev, $delayEnStdDev, $dzStdDev, $zLow, $zHigh, $dtCut, $boolESmear)
.q
EOF

fi

#==============================================
# Make list of cut parameters versus time
echo ========= Running CutsVsTime =========

if [ $4 -eq 1 ]
then

root -l -b <<EOF 
.L Calculate/CutsVsTime.C+
CutsVsTime($dtCut, $timeBin)
.q
EOF

fi

#==============================================
# Run RnPoVsTime
echo ========= Running RnPoVsTime =========

if [ $5 -eq 1 ]
then

root -l -b <<EOF 
.L Calculate/RnPoVsTime.C+
RnPoVsTime($promptPSDStdDev, $delayPSDStdDev, $promptEnStdDev, $delayEnStdDev, $dzStdDev, $zLow, $zHigh, $timeBin, $dtCut, $boolESmear)
.q
EOF

fi

#==============================================
# Run RnPoColVsTime and RnPoRunVsTime
echo ========= Running RnPoCellVsTime =========

if [ $6 -eq 1 ]
then

root -l -b <<EOF 
.L Calculate/RnPoCellVsTime.C+
RnPoCellVsTime($promptPSDStdDev, $delayPSDStdDev, $promptEnStdDev, $delayEnStdDev, $dzStdDev, $zLow, $zHigh, $timeBinCell, $dtCut, $boolESmear)
.q
EOF

fi


#==============================================
# Run RnPoColVsTime and RnPoRunVsTime

#if [ $7 -eq 1 ]
#then

#echo ======= Running RnPoColVsTime =======

#root -l -b <<EOF 
#.L Calculate/RnPoColVsTime.C+
#RnPoColVsTime($p_lowPSD,$d_lowPSD,$p_lowE,$d_lowE,$timeBin,$dtFit)
#.q
#EOF

#echo ======= Running RnPoRowVsTime =======

#root -l -b <<EOF 
#.L Calculate/RnPoRowVsTime.C+
#RnPoRowVsTime($p_lowPSD,$d_lowPSD,$p_lowE,$d_lowE,$timeBin,$dtFit)
#.q
#EOF

#fi

