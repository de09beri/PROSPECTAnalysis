#!/bin/bash

# ./RnPoPlot.sh $1 $2 
# $1: Run PlotRnPoVsCell.C
# $2: Run PlotRnPoVsTime.C

#==============================================
# Run RnPoVsCell.C 
echo ======= Running PlotRnPoVsCell =======

if [ $1 -eq 1 ]
then

root -l -b <<EOF
.L Plot/PlotRnPoVsCell.C+
PlotRnPoVsCell()
.q
EOF

fi

#==============================================
# Run RnPoVsTime.C 
echo ======= Running PlotRnPoVsTime =======

if [ $2 -eq 1 ]
then

root -l -b <<EOF
.L Plot/PlotRnPoVsTime.C+
PlotRnPoVsTime()
.q
EOF

fi

