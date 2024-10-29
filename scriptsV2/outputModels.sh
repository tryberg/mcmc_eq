#!/bin/bash

# write out the final models for use elsewhere

inp=resmcnx.dat
ls $inp
test -f $inp || echo "did you run dispe yet?"
test -f $inp || exit

# P-wave std model 1:
awk '{if ($1=="STAN") print $3-$4, $2}' resmcnx.dat |\
 awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' > sd1_p.mod

# P-wave std model 1:
awk '{if ($1=="STAN") print $3+$4, $2}' resmcnx.dat |\
 awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' > sd2_p.mod

# P-wave models:
awk '{if ($1=="STAN") print $3, $2}' resmcnx.dat |\
 awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' > newAvg_p.mod


#awk '{if ($1=="STAN") print $3, $2}' resmcns.dat |\
# awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' > min_p.mod

awk '{if ($1=="STAN") print $11, $2}' resmcnx.dat |\
 awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' > maxProb_p.mod

ls *p.mod
############ S-wave

# S-wave std model 1:
awk '{if ($1=="STAN") {f=$3/$5; x1=$4/$3; x2=$6/$5; df=f*sqrt(x1*x1+x2*x2); print f+df, $2}}' resmcnx.dat |\
 awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' > sd1_s.mod

# S-wave std model 1:
awk '{if ($1=="STAN") {f=$3/$5; x1=$4/$3; x2=$6/$5; df=f*sqrt(x1*x1+x2*x2); print f-df, $2}}' resmcnx.dat |\
 awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' > sd2_s.mod

#awk '{if ($1=="STAN") print $3/$5, $2}' resmcnx.dat  |\
# awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' > classicAvg_s.mod
#psxy classicAvg_s.mod -JX -R -W5/255/200/200 -O -K -N  >> $output

awk '{if ($1=="STAN") print $3/$5, $2}' resmcnx.dat  |\
 awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' > newAvg_s.mod

#awk '{if ($1=="STAN") print $3/$5, $2}' resmcns.dat |\
# awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' > min_s.mod

awk '{if ($1=="STAN") print $11/$12, $2}' resmcnx.dat |\
 awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' > maxProb_s.mod

ls *s.mod

###### write out reloc file for other plots
egrep EZ resmcnx.dat | awk '{printf "%8.3f %8.3f %8.3f %015.3f %7.3f %f %f %f\n",$3,$4,$5,$9,$10,$6,$7,$8}' > eqs.reloc.xyz
ls eqs.reloc.xyz
