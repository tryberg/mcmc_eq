#!/bin/bash

# write out the final models for use elsewhere

inp=resmcnx.dat
ls $inp
test -f $inp || echo "did you run disp_m_ yet?"
test -f $inp || exit

# P-wave models: classic mean
awk '{if ($1=="STAN") print $3, $2}' $inp |\
 awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' > newAvg_p.mod

# P-wave std model 1: classic mean
awk '{if ($1=="STAN") print $3-$4, $2}' $inp |\
 awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' > sd1_p.mod

# P-wave std model 2: classic mean
awk '{if ($1=="STAN") print $3+$4, $2}' $inp |\
 awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' > sd2_p.mod

# P-wave models: excess mean
awk '{if ($1=="STAN") print $7, $2}' $inp |\
 awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' > newAvg_p2.mod

# P-wave std model 1: excess mean
awk '{if ($1=="STAN") print $7-$8, $2}' $inp |\
 awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' > sd1_p2.mod

# P-wave std model 2: excess mean
awk '{if ($1=="STAN") print $7+$8, $2}' $inp |\
 awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' > sd2_p2.mod


#awk '{if ($1=="STAN") print $3, $2}' resmcns.dat |\
# awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' > min_p.mod

awk '{if ($1=="STAN") print $11, $2}' $inp |\
 awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' > maxProb_p.mod

ls *p.mod
############ S-wave
# S-wave models: classic mean
awk '{if ($1=="STAN") print $3/$5, $2}' $inp  |\
 awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' > newAvg_s.mod

# S-wave std model 1: classic mean
awk '{if ($1=="STAN") {f=$3/$5; x1=$4/$3; x2=$6/$5; df=f*sqrt(x1*x1+x2*x2); print f+df, $2}}' $inp |\
 awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' > sd1_s.mod

# S-wave std model 2: classic mean
awk '{if ($1=="STAN") {f=$3/$5; x1=$4/$3; x2=$6/$5; df=f*sqrt(x1*x1+x2*x2); print f-df, $2}}' $inp |\
 awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' > sd2_s.mod

# S-wave models: excess mean
awk '{if ($1=="STAN") print $7/$9, $2}' $inp  |\
 awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' > newAvg_s2mod

# S-wave std model 1: excess mean
awk '{if ($1=="STAN") {f=$7/$9; x1=$8/$7; x2=$10/$9; df=f*sqrt(x1*x1+x2*x2); print f+df, $2}}' $inp |\
 awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' > sd1_s2.mod

# S-wave std model 2: excess mean
awk '{if ($1=="STAN") {f=$7/$9; x1=$8/$7; x2=$10/$9; df=f*sqrt(x1*x1+x2*x2); print f-df, $2}}' $inp |\
 awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' > sd2_s2.mod

#awk '{if ($1=="STAN") print $3/$5, $2}' $inp  |\
# awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' > classicAvg_s.mod
#psxy classicAvg_s.mod -JX -R -W5/255/200/200 -O -K -N  >> $output

#awk '{if ($1=="STAN") print $3/$5, $2}' resmcns.dat |\
# awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' > min_s.mod

awk '{if ($1=="STAN") print $11/$12, $2}' $inp |\
 awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' > maxProb_s.mod

ls *s.mod

###### write out reloc file for other plots
# output new catalog: {'X','Y','Z','OT','dOT','ex','ey','ez'}
cat="quakes_mcmc.dat"
echo "# {'EVID','X','Y','Z','OT','dOT','ex','ey','ez'} (%8.3f %8.3f %8.3f %015.3f %7.3f %f %f %f\n)" > "$cat"
egrep EZ $inp | awk '{printf "%03d %8.3f %8.3f %8.3f %015.3f %7.3f %f %f %f\n",$2,$3,$4,$5,$9,$10,$6,$7,$8}' >> "$cat"
echo "new output catalog saved:"
ls "$PWD/$cat"

##
picks=picks.mcmc
ls $picks
test -f $picks || exit

awk '{if ($1=="RES") print $2, $3, $4}' "$inp" > t1
awk '{if ($1!="#" && $2!="NA") print $1, $2, $4, $5, $6}' "$picks" | sort | uniq > rec.dat
paste t1 rec.dat > recdata
stac=staCors_mcmc.dat
awk '{print $4,$2,$3}' recdata > "$stac"
echo "station corrections saved:"
ls "$PWD/$stac"

ls stations.dat
test -f stations.dat || exit
paste stations.dat staCors_mcmc.dat |\
awk '{printf "%5s %3d %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",$10,$1,$2,$3,$4,$11,$12,$7,$8}' > stations.out
ls stations.out
