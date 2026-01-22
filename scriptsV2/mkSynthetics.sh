#!/bin/bash

cfg=config_fwd.dat
pick_file=picks # can use output of pha2mcmc here
synth_model=model.inp # must have # of layers = nz in cfg

ls $cfg
test -f $cfg || exit
ls $pick_file
test -f $pick_file || exit
ls $synth_model
test -f $synth_model || exit
ls quakes.dat
test -f quakes.dat || exit
ls stations.dat
test -f stations.dat || exit

nz=`sed -n "4p" $cfg | awk '{print $1}'`
nl=`cat $synth_model | wc -l`

[[ $nz -eq $nl ]] || { echo "Error: $nz != $nl"; exit 1; }

# construct input file
awk '{print "STAN", $1, $2, 0, $3, 0, $2, 0, $3, 0, $2, $3, 0.01}' $synth_model > res.dat
awk '{print "EQ",$1, $2, $3, $4,0,0,0,0,0,0,0}' quakes.dat >> res.dat
awk '{print "EZ",$1, $2, $3, $4,0,0,0,0,0,0,0}' quakes.dat >> res.dat
awk '{print "RES", $1, $5, $6, 0, 0}' stations.dat >> res.dat
echo "NOISE 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1"  >> res.dat

# forward
fw $cfg res.dat $pick_file > predictions

paste $pick_file predictions | awk '{if ($1=="#") {print $1, $2, $3, $4, $5;} else {print $1, $2, $3, $4, $5, $6, $14, $8}}' > synths_wo_noise

# add gaussian noise: P0 0.05 P1 0.1 P2 0.15 P3 0.2 S0 0.175 S1 0.225 S2 0.275 S3 0.325
#rms=0.10
#JP: reduce noise level for volcanoes (see also disp_compare.sh)
# add gaussian noise: P0 0.015 P1 0.030 P2 0.045 P3 0.060 S0 0.0525 S1 0.0675 S2 0.0825 S3 0.0975
rms=0.03
rseed=33
rseed=$RANDOM
awk 'BEGIN {srand(1*"'$rseed'")}{do {x1=2.0*rand()-1; x2=2.0*rand()-1; s=x1*x1+x2*x2;} while ((s>1) || (s==0)); rnd=x1*sqrt(-2.0*log(s)/s); if ($1=="#") {print $0;} else {f=0; if ($3=="S") {f=1;} print $1, $2, $3, $4, $5, $6, $7+rnd*"'$rms'"*(($8+1+2.5*f)/4)*2, $8}}' synths_wo_noise >  synths_with_noise

# P and S equal weight:
#awk 'BEGIN {srand(1*"'$rseed'")}{do {x1=2.0*rand()-1; x2=2.0*rand()-1; s=x1*x1+x2*x2;} while ((s>1) || (s==0)); rnd=x1*sqrt(-2.0*log(s)/s); if ($1=="#") {print $0;} else {print $1, $2, $3, $4, $5, $6, $7+rnd*"'$rms'"*(($8+1)/4)*2, $8}}' synths_wo_noise >  synths_with_noise

# check
fw $cfg res.dat synths_wo_noise > t1
fw $cfg res.dat synths_with_noise > t1
#cp synths_with_noise picks
awk '{if ($1=="#") print $0; else printf "%4s %03d %1s %8.3f %8.3f %8.3f %8.3f %d\n",$1,$2,$3,$4,$5,$6,$7,$8}' synths_with_noise > picks.mcmc

# clean-up
#rm synths_with_noise synths_wo_noise 
#rm res.dat t1 predictions

