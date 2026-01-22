#!/bin/bash

cfg=config_eqx.dat
pick_file=picks.mcmc # can use output of pha2mcmc here
prior_model=model.inp # must have # of layers = nz in cfg

ls $cfg
test -f $cfg || exit
ls $pick_file
test -f $pick_file || exit
#ls quakes.dat
#test -f quakes.dat || exit
#ls stations.dat
#test -f stations.dat || exit


d=$(awk 'NR==1 {print $1}' "$cfg")
nz=$(awk 'NR==4 {print $1}' "$cfg")
zmin=$(awk 'NR==7 {print $1}' "$cfg")
zmax=$(echo "$zmin $nz $d" | awk '{print $1+$2*$3-$3}')

echo $zmin $zmax $d $nz
getGeneric1Dmodel.sh $zmin $zmax $d | awk '{print $2,$1,1.730}' > $prior_model
ls $prior_model
test -f $prior_model || exit

nl=`cat $prior_model | wc -l`

[[ $nz -eq $nl ]] || { echo "Error: $nz != $nl"; exit 1; }

# construct input file
awk '{print "STAN", $1, $2, 0, $3, 0, $2, 0, $3, 0, $2, $3, 0.01}' $prior_model > res.dat
# do these matter or can I put zeros here?
awk '{print "EQ",$1, $2, $3, $4,0,0,0,0,0,0,0}' quakes.dat >> res.dat
awk '{print "EZ",$1, $2, $3, $4,0,0,0,0,0,0,0}' quakes.dat >> res.dat
# to load prior corrections only, not needed? put zeros instead?
awk '{print "RES", $1, $5, $6, 0, 0}' stations.dat >> res.dat
echo "NOISE 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1"  >> res.dat
mv res.dat model.dat
ls model.dat
