#!/bin/bash

# V2 converted from csh to bash by chatcpt and J. Pesicek, Oct 29, 2024
# ported to gmt6 Oct 30, 2024

rm gmt.*
gmt set MEASURE_UNIT INCH
gmt set FONT_ANNOT_PRIMARY 8 
gmt set HEADER_FONT_SIZE 8
gmt set HEADER_OFFSET 0.5c
gmt set LABEL_FONT_SIZE 8
gmt set COLOR_NAN 0/0/100
export LC_NUMERIC=C.UTF-8

cfg=config_eqx.dat
ls "$cfg"
[ -f "$cfg" ] || exit

exe=analyse_eq
which "$exe"
[ -z "$(whereis -q -b "$exe")" ] && echo "can't find $exe" && exit

eq=$(awk '{if (NR==30) print $1}' "$cfg") # needs to be past burn-in phase!
vv=$(awk '{if (NR==30) print $2}' "$cfg")
tot=$(echo "$vv $eq" | awk '{print $1+$2}')

str="need (1) start model input number post burn-in, b/t $eq and $tot, (2) % of models (90% default) to include, and (3) % of models close to boundaries to exclude (3% default)"
echo reading "$#" inputs: "$1" "$2" "$3"
[ "$#" -eq 3 ] || echo "$str"
[ "$#" -eq 3 ] || exit

if (( $1 > $eq && $1 < $tot )); then
    bi=$1
else
    echo "error: bad 1st input value"
    echo "$str"
    exit
fi

if (( $2 >= 0 && $2 <= 100 )); then
    prms=$2
else
    echo "error: bad 2nd input value (0-100%)"
    echo "$str"
    exit
fi

if (( $3 >= 0 && $3 <= 100 )); then
    p=$3
else
    echo "error: bad 3rd input value (0-100%)"
    echo "$str"
    exit
fi

output="xselect_${bi}_${prms}_${p}"
################################################################################################
# set p=3.0 # percentage of range value to detect proximity to technical boundary
# set prms=90 # percentage of RMS
# set bi=200000

f=1
awk -v f="$f" '{if (NR==1) print $1/f; if ((NR==2) || (NR==3) || (NR==4)) print $1*f; if (NR>4) print $0}' "$cfg" > configa.dat

d=$(awk '{if (NR==1) print $1}' configa.dat)
nx=$(awk '{if (NR==4) print $1}' configa.dat)
zmin=$(awk '{if (NR==7) print $1}' configa.dat)
zmax=$(echo "$zmin $nx $d" | awk '{print $1+$2*$3}')

vmin=$(awk '{if (NR==9) print $1}' configa.dat)
vmax=$(awk '{if (NR==10) print $1}' configa.dat)

v1=$(echo "$vmax $vmin" | awk -v p="$p" '{print $2+($1-$2)*p/100.0}')
v2=$(echo "$vmax $vmin" | awk -v p="$p" '{print $1-($1-$2)*p/100.0}')

vmins=$(awk '{if (NR==11) print $1}' configa.dat)
vmaxs=$(awk '{if (NR==12) print $1}' configa.dat)
vs1=$(echo "$vmaxs $vmins" | awk -v p="$p" '{print $2+($1-$2)*p/100.0}')
vs2=$(echo "$vmaxs $vmins" | awk -v p="$p" '{print $1-($1-$2)*p/100.0}')

vminvs=0.25
vmaxvs=6.93

touch tmp; rm tmp; touch tmp

for n in rjx-*.out; do
  echo "$n"
  awk -v bi="$bi" '($3 > 1 * bi) && ($3 < 300000) {print $0}' "$n" | 
  awk '$2 != "BF" {print $0}' | 
  awk '$1 != "cnt" {print $0}' >> tmp
done

# Remove models which are close to technical boundaries
awk '{if ($1=="mod") {f=0; for (i=0; i<$4; i++) {b=14+i*3+1; c=14+i*3+2; if ($b<1*"'$v1'") {f=f+1;} if ($b>1*"'$v2'") {f=f+1;} if ($c<1*"'$vs1'") {f=f+1;} if ($c>1*"'$vs2'") {f=f+1;}}} if (f==0) print $0}' tmp > tmpy

tmax=$(awk '{print $5}' tmpy | gmt pshistogram -R0/20/0/110 -W0.00001 -Z1 -Q -G0 -K -O -V -Io | 
       awk -v prms="$prms" '$2 >= 1 * prms {printf "%f\n", $1}' | head -1)

awk -v tmax="$tmax" '$5 < (1.0 * tmax) {print $0}' tmpy > tmpx

# Count post burn-in models
all=$(egrep mod tmp | wc -l)
allp=$(egrep mod tmpy | wc -l)
allr=$(egrep mod tmpx | wc -l)

dv=0.01
dvpvs=0.01

$exe configa.dat tmpx "$dv" "$dvpvs" > resmcnx.dat

# ------ PPPPPPP
awk '$1 == "BINP" {print $2, $3, $4}' resmcnx.dat | 
awk '$3 > 0 {print $0}' | 
gmt xyz2grd -I"$dv/$d" -R"$vmin/$vmax/$zmin/$zmax" -Ghistp.grd -V

gmt grd2cpt -Cseis -Z -D -I histp.grd > colv.cpt

# Model
gmt psbasemap -JX1.8/-6 -R"$vmin/$vmax/$zmin/$zmax" -B1f0.5:"Vp":/10f5g10000SeWN -K -Y1.6 -X0.4 > "$output.ps"
gmt grdimage histp.grd -JX -R -Ccolv.cpt -K -O >> "$output.ps"
echo "$vmin" "$zmax" "$bi" "$p" "$prms" "$all" "$allp" "$allr" "$(pwd)" | 
awk '{print $1, $2 * 1.15, 10, 0, 0, 0, "BI", $3, "P", $4, "PRMS", $5, $6, "/", $7, "/", $8, $9}' | 
gmt pstext -JX -R -K -O -N -F+jBL >> "$output.ps"

# hist count
gmt psbasemap -JX -R -B1f0.5:"Vp":/10f5g10000SewN -K -X2.05 -O >> "${output}.ps"

# modified
awk '{if ($1=="STAN") print $7, $2}' resmcnx.dat | awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' | gmt psxy -JX -R -W1.5,red -O -K -N >> "${output}.ps"
awk '{if ($1=="STAN") print $7-$8, $2}' resmcnx.dat | awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' | gmt psxy -JX -R -W,gray -O -K -N >> "${output}.ps"
awk '{if ($1=="STAN") print $7+$8, $2}' resmcnx.dat | awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' | gmt psxy -JX -R -W,gray -O -K -N >> "${output}.ps"

# boundary
max=$(awk '{if ($1=="STAN") print $13*1.1}' resmcnx.dat | sort -n | tail -1 | awk '{x=$1; if ($1==0) x=1; print x}')
# set max = 1 (Uncomment the following line if you want to force max to 1)
# max=1
gmt psbasemap -JX1/-6 -R0/$max/$zmin/$zmax -B100f5:"Boundary":/10f5g10000SewN -K -X2.05 -O >> "${output}.ps"
awk '{if ($1=="STAN") print $13, $2}' resmcnx.dat | awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' | gmt psxy -JX -R -W0.75,black -O -K -N >> "${output}.ps"

# ------- SSSSSSS
awk '{if ($1=="BINV") print $2, $3, $4}' resmcnx.dat | awk '{if ($3>0) print $0}' | gmt xyz2grd -I$dvpvs/$d -R$vmins/$vmaxs/$zmin/$zmax -Ghists.grd -V
# hmax=$(grdinfo hists.grd | egrep z_max | awk '{print $5}')

# gmt makecpt -Cseis -Z -T0/$hmax/1 -D -I >! colv.cpt
gmt grd2cpt -Cseis -Z -D -I histp.grd > colv.cpt

gmt psbasemap -JX1.8/-6 -R$vmins/$vmaxs/$zmin/$zmax -B1f0.5:"Vp/Vs":/10f5g1000SewN -K -O -X1.2 >> "${output}.ps"
gmt grdimage hists.grd -JX -R -Ccolv.cpt -K -O >> "${output}.ps"

# hist count
gmt psbasemap -JX -R -B1f0.5g1.732:"Vp/Vs":/10f5g1000SewN -K -X2.05 -O >> "${output}.ps"

# modified

awk '{if ($1=="STAN") print $9, $2}' resmcnx.dat | awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' | gmt psxy -JX -R -W1.5,red -O -K -N >> "${output}.ps"
awk '{if ($1=="STAN") print $9-$10, $2}' resmcnx.dat | awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' | gmt psxy -JX -R -W,gray -O -K -N >> "${output}.ps"
awk '{if ($1=="STAN") print $9+$10, $2}' resmcnx.dat | awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' | gmt psxy -JX -R -W,gray -O -K -N >> "${output}.ps"
#awk '{if ($1=="STAN") print $5, $2}' resmcnbf.dat | awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' | psxy -JX -R -W5/blue -O -K -N >> "${output}.ps"

gmt psbasemap -JX1.5/-6 -R$vminvs/$vmaxvs/$zmin/$zmax -B1f0.5:"Vs":/10f5g10000SEwN -K -O -X2.05 >> "${output}.ps"

awk '{if ($1=="STAN") print $7/$9, $2}' resmcnx.dat | awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' | gmt psxy -JX -R -W1.5,red -O -K -N >> "${output}.ps"

awk '{if ($1=="STAN") {f=$7/$9; x1=$8/$7; x2=$10/$9; df=f*sqrt(x1*x1+x2*x2); print f+df, $2}}' resmcnx.dat | awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' | gmt psxy -JX -R -W,gray -O -K -N >> "${output}.ps"
awk '{if ($1=="STAN") {f=$7/$9; x1=$8/$7; x2=$10/$9; df=f*sqrt(x1*x1+x2*x2); print f-df, $2}}' resmcnx.dat | awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' | gmt psxy -JX -R -W,gray -O -K -N >> "${output}.ps"

#awk '{if ($1=="STAN") print $7/$9, $2}' resmcnbf.dat | awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' | psxy -JX -R -W5/blue -O -K -N >> "${output}.ps"
#awk '{if ($1=="STAN") print $11/$12, $2}' resmcnx.dat | awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' | psxy -JX -R -W3/green -O -K -N >> "${output}.ps"

echo 0 0 | gmt psxy -JX -R -B0 -Sc0.001 -O >> "${output}.ps"

rm configa.dat
gmt psconvert -Tg "${output}.ps" -A
ls "$PWD/${output}.ps"
ls "$PWD/${output}.png"

[[ "$(uname)" == "Darwin" ]] && open "${output}.png"
