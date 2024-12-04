#!/bin/bash

# V2 converted from csh to bash by chatcpt and J. Pesicek, Oct 29, 2024
# ported to gmt6 Oct 30, 2024

cfg="config_eqx.dat"
if [ ! -f "$cfg" ]
then
ls $cfg
exit
fi

z0=$(awk 'NR==7 {print $1}' "$cfg")
eq=$(awk '{if (NR==30) print $1}' "$cfg") # needs to be past burn-in phase!
vv=$(awk '{if (NR==30) print $2}' "$cfg")
tot=$(echo "$vv $eq" | awk '{print $1+$2}')

str="<start model number post burn-in, b/t $eq and $tot> [<event number[s] to plot (default#: 200)>] [<window half length (5 km)>]"
me=`basename "$0"`

if [ "$#" -lt 1 ]
then
echo "${me}: $str"
echo "optionally plot multiple event numbers with quotes for 2nd input"
exit
else
echo post burn-in input number: $1
echo event ID[s]: $2
fi

if (( $1 > $eq && $1 < $tot )); then
    bi=$1
else
    echo "error: bad 1st input value"
    echo "$str"
    exit
fi

[ "$#" -gt 1 ] && eqn0=$2
[ "$#" -gt 1 ] || eqn0=200
echo "using quake number $eqn0"

f="resmcnx.dat"
ls "$f" "$cfg"
[ -f "$f" ] || exit

rm gmt.*
gmt set MEASURE_UNIT INCH
gmt set HEADER_FONT_SIZE 10
gmt set FONT_ANNOT_PRIMARY 10
gmt set HEADER_OFFSET 0.5c
gmt set LABEL_FONT_SIZE 10
gmt set COLOR_NAN  200/200/200

export LC_NUMERIC=C.UTF-8

if [ "$#" -lt 3 ]
then
window=5 # radius around mean
else
window=$3
fi
echo "window length is 2*${window} km"

for eqn in $eqn0; do

    echo "plotting quake: $eqn"

    output="x.ps"

    xy=$(awk '{if (($2=="'"$eqn"'") && ($1=="EZ")) printf"%f/%f/%f/%f\n", $3-"'"$window"'", $3+"'"$window"'", $4-"'"$window"'", $4+"'"$window"'"}' "$f")
    xz=$(awk '{if (($2=="'"$eqn"'") && ($1=="EZ")) printf"%f/%f/%f/%f\n", $3-"'"$window"'", $3+"'"$window"'", $5-"'"$window"'", $5+"'"$window"'"}' "$f")
    x0=$(awk '{if (($2=="'"$eqn"'") && ($1=="EZ")) printf"%f\n", $3-"'"$window"'"}' "$f")
    x1=$(awk '{if (($2=="'"$eqn"'") && ($1=="EZ")) printf"%f\n", $3+"'"$window"'"}' "$f")

    dx=0.25
    dy=0.25
    dz=0.25

    awk '{if (($3>(1*"'$bi'")) && ($4=="'"$eqn"'")) print $0}' tmpx | grep EQ > t77

    # x-y
    gmt psbasemap -JX5 -R$xy -B2f1:"X [km]":/2f1:"Y [km]":swEN -K -P -Y5.75 > "$output"
    awk '{print $6, $7}' t77 | \
    awk '{print int($1/"'"$dx"'"), int($2/"'"$dy"'")}' | \
    sort -n | \
    awk '{if (($1!=xold) || ($2!=yold)) {print xold*"'"$dx"'", yold*"'"$dy"'", s; s=0; xold=$1; yold=$2} else {s=s+1}}' | \
    tail -n +2 | gmt xyz2grd -Gtmpxy.grd -R -I"$dx/$dy" -V -F

    gmt grd2cpt -Chot -Z -D tmpxy.grd > tmp.cpt

    gmt grdimage tmpxy.grd -R -B0 -JX -Ctmp.cpt -K -O >> "$output"
    awk '{if (($2=="'"$eqn"'") && ($1=="EQ")) print $3, $4}' "$f" | \
    gmt psxy -JX -R -Sc0.075 -Glightblue -K -O -m >> "$output"
    awk '{if (($2=="'"$eqn"'") && ($1=="EQ")) print $3, $4, $6, $7}' "$f" | \
    gmt psxy -JX -R -Sc0.03 -Gblue -Exy0.01 -K -O -m >> "$output"

    # x-z
    gmt psbasemap -JX5/-5 -R$xz -B2f1:"X [km] EQ $eqn":/2f1:"Z [km]":SwEn -K -Y-5.0 -O >> "$output"

    awk '{print $6, $8-"'"$z0"'"}' t77 | \
    awk '{print int($1/"'"$dx"'"), int($2/"'"$dy"'")}' | \
    sort -n | \
    awk '{if (($1!=xold) || ($2!=yold)) {print xold*"'"$dx"'", yold*"'"$dy"'"+"'"$z0"'", s; s=0; xold=$1; yold=$2} else {s=s+1}}' | \
    tail -n +2 | gmt xyz2grd -Gtmpxy.grd -R -I"$dx/$dy" -V -F

    gmt grd2cpt -Chot -Z -D tmpxy.grd > tmp.cpt
    gmt grdimage tmpxy.grd -R -B0 -JX -Ctmp.cpt -K -O >> "$output"
    awk '{if (($2=="'"$eqn"'") && ($1=="EQ")) print $3, $5}' "$f" | gmt psxy -JX -R -Sc0.1 -Gwhite -K -O -m >> "$output"
    awk '{if (($2=="'"$eqn"'") && ($1=="EQ")) print $3, $5, $6, $8}' "$f" | gmt psxy -JX -R -Sc0.05 -Gblue -Exy0.01 -K -O -m >> "$output"
    awk '{if (($2=="'"$eqn"'") && ($1=="EZ")) print $3, $5}' "$f" | gmt psxy -JX -R -Sc0.1 -Gwhite -K -O -m >> "$output"
    awk '{if (($2=="'"$eqn"'") && ($1=="EZ")) print $3, $5, $6, $8}' "$f" | gmt psxy -JX -R -Sc0.05 -Ggreen -Exy0.01 -K -O -m >> "$output"
    
    awk '{if (($2=="'"$eqn"'") && ($1=="EM")) print $3, $5}' "$f" | gmt psxy -JX -R -Sc0.1 -Gwhite -K -O -m >> "$output"
    awk '{if (($2=="'"$eqn"'") && ($1=="EM")) print $3, $5}' "$f" | gmt psxy -JX -R -Sc0.05 -Gred -K -O -m >> "$output"
    
    echo "$x0" "$x1" "$z0" | awk '{print $1, $3; print $2, $3; print ">" }' | gmt psxy -JX -R -W -K -O -m >> "$output"
    
    echo 0 0 | gmt psxy -JX -R -B0 -Sc0.001 -O >> "$output"
    mv "$output" "loc_eq${eqn}.ps"
    gmt psconvert -Tg "loc_eq${eqn}.ps" -A
    ls "$PWD/loc_eq${eqn}.p"*
    [[ "$(uname)" == "Darwin" ]] && open "loc_eq${eqn}.png"
    
done

rm t77
