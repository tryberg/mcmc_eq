#!/bin/bash

# V2 converted from csh to bash by chatcpt and J. Pesicek, Oct 29, 2024

rm .gmt*
gmtset MEASURE_UNIT INCH
gmtset ANOT_FONT_SIZE 10
gmtset HEADER_FONT_SIZE 10
gmtset HEADER_OFFSET 0.5c
gmtset LABEL_FONT_SIZE 10
gmtset COLOR_NAN  200/200/200

export LC_NUMERIC=C.UTF-8

f="resmcnx.dat"
ls "$f"
[ -f "$f" ] || exit

window=5

cfg="config_eqx.dat"
ls "$cfg"
[ -f "$cfg" ] || exit
z0=$(awk 'NR==7 {print $1}' "$cfg")

burnin=200000
echo "using $burnin for burnin"

for eqn in 200; do

    echo "plotting quake: $eqn"

    output="x.ps"

    xy=$(awk '{if (($2=="'"$eqn"'") && ($1=="EZ")) printf"%f/%f/%f/%f\n", $3-"'"$window"'", $3+"'"$window"'", $4-"'"$window"'", $4+"'"$window"'"}' "$f")
    xz=$(awk '{if (($2=="'"$eqn"'") && ($1=="EZ")) printf"%f/%f/%f/%f\n", $3-"'"$window"'", $3+"'"$window"'", $5-"'"$window"'", $5+"'"$window"'"}' "$f")
    x0=$(awk '{if (($2=="'"$eqn"'") && ($1=="EZ")) printf"%f\n", $3-"'"$window"'"}' "$f")
    x1=$(awk '{if (($2=="'"$eqn"'") && ($1=="EZ")) printf"%f\n", $3+"'"$window"'"}' "$f")

    dx=0.25
    dy=0.25
    dz=0.25

    awk '{if (($3>(1*"'$burnin'")) && ($4=="'"$eqn"'")) print $0}' tmpx | grep EQ > t77

    # x-y
    psbasemap -JX5 -R$xy -B2f1:"X [km]":/2f1:"Y [km]":swEN -K -P -Y5.75 > "$output"
    awk '{print $6, $7}' t77 | \
    awk '{print int($1/"'"$dx"'"), int($2/"'"$dy"'")}' | \
    sort -n | \
    awk '{if (($1!=xold) || ($2!=yold)) {print xold*"'"$dx"'", yold*"'"$dy"'", s; s=0; xold=$1; yold=$2} else {s=s+1}}' | \
    tail -n +2 | xyz2grd -Gtmpxy.grd -R -I"$dx/$dy" -V -F

    grd2cpt -Chot -Z -D tmpxy.grd > tmp.cpt

    grdimage tmpxy.grd -R -B0 -JX -Ctmp.cpt -K -O >> "$output"
    awk '{if (($2=="'"$eqn"'") && ($1=="EQ")) print $3, $4}' "$f" | \
    psxy -JX -R -Sc0.075 -Glightblue -K -O -M >> "$output"
    awk '{if (($2=="'"$eqn"'") && ($1=="EQ")) print $3, $4, $6, $7}' "$f" | \
    psxy -JX -R -Sc0.03 -Gblue -Exy0.01 -K -O -M >> "$output"

    # x-z
    psbasemap -JX5/-5 -R$xz -B2f1:"X [km] EQ $eqn":/2f1:"Z [km]":SwEn -K -Y-5.0 -O >> "$output"

    awk '{print $6, $8-"'"$z0"'"}' t77 | \
    awk '{print int($1/"'"$dx"'"), int($2/"'"$dy"'")}' | \
    sort -n | \
    awk '{if (($1!=xold) || ($2!=yold)) {print xold*"'"$dx"'", yold*"'"$dy"'"+"'"$z0"'", s; s=0; xold=$1; yold=$2} else {s=s+1}}' | \
    tail -n +2 | xyz2grd -Gtmpxy.grd -R -I"$dx/$dy" -V -F

    grd2cpt -Chot -Z -D tmpxy.grd > tmp.cpt
    grdimage tmpxy.grd -R -B0 -JX -Ctmp.cpt -K -O >> "$output"
    awk '{if (($2=="'"$eqn"'") && ($1=="EQ")) print $3, $5}' "$f" | psxy -JX -R -Sc0.1 -Gwhite -K -O -M >> "$output"
    awk '{if (($2=="'"$eqn"'") && ($1=="EQ")) print $3, $5, $6, $8}' "$f" | psxy -JX -R -Sc0.05 -Gblue -Exy0.01 -K -O -M >> "$output"
    awk '{if (($2=="'"$eqn"'") && ($1=="EZ")) print $3, $5}' "$f" | psxy -JX -R -Sc0.1 -Gwhite -K -O -M >> "$output"
    awk '{if (($2=="'"$eqn"'") && ($1=="EZ")) print $3, $5, $6, $8}' "$f" | psxy -JX -R -Sc0.05 -Ggreen -Exy0.01 -K -O -M >> "$output"
    
    awk '{if (($2=="'"$eqn"'") && ($1=="EM")) print $3, $5}' "$f" | psxy -JX -R -Sc0.1 -Gwhite -K -O -M >> "$output"
    awk '{if (($2=="'"$eqn"'") && ($1=="EM")) print $3, $5}' "$f" | psxy -JX -R -Sc0.05 -Gred -K -O -M >> "$output"
    
    echo "$x0" "$x1" "$z0" | awk '{print $1, $3; print $2, $3; print ">" }' | psxy -JX -R -W1 -K -O -M >> "$output"
    
    echo 0 0 | psxy -JX -R -B0 -Sc0.001 -O >> "$output"
    mv "$output" "eq${eqn}.ps"
    ps2raster -Tg "eq${eqn}.ps" -A
    ls "$PWD/eq${eqn}.p"*
    if [[ "$(uname)" == "Darwin" ]]; then open "eq${eqn}.png"; fi
    
done

rm t77
