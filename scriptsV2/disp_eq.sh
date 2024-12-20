#!/bin/bash

# V2 converted from csh to bash by chatcpt and J. Pesicek, Oct 29, 2024
# ported to gmt6 Oct 30, 2024

rm -f gmt.*
gmt set MEASURE_UNIT INCH
gmt set FONT_ANNOT_PRIMARY 10
gmt set HEADER_FONT_SIZE 10
gmt set HEADER_OFFSET 0.5c
gmt set LABEL_FONT_SIZE 10
gmt set COLOR_NAN  200/200/200

cfg=config_eqx.dat
ls "$cfg"
[[ -f "$cfg" ]] || exit

res=resmcnx.dat
ls "$res"
[[ -f "$res" ]] || exit

echo "reading $# args..."
[[ $# -eq 1 ]] || { echo "need mcmcPhaseFile (picks) as input arg"; exit; }

picks=$1
ls "$picks"
[[ -f "$picks" ]] || exit

quakes=quakes.dat # if prior locations file exists
ls $quakes

d=$(awk 'NR==1 {print $1}' "$cfg")
nz=$(awk 'NR==4 {print $1}' "$cfg")
zmin=$(awk 'NR==7 {print $1}' "$cfg")
zmax=$(echo "$zmin $nz $d" | awk '{print $1+$2*$3}')

nx=$(awk 'NR==2 {print $1}' "$cfg")
xmin=$(awk 'NR==5 {print $1}' "$cfg")
xmax=$(echo "$xmin $nx $d" | awk '{print $1+$2*$3}')

ny=$(awk 'NR==3 {print $1}' "$cfg")
ymin=$(awk 'NR==6 {print $1}' "$cfg")
ymax=$(echo "$ymin $ny $d" | awk '{print $1+$2*$3}')

# set map scaling (JP)
xyr=$(echo "$xmin $xmax $ymin $ymax" | awk '{print ($2 - $1)/($4 - $3)}')
ys=4.5
xs=$(echo "$ys $xyr" | awk '{print $1*$2}')
zs=$(echo "$ys 2" | awk '{print $1/$2}')

awk '{if ($1=="RES") print $2, $3, $4}' "$res" > t1

awk '{if ($1!="#" && $2!="NA") print $1, $2, $4, $5, $6}' "$picks" | sort | uniq > rec.dat

egrep EZ "$res" > resmcna.tmp

# x-y
gmt psbasemap -JX"$xs/$ys" -R"$xmin/$xmax/$ymin/$ymax" -Baf -BNWse -K -X0.7 -Y2.9 > eq.ps

paste t1 rec.dat > recdata
stac=staCors_mcmc.dat
awk '{print $4,$2,$3}' recdata > "$stac"
echo "station corrections saved:"
ls "$PWD/$stac"

awk '{print $6, $7}' recdata | gmt psxy -JX -R -St0.15 -Glightblue -K -O >> eq.ps
awk '{if ($2>0) print $6, $7, $2*0.8}' recdata | gmt psxy -JX -R -Sc -K -O >> eq.ps
awk '{if ($2<0) print $6, $7, -$2*0.8}' recdata | gmt psxy -JX -R -Sx  -K -O >> eq.ps

awk '{print $3, $4}' resmcna.tmp | gmt psxy -JX -R -Sc0.04 -Gred -K -O >> eq.ps
awk '{print $3, $4, $6/2.0, $7/2.0}' resmcna.tmp | gmt psxy -JX -R -Sc0.001 -Exy0.01 -K -O >> eq.ps

if [[ -f $quakes ]]; then
    awk '{print $2, $3}' $quakes | gmt psxy -JX -R -Sc0.01 -Gblue -K -O >> eq.ps
fi

# z-y
gmt psbasemap -JX"$zs/$ys" -R"$zmin/$zmax/$ymin/$ymax" -BsENw -Bxafg1000 -Byaf+l"$PWD" -K -O -X"$xs" >> eq.ps
awk '{print $8, $7}' recdata | gmt psxy -JX -R -St0.15 -Glightblue -K -O -N >> eq.ps
awk '{print $5, $4}' resmcna.tmp | gmt psxy -JX -R -Sc0.04 -Gred -K -O >> eq.ps
awk '{print $5, $4, $8/2.0, $7/2.0}' resmcna.tmp | gmt psxy -JX -R -Sc0.001 -Exy0.01 -K -O -N >> eq.ps

awk '{if ($5<'$zmin') print $0}' resmcna.tmp > zbad.tmp
na=$(cat zbad.tmp | wc -l)
if [ "$na" -gt 0 ]; then
    echo "WARNING: $na events above model plotted in green:"
    cat zbad.tmp
fi
awk '{if ($5<'$zmin') print $5, $4}' resmcna.tmp | gmt psxy -JX -R -Sc0.04 -Ggreen -K -O -N >> eq.ps

if [ -f $quakes ]; then  # known locations file for comparison (ID,X,Y,Z,OT,0)
    awk '{print $4, $3}' $quakes | gmt psxy -JX -R -Sc0.01 -Gblue -K -O >> eq.ps
fi

# x-z
gmt psbasemap -JX"$xs/-$zs" -R"$xmin/$xmax/$zmin/$zmax" -BWsne -Byafg1000 -Bxaf -K -O -X-"$xs" -Y-"$zs" >> eq.ps
awk '{print $6, $8}' recdata | gmt psxy -JX -R -St0.15 -Glightblue -K -O -N >> eq.ps
awk '{print $3, $5}' resmcna.tmp | gmt psxy -JX -R -Sc0.04 -Gred -K -O >> eq.ps
awk '{print $3, $5, $6/2.0, $8/2.0}' resmcna.tmp | gmt psxy -JX -R -Sc0.001 -Exy0.01 -K -O -N >> eq.ps

awk '{if ($5<'$zmin') print $3, $5}' resmcna.tmp | gmt psxy -JX -R -Sc0.04 -Ggreen -K -O -N >> eq.ps

if [ -f $quakes ]; then  # known locations file for comparison (ID,X,Y,Z,OT,0)
    awk '{print $2, $4}' $quakes | gmt psxy -JX -R -Sc0.01 -Gblue -K -O >> eq.ps
fi

echo 0 0 | gmt psxy -JX -R -B0 -Sc0.001 -O >> eq.ps

gmt psconvert -Tg eq.ps -A
ls $PWD/eq.ps 
ls $PWD/eq.png
[[ "$(uname)" == "Darwin" ]] && open eq.png

# output new catalog: {'X','Y','Z','OT','dOT','ex','ey','ez'}
cat="quakes_mcmc.dat"
echo "# {'X','Y','Z','OT','dOT','ex','ey','ez'} (%8.3f %8.3f %8.3f %015.3f %7.3f %f %f %f\n)" > "$cat"
awk '{printf "%8.3f %8.3f %8.3f %015.3f %7.3f %f %f %f\n",$3,$4,$5,$9,$10,$6,$7,$8}' resmcna.tmp >> "$cat"
echo "new output catalog saved:"
ls "$PWD/$cat"

rm t1 rec.dat resmcna.tmp
