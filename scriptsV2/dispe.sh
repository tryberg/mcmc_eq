#!/bin/bash

# V2 converted from csh to bash by chatcpt and J. Pesicek, Oct 29, 2024

echo "Checking for GMT4:"
whereis -b pshistogram
if [[ -z "$(whereis -q -b pshistogram)" ]]; then
    echo "install GMT4?"
    exit
fi

rm -f .gmt*
gmtset MEASURE_UNIT INCH
gmtset ANOT_FONT_SIZE 12
gmtset HEADER_FONT_SIZE 12
gmtset HEADER_OFFSET 0.5c
gmtset LABEL_FONT_SIZE 16
gmtset COLOR_NAN 200

export LC_NUMERIC=C.UTF-8

cfg="config_eqx.dat"
ls "$cfg"
[[ ! -f "$cfg" ]] && exit

eq=$(awk 'NR==30 {print $1}' "$cfg")
v2=$(awk 'NR==30 {print $2}' "$cfg")
bi=$(echo "$v2 $eq" | awk '{print ($1+$2)-$2}')
deci=$(awk 'NR==31 {print $1}' "$cfg")

touch tmp; rm tmp; touch tmp
for n in rjx-*.out; do
    echo "$n"
    awk '$1=="mod" {print $0}' "$n" >> tmp
done

awk -v eq="$eq" '$3 > (1.0 * eq) {print $0}' tmp > tmpeq

awk -v bi="$bi" '$3 > (1.0 * bi) {print $0}' tmp > tmppbi
tmax=$(awk '{print $5}' tmppbi | pshistogram -R0/10/0/110 -W0.00001 -Z1 -Q -G0 -K -O -V -Io | awk '$2 > 50 {printf "%f\n", $1}' | head -1)
awk -v tmax="$tmax" '$5 < (1.0 * tmax) {print $0}' tmppbi > tmppbibf

# Misfit
awk '{print $5}' tmp | pshistogram -R0/20/0/1 -F -W0.001 -IO -Z0 > mall.dat
awk '{print $5}' tmpeq | pshistogram -R0/20/0/1 -F -W0.001 -IO -Z0 > malleq.dat
awk '{print $5}' tmppbi | pshistogram -R0/20/0/1 -F -W0.001 -IO -Z0 > mallpbi.dat
awk '{print $5}' tmppbibf | pshistogram -R0/20/0/1 -F -W0.001 -IO -Z0 > mallpbibf.dat
awk '{print $5}' tmppbi | pshistogram -R0/120/0/110 -W0.001 -Z1 -Q -G0 -K -O -V -Io > mallpbic.dat

# Dim
awk '{print $4}' tmp | pshistogram -R0/500/0/1 -F -W1 -IO -Z0 > dall.dat
awk '{print $4}' tmpeq | pshistogram -R0/500/0/1 -F -W1 -IO -Z0 > dalleq.dat
awk '{print $4}' tmppbi | pshistogram -R0/500/0/1 -F -W1 -IO -Z0 > dallpbi.dat
awk '{print $4}' tmppbibf | pshistogram -R0/500/0/1 -F -W1 -IO -Z0 > dallpbibf.dat

# Misfit
xmin=0
xmax=$(echo "$eq $v2" | awk '{print ($1+$2)*1.1}')

ymin=0.02
ymax=20.0
lymin=$(echo "$ymin" | awk '{print log($1)}')
lymax=$(echo "$ymax" | awk '{print log($1)}')

psbasemap -JX5/3 -R$xmin/$xmax/$lymin/$lymax -B0n -K -P -Y2.0 > evo.ps

dx=0.03

awk '{print $3,log($5)}' tmp | awk '{print int($1/"'$deci'"), int($2/"'$dx'")}' | sort -n | awk '{if (($1!=xold) || ($2!=yold)) {print xold*"'$deci'", yold*"'$dx'", s; s=0; xold=$1; yold=$2} else {s=s+1}}' | tail +2 | xyz2grd -Gtmp.grd -R$xmin/$xmax/$lymin/$lymax -I$deci/$dx -V

zmax=$(grdinfo tmp.grd | grep z_max | awk '{print $5}')
makecpt -Chot -Z -T0/$zmax/1 -D > tmp.cpt
grdimage tmp.grd -R -B0 -JX -Ctmp.cpt -K -O >> evo.ps
echo "$eq $lymin $eq $lymax" | awk '{print $1, $2; print $3, $4; print ">"}' | psxy -JX -R -B -W1 -K -O -M">" >> evo.ps
echo "$bi $lymin $bi $lymax" | awk '{print $1, $2; print $3, $4; print ">"}' | psxy -JX -R -B -W1 -K -O -M">" >> evo.ps
psbasemap -JX5/3l -R$xmin/$xmax/$ymin/$ymax -B50000f10000:"evaluated models $(pwd)":/2f3:"Misfit [s]":nSW -K -O >> evo.ps

psbasemap -JX-1/3l -R0/1/$ymin/$ymax -B0/2f3nSEw -K -O -X5 >> evo.ps

max=$(sort -n -k 2 malleq.dat | tail -1 | awk '{print $2}')
awk -v max="$max" '{if ($1 > 0) print $2 / max, $1}' malleq.dat | psxy -JX -R -K -O -Gblack -L >> evo.ps
awk -v max="$max" '{if ($1 > 0) print $2 / max, $1}' mallpbi.dat | psxy -JX -R -K -O -Gblue -L >> evo.ps
awk -v max="$max" '{if ($1 > 0) print $2 / max, $1}' mallpbibf.dat | psxy -JX -R -K -O -Gred -L >> evo.ps
awk '{if ($1 > 0) print $2 / 100.0, $1}' mallpbic.dat | psxy -JX -R -K -O -W1/green >> evo.ps

# Dim
ymin=0
ymax=30
psbasemap -JX5/2 -R$xmin/$xmax/$ymin/$ymax -B0n -K -P -X-5 -Y4 -O >> evo.ps

awk '{print $3,$4}' tmp | awk '{print int($1/"'$deci'"), $2}' | sort -n | awk '{if (($1!=xold) || ($2!=yold)) {print xold*"'$deci'", yold, s; s=0; xold=$1; yold=$2} else {s=s+1}}' | tail +2 | xyz2grd -Gtmp.grd -R$xmin/$xmax/$ymin/$ymax -I$deci/1 -V

zmax=$(grdinfo tmp.grd | egrep z_max | awk '{print $5}')
makecpt -Chot -Z -T0/$zmax/1 -D > tmp.cpt
grdimage tmp.grd -R -B0 -JX -Ctmp.cpt -K -O >> evo.ps
echo "$eq $ymin $eq $ymax" | awk '{print $1, $2; print $3, $4; print ">"}' | psxy -JX -R -B0 -W1 -K -O -M">" >> evo.ps
echo "$bi $ymin $bi $ymax" | awk '{print $1, $2; print $3, $4; print ">"}' | psxy -JX -R -B0 -W1 -K -O -M">" >> evo.ps

psbasemap -JX5/2 -R$xmin/$xmax/$ymin/$ymax -B50000f10000:"evaluated models":/10f5:"Dim":nSW -K -O >> evo.ps

psbasemap -JX-1/2 -R0/1/$ymin/$ymax -B0/10f5nSEw -K -O -X5 >> evo.ps

max=$(sort -n -k 2 dalleq.dat | tail -1 | awk '{print $2}')
awk -v max="$max" '{print $2/max, $1}' dalleq.dat | psxy -JX -R -K -O -Gblack -L >> evo.ps
awk -v max="$max" '{print $2/max, $1}' dallpbi.dat | psxy -JX -R -K -O -Gblue -L >> evo.ps
awk -v max="$max" '{print $2/max, $1}' dallpbibf.dat | psxy -JX -R -K -O -Gred -L >> evo.ps

echo 0 0 | psxy -JX -R -B0 -Sc0.001 -O >> evo.ps
rm dall*dat mall*dat tmp*

ps2raster -Tg evo.ps -Au
ls "$PWD/evo.ps"
ls "$PWD/evo.png"

[[ "$(uname)" == "Darwin" ]] && open evo.png

