#!/bin/bash

# V2 converted from csh to bash by chatcpt and J. Pesicek, Oct 29, 2024
# ported to gmt6 Oct 30, 2024

echo "GMT version:"
gmt --version

rm -f gmt.*
gmt set MEASURE_UNIT INCH
gmt set FONT_ANNOT_PRIMARY 12
gmt set HEADER_FONT_SIZE 12
gmt set HEADER_OFFSET 0.5c
gmt set LABEL_FONT_SIZE 16
gmt set COLOR_NAN 200

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
tmax=$(awk '{print $5}' tmppbi | gmt pshistogram -R0/10/0/110 -W0.00001 -Z1 -Q -G0 -K -O -V -Io | awk '$2 > 50 {printf "%f\n", $1}' | head -1)
awk -v tmax="$tmax" '$5 < (1.0 * tmax) {print $0}' tmppbi > tmppbibf

# Misfit
awk '{print $5}' tmp | gmt pshistogram -R0/20/0/1 -F -W0.001 -IO -Z0 > mall.dat
awk '{print $5}' tmpeq | gmt pshistogram -R0/20/0/1 -F -W0.001 -IO -Z0 > malleq.dat
awk '{print $5}' tmppbi | gmt pshistogram -R0/20/0/1 -F -W0.001 -IO -Z0 > mallpbi.dat
awk '{print $5}' tmppbibf | gmt pshistogram -R0/20/0/1 -F -W0.001 -IO -Z0 > mallpbibf.dat
awk '{print $5}' tmppbi | gmt pshistogram -R0/120/0/110 -W0.001 -Z1 -Q -G0 -K -O -V -Io > mallpbic.dat

# Dim
awk '{print $4}' tmp | gmt pshistogram -R0/500/0/1 -F -W1 -IO -Z0 > dall.dat
awk '{print $4}' tmpeq | gmt pshistogram -R0/500/0/1 -F -W1 -IO -Z0 > dalleq.dat
awk '{print $4}' tmppbi | gmt pshistogram -R0/500/0/1 -F -W1 -IO -Z0 > dallpbi.dat
awk '{print $4}' tmppbibf | gmt pshistogram -R0/500/0/1 -F -W1 -IO -Z0 > dallpbibf.dat

# Misfit
xmin=0
xmax=$(echo "$eq $v2" | awk '{print ($1+$2)*1.1}')

ymin=0.01
ymax=20.0
lymin=$(echo "$ymin" | awk '{print log($1)}')
lymax=$(echo "$ymax" | awk '{print log($1)}')

gmt psbasemap -JX5/3 -R$xmin/$xmax/$lymin/$lymax -B0n -K -P -Y2.0 -U"`basename $PWD`" > evo.ps

dx=0.03

awk '{if ($5!="nan") print $3,log($5),1}' tmp | gmt blockmean -R$xmin/$xmax/$lymin/$lymax -I$deci/$dx -Sw | awk '{print $1, $2, $3-1}' | gmt xyz2grd -Gtmp.grd -R$xmin/$xmax/$lymin/$lymax -I$deci/$dx -V

zmax=$(gmt grdinfo tmp.grd | grep v_max | awk '{print $5}')
gmt makecpt -Chot -Z -T0/$zmax/1 -D > tmp.cpt
gmt grdimage tmp.grd -R -B0 -JX -Ctmp.cpt -K -O >> evo.ps
echo "$eq $lymin $eq $lymax" | awk '{print $1, $2; print $3, $4; print "#"}' | gmt psxy -JX -R -B -W -K -O >> evo.ps
echo "$bi $lymin $bi $lymax" | awk '{print $1, $2; print $3, $4; print "#"}' | gmt psxy -JX -R -B -W -K -O >> evo.ps
gmt psbasemap -JX5/3l -R$xmin/$xmax/$ymin/$ymax -Baf:"evaluated models":/2f3:"Misfit [s]":nSW -K -O >> evo.ps

gmt psbasemap -JX-1/3l -R0/1/$ymin/$ymax -B0/2f3nSEw -K -O -X5 >> evo.ps

max=$(sort -n -k 2 malleq.dat | tail -1 | awk '{print $2}')
awk -v max="$max" '{if ($1 > 0) print $2 / max, $1}' malleq.dat | gmt psxy -JX -R -K -O -Gblack -L >> evo.ps
awk -v max="$max" '{if ($1 > 0) print $2 / max, $1}' mallpbi.dat | gmt psxy -JX -R -K -O -Gblue -L >> evo.ps
awk -v max="$max" '{if ($1 > 0) print $2 / max, $1}' mallpbibf.dat | gmt psxy -JX -R -K -O -Gred -L >> evo.ps
awk '{if ($1 > 0) print $2 / 100.0, $1}' mallpbic.dat | gmt psxy -JX -R -K -O -W1,green >> evo.ps

# Dim
ymin=0
ymax=30
gmt psbasemap -JX5/2 -R$xmin/$xmax/$ymin/$ymax -B0n -K -P -X-5 -Y4 -O >> evo.ps

awk '{print $3,$4,1}' tmp | gmt blockmean -R$xmin/$xmax/$ymin/$ymax -I$deci/1 -Sw | awk '{print $1, $2, $3-1}' | gmt xyz2grd -Gtmp.grd -R$xmin/$xmax/$ymin/$ymax -I$deci/1 -V

zmax=$(gmt grdinfo tmp.grd | egrep v_max | awk '{print $5}')
gmt makecpt -Chot -Z -T0/$zmax/1 -D > tmp.cpt
gmt grdimage tmp.grd -R -B0 -JX -Ctmp.cpt -K -O >> evo.ps
echo "$eq $ymin $eq $ymax" | awk '{print $1, $2; print $3, $4; print "#"}' | gmt psxy -JX -R -B0 -W -K -O  >> evo.ps
echo "$bi $ymin $bi $ymax" | awk '{print $1, $2; print $3, $4; print "#"}' | gmt psxy -JX -R -B0 -W -K -O  >> evo.ps

gmt psbasemap -JX5/2 -R$xmin/$xmax/$ymin/$ymax -B50000f10000:"evaluated models":/10f5:"Dim":nSW -K -O >> evo.ps

gmt psbasemap -JX-1/2 -R0/1/$ymin/$ymax -B0/10f5nSEw -K -O -X5 >> evo.ps

max=$(sort -n -k 2 dalleq.dat | tail -1 | awk '{print $2}')
awk -v max="$max" '{print $2/max, $1}' dalleq.dat | gmt psxy -JX -R -K -O -Gblack -L >> evo.ps
awk -v max="$max" '{print $2/max, $1}' dallpbi.dat | gmt psxy -JX -R -K -O -Gblue -L >> evo.ps
awk -v max="$max" '{print $2/max, $1}' dallpbibf.dat | gmt psxy -JX -R -K -O -Gred -L >> evo.ps

echo 0 0 | gmt psxy -JX -R -B0 -Sc0.001 -O >> evo.ps
rm dall*dat mall*dat tmp*

gmt psconvert -Tg evo.ps -A
ls "$PWD/evo.ps"
ls "$PWD/evo.png"

[[ "$(uname)" == "Darwin" ]] && open evo.png

