#!/bin/bash

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

z0=$(awk 'NR==7 {print $1}' "$cfg")
eq=$(awk '{if (NR==30) print $1}' "$cfg") # needs to be past burn-in phase!
vv=$(awk '{if (NR==30) print $2}' "$cfg")
tot=$(echo "$vv $eq" | awk '{print $1+$2}')
d=$(awk 'NR==1 {print $1}' "$cfg")

str="<start model number post burn-in, b/t $eq and $tot> [<event number to plot (default#: 100)>]"
me=`basename "$0"`

if [ "$#" -lt 1 ]
then
echo "${me}: $str"
echo "plots one quake"
exit
else
echo post burn-in input number: $1
echo event ID: $2
fi

if (( $1 > $eq && $1 < $tot )); then
    bi=$1
else
    echo "error: bad 1st input value"
    echo "$str"
    exit
fi

[ "$#" -gt 1 ] && eqn0=$2
[ "$#" -gt 1 ] || eqn0=100
echo "using quake number $eqn0"

f="resmcnx.dat"
ls "$f" "$cfg"
[ -f "$f" ] || exit
[[ -f  tmpx ]] || { echo "need to run disp_m_average_sl.sh first" ; exit 1 ; }

eqn=$eqn0
echo "plotting quake: $eqn"

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

awk '{if (($3>(1*"'$bi'")) && ($4=="'"$eqn"'")) print $0}' tmpx | grep EQ > t77

dx=`echo $d | awk '{print $1+$1*0.5}'`
dx=0.5
dy=$dx
dz=$dx
echo "grid increment: $dx km"

# x-y
gmt psbasemap -JX"$xs/$ys" -R"$xmin/$xmax/$ymin/$ymax" -Bxaf+l"X [km]" -Byaf+l"Y [km]" -BNWse -K -X0.7 -Y2.9 > eq$eqn.ps

awk '{print $6, $7}' t77 | \
awk '{print int($1/"'"$dx"'"), int($2/"'"$dy"'")}' | \
sort -n | \
awk '{if (($1!=xold) || ($2!=yold)) {print xold*"'"$dx"'", yold*"'"$dy"'", s; s=0; xold=$1; yold=$2} else {s=s+1}}' | \
tail -n +2 | gmt xyz2grd -Gtmpxy.grd -R -I"$dx/$dy" -V -F

gmt grd2cpt -Chot -Z -D tmpxy.grd > tmp.cpt

gmt grdimage tmpxy.grd -R -B0 -JX -Ctmp.cpt -K -O >> "eq$eqn.ps"
awk '{if (($2=="'"$eqn"'") && ($1=="EZ")) print $3, $4}' "$f" | gmt psxy -JX -R -Sc0.075 -Glightblue -K -O -m >> "eq$eqn.ps"
awk '{if (($2=="'"$eqn"'") && ($1=="EZ")) print $3, $4, $6, $7}' "$f" | gmt psxy -JX -R -Sc0.03 -Gblue -Exy0.01 -K -O -m >> "eq$eqn.ps"

# z-y
gmt psbasemap -JX"$zs/$ys" -R"$zmin/$zmax/$ymin/$ymax" -BseNw -Bxafg1000+l"Z [km]" -Byaf -K -O -X"$xs" >> eq$eqn.ps

awk '{print $8, $7}' t77 | awk '{print int($1/"'"$dz"'"), int($2/"'"$dy"'")}' | sort -n | \
awk '{if (($1!=xold) || ($2!=yold)) {print xold*"'"$dz"'", yold*"'"$dy"'", s; s=0; xold=$1; yold=$2} else {s=s+1}}' | \
tail -n +2 | gmt xyz2grd -Gtmpxy.grd -R -I"$dz/$dy" -V -F

gmt grd2cpt -Chot -Z -D tmpxy.grd > tmp.cpt
gmt grdimage tmpxy.grd -R -B0 -JX -Ctmp.cpt -K -O >> "eq$eqn.ps"
awk '{if (($2=="'"$eqn"'") && ($1=="EZ")) print $5, $4}' "$f" | gmt psxy -JX -R -Sc0.1 -Glightblue -K -O -m >> "eq$eqn.ps"
awk '{if (($2=="'"$eqn"'") && ($1=="EZ")) print $5, $4, $8, $7}' "$f" | gmt psxy -JX -R -Sc0.05 -Gblue -Exy0.01 -K -O -m >> "eq$eqn.ps"

# x-z
gmt psbasemap -JX"$xs/-$zs" -R"$xmin/$xmax/$zmin/$zmax" -BWsne -Byafg1000+l"Z [km]" -Bxaf -K -O -X-"$xs" -Y-"$zs" >> eq$eqn.ps

awk '{print $6, $8}' t77 | awk '{print int($1/"'"$dx"'"), int($2/"'"$dz"'")}' | sort -n | \
awk '{if (($1!=xold) || ($2!=yold)) {print xold*"'"$dx"'", yold*"'"$dz"'", s; s=0; xold=$1; yold=$2} else {s=s+1}}' | \
tail -n +2 | gmt xyz2grd -Gtmpxy.grd -R -I"$dx/$dz" -V -F

gmt grd2cpt -Chot -Z -D tmpxy.grd > tmp.cpt
gmt grdimage tmpxy.grd -R -B0 -JX -Ctmp.cpt -K -O >> "eq$eqn.ps"
awk '{if (($2=="'"$eqn"'") && ($1=="EZ")) print $3, $5}' "$f" | gmt psxy -JX -R -Sc0.1 -Glightblue -K -O -m >> "eq$eqn.ps"
awk '{if (($2=="'"$eqn"'") && ($1=="EZ")) print $3, $5, $6, $8}' "$f" | gmt psxy -JX -R -Sc0.05 -Gblue -Exy0.01 -K -O -m >> "eq$eqn.ps"

echo 0 0 | gmt psxy -JX -R -Sc0.001 -O -U"`basename $PWD`" -Y1c -X1c >> eq$eqn.ps

gmt psconvert -Tg eq$eqn.ps -A
ls $PWD/eq$eqn.ps 
ls $PWD/eq$eqn.png
[[ "$(uname)" == "Darwin" ]] && open eq$eqn.png

