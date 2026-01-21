#!/bin/bash

cfg=config_eqx.dat
ls "$cfg"
[[ -f "$cfg" ]] || exit

res=resmcnx.dat
ls "$res"
[[ -f "$res" ]] || exit
rm -f gmt.*

z0=$(awk 'NR==7 {print $1}' "$cfg")
eq=$(awk '{if (NR==30) print $1}' "$cfg") # needs to be past burn-in phase!
vv=$(awk '{if (NR==30) print $2}' "$cfg")
tot=$(echo "$vv $eq" | awk '{print $1+$2}')

str="<start model number post burn-in, b/t $eq and $tot> [<event number to plot (100)>] [<number of chains to show (20)>]" 
me=`basename "$0"`

if [ "$#" -lt 1 ]
then
echo "${me}: $str"
exit
else
echo post burn-in input number: $1
fi

if (( $1 > $eq && $1 < $tot )); then
    bi=$1
else
    echo "error: bad 1st input value"
    echo "$str"
    exit
fi

gmt set MEASURE_UNIT INCH
gmt set FONT_ANNOT_PRIMARY 10
gmt set HEADER_FONT_SIZE 10
gmt set HEADER_OFFSET 0.5c
gmt set LABEL_FONT_SIZE 10
gmt set COLOR_NAN 200/200/200

[ "$#" -gt 1 ] && eq=$2
[ "$#" -gt 1 ] || eq=100
[ "$#" -gt 2 ] && nc=$3
[ "$#" -gt 2 ] || nc=20 

bi=$1
echo "using quake number: $eq"
output="loc_evo${eq}"
echo "using number of chains: $nc"

p="."

d=$(awk '{if (NR==1) print $1}' "$p/$cfg")
nx=$(awk '{if (NR==2) print $1}' "$p/$cfg")
ny=$(awk '{if (NR==3) print $1}' "$p/$cfg")
nz=$(awk '{if (NR==4) print $1}' "$p/$cfg")
x0=$(awk '{if (NR==5) print $1}' "$p/$cfg")
y0=$(awk '{if (NR==6) print $1}' "$p/$cfg")
z0=$(awk '{if (NR==7) print $1}' "$p/$cfg")
x1=$(echo "$x0 $nx $d" | awk '{print $1+($2-1)*$3}')
y1=$(echo "$y0 $ny $d" | awk '{print $1+($2-1)*$3}')
z1=$(echo "$z0 $nz $d" | awk '{print $1+($2-1)*$3}')

echo "$x0" "$x1" "$y0" "$y1" "$z0" "$z1"

gmt psbasemap -JX6 -R"$x0/$x1/$y0/$y1" -BNWes -Bxaf+l"X [km]" -Byaf+l"Y [km]" -K -P -Y4.5 > "${output}.ps"

echo XY
for n in $(seq 1 $nc); do
    nn=`echo $n | awk '{printf "%03d\n",$1}'`
    f="$p/rjx-$nn"*.out
    ls $f
    cat $f | awk '{if (($4=="'"$eq"'") && ($1=="EQ")) print $6, $7; if ($1=="sta") print ">"}' | gmt psxy -JX -R -W  -K -O -m">" >> "${output}.ps"
    cat $f | awk '{if (($4=="'"$eq"'") && ($1=="EQ")) print $6, $7}' | gmt psxy -JX -R -Sc0.01 -G0 -K -O -m">" >> "${output}.ps"
    cat $f | awk '{if (($4=="'"$eq"'") && ($1=="EQ")) print $6, $7}' | head -1 | gmt psxy -JX -R -Sc0.05 -Ggreen -K -O -m">" >> "${output}.ps"
done

cat "$p/$res" | awk '{if (($2=="'"$eq"'") && ($1=="EQ")) print $3, $4}' | gmt psxy -JX -R -Sc0.075 -Glightblue -K -O -m">" >> "${output}.ps"
cat "$p/$res" | awk '{if (($2=="'"$eq"'") && ($1=="EQ")) print $3, $4}' | gmt psxy -JX -R -Sc0.03 -Gblue -K -O -m">" >> "${output}.ps"

echo $x0 $y1 "Event: "$eq"" | gmt pstext -JX -R -K -O -N -F+jTL -D.1/-.1 >> "${output}.ps"

gmt psbasemap -JX6/-3 -R"$x0/$x1/$z0/$z1" -BSWen -Bxaf+l"X [km]" -Byaf+l"Z [km]" -K -P -Y-3 -O >> "${output}.ps"

echo XZ
for n in $(seq 1 $nc); do
    nn=`echo $n | awk '{printf "%03d\n",$1}'`
    f="$p/rjx-$nn"*.out
    ls $f
    cat $f | awk '{if (($4=="'"$eq"'") && ($1=="EQ")) print $6, $8; if ($1=="sta") print ">"}' | gmt psxy -JX -R -W -K -O -m">" >> "${output}.ps"
    cat $f | awk '{if (($4=="'"$eq"'") && ($1=="EQ")) print $6, $8}' | gmt psxy -JX -R -Sc0.01 -G0 -K -O -m">" >> "${output}.ps"
    cat $f | awk '{if (($4=="'"$eq"'") && ($1=="EQ")) print $6, $8}' | head -1 | gmt psxy -JX -R -Sc0.05 -Ggreen -K -O -m">" >> "${output}.ps"
done

cat "$p/$res" | awk '{if (($2=="'"$eq"'") && ($1=="EQ")) print $3, $5}' | gmt psxy -JX -R -Sc0.075 -Glightblue -K -O -m">" >> "${output}.ps"
cat "$p/$res" | awk '{if (($2=="'"$eq"'") && ($1=="EQ")) print $3, $5}' | gmt psxy -JX -R -Sc0.03 -Gblue -K -O -m">" >> "${output}.ps"

echo 0 0 | gmt psxy -JX -R -B0 -Sc0.001 -O >> "${output}.ps"

gmt psconvert -Tg "${output}.ps" -A
ls "$PWD/${output}.ps"

[[ "$(uname)" == "Darwin" ]] && open $output.png
