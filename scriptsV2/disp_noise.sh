#!/bin/bash

# V2 converted from csh to bash by chatcpt and J. Pesicek, Oct 29, 2024
# ported to gmt6 Oct 30, 2024

rm gmt.*
gmt set MEASURE_UNIT INCH
gmt set FONT_ANNOT_PRIMARY 10
gmt set HEADER_FONT_SIZE 10
gmt set HEADER_OFFSET 0.5c
gmt set LABEL_FONT_SIZE 10
gmt set COLOR_NAN  200/200/200

# noise assumed to be in tmpx
ls tmpx
if [ ! -f tmpx ]; then exit; fi
egrep mod tmpx | awk '{if ($3>0) print $5, $6, $7, $8, $9, $10, $11, $12, $13}' > t1

bw=0.001
ymax=`gmt pshistogram -I t1 -T$bw -o3`
xmax=`awk '{for(i=1;i<=NF;i++){v=($i<0)?-$i:$i; if(v>m)m=v}} END{print m}' t1 | gmt info -C -I0.25 -o1`
xmax=$(awk -v a="$xmax" -v b="1" 'BEGIN{print (a<b)?a:b}') # use reasonable max if less than 1, otherwise 1

gmt psbasemap -JX8/6 -R0/$xmax/0/$ymax -Bxafg0.05:"Seconds" -Byaf -BnSWe -K > noise.ps

awk '{print $1}' t1 | gmt pshistogram -JX -R -T$bw  -Z0 -G0 -K -O -V -F >> noise.ps

awk '{print $2}' t1 | gmt pshistogram -JX -R -T$bw  -Z0 -G0/0/255 -K -O -V -F >> noise.ps
awk '{print $3}' t1 | gmt pshistogram -JX -R -T$bw  -Z0 -G80/80/255 -K -O -V -F >> noise.ps
awk '{print $4}' t1 | gmt pshistogram -JX -R -T$bw  -Z0 -G150/150/255 -K -O -V -F >> noise.ps
awk '{print $5}' t1 | gmt pshistogram -JX -R -T$bw  -Z0 -G200/200/255 -K -O -V -F >> noise.ps

awk '{print $6}' t1 | gmt pshistogram -JX -R -T$bw  -Z0 -G255/0/0 -K -O -V -F >> noise.ps
awk '{print $7}' t1 | gmt pshistogram -JX -R -T$bw  -Z0 -G255/80/80 -K -O -V -F >> noise.ps
awk '{print $8}' t1 | gmt pshistogram -JX -R -T$bw  -Z0 -G255/150/150 -K -O -V -F >> noise.ps
awk '{print $9}' t1 | gmt pshistogram -JX -R -T$bw  -Z0 -G255/200/200 -K -O -V -F >> noise.ps

echo 0 0 | gmt psxy -JX -R -B0 -Sc0.001 -O >> noise.ps

gmt psconvert -Tg noise.ps -A
ls "$PWD/noise.p"*

[[ "$(uname)" == "Darwin" ]] && open noise.png
