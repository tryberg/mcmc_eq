#!/bin/bash

# see also disp_error in scripts
# used to compare to synthetic inputs or previous/original results
# input files are outputs from mcmc_eq and pha2mcmc.m (or other..)

rm -f gmt.*
gmt set MEASURE_UNIT INCH
gmt set FONT_ANNOT_PRIMARY 8 
gmt set HEADER_FONT_SIZE 10
gmt set HEADER_OFFSET 0.5c
gmt set LABEL_FONT_SIZE 10
gmt set COLOR_NAN  200/200/200

res=resmcnx.dat
ls "$res"
[[ -f "$res" ]] || exit

egrep EZ $res > t1

quakes=quakes.dat # if prior locations file exists
ls $quakes
[[ -f "$quakes" ]] || exit

ls stations.dat
[[ -f "stations.dat" ]] || exit

ls model.inp
[[ -f "model.inp" ]] || echo "No input model to plot (model.inp), skip"

ls tmpx 
[[ -f  tmpx ]] || { echo "need to run disp_m_average_sl.sh first" ; exit 1 ; }

out=compare

paste t1 $quakes | awk '{print $3-$14,$4-$15,$5-$16}' > t3
xmax=`awk '{for(i=1;i<=3;i++){v=($i<0)?-$i:$i; if(v>m)m=v}} END{print m}' t3 | gmt info -C -I1 -o1`
bw=0.05
ymax=`gmt pshistogram t3 -i0 -T$bw -I -o3 | gmt info -I10 -C -o1`
echo $ymax $xmax

echo EQ x
gmt psbasemap -JX2/2 -R-$xmax/$xmax/0/$ymax -Bxa+l"Quake_x dx [km]" -Bya -BSWne -K -P -Y8 -X0.5 > $out.ps
paste t1 $quakes | awk '{print $3-$14}' > t2
gmt pshistogram t2 -JX -R -T$bw -Z0 -G0 -K -O -V -F >> $out.ps

m=$(awk '{i++; s+=$1;} END {printf "%5.2f\n", s/i;}' t2)
s=$(awk -v mean="$m" '{i++; s+=($1-mean)*($1-mean);} END {printf "%5.2f\n", sqrt(s/(i-1));}' t2)

echo -$xmax $ymax  "m/s="$m"+-"$s" km" | gmt pstext -JX -R -K -O -N -F+jTL -D0.1c/-0.1c -Gwhite >> $out.ps

echo EQ y
gmt psbasemap -JX2/2 -R-$xmax/$xmax/0/$ymax -Bxa+l"Quake_y dy [km]" -Bya -BSwen -K -O -X2.35 >> $out.ps
paste t1 $quakes | awk '{print $4-$15}' > t2
gmt pshistogram t2 -JX -R -T$bw -Z0 -G0 -K -O -V -F >> $out.ps

m=$(awk '{i++; s+=$1;} END {printf "%5.2f\n", s/i;}' t2)
s=$(awk -v mean="$m" '{i++; s+=($1-mean)*($1-mean);} END {printf "%5.2f\n", sqrt(s/(i-1));}' t2)

echo -$xmax $ymax  "m/s="$m"+-"$s" km" | gmt pstext -JX -R -K -O -N -F+jTL -D0.1c/-0.1c -Gwhite >> $out.ps

echo 0 0 | gmt psxy -JX -R -B0:"`pwd`":/0N -Sc0.001 -O -K >> $out.ps

echo EQ z
gmt psbasemap -JX2/2 -R-$xmax/$xmax/0/$ymax -Bxa+l"shallower    dz [km]     deeper" -Bya -BSwen -K -O -X2.35 >> $out.ps
paste t1 $quakes | awk '{print $5-$16}' > t2
gmt pshistogram t2 -JX -R -T$bw -Z0 -G0 -K -O -V -F >> $out.ps

m=$(awk '{i++; s+=$1;} END {printf "%5.2f\n", s/i;}' t2)
s=$(awk -v mean="$m" '{i++; s+=($1-mean)*($1-mean);} END {printf "%5.2f\n", sqrt(s/(i-1));}' t2)

echo -$xmax $ymax  "m/s="$m"+-"$s" km" | gmt pstext -JX -R -K -O -N -F+jTL -D0.1c/-0.1c -Gwhite >> $out.ps

echo origin time
paste t1 $quakes | awk '{print $9+$10-$17}' > t2
bw=0.025
ymax=`gmt pshistogram t2 -i0 -T$bw -I -o3 | gmt info -I10 -C -o1`
#xmax=`gmt info t2 -I0.5 -C -o1a`
xmax=1
gmt psbasemap -JX2/2 -R-$xmax/$xmax/0/$ymax -Bxa+l"Origin time dt [s]" -Bya -BSWen -K -O -Y-3.0 -X-4.7 >> $out.ps
paste t1 $quakes | awk '{print $9+$10-$17}' | gmt pshistogram -JX -R -T$bw -Z0 -G0 -K -O -V -F >> $out.ps

m=$(awk '{i++; s+=$1;} END {printf "%5.2f\n", s/i;}' t2)
s=$(awk -v mean="$m" '{i++; s+=($1-mean)*($1-mean);} END {printf "%5.2f\n", sqrt(s/(i-1));}' t2)

echo -$xmax $ymax  "m/s="$m"+-"$s" s" | gmt pstext -JX -R -K -O -N -F+jTL -D0.1c/-0.1c -Gwhite >> $out.ps

echo data noise
bw=0.0005
ymax=`egrep mod tmpx | awk '{print $6}' | gmt pshistogram -T$bw -I -o3 | gmt info -I10 -C -o1`
egrep mod tmpx | awk '{print $6,$7,$8,$9,$10,$11,$12,$13}' > t3
xmax=`awk '{for(i=1;i<=NF;i++){v=($i<0)?-$i:$i; if(v>m)m=v}} END{print m}' t3 | gmt info -C -I0.2 -o1`
xmax=$(awk -v a="$xmax" -v b="1.5" 'BEGIN{print (a<b)?a:b}')

gmt psbasemap -JX2/2 -R0.0/$xmax/0/$ymax -Bxa+l"Data noise [s]" -Bya -BSwen -K -O -X2.35 >> $out.ps

# true values
echo 0.025 | awk '{print $1, 0; print $1, '$ymax'}' | gmt psxy -JX -R -W2,200/200/255 -O -K -N >> $out.ps
echo 0.050 | awk '{print $1, 0; print $1, '$ymax'}' | gmt psxy -JX -R -W2,200/200/255 -O -K -N >> $out.ps
echo 0.075 | awk '{print $1, 0; print $1, '$ymax'}' | gmt psxy -JX -R -W2,200/200/255 -O -K -N >> $out.ps
echo 0.100 | awk '{print $1, 0; print $1, '$ymax'}' | gmt psxy -JX -R -W2,200/200/255 -O -K -N >> $out.ps

echo 0.0875 | awk '{print $1, 0; print $1, '$ymax'}' | gmt psxy -JX -R -W2,255/200/200 -O -K -N >> $out.ps
echo 0.1125 | awk '{print $1, 0; print $1, '$ymax'}' | gmt psxy -JX -R -W2,255/200/200 -O -K -N >> $out.ps
echo 0.1375 | awk '{print $1, 0; print $1, '$ymax'}' | gmt psxy -JX -R -W2,255/200/200 -O -K -N >> $out.ps
echo 0.1625 | awk '{print $1, 0; print $1, '$ymax'}' | gmt psxy -JX -R -W2,255/200/200 -O -K -N >> $out.ps

echo "start noise"
s=$(egrep mod tmpx | awk '{print $6}' | gmt math STDIN MEAN = /dev/stdout | head -1)
echo $s | awk '{print $1, 0; print $1, '$ymax'}' | gmt psxy -JX -R -W1,0/0/255,dotted -O -K -N >> $out.ps
s=$(egrep mod tmpx | awk '{print $7}' | gmt math STDIN MEAN = /dev/stdout | head -1)
echo $s | awk '{print $1, 0; print $1, '$ymax'}' | gmt psxy -JX -R -W1,42/42/255,dotted -O -K -N >> $out.ps
s=$(egrep mod tmpx | awk '{print $8}' | gmt math STDIN MEAN = /dev/stdout | head -1)
echo $s | awk '{print $1, 0; print $1, '$ymax'}' | gmt psxy -JX -R -W1,85/85/255,dotted -O -K -N >> $out.ps
s=$(egrep mod tmpx | awk '{print $9}' | gmt math STDIN MEAN = /dev/stdout | head -1)
echo $s | awk '{print $1, 0; print $1, '$ymax'}' | gmt psxy -JX -R -W1,128/128/255,dotted -O -K -N >> $out.ps
egrep mod tmpx | awk '{print $6}' | gmt pshistogram -JX -R -T$bw -Z0 -G0/0/255 -K -O -F -V >> $out.ps
egrep mod tmpx | awk '{print $7}' | gmt pshistogram -JX -R -T$bw -Z0 -G42/42/255 -K -O -F -V >> $out.ps
egrep mod tmpx | awk '{print $8}' | gmt pshistogram -JX -R -T$bw -Z0 -G85/85/255 -K -O -F -V >> $out.ps
egrep mod tmpx | awk '{print $9}' | gmt pshistogram -JX -R -T$bw -Z0 -G128/128/255 -K -O -F -V >> $out.ps

s=$(egrep mod tmpx | awk '{print $10}' | gmt math STDIN MEAN = /dev/stdout | head -1)
echo $s | awk '{print $1, 0; print $1, '$ymax'}' | gmt psxy -JX -R -W1,255/0/0,dotted -O -K -N >> $out.ps
s=$(egrep mod tmpx | awk '{print $11}' | gmt math STDIN MEAN = /dev/stdout | head -1)
echo $s | awk '{print $1, 0; print $1, '$ymax'}' | gmt psxy -JX -R -W1,255/42/42,dotted -O -K -N >> $out.ps
s=$(egrep mod tmpx | awk '{print $12}' | gmt math STDIN MEAN = /dev/stdout | head -1)
echo $s | awk '{print $1, 0; print $1, '$ymax'}' | gmt psxy -JX -R -W1,255/85/85,dotted -O -K -N >> $out.ps
s=$(egrep mod tmpx | awk '{print $13}' | gmt math STDIN MEAN = /dev/stdout | head -1)
echo $s | awk '{print $1, 0; print $1, '$ymax'}' | gmt psxy -JX -R -W1,255/128/128,dotted -O -K -N >> $out.ps
egrep mod tmpx | awk '{print $10}' | gmt pshistogram -JX -R -T$bw -Z0 -G255/0/0 -K -O -F -V >> $out.ps
egrep mod tmpx | awk '{print $11}' | gmt pshistogram -JX -R -T$bw -Z0 -G255/42/42 -K -O -F -V >> $out.ps
egrep mod tmpx | awk '{print $12}' | gmt pshistogram -JX -R -T$bw -Z0 -G255/85/85 -K -O -F -V >> $out.ps
egrep mod tmpx | awk '{print $13}' | gmt pshistogram -JX -R -T$bw -Z0 -G255/128/128 -K -O -F -V >> $out.ps
echo "end noise"

# station correction P
egrep RES $res > t1
paste t1 stations.dat | awk '{print $3-$11}' > t2
bw=0.01
xmax=`gmt info t2 -C -I1 -o1`
test "$xmax" -eq "0" && xmax=1
ymax=`gmt pshistogram t2 -T$bw -I -o3 | gmt info -I5 -C -o1`

gmt psbasemap -JX2/2 -R-$xmax/$xmax/0/$ymax -Bxa+l"Station correction P dt [s]" -Bya -BSWen -K -O -Y-3.0 -X-2.35 >> $out.ps
gmt pshistogram t2 -JX -R -T$bw -Z0 -G0 -K -O -V -F >> $out.ps

m=$(awk '{i++; s+=$1;} END {printf "%5.2f\n", s/i;}' t2)
s=$(awk -v mean="$m" '{i++; s+=($1-mean)*($1-mean);} END {printf "%5.2f\n", sqrt(s/(i-1));}' t2)

echo -$xmax $ymax  "m/s="$m"+-"$s" s" | gmt pstext -JX -R -K -O -N -F+jTL -D0.1c/-0.1c -Gwhite  >> $out.ps

# station correction S
paste t1 stations.dat | awk '{print $4-$12}' > t2
xmax=`gmt info t2 -C -I1 -o1`
test "$xmax" -eq "0" && xmax=1
ymax=`gmt pshistogram t2 -T$bw -I -o3 | gmt info -I5 -C -o1`
gmt psbasemap -JX2/2 -R-$xmax/$xmax/0/$ymax -Bxa+l"Station correction S dt [s]" -Bya -BSwen -K -O -X2.35 >> $out.ps
gmt pshistogram t2 -JX -R -T$bw -Z0 -G0 -K -O -V -F >> $out.ps

m=$(awk '{i++; s+=$1;} END {printf "%5.2f\n", s/i;}' t2)
s=$(awk -v mean="$m" '{i++; s+=($1-mean)*($1-mean);} END {printf "%5.2f\n", sqrt(s/(i-1));}' t2)

echo -$xmax $ymax  "m/s="$m"+-"$s" s" | gmt pstext -JX -R -K -O -N -F+jTL -D0.1c/-0.1c -Gwhite  >> $out.ps

# Vp
gmt psbasemap -JX1.2/-5 -R3/8/-5/32 -B1f0.5:"Vp [km/s]":/10f5g1000:"Depth [km]":Swen -K -O -X2.25 -Y0 >> $out.ps
egrep EZ $res | awk '{print 4, $5}' | gmt psxy -JX -R -Sc0.05 -Ggreen -K -O >> $out.ps
awk '{print 3.5,$4}' $quakes | gmt psxy -JX -R -Sc0.05 -Gdarkgreen -K -O >> $out.ps

test -f model.inp && awk '{print $2, $1}' model.inp | gmt psxy -JX -R -W1,magenta -K -O >> $out.ps
test -f vp_generic.xy && gmt psxy vp_generic.xy -JX -R -W1,blue -K -O >> $out.ps
awk '{if ($1=="STAN") print $7, $2}' $res | awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' | gmt psxy -JX -R -W1,red -O -K -N  >> $out.ps
awk '{if ($1=="STAN") print $7-$8, $2}' $res | awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' | gmt psxy -JX -R -W,gray -O -K -N  >> $out.ps
awk '{if ($1=="STAN") print $7+$8, $2}' $res | awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' | gmt psxy -JX -R -W,gray -O -K -N  >> $out.ps

# Vp/Vs
gmt psbasemap -JX1.2/-5 -R1.501/2/-5/32 -B0.2f0.1:"Vp/Vs":/10f5g1000:"Depth [km]":SwEn -K -O -X1.35  >> $out.ps
test -f model.inp && awk '{print $3, $1}' model.inp | gmt psxy -JX -R -W1,blue -K -O >> $out.ps
awk '{if ($1=="STAN") print $9, $2}' $res | awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' | gmt psxy -JX -R -W1,red -O -K  >> $out.ps
awk '{if ($1=="STAN") print $9-$10, $2}' $res | awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' | gmt psxy -JX -R -W,gray -O -K  >> $out.ps
awk '{if ($1=="STAN") print $9+$10, $2}' $res | awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' | gmt psxy -JX -R -W,gray -O -K  >> $out.ps

echo 0 0 | gmt psxy -JX -R -B0 -Sc0.001 -O >> $out.ps

gmt psconvert -Tg $out.ps -A
ls "$PWD/$out.p"*

[[ "$(uname)" == "Darwin" ]] && open $out.png

