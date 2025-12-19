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

# EQ x
gmt psbasemap -JX2/2 -R-2/2/0/100 -Bxa+l"Quake_x dx [km]" -Bya -BSWne -K -P -Y8 -X0.5 > $out.ps
paste t1 $quakes | awk '{print $3-$14}' | gmt pshistogram -JX -R -W0.1 -Z0 -G0 -K -O -V -F >> $out.ps
paste t1 $quakes | awk '{print $3-$14}' > t2

m=$(awk '{i++; s+=$1;} END {printf "%5.2f\n", s/i;}' t2)
s=$(awk -v mean="$m" '{i++; s+=($1-mean)*($1-mean);} END {printf "%5.2f\n", sqrt(s/(i-1));}' t2)

echo -2 90  "m/s="$m"+-"$s" km" | gmt pstext -JX -R -K -O -N -F+jBL >> $out.ps

# EQ y
gmt psbasemap -JX2/2 -R-2/2/0/100 -Bxa+l"Quake_y dy [km]" -Bya -BSwen -K -O -X2.35 >> $out.ps
paste t1 $quakes | awk '{print $4-$15}' | gmt pshistogram -JX -R -W0.1 -Z0 -G0 -K -O -V -F >> $out.ps
paste t1 $quakes | awk '{print $4-$15}' > t2

m=$(awk '{i++; s+=$1;} END {printf "%5.2f\n", s/i;}' t2)
s=$(awk -v mean="$m" '{i++; s+=($1-mean)*($1-mean);} END {printf "%5.2f\n", sqrt(s/(i-1));}' t2)

echo -2 90  "m/s="$m"+-"$s" km" | gmt pstext -JX -R -K -O -N -F+jBL >> $out.ps

echo 0 0 | gmt psxy -JX -R -B0:"`pwd`":/0N -Sc0.001 -O -K >> $out.ps

# EQ z
gmt psbasemap -JX2/2 -R-2/2/0/100 -Bxa+l"shallower    dz [km]     deeper" -Bya -BSwen -K -O -X2.35 >> $out.ps
paste t1 $quakes | awk '{print $5-$16}' | gmt pshistogram -JX -R -W0.1 -Z0 -G0 -K -O -V -F >> $out.ps
paste t1 $quakes | awk '{print $5-$16}' > t2

m=$(awk '{i++; s+=$1;} END {printf "%5.2f\n", s/i;}' t2)
s=$(awk -v mean="$m" '{i++; s+=($1-mean)*($1-mean);} END {printf "%5.2f\n", sqrt(s/(i-1));}' t2)

echo -2 90  "m/s="$m"+-"$s" km" | gmt pstext -JX -R -K -O -N -F+jBL >> $out.ps

#origin time
gmt psbasemap -JX2/2 -R-2/2/0/100 -Bxa+l"Origin time dt [s]" -Bya -BSWen -K -O -Y-3.0 -X-4.7 >> $out.ps
paste t1 $quakes | awk '{print $9+$10-$17}' | gmt pshistogram -JX -R -W0.025 -Z0 -G0 -K -O -V -F >> $out.ps
paste t1 $quakes | awk '{print $9+$10-$17}' > t2

m=$(awk '{i++; s+=$1;} END {printf "%5.2f\n", s/i;}' t2)
s=$(awk -v mean="$m" '{i++; s+=($1-mean)*($1-mean);} END {printf "%5.2f\n", sqrt(s/(i-1));}' t2)

echo -2 90  "m/s="$m"+-"$s" s" | gmt pstext -JX -R -K -O -N -F+jBL >> $out.ps

# datat noise
gmt psbasemap -JX2/2 -R0.0/0.15/0/1000 -Bxa+l"Data noise [s]" -Bya -BSwen -K -O -X2.35 >> $out.ps

# true values
echo 0.015 | awk '{print $1, 0; print $1, 1000}' | gmt psxy -JX -R -W2,200/200/255 -O -K -N >> $out.ps
echo 0.030 | awk '{print $1, 0; print $1, 1000}' | gmt psxy -JX -R -W2,200/200/255 -O -K -N >> $out.ps
echo 0.045 | awk '{print $1, 0; print $1, 1000}' | gmt psxy -JX -R -W2,200/200/255 -O -K -N >> $out.ps
echo 0.060 | awk '{print $1, 0; print $1, 1000}' | gmt psxy -JX -R -W2,200/200/255 -O -K -N >> $out.ps

echo 0.0525 | awk '{print $1, 0; print $1, 1000}' | gmt psxy -JX -R -W2,255/200/200 -O -K -N >> $out.ps
echo 0.0675 | awk '{print $1, 0; print $1, 1000}' | gmt psxy -JX -R -W2,255/200/200 -O -K -N >> $out.ps
echo 0.0825 | awk '{print $1, 0; print $1, 1000}' | gmt psxy -JX -R -W2,255/200/200 -O -K -N >> $out.ps
echo 0.0975 | awk '{print $1, 0; print $1, 1000}' | gmt psxy -JX -R -W2,255/200/200 -O -K -N >> $out.ps

echo "start noise"
egrep mod tmpx | awk '{print $6}' | gmt pshistogram -JX -R -W0.0005 -Z0 -G0/0/255 -K -O -F -V >> $out.ps
egrep mod tmpx | awk '{print $7}' | gmt pshistogram -JX -R -W0.0005 -Z0 -G42/42/255 -K -O -F -V >> $out.ps
egrep mod tmpx | awk '{print $8}' | gmt pshistogram -JX -R -W0.0005 -Z0 -G85/85/255 -K -O -F -V >> $out.ps
egrep mod tmpx | awk '{print $9}' | gmt pshistogram -JX -R -W0.0005 -Z0 -G128/128/255 -K -O -F -V >> $out.ps

egrep mod tmpx | awk '{print $10}' | gmt pshistogram -JX -R -W0.0005 -Z0 -G255/0/0 -K -O -F -V >> $out.ps
egrep mod tmpx | awk '{print $11}' | gmt pshistogram -JX -R -W0.0005 -Z0 -G255/42/42 -K -O -F -V >> $out.ps
egrep mod tmpx | awk '{print $12}' | gmt pshistogram -JX -R -W0.0005 -Z0 -G255/85/85 -K -O -F -V >> $out.ps
egrep mod tmpx | awk '{print $13}' | gmt pshistogram -JX -R -W0.0005 -Z0 -G255/128/128 -K -O -F -V >> $out.ps
echo "end noise"

# station correction P
egrep RES $res > t1
gmt psbasemap -JX2/2 -R-0.5/0.5/0/15 -Bxa+l"Station correction P dt [s]" -Bya -BSWen -K -O -Y-3.0 -X-2.35 >> $out.ps
paste t1 stations.dat | awk '{print $3-$11}' | gmt pshistogram -JX -R -W0.01 -Z0 -G0 -K -O -V -F >> $out.ps
paste t1 stations.dat | awk '{print $3-$11}' > t2

m=$(awk '{i++; s+=$1;} END {printf "%5.2f\n", s/i;}' t2)
s=$(awk -v mean="$m" '{i++; s+=($1-mean)*($1-mean);} END {printf "%5.2f\n", sqrt(s/(i-1));}' t2)

echo -0.5 13  "m/s="$m"+-"$s" s" | gmt pstext -JX -R -K -O -N -F+jBL >> $out.ps

# station correction S
gmt psbasemap -JX2/2 -R-0.5/0.5/0/15 -Bxa+l"Station correction S dt [s]" -Bya -BSwen -K -O -X2.35 >> $out.ps
paste t1 stations.dat | awk '{print $4-$12}' | gmt pshistogram -JX -R -W0.01 -Z0 -G0 -K -O -V -F >> $out.ps
paste t1 stations.dat | awk '{print $4-$12}' > t2

m=$(awk '{i++; s+=$1;} END {printf "%5.2f\n", s/i;}' t2)
s=$(awk -v mean="$m" '{i++; s+=($1-mean)*($1-mean);} END {printf "%5.2f\n", sqrt(s/(i-1));}' t2)

echo -0.5 13  "m/s="$m"+-"$s" s" | gmt pstext -JX -R -K -O -N -F+jBL >> $out.ps

# Vp
gmt psbasemap -JX1.2/-5 -R3/8/-5/35 -B1f0.5:"Vp [km/s]":/10f5g1000:"Depth [km]":Swen -K -O -X2.25 -Y0 >> $out.ps
egrep EZ $res | awk '{print 4, $5}' | gmt psxy -JX -R -Sc0.05 -Ggreen -K -O >> $out.ps
awk '{print 3.5,$4}' $quakes | gmt psxy -JX -R -Sc0.05 -Gdarkgreen -K -O >> $out.ps

test -f model.inp && awk '{print $2, $1}' model.inp | gmt psxy -JX -R -W1,blue -K -O -N >> $out.ps
awk '{if ($1=="STAN") print $7, $2}' $res | awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' | gmt psxy -JX -R -W1,red -O -K -N  >> $out.ps
awk '{if ($1=="STAN") print $7-$8, $2}' $res | awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' | gmt psxy -JX -R -W,gray -O -K -N  >> $out.ps
awk '{if ($1=="STAN") print $7+$8, $2}' $res | awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' | gmt psxy -JX -R -W,gray -O -K -N  >> $out.ps

# Vp/Vs
gmt psbasemap -JX1.2/-5 -R1.501/2/-5/35 -B0.2f0.1:"Vp/Vs":/10f5g1000:"Depth [km]":SwEn -K -O -X1.35  >> $out.ps
test -f model.inp && awk '{print $3, $1}' model.inp | gmt psxy -JX -R -W1,blue -K -O >> $out.ps
awk '{if ($1=="STAN") print $9, $2}' $res | awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' | gmt psxy -JX -R -W1,red -O -K  >> $out.ps
awk '{if ($1=="STAN") print $9-$10, $2}' $res | awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' | gmt psxy -JX -R -W,gray -O -K  >> $out.ps
awk '{if ($1=="STAN") print $9+$10, $2}' $res | awk '{if (NR==1) {v0=$1;} print v0, $2; print $1, $2; v0=$1;}' | gmt psxy -JX -R -W,gray -O -K  >> $out.ps

echo 0 0 | gmt psxy -JX -R -B0 -Sc0.001 -O >> $out.ps

gmt psconvert -Tg $out.ps -A
ls "$PWD/$out.p"*

[[ "$(uname)" == "Darwin" ]] && open $out.png

