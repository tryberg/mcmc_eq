#!/bin/bash

rm -f gmt.*
gmt set MEASURE_UNIT INCH
gmt set FONT_ANNOT_PRIMARY 10
gmt set HEADER_FONT_SIZE 10
gmt set HEADER_OFFSET 0.5c
gmt set LABEL_FONT_SIZE 10
gmt set COLOR_NAN 200/200/200
export LC_NUMERIC=C.UTF-8

exe="fw"
which "$exe"
if ! command -v "$exe" &> /dev/null; then
    echo "Can't find $exe"
    exit 1
fi

cfg="config_eqx.dat"
ls "$cfg"
if [[ ! -f "$cfg" ]]; then
    exit 1
fi

res="resmcnx.dat"
ls "$res"
if [[ ! -f "$res" ]]; then
    exit 1
fi

z0=$(grep Z0 "$cfg" | awk '{print $1}')
nz=$(grep NZ "$cfg" | awk '{print $1}')
dx=$(awk 'NR==1 {print $1}' "$cfg")

z1=$(echo "$z0 $nz $dx" | awk '{print $1+($2-2)*$3}')

echo "Reading $# inputs..."
if [[ $# -ne 1 ]]; then
    echo "Need mcmcPhaseFile (picks) as input arg"
    exit 1
fi

picks="$1"
ls "$picks"
if [[ ! -f "$picks" ]]; then
    exit 1
fi

egrep STAN "$res" > res.dat
egrep EQ "$res" | awk -v z0="$z0" -v z1="$z1" '{z=$5; if (z<1.0*z0) z=z0; if (z>1.0*z1) z=z1; print $1, $2, $3, $4, z, $6, $7, $8, $9, $10, $11}' >> res.dat
egrep EZ "$res" | awk -v z0="$z0" -v z1="$z1" '{z=$5; if (z<1.0*z0) z=z0; if (z>1.0*z1) z=z1; print $1, $2, $3, $4, z, $6, $7, $8, $9, $10, $11}' >> res.dat
egrep RES "$res" >> res.dat
egrep NOISE "$res" | awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17}' >> res.dat

# Prepare pick file with that quake
cat "$picks" > test

"$exe" "$cfg" res.dat test > predictions

awk '{print $0}' test > t0

paste predictions t0 > residuals.dat

awk '{if ($1!="EVENT") print $1}' residuals.dat > t
#~/bin/rms.linux t

x2=`grep -v "EVENT" residuals.dat | gmt info -i1 -C -o1 -I5`
gmt psbasemap -JX7/5 -R0/$x2/-2/15 -Bxaf+l"epi dist [km] quake" -Byafg100+l"t - epi/8.0 [s]" -BSWen -K -P > msftp.ps

awk '{if ($7=="P") print $2, $5-$4-$2/8}' residuals.dat | gmt psxy -JX -R -Sc0.025 -Gblue -K -O -V >> msftp.ps
awk '{if ($7=="P") print $2, $6-$2/8}' residuals.dat | gmt psxy -JX -R -S+0.05 -K -O -V >> msftp.ps
awk '{if ($7=="S") print $2, $5-$4-$2/8}' residuals.dat | gmt psxy -JX -R -Sc0.025 -Gred -K -O -V >> msftp.ps
awk '{if ($7=="S") print $2, $6-$2/8}' residuals.dat | gmt psxy -JX -R -S+0.05 -K -O -V >> msftp.ps

xmax=`gmt info t -C -I1 -o1`
gmt psbasemap -JX7/4 -R-$xmax/$xmax/1/500 -Bxafg1000+l"misfit dt [s]" -Byaf -BSWen -K -O -Y6 >> msftp.ps

awk '{if (($1!="EVENT") && ($7=="P") && ($15!=10)) print $1}' residuals.dat | gmt pshistogram -JX -R -B0 -W0.01 -L1,blue -K -O -F -S >> msftp.ps
awk '{if (($1!="EVENT") && ($7=="S") && ($15!=10)) print $1}' residuals.dat | gmt pshistogram -JX -R -B0 -W0.01 -L1,red -K -O -F -S >> msftp.ps

awk '{if ($7=="P") print $1}' residuals.dat > t1
prms=`awk '{t=t+$1*$1; n=n+1; printf"%5.3f\n", sqrt(t/n);}' t1 | tail -1`
awk '{if ($7=="S") print $1}' residuals.dat > t1
srms=`awk '{t=t+$1*$1; n=n+1; printf"%5.3f\n", sqrt(t/n);}' t1 | tail -1`
noq=`egrep EVENT residuals.dat  | wc | awk '{print $1}'`
echo $prms $srms $noq

echo $xmax 450 "P_misfit="$prms" s" | gmt pstext -JX -R -B0 -K -Wblue  -O -V -F+f10p+jTR -D-.1c >> msftp.ps
echo $xmax 400 "S_misfit="$srms" s" | gmt pstext -JX -R -B0 -K -Wred   -O -V -F+f10p+jTR -D-.1c >> msftp.ps
echo $xmax 350 ""$noq" quakes"      | gmt pstext -JX -R -B0 -K -Wblack -O -V -F+f10p+jTR -D-.1c >> msftp.ps

# Pick classes
# Uncomment the following lines if needed
#awk '{if (($1!="EVENT") && ($7=="P") && ($15==0)) print $1}' residuals.dat | pshistogram -JX -R -B0 -W0.01 -L1/0/0/255 -K -O -F -S >> msftp.ps
#awk '{if (($1!="EVENT") && ($7=="P") && ($15==1)) print $1}' residuals.dat | pshistogram -JX -R -B0 -W0.01 -L1/42/42/255 -K -O -F -S >> msftp.ps
#awk '{if (($1!="EVENT") && ($7=="P") && ($15==2)) print $1}' residuals.dat | pshistogram -JX -R -B0 -W0.01 -L1/85/85/255 -K -O -F -S >> msftp.ps
#awk '{if (($1!="EVENT") && ($7=="P") && ($15==3)) print $1}' residuals.dat | pshistogram -JX -R -B0 -W0.01 -L1/128/128/255 -K -O -F -S >> msftp.ps

echo 0 0 | gmt psxy -JX -R -B0 -Sc0.001 -O >> msftp.ps

rm -f t residuals.dat predictions res.dat

gmt ps2raster -Tg msftp.ps -Au
ls "$PWD/msftp.p"*

if [[ "$(uname)" == "Darwin" ]]; then
    open msftp.png
fi
