#!/usr/bin/env bash
#
# disp_staCorMap.sh
#
# Plots station correction maps (P-wave and S-wave, separately) from a
# VELEST/MCMC-style residual file, using GMT 6 modern mode. Corrections
# are shown as color-filled circles on a diverging red-blue scale
# (red = negative, blue = positive), sized by magnitude.
#
# Usage:
#   ./disp_staCorMap.sh mcmcPhaseFile

set -euo pipefail

cfg=config_eqx.dat
res=resmcnx.dat

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 mcmcPhaseFile (picks)" >&2
    exit 1
fi

picks="$1"

for f in "$cfg" "$res" "$picks"; do
    if [[ ! -f "$f" ]]; then
        echo "Error: required file not found: $f" >&2
        exit 1
    fi
done

# --- Read map extent from config file ---
d=$(awk 'NR==1 {print $1}' "$cfg")
nx=$(awk 'NR==2 {print $1}' "$cfg")
xmin=$(awk 'NR==5 {print $1}' "$cfg")
xmax=$(echo "$xmin $nx $d" | awk '{print $1+$2*$3}')

ny=$(awk 'NR==3 {print $1}' "$cfg")
ymin=$(awk 'NR==6 {print $1}' "$cfg")
ymax=$(echo "$ymin $ny $d" | awk '{print $1+$2*$3}')

# --- Map scaling ---
xyr=$(echo "$xmin $xmax $ymin $ymax" | awk '{print ($2 - $1)/($4 - $3)}')
ys=4.5
xs=$(echo "$ys $xyr" | awk '{print $1*$2}')

# --- Build working data: station corrections (P, S) + station x,y ---
tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

t1="$tmpdir/t1.txt"
rec="$tmpdir/rec.dat"
recdata="$tmpdir/recdata.txt"

awk '{if ($1=="RES") print $2, $3, $4}' "$res" > "$t1"
awk '{if ($1!="#" && $2!="NA") print $1, $2, $4, $5, $6}' "$picks" | sort -u > "$rec"
paste "$t1" "$rec" > "$recdata"

# recdata columns: 1=stationID(t1) 2=P_corr 3=S_corr 4=stationID(rec) 5=? 6=X 7=Y

# --- Determine a shared symmetric color scale across BOTH P and S, so
#     the two maps are directly comparable ---
absmax=$(awk '{
    p = ($2<0)?-$2:$2
    s = ($3<0)?-$3:$3
    if (p>m) m=p
    if (s>m) m=s
} END { print m }' "$recdata")

# pad slightly and guard against a degenerate all-zero case
cptmax=$(echo "$absmax" | awk '{v=$1*1.1; if (v<=0) v=1; print v}')

cpt="$tmpdir/corr.cpt"
gmt makecpt -Cpolar -T-"${cptmax}"/"${cptmax}"/0.01 -D -Z > "$cpt"

# --- Function to plot one map (P or S) ---
# args: column_number  label  output_basename
plot_map() {
    local col="$1"
    local label="$2"
    local outname="$3"

    gmt begin "$outname" png
        gmt basemap -JX"${xs}i/${ys}i" -R"${xmin}/${xmax}/${ymin}/${ymax}" \
            -Bxaf+l"X [km]" -Byaf+l"Y [km]" -BNWse

        # Skip plotting circles entirely if every value in this column is zero
        nonzero=$(awk -v c="$col" '{v=$c; if (v<0) v=-v; if (v>0) {print "yes"; exit}}' "$recdata")
        if [[ -n "$nonzero" ]]; then
            awk -v c="$col" '{
                v = $c
                av = (v<0)?-v:v
                size = av*0.8
                if (size < 0.02) size = 0.02
                print $6, $7, v, size
            }' "$recdata" | sort -k4,4 -rn | gmt plot -Sc -C"$cpt" -W0.5p,black
        else
            echo "  (all ${label}-wave values are zero -- skipping circles)"
        fi

        # Station name labels (column 4 of recdata), offset from each point
        awk '{ print $6, $7, $4 }' "$recdata" | gmt text -F+f7p,Helvetica+jLM -D0.15c/0c

        gmt colorbar -C"$cpt" -Dx"$(echo "$xs" | awk '{print $1/2}')i/-0.6i+w${xs}i/0.15i+h+jTC" -Bxaf+l"${label} correction"

        gmt text -F+f10p,Helvetica-Bold+jTL -D0.2c/-0.2c << EOF
${xmin} ${ymax} ${label}-wave Station Corrections
EOF
    gmt end
}

plot_map 2 "P" "staCorMap_P"
plot_map 3 "S" "staCorMap_S"

echo "Wrote staCorMap_P.png and staCorMap_S.png"

if [[ "$(uname)" == "Darwin" ]]; then
    open staCorMap_P.png staCorMap_S.png
fi
