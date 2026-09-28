#!/bin/bash

# disp_eq_z_ellipse.sh
# Based on disp_eq_z.sh (V2, gmt6). Adds, on top of the windowed MCMC posterior-
# sample density map for each event, the FORMAL spatial error ellipse derived from
# the new gmatrix covariance analysis:
#     Cm = (Ge^T Cd^-1 Ge)^-1     (per-event 4x4 [x y z ot] location covariance)
# The ellipse uses the off-diagonal covariance, so it is correctly oriented and
# shows the x-y / x-z trade-off -- information the axis-aligned marginal std-devs
# in resmcnx.dat do not carry.
#
# Requires (in the run directory): config_eqx.dat, resmcnx.dat, tmpx, and a
# loc-mode gmatrix output. If the gmatrix output is absent, this script will try
# to build it, needing the chain file and pick file (GMATRIX_CHAIN/GMATRIX_PICKS).
#
# J. Pesicek / Kiro, 2026

cfg="config_eqx.dat"
if [ ! -f "$cfg" ]
then
ls $cfg
exit
fi

# --- locate the helper + gmatrix binary relative to this script -------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC_DIR="$(cd "$SCRIPT_DIR/../src" 2>/dev/null && pwd)"
ELLIPSE_PY="$SRC_DIR/gmatrix_event_ellipse.py"
GMATRIX_BIN="$SRC_DIR/gmatrix"

# loc-mode gmatrix output (override with env GMATRIX_LOC); chain/picks for auto-build
GMATRIX_LOC="${GMATRIX_LOC:-gmatrix_loc.out}"
GMATRIX_CHAIN="${GMATRIX_CHAIN:-}"
GMATRIX_PICKS="${GMATRIX_PICKS:-picks.mcmc}"
GMATRIX_TAG="${GMATRIX_TAG:-bat}"
# error-ellipse confidence (2 DOF). Override with env: CONF=0.68 etc. or NSIGMA=1
CONF="${CONF:-0.95}"
NSIGMA="${NSIGMA:-}"

z0=$(awk 'NR==7 {print $1}' "$cfg")
eq=$(awk '{if (NR==30) print $1}' "$cfg") # needs to be past burn-in phase!
vv=$(awk '{if (NR==30) print $2}' "$cfg")
tot=$(echo "$vv $eq" | awk '{print $1+$2}')

str="<start model number post burn-in, b/t $eq and $tot> [<event number[s] to plot (default#: 100)>] [<window half length (5 km)>]"
me=`basename "$0"`

if [ "$#" -lt 1 ]
then
echo "${me}: $str"
echo "optionally plot multiple event numbers with quotes for 2nd input"
exit
else
echo post burn-in input number: $1
echo event ID[s]: $2
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

# --- ensure a loc-mode gmatrix output exists --------------------------------
if [ ! -f "$GMATRIX_LOC" ]; then
    echo "loc-mode gmatrix output '$GMATRIX_LOC' not found; attempting to build it"
    if [ -z "$GMATRIX_CHAIN" ]; then
        echo "  set GMATRIX_CHAIN=<chain file> (e.g. rjx-000.out) to auto-build,"
        echo "  or provide $GMATRIX_LOC yourself:"
        echo "     $GMATRIX_BIN $cfg <chain> $GMATRIX_PICKS loc $GMATRIX_TAG $GMATRIX_LOC"
        exit 1
    fi
    if [ ! -x "$GMATRIX_BIN" ]; then
        echo "  gmatrix binary not found at $GMATRIX_BIN (build with: make -C $SRC_DIR gmatrix)"
        exit 1
    fi
    "$GMATRIX_BIN" "$cfg" "$GMATRIX_CHAIN" "$GMATRIX_PICKS" loc "$GMATRIX_TAG" "$GMATRIX_LOC" \
        || { echo "gmatrix loc-mode build failed"; exit 1; }
fi

# confidence flag for the helper
if [ -n "$NSIGMA" ]; then
    ELL_ARGS="--nsigma $NSIGMA"
    ell_label="${NSIGMA}-sigma"
else
    ELL_ARGS="--conf $CONF"
    ell_label="$(echo "$CONF" | awk '{printf "%g%%", $1*100}')"
fi

rm -f gmt.*
gmt set MEASURE_UNIT INCH
gmt set HEADER_FONT_SIZE 10
gmt set FONT_ANNOT_PRIMARY 10
gmt set HEADER_OFFSET 0.5c
gmt set LABEL_FONT_SIZE 10
gmt set COLOR_NAN  200/200/200

export LC_NUMERIC=C.UTF-8

if [ "$#" -lt 3 ]
then
window=5 # radius around mean
else
window=$3
fi
echo "window length is 2*${window} km"

for eqn in $eqn0; do

    echo "plotting quake: $eqn"

    output="x.ps"

    xy=$(awk '{if (($2=="'"$eqn"'") && ($1=="EZ")) printf"%f/%f/%f/%f\n", $3-"'"$window"'", $3+"'"$window"'", $4-"'"$window"'", $4+"'"$window"'"}' "$f")
    xz=$(awk '{if (($2=="'"$eqn"'") && ($1=="EZ")) printf"%f/%f/%f/%f\n", $3-"'"$window"'", $3+"'"$window"'", $5-"'"$window"'", $5+"'"$window"'"}' "$f")
    x0=$(awk '{if (($2=="'"$eqn"'") && ($1=="EZ")) printf"%f\n", $3-"'"$window"'"}' "$f")
    x1=$(awk '{if (($2=="'"$eqn"'") && ($1=="EZ")) printf"%f\n", $3+"'"$window"'"}' "$f")

    dx=0.25
    dy=0.25
    dz=0.25

    awk '{if (($3>(1*"'$bi'")) && ($4=="'"$eqn"'")) print $0}' tmpx | grep EQ > t77

    # --- formal error ellipses from the new gmatrix analysis ----------------
    ell_xy=$(python3 "$ELLIPSE_PY" "$GMATRIX_LOC" "$eqn" --plane xy $ELL_ARGS --print-cov 2>ell.log)
    ell_xz=$(python3 "$ELLIPSE_PY" "$GMATRIX_LOC" "$eqn" --plane xz $ELL_ARGS 2>>ell.log)
    cat ell.log
    echo "ellipse xy (x y az major minor): $ell_xy"
    echo "ellipse xz (x z az major minor): $ell_xz"

    # x-y
    gmt psbasemap -JX5 -R$xy -B2f1:"X [km]":/2f1:"Y [km]":swEN -K -P -Y5.75 > "$output"
    awk '{print $6, $7}' t77 | \
    awk '{print int($1/"'"$dx"'"), int($2/"'"$dy"'")}' | \
    sort -n | \
    awk '{if (($1!=xold) || ($2!=yold)) {print xold*"'"$dx"'", yold*"'"$dy"'", s; s=0; xold=$1; yold=$2} else {s=s+1}}' | \
    tail -n +2 | gmt xyz2grd -Gtmpxy.grd -R -I"$dx/$dy" -V -F

    gmt grd2cpt -Chot -Z -D tmpxy.grd > tmp.cpt

    gmt grdimage tmpxy.grd -R -B0 -JX -Ctmp.cpt -K -O >> "$output"
    awk '{if (($2=="'"$eqn"'") && ($1=="EQ")) print $3, $4, $6, $7}' "$f" | \
    gmt psxy -JX -R -Sc0.075 -Gblue -W.5p,white -Exy+p0.5p,white -K -O -m >> "$output"

    # NEW: formal covariance error ellipse (from gmatrix loc-mode analysis)
    if [ -n "$ell_xy" ]; then
        echo "$ell_xy" | gmt psxy -JX -R -SE -W1.5p,cyan -K -O -N >> "$output"
    fi

awk '{print $6}' t77 > tjp
m=$(awk '{i++; s+=$1;} END {printf "%5.2f\n", s/i;}' tjp)
s=$(awk -v mean="$m" '{i++; s+=($1-mean)*($1-mean);} END {printf "%5.2f\n", sqrt(s/(i-1));}' tjp)
echo "X = $m +/- $s km" | gmt pstext -JX -R -K -O -N -F+cTL+jTL -D0.1i/-0.1i -Gwhite >> "$output"
awk '{print $7}' t77 > tjp
m=$(awk '{i++; s+=$1;} END {printf "%5.2f\n", s/i;}' tjp)
s=$(awk -v mean="$m" '{i++; s+=($1-mean)*($1-mean);} END {printf "%5.2f\n", sqrt(s/(i-1));}' tjp)
echo "Y = $m +/- $s km" | gmt pstext -JX -R -K -O -N -F+cTL+jTL -D0.1i/-0.3i -Gwhite >> "$output"
# ellipse legend
echo "cyan: ${ell_label} error ellipse (gmatrix Cm)" | \
  gmt pstext -JX -R -K -O -N -F+cTR+jTR+f8p,Helvetica,cyan4 -D-0.1i/-0.1i -Gwhite >> "$output"

    # x-z
    gmt psbasemap -JX5/-5 -R$xz -B2f1:"X [km] EQ $eqn":/2f1:"Z [km]":SwEn -K -Y-5.0 -O >> "$output"

    awk '{print $6, $8-"'"$z0"'"}' t77 | \
    awk '{print int($1/"'"$dx"'"), int($2/"'"$dy"'")}' | \
    sort -n | \
    awk '{if (($1!=xold) || ($2!=yold)) {print xold*"'"$dx"'", yold*"'"$dy"'"+"'"$z0"'", s; s=0; xold=$1; yold=$2} else {s=s+1}}' | \
    tail -n +2 | gmt xyz2grd -Gtmpxy.grd -R -I"$dx/$dy" -V -F

    gmt grd2cpt -Chot -Z -D tmpxy.grd > tmp.cpt
    gmt grdimage tmpxy.grd -R -B0 -JX -Ctmp.cpt -K -O >> "$output"
    awk '{if (($2=="'"$eqn"'") && ($1=="EQ")) print $3, $5, $6, $8}' "$f" | gmt psxy -JX -R -Sc0.075 -Gblue -W.5p,white -Exy+p0.5p,white -K -O -m >> "$output"
    awk '{if (($2=="'"$eqn"'") && ($1=="EZ")) print $3, $5, $6, $8}' "$f" | gmt psxy -JX -R -Sc0.075 -Ggreen -W.5p,white -Exy+p0.5p,white -K -O -m >> "$output"

    # NEW: formal covariance error ellipse (x-z plane)
    if [ -n "$ell_xz" ]; then
        echo "$ell_xz" | gmt psxy -JX -R -SE -W1.5p,cyan -K -O -N >> "$output"
    fi

    echo "$x0" "$x1" "$z0" | awk '{print $1, $3; print $2, $3; print ">" }' | gmt psxy -JX -R -W -K -O -m >> "$output"

awk '{print $8}' t77 > tjp
m=$(awk '{i++; s+=$1;} END {printf "%5.2f\n", s/i;}' tjp)
s=$(awk -v mean="$m" '{i++; s+=($1-mean)*($1-mean);} END {printf "%5.2f\n", sqrt(s/(i-1));}' tjp)
echo "Z = $m +/- $s km" | gmt pstext -JX -R -K -O -N -F+cTL+jTL -D0.1i/-0.1i -Gwhite >> "$output"

    echo 0 0 | gmt psxy -JX -R -B0 -Sc0.001 -O >> "$output"
    mv "$output" "loc_eq_ellipse_${eqn}.ps"
    gmt psconvert -Tg "loc_eq_ellipse_${eqn}.ps" -A
    ls "$PWD/loc_eq_ellipse_${eqn}.p"*
    [[ "$(uname)" == "Darwin" ]] && open "loc_eq_ellipse_${eqn}.png"

done

rm -f t77 tjp tmp.cpt tmpxy.grd ell.log
