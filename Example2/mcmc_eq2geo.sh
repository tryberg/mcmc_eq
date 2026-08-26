#!/usr/bin/env bash
#
# mcmc2geo.sh
#
# Add geographic coordinates (lon, lat) to a "quakes_mcmc.dat" file (EVID +
# local Cartesian X,Y,Z + other columns), using the same azimuthal-
# equidistant projection convention as pha2mcmc_eq_picks.sh:
#
#   forward (lon,lat -> X,Y km):
#     gmt mapproject lonlat.txt -Je<origin_lon>/<origin_lat>/1:1 -Fk -C -Rg
#
#   inverse (X,Y km -> lon,lat), used here:
#     gmt mapproject xy.txt -Je<origin_lon>/<origin_lat>/1:1 -Fk -C -Rg -I
#
# Input file format (quakes_mcmc.dat), one comment header line then rows:
#   # {'EVID','X','Y','Z','OT','dOT','ex','ey','ez'} (%8.3f %8.3f %8.3f %015.3f %7.3f %f %f %f\n)
#   EVID     X        Y        Z        OT                  dOT      ex       ey       ez
#
# Output format: same columns as the input, in the same order, but
# reformatted to fixed-width (values unchanged, just consistently
# aligned) with two new columns -- Lon and Lat -- appended at the end.
#
# Usage:
#   ./mcmc2geo.sh quakes_mcmc.dat origin_lon origin_lat > quakes_geo.dat
#
# Requires: gmt (mapproject)

set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 quakes_mcmc.dat origin_lon origin_lat" >&2
    exit 1
fi

infile="$1"
origin_lon="$2"
origin_lat="$3"

if [[ ! -f "$infile" ]]; then
    echo "Error: file '$infile' not found" >&2
    exit 1
fi

if ! command -v gmt &> /dev/null; then
    echo "Error: gmt not found in PATH" >&2
    exit 1
fi

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

xy_file="$tmpdir/xy.txt"
lonlat_file="$tmpdir/lonlat.txt"
data_file="$tmpdir/data.txt"

# --- 1. Extract X,Y (for projection) and keep the full original data
#        lines untouched, skipping the comment header line ---
awk '$1 != "#" { print $2, $3 }' "$infile" > "$xy_file"
grep -v '^#' "$infile" > "$data_file"

# --- 2. Inverse-project local km X,Y -> lon,lat ---
#     -I requests the inverse transform (Cartesian -> geographic)
gmt mapproject "$xy_file" -Je"${origin_lon}"/"${origin_lat}"/1:1 -Fk -C -R"g" -I > "$lonlat_file"

# --- 3. Output: original header + new columns noted, then each row
#        reformatted with fixed-width columns (original values unchanged,
#        just consistently aligned) with lon/lat appended ---
grep '^#' "$infile"
echo "# appended columns: Lon Lat (%9.4f %9.4f)"
paste "$data_file" "$lonlat_file" | awk '{
    evid=$1; x=$2; y=$3; z=$4; ot=$5; dot=$6; ex=$7; ey=$8; ez=$9
    lon=$10; lat=$11
    printf "%03d %8.3f %8.3f %8.3f %015.3f %7.3f %8.6f %8.6f %8.6f %9.4f %9.4f\n", \
        evid, x, y, z, ot, dot, ex, ey, ez, lon, lat
}'
