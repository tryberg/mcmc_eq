#!/usr/bin/env bash
#
# pha2mcmc_eq_picks.sh
#
# Convert a hypoDD-style phase.dat file + hypoDD station file into the
# "picks.mcmc" format, matching pha2mcmc.m's station handling:
#
#   # event_num numP numS reftime(YYYYMMDDHHMMSS.sss)
#   STA index PHASE   X        Y        Z    traveltime  qual
#
# Usage:
#   ./pha2mcmc_eq_picks.sh phase.dat station.dat origin_lon origin_lat > picks.dat
#
# Inputs:
#   phase.dat   : hypoDD-style phase file. Event header lines:
#                   # YR MO DY HR MN SEC LAT LON DEPTH ...
#                 Pick lines:
#                   STA  traveltime  weight  PHASE
#
#   station.dat : hypoDD-style station file, one station per line:
#                   STA  LAT  LON  [ELEV]
#                 (ELEV assumed meters, positive up; optional, defaults to 0)
#
#   origin_lon/lat : local Cartesian origin (decimal degrees)
#
# Station filtering & indexing (matches pha2mcmc.m):
#   - Only stations that actually appear in phase.dat are kept; any station
#     in station.dat with zero picks is dropped entirely (never indexed).
#   - The kept stations are then sorted ALPHABETICALLY by code, and the
#     station INDEX = 0-based row position in that filtered+sorted list
#     (i.e. matches MATLAB's `intersect`, which sorts its output by default).
#     This index is NOT the line number in the original station.dat.
#
# Output fields per pick:
#   X, Y  = station position (km) relative to origin, via GMT azimuthal
#           equidistant projection (X=east, Y=north)
#   Z     = -elevation_km  (depth-positive-down convention)
#   qual  = weight bin, from continuous weight w:
#             0.75 < w <= 1.00 -> 0
#             0.50 < w <= 0.75 -> 1
#             0.25 < w <= 0.50 -> 2
#             0.00 <= w <= 0.25 -> 3
#             w < 0            -> 0 (warning printed to stderr)
#             w > 1            -> error
#
# Picks within each event are sorted: all P picks first (alphabetical by
# station), then all S picks (alphabetical by station) — matches
# sortrows(pha,[4,1]) in pha2mcmc.m.
#
# Requires: gmt (mapproject)

set -euo pipefail

if [[ $# -ne 4 ]]; then
    echo "Usage: $0 phase.dat station.dat origin_lon origin_lat" >&2
    exit 1
fi

phase_file="$1"
station_file="$2"
origin_lon="$3"
origin_lat="$4"

for f in "$phase_file" "$station_file"; do
    if [[ ! -f "$f" ]]; then
        echo "Error: file '$f' not found" >&2
        exit 1
    fi
done

if ! command -v gmt &> /dev/null; then
    echo "Error: gmt not found in PATH" >&2
    exit 1
fi

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

used_stations="$tmpdir/used_stations.txt"
filtered_sorted_stations="$tmpdir/filtered_sorted_stations.txt"
sta_lonlat="$tmpdir/sta_lonlat.txt"
sta_xy="$tmpdir/sta_xy.txt"
sta_lookup="$tmpdir/sta_lookup.txt"

# --- 1. Find every station code actually used in phase.dat (pick lines only) ---
awk '$1 != "#" { print $1 }' "$phase_file" | sort -u > "$used_stations"

# --- 2. Filter station.dat to only used stations, then sort alphabetically by code ---
#     (mirrors MATLAB's intersect(), which filters AND sorts by default)
awk 'NR==FNR { used[$1]=1; next } ($1 in used)' "$used_stations" "$station_file" \
    | sort -k1,1 > "$filtered_sorted_stations"

# Sanity check: warn about any used station missing from station.dat
awk 'NR==FNR { sta[$1]=1; next } { if (!($1 in sta)) print $1 }' \
    <(awk '{print $1}' "$filtered_sorted_stations") "$used_stations" > "$tmpdir/missing.txt"
if [[ -s "$tmpdir/missing.txt" ]]; then
    echo "Error: phase.dat references station(s) not found in station.dat:" >&2
    cat "$tmpdir/missing.txt" >&2
    exit 1
fi

# --- 3. Project filtered+sorted station lon/lat -> local km X,Y ---
awk '{ print $3, $2 }' "$filtered_sorted_stations" > "$sta_lonlat"   # LON LAT
gmt mapproject "$sta_lonlat" -Je"${origin_lon}"/"${origin_lat}"/1:1 -Fk -C -R"g" > "$sta_xy"

# --- 4. Build station lookup: CODE INDEX X Y Z ---
#     Index = 0-based row position in filtered_sorted_stations (NOT line number
#     in the original station.dat) — matches pha2mcmc.m's sii-1 after intersect.
paste "$filtered_sorted_stations" "$sta_xy" | awk '
{
    code = $1
    elev_m = (NF == 6) ? $4 : 0    # STA LAT LON [ELEV] X Y
    x = $(NF-1)
    y = $NF
    z = -elev_m / 1000.0           # km, depth-positive-down
    printf "%s %d %.4f %.4f %.4f\n", code, NR-1, x, y, z
}
' > "$sta_lookup"

# --- 5. Process phase.dat, event-block aware, emit new format ---
awk -v lookup_file="$sta_lookup" '
BEGIN {
    while ((getline line < lookup_file) > 0) {
        split(line, a, " ")
        idx[a[1]] = a[2]
        sx[a[1]]  = a[3]
        sy[a[1]]  = a[4]
        sz[a[1]]  = a[5]
    }
    event_num = -1
    have_event = 0
}

function qual_from_wgt(w,   q) {
    if (w <= 1.0 && w > 0.75)      q = 0
    else if (w <= 0.75 && w > 0.5) q = 1
    else if (w <= 0.5 && w > 0.25) q = 2
    else if (w <= 0.25 && w >= 0)  q = 3
    else if (w < 0) {
        printf "warning: neg pick wgt %s set to qual 0: flagged or error?\n", w > "/dev/stderr"
        q = 0
    } else {
        printf "error: bad weight %s\n", w > "/dev/stderr"
        exit 1
    }
    return q
}

function flush_event() {
    if (!have_event) return
    printf "# %d %d %d %s\n", event_num, nP, nS, reftime

    nPk = 0
    for (i = 1; i <= nlines; i++) {
        if (phase_arr[i] == "P") { nPk++; p_sta[nPk] = sta_arr[i]; p_tt[nPk] = tt_arr[i]; p_wgt[nPk] = wgt_arr[i] }
    }
    nSk = 0
    for (i = 1; i <= nlines; i++) {
        if (phase_arr[i] == "S") { nSk++; s_sta[nSk] = sta_arr[i]; s_tt[nSk] = tt_arr[i]; s_wgt[nSk] = wgt_arr[i] }
    }

    for (i = 2; i <= nPk; i++) {
        key_s = p_sta[i]; key_t = p_tt[i]; key_w = p_wgt[i]; j = i - 1
        while (j >= 1 && p_sta[j] > key_s) {
            p_sta[j+1] = p_sta[j]; p_tt[j+1] = p_tt[j]; p_wgt[j+1] = p_wgt[j]; j--
        }
        p_sta[j+1] = key_s; p_tt[j+1] = key_t; p_wgt[j+1] = key_w
    }
    for (i = 2; i <= nSk; i++) {
        key_s = s_sta[i]; key_t = s_tt[i]; key_w = s_wgt[i]; j = i - 1
        while (j >= 1 && s_sta[j] > key_s) {
            s_sta[j+1] = s_sta[j]; s_tt[j+1] = s_tt[j]; s_wgt[j+1] = s_wgt[j]; j--
        }
        s_sta[j+1] = key_s; s_tt[j+1] = key_t; s_wgt[j+1] = key_w
    }

    for (i = 1; i <= nPk; i++) {
        code = p_sta[i]
        q = qual_from_wgt(p_wgt[i])
        printf "%s %03d P %8.3f %8.3f %8.3f %8.3f %d\n", code, idx[code], sx[code], sy[code], sz[code], p_tt[i], q
    }
    for (i = 1; i <= nSk; i++) {
        code = s_sta[i]
        q = qual_from_wgt(s_wgt[i])
        printf "%s %03d S %8.3f %8.3f %8.3f %8.3f %d\n", code, idx[code], sx[code], sy[code], sz[code], s_tt[i], q
    }
}

$1 == "#" {
    flush_event()

    event_num++
    yr = $2; mo = $3; dy = $4; hr = $5; mn = $6; sec = $7
    reftime = sprintf("%04d%02d%02d%02d%02d%06.3f", yr, mo, dy, hr, mn, sec)

    nlines = 0
    nP = 0
    nS = 0
    have_event = 1
    next
}

{
    nlines++
    sta_arr[nlines]   = $1
    tt_arr[nlines]    = $2
    wgt_arr[nlines]   = $3
    phase_arr[nlines] = $4
    if ($4 == "P") nP++
    else if ($4 == "S") nS++
}

END {
    flush_event()
}
' "$phase_file"
