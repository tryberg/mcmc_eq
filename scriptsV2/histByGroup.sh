#!/bin/bash
# histByGroup.sh
#
# Reads a data file with whitespace-separated columns, e.g.:
#   RES P.   225499 0 0.115835 0.165560 0.379709
# and creates histograms of column 6, one panel per unique value in
# column 4, tiled onto pages of up to (rows x cols) panels each.
#
# Also reads a station lookup file mapping station number -> station
# name, e.g.:
#   0 0.143 0.349	AKRO 000 -7.747 -18.374 -0.179
# (column 1 = station number, column 4 = station name) and uses the
# station name in each panel's title instead of the raw number.
#
# Usage: ./histByGroup.sh datafile.txt stationfile.txt [outdir]

infile=$1
stafile=$2
outdir=${3:-.}
cols=4
rows=2
maxPerPage=$(( cols * rows ))
panelW=6
panelH=5
panelMargin=0.3

synStaCorFile=synStaCors.dat #JP add

if [ -z "$infile" ] || [ -z "$stafile" ]; then
    echo "Usage: $0 datafile.txt stationfile.txt [outdir]"
    exit 1
fi

if [ ! -f "$infile" ]; then
    echo "Error: file not found: $infile"
    exit 1
fi

if [ ! -f "$stafile" ]; then
    echo "Error: file not found: $stafile"
    exit 1
fi

mkdir -p "$outdir"

# ---- Station-number -> station-name lookup (Bash 3.2 compatible: no
# associative arrays -- looks up directly from the file each time) ----
getStaName() {
    local num="$1"
    local result
    result=$(awk -v n="$num" '$1 == n {print $4; exit}' "$stafile")
    if [ -z "$result" ]; then
        echo "$num"
    else
        echo "$result"
    fi
}

# Get unique values from column 4, into an array
groups=($(awk '{print $4}' "$infile" | sort -n -u))
nGroups=${#groups[@]}

nPages=$(( (nGroups + maxPerPage - 1) / maxPerPage ))

echo "Found $nGroups groups -> $nPages page(s)"

binWidth=0.002 #synthetics
#binWidth=0.05

# ---- Compute a shared X range across ALL groups (columns 6 AND 7) ----
read xmin xmax <<< $(awk '{print $6"\n"$7}' "$infile" | sort -n | awk '
    NR==1 {min=$1}
    {max=$1}
    END {print min, max}
')

# Pad the x range slightly and snap to the bin width so bins align nicely
xmax=`awk '{for(i=6;i<=NF;i++){v=($i<0)?-$i:$i; if(v>m)m=v}} END{print m}' $infile | gmt info -C -I0.25 -o1`
#xmax=$(awk -v a="$xmax" -v b="1" 'BEGIN{print (a<b)?a:b}') # use reasonable max if less than 1, otherwise 1
xmin=-${xmax}
echo "Shared X range: $xmin to $xmax"

# ---- Compute a shared Y range (max bin count) across ALL groups and BOTH columns (6 and 7) ----
# Done with pure awk (no gmt call) to avoid invoking the histogram module
# outside of a gmt begin/end session, which can trigger module-resolution
# errors when called repeatedly in a loop.
ymax=0
for g in "${groups[@]}"; do
    tmpfile6="${outdir}/tmp_yscan6_${g}.txt"
    tmpfile7="${outdir}/tmp_yscan7_${g}.txt"
    awk -v grp="$g" '$4 == grp {print $6}' "$infile" > "$tmpfile6"
    awk -v grp="$g" '$4 == grp {print $7}' "$infile" > "$tmpfile7"

    for tmpfile in "$tmpfile6" "$tmpfile7"; do
        n=$(wc -l < "$tmpfile")
        if [ "$n" -gt 0 ]; then
            thisMax=$(awk -v bw="$binWidth" -v xmin="$xmin" '
                { bin = int(($1 - xmin) / bw); count[bin]++ }
                END {
                    m = 0
                    for (b in count) if (count[b] > m) m = count[b]
                    print m
                }
            ' "$tmpfile")
            ymax=$(echo "$ymax $thisMax" | awk '{print ($1>$2)?$1:$2}')
        fi
    done
    rm -f "$tmpfile6" "$tmpfile7"
done

# pad y max a bit so tallest bar isn't flush with the top of the panel
ymax=$(echo "$ymax" | awk '{print $1*1.1}')

echo "Shared Y range: 0 to $ymax"

sharedR="${xmin}/${xmax}/0/${ymax}"

# Legend spec file (same for every panel) -- small color-swatch legend
# mapping column 6 (P-wave) and column 7 (S-wave) to their fill colors.
legendFile="${outdir}/tmp_legend.txt"
cat << EOF > "$legendFile"
S 0.1c s 0.3c grey 0.0p 0.4c P-wave
S 0.1c s 0.3c steelblue 0.0p 0.4c S-wave
EOF

for (( p=0; p<nPages; p++ )); do

    pageOut="${outdir}/staCorHist_pg_$((p+1))"
    startIdx=$(( p * maxPerPage ))
    endIdx=$(( startIdx + maxPerPage - 1 ))
    if [ $endIdx -ge $nGroups ]; then endIdx=$(( nGroups - 1 )); fi

    echo "Page $((p+1)): groups ${groups[$startIdx]}..${groups[$endIdx]}"

    gmt begin "$pageOut" png
        gmt set MAP_TITLE_OFFSET 1.5p
        gmt subplot begin ${rows}x${cols} -Fs${panelW}c/${panelH}c -M${panelMargin}c

        panel=0
        for (( gi=startIdx; gi<=endIdx; gi++ )); do
            g=${groups[$gi]}
            tmpfile6="${outdir}/tmp_group6_${g}.txt"
            tmpfile7="${outdir}/tmp_group7_${g}.txt"
            statsFile="${outdir}/tmp_stats_${g}.txt"

            awk -v grp="$g" '$4 == grp {print $6}' "$infile" > "$tmpfile6"
            awk -v grp="$g" '$4 == grp {print $7}' "$infile" > "$tmpfile7"
            n=$(wc -l < "$tmpfile6")

            # Mean and standard deviation for both columns
            read mean6 std6 <<< $(awk '
                { s+=$1; ss+=$1*$1; n++ }
                END {
                    if (n>0) {
                        m = s/n
                        v = ss/n - m*m
                        if (v<0) v=0
                        printf "%.4f %.4f", m, sqrt(v)
                    } else {
                        print "0 0"
                    }
                }
            ' "$tmpfile6")

            read mean7 std7 <<< $(awk '
                { s+=$1; ss+=$1*$1; n++ }
                END {
                    if (n>0) {
                        m = s/n
                        v = ss/n - m*m
                        if (v<0) v=0
                        printf "%.4f %.4f", m, sqrt(v)
                    } else {
                        print "0 0"
                    }
                }
            ' "$tmpfile7")

            row=$(( panel / cols ))
            col=$(( panel % cols ))

            gmt subplot set $row,$col

            staLabel=$(getStaName "$g")

            if [ "$n" -eq 0 ]; then
                echo "  no data for group $g, leaving panel blank"
            else
		#JP addition to plot syn sta cor if present:
		if [ -f "$synStaCorFile" ]; then
		  psyncsc=`grep $staLabel $synStaCorFile | awk '{print $5}'`	
		  echo "$psyncsc 0" > ptmp
		  echo "$psyncsc $ymax" >> ptmp
		  gmt psxy ptmp -JX?/? -R$sharedR -W0.5p,grey
		  ssyncsc=`grep $staLabel $synStaCorFile | awk '{print $6}'`	
		  echo "$ssyncsc 0" > stmp
		  echo "$ssyncsc $ymax" >> stmp
		  gmt psxy stmp -W0.5p,steelblue
		fi

                # Column 7 plotted first (background layer), fully opaque
                gmt histogram "$tmpfile7" -JX?/? -R${sharedR} -T${binWidth} -B+t"Station ${staLabel}" -Bxaf -Byaf -Gsteelblue 

                # Column 6 plotted on top (foreground layer), semi-transparent
                # grey so both distributions remain visible where they overlap
                gmt histogram "$tmpfile6" -R${sharedR} -T${binWidth} -Ggrey -t30

                # Small legend, upper-left corner of this panel
                gmt legend "$legendFile" -DjTL+w2.2c+o0.1c/0.1c -F+p0.5p+gwhite

                # Mean/std text box, upper-right corner of this panel
                # (positioned via actual data coordinates near xmax/ymax,
                # since this panel's -R is sharedR)
                xpos=$(echo "$xmax" | awk '{print $1*0.98}')
                ypos1=$(echo "$ymax" | awk '{print $1*0.95}')
                ypos2=$(echo "$ymax" | awk '{print $1*0.86}')
                cat << EOF > "$statsFile"
$xpos $ypos1 P: ${mean6} +/- ${std6}
$xpos $ypos2 S: ${mean7} +/- ${std7}
EOF
                gmt text "$statsFile" -R${sharedR} -F+f7p,Helvetica,black+jTR -Gwhite -W0.5p
            fi

            rm -f "$tmpfile6" "$tmpfile7" "$statsFile"
            panel=$(( panel + 1 ))
        done

        gmt subplot end

        # Small page indicator, bottom-right corner of the physical page
        pageW=$(echo "$cols" | awk -v s="$panelW" -v m="$panelMargin" '{print $1*s + ($1-1)*m}')
        pageH=$(echo "$rows" | awk -v s="$panelH" -v m="$panelMargin" '{print $1*s + ($1-1)*m}')
#        echo "$pageW 0 Page $((p+1))/${nPages}" | gmt text -R0/${pageW}/0/${pageH} -Jx1c -F+f8p,Helvetica,black+jBR -D-0.1c/0.1c -N
    gmt end
    ls $outdir/$pageOut.png
    [[ "$(uname)" == "Darwin" ]] && open $pageOut.png

done

rm -f "$legendFile"
