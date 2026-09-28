#!/usr/bin/env bash
#
# run_gmatrix.sh
#
# One-shot driver for the gmatrix covariance / resolution analysis of an
# mcmc_eq output chain. It:
#   1. builds the `gmatrix` binary if it is missing (via the src Makefile),
#   2. runs gmatrix in GLOBAL 1c mode to emit the sparse sensitivity matrix,
#   3. runs gmatrix_cov_res.py to compute Cm and R and write ALL figures,
#   4. (optional) runs the loc-mode + GMT per-event error-ellipse plot.
#
# Usage:
#   run_gmatrix.sh <chain_file> [pick_file] [config] [outdir]
#
#   chain_file : an mcmc_eq output chain (e.g. rjx-009_4621867.out)   [required]
#   pick_file  : the pick file used for the run   (default: picks.mcmc)
#   config     : the config file                  (default: config_eqx.dat)
#   outdir     : where to write gmatrix.out + figures (default: dir of chain)
#
# Env overrides:
#   TAG=bat|mod        model to extract from the chain        (default bat)
#   FDSTEP=0.01        finite-difference step (km) for x,y,z
#   CWEIGHT=1e6        zero-mean station-correction constraint weight
#   DVP=0.02 DVPVS=0.01 velocity finite-difference steps
#   RCOND=1e-8         SVD truncation ratio for cov_res.py
#   DPI=130            figure resolution
#   SAVE_NPZ=1         also dump Cm/R/spectrum to a .npz
#   DO_ELLIPSE=1       also produce the loc-mode GMT error-ellipse figure
#                      (needs GMT + resmcnx.dat/tmpx in the run dir)
#   PYTHON=python3     python interpreter to use
#
# J. Pesicek / Kiro, 2026

set -euo pipefail

# --- locate this script, the src dir, and the tools ------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC_DIR="$(cd "$SCRIPT_DIR/../src" 2>/dev/null && pwd || true)"
if [[ -z "${SRC_DIR:-}" || ! -d "$SRC_DIR" ]]; then
    echo "error: cannot locate src/ relative to $SCRIPT_DIR" >&2
    exit 1
fi
GMATRIX_BIN="$SRC_DIR/gmatrix"
COV_RES_PY="$SRC_DIR/gmatrix_cov_res.py"
ELLIPSE_SH="$SCRIPT_DIR/disp_eq_z_ellipse.sh"
PYTHON="${PYTHON:-python3}"

# --- arguments -------------------------------------------------------------
if [[ $# -lt 1 ]]; then
    grep '^#' "$0" | sed -n '2,40p' | sed 's/^# \{0,1\}//'
    exit 1
fi
CHAIN="$1"
PICKS="${2:-picks.mcmc}"
CONFIG="${3:-config_eqx.dat}"
OUTDIR="${4:-$(dirname "$CHAIN")}"

TAG="${TAG:-bat}"
FDSTEP="${FDSTEP:-0.01}"
CWEIGHT="${CWEIGHT:-1e6}"
DVP="${DVP:-0.02}"
DVPVS="${DVPVS:-0.01}"
RCOND="${RCOND:-1e-8}"
DPI="${DPI:-130}"

# --- sanity checks ---------------------------------------------------------
for f in "$CHAIN" "$PICKS" "$CONFIG"; do
    if [[ ! -f "$f" ]]; then
        echo "error: required file not found: $f" >&2
        exit 1
    fi
done
mkdir -p "$OUTDIR"

# --- 1. build gmatrix if needed --------------------------------------------
if [[ ! -x "$GMATRIX_BIN" ]]; then
    echo ">> gmatrix binary not found; building via make ..."
    ( cd "$SRC_DIR" && make gmatrix )
fi

# --- 2. run gmatrix in GLOBAL mode -----------------------------------------
# default output name embeds the chain id (gmatrix picks it automatically when
# the 6th arg is omitted). We pass an explicit path so it lands in OUTDIR.
base="$(basename "$CHAIN")"
cid="${base#rjx-}"; cid="${cid%%_*}"; cid="${cid%%.*}"
GOUT="$OUTDIR/gmatrix_${cid}.out"

echo ">> running gmatrix (global) on $CHAIN"
echo "   config=$CONFIG picks=$PICKS tag=$TAG -> $GOUT"
"$GMATRIX_BIN" "$CONFIG" "$CHAIN" "$PICKS" global "$TAG" "$GOUT" \
    "$FDSTEP" "$CWEIGHT" "$DVP" "$DVPVS"

# --- 3. covariance + resolution figures ------------------------------------
echo ">> computing Cm & R and writing figures with gmatrix_cov_res.py"
covres_args=( "$GOUT" --outdir "$OUTDIR" --rcond "$RCOND" --dpi "$DPI" )
if [[ "${SAVE_NPZ:-0}" == "1" ]]; then
    covres_args+=( --save-npz )
fi
"$PYTHON" "$COV_RES_PY" "${covres_args[@]}"

# --- 4. optional: loc-mode GMT per-event error-ellipse figure --------------
if [[ "${DO_ELLIPSE:-0}" == "1" ]]; then
    if [[ -x "$ELLIPSE_SH" || -f "$ELLIPSE_SH" ]]; then
        echo ">> producing loc-mode GMT error-ellipse figure"
        # the ellipse script auto-builds its loc-mode gmatrix output if given the
        # chain + picks via env; run it from OUTDIR so it finds config/resmcnx/tmpx.
        ( cd "$OUTDIR" && \
          GMATRIX_CHAIN="$(cd "$(dirname "$CHAIN")" && pwd)/$(basename "$CHAIN")" \
          GMATRIX_PICKS="$(cd "$(dirname "$PICKS")" && pwd)/$(basename "$PICKS")" \
          GMATRIX_TAG="$TAG" \
          bash "$ELLIPSE_SH" ) || \
          echo "   (ellipse step failed -- needs GMT + resmcnx.dat/tmpx in $OUTDIR)"
    else
        echo "   (ellipse script not found at $ELLIPSE_SH -- skipping)"
    fi
fi

echo ">> done. outputs in $OUTDIR"
ls -1 "$OUTDIR"/gmatrix_"${cid}"*.png 2>/dev/null || true
