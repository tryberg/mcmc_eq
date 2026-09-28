#!/usr/bin/env python3
"""
gmatrix_cov_res.py -- covariance & resolution from a global gmatrix.out

Reads the GLOBAL-mode output of the `gmatrix` tool (records: header `#`,
`D` data-row meta, `G` sparse triplets, `C` constraint-row meta), forms the
weighted normal matrix G^T W G, inverts it with a *truncated SVD*, and produces:

    posterior model covariance     Cm = (G^T W G)^-1        (Tarantola)
    model resolution matrix        R  = Cm Gd^T Wd Gd        (Menke)

Then it writes a set of figures visualizing Cm and R (full matrices, the
model-parameter correlation matrix, focused sub-views on the velocity and
station-correction blocks, the resolution diagonal, and the singular-value
spectrum used for truncation).

Usage:
    gmatrix_cov_res.py gmatrix.out [--outdir DIR] [--rcond RCOND]
                                   [--dpi DPI] [--no-show]

The column layout (loc / vp / vpvs / pres / sres blocks) is read automatically
from the `# colmap` header line, so block sub-views are labelled correctly.

Requires: numpy, scipy, matplotlib.
"""

import argparse
import os
import re
import sys

import numpy as np
from scipy.sparse import coo_matrix

import matplotlib
matplotlib.use("Agg")            # default: headless-safe; overridden by --show handling below
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------
class GMatrix:
    """Container for the parsed global gmatrix.out."""

    def __init__(self):
        self.ncols = None
        self.ndata_rows = None
        self.colmap = {}          # name -> (lo, hi) inclusive column range
        self.rows = []            # triplet row indices
        self.cols = []            # triplet col indices
        self.vals = []            # triplet values
        self.sigma = {}           # data row -> pick noise (sigma)
        self.is_data = {}         # row -> True (data) / False (constraint)
        self.dims = {}            # ne, dim, nos, ... from header
        self.cweight = None       # station-correction constraint weight (from header)
        self.loc_per_event = 4    # location params per event (4=x,y,z,ot; 1=ot only if loc fixed)
        self.fixed = {}           # e.g. {"vel":1,"loc":0,"res":0} from the header, if present


def _parse_colmap(line):
    """Parse `# colmap loc 0..879 ... vp 880..894 ...` into {name:(lo,hi)}."""
    out = {}
    for name, lo, hi in re.findall(r"([a-zA-Z]+)\s+(\d+)\.\.(\d+)", line):
        out[name] = (int(lo), int(hi))
    return out


def parse_gmatrix(path):
    g = GMatrix()
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            p = line.split()
            tag = p[0]
            if tag == "#":
                if "ncols" in p:
                    g.ncols = int(p[p.index("ncols") + 1])
                if "ndata_rows" in p:
                    g.ndata_rows = int(p[p.index("ndata_rows") + 1])
                if "colmap" in line:
                    g.colmap = _parse_colmap(line)
                    # loc stride: "loc A..B (N/event: ...)"
                    mstride = re.search(r"\((\d+)/event", line)
                    if mstride:
                        g.loc_per_event = int(mstride.group(1))
                # fixed-class flags: "# fixed  vp 0  vpvs 1  loc 0  res 0 ..."
                if p[1:2] == ["fixed"]:
                    for key in ("vp", "vpvs", "vel", "loc", "res"):
                        if key in p:
                            try:
                                g.fixed[key] = int(p[p.index(key) + 1])
                            except (ValueError, IndexError):
                                pass
                # dims header: "# dims ne 220  dim 15  nos 130 ..."
                if p[1:2] == ["dims"]:
                    for key in ("ne", "dim", "nos", "eikonal", "TRIA"):
                        if key in p:
                            g.dims[key] = int(p[p.index(key) + 1])
                # fd header carries the constraint weight: "... cweight 1.000e+06"
                if "cweight" in p:
                    try:
                        g.cweight = float(p[p.index("cweight") + 1])
                    except (ValueError, IndexError):
                        pass
            elif tag == "D":
                # D <row> <P|S> <evt> <st> <cl> <sigma> <resid>
                r = int(p[1])
                g.sigma[r] = float(p[6])
                g.is_data[r] = True
            elif tag == "C":
                # C <row> <meanP|meanS> <weight>
                g.is_data[int(p[1])] = False
            elif tag == "G":
                # G <row> <col> <value>
                g.rows.append(int(p[1]))
                g.cols.append(int(p[2]))
                g.vals.append(float(p[3]))
    if g.ncols is None:
        # fall back to inferring from triplets
        g.ncols = max(g.cols) + 1
    return g


# ---------------------------------------------------------------------------
# Linear algebra: Cm and R
# ---------------------------------------------------------------------------
def build_sparse_G(g):
    nrows = max(g.rows) + 1
    G = coo_matrix((g.vals, (g.rows, g.cols)), shape=(nrows, g.ncols)).tocsr()
    return G, nrows


def compute_cov_res(g, rcond=1e-8, verbose=True):
    """Return (Cm, R, S, ntrunc, data_rows) using a truncated SVD of GtWG.

    Row weighting: sqrt(1/sigma^2)=1/sigma for data rows; constraint rows keep
    weight 1 (their strength `cweight` is already baked into the G entries).
    """
    G, nrows = build_sparse_G(g)

    # row weight vector w such that (sqrt(W) G) has the correct scaling
    w = np.ones(nrows)
    for r, s in g.sigma.items():
        w[r] = 1.0 / s                          # sqrt(1/sigma^2)
    Gw = G.multiply(w[:, None]).tocsr()

    GtG = (Gw.T @ Gw).toarray()                 # (ncols x ncols) dense normal matrix

    if verbose:
        print(f"  normal matrix GtWG: {GtG.shape}, "
              f"symmetric residual {np.abs(GtG - GtG.T).max():.2e}")

    # symmetric eigendecomposition (GtWG is SPD up to null space) is more stable
    # and cheaper than a full SVD for a symmetric matrix.
    evals, evecs = np.linalg.eigh(GtG)
    # eigh returns ascending; flip to descending for a singular-value-like spectrum
    order = np.argsort(evals)[::-1]
    S = evals[order]
    V = evecs[:, order]

    # --- truncation reference -------------------------------------------------
    # The zero-mean station-correction constraint rows inject `cweight` entries,
    # so they contribute eigenvalues ~ cweight^2 * (#columns constrained) that are
    # orders of magnitude larger than any data-resolved direction. Referencing the
    # truncation tolerance to the global max eigenvalue would then wrongly discard
    # every real data direction. We therefore reference `rcond` to the largest
    # *data* eigenvalue (i.e. exclude the near-hard constraint modes), which is the
    # physically meaningful scale for "well- vs weakly-resolved".
    constraint_floor = None
    if g.cweight and g.cweight > 0:
        # constraint eigenvalues are O(cweight^2); use a safe fraction as the cut
        constraint_floor = (g.cweight ** 2)
    if constraint_floor is not None:
        data_mask = S < constraint_floor * 0.5
        smax = S[data_mask].max() if np.any(data_mask) else S.max()
        n_constraint = int((~data_mask).sum())
    else:
        smax = S.max()
        n_constraint = 0

    tol = smax * rcond
    keep = S > tol                              # keep constraint modes too (S huge)
    ntrunc = int((~keep).sum())
    Sinv = np.zeros_like(S)
    Sinv[keep] = 1.0 / S[keep]                  # guarded: no division on dropped modes

    # (GtWG)^-1 via truncated pseudo-inverse:  V diag(Sinv) V^T
    # NOTE: numpy 2.x can raise spurious "divide by zero / overflow in matmul"
    # FP-flag warnings from the BLAS path even when the result is fully finite
    # (subnormal intermediates). Sinv is guarded above, so we assert finiteness.
    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
        GtG_inv = (V * Sinv) @ V.T
    if not np.all(np.isfinite(GtG_inv)):
        raise FloatingPointError(
            "non-finite covariance after truncation; try a larger --rcond")
    Cm = GtG_inv

    if verbose:
        print(f"  eigen-spectrum: max {S.max():.3e}, min {S.min():.3e}")
        if n_constraint:
            print(f"  constraint modes (~cweight^2): {n_constraint}; "
                  f"largest data eigenvalue {smax:.3e}")
        print(f"  truncated {ntrunc}/{len(S)} modes below tol {tol:.3e} "
              f"(rcond={rcond:g})")

    # resolution: use DATA rows only for the mapping term
    data_rows = sorted(r for r in range(nrows) if g.is_data.get(r, False))
    Gd = G[data_rows, :]
    Wd = np.array([1.0 / g.sigma[r] ** 2 for r in data_rows])
    # Gg = (GtWG)^-1 Gd^T Wd ;  R = Gg Gd
    GdT_Wd = Gd.T.multiply(Wd).tocsr()          # ncols x ndata
    Gg = GtG_inv @ GdT_Wd                       # ncols x ndata (dense)
    R = Gg @ Gd                                 # ncols x ncols (dense)
    R = np.asarray(R)

    return Cm, R, S, ntrunc, data_rows, smax


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def _diverging(ax, M, title, vmax=None, cmap="RdBu_r"):
    if vmax is None:
        vmax = np.nanmax(np.abs(M))
        if vmax == 0 or not np.isfinite(vmax):
            vmax = 1.0
    norm = TwoSlopeNorm(vcenter=0.0, vmin=-vmax, vmax=vmax)
    im = ax.imshow(M, cmap=cmap, norm=norm, aspect="auto", interpolation="nearest")
    ax.set_title(title, fontsize=10)
    return im


def _block_boundaries(colmap):
    """Return list of (name, lo, hi) sorted by lo for drawing block guides."""
    blocks = [(name, lo, hi) for name, (lo, hi) in colmap.items()]
    blocks.sort(key=lambda t: t[1])
    return blocks


def _draw_block_guides(ax, blocks, n):
    for _, lo, hi in blocks:
        for edge in (lo, hi + 1):
            if 0 < edge < n:
                ax.axhline(edge - 0.5, color="k", lw=0.4, alpha=0.3)
                ax.axvline(edge - 0.5, color="k", lw=0.4, alpha=0.3)


def correlation_from_cov(Cm):
    d = np.sqrt(np.clip(np.diag(Cm), 0, None))
    with np.errstate(divide="ignore", invalid="ignore"):
        Dinv = np.where(d > 0, 1.0 / d, 0.0)
    Corr = (Cm * Dinv[:, None]) * Dinv[None, :]
    np.fill_diagonal(Corr, 1.0)
    return Corr


def fig_full_matrices(Cm, R, colmap, outpath, label=""):
    blocks = _block_boundaries(colmap)
    n = Cm.shape[0]
    # width kept <= ~1950px at the default 130 dpi so the 3-panel figure stays under
    # common multi-image viewer limits (2000px/side).
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.2))

    # covariance (log-abs to cope with huge dynamic range across blocks)
    Cabs = np.abs(Cm)
    floor = Cabs[Cabs > 0].min() if np.any(Cabs > 0) else 1e-30
    im0 = axes[0].imshow(np.log10(Cabs + floor), cmap="viridis",
                         aspect="auto", interpolation="nearest")
    axes[0].set_title("Covariance  log10|Cm|", fontsize=10)
    fig.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04)

    # correlation matrix (bounded [-1,1], most interpretable full view)
    Corr = correlation_from_cov(Cm)
    im1 = _diverging(axes[1], Corr, "Correlation matrix", vmax=1.0)
    fig.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)

    # resolution
    im2 = _diverging(axes[2], R, "Resolution R", vmax=1.0)
    fig.colorbar(im2, ax=axes[2], fraction=0.046, pad=0.04)

    for ax in axes:
        _draw_block_guides(ax, blocks, n)
        ax.set_xlabel("model parameter (column)")
    axes[0].set_ylabel("model parameter (row)")

    # annotate block names at the centre of each block along the x-axis via ticks
    tick_pos = [(lo + hi) / 2 for _, lo, hi in blocks]
    tick_lbl = [name for name, _, _ in blocks]
    for ax in axes:
        ax.set_xticks(tick_pos)
        ax.set_xticklabels(tick_lbl, fontsize=7)

    fig.suptitle(_titled("Global model covariance & resolution", label), fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(outpath, dpi=fig.dpi)
    plt.close(fig)
    return outpath


def fig_block_view(Cm, R, colmap, block_names, title, outpath):
    """Focused sub-view: correlation + resolution restricted to given blocks."""
    idx = []
    labels_at = []
    for name in block_names:
        if name not in colmap:
            continue
        lo, hi = colmap[name]
        start = len(idx)
        idx.extend(range(lo, hi + 1))
        labels_at.append((name, start, len(idx)))
    if not idx:
        return None
    idx = np.array(idx)
    sub_C = Cm[np.ix_(idx, idx)]
    sub_R = R[np.ix_(idx, idx)]
    Corr = correlation_from_cov(sub_C)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.4))
    im0 = _diverging(axes[0], Corr, "Correlation", vmax=1.0)
    fig.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04)
    im1 = _diverging(axes[1], sub_R, "Resolution R", vmax=1.0)
    fig.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)

    n = len(idx)
    for ax in axes:
        for _, s, e in labels_at:
            for edge in (s, e):
                if 0 < edge < n:
                    ax.axhline(edge - 0.5, color="k", lw=0.5, alpha=0.4)
                    ax.axvline(edge - 0.5, color="k", lw=0.5, alpha=0.4)
        # place block labels at tick centres
        ticks = [(s + e) / 2 - 0.5 for _, s, e in labels_at]
        names = [nm for nm, _, _ in labels_at]
        ax.set_xticks(ticks)
        ax.set_xticklabels(names, fontsize=8)
        ax.set_yticks(ticks)
        ax.set_yticklabels(names, fontsize=8)

    fig.suptitle(title, fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(outpath, dpi=fig.dpi)
    plt.close(fig)
    return outpath


def _depth_columns(colmap, per_event):
    """Return the model-column indices of every event's hypocentre depth (z).

    In mcmc_eq the location block is stored interleaved per event as
    [x y z ot] (per_event == 4) or [ot] only (per_event == 1, locations fixed).
    Depth is component 2 of each length-`per_event` group. If locations were held
    fixed (per_event < 3) there is no depth degree of freedom and we return [].
    """
    if "loc" not in colmap or per_event < 3:
        return []
    lo, hi = colmap["loc"]
    return [c for c in range(lo, hi + 1) if (c - lo) % per_event == 2]


def fig_depth_cross(Cm, colmap, col_blocks, title, outpath, per_event=4):
    """Cross-correlation between every event's hypocentre depth (z) and another
    block (velocity nodes, or station corrections).

    Rows of the heatmap are events (their depth column), columns are the target
    block parameters (e.g. vp nodes, or pres/sres station corrections). This is the
    depth analogue of `fig_cross_block` and exposes the classic depth <-> velocity
    and depth <-> station-correction trade-offs. Left panel = the full
    cross-correlation heatmap; right panel = per-event max |correlation| against
    each target block, a compact summary of how strongly each event's depth trades
    off with that block.
    """
    ridx = _depth_columns(colmap, per_event)
    if not ridx:
        return None
    cols = [b for b in col_blocks if b in colmap]
    if not cols:
        return None

    # build column index list and per-block spans for labelling
    cidx, spans = [], []
    for b in cols:
        lo, hi = colmap[b]
        start = len(cidx)
        cidx.extend(range(lo, hi + 1))
        spans.append((b, start, len(cidx)))

    ridx = np.array(ridx)
    cidx_arr = np.array(cidx)
    ne = len(ridx)

    d = np.sqrt(np.clip(np.diag(Cm), 0, None))
    with np.errstate(divide="ignore", invalid="ignore"):
        dinv = np.where(d > 0, 1.0 / d, 0.0)
    sub = Cm[np.ix_(ridx, cidx_arr)] * dinv[ridx][:, None] * dinv[cidx_arr][None, :]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.4),
                             gridspec_kw={"width_ratios": [2.4, 1]})

    # ---- left: cross-correlation heatmap (events x target columns) ----
    vmax = 1.0
    im = axes[0].imshow(sub, cmap="RdBu_r",
                        norm=TwoSlopeNorm(vcenter=0.0, vmin=-vmax, vmax=vmax),
                        aspect="auto", interpolation="nearest")
    fig.colorbar(im, ax=axes[0], fraction=0.046, pad=0.04, label="correlation")
    axes[0].set_ylabel("event (hypocentre depth, z)")
    axes[0].set_xlabel("model parameter column")
    axes[0].set_title("cross-correlation", fontsize=10)
    # column block guides + labels (placed just inside the top edge so they
    # never collide with the panel title)
    for b, s, e in spans:
        if s > 0:
            axes[0].axvline(s - 0.5, color="k", lw=0.8, alpha=0.5)
        axes[0].text((s + e) / 2 - 0.5, 0.01 * ne, b, ha="center", va="top",
                     fontsize=9)

    # ---- right: per-event max |corr| against each target block ----
    y = np.arange(ne)
    for b, s, e in spans:
        maxc = np.abs(sub[:, s:e]).max(axis=1) if e > s else np.zeros(ne)
        axes[1].plot(maxc, y, lw=0.6, label=b)
    axes[1].set_ylim(ne - 0.5, -0.5)          # match imshow row order (0 at top)
    axes[1].set_ylabel("event")
    axes[1].set_xlabel("max |correlation|")
    axes[1].set_xlim(0, 1)
    axes[1].axvline(0.3, color="0.5", ls="--", lw=0.8)
    axes[1].legend(fontsize=8, loc="lower right")
    axes[1].set_title("peak coupling per event", fontsize=10)

    fig.suptitle(title, fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(outpath, dpi=fig.dpi)
    plt.close(fig)
    return outpath


def fig_cross_block(Cm, colmap, row_block, col_blocks, title, outpath):
    """Off-diagonal correlation between one block (rows) and others (columns).

    Designed to expose the velocity <-> station-correction trade-off: rows are the
    vp columns, columns are pres/sres. Left panel = the cross-correlation heatmap
    (rows = vp node, columns = station corrections, split by phase); right panel =
    per-vp-node max |correlation| against P and S corrections, a compact summary of
    how strongly each velocity node trades off with the statics.
    """
    if row_block not in colmap:
        return None
    cols = [b for b in col_blocks if b in colmap]
    if not cols:
        return None
    rlo, rhi = colmap[row_block]
    ridx = list(range(rlo, rhi + 1))

    # build column index list and per-block spans for labelling
    cidx, spans = [], []
    for b in cols:
        lo, hi = colmap[b]
        start = len(cidx)
        cidx.extend(range(lo, hi + 1))
        spans.append((b, start, len(cidx)))

    d = np.sqrt(np.clip(np.diag(Cm), 0, None))
    with np.errstate(divide="ignore", invalid="ignore"):
        dinv = np.where(d > 0, 1.0 / d, 0.0)
    sub = Cm[np.ix_(ridx, cidx)] * dinv[ridx][:, None] * dinv[cidx][None, :]

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8),
                             gridspec_kw={"width_ratios": [2.4, 1]})

    # ---- left: cross-correlation heatmap ----
    vmax = 1.0
    im = axes[0].imshow(sub, cmap="RdBu_r",
                        norm=TwoSlopeNorm(vcenter=0.0, vmin=-vmax, vmax=vmax),
                        aspect="auto", interpolation="nearest")
    fig.colorbar(im, ax=axes[0], fraction=0.046, pad=0.04, label="correlation")
    axes[0].set_yticks(range(len(ridx)))
    axes[0].set_yticklabels([f"{row_block}{k}" for k in range(len(ridx))], fontsize=8)
    axes[0].set_ylabel("velocity node")
    # column block guides + labels
    for b, s, e in spans:
        if s > 0:
            axes[0].axvline(s - 0.5, color="k", lw=0.8, alpha=0.5)
        axes[0].text((s + e) / 2 - 0.5, -0.6, b, ha="center", va="bottom", fontsize=9)
    axes[0].set_xlabel("station-correction column")
    axes[0].set_title("cross-correlation", fontsize=10)

    # ---- right: per-node max |corr| against each phase ----
    width = 0.38
    y = np.arange(len(ridx))
    handles = []
    for i, (b, s, e) in enumerate(spans):
        maxc = np.abs(sub[:, s:e]).max(axis=1) if e > s else np.zeros(len(ridx))
        h = axes[1].barh(y + (i - (len(spans)-1)/2.0) * width, maxc, height=width,
                         label=b)
        handles.append(h)
    axes[1].set_yticks(y)
    axes[1].set_yticklabels([f"{row_block}{k}" for k in range(len(ridx))], fontsize=8)
    axes[1].invert_yaxis()
    axes[1].set_xlabel("max |correlation|")
    axes[1].set_xlim(0, 1)
    axes[1].axvline(0.3, color="0.5", ls="--", lw=0.8)
    axes[1].legend(fontsize=8, loc="lower right")
    axes[1].set_title("peak coupling per node", fontsize=10)

    fig.suptitle(title, fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(outpath, dpi=fig.dpi)
    plt.close(fig)
    return outpath


def fig_location_no_ot(Cm, R, colmap, outpath, per_event=4, label=""):
    """Location sub-view over hypocentre x,y,z only, dropping origin-time columns.

    In mcmc_eq the origin time is NOT a sampled free parameter: forward.c sets it
    analytically to the mean travel-time residual (m->origin = -sum), and gmatrix
    encodes its column as a hard-coded analytic +1. So the origin-time rows/cols in
    Cm/R are not a genuine sampled degree of freedom. This view keeps only the
    x,y,z columns of every event's [x y z ot] group. Because we simply omit the ot
    rows/cols from the already-inverted Cm, the x,y,z (co)variances shown here are
    *marginalized over* origin time -- i.e. they reflect the location uncertainty
    with origin time free to absorb error, matching what the MCMC posterior samples.
    """
    if "loc" not in colmap:
        return None
    if per_event < 4:
        # locations were held fixed (loc block is origin-time only) -> nothing to do
        return None
    lo, hi = colmap["loc"]
    # keep components 0,1,2 (x,y,z) of each length-`per_event` group; drop comp 3 (ot)
    idx = [c for c in range(lo, hi + 1) if (c - lo) % per_event != (per_event - 1)]
    if not idx:
        return None
    idx = np.array(idx)
    sub_C = Cm[np.ix_(idx, idx)]
    sub_R = R[np.ix_(idx, idx)]
    Corr = correlation_from_cov(sub_C)
    ne = len(idx) // 3

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.4))
    im0 = _diverging(axes[0], Corr, "Correlation", vmax=1.0)
    fig.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04)
    im1 = _diverging(axes[1], sub_R, "Resolution R", vmax=1.0)
    fig.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)

    # label the three x/y/z sub-bands (columns are ordered x0,y0,z0,x1,y1,z1,...
    # -> after dropping ot the ordering is still interleaved per event, so we just
    #    annotate with a single "x,y,z per event (ot dropped)" note)
    n = len(idx)
    for ax in axes:
        ax.set_xticks([n / 2.0 - 0.5])
        ax.set_xticklabels([f"loc x,y,z  ({ne} events, ot dropped)"], fontsize=8)
        ax.set_yticks([n / 2.0 - 0.5])
        ax.set_yticklabels(["loc x,y,z"], fontsize=8, rotation=90, va="center")

    fig.suptitle(_titled("Location block, origin-time excluded (x, y, z only)", label), fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(outpath, dpi=fig.dpi)
    plt.close(fig)
    return outpath


def fig_resolution_diagonal(R, colmap, outpath, label=""):
    diag = np.diag(R)
    blocks = _block_boundaries(colmap)
    n = len(diag)

    fig, ax = plt.subplots(figsize=(13, 4.2))
    ax.plot(np.arange(n), diag, lw=0.6, color="0.3")
    ax.scatter(np.arange(n), diag, s=3, c=diag, cmap="viridis", vmin=0, vmax=1)
    ax.axhline(1.0, color="g", lw=0.8, ls="--", alpha=0.7, label="perfect (=1)")
    ax.axhline(0.0, color="r", lw=0.8, ls="--", alpha=0.5)

    # shade / label blocks
    ymax = max(1.05, np.nanmax(diag) * 1.05)
    for i, (name, lo, hi) in enumerate(blocks):
        ax.axvspan(lo - 0.5, hi + 0.5, color="C%d" % (i % 10), alpha=0.06)
        ax.text((lo + hi) / 2, ymax * 0.97, name, ha="center", va="top",
                fontsize=8)
        ax.axvline(hi + 0.5, color="k", lw=0.4, alpha=0.3)

    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(min(-0.05, np.nanmin(diag) - 0.05), ymax)
    ax.set_xlabel("model parameter index")
    ax.set_ylabel("resolution diagonal  R_ii")
    ax.set_title(_titled("Resolution diagonal (1 = fully resolved, 0 = unresolved)", label))
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(outpath, dpi=fig.dpi)
    plt.close(fig)
    return outpath


def fig_singular_spectrum(S, ntrunc, rcond, data_smax, outpath, label=""):
    S = np.asarray(S, dtype=float)
    n = len(S)
    smax = S.max()
    tol = data_smax * rcond
    fig, ax = plt.subplots(figsize=(8, 5))
    # clip tiny/negative eigenvalues for log display
    S_disp = np.clip(S, smax * 1e-20, None)
    ax.semilogy(np.arange(n), S_disp, ".", ms=3, color="0.2")
    ax.axhline(tol, color="r", ls="--", lw=1.0,
               label=f"truncation tol = data_smax·{rcond:g}")
    if data_smax < smax:
        ax.axhline(data_smax, color="C0", ls=":", lw=1.0,
                   label="largest data eigenvalue")
    if ntrunc > 0:
        ax.axvspan(n - ntrunc - 0.5, n - 0.5, color="r", alpha=0.10,
                   label=f"{ntrunc} truncated modes")
    ax.set_xlabel("mode index (descending eigenvalue)")
    ax.set_ylabel("eigenvalue of G$^T$WG")
    ax.set_title(_titled("Spectrum of the weighted normal matrix", label))
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outpath, dpi=fig.dpi)
    plt.close(fig)
    return outpath


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def _resolve_label(explicit, inp_path):
    """Return a chain/run label for figure titles.

    Priority: explicit --label > chain id parsed from an rjx-* file in the input
    directory (e.g. rjx-009_4621867.out -> '009') > input filename stem.
    """
    if explicit:
        return explicit
    # 1) chain id embedded in the gmatrix output filename: gmatrix_<id>.out
    stem = os.path.splitext(os.path.basename(inp_path))[0]
    mgm = re.match(r"gmatrix[_-](\w+?)(?:_.*)?$", stem)
    if mgm and mgm.group(1) not in ("global", "loc"):
        return mgm.group(1)
    # 2) otherwise, a chain id parsed from an rjx-* file in the input directory
    d = os.path.dirname(os.path.abspath(inp_path)) or "."
    try:
        rjx = [f for f in os.listdir(d) if f.startswith("rjx-")]
    except OSError:
        rjx = []
    for f in sorted(rjx):
        mtag = re.search(r"rjx-(\w+?)(?:_|\.)", f)
        if mtag:
            return mtag.group(1)
    return ""


def _titled(title, label):
    """Prefix a figure title with the chain/run label, if any."""
    return f"[{label}]  {title}" if label else title


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("gmatrix_out", help="path to global-mode gmatrix.out")
    ap.add_argument("--outdir", default=None,
                    help="directory for figures (default: alongside input)")
    ap.add_argument("--rcond", type=float, default=1e-8,
                    help="SVD truncation ratio relative to largest eigenvalue")
    ap.add_argument("--dpi", type=int, default=130, help="figure DPI")
    ap.add_argument("--save-npz", action="store_true",
                    help="also save Cm, R, spectrum to a .npz")
    ap.add_argument("--label", default=None,
                    help="label added to every figure title (e.g. a chain id like "
                         "'009') for comparing runs. If omitted, auto-detected from "
                         "an rjx-* file in the input directory when possible.")
    args = ap.parse_args(argv)

    inp = args.gmatrix_out
    if not os.path.isfile(inp):
        ap.error(f"file not found: {inp}")
    outdir = args.outdir or os.path.dirname(os.path.abspath(inp)) or "."
    os.makedirs(outdir, exist_ok=True)
    base = os.path.splitext(os.path.basename(inp))[0]
    plt.rcParams["figure.dpi"] = args.dpi

    label = _resolve_label(args.label, inp)
    if label:
        print(f"figure label: {label}")
        # put the chain/run label into the output filename prefix so figures from
        # different chains don't overwrite each other (e.g. gmatrix_global_009_full.png)
        safe = re.sub(r"[^A-Za-z0-9._-]", "_", label)
        if safe not in base.split("_"):
            base = f"{base}_{safe}"

    print(f"parsing {inp} ...")
    g = parse_gmatrix(inp)
    if not g.colmap:
        ap.error("no `# colmap` header found -- is this a GLOBAL-mode gmatrix.out? "
                 "(loc-mode output is per-event and not handled by this script)")
    print(f"  ncols={g.ncols}  data_rows={g.ndata_rows}  "
          f"triplets={len(g.vals)}  blocks={list(g.colmap)}")

    print("computing covariance & resolution (truncated SVD) ...")
    Cm, R, S, ntrunc, data_rows, data_smax = compute_cov_res(g, rcond=args.rcond)

    produced = []
    produced.append(fig_full_matrices(
        Cm, R, g.colmap, os.path.join(outdir, f"{base}_full.png"), label=label))
    produced.append(fig_resolution_diagonal(
        R, g.colmap, os.path.join(outdir, f"{base}_resdiag.png"), label=label))
    produced.append(fig_singular_spectrum(
        S, ntrunc, args.rcond, data_smax,
        os.path.join(outdir, f"{base}_spectrum.png"), label=label))

    # focused block sub-views (only if those blocks exist)
    vel_blocks = [b for b in ("vp", "vpvs") if b in g.colmap]
    vel_title = {("vp", "vpvs"): "Velocity block (vp, vp/vs)",
                 ("vp",): "Velocity block (vp only; vp/vs fixed)",
                 ("vpvs",): "Velocity block (vp/vs only; vp fixed)"}.get(
                     tuple(vel_blocks), "Velocity block")
    vel = fig_block_view(Cm, R, g.colmap, vel_blocks, _titled(vel_title, label),
                         os.path.join(outdir, f"{base}_velocity.png")) if vel_blocks else None
    if vel:
        produced.append(vel)
    sta = fig_block_view(Cm, R, g.colmap, ["pres", "sres"],
                         _titled("Station-correction block (P res, S res)", label),
                         os.path.join(outdir, f"{base}_stacorr.png"))
    if sta:
        produced.append(sta)

    # vp <-> station-correction cross-coupling (the velocity/static trade-off)
    vpsc = fig_cross_block(
        Cm, g.colmap, "vp", ["pres", "sres"],
        _titled("vp \u2194 station-correction trade-off", label),
        os.path.join(outdir, f"{base}_vp_statcor.png"))
    if vpsc:
        produced.append(vpsc)

    # quake depth (z) <-> velocity cross-coupling (the depth/velocity trade-off)
    depth_vel_blocks = [b for b in ("vp", "vpvs") if b in g.colmap]
    if depth_vel_blocks:
        dvel = fig_depth_cross(
            Cm, g.colmap, depth_vel_blocks,
            _titled("quake depth \u2194 velocity trade-off", label),
            os.path.join(outdir, f"{base}_depth_vel.png"),
            per_event=g.loc_per_event)
        if dvel:
            produced.append(dvel)

    # quake depth (z) <-> station-correction cross-coupling
    dsc = fig_depth_cross(
        Cm, g.colmap, ["pres", "sres"],
        _titled("quake depth \u2194 station-correction trade-off", label),
        os.path.join(outdir, f"{base}_depth_statcor.png"),
        per_event=g.loc_per_event)
    if dsc:
        produced.append(dsc)

    # location block only (hypocentre x,y,z + origin-time for every event)
    loc_title = ("Location block (origin-time only; quake locations fixed)"
                 if g.loc_per_event < 4
                 else "Location block (x, y, z, origin-time per event)")
    loc = fig_block_view(Cm, R, g.colmap, ["loc"], _titled(loc_title, label),
                         os.path.join(outdir, f"{base}_location.png"))
    if loc:
        produced.append(loc)

    # location block WITHOUT origin time (ot is analytic in mcmc_eq, not sampled)
    loc_no_ot = fig_location_no_ot(
        Cm, R, g.colmap, os.path.join(outdir, f"{base}_location_noOT.png"),
        per_event=g.loc_per_event, label=label)
    if loc_no_ot:
        produced.append(loc_no_ot)

    # everything except location, together (velocity + station corrections),
    # so the cross-block vp / vp-vs / pres / sres coupling is visible in one view
    nonloc_blocks = [b for b in ("vp", "vpvs", "pres", "sres") if b in g.colmap]
    nonloc = fig_block_view(Cm, R, g.colmap, nonloc_blocks,
                            _titled("Non-location parameters (velocity + station corrections)", label),
                            os.path.join(outdir, f"{base}_nonloc.png"))
    if nonloc:
        produced.append(nonloc)

    if args.save_npz:
        npz = os.path.join(outdir, f"{base}_cov_res.npz")
        np.savez_compressed(npz, Cm=Cm, R=R, spectrum=S,
                            colmap=np.array(list(g.colmap.items()), dtype=object))
        print(f"  saved arrays -> {npz}")

    print("\nfigures written:")
    for p in produced:
        print(f"  {p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
