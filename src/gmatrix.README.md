# gmatrix — sensitivity matrix G for covariance and resolution

`gmatrix` builds the linearized sensitivity (Jacobian) matrix

    G[i][j] = d(predicted_time_i) / d(model_param_j)

evaluated at a model taken from an `mcmc_eq` output chain. Derivatives are computed
by finite differences **through the same forward code the inversion uses**
(`forward.c`: `cal_fit_newx` / `traveltimet` / `setup_table_new`), so G is exactly
consistent with the travel-time tables the model produced.

Once you have G you can form the model **covariance** and **resolution** matrices
downstream (see recipes below).

## Build

    make gmatrix

## Usage

    gmatrix config.dat chain_file pick_file [mode] [tag] [outfile] [fdstep] [cweight] [dvp] [dvpvs]

- `mode`    : `loc` (per-event location G, default) or `global` (full 1c matrix)
- `tag`     : `bat` (best-fit model, default) or `mod` (last sampled model)
- `outfile` : output path (default `gmatrix.out`)
- `fdstep`  : central-difference step in km for x,y,z (default 0.01)
- `cweight` : weight of zero-mean station-correction constraint rows (default 1e6)
- `dvp`     : central-difference step in km/s for vp columns (default 0.02)
- `dvpvs`   : central-difference step for vp/vs columns (default 0.01)

The config's `eikonal`, `TRIA`, `scor_flag`, and `reference_station` are read
automatically and must match how the chain was produced.

## mode = loc (per-event location G)

Four parameters per event `[x, y, z, origin_time]`. Each output row is one pick:

    <evt> <P|S> st <id> cl <class> sigma <s>  G <dTp/dx> <dTp/dy> <dTp/dz> <dTp/dot>

For each event, stack its rows into a matrix `Ge` (npick x 4) and its pick noise
into `sigma`. Then:

    Cd^-1 = diag(1/sigma^2)
    Cm    = (Ge^T Cd^-1 Ge)^-1            # 4x4 location+origin-time covariance

The diagonal of `Cm` gives variances of x, y, z, origin-time; off-diagonals give
trade-offs (e.g. depth vs origin-time).

## mode = global (full 1c matrix)

Parameter vector (column order, fixed):

    [ (x,y,z,ot)_e for e=0..ne-1 | vp_k k=0..dim-1 | vpvs_k k=0..dim-1
      | pres_s s=0..nos-1 | sres_s s=0..nos-1 ]

The manifest header lists the exact column ranges, e.g.

    # colmap loc 0..879 (4/event: x y z ot)  vp 880..894  vpvs 895..909  pres 910..1039  sres 1040..1169

Output records:

    D <row> <P|S> <evt> <st> <cl> <sigma> <resid>   # one per data row (pick)
    G <row> <col> <value>                           # sparse Jacobian triplets
    C <row> <meanP|meanS> <weight>                  # constraint pseudo-row meta

- Location columns: central FD, block-diagonal per event; origin-time column = +1.
- Station-correction columns: analytic +1 (`pres_s` for a P pick at station s;
  `sres_s` for an S pick). `d(t)/d(vpvs)` is exactly 0 for P picks (S-only).
- Velocity columns: central FD. `vp_k` rebuilds both P and S tables (Vs=Vp/(Vp/Vs));
  `vpvs_k` rebuilds S only.

### Fixed parameters (JHD / input-from-file runs)

If the run held some parameter classes fixed (config line 34 = `3 <switch>`, i.e.
`aflag==3` "input from model file"), those classes are **not free parameters** and
`gmatrix` omits their columns so the covariance/resolution reflect the actual
posterior. The switch letters map as:

- `V` velocity fixed  -> **vp and vpvs columns dropped**
- `Q` quake locations fixed -> **x,y,z columns dropped**; the per-event origin-time
  column is kept (mcmc_eq still sets it analytically), so the `loc` block becomes
  1 column/event instead of 4.
- `R` station corrections fixed -> **pres/sres columns dropped** (and the zero-mean
  constraint pseudo-rows are skipped, since there are no correction columns).

Examples: `3 V` (Example3) drops velocity; `3 QV` drops velocity and locations.
The manifest header records this so downstream tools adapt automatically:

    # colmap loc 0..879 (4/event: x y z ot)  pres 880..1009  sres 1010..1139   # 3 V
    # colmap loc 0..219 (1/event: ot)  pres 220..349  sres 350..479             # 3 QV
    # fixed  vel 1  loc 0  res 0   (1 = held fixed / not a free parameter)

`gmatrix_cov_res.py` reads the `colmap`/`fixed` header, so it only draws figures for
blocks that exist (e.g. no `*_velocity.png` when velocity is fixed) and skips the
`*_location_noOT.png` figure when locations are fixed (the loc block is already
origin-time only). For loc-mode with fixed quake locations, `gmatrix` prints a
warning because per-event location covariance is not meaningful there.

### Station-correction constraint handling (config-driven)

The MCMC does not sample corrections freely; `gmatrix` mirrors what the run did:

- `scor_flag == 0` (zero-mean): all columns kept, plus two **constraint pseudo-rows**
  (`meanP`, `meanS`) each with entries = `cweight` on every `pres` (resp. `sres`)
  column and a target residual of 0. This removes the exact null space
  (adding a constant to all corrections) that would otherwise make `G^T G` singular,
  while keeping the covariance in the **same gauge as the posterior samples**.
- `scor_flag == 1` (P reference fixed): the reference station's `pres` column is dropped.
- `scor_flag == 2` (P & S reference fixed): reference `pres` and `sres` columns dropped.
- `scor_flag == -1` (invert P only): `sres` columns omitted.
- `scor_flag == -2` (invert S only): `pres` columns omitted.

`cweight` sets the strength: large (default 1e6) = near-hard constraint; set it to
`1/sigma_prior` for a soft Bayesian prior instead.

## Downstream: covariance and resolution

Assemble the sparse `G` (data rows + constraint rows), the data weight
`Wd = diag(1/sigma^2)` over data rows, and the constraint weight over pseudo-rows.
Let `W` be the combined row-weight (data + constraint). Then, using
`Gw = sqrt(W) G` and `dw = sqrt(W) d`:

Covariance (posterior model covariance, Tarantola):

    Cm = (G^T W G)^-1

Resolution (Menke model resolution matrix):

    G^-g = (G^T W G)^-1 G_data^T Wd        # generalized inverse (data part only)
    R    = G^-g G_data

`R` should be close to identity where parameters are well resolved; row/column
spread shows smearing between parameters. Use only the **data** rows in `G_data`
for the `G^-g G_data` product — the constraint rows regularize the inverse
`(G^T W G)^-1` but are not "data" being mapped.

Because a full 1c problem generally has weakly-resolved directions (especially deep
velocity layers with little ray coverage), form the inverse with a **truncated SVD**
of `G^T W G` (or of `Gw` directly) and inspect the singular-value spectrum. The
constraint rows fix the *known* exact null spaces (station-correction mean); SVD
truncation handles the remaining weakly-resolved ones. Report which singular values
were truncated alongside `Cm`/`R`.

### Standalone visualization script

For global-mode output, use the ready-made script `gmatrix_cov_res.py` (in `src/`)
instead of the sketch below. It parses `gmatrix.out`, forms `Cm` and `R` with a
truncated SVD, and writes a set of figures:

    python3 src/gmatrix_cov_res.py Example/gmatrix_global.out
    # options: --outdir DIR  --rcond 1e-8  --dpi 130  --save-npz

Figures written next to the input (prefix = input basename):

- `*_full.png`     : log10|Cm|, the model correlation matrix, and R (full 1c),
                     with block guides for loc / vp / vpvs / pres / sres.
- `*_resdiag.png`  : the resolution diagonal R_ii per parameter (1 = resolved).
- `*_spectrum.png` : eigenvalue spectrum of G^T W G with the truncation level.
- `*_velocity.png` : correlation + resolution restricted to the vp / vpvs block.
- `*_stacorr.png`  : correlation + resolution restricted to the pres / sres block.

`--save-npz` also dumps `Cm`, `R`, and the spectrum to a compressed `.npz`.

**Truncation caveat (important).** With `scor_flag == 0` the zero-mean constraint
rows carry weight `cweight` (default 1e6), so the normal matrix `G^T W G` picks up
a few eigenvalues of order `cweight^2` (~1e12-1e14) that dwarf every data-resolved
direction. If you truncate relative to the global largest eigenvalue you will throw
away all the real data modes. The script handles this by referencing the truncation
tolerance to the largest *data* eigenvalue (constraint modes ~`cweight^2` are
detected and excluded from that reference), so `--rcond` behaves as expected. On the
Example chain this truncates only ~6 genuinely weak modes, and `trace(R)` approaches
`ncols`.

### Minimal Python sketch

```python
import numpy as np
from scipy.sparse import coo_matrix

rows, cols, vals = [], [], []
sigma = {}; is_data = {}
ncols = None
for line in open("gmatrix.out"):
    p = line.split()
    if not p: continue
    if p[0] == "#" and "ncols" in line:
        ncols = int(p[p.index("ncols")+1])
    elif p[0] == "D":
        r = int(p[1]); sigma[r] = float(p[6]); is_data[r] = True
    elif p[0] == "C":
        is_data[int(p[1])] = False        # constraint row
    elif p[0] == "G":
        rows.append(int(p[1])); cols.append(int(p[2])); vals.append(float(p[3]))

nrows = max(rows) + 1
G = coo_matrix((vals, (rows, cols)), shape=(nrows, ncols)).tocsr()

# row weights: 1/sigma for data rows, 1 for constraint rows (weight already in G value)
w = np.ones(nrows)
for r, s in sigma.items():
    w[r] = 1.0 / s
Gw = G.multiply(w[:, None]).tocsr()

GtG = (Gw.T @ Gw).toarray()
# truncated-SVD pseudo-inverse
U, S, Vt = np.linalg.svd(GtG)
tol = S.max() * 1e-8
Sinv = np.array([1/x if x > tol else 0 for x in S])
GtG_inv = (Vt.T * Sinv) @ U.T
Cm = GtG_inv                                  # posterior covariance

# resolution: use data rows only for the mapping term
data_rows = [r for r in range(nrows) if is_data.get(r, False)]
Gd = G[data_rows, :]
Wd = np.array([1.0/sigma[r]**2 for r in data_rows])
Gg = GtG_inv @ (Gd.T.multiply(Wd)).tocsr()    # generalized inverse
R  = (Gg @ Gd).toarray()                      # model resolution matrix
```

## Validation notes

- `d(t)/d(vp) < 0` everywhere (faster medium -> earlier arrival); confirmed on the
  Example chain (all entries negative).
- `d(t)/d(vp/vs) > 0` for S picks, exactly 0 for P picks; confirmed.
- Finite-difference derivatives are stable to ~4 significant figures across step
  sizes (x,y,z at 0.005/0.01/0.02 km; vp at 0.02/0.04 km/s), so the step choice is
  not sensitive.
- Origin-time and station-correction columns are analytic (+1), not differenced.
