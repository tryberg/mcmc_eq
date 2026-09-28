#!/usr/bin/env python3
"""
gmatrix_event_ellipse.py -- formal location error ellipses from loc-mode gmatrix.out

Reads the per-event location Jacobian written by `gmatrix ... loc ...` and, for a
requested event, forms the 4x4 location+origin-time covariance

    Cm = (Ge^T Cd^-1 Ge)^-1 ,   Cd^-1 = diag(1/sigma^2)

then emits GMT `psxy -SE` ellipse parameters for the requested projection plane.
Unlike the axis-aligned marginal std-devs stored in resmcnx.dat, this ellipse uses
the *off-diagonal* covariance, so it is correctly oriented (tilted) and shows the
x-y (or x-z) parameter trade-off.

Output (one line, GMT -SE order):

    <cx> <cy> <azimuth_deg> <major_diam> <minor_diam>

where azimuth is degrees clockwise from +Y (north/up) and the axes are FULL
diameters in the same km units as the plot, scaled to the requested confidence.

Usage:
    gmatrix_event_ellipse.py gmatrix_loc.out EVENT [--plane xy|xz]
                             [--conf 0.95 | --nsigma N] [--center-from-cov]

By default the center is the event hypocentre from the `# EVENT` header. With
`--center-from-cov` no center shift is applied (still uses hypocentre); the flag is
reserved for future use. Prints nothing (empty) and exits 2 if the event has too
few picks to invert.
"""

import argparse
import math
import sys

import numpy as np
from scipy.stats import chi2


PLANE_IDX = {"xy": (0, 1), "xz": (0, 2)}


def parse_event(path, event):
    """Return (hypo(x,y,z), Ge (npick x 4), sigma (npick,)) for the given event."""
    hypo = None
    G = []
    sig = []
    in_evt = False
    target = f"# EVENT {event} "
    with open(path) as fh:
        for line in fh:
            if line.startswith("# EVENT "):
                if in_evt:            # reached the next event: stop
                    break
                if line.startswith(target):
                    in_evt = True
                    p = line.split()
                    # ... hypo <x> <y> <z>
                    hi = p.index("hypo")
                    hypo = (float(p[hi + 1]), float(p[hi + 2]), float(p[hi + 3]))
                continue
            if not in_evt:
                continue
            if line.startswith("#"):
                continue
            p = line.split()
            if len(p) < 12:
                continue
            # <evt> <P|S> st <id> cl <class> sigma <s>  G <gx> <gy> <gz> <got>
            try:
                gi = p.index("G")
                si = p.index("sigma")
            except ValueError:
                continue
            sigma = float(p[si + 1])
            gx, gy, gz, got = (float(p[gi + 1]), float(p[gi + 2]),
                               float(p[gi + 3]), float(p[gi + 4]))
            G.append([gx, gy, gz, got])
            sig.append(sigma)
    if hypo is None:
        raise SystemExit(f"event {event} not found in {path}")
    return hypo, np.array(G), np.array(sig)


def event_covariance(Ge, sigma):
    """4x4 location+ot covariance Cm = (Ge^T Cd^-1 Ge)^-1, or None if singular."""
    if Ge.shape[0] < 4:
        return None
    Cd_inv = 1.0 / (sigma ** 2)
    N = Ge.T @ (Ge * Cd_inv[:, None])         # 4x4 normal matrix
    try:
        Cm = np.linalg.inv(N)
    except np.linalg.LinAlgError:
        return None
    return Cm


def ellipse_from_cov2(C2, scale):
    """2x2 covariance -> (azimuth_deg_from_+Y, major_diam, minor_diam) * scale.

    Eigen-decompose C2; semi-axis lengths = scale*sqrt(eigenvalue). Azimuth is the
    orientation of the MAJOR axis measured clockwise from +Y (GMT -SE convention:
    angle CCW from +X is converted to CW-from-north = 90 - angle_ccw_from_x).
    """
    evals, evecs = np.linalg.eigh(C2)          # ascending
    # major axis = largest eigenvalue
    lam_major, lam_minor = evals[1], evals[0]
    vmaj = evecs[:, 1]
    # angle of major axis CCW from +X (plot axis 0)
    ang_ccw_x = math.degrees(math.atan2(vmaj[1], vmaj[0]))
    # GMT -SE azimuth = degrees CW from +Y  ->  90 - ang_ccw_x
    azimuth = 90.0 - ang_ccw_x
    major_diam = 2.0 * scale * math.sqrt(max(lam_major, 0.0))
    minor_diam = 2.0 * scale * math.sqrt(max(lam_minor, 0.0))
    return azimuth, major_diam, minor_diam


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("gmatrix_loc", help="loc-mode gmatrix.out")
    ap.add_argument("event", type=int, help="event number to plot")
    ap.add_argument("--plane", choices=["xy", "xz"], default="xy")
    grp = ap.add_mutually_exclusive_group()
    grp.add_argument("--conf", type=float, default=None,
                     help="confidence level in (0,1), e.g. 0.95 (2 DOF chi-square)")
    grp.add_argument("--nsigma", type=float, default=None,
                     help="axis scale in sigma units (e.g. 1, 2)")
    ap.add_argument("--print-cov", action="store_true",
                    help="print the full 4x4 covariance to stderr for inspection")
    args = ap.parse_args(argv)

    hypo, Ge, sigma = parse_event(args.gmatrix_loc, args.event)
    Cm = event_covariance(Ge, sigma)
    if Cm is None:
        sys.stderr.write(f"event {args.event}: too few picks / singular; no ellipse\n")
        return 2

    if args.print_cov:
        np.set_printoptions(precision=5, suppress=True)
        sys.stderr.write(f"event {args.event} 4x4 Cm [x y z ot]:\n{Cm}\n")
        sx, sy, sz = (math.sqrt(Cm[0, 0]), math.sqrt(Cm[1, 1]), math.sqrt(Cm[2, 2]))
        sys.stderr.write(f"  formal 1-sigma: sx={sx:.3f} sy={sy:.3f} sz={sz:.3f} km\n")

    # confidence scaling for a 2-DOF error ellipse
    if args.nsigma is not None:
        scale = args.nsigma
    else:
        conf = args.conf if args.conf is not None else 0.95
        scale = math.sqrt(chi2.ppf(conf, df=2))

    i, j = PLANE_IDX[args.plane]
    C2 = np.array([[Cm[i, i], Cm[i, j]],
                   [Cm[j, i], Cm[j, j]]])
    azimuth, major, minor = ellipse_from_cov2(C2, scale)

    cx = hypo[i]
    cy = hypo[j]
    # GMT -SE line: x y azimuth major_diameter minor_diameter
    print(f"{cx:.5f} {cy:.5f} {azimuth:.3f} {major:.5f} {minor:.5f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
