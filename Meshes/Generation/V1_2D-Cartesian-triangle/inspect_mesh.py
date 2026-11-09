#!/usr/bin/env python3
"""
inspect_mesh.py  -  check a mesh written by mesh2d.py (Gmsh 2.2 ASCII, tri + quad).

Checks
  nodes        non-finite values, coincident nodes, unused nodes
  cells        repeated node ids, zero / negative area, near-zero-area slivers,
               clockwise cells, duplicate cells, very small angles
  connectivity edges shared by more than two cells, boundary nodes that are not on a
               single closed loop, zero-area boundary loops (hanging nodes),
               disconnected pieces
  area         sum of cell areas against (domain box - body)

Usage
  python3 inspect_mesh.py mesh.msh
  python3 inspect_mesh.py mesh.msh --plot            # writes <mesh>_inspect.png
  python3 inspect_mesh.py mesh.msh --area-tol 1e-10  # zero-volume threshold (relative to median area)

Exit code: 0 = no errors, 1 = errors found (warnings do not change the exit code).
"""

import argparse
import sys

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree


def read_msh(fname):
    nodes = {}
    tris, quads, lines = [], [], []
    with open(fname) as f:
        lns = f.read().split("\n")
    k = 0
    while k < len(lns):
        s = lns[k].strip()
        if s == "$Nodes":
            n = int(lns[k + 1])
            for r in range(n):
                p = lns[k + 2 + r].split()
                nodes[int(p[0])] = (float(p[1]), float(p[2]))
            k += n + 2
        elif s == "$Elements":
            n = int(lns[k + 1])
            for r in range(n):
                p = lns[k + 2 + r].split()
                et = int(p[1])
                nt = int(p[2])
                ids = [int(v) for v in p[3 + nt:]]
                if et == 2:
                    tris.append(ids)
                elif et == 3:
                    quads.append(ids)
                elif et == 1:
                    lines.append(ids + [int(p[3])])
            k += n + 2
        k += 1
    if not nodes:
        sys.exit("error: no nodes found in %s" % fname)
    ids = np.array(sorted(nodes))
    xy = np.array([nodes[i] for i in ids])
    lut = {i: n for n, i in enumerate(ids)}
    tris = np.array([[lut[v] for v in t] for t in tris], dtype=int).reshape(-1, 3)
    quads = np.array([[lut[v] for v in q] for q in quads], dtype=int).reshape(-1, 4)
    lines = np.array([[lut[a], lut[b], t] for a, b, t in lines], dtype=int).reshape(-1, 3)
    return xy, tris, quads, lines


def poly_area(P):
    x = P[:, :, 0]
    y = P[:, :, 1]
    return 0.5 * np.sum(x * np.roll(y, -1, axis=1) - np.roll(x, -1, axis=1) * y, axis=1)


def min_angle_deg(P):
    a = np.roll(P, 1, axis=1) - P
    b = np.roll(P, -1, axis=1) - P
    la = np.linalg.norm(a, axis=2)
    lb = np.linalg.norm(b, axis=2)
    with np.errstate(invalid="ignore", divide="ignore"):
        c = np.sum(a * b, axis=2) / (la * lb)
    c = np.clip(c, -1, 1)
    ang = np.degrees(np.arccos(c))
    ang = np.where(np.isfinite(ang), ang, 0.0)
    return ang.min(axis=1)


class Report:
    def __init__(self):
        self.errors = 0
        self.warnings = 0

    def ok(self, msg):
        print("  [ ok ] " + msg)

    def warn(self, msg):
        self.warnings += 1
        print("  [warn] " + msg)

    def err(self, msg):
        self.errors += 1
        print("  [FAIL] " + msg)


def examples(ids, cent, n=5):
    out = []
    for c in ids[:n]:
        out.append("cell %d at (%.5g, %.5g)" % (c + 1, cent[c, 0], cent[c, 1]))
    return "; ".join(out) + (" ..." if len(ids) > n else "")


def main():
    ap = argparse.ArgumentParser(description="Inspect a 2D mesh for degenerate / zero-area cells")
    ap.add_argument("mesh", help="mesh file (.msh, Gmsh 2.2 ASCII)")
    ap.add_argument("--area-tol", type=float, default=1e-12,
                    help="area <= tol * median area counts as zero volume (default 1e-12)")
    ap.add_argument("--sliver", type=float, default=1e-3,
                    help="area < sliver * median area is reported as a sliver (default 1e-3)")
    ap.add_argument("--min-angle", type=float, default=1.0,
                    help="cells with an interior angle below this (deg) are reported (default 1)")
    ap.add_argument("--bodies", type=int, default=1,
                    help="expected number of bodies (interior boundary loops) (default 1)")
    ap.add_argument("--plot", action="store_true", help="write <mesh>_inspect.png with bad cells in red")
    args = ap.parse_args()

    xy, tris, quads, lines = read_msh(args.mesh)
    nt, nq = len(tris), len(quads)
    print("mesh: %s" % args.mesh)
    print("  nodes %d, triangles %d, quads %d, boundary lines %d" % (len(xy), nt, nq, len(lines)))
    R = Report()

    # ---- nodes ---------------------------------------------------------------
    print("nodes")
    if not np.all(np.isfinite(xy)):
        R.err("non-finite node coordinates: %d" % (~np.isfinite(xy).all(axis=1)).sum())
    else:
        R.ok("all coordinates finite")
    ext = np.ptp(xy, axis=0).max()
    pairs = cKDTree(xy).query_pairs(1e-10 * ext, output_type="ndarray")
    if len(pairs):
        R.err("%d pairs of coincident nodes (e.g. nodes %d and %d)" % (len(pairs), pairs[0, 0] + 1, pairs[0, 1] + 1))
    else:
        R.ok("no coincident nodes")
    used = np.zeros(len(xy), dtype=bool)
    used[tris.ravel()] = True
    used[quads.ravel()] = True
    if (~used).any():
        R.warn("%d unused nodes" % (~used).sum())
    else:
        R.ok("no unused nodes")

    # ---- cells ---------------------------------------------------------------
    print("cells")
    cells = [("tri", tris), ("quad", quads)]
    allc = []
    for name, c in cells:
        for row in c:
            allc.append(row)
    area = np.concatenate([poly_area(xy[c]) if len(c) else np.zeros(0) for _, c in cells])
    ang = np.concatenate([min_angle_deg(xy[c]) if len(c) else np.zeros(0) for _, c in cells])
    cent = np.concatenate([xy[c].mean(axis=1) if len(c) else np.zeros((0, 2)) for _, c in cells])
    ncell = len(area)
    bad = np.zeros(ncell, dtype=bool)

    if ncell == 0:
        R.err("mesh has no cells")
    else:
        med = np.median(np.abs(area))
        rep = np.array([len(set(r)) != len(r) for r in allc])
        if rep.any():
            R.err("%d cells repeat a node: %s" % (rep.sum(), examples(np.nonzero(rep)[0], cent)))
            bad |= rep
        else:
            R.ok("no cells with repeated nodes")

        zero = np.abs(area) <= args.area_tol * med
        if zero.any():
            R.err("%d zero-volume cells: %s" % (zero.sum(), examples(np.nonzero(zero)[0], cent)))
            bad |= zero
        else:
            R.ok("no zero-volume cells (|area| <= %.1e x median)" % args.area_tol)

        neg = area < -args.area_tol * med
        if neg.any():
            R.err("%d clockwise / negative-area cells: %s" % (neg.sum(), examples(np.nonzero(neg)[0], cent)))
            bad |= neg
        else:
            R.ok("all cells counter-clockwise")

        sl = (~zero) & (np.abs(area) < args.sliver * med)
        if sl.any():
            R.warn("%d sliver cells (area < %.0e x median): %s" % (sl.sum(), args.sliver, examples(np.nonzero(sl)[0], cent)))
            bad |= sl
        else:
            R.ok("no sliver cells (area < %.0e x median)" % args.sliver)

        sm = (~zero) & (ang < args.min_angle)
        if sm.any():
            R.warn("%d cells with an interior angle below %.1f deg (min %.3g): %s"
                   % (sm.sum(), args.min_angle, ang.min(), examples(np.nonzero(sm)[0], cent)))
            bad |= sm
        else:
            R.ok("no interior angle below %.1f deg" % args.min_angle)

        seen = {}
        dup = 0
        for r in allc:
            kk = tuple(sorted(int(v) for v in r))
            seen[kk] = seen.get(kk, 0) + 1
        dup = sum(v - 1 for v in seen.values() if v > 1)
        if dup:
            R.err("%d duplicate cells" % dup)
        else:
            R.ok("no duplicate cells")
        print("  area: min %.4g, median %.4g, max %.4g; smallest interior angle %.3g deg"
              % (np.abs(area).min(), med, np.abs(area).max(), ang.min()))

    # ---- connectivity ---------------------------------------------------------
    print("connectivity")
    E = []
    owner = []
    off = 0
    for name, c in cells:
        if len(c):
            k = c.shape[1]
            for s in range(k):
                E.append(np.column_stack((c[:, s], c[:, (s + 1) % k])))
                owner.append(np.arange(len(c)) + off)
        off += len(c)
    E = np.vstack(E)
    owner = np.concatenate(owner)
    ek = np.sort(E, axis=1)
    uk, inv, cnt = np.unique(ek, axis=0, return_inverse=True, return_counts=True)
    inv = inv.ravel()
    nm = cnt > 2
    if nm.any():
        R.err("%d edges shared by more than two cells (e.g. nodes %d-%d)" % (nm.sum(), uk[nm][0, 0] + 1, uk[nm][0, 1] + 1))
        for e in np.nonzero(nm[inv])[0]:
            bad[owner[e]] = True
    else:
        R.ok("no edge is shared by more than two cells")

    bmask = cnt[inv] == 1
    BE = E[bmask]
    if len(BE) == 0:
        R.err("no boundary edges: mesh has no outer boundary")
    else:
        deg = np.bincount(BE.ravel(), minlength=len(xy))
        bnodes = np.nonzero(deg)[0]
        wrong = bnodes[deg[bnodes] != 2]
        if len(wrong):
            R.err("%d boundary nodes are not on a single closed loop (degree != 2), e.g. node %d at (%.5g, %.5g)"
                  % (len(wrong), wrong[0] + 1, xy[wrong[0], 0], xy[wrong[0], 1]))
        else:
            R.ok("every boundary node lies on a closed loop")

        n = len(xy)
        A = coo_matrix((np.ones(len(BE)), (BE[:, 0], BE[:, 1])), shape=(n, n))
        nc, lab = connected_components(A, directed=False)
        loopA = {}
        for a, b in BE:
            L = lab[a]
            loopA[L] = loopA.get(L, 0.0) + 0.5 * (xy[a, 0] * xy[b, 1] - xy[b, 0] * xy[a, 1])
        loops = [L for L in loopA if deg[np.nonzero(lab == L)[0][0]] > 0]
        scale = np.ptp(xy, axis=0).max() ** 2
        flat = [L for L in loops if abs(loopA[L]) < 1e-12 * scale]
        if flat:
            R.err("%d zero-area boundary loops (hanging nodes / T-junctions)" % len(flat))
            for L in flat:
                for a, b in BE:
                    if lab[a] == L:
                        print("         near (%.5g, %.5g)" % tuple(xy[a]))
                        break
        else:
            R.ok("no zero-area boundary loops (no hanging nodes)")
        print("  boundary loops: %d, boundary edges: %d" % (len(loops), len(BE)))

    # disconnected pieces (cells connected through shared edges)
    first = {}
    pairs = []
    for e in range(len(E)):
        k = inv[e]
        if k in first:
            pairs.append((first[k], owner[e]))
        else:
            first[k] = owner[e]
    if ncell:
        if pairs:
            pa = np.array(pairs)
            Gc = coo_matrix((np.ones(len(pa)), (pa[:, 0], pa[:, 1])), shape=(ncell, ncell))
            npieces, _ = connected_components(Gc, directed=False)
        else:
            npieces = ncell
        if npieces > 1:
            R.warn("mesh has %d disconnected pieces" % npieces)
        else:
            R.ok("mesh is one connected piece")

    # ---- area conservation ------------------------------------------------------
    print("area")
    if len(BE):
        tol = 1e-9 * ext
        x0, x1 = xy[:, 0].min(), xy[:, 0].max()
        y0, y1 = xy[:, 1].min(), xy[:, 1].max()
        mid = 0.5 * (xy[BE[:, 0]] + xy[BE[:, 1]])
        on_dom = ((np.abs(mid[:, 0] - x0) < tol) | (np.abs(mid[:, 0] - x1) < tol)
                  | (np.abs(mid[:, 1] - y0) < tol) | (np.abs(mid[:, 1] - y1) < tol))
        body = BE[~on_dom]
        n_body_loops = len(np.unique(lab[body[:, 0]])) if len(body) else 0
        if n_body_loops != args.bodies:
            R.err("found %d closed interior boundary loop(s), expected %d (a missing cell shows up as an extra loop)"
                  % (n_body_loops, args.bodies))
        else:
            R.ok("%d interior boundary loop(s), as expected" % n_body_loops)
        sa = 0.5 * np.sum(xy[body[:, 0], 0] * xy[body[:, 1], 1] - xy[body[:, 1], 0] * xy[body[:, 0], 1])
        body_area = -sa                       # body loops run clockwise (fluid on the left)
        expect = (x1 - x0) * (y1 - y0) - body_area
        got = area.sum()
        rel = abs(got - expect)/expect
        msg = "cell area %.8g, domain box - body %.8g (relative difference %.2e)" % (got, expect, rel)
        if rel > 1e-8:
            R.err(msg + "  -> gaps or overlaps in the mesh")
        else:
            R.ok(msg)

    print("")
    print("result: %d error(s), %d warning(s)" % (R.errors, R.warnings))

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.collections import PolyCollection
        polys = [xy[c] for _, cs in cells for c in cs]
        fc = ["#ff4040" if b else "none" for b in bad]
        fig, ax = plt.subplots(1, 2, figsize=(13, 6))
        for a in ax:
            a.add_collection(PolyCollection(polys, facecolors=fc, edgecolors="k", linewidths=0.25))
            a.set_aspect("equal")
        ax[0].set_xlim(xy[:, 0].min(), xy[:, 0].max())
        ax[0].set_ylim(xy[:, 1].min(), xy[:, 1].max())
        ax[0].set_title("full domain (flagged cells in red)")
        if bad.any():
            c = cent[bad]
            w = max(np.ptp(c[:, 0]), np.ptp(c[:, 1]), 20 * np.median(np.sqrt(np.abs(area)))) * 0.75
            m = c.mean(axis=0)
            ax[1].set_xlim(m[0] - w, m[0] + w)
            ax[1].set_ylim(m[1] - w, m[1] + w)
            ax[1].set_title("around the flagged cells")
        else:
            ax[1].set_title("no flagged cells")
        out = args.mesh.rsplit(".", 1)[0] + "_inspect.png"
        fig.tight_layout()
        fig.savefig(out, dpi=130)
        print("plot written: %s" % out)

    sys.exit(1 if R.errors else 0)


if __name__ == "__main__":
    main()
