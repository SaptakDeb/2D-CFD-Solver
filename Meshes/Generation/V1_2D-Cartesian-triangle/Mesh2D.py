#!/usr/bin/env python3
"""
mesh2d.py  -  simple 2D mesher: Cartesian grid + triangles around the geometry.

  * The domain is covered by a uniform Cartesian grid of quad cells.
  * Cells fully inside the body are removed.
  * Cells cut by the body surface are clipped to the fluid region and split into
    triangles (constrained Delaunay), so the surface is represented exactly.
  * No inflation layers, no local refinement. Cell quality is not controlled
    (use inspect_mesh.py to look for degenerate / zero-area cells).

Geometry (one closed body):
  --geometry FILE      text file with two columns "x y" (header lines are ignored,
                       works with Selig/Lednicer airfoil .dat files)
  --shape NAME         circle | ellipse | rectangle | naca     (built-in shapes)

Output (basename given by --out):
  <out>.msh            Gmsh 2.2 ASCII (quads + triangles, boundary lines tagged
                       1 = body, 2 = farfield, cells tagged 3 = fluid)
  <out>_nodes.txt      x y                      (1-based ids = row numbers)
  <out>_tri.txt        n1 n2 n3                 (counter-clockwise, 1-based)
  <out>_quad.txt       n1 n2 n3 n4              (counter-clockwise, 1-based)
  <out>_bnd.txt        n1 n2 tag                (1 = body, 2 = farfield)
  <out>.png            preview (with --plot)

Requires: python3, numpy, scipy, shapely >= 2.1, matplotlib (only for --plot)
  pip install numpy scipy "shapely>=2.1" matplotlib

Examples:
  python3 mesh2d.py --shape circle --radius 0.5 --h 0.05 --out cyl --plot
  python3 mesh2d.py --shape naca --naca 2412 --aoa 5 --h 0.02 --out wing
  python3 mesh2d.py --geometry mybody.dat --xmin -5 --xmax 15 --ymin -5 --ymax 5 --h 0.05
"""

import argparse
import sys

import numpy as np
import shapely
from shapely.geometry import Polygon
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree


# ----------------------------------------------------------------------------
#  Geometry
# ----------------------------------------------------------------------------

def load_points(fname):
    pts = []
    with open(fname) as f:
        for line in f:
            parts = line.replace(",", " ").replace(";", " ").split()
            if len(parts) < 2:
                continue
            try:
                pts.append((float(parts[0]), float(parts[1])))
            except ValueError:
                continue
    if len(pts) < 3:
        sys.exit("error: geometry file needs at least 3 numeric 'x y' rows")
    return np.array(pts)


def shape_circle(r, n):
    t = 2 * np.pi * np.arange(n) / n
    return np.column_stack((r * np.cos(t), r * np.sin(t)))


def shape_ellipse(a, b, n):
    t = 2 * np.pi * np.arange(n) / n
    return np.column_stack((a * np.cos(t), b * np.sin(t)))


def shape_rectangle(length, thick):
    L, H = length / 2, thick / 2
    return np.array([[L, -H], [L, H], [-L, H], [-L, -H]])


def shape_naca4(code, chord, n):
    m = int(code[0]) / 100.0
    p = int(code[1]) / 10.0
    t = int(code[2:4]) / 100.0
    beta = np.linspace(0, np.pi, n)
    x = 0.5 * (1 - np.cos(beta))
    yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x**2
                  + 0.2843 * x**3 - 0.1036 * x**4)
    if m > 0 and p > 0:
        yc = np.where(x < p, m / p**2 * (2 * p * x - x**2),
                      m / (1 - p)**2 * ((1 - 2 * p) + 2 * p * x - x**2))
        dyc = np.where(x < p, 2 * m / p**2 * (p - x),
                       2 * m / (1 - p)**2 * (p - x))
    else:
        yc = np.zeros_like(x)
        dyc = np.zeros_like(x)
    th = np.arctan(dyc)
    xu, yu = x - yt * np.sin(th), yc + yt * np.cos(th)
    xl, yl = x + yt * np.sin(th), yc - yt * np.cos(th)
    # trailing edge -> upper surface -> leading edge -> lower surface (counter-clockwise)
    X = np.concatenate((xu[::-1], xl[1:-1]))
    Y = np.concatenate((yu[::-1], yl[1:-1]))
    return np.column_stack((X - 0.5, Y)) * chord


def clean_points(pts):
    keep = [0]
    for k in range(1, len(pts)):
        if np.hypot(*(pts[k] - pts[keep[-1]])) > 1e-12:
            keep.append(k)
    pts = pts[keep]
    if np.hypot(*(pts[0] - pts[-1])) < 1e-12:
        pts = pts[:-1]
    return pts


def rotate(pts, aoa_deg, c):
    a = np.radians(aoa_deg)
    d = pts - c
    x = c[0] + d[:, 0] * np.cos(a) + d[:, 1] * np.sin(a)
    y = c[1] - d[:, 0] * np.sin(a) + d[:, 1] * np.cos(a)
    return np.column_stack((x, y))


def make_body(args):
    if args.geometry:
        pts = load_points(args.geometry)
        c = np.array(args.center) if args.center else 0.5 * (pts.min(0) + pts.max(0))
    else:
        c = np.array(args.center if args.center else (0.0, 0.0))
        if args.shape == "circle":
            pts = shape_circle(args.radius, args.npts)
        elif args.shape == "ellipse":
            pts = shape_ellipse(0.5 * args.length, 0.5 * args.thickness, args.npts)
        elif args.shape == "rectangle":
            pts = shape_rectangle(args.length, args.thickness)
        elif args.shape == "naca":
            pts = shape_naca4(args.naca, args.length, args.npts // 2)
        else:
            sys.exit("error: give --geometry FILE or --shape NAME")
        pts = pts + c
    pts = clean_points(pts)
    if abs(args.aoa) > 0:
        pts = rotate(pts, args.aoa, c)
    poly = Polygon(pts)
    if not poly.is_valid:
        print("warning: body polygon is not valid (self-intersections); applying buffer(0)")
        poly = poly.buffer(0)
        if poly.geom_type != "Polygon":
            sys.exit("error: body polygon could not be repaired")
    return poly


# ----------------------------------------------------------------------------
#  Mesh construction
# ----------------------------------------------------------------------------

def build_mesh(body, xmin, xmax, ymin, ymax, h, sliver_tol):
    bx0, by0, bx1, by1 = body.bounds
    if bx0 <= xmin or bx1 >= xmax or by0 <= ymin or by1 >= ymax:
        sys.exit("error: the body must lie strictly inside the domain")

    nx = max(1, int(round((xmax - xmin) / h)))
    ny = max(1, int(round((ymax - ymin) / h)))
    xs = np.linspace(xmin, xmax, nx + 1)
    ys = np.linspace(ymin, ymax, ny + 1)
    hx = xs[1] - xs[0]
    hy = ys[1] - ys[0]
    cell_area = hx * hy
    tol = 1e-9 * min(hx, hy)

    fluid_quad = np.ones((nx, ny), dtype=bool)
    cut_list = []

    # only cells near the body need the geometric tests
    i0 = max(int(np.floor((bx0 - xmin) / hx)) - 1, 0)
    i1 = min(int(np.ceil((bx1 - xmin) / hx)) + 1, nx)
    j0 = max(int(np.floor((by0 - ymin) / hy)) - 1, 0)
    j1 = min(int(np.ceil((by1 - ymin) / hy)) + 1, ny)

    ii, jj = np.meshgrid(np.arange(i0, i1), np.arange(j0, j1), indexing="ij")
    ii = ii.ravel()
    jj = jj.ravel()
    boxes = shapely.box(xs[ii], ys[jj], xs[ii + 1], ys[jj + 1])
    shapely.prepare(body)
    hit = shapely.intersects(boxes, body)
    inside = shapely.contains(body, boxes)

    fluid_quad[ii[inside], jj[inside]] = False

    # extra nodes created in cut cells
    new_pts = []          # coordinates
    cut_tris = []         # (3 entries, each ('g', gridnode) or ('n', newindex))

    def gid(i, j):
        return i * (ny + 1) + j

    n_cut = 0
    n_touch = 0
    for k in np.nonzero(hit & ~inside)[0]:
        i, j = int(ii[k]), int(jj[k])
        box = boxes[k]
        fl = box.difference(body)
        parts = [g for g in getattr(fl, "geoms", [fl])
                 if g.geom_type == "Polygon" and g.area > sliver_tol * cell_area]
        if not parts:
            fluid_quad[i, j] = False
            continue
        if len(parts) == 1 and abs(parts[0].area - cell_area) <= 1e-12 * cell_area \
                and len(parts[0].interiors) == 0:
            n_touch += 1                     # body only touches the cell: keep the quad
            continue
        fluid_quad[i, j] = False
        n_cut += 1
        corners = np.array([[xs[i], ys[j]], [xs[i + 1], ys[j]],
                            [xs[i + 1], ys[j + 1]], [xs[i], ys[j + 1]]])
        cids = [gid(i, j), gid(i + 1, j), gid(i + 1, j + 1), gid(i, j + 1)]
        for part in parts:
            tris = shapely.constrained_delaunay_triangles(part)
            for tri in tris.geoms:
                c = np.array(tri.exterior.coords)[:3]
                entry = []
                for p in c:
                    d = np.hypot(corners[:, 0] - p[0], corners[:, 1] - p[1])
                    m = int(np.argmin(d))
                    if d[m] < tol:
                        entry.append(("g", cids[m]))
                    else:
                        new_pts.append(p)
                        entry.append(("n", len(new_pts) - 1))
                cut_tris.append(entry)

    # ---- nodes: grid nodes + merged new nodes -----------------------------
    GX, GY = np.meshgrid(xs, ys, indexing="ij")
    grid_xy = np.column_stack((GX.ravel(), GY.ravel()))
    n_grid = len(grid_xy)

    if new_pts:
        new_pts = np.array(new_pts)
        tree = cKDTree(new_pts)
        pairs = tree.query_pairs(tol, output_type="ndarray")
        npn = len(new_pts)
        if len(pairs):
            A = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(npn, npn))
            _, lab = connected_components(A, directed=False)
        else:
            lab = np.arange(npn)
        nlab = lab.max() + 1
        rep = np.zeros((nlab, 2))
        cnt = np.zeros(nlab)
        np.add.at(rep, lab, new_pts)
        np.add.at(cnt, lab, 1)
        rep /= cnt[:, None]
    else:
        lab = np.zeros(0, dtype=int)
        rep = np.zeros((0, 2))

    all_xy = np.vstack((grid_xy, rep))

    def resolve(entry):
        return entry[1] if entry[0] == "g" else n_grid + int(lab[entry[1]])

    tris = np.array([[resolve(e) for e in t] for t in cut_tris], dtype=int).reshape(-1, 3)

    qi, qj = np.nonzero(fluid_quad)
    quads = np.column_stack((qi * (ny + 1) + qj, (qi + 1) * (ny + 1) + qj,
                             (qi + 1) * (ny + 1) + qj + 1, qi * (ny + 1) + qj + 1))

    # counter-clockwise triangles
    if len(tris):
        P = all_xy[tris]
        a = 0.5 * ((P[:, 1, 0] - P[:, 0, 0]) * (P[:, 2, 1] - P[:, 0, 1])
                   - (P[:, 2, 0] - P[:, 0, 0]) * (P[:, 1, 1] - P[:, 0, 1]))
        flip = a < 0
        tris[flip] = tris[flip][:, [0, 2, 1]]

    # drop unused nodes
    used = np.zeros(len(all_xy), dtype=bool)
    used[quads.ravel()] = True
    used[tris.ravel()] = True
    remap = -np.ones(len(all_xy), dtype=int)
    remap[used] = np.arange(used.sum())
    nodes = all_xy[used]
    quads = remap[quads]
    tris = remap[tris]

    info = dict(nx=nx, ny=ny, hx=hx, hy=hy, n_cut=n_cut, n_touch=n_touch)
    return nodes, tris, quads, info


def boundary_edges(nodes, tris, quads, xmin, xmax, ymin, ymax, tol):
    e = []
    if len(tris):
        e.append(np.vstack((tris[:, [0, 1]], tris[:, [1, 2]], tris[:, [2, 0]])))
    if len(quads):
        e.append(np.vstack((quads[:, [0, 1]], quads[:, [1, 2]], quads[:, [2, 3]], quads[:, [3, 0]])))
    e = np.vstack(e)
    key = np.sort(e, axis=1)
    _, inv, cnt = np.unique(key, axis=0, return_inverse=True, return_counts=True)
    inv = inv.ravel()
    bmask = cnt[inv] == 1
    be = e[bmask]                          # directed, fluid on the left
    mid = 0.5 * (nodes[be[:, 0]] + nodes[be[:, 1]])
    on_dom = ((np.abs(mid[:, 0] - xmin) < tol) | (np.abs(mid[:, 0] - xmax) < tol)
              | (np.abs(mid[:, 1] - ymin) < tol) | (np.abs(mid[:, 1] - ymax) < tol))
    tag = np.where(on_dom, 2, 1)
    return be, tag


# ----------------------------------------------------------------------------
#  Output
# ----------------------------------------------------------------------------

def write_msh(fname, nodes, tris, quads, be, tag):
    with open(fname, "w") as f:
        f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        f.write('$PhysicalNames\n3\n1 1 "body"\n1 2 "farfield"\n2 3 "fluid"\n$EndPhysicalNames\n')
        f.write("$Nodes\n%d\n" % len(nodes))
        for k, p in enumerate(nodes, 1):
            f.write("%d %.16g %.16g 0\n" % (k, p[0], p[1]))
        f.write("$EndNodes\n")
        f.write("$Elements\n%d\n" % (len(be) + len(tris) + len(quads)))
        eid = 1
        for (a, b), t in zip(be, tag):
            f.write("%d 1 2 %d %d %d %d\n" % (eid, t, t, a + 1, b + 1))
            eid += 1
        for t in tris:
            f.write("%d 2 2 3 3 %d %d %d\n" % (eid, t[0] + 1, t[1] + 1, t[2] + 1))
            eid += 1
        for q in quads:
            f.write("%d 3 2 3 3 %d %d %d %d\n" % (eid, q[0] + 1, q[1] + 1, q[2] + 1, q[3] + 1))
            eid += 1
        f.write("$EndElements\n")


def write_txt(base, nodes, tris, quads, be, tag):
    np.savetxt(base + "_nodes.txt", nodes, fmt="%.16g")
    np.savetxt(base + "_tri.txt", tris + 1, fmt="%d")
    np.savetxt(base + "_quad.txt", quads + 1, fmt="%d")
    np.savetxt(base + "_bnd.txt", np.column_stack((be + 1, tag)), fmt="%d")


def plot_mesh(fname, nodes, tris, quads, body):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection

    fig, ax = plt.subplots(1, 2, figsize=(13, 6))
    polys = [nodes[c] for c in quads] + [nodes[c] for c in tris]
    cols = ["none"] * len(quads) + ["#f4c542"] * len(tris)
    bx0, by0, bx1, by1 = body.bounds
    w = max(bx1 - bx0, by1 - by0)
    for a, lim in zip(ax, [None, (bx0 - 0.6 * w, bx1 + 0.6 * w, by0 - 0.6 * w, by1 + 0.6 * w)]):
        a.add_collection(PolyCollection(polys, facecolors=cols, edgecolors="k", linewidths=0.3))
        x, y = body.exterior.xy
        a.plot(x, y, "r-", lw=1.2)
        a.set_aspect("equal")
        if lim is None:
            a.set_xlim(nodes[:, 0].min(), nodes[:, 0].max())
            a.set_ylim(nodes[:, 1].min(), nodes[:, 1].max())
            a.set_title("full domain")
        else:
            a.set_xlim(lim[0], lim[1])
            a.set_ylim(lim[2], lim[3])
            a.set_title("near the body (triangles in yellow)")
    fig.tight_layout()
    fig.savefig(fname, dpi=130)


# ----------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="Cartesian + cut-cell triangle mesher for a 2D geometry")
    ap.add_argument("--geometry", help="text file with 'x y' rows (closed or open polygon)")
    ap.add_argument("--shape", choices=["circle", "ellipse", "rectangle", "naca"])
    ap.add_argument("--radius", type=float, default=0.5, help="circle radius")
    ap.add_argument("--length", type=float, default=1.0, help="ellipse/rectangle length, NACA chord")
    ap.add_argument("--thickness", type=float, default=0.7, help="ellipse/rectangle thickness")
    ap.add_argument("--naca", default="0012", help="NACA 4-digit code")
    ap.add_argument("--npts", type=int, default=400, help="points on curved built-in shapes")
    ap.add_argument("--center", type=float, nargs=2, help="body centre / rotation pivot")
    ap.add_argument("--aoa", type=float, default=0.0, help="rotate the body, deg (positive = nose up)")
    ap.add_argument("--xmin", type=float, default=-10.0)
    ap.add_argument("--xmax", type=float, default=20.0)
    ap.add_argument("--ymin", type=float, default=-10.0)
    ap.add_argument("--ymax", type=float, default=10.0)
    ap.add_argument("--h", type=float, default=0.1, help="Cartesian cell size")
    ap.add_argument("--sliver-tol", type=float, default=1e-10,
                    help="fluid pieces smaller than this fraction of a cell are dropped")
    ap.add_argument("--out", default="mesh", help="output basename")
    ap.add_argument("--plot", action="store_true", help="write <out>.png")
    args = ap.parse_args()

    if not args.geometry and not args.shape:
        ap.error("give --geometry FILE or --shape NAME")

    body = make_body(args)
    print("body: %d vertices, area %.6g, bounds %s" % (len(body.exterior.coords) - 1, body.area,
                                                       tuple(round(v, 4) for v in body.bounds)))

    nodes, tris, quads, info = build_mesh(body, args.xmin, args.xmax, args.ymin, args.ymax,
                                          args.h, args.sliver_tol)
    tol = 1e-9 * min(info["hx"], info["hy"])
    be, tag = boundary_edges(nodes, tris, quads, args.xmin, args.xmax, args.ymin, args.ymax,
                             max(tol, 1e-9))

    write_msh(args.out + ".msh", nodes, tris, quads, be, tag)
    write_txt(args.out, nodes, tris, quads, be, tag)
    if args.plot:
        plot_mesh(args.out + ".png", nodes, tris, quads, body)

    print("grid: %d x %d, hx = %.5g, hy = %.5g" % (info["nx"], info["ny"], info["hx"], info["hy"]))
    print("nodes: %d   quads: %d   triangles: %d   cut cells: %d" % (len(nodes), len(quads), len(tris), info["n_cut"]))
    print("boundary edges: body %d, farfield %d" % ((tag == 1).sum(), (tag == 2).sum()))
    print("written: %s.msh, %s_nodes.txt, %s_tri.txt, %s_quad.txt, %s_bnd.txt" % ((args.out,) * 5))
    print("next: python3 inspect_mesh.py %s.msh" % args.out)


if __name__ == "__main__":
    main()
