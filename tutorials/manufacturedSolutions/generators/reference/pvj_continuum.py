"""Continuum Purkinje-myocardium junction on a mesh: NumPy only.

A terminal at x_t exchanges with the tissue through a radial kernel of support radius R,

    W      = int_{Omega cap B} w dx                  (normalisation)
    <f>    = int_{Omega cap B} w f dx / W            (gathered value of a field f)
    S(x)   = I w(|x - x_t|) / W                      (deposited source, current I)

over the part of the ball B = {|x - x_t| <= R} inside the true tissue Omega.  Kernels
(w(r) = phi(r/R) for r <= R, zero beyond):

    uniform   1                       linear    1 - s
    gaussian  exp(-9 s^2/2)           cosine    (1 + cos(pi s))/2

This module evaluates the continuum quantities of a terminal and their integrals over the
cells of a mesh, with the algorithm of the cardiacFoam class pvjContinuumReference
(src/verificationModels/coupledVerification/pvjContinuumReference.{H,C}), so that a
solver's discrete junction map can be tested against the continuum coupling.  It needs
NumPy only.

Functions
---------
support_integral          W and the numerator of <f>: int_{Omega cap B} w f dx for a ball
                          cut by nothing (interior), by a flat wall (plane) or by a sphere
                          with the tissue outside it (cavitySphere)
cell_integrals            int_{K cap B} w dx and int_{K cap B} w f dx for each mesh cell K
junction_source           the cell average I omega_K / (W |K|) of the continuum source
cell_geometry, cell_volumes, outward_faces
                          cell centre, volume and face centres as OpenFOAM computes them; the
                          faces of a cell with inward-pointing faces reversed
cells_from_hexahedra, cells_from_tetrahedra, cells_from_polymesh
                          build the cell input from vertex coordinates or a polyhedral mesh
box_support, shell_support
                          the support kind of a terminal on a box or in a sealed spherical
                          shell, as the cardiacFoam verifier chooses it

Mesh input
----------
A mesh is a list `cells`; a cell is a list of faces; a face is an array (n_vertices, 3) of
vertex coordinates in cyclic order.  Cells may be any polyhedra whose faces are planar or
nearly so: hexahedra, tetrahedra, prisms, general polyhedra.  Nothing else is needed:
connectivity is not used, so a vertex shared by two cells is simply given twice.  The helpers
build this input

    cells_from_hexahedra(v)    v (n, 8, 3), vertices in OpenFOAM/blockMesh hex order
                               (0-3 the bottom quadrilateral, 4-7 the top one above it)
    cells_from_tetrahedra(v)   v (n, 4, 3), any vertex order
    cells_from_polymesh(points, faces, owner, neighbour)
                               a polyMesh: points (n, 3), faces lists of vertex indices,
                               owner and neighbour cell of each face (neighbour covers the
                               internal faces only)

The vertex order of a face fixes the order of the vertices of the tetrahedra
(cell centre, face centre, p_i, p_i+1) it generates, and the collapsed Gauss rule is not
symmetric in the vertex order of a tetrahedron.  Reversing a face therefore changes the
integrals of a cell that the sphere cuts at the level of the quadrature error (up to about 1e-7
of a weakly cut cell and 2e-9 of the sum at cut size R/8, 1e-11 at R/16 in the tests of this
module); a wholly inside cell changes at rounding level.  cells_from_polymesh keeps each stored
face as it is in both cells it bounds, which is what OpenFOAM does, so on a polyMesh this module
reproduces pvjContinuumReference to rounding (1e-14).  cells_from_hexahedra and
cells_from_tetrahedra orient every face out of its cell.  The orientation does not enter the
cell geometry.

A tetrahedron whose sub-tetrahedra have exactly equal octahedron diagonals (the Kuhn
tetrahedra of an unperturbed grid, for one) is split by the last bit of a distance, so two
implementations, or two compilers, may split it differently and then differ by the quadrature
error of the cut tetrahedra (1e-6 of a weakly cut cell in the tests).  Meshes of real
geometry do not tie.

The tissue boundary is carried by the mesh: cell_integrals integrates over the cells given,
so the part of the ball outside the mesh contributes nothing.  The support kind describes
the TRUE tissue and only enters W: the mesh of a curved wall is polyhedral, and the
difference between the two is what a junction verification measures.

Cell geometry follows OpenFOAM: the face centre is the area-weighted centroid of the
triangles (face centroid estimate, edge) and the cell centre the volume-weighted centroid of
the pyramids over the faces about the mean of the face centres.  Each cell is split into the
tetrahedra (cell centre, face centre, edge) and each tetrahedron is integrated as below.

Algorithm of cell_integrals (every number is a parameter of pvjContinuumReference and is
documented there; the sizes are fractions of R)
  a tetrahedron that cannot meet the ball is dropped;
  one inside the ball is integrated whole with a collapsed (Duffy) Gauss rule of 4 x 4 x 4
  points once its longest edge is at most inside_size = 1/4 (halved otherwise), and, for the
  linear kernel, a tetrahedron that contains the kernel's cone tip (the terminal) is halved
  until its longest edge is at most tip_size = 1/64;
  one cut by the sphere |x - x_t| = R is halved (four corner tetrahedra and the octahedron
  between them, split about its shortest diagonal) until its longest edge is at most
  cut_size (1/96 for uniform and gaussian, which jump at R, 1/48 for linear, 1/24 for
  cosine), then clipped by the linear interpolant of |x - x_t| - R and integrated with a
  2 x 2 x 2 collapsed Gauss rule on the pieces, the kernel formula continued past R.
support_integral integrates in spherical coordinates about x_t: Gauss-Legendre with 48 points
in r and in the polar cosine mu, the periodic trapezoid rule with 96 points in phi.

Conventions: points are arrays (N, 3); a field f is a function of such an array returning an
array (N,).  Lengths are in the units of the mesh.

Not covered (neither by this module nor by the C++ class): a support cut by two walls, by a
sphere with the tissue inside it (a terminal on the outer sphere of a shell) or by a curved wall
that is not a sphere; the mesh must hold the whole of Omega cap B.
Overlapping supports are superposed by summing the sources of the terminals.
"""

import collections

import numpy

# module: pvj_continuum
# Generated by coupled_1d3d_mms.py from mms_common/pvj_continuum.py

KERNELS = ("uniform", "linear", "gaussian", "cosine")
SUPPORT_KINDS = ("interior", "plane", "cavitySphere")

# Defaults of pvjContinuumReference (referenceRadialPoints, referencePolarPoints,
# referenceAzimuthalPoints, referenceCutSize, referenceInsideSize, referenceTipSize)
N_RADIAL = 48
N_POLAR = 48
N_AZIMUTHAL = 96
INSIDE_SIZE = 0.25
TIP_SIZE = 1.0/64.0
CUT_SIZE = {"cosine": 1.0/24.0, "linear": 1.0/48.0, "uniform": 1.0/96.0, "gaussian": 1.0/96.0}

VSMALL = 1e-300
ROOTVSMALL = 1e-150

CellJunction = collections.namedtuple(
    "CellJunction", "cells omega omega_f volumes normalisation source")


def _norm(v):
    """Euclidean length along the last axis."""
    return numpy.sqrt(numpy.einsum("...i,...i->...", v, v))


def _check_kernel(kernel):
    if kernel not in KERNELS:
        raise ValueError("no continuum definition of kernel %r; defined: %s"
                         % (kernel, " ".join(KERNELS)))


def kernel_profile(kernel, s):
    """phi(s), s = r/R: the kernel formula, continued past s = 1."""
    _check_kernel(kernel)
    s = numpy.asarray(s, dtype=float)
    if kernel == "linear":
        return 1.0 - s
    if kernel == "gaussian":
        return numpy.exp(-4.5*s*s)
    if kernel == "cosine":
        return 0.5*(1.0 + numpy.cos(numpy.pi*s))
    return numpy.ones_like(s)


def kernel_weight(kernel, r, radius):
    """w(r): the kernel inside the support radius, zero beyond it."""
    r = numpy.asarray(r, dtype=float)
    return numpy.where(r <= radius, kernel_profile(kernel, r/radius), 0.0)


def _unit_gauss_rule(n):
    """Gauss-Legendre nodes and weights on [0, 1]."""
    x, w = numpy.polynomial.legendre.leggauss(n)
    return 0.5*(x + 1.0), 0.5*w


# ---- support integral ---------------------------------------------------------------

def support_integral(centre, radius, kernel, f=None, kind="interior", normal=None,
                     wall_radius=None, n_radial=N_RADIAL, n_polar=N_POLAR,
                     n_azimuthal=N_AZIMUTHAL):
    """int_{Omega cap B} w f dx, f = 1 when None.

    kind     interior      B lies inside Omega
             plane         the terminal is on a flat wall, Omega the half space on the side
                           of `normal` (the unit normal into the tissue): mu >= 0
             cavitySphere  the terminal is on a sphere of radius wall_radius with the tissue
                           outside it, `normal` the outward radius at the terminal:
                           mu >= -r/(2 wall_radius), with mu the cosine of the polar angle
                           about `normal`
    """
    _check_kernel(kernel)
    if kind not in SUPPORT_KINDS:
        raise ValueError("support kind %r not in %s" % (kind, ", ".join(SUPPORT_KINDS)))
    centre = numpy.asarray(centre, dtype=float)
    e3 = numpy.array([0.0, 0.0, 1.0])
    if kind != "interior":
        if normal is None:
            raise ValueError("support kind %s needs a normal" % kind)
        normal = numpy.asarray(normal, dtype=float)
        e3 = normal/numpy.linalg.norm(normal)
    if kind == "cavitySphere" and not (wall_radius and wall_radius > 0):
        raise ValueError("support kind cavitySphere needs a positive wall_radius")
    trial = numpy.array([1.0, 0.0, 0.0]) if abs(e3[0]) < 0.9 else numpy.array([0.0, 1.0, 0.0])
    e1 = trial - (trial @ e3)*e3
    e1 /= numpy.linalg.norm(e1)
    e2 = numpy.cross(e3, e1)

    xr, wr = _unit_gauss_rule(n_radial)
    xm, wm = _unit_gauss_rule(n_polar)
    d_phi = 2.0*numpy.pi/n_azimuthal
    phi = (numpy.arange(n_azimuthal) + 0.5)*d_phi

    r = radius*xr
    radial_weight = radius*wr*r*r*kernel_weight(kernel, r, radius)          # (nr,)
    if kind == "plane":
        mu_min = numpy.zeros_like(r)
    elif kind == "cavitySphere":
        mu_min = -r/(2.0*wall_radius)
    else:
        mu_min = -numpy.ones_like(r)
    mu_span = 1.0 - mu_min
    mu = mu_min[:, None] + mu_span[:, None]*xm[None, :]                       # (nr, nm)
    sin_theta = numpy.sqrt(numpy.maximum(1.0 - mu*mu, 0.0))

    cos_phi, sin_phi = numpy.cos(phi), numpy.sin(phi)
    direction = (mu[..., None, None]*e3
                 + sin_theta[..., None, None]
                 *(cos_phi[:, None]*e1 + sin_phi[:, None]*e2))                # (nr, nm, na, 3)
    points = centre + r[:, None, None, None]*direction
    if f is None:
        sum_phi = numpy.full(mu.shape, float(n_azimuthal))
    else:
        values = numpy.asarray(f(points.reshape(-1, 3)), dtype=float).reshape(points.shape[:-1])
        sum_phi = values.sum(axis=-1)
    sum_mu = (mu_span[:, None]*wm[None, :]*d_phi*sum_phi).sum(axis=1)
    return float(numpy.sum(radial_weight*sum_mu))


# ---- cell geometry --------------------------------------------------------------------

def _face_centre_and_area(p):
    """Centre and area vector of a face (n, 3), as primitiveMeshTools::makeFaceCentresAndAreas."""
    n = len(p)
    if n == 3:
        return (1.0/3.0)*(p[0] + p[1] + p[2]), 0.5*numpy.cross(p[1] - p[0], p[2] - p[0])
    estimate = p.sum(axis=0)/n
    following = numpy.roll(p, -1, axis=0)
    cross = numpy.cross(following - p, estimate - p)
    area = numpy.linalg.norm(cross, axis=1)
    sum_a = area.sum()
    if sum_a < ROOTVSMALL:
        return estimate, numpy.zeros(3)
    centre = (1.0/3.0)*((p + following + estimate)*area[:, None]).sum(axis=0)/sum_a
    return centre, 0.5*cross.sum(axis=0)


def cell_geometry(faces):
    """(cell centre, cell volume, face centres (nf, 3)) of a cell, as
    primitiveMeshTools::makeCellCentresAndVols.  Each face pyramid counts with the sign of an
    outward face (OpenFOAM flips the face for the neighbour), which assumes the mean of the
    face centres lies inside the cell, as it does for any cell that is not badly distorted."""
    pairs = [_face_centre_and_area(numpy.asarray(f, dtype=float)) for f in faces]
    centres = numpy.array([c for c, _ in pairs])
    areas = numpy.array([a for _, a in pairs])
    estimate = centres.mean(axis=0)
    pyramid3 = numpy.abs(numpy.einsum("ij,ij->i", areas, centres - estimate))
    pyramid_centre = 0.75*centres + 0.25*estimate
    total = pyramid3.sum()
    if abs(total) > VSMALL:
        centre = (pyramid3[:, None]*pyramid_centre).sum(axis=0)/total
    else:
        centre = estimate
    return centre, total/3.0, centres


def cell_volumes(cells):
    """Volume |K| of every cell, as OpenFOAM computes it (mesh.V())."""
    return numpy.array([cell_geometry(c)[1] for c in cells])


def outward_faces(faces):
    """The faces of a cell with their vertex cycles reversed where the area vector points into
    the cell (relative to the mean of the face centres)."""
    faces = [numpy.asarray(f, dtype=float) for f in faces]
    pairs = [_face_centre_and_area(f) for f in faces]
    estimate = numpy.mean([c for c, _ in pairs], axis=0)
    return [f if a @ (c - estimate) >= 0.0 else f[::-1] for f, (c, a) in zip(faces, pairs)]


_HEX_FACES = ((0, 4, 7, 3), (1, 2, 6, 5), (0, 1, 5, 4), (3, 7, 6, 2), (0, 3, 2, 1), (4, 5, 6, 7))
_TET_FACES = ((1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1))


def cells_from_hexahedra(vertices):
    """Cells of hexahedra v (n, 8, 3) in OpenFOAM vertex order."""
    vertices = numpy.asarray(vertices, dtype=float)
    return [outward_faces([v[list(f)] for f in _HEX_FACES]) for v in vertices]


def cells_from_tetrahedra(vertices):
    """Cells of tetrahedra v (n, 4, 3), any vertex order."""
    vertices = numpy.asarray(vertices, dtype=float)
    return [outward_faces([v[list(f)] for f in _TET_FACES]) for v in vertices]


def cells_from_polymesh(points, faces, owner, neighbour):
    """Cells of a polyhedral mesh: face i joins cell owner[i] to neighbour[i] (internal
    faces, i < len(neighbour)) or bounds cell owner[i].  Each face keeps its stored vertex
    order in both cells, as OpenFOAM's pvjContinuumReference iterates it."""
    points = numpy.asarray(points, dtype=float)
    n_cells = max(max(owner), max(neighbour) if len(neighbour) else -1) + 1
    cells = [[] for _ in range(n_cells)]
    for i, face in enumerate(faces):
        cells[owner[i]].append(points[list(face)])
        if i < len(neighbour):
            cells[neighbour[i]].append(points[list(face)])
    return cells


# ---- cell integrals -------------------------------------------------------------------

def _duffy_rule(n):
    """Barycentric offsets (l1, l2, l3), shape (n^3, 3), of the collapsed Gauss rule of n^3
    points on the unit simplex and its weights: pvjContinuumReference::quadratureTet."""
    x, w = _unit_gauss_rule(n)
    u = x[:, None, None]
    v = x[None, :, None]
    k = x[None, None, :]
    l1 = numpy.broadcast_to(u, (n, n, n))
    l2 = numpy.broadcast_to((1.0 - u)*v, (n, n, n))
    l3 = (1.0 - u)*(1.0 - v)*k
    jac = (w[:, None, None]*w[None, :, None]*w[None, None, :]*(1.0 - u)**2*(1.0 - v))
    coefficients = numpy.stack([a.reshape(-1) for a in (l1, l2, l3)], axis=1)      # (n^3, 3)
    return coefficients, jac.reshape(-1)


_RULES = {4: _duffy_rule(4), 2: _duffy_rule(2)}
_CHUNK = 20000
_BLOCK = 8192


class _Integrator:
    """Integrals of w and w f over tetrahedra, accumulated per owning cell."""

    def __init__(self, terminal, radius, kernel, f, n_cells, cut_size, inside_size, tip_size):
        self.t = numpy.asarray(terminal, dtype=float)
        self.R = float(radius)
        self.kernel = kernel
        self.f = f
        self.cut = self.R*cut_size
        self.inside = self.R*inside_size
        self.tip = self.R*tip_size
        self.omega = numpy.zeros(n_cells)
        self.omega_f = numpy.zeros(n_cells)

    def _distance(self, p):
        return _norm(p - self.t)

    def _quadrature(self, tets, owner, sign, clipped):
        """Collapsed Gauss rule on each tet (N, 4, 3), the kernel extended past R."""
        coefficients, jac = _RULES[2 if clipped else 4]
        for lo in range(0, len(tets), _CHUNK):
            part = tets[lo:lo + _CHUNK]
            who = owner[lo:lo + _CHUNK]
            a = part[:, 0]
            edges = part[:, 1:] - a[:, None, :]                       # ab, ac, ad
            det = numpy.abs(numpy.einsum("ni,ni->n", edges[:, 0], numpy.cross(edges[:, 1], edges[:, 2])))
            det = numpy.where(det > VSMALL, det, 0.0)
            p = a[:, None, :] + numpy.matmul(coefficients, edges)         # (N, n^3, 3)
            wp = jac[None, :]*kernel_profile(self.kernel, self._distance(p)/self.R)
            self.omega += numpy.bincount(who, sign*det*wp.sum(axis=1), len(self.omega))
            if self.f is not None:
                values = numpy.asarray(self.f(p.reshape(-1, 3)), dtype=float).reshape(wp.shape)
                self.omega_f += numpy.bincount(
                    who, sign*det*(wp*values).sum(axis=1), len(self.omega_f))

    @staticmethod
    def _zero_on_edge(p, q, phi_p, phi_q):
        s = (phi_p/(phi_p - phi_q))[:, None]
        return p + s*(q - p)

    def _clip(self, tets, owner):
        """The part of tets (all cut by the sphere) inside it, as pvjContinuumReference::clipTet."""
        phi = self._distance(tets) - self.R
        inside = phi < 0
        n_in = inside.sum(axis=1)

        whole = n_in == 4
        if whole.any():
            self._quadrature(tets[whole], owner[whole], 1.0, True)

        others = numpy.array([[i for i in range(4) if i != lone] for lone in range(4)])
        for count in (1, 3):
            sel = n_in == count
            if not sel.any():
                continue
            vt, ph, ow = tets[sel], phi[sel], owner[sel]
            # the vertex alone on its side of the cut: the last one on that side
            on_side = (ph < 0) == (count == 1)
            lone = 3 - numpy.argmax(on_side[:, ::-1], axis=1)
            r = numpy.arange(len(vt))
            others_of = others[lone]                                            # (m, 3)
            cut = [self._zero_on_edge(vt[r, lone], vt[r, others_of[:, i]],
                                      ph[r, lone], ph[r, others_of[:, i]]) for i in range(3)]
            corner = numpy.stack([vt[r, lone], *cut], axis=1)
            if count == 1:
                self._quadrature(corner, ow, 1.0, True)
            else:
                self._quadrature(vt, ow, 1.0, True)
                self._quadrature(corner, ow, -1.0, True)

        sel = n_in == 2
        if sel.any():
            vt, ph, ow = tets[sel], phi[sel], owner[sel]
            order = numpy.argsort(ph >= 0, axis=1, kind="stable")                  # inside first
            r = numpy.arange(len(vt))
            i0, i1, o0, o1 = (order[:, k] for k in range(4))
            a0, b0 = vt[r, i0], vt[r, i1]
            a1 = self._zero_on_edge(a0, vt[r, o0], ph[r, i0], ph[r, o0])
            a2 = self._zero_on_edge(a0, vt[r, o1], ph[r, i0], ph[r, o1])
            b1 = self._zero_on_edge(b0, vt[r, o0], ph[r, i1], ph[r, o0])
            b2 = self._zero_on_edge(b0, vt[r, o1], ph[r, i1], ph[r, o1])
            for corners in ((a0, a1, a2, b2), (a0, a1, b1, b2), (a0, b0, b1, b2)):
                self._quadrature(numpy.stack(corners, axis=1), ow, 1.0, True)

    def run(self, tets, owner):
        """Integrate the tetrahedra: drop, integrate or clip, or halve and go on with the
        children.  Blocks of _BLOCK tetrahedra keep the arrays small."""
        linear = self.kernel == "linear"
        stack = [(tets[i:i + _BLOCK], owner[i:i + _BLOCK]) for i in range(0, len(tets), _BLOCK)]
        while stack:
            tets, owner = stack.pop()
            a, b, c, d = tets[:, 0], tets[:, 1], tets[:, 2], tets[:, 3]
            g = 0.25*(a + b + c + d)
            rad = numpy.maximum(numpy.maximum(_norm(a - g), _norm(b - g)),
                             numpy.maximum(_norm(c - g), _norm(d - g)))
            dg = self._distance(g)
            keep = ~(dg - rad >= self.R)
            tets, owner, rad, dg = tets[keep], owner[keep], rad[keep], dg[keep]
            if not len(tets):
                continue
            a, b, c, d = tets[:, 0], tets[:, 1], tets[:, 2], tets[:, 3]
            d_max = numpy.max(self._distance(tets), axis=1)
            size = numpy.maximum(numpy.maximum(numpy.maximum(_norm(b - a), _norm(c - a)),
                                         numpy.maximum(_norm(d - a), _norm(c - b))),
                              numpy.maximum(_norm(d - b), _norm(d - c)))
            inside = d_max <= self.R
            split = numpy.where(inside,
                             (size > self.inside) | (linear & (dg < rad) & (size > self.tip)),
                             size > self.cut)
            whole = ~split & inside
            if whole.any():
                self._quadrature(tets[whole], owner[whole], 1.0, False)
            clipped = ~split & ~inside
            if clipped.any():
                self._clip(tets[clipped], owner[clipped])
            if split.any():
                children, who = _halve(tets[split]), numpy.repeat(owner[split], 8)
                stack += [(children[i:i + _BLOCK], who[i:i + _BLOCK])
                          for i in range(0, len(children), _BLOCK)]


def _halve(tets):
    """The 8 children of each tet (N, 4, 3): four corners and the octahedron split about its
    shortest diagonal, as pvjContinuumReference::integrateTet."""
    a, b, c, d = tets[:, 0], tets[:, 1], tets[:, 2], tets[:, 3]
    m01, m02, m03 = 0.5*(a + b), 0.5*(a + c), 0.5*(a + d)
    m12, m13, m23 = 0.5*(b + c), 0.5*(b + d), 0.5*(c + d)
    corners = [numpy.stack(t, axis=1) for t in
               ((a, m01, m02, m03), (b, m01, m12, m13), (c, m02, m12, m23), (d, m03, m13, m23))]
    d0 = _norm(m01 - m23)
    d1 = _norm(m02 - m13)
    d2 = _norm(m03 - m12)
    first = (d0 <= d1) & (d0 <= d2)
    second = ~first & (d2 < d1)
    pick = lambda x0, x1, x2: numpy.where(first[:, None], x0, numpy.where(second[:, None], x1, x2))
    p = pick(m01, m03, m02)
    q = pick(m23, m12, m13)
    ring = [pick(m02, m01, m01), pick(m03, m02, m03), pick(m13, m23, m23), pick(m12, m13, m12)]
    octa = [numpy.stack((p, q, ring[k], ring[(k + 1) % 4]), axis=1) for k in range(4)]
    children = numpy.stack(corners + octa, axis=1)                                # (N, 8, 4, 3)
    return children.reshape(-1, 4, 3)


def _cell_tets(faces, centre, face_centres):
    """Tetrahedra (cell centre, face centre, edge) of one cell, shape (m, 4, 3)."""
    out = []
    for face, fc in zip(faces, face_centres):
        p = numpy.asarray(face, dtype=float)
        following = numpy.roll(p, -1, axis=0)
        n = len(p)
        out.append(numpy.stack([numpy.broadcast_to(centre, (n, 3)), numpy.broadcast_to(fc, (n, 3)),
                             p, following], axis=1))
    return numpy.concatenate(out)


def cell_integrals(cells, terminal, radius, kernel="linear", f=None, cut_size=None,
                   inside_size=INSIDE_SIZE, tip_size=TIP_SIZE):
    """int_{K cap B} w dx and int_{K cap B} w f dx for each cell K of `cells`.

    Returns (index, omega, omega_f): the positions in `cells` of the cells that meet the ball
    (omega or omega_f nonzero), in increasing order, and the two integrals over each; omega_f is
    zero when f is None.  cut_size, inside_size and tip_size are fractions of the radius
    (see the module text for the defaults).
    """
    _check_kernel(kernel)
    if radius <= 0:
        raise ValueError("support radius must be positive; found %r" % (radius,))
    terminal = numpy.asarray(terminal, dtype=float)
    if cut_size is None:
        cut_size = CUT_SIZE[kernel]
    tet_blocks, owners, candidates = [], [], []
    for i, faces in enumerate(cells):
        corners = numpy.concatenate([numpy.asarray(fc, dtype=float) for fc in faces])
        gap = numpy.maximum(numpy.maximum(corners.min(axis=0) - terminal, terminal - corners.max(axis=0)),
                         0.0)
        if numpy.linalg.norm(gap) >= radius:       # the ball does not reach the bounding box
            continue
        centre, _, face_centres = cell_geometry(faces)
        rad = numpy.max(numpy.linalg.norm(corners - centre, axis=1))
        if numpy.linalg.norm(centre - terminal) - rad >= radius:
            continue
        tets = _cell_tets(faces, centre, face_centres)
        tet_blocks.append(tets)
        owners.append(numpy.full(len(tets), len(candidates), dtype=numpy.intp))
        candidates.append(i)
    index = numpy.array(candidates, dtype=numpy.intp)
    if not candidates:
        return index, numpy.zeros(0), numpy.zeros(0)
    work = _Integrator(terminal, radius, kernel, f, len(candidates), cut_size, inside_size,
                       tip_size)
    work.run(numpy.concatenate(tet_blocks), numpy.concatenate(owners))
    hit = (work.omega != 0) | (work.omega_f != 0)
    return index[hit], work.omega[hit], work.omega_f[hit]


def junction_source(cells, terminal, radius, kernel="linear", support="interior", current=1.0,
                    f=None, normal=None, wall_radius=None, **options):
    """The continuum junction of one terminal on a mesh.

    Returns CellJunction(cells, omega, omega_f, volumes, normalisation, source): the cells
    that meet the ball, int_{K cap B} w and int_{K cap B} w f over each, their volumes |K|,
    W = support_integral over the true support (`support` is its kind, with `normal` and
    `wall_radius` for plane and cavitySphere) and the cell average of the continuum source
    S = I w/W, that is current*omega/(W*|K|).  The gathered value of f is
    sum(omega_f)/W when the mesh holds the whole support.
    """
    index, omega, omega_f = cell_integrals(cells, terminal, radius, kernel, f, **options)
    volumes = numpy.array([cell_geometry(cells[i])[1] for i in index])
    normalisation = support_integral(terminal, radius, kernel, None, support, normal, wall_radius)
    source = current*omega/(normalisation*volumes)
    return CellJunction(index, omega, omega_f, volumes, normalisation, source)


# ---- the support of a terminal ------------------------------------------------------------

def box_support(terminal, lo, hi, radius, tolerance=None):
    """(kind, normal) of the support of a terminal on or inside the box [lo, hi]^3, as
    coupled1D3DMonodomainVerifier::supportGeometry: plane with the inward normal when the
    terminal lies on one face, interior when the ball touches no face; a ball cut by a face
    the terminal is not on, or a terminal on an edge or a corner, is an error."""
    terminal, lo, hi = (numpy.asarray(v, dtype=float) for v in (terminal, lo, hi))
    tol = 1e-9*max(radius, 1.0) if tolerance is None else tolerance
    kind, normal, n_on = "interior", numpy.zeros(3), 0
    for d in range(3):
        to_min, to_max = terminal[d] - lo[d], hi[d] - terminal[d]
        if abs(to_min) <= tol or abs(to_max) <= tol:
            n_on += 1
            kind, normal = "plane", numpy.zeros(3)
            normal[d] = 1.0 if abs(to_min) <= tol else -1.0
            if (to_max if abs(to_min) <= tol else to_min) < radius:
                raise ValueError("support of the terminal at %s reaches the opposite wall"
                                 % (terminal,))
        elif to_min < radius or to_max < radius:
            raise ValueError("support of the terminal at %s is cut by a wall the terminal is "
                             "not on; the continuum reference does not cover that geometry"
                             % (terminal,))
    if n_on > 1:
        raise ValueError("terminal at %s is on an edge or corner of the box" % (terminal,))
    return kind, normal


def shell_support(terminal, r1, r2, radius, tolerance=None):
    """(kind, normal, wall_radius) of the support of a terminal in the sealed shell
    r1 < |x| < r2: cavitySphere (normal the outward radius, wall_radius r1) when the terminal
    is on the inner sphere, interior when the ball touches neither sphere; a terminal on the
    outer sphere, or a ball that meets a sphere the terminal is not on, is an error."""
    terminal = numpy.asarray(terminal, dtype=float)
    tol = 1e-9*max(radius, 1.0) if tolerance is None else tolerance
    rt = numpy.linalg.norm(terminal)
    if abs(rt - r1) <= tol:
        if r1 + radius > r2:
            raise ValueError("support of the terminal at %s reaches the outer sphere: "
                             "r1 + R = %g > r2 = %g" % (terminal, r1 + radius, r2))
        return "cavitySphere", terminal/rt, r1
    if rt - radius < r1 or rt + radius > r2:
        raise ValueError("terminal at %s is neither on the inner sphere nor at least R = %g "
                         "from both spheres" % (terminal, radius))
    return "interior", numpy.zeros(3), 0.0
