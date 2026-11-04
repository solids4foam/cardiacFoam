#!/usr/bin/env python3
"""Exact source-free eikonal solution in a sealed sphere around a held ball.

With an isotropic M = m I and a field of the radius alone, the eikonal operator
of the solver

    -div(M grad psi) + c0 sqrt(grad psi . M grad psi) = 1

reduces to -m (psi'' + 2 psi'/r) + V psi' = 1, V = c0 sqrt(m), with psi(R) = Tk on
the ball and psi'(Ro) = 0 on the sealed wall. With L = Ro - R, s = r/L, a = R/L,
b = Ro/L, Pe = c0 L/sqrt(m) and w = V psi' this is the linear equation

    (w' + 2 w/s)/Pe = w - 1,   w(b) = 0,

so that, with q(s) = s^2 + 2 s/Pe + 2/Pe^2,

    w(s) = (q(s) - q(b) exp(-Pe (b - s)))/s^2,
    psi(s) = Tk + (L/V) (W(s) - W(a)),
    W(s) = s + 2 log(s)/Pe - 2/(Pe^2 s)
           - q(b) exp(-Pe b) (Pe Ei(Pe s) - exp(Pe s)/s).

The solution extends below R and is singular at r = 0.

A constant anisotropic M = m Q diag(a^2) Q^T keeps the solution exact in the
coordinate xi = diag(1/a) Q^T x, with r = |xi|: the domain and the ball are
ellipsoids, and the sealed condition n.M grad(psi) = 0 maps to psi'(Ro) = 0.

The closed form needs Pe (a + 1) below about 700, where exp(Pe s) still fits a double.
The C++ reference manufacturedRadialEikonalReference.H evaluates the same solution, and
test_radial_eikonal_mms.py checks its samples.

    python radial_eikonal_mms.py --PeL 46 --RoverL 0.25
"""

import argparse

import numpy as np
from scipy.special import expi


def q(s, Pe):
    return s*s + 2*s/Pe + 2/Pe**2


def slope(s, Pe, a):
    """w = V dpsi/dr, 1 for a plane front."""
    b = a + 1
    return (q(s, Pe) - q(b, Pe)*np.exp(-Pe*(b - s)))/s**2


def scaled_time(s, Pe, a):
    """V (psi - Tk)/L, zero at s = a."""
    b = a + 1

    def W(t):
        return (
            t + 2*np.log(t)/Pe - 2/(Pe*Pe*t)
          - q(b, Pe)*np.exp(-Pe*b)*(Pe*expi(Pe*t) - np.exp(Pe*t)/t)
        )

    return W(s) - W(a)


def rad_arrival(Pe, a):
    return float(scaled_time(np.array(a + 1.0), Pe, a))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--PeL", type=float, default=46.0)
    ap.add_argument("--RoverL", type=float, default=0.25)
    args = ap.parse_args()
    a = args.RoverL
    if args.PeL*(a + 1) > 700:
        ap.error("--PeL times (R/L + 1) above 700 overflows the closed form")
    print("Pe_L %g, R/L %g: scaled arrival at the wall %.6f, start slope %.4f" % (
        args.PeL, a, rad_arrival(args.PeL, a), slope(np.array(a), args.PeL, a)))


if __name__ == "__main__":
    main()
