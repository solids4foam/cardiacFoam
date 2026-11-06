"""The radial solution satisfies its equation, its wall condition and its scaling."""

import os
import re

import numpy as np
import pytest

import radial_eikonal_mms as rad

HEADER = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "src", "verificationModels",
    "eikonalVerification", "manufacturedRadialEikonalReference.H",
)


@pytest.mark.parametrize("Pe", [5.0, 21.0, 46.0, 94.0])
@pytest.mark.parametrize("a", [0.1, 0.25, 0.5])
def test_equation_residual(Pe, a):
    s = np.linspace(a, a + 1, 20001)[1:-1]
    ds = 1e-5
    w = rad.slope(s, Pe, a)
    dw = (rad.slope(s + ds, Pe, a) - rad.slope(s - ds, Pe, a))/(2*ds)
    assert np.abs((dw + 2*w/s)/Pe - (w - 1)).max() < 1e-6


@pytest.mark.parametrize("Pe", [5.0, 46.0, 94.0])
def test_time_integrates_slope(Pe):
    a = 0.25
    s = np.linspace(a, a + 1, 400001)
    w = rad.slope(s, Pe, a)
    quadrature = np.sum(0.5*(w[1:] + w[:-1])*np.diff(s))
    assert abs(quadrature - rad.scaled_time(np.array(a + 1.0), Pe, a)) < 1e-8


def test_boundary_values():
    Pe, a = 46.0, 0.25
    assert abs(rad.slope(np.array(a + 1.0), Pe, a)) < 1e-14
    assert abs(rad.scaled_time(np.array(a), Pe, a)) < 1e-14


def test_front_limit():
    # Away from the ball and the wall, V dpsi/dr tends to 1 + 2/(Pe s), the curvature slowing
    Pe, a, s = 500.0, 0.25, 0.8
    assert abs(rad.slope(np.array(s), Pe, a) - (1 + 2/(Pe*s))) < 1e-4


def test_header_samples():
    # The self-test samples of the C++ reference come from this solution
    with open(HEADER) as f:
        text = f.read()
    rows = re.findall(r"\{([0-9.]+), ([0-9.]+), ([0-9.]+), ([-0-9.e]+), ([-0-9.e]+)\},", text)
    assert len(rows) == 7
    for row in rows:
        Pe, a, s, w, time = map(float, row)
        assert abs(rad.slope(np.array(s), Pe, a) - w) < 1e-13
        assert abs(rad.scaled_time(np.array(s), Pe, a) - time) < 1e-11


def test_seed_offset_slope_solves_unbounded_radial_equation():
    # pvjMapper::setSeedOffsets integrates s = 1/c0 + 2/(c0^2 rho) + 2/(c0^3 rho^2)
    # of c0 s - s' - 2 s/rho = 1 in the metric distance, from rho = 2/c0
    c0 = 30.0
    rho = np.linspace(2/c0, 5/c0, 7)
    s = 1/c0 + 2/(c0**2*rho) + 2/(c0**3*rho**2)
    ds = -2/(c0**2*rho**2) - 4/(c0**3*rho**3)
    assert np.max(np.abs(c0*s - ds - 2*s/rho - 1)) < 1e-12
    rhoC = 2/c0
    offset = lambda r: (r - rhoC)/c0 + 2*np.log(r/rhoC)/c0**2 - 2*(1/r - 1/rhoC)/c0**3
    assert abs(offset(rhoC)) < 1e-15
    eps = 1e-7
    assert np.max(np.abs((offset(rho + eps) - offset(rho - eps))/(2*eps) - s)) < 1e-7
