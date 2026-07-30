#!/usr/bin/env python3
"""Derive and verify the rotated-anisotropy monodomain MMS source.

This script is an independent symbolic audit of the analytic expression in
manufacturedAnisotropicMonodomainReference.H. It is not used at runtime and
does not evaluate the source with a finite-volume operator.
"""

from __future__ import annotations

import sympy as sp


def main() -> None:
    x, y, z, t, beta = sp.symbols("x y z t beta", real=True)
    kxx, kxy, kxz, kyy, kyz, kzz = sp.symbols(
        "kxx kxy kxz kyy kyz kzz", real=True
    )

    conductivity = sp.Matrix(
        [
            [kxx, kxy, kxz],
            [kxy, kyy, kyz],
            [kxz, kyz, kzz],
        ]
    )
    coordinates = (x, y, z)
    spatial_factor = (
        sp.sin(sp.pi * x) ** 2
        * sp.sin(sp.pi * y) ** 2
        * sp.sin(sp.pi * z) ** 2
    )
    voltage = sp.sqrt(1 + t) * spatial_factor
    voltage_gradient = sp.Matrix(
        [sp.diff(voltage, coordinate) for coordinate in coordinates]
    )
    divergence = sp.expand(
        sum(
            sp.diff((conductivity * voltage_gradient)[i], coordinates[i])
            for i in range(3)
        )
    )
    source = sp.simplify(beta * voltage - divergence)

    sx2 = sp.sin(sp.pi * x) ** 2
    sy2 = sp.sin(sp.pi * y) ** 2
    sz2 = sp.sin(sp.pi * z) ** 2
    dx = sp.pi * sp.sin(2 * sp.pi * x)
    dy = sp.pi * sp.sin(2 * sp.pi * y)
    dz = sp.pi * sp.sin(2 * sp.pi * z)
    dxx = 2 * sp.pi**2 * sp.cos(2 * sp.pi * x)
    dyy = 2 * sp.pi**2 * sp.cos(2 * sp.pi * y)
    dzz = 2 * sp.pi**2 * sp.cos(2 * sp.pi * z)

    implemented_divergence = sp.sqrt(1 + t) * (
        kxx * dxx * sy2 * sz2
        + kyy * sx2 * dyy * sz2
        + kzz * sx2 * sy2 * dzz
        + 2 * kxy * dx * dy * sz2
        + 2 * kxz * dx * sy2 * dz
        + 2 * kyz * sx2 * dy * dz
    )
    implemented_source = sp.sqrt(1 + t) * (
        beta * spatial_factor
        - implemented_divergence / sp.sqrt(1 + t)
    )

    divergence_difference = sp.simplify(
        divergence - implemented_divergence
    )
    source_difference = sp.simplify(source - implemented_source)

    if divergence_difference != 0 or source_difference != 0:
        raise RuntimeError(
            "The symbolic MMS source does not match the C++ expression: "
            f"divergence difference={divergence_difference}, "
            f"source difference={source_difference}"
        )

    print("div(K grad(V)) / sqrt(1+t) =")
    print(sp.simplify(divergence / sp.sqrt(1 + t)))
    print("\nS / sqrt(1+t) =")
    print(sp.simplify(source / sp.sqrt(1 + t)))
    print("\nSymbolic comparison with the C++ expression: PASS")


if __name__ == "__main__":
    main()
