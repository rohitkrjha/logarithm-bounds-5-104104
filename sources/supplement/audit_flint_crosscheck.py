#!/usr/bin/env python3
"""Independent FLINT/Arb root cross-check for the 5.104104 theorem.

The theorem verifier uses SymPy's exact VAS isolation for the real critical
polynomial and exact Bernstein signs on the circle.  This audit reconstructs
the same integer polynomials from the input construction, then asks
python-flint/Arb to certify all complex roots.  It is deliberately a different
root-isolation implementation from the theorem verifier.
"""

from __future__ import annotations

from flint import arb, ctx, fmpz_poly

import verify_log3_bound_5_104104 as cert


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def flint_poly(poly) -> fmpz_poly:
    """Convert a SymPy integer polynomial to FLINT's low-to-high format."""
    return fmpz_poly([int(poly.nth(i)) for i in range(poly.degree() + 1)])


def aq(value) -> arb:
    """Embed a SymPy rational exactly as an Arb ball."""
    return arb(int(value.p)) / int(value.q)


def certified_real_roots(poly):
    roots = flint_poly(poly).complex_roots()
    require(sum(int(multiplicity) for _, multiplicity in roots) == poly.degree(),
            "FLINT did not account for the full polynomial degree")
    return [(root.real, int(multiplicity)) for root, multiplicity in roots
            if root.imag.is_zero()]


ctx.prec = 256

print("FLINT: isolating all roots of H", flush=True)
h_real = certified_real_roots(cert.H)
require(all(multiplicity == 1 for _, multiplicity in h_real),
        "H has a multiple real root")

left_count = sum(m for root, m in h_real if root > arb(0) and root < arb(25))
right_count = sum(m for root, m in h_real if root > arb(25) and root < arb(49))
require(left_count == 41, f"FLINT found {left_count}, not 41, H roots in (0,25)")
require(right_count == 22, f"FLINT found {right_count}, not 22, H roots in (25,49)")

# Each VAS interval exported by the theorem verifier must contain exactly one
# of the independently isolated FLINT real-root balls.
for index, (left, right) in enumerate(cert.ROOT_INTERVALS, start=1):
    contained = sum(
        multiplicity
        for root, multiplicity in h_real
        if root > aq(left) and root < aq(right)
    )
    require(contained == 1, f"VAS interval {index} contains {contained} FLINT roots")

print("FLINT: isolating all roots of the degree-169 circle quotient K", flush=True)
k_real = certified_real_roots(cert.K)
k_interval_count = sum(
    multiplicity
    for root, multiplicity in k_real
    if root > arb(-1) and root < arb(1)
)
require(k_interval_count == 0, "FLINT found a K root in (-1,1)")

print("FLINT: checking all 21 circle-factor polynomials", flush=True)
for index, factor in enumerate(cert.CIRCLE_FACTORS, start=1):
    factor_real = certified_real_roots(factor)
    interval_count = sum(
        multiplicity
        for root, multiplicity in factor_real
        if root > arb(-1) and root < arb(1)
    )
    require(interval_count == 0, f"circle factor {index} has a root in (-1,1)")
    require(factor.eval(0) > 0, f"circle factor {index} has the wrong sign")

require(cert.J.eval(cert.C0) == 0, "the exact linear circle root is absent")
require(cert.K.eval(cert.C0) < 0, "K has the wrong sign at the circle root")

print("PASS: FLINT/Arb independently confirms the decisive root claims")
print(f"H roots in (0,25)/(25,49): {left_count}/{right_count}")
print("each of the 63 VAS intervals contains one independently isolated root")
print("K and every circle factor have no root in (-1,1)")
