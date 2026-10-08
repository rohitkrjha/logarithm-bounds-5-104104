#!/usr/bin/env python3
"""Exact verifier for a computer-assisted irrationality-measure proof.

Claims certified together with the analytic lemmas stated in the accompanying
manuscript:

    |h0 + h1*log(2) + h2*log(3)| > H^(-4.104104)
    for all sufficiently large H, and consequently
    mu(log 3), mu(log(3)/log(2)) <= 5.104104.

The verifier makes no floating-point decisions.  It uses exact integer and
rational arithmetic, exact Vincent--Akritas--Strzebonski real-root isolation,
Bernstein-basis sign certificates, and positive rational series enclosures for
logarithms.  Decimal values are printed only for readability.

Tested with Python 3.13.5 and SymPy 1.14.0.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from decimal import Decimal, getcontext
from fractions import Fraction
from math import comb, gcd
from pathlib import Path
from typing import Iterable, Sequence

import sympy as sp
from sympy.polys import rootisolation as ri
from sympy.polys.domains import ZZ

if not __debug__:
    raise RuntimeError(
        "This verifier relies on assertion checks; rerun Python without -O or PYTHONOPTIMIZE."
    )

START = time.time()

def checkpoint(message: str) -> None:
    print(f"[{time.time() - START:7.2f}s] {message}", flush=True)

# ---------------------------------------------------------------------------
# Exact construction data
# ---------------------------------------------------------------------------

# Each entry is (historical search id, coefficients high-to-low, profile
# (tau_7,tau_5,tau_3,tau_2)).
FACTOR_DATA: list[tuple[int, list[int], tuple[int, int, int, int]]] = [
    (0,   [1, -49], (2, 0, 1, 3)),
    (1,   [1, -25], (0, 2, 1, 4)),
    (2,   [1, 0], (2, 2, 0, 0)),
    (43,  [38553, -19202708, 1088493350, -17941472500, 90075015625], (8, 6, 2, 16)),
    (99,  [463, -2770460, 172374650, -2731137500, 12867859375], (7, 6, 2, 15)),
    (105, [1579, 4041716, -242861150, 3823592500, -18015003125], (8, 5, 1, 16)),
    (124, [11, -38729, 834225, -4501875], (4, 4, 1, 12)),
    (144, [16201, -64316693, 6078999290, -194217490250, 2344523978125, -9457876640625], (9, 7, 2, 20)),
    (148, [3723, -14434700, 867361250, -13655687500, 64339296875], (7, 7, 0, 16)),
    (159, [154073, -706525365, 69673778650, -2292537826250, 27833179828125, -110341894140625], (10, 8, 2, 20)),
    (161, [13627, -789705, 11790625, -52521875], (5, 5, 1, 12)),
    (163, [4259, -19924359, 1880315710, -60526208750, 754056559375, -3152625546875], (9, 7, 1, 20)),
    (178, [8831, -33425731, 2301541270, -49661683750, 434933646875, -1351125234375], (8, 7, 2, 20)),
    (189, [11881, -45475525, 3256944250, -74061846250, 685673078125, -2251875390625], (8, 8, 2, 19)),
    (190, [31809, -138944253, 9336110826, -192597295450, 1574511273125, -4413675765625], (10, 6, 2, 20)),
    (195, [85343, -320906467, 22689947350, -517088363750, 4799711546875, -15763127734375], (9, 8, 3, 19)),
    (198, [21629, -86995629, 6252170386, -142667059850, 1329507230625, -4413675765625], (10, 6, 2, 18)),
    (199, [359787, -1983981727, 196332531150, -6510548598750, 80617138984375, -331025682421875], (10, 8, 3, 20)),
    (201, [14179, -60400263, 4236862910, -96560416750, 918765159375, -3152625546875], (9, 7, 3, 20)),
    (202, [43951, -94114083, 9117300150, -301286483750, 3770282796875, -15763127734375], (9, 8, 2, 19)),
    (208, [162521, -717601717, 70580396250, -2311361666250, 27833179828125, -110341894140625], (10, 8, 2, 19)),
]

# The last fifteen weights are fixed rational numbers FREE_NUM[i]/10^12.
# The first six are the unique solution of the six exact constraints below.
FREE_DEN = 10**12
FREE_NUM = [
    119_873_656,
    878_336_597,
    101_691_070,
    103_130_300,
    219_953_072,
    757_245_952,
    781_996_563,
    232_727_309,
    630_848_726,
    330_420_926,
    232_427_174,
    314_945_936,
    103_903_210,
    14_385_085,
    9_664_407,
]

# Exact rational complex saddle on |z|=26.
Z0_RE = sp.Rational(-23374, 901)
Z0_IM = sp.Rational(1560, 901)
Z0 = Z0_RE + sp.I * Z0_IM
C0 = sp.Rational(-899, 901)
S0 = sp.Rational(60, 901)
assert C0**2 + S0**2 == 1
assert sp.expand(Z0 * sp.conjugate(Z0)) == 26**2

T, C, Y, Z, X = sp.symbols("t c y z X")
POLYS = [sp.Poly.from_list(coeffs, gens=T, domain=ZZ) for _, coeffs, _ in FACTOR_DATA]
PROFILES = [profile for _, _, profile in FACTOR_DATA]
T0 = sp.expand((Z0 - 35) ** 2)
assert all(poly.eval(T0) != 0 for poly in POLYS)


def exact_weights() -> tuple[list[sp.Rational], int, list[int]]:
    """Construct the 21 exact weights and their common denominator."""
    rows: list[list[sp.Rational]] = [
        [sp.Rational(p.degree()) for p in POLYS]
    ]
    for j in range(3):
        rows.append([sp.Rational(profile[j]) for profile in PROFILES])

    log_derivatives = [
        sp.cancel(
            2 * (Z0 - 35) * sp.diff(p.as_expr(), T).subs(T, T0)
            / p.as_expr().subs(T, T0)
        )
        for p in POLYS
    ]
    rows.append([sp.re(value) for value in log_derivatives])
    rows.append([sp.im(value) for value in log_derivatives])

    target = sp.Matrix([
        sp.Rational(3, 2),
        sp.Rational(2),
        sp.Rational(2),
        sp.Rational(1),
        sp.re(1 / Z0 - 1 / (70 - Z0)),
        sp.im(1 / Z0 - 1 / (70 - Z0)),
    ])
    matrix = sp.Matrix(rows)
    assert matrix.rank() == 6
    pivot = matrix[:, :6]
    assert pivot.det() != 0

    free = [sp.Rational(n, FREE_DEN) for n in FREE_NUM]
    first = list(pivot.inv() * (target - matrix[:, 6:] * sp.Matrix(free)))
    alpha = first + free

    assert len(alpha) == len(POLYS) == 21
    assert all(value > 0 for value in alpha)
    assert matrix * sp.Matrix(alpha) == target

    denominator = 1
    for value in alpha:
        denominator = int(sp.ilcm(denominator, value.q))
    numerators = [int(value * denominator) for value in alpha]
    assert all(value > 0 for value in numerators)
    return alpha, denominator, numerators


ALPHA, DEN, WEIGHT_NUM = exact_weights()
checkpoint(f"exact weights reconstructed; common denominator has {len(str(DEN))} digits")

# ---------------------------------------------------------------------------
# Arithmetic profiles and balances
# ---------------------------------------------------------------------------


def valuation(number: int, prime: int) -> int:
    number = abs(int(number))
    if number == 0:
        raise ValueError("valuation of zero is not used")
    result = 0
    while number % prime == 0:
        number //= prime
        result += 1
    return result


def coefficient_profile(poly: sp.Poly) -> tuple[int, int, int, int]:
    scaled = sp.Poly(sp.expand(poly.as_expr().subs(T, (840 * Y - 35) ** 2)), Y, domain=ZZ)
    values = [int(value) for value in scaled.all_coeffs() if value != 0]
    return tuple(min(valuation(value, prime) for value in values) for prime in (7, 5, 3, 2))


COMPUTED_PROFILES = [coefficient_profile(poly) for poly in POLYS]
assert COMPUTED_PROFILES == PROFILES
assert sum(ALPHA[i] * POLYS[i].degree() for i in range(21)) == sp.Rational(3, 2)
assert sum(ALPHA[i] * PROFILES[i][0] for i in range(21)) == 2
assert sum(ALPHA[i] * PROFILES[i][1] for i in range(21)) == 2
assert sum(ALPHA[i] * PROFILES[i][2] for i in range(21)) == 1
TAU2_SUM = sum(ALPHA[i] * PROFILES[i][3] for i in range(21))
CLEARING_ALPHA = 4 - TAU2_SUM
assert CLEARING_ALPHA > 0
CLEARING_NUM = 4 * DEN - sum(WEIGHT_NUM[i] * PROFILES[i][3] for i in range(21))
assert CLEARING_ALPHA == sp.Rational(CLEARING_NUM, DEN)

# Recheck the exact complex stationary equation independently.
LOG_DERIVATIVE_AT_Z0 = sum(
    ALPHA[i]
    * 2
    * (Z0 - 35)
    * sp.diff(POLYS[i].as_expr(), T).subs(T, T0)
    / POLYS[i].as_expr().subs(T, T0)
    for i in range(21)
) - 1 / Z0 + 1 / (70 - Z0)
assert sp.cancel(LOG_DERIVATIVE_AT_Z0) == 0
checkpoint("all coefficient profiles, balances, and the exact saddle equation verified")

# ---------------------------------------------------------------------------
# Rational logarithm enclosures
# ---------------------------------------------------------------------------

SCALE = 10**18
LOG_TERMS = 24


def to_fraction(value) -> Fraction:
    if isinstance(value, Fraction):
        return value
    return Fraction(int(value.numerator), int(value.denominator))


def pow2_fraction(k: int) -> Fraction:
    return Fraction(1 << k, 1) if k >= 0 else Fraction(1, 1 << (-k))


def floor_log2_fraction(value: Fraction) -> int:
    if value <= 0:
        raise ValueError("logarithm input must be positive")
    k = value.numerator.bit_length() - value.denominator.bit_length()
    while value < pow2_fraction(k):
        k -= 1
    while value >= pow2_fraction(k + 1):
        k += 1
    return k


def atanh_log_bounds(z: Fraction, terms: int = LOG_TERMS) -> tuple[Fraction, Fraction]:
    """Bounds 2*atanh(z) by a positive series, for 0 <= z < 1."""
    if not (Fraction(0) <= z < Fraction(1)):
        raise ValueError("atanh argument outside [0,1)")
    if z == 0:
        return Fraction(0), Fraction(0)
    z2 = z * z
    # Horner evaluation of sum_{j=0}^{terms-1} z^(2j)/(2j+1).
    series = Fraction(1, 2 * terms - 1)
    for j in range(terms - 2, -1, -1):
        series = Fraction(1, 2 * j + 1) + z2 * series
    lower = 2 * z * series
    tail = 2 * z * (z2**terms) / ((2 * terms + 1) * (1 - z2))
    return lower, lower + tail


LOG2_LOWER, LOG2_UPPER = atanh_log_bounds(Fraction(1, 3))


def log_bounds(value: Fraction) -> tuple[Fraction, Fraction]:
    """Exact rational lower and upper bounds for log(value), value>0."""
    value = Fraction(value)
    if value <= 0:
        raise ValueError("logarithm input must be positive")
    k = floor_log2_fraction(value)
    scaled = value / pow2_fraction(k)  # in [1,2)
    z = (scaled - 1) / (scaled + 1)
    lower, upper = atanh_log_bounds(z)
    if k >= 0:
        return k * LOG2_LOWER + lower, k * LOG2_UPPER + upper
    return k * LOG2_UPPER + lower, k * LOG2_LOWER + upper


def floor_scaled(value: Fraction) -> int:
    return (value.numerator * SCALE) // value.denominator


def ceil_scaled(value: Fraction) -> int:
    return -((-value.numerator * SCALE) // value.denominator)


def log_fixed(value: Fraction) -> tuple[int, int]:
    """Return outward integer bounds L,U with L/SCALE <= log(value) <= U/SCALE."""
    lower, upper = log_bounds(value)
    return floor_scaled(lower), ceil_scaled(upper)

# ---------------------------------------------------------------------------
# Real critical-point certificate
# ---------------------------------------------------------------------------

PRODUCT_T = sp.Poly(1, T, domain=ZZ)
for poly in POLYS:
    PRODUCT_T *= poly

# Numerator of D*d(log g)/dt.
H = sp.Poly(DEN * PRODUCT_T.as_expr(), T, domain=ZZ)
for i, poly in enumerate(POLYS):
    H += sp.Poly(
        (1225 - T)
        * WEIGHT_NUM[i]
        * sp.diff(poly.as_expr(), T)
        * PRODUCT_T.exquo(poly).as_expr(),
        T,
        domain=ZZ,
    )
H_CONTENT, H = H.primitive()
assert H_CONTENT > 0
assert H.degree() == 85
assert sp.gcd(H, H.diff()).degree() == 0
assert sp.gcd(H, PRODUCT_T).degree() == 0

H_COEFF = [int(value) for value in H.all_coeffs()]
H_HASH = hashlib.sha256(",".join(map(str, H_COEFF)).encode()).hexdigest()
ROOT_INTERVALS = ri.dup_isolate_real_roots_sqf(
    H_COEFF,
    ZZ,
    eps=sp.Rational(1, 10**15),
    inf=sp.Rational(0),
    sup=sp.Rational(49),
    fast=True,
)
assert len(ROOT_INTERVALS) == 63
checkpoint("degree-85 real critical polynomial built and all 63 roots isolated exactly")


def interval_poly_high(coefficients: Sequence[int], a: Fraction, b: Fraction) -> tuple[Fraction, Fraction]:
    """Exact Horner interval evaluation on 0 <= a <= b."""
    if not (0 <= a <= b):
        raise ValueError("this evaluator expects a nonnegative interval")
    lower = upper = Fraction(coefficients[0])
    for coefficient in coefficients[1:]:
        products = (lower * a, lower * b, upper * a, upper * b)
        lower = min(products) + coefficient
        upper = max(products) + coefficient
    return lower, upper


def absolute_interval(lower: Fraction, upper: Fraction) -> tuple[Fraction, Fraction]:
    if lower > 0:
        return lower, upper
    if upper < 0:
        return -upper, -lower
    raise AssertionError("a factor zero lies in a critical-root isolating interval")


def potential_bounds(a: Fraction, b: Fraction) -> tuple[int, int]:
    """Outward fixed-point bounds for log g throughout [a,b]."""
    lower = 0
    upper = 0
    for weight, (_, coefficients, _) in zip(WEIGHT_NUM, FACTOR_DATA):
        p_lower, p_upper = interval_poly_high(coefficients, a, b)
        p_min, p_max = absolute_interval(p_lower, p_upper)
        log_min_lower, _ = log_fixed(p_min)
        _, log_max_upper = log_fixed(p_max)
        lower += weight * log_min_lower
        upper += weight * log_max_upper

    d_min = Fraction(1225) - b
    d_max = Fraction(1225) - a
    _, log_dmax_upper = log_fixed(d_max)
    log_dmin_lower, _ = log_fixed(d_min)
    lower -= DEN * log_dmax_upper
    upper -= DEN * log_dmin_lower
    return lower, upper


REAL_VALUES: list[tuple[Fraction, Fraction, int, int]] = []
for left, right in ROOT_INTERVALS:
    a = to_fraction(left)
    b = to_fraction(right)
    lower, upper = potential_bounds(a, b)
    REAL_VALUES.append((a, b, lower, upper))

REAL_DEN = DEN * SCALE
LEFT_VALUES = [item for item in REAL_VALUES if item[1] < 25]
RIGHT_VALUES = [item for item in REAL_VALUES if item[0] > 25]
assert len(LEFT_VALUES) == 41
assert len(RIGHT_VALUES) == 22
assert len(LEFT_VALUES) + len(RIGHT_VALUES) == 63

M1_LOWER_NUM = max(item[2] for item in LEFT_VALUES)
M1_UPPER_NUM = max(item[3] for item in LEFT_VALUES)
M2_LOWER_NUM = max(item[2] for item in RIGHT_VALUES)
M2_UPPER_NUM = max(item[3] for item in RIGHT_VALUES)
assert M1_LOWER_NUM > M2_UPPER_NUM
checkpoint("all 63 real critical values bounded; the two exponential rates are separated")

# ---------------------------------------------------------------------------
# Complex-circle maximum and exact saddle certificate
# ---------------------------------------------------------------------------


def modulus_squared_on_circle(poly_z: sp.Poly, radius: int) -> sp.Poly:
    """Return |poly_z(r*e^{i theta})|^2 as an integer polynomial in cos(theta)."""
    coefficients = {int(monomial[0]): sp.Rational(value) for monomial, value in poly_z.terms()}
    degrees = sorted(coefficients)
    expression = 0
    for j in degrees:
        expression += coefficients[j] ** 2 * radius ** (2 * j)
    for position, j in enumerate(degrees):
        for k in degrees[position + 1 :]:
            expression += (
                2
                * coefficients[j]
                * coefficients[k]
                * radius ** (j + k)
                * sp.chebyshevt(k - j, C)
            )
    return sp.Poly(sp.expand(expression), C, domain=ZZ)


CIRCLE_FACTORS: list[sp.Poly] = []
for poly in POLYS:
    q_poly = sp.Poly(sp.expand(poly.as_expr().subs(T, (Z - 35) ** 2)), Z, domain=ZZ)
    CIRCLE_FACTORS.append(modulus_squared_on_circle(q_poly, 26))

B_CIRCLE = sp.Poly(676 * (5576 - 3640 * C), C, domain=ZZ)
CIRCLE_PRODUCT = sp.Poly(1, C, domain=ZZ)
for circle_factor in CIRCLE_FACTORS:
    CIRCLE_PRODUCT *= circle_factor

# Numerator of 2*DEN*d/dc log|Phi(26e^{itheta})|.
J = sp.Poly(-DEN * sp.diff(B_CIRCLE.as_expr(), C) * CIRCLE_PRODUCT.as_expr(), C, domain=ZZ)
for i, circle_factor in enumerate(CIRCLE_FACTORS):
    J += sp.Poly(
        B_CIRCLE.as_expr()
        * WEIGHT_NUM[i]
        * sp.diff(circle_factor.as_expr(), C)
        * CIRCLE_PRODUCT.exquo(circle_factor).as_expr(),
        C,
        domain=ZZ,
    )
J_CONTENT, J = J.primitive()
assert J_CONTENT > 0
assert J.degree() == 170
assert J.eval(C0) == 0
assert J.eval(-1) > 0
assert J.eval(1) < 0
assert J.diff().eval(C0) < 0

J_COEFF = [int(value) for value in J.all_coeffs()]
J_HASH = hashlib.sha256(",".join(map(str, J_COEFF)).encode()).hexdigest()
LINEAR_SADDLE = sp.Poly(901 * C + 899, C, domain=ZZ)
K = J.exquo(LINEAR_SADDLE)
assert K.degree() == 169


def bernstein_coefficients_minus_one_one(poly: sp.Poly) -> list[Fraction]:
    """Bernstein coefficients after mapping c=2x-1, x in [0,1]."""
    transformed = sp.Poly(sp.expand(poly.as_expr().subs(C, 2 * X - 1)), X, domain=ZZ)
    degree = transformed.degree()
    monomial = [int(transformed.nth(j)) for j in range(degree + 1)]
    result: list[Fraction] = []
    for k in range(degree + 1):
        value = Fraction(0)
        for j in range(k + 1):
            value += Fraction(monomial[j] * comb(k, j), comb(degree, j))
        result.append(value)
    return result


K_BERNSTEIN = bernstein_coefficients_minus_one_one(K)
assert len(K_BERNSTEIN) == 170
assert all(value < 0 for value in K_BERNSTEIN)

for circle_factor in CIRCLE_FACTORS:
    coefficients = bernstein_coefficients_minus_one_one(circle_factor)
    assert all(value > 0 for value in coefficients)
assert B_CIRCLE.eval(-1) > 0 and B_CIRCLE.eval(1) > 0
checkpoint("circle maximum certified by exact Bernstein-basis signs")

# Exact value of the circle potential at c0.
L_BRACKET_LOWER = 0
L_BRACKET_UPPER = 0
for weight, circle_factor in zip(WEIGHT_NUM, CIRCLE_FACTORS):
    value = to_fraction(circle_factor.eval(C0))
    lower, upper = log_fixed(value)
    L_BRACKET_LOWER += weight * lower
    L_BRACKET_UPPER += weight * upper
b_value = to_fraction(B_CIRCLE.eval(C0))
b_lower, b_upper = log_fixed(b_value)
L_BRACKET_LOWER -= DEN * b_upper
L_BRACKET_UPPER -= DEN * b_lower
# L lies in [L_BRACKET_LOWER,L_BRACKET_UPPER]/(2*DEN*SCALE).

# ---------------------------------------------------------------------------
# Final constants and exact comparison
# ---------------------------------------------------------------------------

LOG2_FIXED_LOWER, LOG2_FIXED_UPPER = log_fixed(Fraction(2))
C_LOWER_NUM = DEN * SCALE + CLEARING_NUM * LOG2_FIXED_LOWER
C_UPPER_NUM = DEN * SCALE + CLEARING_NUM * LOG2_FIXED_UPPER

# lambda < C_upper + L_upper, with denominator 2*DEN*SCALE.
LAMBDA_UPPER_NUM = 2 * C_UPPER_NUM + L_BRACKET_UPPER
# tau_1 > -(C_upper + M1_upper), with denominator DEN*SCALE.
TAU_LOWER_NUM = -(C_UPPER_NUM + M1_UPPER_NUM)
assert TAU_LOWER_NUM > 0

# Target 5.104104 = 638013/125000; homogeneous target is target-1.
TARGET_NUM = 638_013
TARGET_DEN = 125_000
HOMOGENEOUS_TARGET_NUM = TARGET_NUM - TARGET_DEN
# lambda/tau < homogeneous target.  lambda has denominator 2*DEN*SCALE,
# tau has denominator DEN*SCALE.
assert LAMBDA_UPPER_NUM * TARGET_DEN < 2 * TAU_LOWER_NUM * HOMOGENEOUS_TARGET_NUM


def decimal_ratio(numerator: int, denominator: int, digits: int = 30) -> str:
    getcontext().prec = digits
    return str(Decimal(numerator) / Decimal(denominator))


HOMOGENEOUS_UPPER = Fraction(LAMBDA_UPPER_NUM, 2 * TAU_LOWER_NUM)
IRRATIONALITY_UPPER = Fraction(1) + HOMOGENEOUS_UPPER

print("PASS: exact finite certificate verified")
print("active factors =", len(POLYS))
print("weight common denominator digits =", len(str(DEN)))
print("clearing exponent alpha =", decimal_ratio(CLEARING_NUM, DEN, 32))
print("real critical polynomial: degree", H.degree(), "with", len(ROOT_INTERVALS), "roots in [0,49]")
print("critical roots in [0,25] / [25,49] =", len(LEFT_VALUES), "/", len(RIGHT_VALUES))
print("M1 in [",
      decimal_ratio(M1_LOWER_NUM, REAL_DEN, 28), ",",
      decimal_ratio(M1_UPPER_NUM, REAL_DEN, 28), "]")
print("M2 in [",
      decimal_ratio(M2_LOWER_NUM, REAL_DEN, 28), ",",
      decimal_ratio(M2_UPPER_NUM, REAL_DEN, 28), "]")
print("M1-M2 >", decimal_ratio(M1_LOWER_NUM - M2_UPPER_NUM, REAL_DEN, 24))
print("circle derivative numerator: degree", J.degree())
print("J=(901*c+899)*K and all", len(K_BERNSTEIN), "Bernstein coefficients of K are negative")
print("circle maximum c0 = -899/901, z0 = 26*(-899+60*i)/901")
print("L in [",
      decimal_ratio(L_BRACKET_LOWER, 2 * REAL_DEN, 28), ",",
      decimal_ratio(L_BRACKET_UPPER, 2 * REAL_DEN, 28), "]")
print("C in [",
      decimal_ratio(C_LOWER_NUM, REAL_DEN, 28), ",",
      decimal_ratio(C_UPPER_NUM, REAL_DEN, 28), "]")
print("lambda <", decimal_ratio(LAMBDA_UPPER_NUM, 2 * REAL_DEN, 28))
print("tau_1 >", decimal_ratio(TAU_LOWER_NUM, REAL_DEN, 28))
print("lambda/tau_1 <",
      decimal_ratio(HOMOGENEOUS_UPPER.numerator, HOMOGENEOUS_UPPER.denominator, 30))
print("1 + lambda/tau_1 <",
      decimal_ratio(IRRATIONALITY_UPPER.numerator, IRRATIONALITY_UPPER.denominator, 30))
print("certified homogeneous target = 4.104104")
print("certified irrationality-exponent target = 5.104104")
print("Therefore, with the manuscript's analytic lemmas, both stated claims follow")


def rational_string(value: sp.Rational | Fraction) -> str:
    if isinstance(value, Fraction):
        return f"{value.numerator}/{value.denominator}"
    return f"{int(value.p)}/{int(value.q)}"


def interval_json(item) -> dict[str, str]:
    a, b, lower, upper = item
    return {
        "left": rational_string(a),
        "right": rational_string(b),
        "potential_lower_numerator": str(lower),
        "potential_upper_numerator": str(upper),
        "potential_denominator": str(REAL_DEN),
    }


def write_json(path: Path) -> None:
    payload = {
        "claims": {
            "homogeneous_exponent": "4.104104",
            "irrationality_exponents": {
                "log(3)": "5.104104",
                "log(3)/log(2)": "5.104104",
            },
        },
        "status": "exact finite certificate; analytic lemmas are proved in the companion manuscript",
        "factor_data": [
            {"search_id": idx, "coefficients_high_to_low": coeffs, "profile_7_5_3_2": list(profile)}
            for idx, coeffs, profile in FACTOR_DATA
        ],
        "free_weight_denominator": FREE_DEN,
        "free_weight_numerators": FREE_NUM,
        "exact_weights": [rational_string(value) for value in ALPHA],
        "common_denominator": str(DEN),
        "weight_numerators": [str(value) for value in WEIGHT_NUM],
        "clearing_alpha": rational_string(CLEARING_ALPHA),
        "saddle": {
            "z0_real": rational_string(Z0_RE),
            "z0_imag": rational_string(Z0_IM),
            "c0": rational_string(C0),
            "s0": rational_string(S0),
        },
        "real_polynomial": {
            "degree": H.degree(),
            "coefficient_sha256": H_HASH,
            "root_count_0_49": len(ROOT_INTERVALS),
            "root_intervals_and_bounds": [interval_json(item) for item in REAL_VALUES],
        },
        "circle_polynomial": {
            "degree": J.degree(),
            "coefficient_sha256": J_HASH,
            "quotient_degree": K.degree(),
            "negative_bernstein_coefficients": len(K_BERNSTEIN),
        },
        "bounds": {
            "scale": SCALE,
            "M1_lower_numerator": str(M1_LOWER_NUM),
            "M1_upper_numerator": str(M1_UPPER_NUM),
            "M2_lower_numerator": str(M2_LOWER_NUM),
            "M2_upper_numerator": str(M2_UPPER_NUM),
            "real_denominator": str(REAL_DEN),
            "L_lower_numerator": str(L_BRACKET_LOWER),
            "L_upper_numerator": str(L_BRACKET_UPPER),
            "L_denominator": str(2 * REAL_DEN),
            "C_lower_numerator": str(C_LOWER_NUM),
            "C_upper_numerator": str(C_UPPER_NUM),
            "C_denominator": str(REAL_DEN),
            "lambda_upper_numerator": str(LAMBDA_UPPER_NUM),
            "lambda_upper_denominator": str(2 * REAL_DEN),
            "tau_lower_numerator": str(TAU_LOWER_NUM),
            "tau_lower_denominator": str(REAL_DEN),
            "homogeneous_upper": rational_string(HOMOGENEOUS_UPPER),
            "irrationality_upper": rational_string(IRRATIONALITY_UPPER),
        },
    }
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print("wrote certificate JSON:", path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", type=Path, help="write the exact certificate data to this JSON file")
    args = parser.parse_args()
    if args.json:
        write_json(args.json)
