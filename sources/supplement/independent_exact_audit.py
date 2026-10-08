#!/usr/bin/env python3
"""
Independent exact audit for the proposed bound mu(log 3) <= 5.104104.

This checker intentionally does NOT import SymPy or any other CAS.  It uses
only Python's standard library (Fraction and arbitrary-precision integers).

It independently:
  * reconstructs all 21 rational weights by exact 6x6 Gaussian elimination;
  * recomputes all 84 coefficient-content valuation profiles;
  * reconstructs the degree-85 real critical polynomial and its hash;
  * verifies every supplied root interval and every complementary gap using
    exact Descartes transformations (variation 1 / variation 0);
  * recomputes every rational logarithm enclosure and both real maxima;
  * reconstructs all 21 circle-modulus polynomials and the degree-170
    derivative numerator;
  * verifies the exact saddle factor and all Bernstein coefficient signs;
  * recomputes the final exact homogeneous and irrationality-exponent
    comparisons below 4.104104 and 5.104104.

The finite certificate remains logically conditional on the analytic lemmas
proved in the manuscript (the positive Laplace principle, the lcm asymptotic,
and the limit-superior two-form criterion).  No saddle asymptotic or phase
noncancellation argument is used in the publication proof.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from math import comb, gcd
from pathlib import Path

if not __debug__:
    raise RuntimeError(
        "This audit relies on assertion checks; rerun Python without -O or PYTHONOPTIMIZE."
    )


def F(value) -> Fraction:
    if isinstance(value, Fraction):
        return value
    if isinstance(value, int):
        return Fraction(value)
    if isinstance(value, str):
        if "/" in value:
            a, b = value.split("/")
            return Fraction(int(a), int(b))
        return Fraction(int(value))
    return Fraction(value)


def lcm(a: int, b: int) -> int:
    return abs(a // gcd(a, b) * b)


# ---------------------------------------------------------------------------
# Integer/rational polynomial arithmetic, coefficients low-to-high
# ---------------------------------------------------------------------------

def trim(poly):
    out = list(poly)
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def padd(a, b):
    n = max(len(a), len(b))
    return trim([
        (a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
        for i in range(n)
    ])


def psub(a, b):
    n = max(len(a), len(b))
    return trim([
        (a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0)
        for i in range(n)
    ])


def pscale(poly, scalar):
    return trim([scalar * x for x in poly])


def pmul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    out[i + j] += x * y
    return trim(out)


def pderiv(poly):
    return trim([i * poly[i] for i in range(1, len(poly))] or [0])


def peval(poly, x):
    result = 0
    for coefficient in reversed(poly):
        result = result * x + coefficient
    return result


def pcompose(poly, inner):
    result = [0]
    for coefficient in reversed(poly):
        result = padd(pmul(result, inner), [coefficient])
    return result


def pdiv_linear_exact_int(poly, constant: int, linear: int):
    """Exact division by constant + linear*x."""
    remainder = list(poly)
    quotient = [0] * max(1, len(poly) - 1)
    divisor = [constant, linear]
    while len(remainder) >= 2 and any(remainder):
        degree = len(remainder) - 2
        lead = remainder[-1]
        assert lead % linear == 0
        coefficient = lead // linear
        quotient[degree] = coefficient
        remainder[degree] -= coefficient * constant
        remainder[degree + 1] -= coefficient * linear
        remainder = trim(remainder)
    assert all(value == 0 for value in remainder)
    return trim(quotient)


def primitive(poly):
    content = 0
    for value in poly:
        content = gcd(content, abs(value))
    assert content > 0
    return [value // content for value in poly], content


# ---------------------------------------------------------------------------
# Exact rational Gaussian arithmetic
# ---------------------------------------------------------------------------

def cadd(a, b):
    return (a[0] + b[0], a[1] + b[1])


def csub(a, b):
    return (a[0] - b[0], a[1] - b[1])


def cmul(a, b):
    return (a[0] * b[0] - a[1] * b[1],
            a[0] * b[1] + a[1] * b[0])


def cinv(a):
    denominator = a[0] * a[0] + a[1] * a[1]
    assert denominator != 0
    return (a[0] / denominator, -a[1] / denominator)


def cdiv(a, b):
    return cmul(a, cinv(b))


def cscale(a, scalar):
    return (a[0] * scalar, a[1] * scalar)


def cpeval(poly, z):
    result = (Fraction(0), Fraction(0))
    for coefficient in reversed(poly):
        result = cadd(cmul(result, z), (Fraction(coefficient), Fraction(0)))
    return result


# ---------------------------------------------------------------------------
# Exact 6x6 solve
# ---------------------------------------------------------------------------

def solve_linear(matrix, vector):
    n = len(matrix)
    augmented = [
        [F(x) for x in row] + [F(vector[i])]
        for i, row in enumerate(matrix)
    ]
    assert all(len(row) == n + 1 for row in augmented)

    for column in range(n):
        pivot = next(row for row in range(column, n)
                     if augmented[row][column] != 0)
        if pivot != column:
            augmented[column], augmented[pivot] = \
                augmented[pivot], augmented[column]

        pivot_value = augmented[column][column]
        augmented[column] = [x / pivot_value for x in augmented[column]]

        for row in range(n):
            if row == column:
                continue
            factor = augmented[row][column]
            if factor:
                augmented[row] = [
                    augmented[row][j] - factor * augmented[column][j]
                    for j in range(n + 1)
                ]
    return [augmented[i][-1] for i in range(n)]


# ---------------------------------------------------------------------------
# Valuations and rational logarithm enclosures
# ---------------------------------------------------------------------------

def valuation(number: int, prime: int) -> int:
    number = abs(number)
    assert number != 0
    result = 0
    while number % prime == 0:
        number //= prime
        result += 1
    return result


SCALE = 10**18
LOG_TERMS = 24


def pow2_fraction(k: int) -> Fraction:
    return Fraction(1 << k, 1) if k >= 0 else Fraction(1, 1 << (-k))


def floor_log2_fraction(value: Fraction) -> int:
    assert value > 0
    k = value.numerator.bit_length() - value.denominator.bit_length()
    while value < pow2_fraction(k):
        k -= 1
    while value >= pow2_fraction(k + 1):
        k += 1
    return k


def atanh_log_bounds(z: Fraction, terms: int = LOG_TERMS):
    assert 0 <= z < 1
    if z == 0:
        return Fraction(0), Fraction(0)
    z2 = z * z
    series = Fraction(1, 2 * terms - 1)
    for j in range(terms - 2, -1, -1):
        series = Fraction(1, 2 * j + 1) + z2 * series
    lower = 2 * z * series
    tail = 2 * z * (z2 ** terms) / ((2 * terms + 1) * (1 - z2))
    return lower, lower + tail


LOG2_LOWER, LOG2_UPPER = atanh_log_bounds(Fraction(1, 3))


def log_bounds(value: Fraction):
    assert value > 0
    k = floor_log2_fraction(value)
    scaled = value / pow2_fraction(k)
    z = (scaled - 1) / (scaled + 1)
    lower, upper = atanh_log_bounds(z)
    if k >= 0:
        return k * LOG2_LOWER + lower, k * LOG2_UPPER + upper
    return k * LOG2_UPPER + lower, k * LOG2_LOWER + upper


def floor_scaled(value: Fraction) -> int:
    return value.numerator * SCALE // value.denominator


def ceil_scaled(value: Fraction) -> int:
    return -((-value.numerator * SCALE) // value.denominator)


def log_fixed(value: Fraction):
    lower, upper = log_bounds(value)
    return floor_scaled(lower), ceil_scaled(upper)


def interval_poly_high(coefficients, a: Fraction, b: Fraction):
    assert 0 <= a <= b
    lower = upper = Fraction(coefficients[0])
    for coefficient in coefficients[1:]:
        products = (
            lower * a, lower * b, upper * a, upper * b
        )
        lower = min(products) + coefficient
        upper = max(products) + coefficient
    return lower, upper


def absolute_interval(lower: Fraction, upper: Fraction):
    if lower > 0:
        return lower, upper
    if upper < 0:
        return -upper, -lower
    raise AssertionError("factor zero inside a critical-root interval")


# ---------------------------------------------------------------------------
# Descartes interval transform
# ---------------------------------------------------------------------------

def mul_linear(poly, constant: int, linear: int):
    out = [0] * (len(poly) + 1)
    for i, value in enumerate(poly):
        out[i] += value * constant
        out[i + 1] += value * linear
    return out


def descartes_transform(poly, a: Fraction, b: Fraction, binomials):
    """
    Return integer coefficients of a positive scalar multiple of
      (1+x)^n p((a+b*x)/(1+x)).
    Positive roots correspond exactly to roots of p in (a,b).
    """
    A, C = a.numerator, a.denominator
    B, D = b.numerator, b.denominator
    L0 = A * D
    L1 = B * C
    CD = C * D

    n = len(poly) - 1
    result = [poly[n]]
    cd_power = 1
    for k in range(n - 1, -1, -1):
        j = n - k
        result = mul_linear(result, L0, L1)
        cd_power *= CD
        scale = poly[k] * cd_power
        for index, coefficient in enumerate(binomials[j]):
            result[index] += scale * coefficient
    return trim(result)


def sign_variations(poly) -> int:
    signs = [1 if value > 0 else -1 for value in reversed(poly) if value]
    return sum(a != b for a, b in zip(signs, signs[1:]))


# ---------------------------------------------------------------------------
# Circle polynomials and Bernstein coefficients
# ---------------------------------------------------------------------------

def chebyshev_polynomials(nmax: int):
    result = [[1]]
    if nmax >= 1:
        result.append([0, 1])
    for n in range(1, nmax):
        result.append(psub(pscale(pmul([0, 1], result[n]), 2),
                           result[n - 1]))
    return result


def modulus_squared_on_circle(poly_z, radius: int, chebyshev):
    degree = len(poly_z) - 1
    result = [0]
    for j in range(degree + 1):
        result[0] += poly_z[j] * poly_z[j] * radius ** (2 * j)
        for k in range(j + 1, degree + 1):
            result = padd(
                result,
                pscale(
                    chebyshev[k - j],
                    2 * poly_z[j] * poly_z[k] * radius ** (j + k),
                ),
            )
    return trim(result)


def substitute_c_eq_2x_minus_1(poly):
    out = [0] * len(poly)
    for k, value in enumerate(poly):
        if value == 0:
            continue
        for j in range(k + 1):
            out[j] += value * comb(k, j) * (2 ** j) * ((-1) ** (k - j))
    return trim(out)


def bernstein_coefficients(poly):
    transformed = substitute_c_eq_2x_minus_1(poly)
    degree = len(transformed) - 1
    result = []
    for k in range(degree + 1):
        value = Fraction(0)
        for j in range(k + 1):
            value += Fraction(
                transformed[j] * comb(k, j),
                comb(degree, j),
            )
        result.append(value)
    return result


def main(certificate_path: Path) -> None:
    cert = json.loads(certificate_path.read_text(encoding="utf-8"))
    factor_data = cert["factor_data"]
    polynomials = [
        list(reversed(row["coefficients_high_to_low"]))
        for row in factor_data
    ]
    profiles = [
        tuple(row["profile_7_5_3_2"])
        for row in factor_data
    ]
    degrees = [len(poly) - 1 for poly in polynomials]
    assert len(polynomials) == 21

    z0 = (F(cert["saddle"]["z0_real"]), F(cert["saddle"]["z0_imag"]))
    z0_minus_35 = csub(z0, (Fraction(35), Fraction(0)))
    t0 = cmul(z0_minus_35, z0_minus_35)

    ell = []
    for poly in polynomials:
        derivative = pderiv(poly)
        ell.append(cmul(
            z0_minus_35,
            cscale(cdiv(cpeval(derivative, t0), cpeval(poly, t0)), 2),
        ))

    rhs = csub(
        cinv(z0),
        cinv(csub((Fraction(70), Fraction(0)), z0)),
    )
    rows = [
        [F(x) for x in degrees],
        [F(profile[0]) for profile in profiles],
        [F(profile[1]) for profile in profiles],
        [F(profile[2]) for profile in profiles],
        [value[0] for value in ell],
        [value[1] for value in ell],
    ]
    target = [F(3) / 2, F(2), F(2), F(1), rhs[0], rhs[1]]

    free = [
        F(value) / F(cert["free_weight_denominator"])
        for value in cert["free_weight_numerators"]
    ]
    pivot = [row[:6] for row in rows]
    reduced_target = [
        target[i] - sum(rows[i][6 + j] * free[j] for j in range(15))
        for i in range(6)
    ]
    weights = solve_linear(pivot, reduced_target) + free
    assert all(weight > 0 for weight in weights)
    assert [f"{x.numerator}/{x.denominator}" for x in weights] == \
        cert["exact_weights"]

    common_denominator = 1
    for weight in weights:
        common_denominator = lcm(common_denominator, weight.denominator)
    weight_numerators = [
        weight.numerator * (common_denominator // weight.denominator)
        for weight in weights
    ]
    assert common_denominator == int(cert["common_denominator"])
    assert weight_numerators == [int(x) for x in cert["weight_numerators"]]

    assert sum(weights[i] * degrees[i] for i in range(21)) == F(3) / 2
    assert sum(weights[i] * profiles[i][0] for i in range(21)) == 2
    assert sum(weights[i] * profiles[i][1] for i in range(21)) == 2
    assert sum(weights[i] * profiles[i][2] for i in range(21)) == 1

    # Coefficient-content profiles.
    t_scaled = [1225, -58800, 705600]
    computed_profiles = []
    for poly in polynomials:
        scaled = pcompose(poly, t_scaled)
        nonzero = [x for x in scaled if x]
        computed_profiles.append(tuple(
            min(valuation(x, prime) for x in nonzero)
            for prime in (7, 5, 3, 2)
        ))
    assert computed_profiles == profiles

    # Prefix/suffix products.
    count = len(polynomials)
    prefix = [[1]]
    for poly in polynomials:
        prefix.append(pmul(prefix[-1], poly))
    suffix = [None] * (count + 1)
    suffix[count] = [1]
    for i in range(count - 1, -1, -1):
        suffix[i] = pmul(polynomials[i], suffix[i + 1])
    product = prefix[count]

    # Real derivative numerator H.
    sum_term = [0]
    for i, poly in enumerate(polynomials):
        product_except = pmul(prefix[i], suffix[i + 1])
        term = pmul(pderiv(poly), product_except)
        sum_term = padd(sum_term, pscale(term, weight_numerators[i]))
    H_raw = padd(
        pscale(product, common_denominator),
        pmul([1225, -1], sum_term),
    )
    H, _ = primitive(H_raw)
    assert len(H) - 1 == 85
    H_hash = hashlib.sha256(
        ",".join(map(str, reversed(H))).encode()
    ).hexdigest()
    assert H_hash == cert["real_polynomial"]["coefficient_sha256"]

    # Exact root coverage by Descartes transformations.
    root_data = cert["real_polynomial"]["root_intervals_and_bounds"]
    intervals = [(F(item["left"]), F(item["right"])) for item in root_data]
    assert len(intervals) == 63
    assert all(F(0) < a < b < F(49) for a, b in intervals)
    assert all(intervals[i][1] < intervals[i + 1][0]
               for i in range(len(intervals) - 1))
    endpoints = {F(0), F(25), F(49)}
    endpoints.update(value for interval in intervals for value in interval)
    assert all(peval(H, endpoint) != 0 for endpoint in endpoints)
    binomials = [
        [comb(j, k) for k in range(j + 1)]
        for j in range(len(H))
    ]
    for a, b in intervals:
        assert sign_variations(descartes_transform(H, a, b, binomials)) == 1
    previous = Fraction(0)
    for a, b in intervals:
        assert sign_variations(
            descartes_transform(H, previous, a, binomials)
        ) == 0
        previous = b
    assert sign_variations(
        descartes_transform(H, previous, Fraction(49), binomials)
    ) == 0

    # Real potential intervals.
    def potential_bounds(a: Fraction, b: Fraction):
        lower = upper = 0
        for weight, row in zip(weight_numerators, factor_data):
            p_lower, p_upper = interval_poly_high(
                row["coefficients_high_to_low"], a, b
            )
            p_min, p_max = absolute_interval(p_lower, p_upper)
            log_min_lower, _ = log_fixed(p_min)
            _, log_max_upper = log_fixed(p_max)
            lower += weight * log_min_lower
            upper += weight * log_max_upper
        d_min = Fraction(1225) - b
        d_max = Fraction(1225) - a
        _, log_dmax_upper = log_fixed(d_max)
        log_dmin_lower, _ = log_fixed(d_min)
        lower -= common_denominator * log_dmax_upper
        upper -= common_denominator * log_dmin_lower
        return lower, upper

    real_values = []
    for interval, item in zip(intervals, root_data):
        lower, upper = potential_bounds(*interval)
        assert lower == int(item["potential_lower_numerator"])
        assert upper == int(item["potential_upper_numerator"])
        real_values.append((interval[0], interval[1], lower, upper))

    left_values = [item for item in real_values if item[1] < 25]
    right_values = [item for item in real_values if item[0] > 25]
    assert len(left_values) == 41 and len(right_values) == 22

    M1_lower = max(item[2] for item in left_values)
    M1_upper = max(item[3] for item in left_values)
    M2_lower = max(item[2] for item in right_values)
    M2_upper = max(item[3] for item in right_values)
    assert M1_lower == int(cert["bounds"]["M1_lower_numerator"])
    assert M1_upper == int(cert["bounds"]["M1_upper_numerator"])
    assert M2_lower == int(cert["bounds"]["M2_lower_numerator"])
    assert M2_upper == int(cert["bounds"]["M2_upper_numerator"])
    assert M1_lower > M2_upper

    # Circle factors.
    q_polynomials = [
        pcompose(poly, [1225, -70, 1])
        for poly in polynomials
    ]
    max_q_degree = max(len(poly) - 1 for poly in q_polynomials)
    chebyshev = chebyshev_polynomials(max_q_degree)
    circle_factors = [
        modulus_squared_on_circle(poly, 26, chebyshev)
        for poly in q_polynomials
    ]
    B_circle = [676 * 5576, -676 * 3640]

    prefix_circle = [[1]]
    for poly in circle_factors:
        prefix_circle.append(pmul(prefix_circle[-1], poly))
    suffix_circle = [None] * (count + 1)
    suffix_circle[count] = [1]
    for i in range(count - 1, -1, -1):
        suffix_circle[i] = pmul(circle_factors[i], suffix_circle[i + 1])
    circle_product = prefix_circle[count]

    circle_sum = [0]
    for i, poly in enumerate(circle_factors):
        product_except = pmul(prefix_circle[i], suffix_circle[i + 1])
        circle_sum = padd(
            circle_sum,
            pscale(
                pmul(pderiv(poly), product_except),
                weight_numerators[i],
            ),
        )
    J_raw = padd(
        pscale(pmul(pderiv(B_circle), circle_product),
               -common_denominator),
        pmul(B_circle, circle_sum),
    )
    J, _ = primitive(J_raw)
    assert len(J) - 1 == 170
    J_hash = hashlib.sha256(
        ",".join(map(str, reversed(J))).encode()
    ).hexdigest()
    assert J_hash == cert["circle_polynomial"]["coefficient_sha256"]

    c0 = F(cert["saddle"]["c0"])
    assert peval(J, c0) == 0
    assert peval(pderiv(J), c0) < 0
    K = pdiv_linear_exact_int(J, 899, 901)
    assert len(K) - 1 == 169
    K_bernstein = bernstein_coefficients(K)
    assert len(K_bernstein) == 170
    assert all(value < 0 for value in K_bernstein)
    assert all(
        all(value > 0 for value in bernstein_coefficients(poly))
        for poly in circle_factors
    )

    # Circle potential L.
    L_lower = L_upper = 0
    for weight, circle_factor in zip(weight_numerators, circle_factors):
        value = peval(circle_factor, c0)
        assert value > 0
        lower, upper = log_fixed(value)
        L_lower += weight * lower
        L_upper += weight * upper
    b_value = peval(B_circle, c0)
    b_lower, b_upper = log_fixed(b_value)
    L_lower -= common_denominator * b_upper
    L_upper -= common_denominator * b_lower
    assert L_lower == int(cert["bounds"]["L_lower_numerator"])
    assert L_upper == int(cert["bounds"]["L_upper_numerator"])

    # Final exact bound.
    clearing_numerator = (
        4 * common_denominator
        - sum(weight_numerators[i] * profiles[i][3] for i in range(21))
    )
    assert Fraction(clearing_numerator, common_denominator) == \
        F(cert["clearing_alpha"])

    log2_lower, log2_upper = log_fixed(Fraction(2))
    C_lower = common_denominator * SCALE + clearing_numerator * log2_lower
    C_upper = common_denominator * SCALE + clearing_numerator * log2_upper
    assert C_lower == int(cert["bounds"]["C_lower_numerator"])
    assert C_upper == int(cert["bounds"]["C_upper_numerator"])

    lambda_upper = 2 * C_upper + L_upper
    tau_lower = -(C_upper + M1_upper)
    assert lambda_upper == int(cert["bounds"]["lambda_upper_numerator"])
    assert tau_lower == int(cert["bounds"]["tau_lower_numerator"])
    assert tau_lower > 0

    target_num = 638_013
    target_den = 125_000
    homogeneous_num = target_num - target_den
    assert lambda_upper * target_den < \
        2 * tau_lower * homogeneous_num

    homogeneous_upper = Fraction(lambda_upper, 2 * tau_lower)
    irrationality_upper = Fraction(1) + homogeneous_upper
    assert f"{homogeneous_upper.numerator}/{homogeneous_upper.denominator}" == \
        cert["bounds"]["homogeneous_upper"]
    assert f"{irrationality_upper.numerator}/{irrationality_upper.denominator}" == \
        cert["bounds"]["irrationality_upper"]

    print("PASS: independent standard-library exact audit")
    print("CAS dependencies: none")
    print("weights reconstructed:", len(weights))
    print("common denominator digits:", len(str(common_denominator)))
    print("valuation profiles recomputed:", len(computed_profiles) * 4)
    print("H degree/hash:", len(H) - 1, H_hash)
    print("Descartes root intervals/gaps:", len(intervals), len(intervals) + 1)
    print("real roots in (0,25)/(25,49):", len(left_values), len(right_values))
    print("strict real-rate separation:", M1_lower > M2_upper)
    print("J degree/hash:", len(J) - 1, J_hash)
    print("negative Bernstein coefficients of K:", len(K_bernstein))
    print("exact final cross-product: PASS")
    print("homogeneous target = 4.104104")
    print("irrationality-exponent target = 5.104104")
    print("STATUS: finite certificate independently reproduced;")
    print("        analytic lemmas remain outside this executable audit.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--certificate",
        type=Path,
        default=Path(__file__).with_name("certificate_log3_5_104104.json"),
    )
    args = parser.parse_args()
    main(args.certificate)
