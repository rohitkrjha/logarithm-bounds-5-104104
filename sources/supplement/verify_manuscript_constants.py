#!/usr/bin/env python3
"""Dependency-free consistency checks between the paper and exact JSON data."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path


HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument(
    "--manuscript",
    type=Path,
    help="path to the submitted manuscript.tex (default: one directory above this script)",
)
args = parser.parse_args()
manuscript_path = args.manuscript or (HERE.parent / "manuscript.tex")
if not manuscript_path.is_file():
    raise RuntimeError(
        f"manuscript source not found at {manuscript_path}; "
        "supply it with --manuscript /path/to/manuscript.tex"
    )
MANUSCRIPT = manuscript_path.read_text(encoding="utf-8")
COMPACT_MANUSCRIPT = "".join(MANUSCRIPT.split())
DATA = json.loads(
    (HERE / "certificate_log3_5_104104.json").read_text(encoding="utf-8")
)
BOUNDS = DATA["bounds"]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def decimal_fraction(text: str) -> Fraction:
    sign = -1 if text.startswith("-") else 1
    unsigned = text.lstrip("+-")
    whole, fractional = unsigned.split(".")
    return sign * Fraction(
        int(whole + fractional),
        10 ** len(fractional),
    )


def certified(numerator_key: str, denominator_key: str) -> Fraction:
    return Fraction(
        int(BOUNDS[numerator_key]),
        int(BOUNDS[denominator_key]),
    )


display = {
    "M1_lower": "-2.722957495668415052959634",
    "M1_upper": "-2.722957495668409711051756",
    "M2_lower": "-2.722977562389619550059089",
    "M2_upper": "-2.722977562389619502813003",
    "L_lower": "4.537901075598406403746477",
    "L_upper": "4.537901075598406404986681",
    "C_lower": "1.300404222911437365598758",
    "C_upper": "1.300404222911437366032151",
    "lambda_upper": "5.838305298509843771018831",
    "tau_lower": "1.422553272756972345019606",
    "homogeneous_upper": "4.104103101316511697013643",
}

for value in display.values():
    require(value in MANUSCRIPT, f"manuscript is missing displayed value {value}")

for index, factor in enumerate(DATA["factor_data"], start=1):
    coefficient_text = "[" + ", ".join(
        str(value) for value in factor["coefficients_high_to_low"]
    ) + "]"
    profile_text = "&".join(str(value) for value in factor["profile_7_5_3_2"])
    row = (
        f"{index}&{factor['search_id']}&{coefficient_text}&"
        f"{profile_text}\\\\"
    )
    require(row in MANUSCRIPT, f"factor-table row {index} does not match JSON")

for numerator in DATA["free_weight_numerators"]:
    require(str(numerator) in MANUSCRIPT,
            f"free-weight numerator {numerator} is absent from the manuscript")

free_numerators = DATA["free_weight_numerators"]
for offset in range(0, len(free_numerators), 4):
    keys = list(range(7 + offset, 7 + min(offset + 4, len(free_numerators))))
    values = free_numerators[offset:offset + 4]
    header = "$k$&" + "&".join(str(value) for value in keys)
    row = "$u_k$&" + "&".join(str(value) for value in values)
    header_position = COMPACT_MANUSCRIPT.find(header)
    row_position = COMPACT_MANUSCRIPT.find(row, header_position)
    require(header_position >= 0 and
            header_position < row_position < header_position + 100,
            f"free-weight table block beginning at k={keys[0]} does not match JSON")

clearing_display_text = "0.43339168265643340846553249867665"
require(clearing_display_text in MANUSCRIPT,
        "manuscript is missing the displayed 2-adic clearing exponent")
clearing_display = decimal_fraction(clearing_display_text)
clearing_exact = Fraction(DATA["clearing_alpha"])
clearing_unit = Fraction(1, 10 ** len(clearing_display_text.split(".")[1]))
require(clearing_display < clearing_exact < clearing_display + clearing_unit,
        "displayed 2-adic clearing exponent is not a valid truncation")

require(
    decimal_fraction(display["M1_lower"])
    < certified("M1_lower_numerator", "real_denominator"),
    "M1 lower decimal is not outward",
)
require(
    decimal_fraction(display["M1_upper"])
    > certified("M1_upper_numerator", "real_denominator"),
    "M1 upper decimal is not outward",
)
require(
    decimal_fraction(display["M2_lower"])
    < certified("M2_lower_numerator", "real_denominator"),
    "M2 lower decimal is not outward",
)
require(
    decimal_fraction(display["M2_upper"])
    > certified("M2_upper_numerator", "real_denominator"),
    "M2 upper decimal is not outward",
)
require(
    decimal_fraction(display["L_lower"])
    < certified("L_lower_numerator", "L_denominator"),
    "L lower decimal is not outward",
)
require(
    decimal_fraction(display["L_upper"])
    > certified("L_upper_numerator", "L_denominator"),
    "L upper decimal is not outward",
)
require(
    decimal_fraction(display["C_lower"])
    < certified("C_lower_numerator", "C_denominator"),
    "C lower decimal is not outward",
)
require(
    decimal_fraction(display["C_upper"])
    > certified("C_upper_numerator", "C_denominator"),
    "C upper decimal is not outward",
)
require(
    decimal_fraction(display["lambda_upper"])
    > certified("lambda_upper_numerator", "lambda_upper_denominator"),
    "lambda upper decimal is not outward",
)
require(
    decimal_fraction(display["tau_lower"])
    < certified("tau_lower_numerator", "tau_lower_denominator"),
    "tau lower decimal is not outward",
)

homogeneous_exact = Fraction(BOUNDS["homogeneous_upper"])
homogeneous_display = decimal_fraction(display["homogeneous_upper"])
target = Fraction(513_013, 125_000)
require(homogeneous_exact < homogeneous_display < target,
        "displayed homogeneous comparison is invalid")

real_hash = DATA["real_polynomial"]["coefficient_sha256"]
circle_hash = DATA["circle_polynomial"]["coefficient_sha256"]
require(real_hash in MANUSCRIPT, "real-polynomial hash mismatch")
require(circle_hash in MANUSCRIPT, "circle-polynomial hash mismatch")
require(DATA["real_polynomial"]["degree"] == 85, "unexpected degree of H")
require(DATA["real_polynomial"]["root_count_0_49"] == 63,
        "unexpected root count for H")
require(DATA["circle_polynomial"]["degree"] == 170, "unexpected degree of J")
require(DATA["circle_polynomial"]["quotient_degree"] == 169,
        "unexpected degree of K")
require(DATA["circle_polynomial"]["negative_bernstein_coefficients"] == 170,
        "unexpected Bernstein count")
require(DATA["claims"]["homogeneous_exponent"] == "4.104104",
        "JSON homogeneous target mismatch")
require(DATA["claims"]["irrationality_exponents"]["log(3)"] == "5.104104",
        "JSON log(3) target mismatch")
require(
    DATA["claims"]["irrationality_exponents"]["log(3)/log(2)"] == "5.104104",
    "JSON logarithm-ratio target mismatch",
)

print("PASS: manuscript decimals are outward and match the exact certificate")
print("PASS: factor rows and free-weight table agree with the exact certificate")
print("PASS: theorem targets, degrees, root counts, Bernstein count, and hashes agree")
