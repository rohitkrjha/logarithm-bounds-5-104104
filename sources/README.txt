Improved irrationality-exponent bounds for log 3 and log₂ 3
=========================================================

Author: Rohit Kumar Jha
Manuscript/submission date: 9 August 2026
Archive version: v1

This is the author's manuscript submitted to the Journal of Number Theory
on 9 August 2026. In this public copy, the author block contains only the
author's name: the affiliation, address, and institutional correspondence
details have been removed. No ORCID is supplied. The approximate-weight
table uses spacing rather than equality signs between its index and weight
columns. The elementary Laplace lemma explicitly requires a nondegenerate
interval, as satisfied by both intervals used in the proof. The main results,
proof construction, and exact computational supplement are unchanged.
The submission date is
distinct from the date of public posting.

RESULT

  H = max(|h₁|, |h₂|)

  |h₀ + h₁ log 2 + h₂ log 3| > H^(−4.104104)

for all sufficiently large H and integer coefficients with (h₁,h₂) ≠ (0,0).
Consequently,

  μ(log 3) ≤ 5.104104,
  μ(log 3 / log 2) ≤ 5.104104.

FILES

  logarithm-bounds-5-104104-v1.pdf
      The 21-page public preprint, with the author details above.
  logarithm-bounds-5-104104-v1-sources.zip
      Manuscript sources, exact certificate, all three checkers, pinned
      requirements, reference outputs, and source-provenance records.
  README.txt
      These reproduction instructions.
  RIGHTS.txt
      Component-specific licensing.
  SHA256SUMS.txt
      SHA-256 hashes of the other four deposited files.

MANUSCRIPT BUILD

Extract the source archive into an empty directory. With a TeX distribution
containing elsarticle and latexmk installed, run:

  latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex

main.tex is a build entry point that inputs manuscript.tex.
The bibliography is contained in manuscript.tex.

EXACT VERIFICATION

Tested with Python 3.12.4, SymPy 1.14.0, and python-flint 0.9.0. From the
extracted source directory, create an isolated Python environment:

  python3.12 -m venv .venv
  . .venv/bin/activate
  python -m pip install -r supplement/requirements-audit.txt
  cd supplement
  python verify_log3_bound_5_104104.py --json reproduced_certificate.json
  cmp certificate_log3_5_104104.json reproduced_certificate.json
  python verify_manuscript_constants.py
  python independent_exact_audit.py
  python audit_flint_crosscheck.py

The primary check needs SymPy; the second exact implementation uses only
Python's standard library. The third check uses FLINT/Arb. Once installed,
the checks themselves require no network access. Do not use python -O or
set PYTHONOPTIMIZE: the exact checkers deliberately refuse optimized Python.
Expected outputs and detailed instructions are in supplement/README.md.

VERIFICATION SCOPE

The executable checks certify the finite algebraic and rational-arithmetic
claims: 21 exact weights and local content profiles, 63 isolated real
critical roots (split 41/22), exact logarithm enclosures, Bernstein signs
for the circle bound, and the final rational comparison. The regenerated
certificate is compared byte-for-byte with the submitted certificate.

The analytic and Diophantine arguments connecting these finite facts to
the theorem are proved in the manuscript. This package is not an
end-to-end Lean formalization; incomplete Lean development notes from the
older working directory are not included.

INTEGRITY AND RIGHTS

From the directory containing the five deposited files, on macOS:

  shasum -a 256 -c SHA256SUMS.txt

On Linux:

  sha256sum -c SHA256SUMS.txt

The extracted sources contain a complete root SHA256SUMS.txt and the
original supplement/SHA256SUMS.txt. Check the latter from supplement.
The exploratory transcript is retained because the manuscript cites it
as provenance; its numerical design value is not the proved bound.
See RIGHTS.txt for CC BY 4.0 manuscript/documentation and Apache 2.0 code.
