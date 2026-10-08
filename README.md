# Improved irrationality-exponent bounds for log 3 and log₂ 3

Rohit Kumar Jha · Manuscript submitted to the Journal of Number Theory on **9 August 2026**

[Read the paper](release/logarithm-bounds-5-104104-v1.pdf) · [LaTeX source](sources/manuscript.tex) · [Exact verification package](sources/supplement) · [Citation](CITATION.cff)

## Results

For integer coefficients with (h₁,h₂) ≠ (0,0), and all sufficiently large
H = max(|h₁|,|h₂|),

$$
|h_0+h_1\log 2+h_2\log 3|>H^{-4.104104}.
$$

Consequently,

$$
\mu(\log 3)\leq 5.104104,
\qquad
\mu\!\left(\frac{\log 3}{\log 2}\right)\leq 5.104104.
$$

This repository preserves the mathematics of the August 2026 manuscript and
its exact certificate. The public copy uses a name-only author block, without
affiliation, institutional contact details, or ORCID.
The submission date records the manuscript's history; it is not a backdated
Git commit or public-release timestamp.

## Reproduction

From the repository root:

    python3 scripts/check_integrity.py
    cd sources
    latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex

To reproduce the computational certificate, start from the repository root
with Python 3.12 available:

    python3.12 -m venv .venv
    . .venv/bin/activate
    python -m pip install -r sources/supplement/requirements-audit.txt
    cd sources/supplement
    python verify_log3_bound_5_104104.py --json reproduced_certificate.json
    cmp certificate_log3_5_104104.json reproduced_certificate.json
    python verify_manuscript_constants.py
    python independent_exact_audit.py
    python audit_flint_crosscheck.py

The dependencies are pinned to SymPy 1.14.0 and python-flint 0.9.0.
The independent exact checker needs only Python's standard library.
Run without Python optimization flags. Once dependencies are installed,
the computations are offline.

## Verification

The primary verifier, independent exact checker, and FLINT/Arb root audit
were rerun successfully. The regenerated certificate was byte-identical
to the frozen certificate. The manuscript-constant checks passed, and a
fresh TeX build of the public copy completed without warnings.

The software checks the finite polynomial, root, sign, interval, and rational
comparison claims. The analytic and Diophantine lemmas connecting those claims
to the theorem are in the manuscript. This is not an end-to-end Lean proof.
GitHub Actions checks archive integrity and manuscript/certificate consistency,
not a fresh run of all computational audits.

## Files and citation

The release directory contains the five files for the corresponding Zenodo
preprint deposit. The sources directory contains the exact unpacked source
archive. The submitted supplement and proof constants are preserved. The
public edits remove the author-contact block and replace equality-sign
separators in the approximate-weight table with spacing. Private journal-portal records,
cover letters, and the old author-contact block are excluded.

Use [CITATION.cff](CITATION.cff) or GitHub's **Cite this repository** control.
The Zenodo DOI will be added when available. No acceptance, journal DOI,
volume, issue, or journal publication date is claimed.

## Rights

Manuscript, certificate data, and documentation: **CC BY 4.0**.
Original verification/audit code and build configuration: **Apache 2.0**.
Third-party software retains its own licenses.
See [LICENSE](LICENSE) and [RIGHTS.txt](sources/RIGHTS.txt).
