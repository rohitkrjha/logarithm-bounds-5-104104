# Supplement for the 5.104104 logarithm bound

This directory accompanies “Improved irrationality-exponent bounds for log 3
and log_2 3.” Together with the analytic and Diophantine lemmas proved
in the manuscript, the exact finite certificate establishes

```text
|h0 + h1 log(2) + h2 log(3)| > H^(-4.104104)
mu(log(3)) <= 5.104104
mu(log(3)/log(2)) <= 5.104104.
```

The software certifies the finite polynomial, root, interval, and rational
comparison claims. The bridge from those claims to the theorem is proved in
the manuscript.

## Files

- `verify_log3_bound_5_104104.py`: primary exact verifier.
- `certificate_log3_5_104104.json`: machine-readable exact certificate.
- `verify_manuscript_constants.py`: dependency-free consistency checker for
  the outward decimal bounds, factor table, targets, hashes, and JSON data.
- `independent_exact_audit.py`: CAS-free, standard-library reconstruction of
  the finite certificate, including exact Descartes root coverage.
- `audit_flint_crosscheck.py`: optional independent FLINT/Arb root audit.
- `EXPECTED_INDEPENDENT_AUDIT_OUTPUT.txt`: reference output of the CAS-free
  audit.
- `EXPECTED_VERIFIER_OUTPUT.txt`: stable theorem-facing output of the primary
  verifier (timing checkpoints and the chosen output path are omitted).
- `EXPECTED_FLINT_AUDIT_OUTPUT.txt`: stable concluding output of the FLINT/Arb
  audit (the imported primary-verifier output and progress lines are omitted).
- `NUMERICAL_TO_EXACT_HISTORY.md`: retrospective numerical-to-exact history.
- `exploratory_original_5_10407658.txt`: surviving numerical design output;
  it is provenance only and is not part of the proof.
- `requirements.txt`: primary verifier dependency.
- `requirements-audit.txt`: optional audit dependencies.
- `SHA256SUMS.txt`: hashes for the frozen supplement.

## Primary verification

Tested with SymPy 1.14.0 under Python 3.10.14 and Python 3.12.4. The supplied
reference package was also generated under Python 3.13.5.

```bash
python -m pip install -r requirements.txt
python verify_log3_bound_5_104104.py --json reproduced_certificate.json
cmp certificate_log3_5_104104.json reproduced_certificate.json
python verify_manuscript_constants.py
```

The last command compares the exact JSON with the submitted manuscript source.
It expects `manuscript.tex` one directory above this supplementary directory.
If the files were downloaded elsewhere, use
`python verify_manuscript_constants.py --manuscript /path/to/manuscript.tex`.

The verifier uses exact integers and rationals, exact VAS real-root
isolation, exact Bernstein coefficients, and positive rational logarithm
series with explicit tails. Decimal text is diagnostic only. The program
refuses to run under `python -O` or with `PYTHONOPTIMIZE`, because those modes
would disable its assertion checks.

## CAS-free independent reproduction

```bash
python independent_exact_audit.py
```

This second exact implementation imports no CAS or third-party package. It
independently rebuilds the weights, coefficient profiles, degree-85 and
degree-170 polynomials, root coverage, interval bounds, Bernstein signs, and
final integer comparison. For the real critical polynomial it verifies one
root in each of the 63 exported intervals and no roots in all 64 complementary
gaps by exact Descartes transformations. It also refuses optimized Python.

## Additional FLINT/Arb root audit

```bash
python -m pip install -r requirements-audit.txt
python audit_flint_crosscheck.py
```

This third program loads the exact integer polynomials reconstructed by the
primary verifier and passes them to FLINT/Arb. It independently confirms the
41/22 real-root split, one root in each of the 63 exported rational intervals,
and the absence of roots of the circle quotient and circle factors on the
forbidden interval.

## Proof scope

The publication proof uses only Cauchy's limit-superior bound for the common
coefficient. It does not require a conjugate-saddle asymptotic, a phase
nonperiodicity argument, or a selected noncancelling subsequence. No
optimality is claimed for the exploratory factor pool.
