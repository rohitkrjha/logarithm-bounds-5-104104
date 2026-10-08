# From the numerical 5.10407658 design to the exact 5.104104 theorem

## 1. Original numerical design

The surviving exploratory LP transcript gave

```text
mu = 5.104076581512609
C  = 1.3003307431474953
M  = -2.722890642820573
L  = 4.537964026899708
```

with 21 active factors. This was a discretized semi-infinite linear-programming result, not a theorem.

## 2. The original coefficient-growth concern

At the numerically optimized contour radius, the modulus did not have a single
real dominant saddle. It had two conjugate nonreal maxima. The original
proof-candidate route sought a full coefficient-growth limit and therefore
had to address possible cancellation between their leading contributions.
The final paper instead proves and uses a limit-superior two-form criterion;
for that argument, Cauchy's upper bound is exactly what is required.

## 3. Rational saddle repair

A rational point on the circle was selected:

\[
z_0=26\frac{-899+60i}{901}.
\]

Two exact linear constraints were imposed on the weights so that the full complex logarithmic derivative vanished at `z0`. Reoptimizing under this constraint barely changed the design.

## 4. Deliberate separation of the real rates

The numerical optimization was then constrained so that the maximum on `[25,49]` lay at least `2e-5` below the maximum on `[0,25]`. This gives a robust exact verification of the ratio condition in Hata's two-form lemma.

## 5. Exact rationalization

The final fifteen small weights were fixed with denominator `10^12`. The first six weights were solved exactly from:

- weighted degree `3/2`;
- exact weighted `7`-, `5`-, and `3`-adic profile balances;
- real and imaginary parts of the exact complex saddle equation.

The resulting common denominator has 1061 digits. This is asymptotically harmless: the sequence is taken over multiples of twice that denominator.

## 6. Exact real certificate

The rationalized weights produce a degree-85 critical polynomial. Exact VAS isolation finds 63 roots in `[0,49]`. Rational interval arithmetic gives

```text
M1 = -2.72295749566841...
M2 = -2.72297756238961...
M1-M2 > 2.0066721204e-5.
```

## 7. Exact circle certificate

At radius 26, the degree-170 derivative numerator factors as

\[
J(c)=(901c+899)K(c).
\]

Every one of the 170 Bernstein coefficients of `K(2x-1)` is negative; every Bernstein coefficient of every circle factor is positive. Thus the only circle maxima are the conjugate points corresponding to `c=-899/901`, and they are nondegenerate.

## 8. Historical noncancellation check and final simplification

The original package also checked, using the Gaussian prime 1+2i, that the
relative saddle phase is not a root of unity. This can produce a
positive-density noncancelling subsequence and a full growth limit. That
additional argument is not used in the JNT manuscript. The exact circle
maximum and Cauchy's inequality give the limit-superior bound directly, while
the distinct real rates supply the cancellation control needed by the
self-contained two-form lemma.

## 9. Final exact bound

The exact certificate gives

```text
C       < 1.3004042229114373661
L       < 4.5379010755984064050
lambda  < 5.8383052985098437711
tau_1   > 1.4225532727569723450
1+lambda/tau_1 < 5.1041031013165116971.
```

The theorem-facing outward-rounded value is therefore

\[
|h_0+h_1\log2+h_2\log3|>H^{-4.104104}
\]

for all sufficiently large \(H\), where
\(H=\max(|h_1|,|h_2|)\), and consequently

\[
\mu(\log3)\le 5.104104,
\qquad
\mu(\log3/\log2)\le 5.104104.
\]
