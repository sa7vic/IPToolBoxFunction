# `normxcorr2` — Normalised 2-D Cross-Correlation

**Scilab Image Processing Toolbox** · FOSSEE Summer Fellowship 2026 Screening Task  
Reimplementation of Octave `image` package · [`octave-forge/image: normxcorr2`](https://octave.sourceforge.io/image/function/normxcorr2.html)

---

## Table of Contents

1. [Overview](#overview)
2. [Mathematical Background](#mathematical-background)
3. [Files in This Folder](#files-in-this-folder)
4. [Calling Sequence](#calling-sequence)
5. [Parameters](#parameters)
6. [Return Value](#return-value)
7. [Error Conditions](#error-conditions)
8. [Implementation Notes](#implementation-notes)
9. [Test Cases](#test-cases)
10. [Running the Tests](#running-the-tests)
11. [Comparison with Octave](#comparison-with-octave)
12. [License](#license)

---

## Overview

`normxcorr2(T, A)` computes the **normalised cross-correlation** between a template matrix `T` and a search image `A`. The output is a correlation map `C` with values in `[-1, 1]`:

- `C(i,j) = 1` — perfect match of `T` at position `(i,j)` in `A`  
- `C(i,j) = 0` — no linear correlation  
- `C(i,j) = -1` — perfect inverse match  

Unlike raw cross-correlation, the normalised form is **invariant to additive and multiplicative intensity offsets**, making it robust to illumination changes between the template and the image.

---

## Mathematical Background

Let `T0 = T - mean(T)` be the zero-mean template. For each output position `(i,j)`, let `W_ij` denote the `m×n` window of `A` (zero-padded outside boundaries):

```
C(i,j) =     Σ  T0(p,q) · [W_ij(p,q) - mean(W_ij)]
          ─────────────────────────────────────────────────────
          ||T0|| · ||W_ij - mean(W_ij)||
```

This is the **Pearson correlation coefficient** between the zero-mean template and each local zero-mean image patch (Lewis 1995).

### Output Size

The output follows full cross-correlation convention:

```
rows(C) = rows(A) + rows(T) - 1
cols(C) = cols(A) + cols(T) - 1
```

---

## Files in This Folder

| File | Purpose |
|------|---------|
| `normxcorr2.sci` | Function implementation (FFT-based, O(N log N)) |
| `normxcorr2_test.sce` | Full test suite (24 tests, 6 sections + brute-force reference) |
| `README.md` | This document |

---

## Calling Sequence

```scilab
C = normxcorr2(T, A)
```

---

## Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `T` | Real 2-D numeric matrix | Template to search for. All real numeric types accepted. If `T` has zero variance (constant), `C` is all zeros. |
| `A` | Real 2-D numeric matrix | Image to search within. Any real numeric type. |

---

## Return Value

| Variable | Type | Description |
|----------|------|-------------|
| `C` | `double` matrix, size `(rows(A)+rows(T)-1) × (cols(A)+cols(T)-1)` | Normalised correlation map. All values clamped to `[-1, 1]`. |

---

## Error Conditions

| Error ID | Trigger |
|----------|---------|
| `normxcorr2:InvalidNumInputs` | Called with ≠ 2 arguments |
| `normxcorr2:ComplexInput` | `T` or `A` contains imaginary values |
| `normxcorr2:NotMatrix` | `T` or `A` has `ndims ≠ 2` |
| `normxcorr2:EmptyInput` | `T` or `A` has zero elements |

---

## Implementation Notes

### Critical Scilab 6.1.0 Compatibility Issues

Your code had two bugs that would crash on 6.1.0:

#### Bug 1: `size(x, "ndims")` is not valid in Scilab 6.1.0

```scilab
size(T, "ndims")    // WRONG — crashes on Scilab 6.1.0
ndims(T)            // CORRECT
```

`size()` only accepts `"r"`, `"c"`, or `"*"` as string arguments.

#### Bug 2: `fft(A, P, Q)` and `ifft()` do not exist in Scilab 6.1.0

Scilab 6.1.0's `fft` function does not accept a zero-padded 2-D calling form. The signature `fft(matrix, rows, cols)` is **not available**. Similarly, `ifft()` as a standalone function **does not exist** in 6.1.0.

The correct approach for 2-D FFT with zero-padding in Scilab 6.1.0:

```scilab
// Step 1: Zero-pad to (P x Q)
Xpad = zeros(P, Q);
Xpad(1:r, 1:c) = X;

// Step 2: 1-D FFT along each row
for i = 1:P
    Frow(i, :) = fft(Xpad(i, :));
end

// Step 3: 1-D FFT along each column
for j = 1:Q
    F(:, j) = fft(Frow(:, j));
end
```

For the inverse 2-D FFT, we use the identity `IFFT2(F) = conj(FFT2(conj(F))) / (P*Q)`:

```scilab
function X = _ifft2(F)
    [P, Q] = size(F);
    Fc = conj(F);
    // ... two-pass forward FFT on conj(F) ...
    X = conj(result) / (P * Q);
endfunction
```

This is mathematically exact and works on any Scilab 6.x version.

### FFT-Based Local Statistics

The local window mean and variance needed for normalisation are computed entirely in the frequency domain using the **integral image trick** (Lewis 1995). Convolving with a `ones(m,n)` kernel via FFT gives the running sum over every `m×n` window simultaneously:

```scilab
sumA  = real(_ifft2(FA   .* Fones))   // local sum of A
sumA2 = real(_ifft2(FA2  .* Fones))   // local sum of A^2

varA = sumA2 - (sumA .^ 2) / area     // local variance
stdA = sqrt(max(varA, 0))             // clamp to avoid sqrt(-eps) noise
```

This achieves `O(N log N)` complexity vs `O(N · m · n)` for a direct loop.

### Numerical Clamping

Floating-point arithmetic can produce values like `1.0000000000002` due to rounding in the FFT. The output is clamped to `[-1, 1]` to prevent such overshoot:

```scilab
C = min(max(C, -1), 1);
```

---

## Test Cases

The suite contains **24 test cases** across 6 sections, plus a pure Scilab brute-force reference implementation used as ground truth.

---

### Section 1 — Output Range

**TC-01–02 · Random inputs: all values in [-1, 1]**

```scilab
A = rand(10, 12);
T = rand(3, 4);
C = normxcorr2(T, A);
// min(C) >= -1  and  max(C) <= 1
```

---

### Section 2 — Output Size

**TC-04 · Size formula**

```scilab
A = rand(6, 7);  T = rand(2, 3);
C = normxcorr2(T, A);
size(C)   // Expected: [7, 9]  = [6+2-1, 7+3-1]
```

**TC-06 · Scalar template**

```scilab
C = normxcorr2([2], rand(5,6));
size(C)   // Expected: [5, 6]  (1x1 template → same size as A)
```

---

### Section 3 — Degenerate Inputs

**TC-07 · Constant template → all zeros**

```scilab
normxcorr2(ones(3,3)*5, rand(8,8))
// A constant template has zero variance => no correlation possible.
// Expected: C is all zeros.
```

**TC-08 · Constant image → all zeros**

```scilab
normxcorr2(rand(3,3), ones(8,8)*7)
// Every window has zero std => denominator = 0 everywhere.
// Expected: C is all zeros.
```

---

### Section 4 — Known Exact Values (FFT vs Brute Force)

**TC-12 · Template = subimage of A → C = 1 at match location**

```scilab
A = [0 0 0 0 0;
     0 1 2 0 0;
     0 3 4 0 0;
     0 0 0 0 0];
T = [1 2; 3 4];      // exact subimage
C = normxcorr2(T, A);
C(3, 3)   // Expected: 1.0  (perfect positive match)
```

**TC-13 · Negated template → C = -1 at match location**

```scilab
C = normxcorr2(-T, A);
C(3, 3)   // Expected: -1.0  (perfect inverse match)
```

**TC-10–11 · FFT result matches brute-force reference**

```scilab
// A 5x5 test image and 3x3 template
max(abs(normxcorr2(T,A) - brute_force_normxcorr2(T,A)))
// Expected: < 1e-8
```

---

### Section 5 — Symmetry and Invariance Properties

**TC-14 · Scale invariance: `normxcorr2(T, k*A) == normxcorr2(T, A)`**

```scilab
max(abs(normxcorr2(T, A) - normxcorr2(T, 5*A)))   // Expected: ~0
```

**TC-15 · Offset invariance: `normxcorr2(T, A+c) == normxcorr2(T, A)`**

```scilab
max(abs(normxcorr2(T, A) - normxcorr2(T, A+100)))  // Expected: ~0
```

**TC-16 · Template scale: `normxcorr2(k*T, A) == normxcorr2(T, A)`**

```scilab
max(abs(normxcorr2(T, A) - normxcorr2(7*T, A)))    // Expected: ~0
```

---

### Section 6 — Error Paths

| TC | Command | Expected |
|----|---------|----------|
| TC-17 | `normxcorr2(rand(3,3))` | Error: 1 argument |
| TC-19 | `normxcorr2([1+%i 2;3 4], rand(5,5))` | Error: complex T |
| TC-21 | `normxcorr2(ones(2,2,2), rand(5,5))` | Error: 3-D T |
| TC-23 | `normxcorr2([], rand(5,5))` | Error: empty T |

---

## Running the Tests

```scilab
exec("normxcorr2_test.sce", -1)
```

Expected final output:

```
----------------------------------------------------------------
RESULT: normxcorr2 test suite PASSED  (24 tests)
----------------------------------------------------------------
```

---

## Comparison with Octave

| Scenario | Octave result | Scilab `normxcorr2` result |
|----------|---------------|---------------------------|
| Constant template | all zeros | all zeros ✓ |
| Constant image | all zeros | all zeros ✓ |
| Template = subimage, C at peak | `1` | `1` ✓ |
| Output size `(M,N)+(m,n)-1` | correct | correct ✓ |
| Scale/offset invariance | invariant | invariant ✓ |
| Complex input | error | error ✓ |

---

## License

```
Copyright (C) 2026 – FOSSEE Summer Fellowship Submission
Released under the GNU General Public License v3 or later.
https://www.gnu.org/licenses/
```

*Submitted as part of the FOSSEE Summer Fellowship 2026 Screening Task:  
**Image Processing Toolbox Development** — IIT Bombay.*