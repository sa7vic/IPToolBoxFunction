# `otsuthresh` — Otsu's Global Binarisation Threshold

**Scilab Image Processing Toolbox** · FOSSEE Summer Fellowship 2026 Screening Task  
Reimplementation of Octave `image` package · [`octave-forge/image: otsuthresh`](https://octave.sourceforge.io/image/function/otsuthresh.html)

---

## Table of Contents

1. [Overview](#overview)
2. [Mathematical Background](#mathematical-background)
3. [Files in This Folder](#files-in-this-folder)
4. [Calling Sequence](#calling-sequence)
5. [Parameters](#parameters)
6. [Return Value](#return-value)
7. [Supported Data Types](#supported-data-types)
8. [Error Conditions](#error-conditions)
9. [Implementation Notes](#implementation-notes)
10. [Test Cases](#test-cases)
11. [Running the Tests](#running-the-tests)
12. [Comparison with Octave](#comparison-with-octave)
13. [License](#license)

---

## Overview

`otsuthresh(x)` computes the globally optimal greyscale threshold using **Otsu's method** (1979). The threshold divides pixels into two classes — background and foreground — such that the **between-class variance is maximised**, or equivalently, the weighted within-class variance is minimised.

The function accepts either a pre-computed histogram vector or a raw 2-D grayscale image, and returns a normalised threshold `t ∈ [0, 1]`.

Typical usage with `im2bw`:

```scilab
exec("otsuthresh.sci", -1);
I  = imread("coins.png");       // requires IPCV toolbox for imread
t  = otsuthresh(I);
BW = I > uint8(t * 255);        // binarise
```

---

## Mathematical Background

Given a normalised histogram `p(i)` for `i = 0, 1, ..., L-1` (L grey levels):

```
omega_0(k) = sum_{i=0}^{k}   p(i)          [background cumulative probability]
omega_1(k) = 1 - omega_0(k)                [foreground cumulative probability]

mu(k)      = sum_{i=0}^{k}   i * p(i)      [cumulative mean]
mu_T       = mu(L-1)                        [total mean]

sigma_b^2(k) = [mu_T * omega_0(k) - mu(k)]^2
               ---------------------------------
               omega_0(k) * [1 - omega_0(k)]
```

The optimal threshold is:

```
k* = argmax_k  sigma_b^2(k)
t  = k* / (L - 1)               [normalised to [0, 1]]
```

In case of ties, the **lowest** k is chosen for deterministic behaviour.

---

## Files in This Folder

| File | Purpose |
|------|---------|
| `otsuthresh.sci` | Function implementation |
| `otsuthresh_test.sce` | Full test suite (25 tests, 5 sections) |
| `README.md` | This document |

---

## Calling Sequence

```scilab
t = otsuthresh(I)        // 2-D grayscale image
t = otsuthresh(counts)   // pre-computed histogram count vector
```

---

## Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `I` | Real 2-D numeric matrix | Grayscale image. Any numeric type. Pixel values are linearly scaled to `[0,1]` internally before histogramming. |
| `counts` | Real numeric vector (1-D) | Histogram bin counts. Length ≥ 2. All values ≥ 0. Sum > 0. |

---

## Return Value

| Variable | Type | Description |
|----------|------|-------------|
| `t` | `double` scalar | Normalised threshold in `[0, 1]`. Multiply by `255` for uint8 threshold, or by `(L-1)` for a bin index. |

---

## Supported Data Types

| Input Type | Notes |
|------------|-------|
| `double` | Direct; no conversion overhead |
| `uint8` | Cast to double; values `[0, 255]` scaled to `[0, 1]` |
| `uint16` | Cast to double; values scaled via min/max |
| `int16` | Signed integers correctly handled via min/max scaling |
| `int8`, `int32`, `uint32` | All cast to `double` safely |
| Histogram vector | Any real type; sum must be > 0 |

---

## Error Conditions

| Error ID | Trigger |
|----------|---------|
| `otsuthresh:InvalidNumInputs` | Called with ≠ 1 argument |
| `otsuthresh:ComplexInput` | Input contains imaginary component |
| `otsuthresh:NotVectorOrMatrix` | Input has `ndims > 2` (e.g. 3-D array) |
| `otsuthresh:TooFewBins` | Histogram vector has fewer than 2 elements |
| `otsuthresh:NegativeCounts` | Any histogram bin is negative |
| `otsuthresh:ZeroHistogram` | Sum of histogram is ≤ 0 |

---

## Implementation Notes

### Scilab 6.1.0 Compatibility: `ndims()` vs `size(x, "ndims")`

A critical Scilab 6.1.0 portability issue: `size(x, "ndims")` is **not a valid call** in Scilab 6.1.0. The `size()` function only accepts `"r"`, `"c"`, or `"*"` as string arguments. To query the number of dimensions, the correct call is:

```scilab
ndims(x)    // CORRECT in Scilab 6.1.0
size(x, "ndims")    // WRONG - runtime error in Scilab 6.1.0
```

This implementation uses `ndims(x)` throughout.

### Vectorised Otsu — No for-loop

The between-class variance is computed entirely with vectorised operations using `cumsum`. This is faster and more readable than a per-bin loop:

```scilab
p      = counts / total;          // probabilities
omega  = cumsum(p);               // cumulative class prob
mu     = cumsum(p .* levels);     // cumulative mean
num_b2 = (mu_T .* omega - mu).^2;
den_b2 = omega .* (1 - omega);
sigma_b2(valid) = num_b2(valid) ./ den_b2(valid);
```

Degenerate bins (where `omega = 0` or `omega = 1`) are assigned `-%inf` so they cannot win the `max`, without any conditional branching.

### Image Mode: Histogram Recycling

Image mode builds a 256-bin histogram from the scaled pixel values, then calls `otsuthresh(counts)` recursively. This ensures complete consistency between the two calling modes and avoids code duplication.

---

## Test Cases

The test suite (`otsuthresh_test.sce`) contains **25 test cases** across 5 sections.

---

### Section 1 — Histogram Mode: Range

**TC-01 · Classic bimodal histogram**

```scilab
counts = zeros(256, 1);
counts(1)   = 1000;     // dark class
counts(256) = 1000;     // bright class
t = otsuthresh(counts);
// Expected: t in [0, 1]
```

Two equal-weight peaks at extremes — Otsu threshold should land near the midpoint.

**TC-04 · Row-vector histogram**

```scilab
counts = zeros(1, 64);
counts(1) = 200; counts(64) = 200;
t = otsuthresh(counts);   // row vector, not column
// Expected: t in [0, 1]  (both orientations must work)
```

---

### Section 2 — Histogram Mode: Known Exact Values

**TC-05 · Two-bin histogram → t = 0**

```scilab
otsuthresh([100; 100])
// Only one candidate threshold: k=1.
// omega_0=0.5, denom=0.25 => max is at k=1.
// t = (1-1)/(2-1) = 0
```

**TC-09 · Normalisation check**

```scilab
counts = [10; 10; 10; 970; 0];  // 5-bin histogram
t = otsuthresh(counts);
// t must be one of {0, 0.25, 0.5, 0.75, 1.0}
// because t = (k-1)/(nbins-1) = (k-1)/4
```

---

### Section 3 — Image Mode

**TC-11 · Uniform image → 0**

```scilab
otsuthresh(ones(8,8) * 7)    // all pixels identical
// Expected: t = 0  (deterministic, defined behaviour)
```

**TC-13 · uint8 bimodal image**

```scilab
I = uint8([zeros(1,128), 255*ones(1,128)]);
t = otsuthresh(I);
// Expected: t in [0, 1]  (overflow-safe uint8 handling)
```

**TC-15 · int16 with negative pixels**

```scilab
I = int16([-100*ones(5,5), 100*ones(5,5)]);
t = otsuthresh(I);
// Expected: t in [0, 1]  (signed int correctly scaled)
```

---

### Section 4 — Image/Histogram Consistency

**TC-18 · Image mode == manually-built histogram**

```scilab
I = [zeros(10,10); ones(10,10)*0.8];
// Build histogram manually ...
t_img  = otsuthresh(I);
t_hist = otsuthresh(h);       // h = manually built 256-bin histogram
// Expected: t_img == t_hist  (to floating-point precision)
```

This is the key consistency test: both calling modes must agree.

---

### Section 5 — Error Paths

| TC | Command | Expected |
|----|---------|----------|
| TC-19 | `otsuthresh()` | Error: wrong number of inputs |
| TC-21 | `otsuthresh([1+%i 2])` | Error: complex input |
| TC-22 | `otsuthresh([100])` | Error: fewer than 2 bins |
| TC-23 | `otsuthresh([-1; 10; 5])` | Error: negative counts |
| TC-24 | `otsuthresh(zeros(10,1))` | Error: zero-sum histogram |
| TC-25 | `otsuthresh(ones(3,3,3))` | Error: 3-D array |

---

## Running the Tests

```scilab
exec("otsuthresh_test.sce", -1)
```

Expected final output:

```
----------------------------------------------------------------
RESULT: otsuthresh test suite PASSED  (25 tests)
----------------------------------------------------------------
```

---

## Comparison with Octave

| Scenario | Octave result | Scilab `otsuthresh` result |
|----------|---------------|---------------------------|
| `otsuthresh([100;100])` | `0` | `0` ✓ |
| Bimodal uint8 image (0 and 255) | `~0.5` | `~0.5` ✓ |
| Uniform image | `0` | `0` ✓ |
| Empty histogram | error | error ✓ |
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