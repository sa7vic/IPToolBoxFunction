# `immse` — Mean Squared Error for Images & Arrays

**Scilab Image Processing Toolbox** · FOSSEE Summer Fellowship 2026 Screening Task  
Reimplementation of Octave `image` package · [`octave-forge/image: immse`](https://octave.sourceforge.io/image/function/immse.html)

---

## Table of Contents

1. [Overview](#overview)
2. [Mathematical Definition](#mathematical-definition)
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

`immse(A, B)` computes the **Mean Squared Error (MSE)** between two numeric arrays or images of identical dimensions.  
MSE is one of the most widely used full-reference image quality metrics. It quantifies the average per-pixel squared difference between a reference image and a distorted (compressed, noisy, or reconstructed) image.

A lower MSE indicates higher similarity between two images. When `A` and `B` are identical, `immse` returns exactly `0`.

---

## Mathematical Definition

Given two arrays **A** and **B** with *N* total elements:

```
MSE = (1/N) · Σᵢ (Aᵢ − Bᵢ)²   for i = 1 … N
```

Both arrays are cast to `double` before the subtraction, preventing silent integer-arithmetic overflow that would otherwise occur with `uint8`, `uint16`, or `int16` pixel data.

---

## Files in This Folder

| File | Purpose |
|------|---------|
| `immse.sci` | Function implementation — source to exec or load |
| `immse_test.sce` | Full test suite (27 tests, 6 sections) |
| `README.md` | This document |

---

## Calling Sequence

```scilab
mse = immse(A, B)
```

---

## Parameters

| Parameter | Type | Size | Description |
|-----------|------|------|-------------|
| `A` | Real numeric | Any | First image or array (reference). Accepts `double`, `single`, `int8`, `int16`, `int32`, `uint8`, `uint16`, `uint32`. |
| `B` | Real numeric | Same as `A` | Second image or array (distorted / reconstructed). Must match `A` in every dimension. |

> **Note on 3-D arrays:** For an RGB image of shape `H × W × 3`, pass the full 3-D matrix directly. `immse` flattens both arrays before computing the mean, so channel weighting is uniform — consistent with Octave's behaviour.

---

## Return Value

| Variable | Type | Description |
|----------|------|-------------|
| `mse` | `double` scalar | Mean squared error. Always `≥ 0`. Returns `0` when `A == B` element-wise. |

---

## Supported Data Types

| Type | Notes |
|------|-------|
| `double` | Native; no conversion overhead |
| `single` | Upcast to `double` before arithmetic |
| `uint8` | Typical 8-bit image; cast prevents `0 − 200 → 56` wrap-around |
| `uint16` | 16-bit medical / HDR images |
| `int8`, `int16`, `int32`, `uint32` | All cast to `double` safely |

Complex-valued arrays are **not** supported and will raise an error.

---

## Error Conditions

| Error ID | Trigger | Message |
|----------|---------|---------|
| `immse:InvalidNumInputs` | Called with ≠ 2 arguments | `immse requires exactly 2 inputs: immse(A, B)` |
| `immse:EmptyInput` | Either `A` or `B` is empty (`[]`) | `A and B must be non-empty arrays` |
| `immse:ComplexNotSupported` | Either input has an imaginary component | `A and B must be real numeric arrays` |
| `immse:SizeMismatch` | `size(A) ≠ size(B)` | Reports actual sizes of both arguments |

---

## Implementation Notes

### Why cast to `double` first?

Scilab (like MATLAB) performs arithmetic on integer types *in the integer domain*, meaning subtraction of `uint8` values saturates at 0 rather than wrapping or going negative. For example:

```scilab
uint8(10) - uint8(200)   // returns uint8(0), NOT -190 or 66
```

This silently corrupts the MSE calculation. By calling `double(A)` and `double(B)` before subtraction, we guarantee exact signed arithmetic across the full dynamic range of any pixel format.

### Why flatten with `D(:)`?

Scilab's `mean()` applied to a matrix computes column-wise means by default (returning a row vector, not a scalar). Flattening to a column vector with `D(:)` ensures `mean()` always returns a single scalar, regardless of input dimensionality — essential for 3-D image arrays.

```scilab
D  = Ad(:) - Bd(:);   // 1-D column vector
mse = mean(D .^ 2);   // guaranteed scalar
```

### Vectorised, no loops

The implementation uses Scilab's native vectorised operators throughout. There are no `for` or `while` loops, making it efficient even for high-resolution images.

---

## Test Cases

The test suite (`immse_test.sce`) contains **27 test cases** across six sections. The most important ones are summarised below.

---

### Section 1 — Basic Correctness

**TC-01 · Identical arrays → 0**

```scilab
A = [1 2; 3 4];
immse(A, A)          // Expected: 0
```

Validates the zero-path: no rounding should occur when inputs are identical.

---

**TC-02 · Known 2×2 example → 1.5**

```scilab
A = [1 2; 3 4];
B = [2 2; 2 2];
// D = A−B = [-1 0; 1 2]
// D² = [1 0; 1 4], mean(D²) = 6/4 = 1.5
immse(A, B)          // Expected: 1.5
```

---

**TC-03 · Zeros vs ones → 1.0**

```scilab
immse(zeros(2,2), ones(2,2))   // Expected: 1.0
```

---

**TC-04 · Scalar inputs**

```scilab
immse(3, 7)    // D=−4, D²=16  → Expected: 16.0
```

---

**TC-05 · Row vector → 1.5**

```scilab
immse([0 1 2 3], [1 1 1 1])
// D=[-1 0 1 2], D²=[1 0 1 4], mean=1.5  → Expected: 1.5
```

---

### Section 2 — Integer-Type Safety

**TC-09 · uint8 — overflow prevention → 350**

```scilab
A8 = uint8([0 10; 20 30]);
B8 = uint8([0  0;  0  0]);
// Without double cast: Scilab integer arithmetic would be wrong.
// D²=[0 100; 400 900], mean=1400/4=350
immse(A8, B8)    // Expected: 350
```

**TC-10 · uint8 near-max values → 48062.5**

```scilab
A8 = uint8([200 250]);
B8 = uint8([ 10   5]);
// D=[190 245], D²=[36100 60025], mean=48062.5
immse(A8, B8)    // Expected: 48062.5
```

**TC-11 · uint16 → 2500000**

```scilab
immse(uint16([1000 2000]), uint16([0 0]))
// D=[1000 2000], D²=[1e6 4e6], mean=2.5e6 → Expected: 2500000
```

---

### Section 3 — Multi-Dimensional Arrays

**TC-13 · 3-D ones vs zeros → 1.0**

```scilab
immse(ones(2,2,3), zeros(2,2,3))
// 12 elements, all D=1, MSE=1  → Expected: 1.0
```

**TC-14 · 3-D RGB-like, mixed channels → 5/3**

```scilab
A3 = cat(3, zeros(2,2), ones(2,2),     zeros(2,2));
B3 = cat(3, zeros(2,2), zeros(2,2), 2*ones(2,2));
// Ch1 all-0, Ch2 all-1, Ch3 all-4 → (4×0 + 4×1 + 4×4)/12 = 20/12 = 5/3
immse(A3, B3)    // Expected: 1.6667 (5/3)
```

---

### Section 4 — Symmetry

**TC-16 · immse(A,B) == immse(B,A)**

```scilab
A = [1 3 5; 7 9 11];
B = [2 2 6; 6 10 10];
immse(A,B) == immse(B,A)   // Must be %T
```

---

### Section 5 — Numerical Precision

**TC-20 · Fractional inputs → 0.065**

```scilab
immse([0.1 0.4], [0.3 0.1])
// D=[-0.2 0.3], D²=[0.04 0.09], mean=0.065
```

---

### Section 6 — Error Paths (all must throw)

| TC | Command | Expected |
|----|---------|----------|
| TC-21 | `immse([1 2 3])` | Error: wrong number of inputs |
| TC-23 | `immse([1 2],[1 2 3])` | Error: size mismatch |
| TC-24 | `immse(ones(2,2),ones(4,1))` | Error: shape mismatch |
| TC-25 | `immse([],[])` | Error: empty input |
| TC-26 | `immse([1+2*%i 3],[1 3])` | Error: complex not supported |

---

## Running the Tests

1. Open Scilab 6.1.0 or later.
2. `cd` to the folder containing `immse.sci` and `immse_test.sce`.
3. Execute:

```scilab
exec("immse_test.sce", -1)
```

Expected final output:

```
--------------------------------------------------------------------
RESULT: immse test suite PASSED  (27 tests)
--------------------------------------------------------------------
```

If any test fails, the suite halts immediately and prints a descriptive `FAIL:` message indicating which test case failed and the actual vs expected values.

---

## Comparison with Octave

The table below confirms functional parity with the Octave `image` package `immse`.

| Scenario | Octave result | Scilab `immse` result |
|----------|---------------|-----------------------|
| `immse([1 2;3 4],[2 2;2 2])` | `1.5000` | `1.5` ✓ |
| `immse(uint8([0 10;20 30]),uint8(zeros(2)))` | `350` | `350` ✓ |
| `immse(ones(2,2,3),zeros(2,2,3))` | `1` | `1` ✓ |
| `immse([1 2],[1 2 3])` | error | error ✓ |
| `immse([],[])` | error | error ✓ |

---

## License

This implementation is released under the **GNU General Public License v3** (or any later version), consistent with the FOSSEE project's open-source licensing policy.

```
Copyright (C) 2026 – FOSSEE Summer Fellowship Submission
This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.
```

---

*Submitted as part of the FOSSEE Summer Fellowship 2026 Screening Task:  
**Image Processing Toolbox Development** — IIT Bombay.*
