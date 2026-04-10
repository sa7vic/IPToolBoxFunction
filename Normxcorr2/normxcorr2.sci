// =============================================================================
// normxcorr2.sci
// =============================================================================
//
// Copyright (C) 2026 - FOSSEE Summer Fellowship Submission
// This file is part of a Scilab Image Processing Toolbox implementation.
// Released under the GNU General Public License v3 or later.
// See <https://www.gnu.org/licenses/> for details.
//
// Function  : normxcorr2
// Purpose   : Normalised 2-D cross-correlation of a template with an image.
// Toolbox   : Image Processing (Scilab reimplementation of Octave 'image' pkg)
// Reference : Octave image package - normxcorr2()
//             https://octave.sourceforge.io/image/function/normxcorr2.html
//             J. P. Lewis, "Fast Normalized Cross-Correlation,"
//             Vision Interface, 1995.
//
// Author    : FOSSEE Screening Task - Image Processing Toolbox Development
// Created   : April 2026
// Scilab    : 6.1.0+  (tested on 6.1.0 64-bit, Windows & Linux)
//
// =============================================================================
//
// CALLING SEQUENCE
//   C = normxcorr2(T, A)
//
// DESCRIPTION
//   normxcorr2 computes the normalised cross-correlation between a template
//   matrix T and an image matrix A.  The output C has values in [-1, 1]:
//
//     C(i,j) =  1  => perfect positive match of T in A at position (i,j)
//     C(i,j) =  0  => no linear correlation
//     C(i,j) = -1  => perfect inverse match
//
//   The output size follows the "full" convolution convention:
//     rows(C) = rows(A) + rows(T) - 1
//     cols(C) = cols(A) + cols(T) - 1
//
//   The normalisation makes the metric invariant to both additive offsets
//   (brightness shifts) and multiplicative scaling of image intensities.
//
// ALGORITHM  (FFT-based, O(N log N))
//   Let T0 = T - mean(T)  (zero-mean template, L2-norm = denT).
//   For each output position (i,j), let W be the m*n image window
//   (zero-padded outside A boundaries):
//
//       C(i,j) = sum(T0 .* (W - mean(W)))
//                / ( denT * ||W - mean(W)|| )
//
//   Numerator:
//     num = conv2(A, rot180(T0), "full")  [via FFT]
//     Because sum(T0) = 0, adding any constant c to A leaves num unchanged:
//       sum((A+c).*T0) = sum(A.*T0) + c*sum(T0) = sum(A.*T0).
//
//   Denominator per window (variance term):
//     varA = sumA2 - sumA^2/area
//     sumA and sumA2 are obtained by convolving A and A^2 with ones(m,n)
//     via FFT.  In exact arithmetic, variance is offset-invariant.
//
// VARIANCE ROBUSTNESS  (key implementation note)
//   Computing varA = sumA2 - sumA^2/area involves catastrophic cancellation
//   when A has a large DC offset (e.g. rand()+100): both terms are ~100^2*area,
//   their difference is small, and FFT rounding noise ~eps*100^2*area can
//   dominate the result.
//
//   WRONG FIX: using eps_rel * sumA2 as threshold.
//     When DC is large, sumA2 is large, so this threshold is large and
//     INCORRECTLY suppresses valid (non-zero) variances.  This breaks
//     offset invariance (TC-15 failure).
//
//   CORRECT FIX: absolute threshold scaled by the image's own magnitude.
//     varA_tol = eps_machine^(2/3) * max(1, max(|Ad|)^2)
//              ~ 3.67e-11 * max(1, max(|Ad|)^2)
//
//     This threshold is:
//       - Safe for constant images:  FFT noise << varA_tol  => zeroed.
//       - Safe for offset images:    true varA >> varA_tol  => kept.
//       - Does NOT grow with DC, only with the image's absolute maximum.
//
// SCILAB 6.1.0 COMPATIBILITY
//   The following do NOT exist in Scilab 6.1.0 and are avoided:
//     - ifft()             use conjugate-FFT identity instead
//     - fft(A, P, Q)       zero-padded 2-D form does not exist
//     - size(x,"ndims")    use ndims(x) instead
//   2-D FFT is implemented via two sequential 1-D passes (row then column),
//   which is mathematically identical to the 2-D DFT.
//
// PARAMETERS
//   T  : 2-D real numeric matrix (template).  Non-empty.
//        Accepted types: double, uint8, uint16, int8, int16, int32, uint32.
//        A zero-variance (constant) template produces C = zeros(...).
//
//   A  : 2-D real numeric matrix (image to search).  Non-empty.
//        Same type constraints as T.
//
// RETURN VALUE
//   C  : double matrix, size (rows(A)+rows(T)-1) x (cols(A)+cols(T)-1).
//        All values clamped to [-1, 1].
//
// ERROR CONDITIONS
//   normxcorr2:InvalidNumInputs  - called with != 2 arguments.
//   normxcorr2:ComplexInput      - T or A contains complex values.
//   normxcorr2:NotMatrix         - T or A is not a 2-D matrix.
//   normxcorr2:EmptyInput        - T or A has zero elements.
//
// =============================================================================

funcprot(0)

// ---------------------------------------------------------------------------
// _fft2pad  -  Zero-padded 2-D FFT via two sequential 1-D passes.
//
// Pads X to size (P x Q) with zeros, then computes the 2-D DFT by
// applying a 1-D FFT across each row, then each column.
// This is the correct approach for Scilab 6.1.0, which has no built-in
// zero-padded 2-D FFT.
// ---------------------------------------------------------------------------
function F = _fft2pad(X, P, Q)
    [r, c] = size(X);
    Xpad   = zeros(P, Q);
    Xpad(1:r, 1:c) = double(X);
    // Row pass: 1-D FFT along each row
    Frow = zeros(P, Q);
    for i = 1:P
        Frow(i, :) = fft(Xpad(i, :));
    end
    // Column pass: 1-D FFT along each column
    F = zeros(P, Q);
    for j = 1:Q
        F(:, j) = fft(Frow(:, j));
    end
endfunction

// ---------------------------------------------------------------------------
// _ifft2  -  Inverse 2-D FFT via the conjugate identity.
//
// Identity: IFFT2(F) = conj( FFT2( conj(F) ) ) / (P*Q)
//
// This avoids ifft(), which does not exist in Scilab 6.1.0.
// ---------------------------------------------------------------------------
function X = _ifft2(F)
    [P, Q] = size(F);
    Fc   = conj(F);
    // Row pass on conjugate
    Frow = zeros(P, Q);
    for i = 1:P
        Frow(i, :) = fft(Fc(i, :));
    end
    // Column pass
    Ftmp = zeros(P, Q);
    for j = 1:Q
        Ftmp(:, j) = fft(Frow(:, j));
    end
    X = conj(Ftmp) / (P * Q);
endfunction

// ---------------------------------------------------------------------------
// Main function
// ---------------------------------------------------------------------------
function C = normxcorr2(T, A)

    // ------------------------------------------------------------------
    // 1. Argument-count guard
    // ------------------------------------------------------------------
    if argn(2) <> 2 then
        error("normxcorr2:InvalidNumInputs: " + ...
              "normxcorr2 requires exactly 2 inputs: normxcorr2(T, A).");
    end

    // ------------------------------------------------------------------
    // 2. Complex-input guard
    // ------------------------------------------------------------------
    if ~isreal(T) | ~isreal(A) then
        error("normxcorr2:ComplexInput: " + ...
              "T and A must be real numeric matrices.");
    end

    // ------------------------------------------------------------------
    // 3. Dimensionality guard
    //    Use ndims(), NOT size(x,"ndims") which is invalid in Scilab 6.1.0.
    // ------------------------------------------------------------------
    if ndims(T) <> 2 | ndims(A) <> 2 then
        error("normxcorr2:NotMatrix: " + ...
              "T and A must be 2-D matrices (grayscale images).");
    end

    // ------------------------------------------------------------------
    // 4. Empty-input guard
    // ------------------------------------------------------------------
    if size(T, "*") == 0 | size(A, "*") == 0 then
        error("normxcorr2:EmptyInput: " + ...
              "T and A must be non-empty.");
    end

    // ------------------------------------------------------------------
    // 5. Setup
    // ------------------------------------------------------------------
    Td = double(T);
    Ad = double(A);

    [m, n] = size(Td);
    [M, N] = size(Ad);

    outH = M + m - 1;
    outW = N + n - 1;

    // ------------------------------------------------------------------
    // 6. Zero-mean template and its L2-norm
    // ------------------------------------------------------------------
    T0   = Td - mean(Td);
    denT = sqrt(sum(T0 .^ 2));

    // Constant template => no meaningful correlation.
    if denT == 0 then
        C = zeros(outH, outW);
        return
    end

    // ------------------------------------------------------------------
    // 7. Numerator: full cross-correlation of A with zero-mean template T0.
    //    Computed as conv2(A, rot180(T0), "full") via FFT product.
    //    Invariant to additive offset in A because sum(T0) = 0.
    // ------------------------------------------------------------------
    rotT0 = T0($ : -1 : 1, $ : -1 : 1);
    FA    = _fft2pad(Ad,    outH, outW);
    FrotT = _fft2pad(rotT0, outH, outW);
    num   = real(_ifft2(FA .* FrotT));

    // ------------------------------------------------------------------
    // 8. Local window statistics of A
    //    Convolving with ones(m,n) computes the running sum over every
    //    m*n window simultaneously (integral image in frequency domain).
    // ------------------------------------------------------------------
    Fones = _fft2pad(ones(m, n), outH, outW);
    FA2   = _fft2pad(Ad .^ 2,    outH, outW);

    sumA  = real(_ifft2(FA  .* Fones));
    sumA2 = real(_ifft2(FA2 .* Fones));

    area = m * n;

    // Local variance: Var(W) = E[W^2] - (E[W])^2
    varA = sumA2 - (sumA .^ 2) / area;

    // ------------------------------------------------------------------
    // 9. Variance robustness: absolute threshold scaled by image magnitude.
    //
    //    FFT rounding noise in varA is of order:
    //        noise ~ eps_machine * area * max(|Ad|)^2
    //
    //    We zero out varA values below:
    //        varA_tol = eps_machine^(2/3) * max(1, max(|Ad|)^2)
    //
    //    This threshold is proportional to max(|Ad|)^2, the same scale as
    //    the noise, but is ~eps^(2/3)/eps = eps^(-1/3) ~ 1.6e5 times
    //    LARGER than the noise floor, giving a safe margin.
    //
    //    Crucially, it does NOT depend on sumA2 or DC offset, so it cannot
    //    accidentally suppress valid variances in high-DC images.
    // ------------------------------------------------------------------
    eps_machine = 2.22e-16;
    varA_tol    = (eps_machine ^ (2/3)) * max(1, max(abs(Ad(:))) ^ 2);

    varA = max(varA, 0);           // clamp rounding negatives to zero
    varA(varA < varA_tol) = 0;     // suppress rounding positives near zero

    stdA  = sqrt(varA);
    denom = denT .* stdA;

    // ------------------------------------------------------------------
    // 10. Normalised correlation
    //     Windows with denom == 0 (constant or fully zero-padded) get C = 0.
    // ------------------------------------------------------------------
    C    = zeros(outH, outW);
    mask = denom > 0;
    C(mask) = num(mask) ./ denom(mask);

    // Clamp to [-1, 1] to remove tiny floating-point overshoot.
    C = min(max(C, -1), 1);

endfunction

// =============================================================================
// END OF normxcorr2.sci
// =============================================================================
