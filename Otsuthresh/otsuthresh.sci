// =============================================================================
// otsuthresh.sci
// =============================================================================
//
// Copyright (C) 2026 - FOSSEE Summer Fellowship Submission
// This file is part of a Scilab Image Processing Toolbox implementation.
// Released under the GNU General Public License v3 or later.
// See <https://www.gnu.org/licenses/> for details.
//
// Function  : otsuthresh
// Purpose   : Compute Otsu's globally optimal binarisation threshold.
// Toolbox   : Image Processing (Scilab reimplementation of Octave 'image' pkg)
// Reference : Octave image package - otsuthresh()
//             https://octave.sourceforge.io/image/function/otsuthresh.html
//             N. Otsu, "A threshold selection method from gray-level
//             histograms," IEEE Trans. SMC, 9(1):62-66, 1979.
//
// Author    : FOSSEE Screening Task - Image Processing Toolbox Development
// Created   : April 2026
// Scilab    : 6.1.0+  (tested on 6.1.0 64-bit, Windows & Linux)
//
// =============================================================================
//
// CALLING SEQUENCE
//   t = otsuthresh(I)        // 2-D grayscale image (any real numeric type)
//   t = otsuthresh(counts)   // pre-computed histogram count vector
//
// DESCRIPTION
//   otsuthresh finds the threshold t in [0,1] that maximises the
//   between-class variance when pixels are split into two groups:
//   "background" (intensity <= t) and "foreground" (intensity > t).
//   This is Otsu's method (1979), which is equivalent to minimising the
//   weighted within-class variance.
//
//   Two calling modes are supported:
//
//   (a) Histogram mode  [vector input]:
//       Pass a 1-D non-negative vector of bin counts of length L.
//       The returned threshold is normalised so that bin 1 maps to 0 and
//       bin L maps to 1.  Typical usage: L = 256 for uint8 images.
//
//   (b) Image mode  [2-D matrix input]:
//       Pass a real 2-D numeric matrix (grayscale image).  Pixel values are
//       linearly scaled to [0, 1] (preserving relative differences), a
//       256-bin histogram is built internally, and mode (a) is applied.
//       Uniform images (maxI == minI) return t = 0 deterministically.
//
// ALGORITHM
//   For each candidate threshold k (bin index 1..L-1):
//
//     omega_0(k) = sum_{i=1}^{k} p(i)          // background class prob
//     omega_1(k) = 1 - omega_0(k)               // foreground class prob
//     mu_0(k)    = sum_{i=1}^{k} i*p(i) / omega_0(k)
//     mu_1(k)    = (mu_T - omega_0*mu_T_partial) / omega_1(k)
//
//   Between-class variance:
//     sigma_b^2(k) = omega_0 * omega_1 * (mu_0 - mu_1)^2
//                  = (mu_T * omega_0 - mu)^2 / (omega_0 * (1 - omega_0))
//
//   The threshold is the bin k that maximises sigma_b^2(k).
//   In case of ties the lowest index is chosen (first maximum).
//
// PARAMETERS
//   x  : real numeric.
//        - Vector (1-D): histogram bin counts.  Length >= 2.
//                        All values must be >= 0 and sum to > 0.
//        - Matrix (2-D): grayscale image.  Any real numeric type
//                        (double, uint8, uint16, int16, etc.).
//
// RETURN VALUE
//   t  : scalar double in [0, 1].
//        Normalised threshold.  Multiply by (L-1) to get a bin index,
//        or by 255 to get a uint8 pixel threshold.
//
// ERROR CONDITIONS
//   otsuthresh:InvalidNumInputs  - called with != 1 argument.
//   otsuthresh:ComplexInput      - input contains complex values.
//   otsuthresh:NotVectorOrMatrix - input is not a vector or 2-D matrix.
//   otsuthresh:TooFewBins        - histogram vector has fewer than 2 bins.
//   otsuthresh:NegativeCounts    - histogram contains negative values.
//   otsuthresh:ZeroHistogram     - histogram sums to zero.
//
// COMPATIBILITY NOTE
//   This file uses only Scilab 6.1.0 built-ins.
//   Specifically avoided:  size(x,"ndims")  [not valid in 6.1.0].
//   Use ndims(x) for dimensionality queries.
//
// =============================================================================

funcprot(0)

function t = otsuthresh(x)

    // ------------------------------------------------------------------
    // 1. Argument-count guard
    // ------------------------------------------------------------------
    if argn(2) <> 1 then
        error("otsuthresh:InvalidNumInputs: " + ...
              "otsuthresh requires exactly 1 input: otsuthresh(x).");
    end

    // ------------------------------------------------------------------
    // 2. Complex-input guard
    // ------------------------------------------------------------------
    if ~isreal(x) then
        error("otsuthresh:ComplexInput: " + ...
              "input must be a real numeric image or histogram counts.");
    end

    // ------------------------------------------------------------------
    // 3. Route: vector => histogram mode / 2-D matrix => image mode
    //    NOTE: use ndims(), NOT size(x,"ndims") which is invalid in 6.1.0
    // ------------------------------------------------------------------
    // Routing rule:
    //   - vector (including scalar 1x1) -> histogram mode.
    //     A 1-bin or 0-bin histogram is rejected by the TooFewBins guard below.
    //     This matches Octave: a scalar input is treated as a degenerate histogram.
    //   - 2-D non-vector matrix           -> image mode.
    if min(size(x)) == 1 then
        // ================================================================
        // PATH A: histogram mode
        // ================================================================
        counts = double(x(:));      // ensure column vector of doubles
        nbins  = length(counts);

        // --- validation ---
        if nbins < 2 then
            error("otsuthresh:TooFewBins: " + ...
                  "histogram vector must have at least 2 bins.");
        end
        if or(counts < 0) then
            error("otsuthresh:NegativeCounts: " + ...
                  "histogram counts must all be >= 0.");
        end
        total = sum(counts);
        if total <= 0 then
            error("otsuthresh:ZeroHistogram: " + ...
                  "histogram sum must be > 0 (all-zero histogram passed).");
        end

        // --- Otsu's algorithm (vectorised, no for-loop) ---
        p     = counts / total;                    // normalised probabilities
        omega = cumsum(p);                         // cumulative class probability
        // bin levels: 0-based index as intensity level
        levels = (0:nbins-1)';                     // column vector
        mu    = cumsum(p .* levels);               // cumulative mean
        mu_T  = mu(nbins);                         // total mean

        // Between-class variance numerator: (mu_T * omega - mu)^2
        num_b2 = (mu_T .* omega - mu) .^ 2;

        // Denominator: omega * (1 - omega)
        denom_b2 = omega .* (1 - omega);

        // Where denominator <= 0 (degenerate classes), set variance to -Inf
        // so those bins can never win the argmax.
        sigma_b2 = zeros(nbins, 1);
        valid = denom_b2 > 0;
        sigma_b2(valid)  = num_b2(valid) ./ denom_b2(valid);
        sigma_b2(~valid) = -%inf;

        // First maximum for tie-breaking stability
        [mx, k] = max(sigma_b2);

        // Normalise bin index to [0, 1]
        t = (k - 1) / (nbins - 1);

    elseif ndims(x) == 2 then
        // ================================================================
        // PATH B: image mode (2-D grayscale matrix)
        // ================================================================
        I    = double(x);
        minI = min(I);
        maxI = max(I);

        // Uniform image: threshold is arbitrary; return 0 for determinism.
        if maxI == minI then
            t = 0;
            return
        end

        // Linearly scale pixel values to [0, 1]
        J = (I - minI) / (maxI - minI);
        J = min(max(J, 0), 1);     // numerical clamp for floating-point edge cases

        // Build a 256-bin histogram via bin-index accumulation
        // Map [0,1] -> bin index in {1,...,256}
        nbins = 256;
        idx   = floor(J * (nbins - 1)) + 1;   // values in {1,...,256}
        idx   = min(idx, nbins);               // guard: prevent index 257 from round-off

        counts = zeros(nbins, 1);
        // Vectorised accumulation using loops avoided by accumarray-style trick:
        // Scilab 6.1.0 does not have accumarray; use a for-loop over unique bins
        // for correctness.  For typical 256-bin histograms this is negligible.
        flat = idx(:);
        for b = 1:nbins
            counts(b) = sum(flat == b);
        end

        // Recurse into histogram mode
        t = otsuthresh(counts);

    else
        error("otsuthresh:NotVectorOrMatrix: " + ...
              "input must be a vector (histogram) or a 2-D image matrix. " + ...
              msprintf("Got an array with %d dimensions.", ndims(x)));
    end

endfunction

// =============================================================================
// END OF otsuthresh.sci
// =============================================================================
