// =============================================================================
// immse.sci
// =============================================================================
//
// Copyright (C) 2026 - FOSSEE Summer Fellowship Submission
// This file is part of a Scilab Image Processing Toolbox implementation.
// Released under the GNU General Public License v3 or later.
// See <https://www.gnu.org/licenses/> for details.
//
// Function  : immse
// Purpose   : Compute the Mean Squared Error (MSE) between two images or
//             numeric arrays of identical dimensions.
// Toolbox   : Image Processing (Scilab reimplementation of Octave 'image' pkg)
// Reference : Octave image package – immse()
//             https://octave.sourceforge.io/image/function/immse.html
//
// Author    : FOSSEE Screening Task – Image Processing Toolbox Development
// Created   : April 2026
// Scilab    : 6.1.0+ (tested on 6.1.0 64-bit, Windows & Linux)
//
// =============================================================================

// ---------------------------------------------------------------------------
// FUNCTION SIGNATURE
// ---------------------------------------------------------------------------
//
//   mse = immse(A, B)
//
// ---------------------------------------------------------------------------
// DESCRIPTION
// ---------------------------------------------------------------------------
//
//   immse computes the Mean Squared Error between two numeric arrays A and B.
//   This is a fundamental metric in image quality assessment, compression
//   analysis, and denoising evaluation.
//
//   The MSE is defined as:
//
//       MSE = (1/N) * sum( (A_i - B_i)^2 )   for i = 1 … N
//
//   where N is the total number of elements.
//
//   Both inputs are cast to double before subtraction to prevent integer
//   overflow (e.g., uint8(0) - uint8(200) would wrap incorrectly otherwise).
//
//   Supports 2D grayscale images, 3D colour images (H×W×C), and arbitrary
//   N-dimensional numeric arrays, matching Octave's immse behaviour.
//
// ---------------------------------------------------------------------------
// PARAMETERS
// ---------------------------------------------------------------------------
//
//   A  –  First input array / image.
//         Type : real numeric (double, int8, int16, int32, uint8, uint16,
//                uint32, single). Complex inputs are rejected.
//         Size : Any – but must match B exactly.
//
//   B  –  Second input array / image.
//         Type : same constraints as A.
//         Size : Must be identical to size(A).
//
// ---------------------------------------------------------------------------
// RETURN VALUE
// ---------------------------------------------------------------------------
//
//   mse  –  Non-negative scalar (double).
//            0  when A and B are identical.
//           >0  otherwise.
//
// ---------------------------------------------------------------------------
// ERROR CONDITIONS
// ---------------------------------------------------------------------------
//
//   immse:InvalidNumInputs      – called with ≠ 2 arguments.
//   immse:EmptyInput            – either A or B has zero elements.
//   immse:ComplexNotSupported   – either A or B is complex-valued.
//   immse:SizeMismatch          – A and B have different sizes.
//
// ---------------------------------------------------------------------------

funcprot(0)  // suppress "redefining function" warnings on repeated exec

function mse = immse(A, B)

    // ------------------------------------------------------------------
    // 1. Argument-count guard
    // ------------------------------------------------------------------
    if argn(2) <> 2 then
        error("immse:InvalidNumInputs: immse requires exactly 2 inputs: immse(A, B).");
    end

    // ------------------------------------------------------------------
    // 2. Empty-input guard
    // ------------------------------------------------------------------
    if size(A, "*") == 0 | size(B, "*") == 0 then
        error("immse:EmptyInput: A and B must be non-empty arrays.");
    end

    // ------------------------------------------------------------------
    // 3. Real-valued guard  (complex images not supported)
    // ------------------------------------------------------------------
    if ~isreal(A) | ~isreal(B) then
        error("immse:ComplexNotSupported: " + ...
              "A and B must be real numeric arrays. " + ...
              "Complex images are not supported.");
    end

    // ------------------------------------------------------------------
    // 4. Size-compatibility guard
    // ------------------------------------------------------------------
    if or(size(A) <> size(B)) then
        error("immse:SizeMismatch: " + ...
              "A and B must have the same size. " + ...
              msprintf("Got size(A)=[%s] vs size(B)=[%s].", ...
                       strcat(string(size(A)), "x"), ...
                       strcat(string(size(B)), "x")));
    end

    // ------------------------------------------------------------------
    // 5. Compute MSE
    //    Cast to double FIRST to avoid integer-arithmetic overflow when
    //    dealing with uint8/uint16/int* typed images.
    // ------------------------------------------------------------------
    Ad = double(A);
    Bd = double(B);

    // Flatten to a column vector so mean() always returns a scalar,
    // regardless of whether the input is 2D, 3D, or higher-dimensional.
    D = Ad(:) - Bd(:);

    mse = mean(D .^ 2);

endfunction

// =============================================================================
// END OF immse.sci
// =============================================================================
