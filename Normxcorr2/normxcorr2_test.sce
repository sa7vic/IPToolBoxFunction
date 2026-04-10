// =============================================================================
// normxcorr2_test.sce  -  Comprehensive Test Suite for normxcorr2()
// =============================================================================
//
// Copyright (C) 2026 - FOSSEE Summer Fellowship Submission
// Released under the GNU General Public License v3 or later.
//
// Description:
//   Standalone test script for normxcorr2() implemented in normxcorr2.sci.
//   Covers output range, output size, degenerate inputs, known exact values,
//   symmetry properties, and all error paths.
//   No external toolbox required.  Runs on Scilab 6.1.0+.
//
// Usage:
//   exec("normxcorr2_test.sce", -1)
//
// Scilab: 6.1.0+
// =============================================================================

funcprot(0)
exec("normxcorr2.sci", -1);

// =============================================================================
// SECTION 0 - Test Harness
// =============================================================================

function s = _s(x)
    if type(x) == 10 then
        if size(x, "*") == 1 then; s = x; else; s = x(1); end
    else
        sv = string(x);
        if size(sv, "*") == 1 then; s = sv; else; s = sv(1); end
    end
endfunction

function s = _fmt(v)
    if size(v, "*") == 1 then
        s = msprintf("%0.10g", double(v));
    else
        sz = size(v);
        s = msprintf("[%dx%d]", sz(1), sz(2));
    end
endfunction

function _printRow(name, expectedStr, actualStr, statusStr)
    mprintf("TEST: %-52s | expected=%-16s | actual=%-16s | PASS=%s\n", ...
            _s(name), _s(expectedStr), _s(actualStr), _s(statusStr));
endfunction

function test_almost(name, actual, expected, tol)
    if argn(2) < 4 then tol = 1e-8; end
    err = abs(double(actual) - double(expected));
    ok  = (err <= tol);
    mprintf("TEST: %-52s | expected=%-16s | actual=%-16s | err=%0.3g tol=%g PASS=%s\n", ...
            _s(name), _s(_fmt(expected)), _s(_fmt(actual)), err, tol, _s(ok));
    if ~ok then
        error(msprintf("FAIL: %s  (err=%g > tol=%g)", _s(name), err, tol));
    end
endfunction

function test_in_range(name, actual, lo, hi)
    ok = (actual >= lo) & (actual <= hi);
    _printRow(name, msprintf("[%g,%g]", lo, hi), _fmt(actual), _s(ok));
    if ~ok then
        error(msprintf("FAIL: %s  value=%g not in [%g, %g]", _s(name), actual, lo, hi));
    end
endfunction

function test_throws(name, cmd)
    ok = %f;
    try; execstr(cmd); catch; ok = %t; end
    _printRow(name, "throw", "throw", _s(ok));
    if ~ok then
        error(msprintf("FAIL: %s  (no error was thrown)", _s(name)));
    end
endfunction

// ---------------------------------------------------------------------------
// Brute-force normxcorr2 (pure Scilab, no FFT) for correctness comparison.
// Computes the "full" normalised cross-correlation position by position.
// Used as ground truth against the FFT implementation.
// ---------------------------------------------------------------------------
function C = _normxcorr2_brute(T, A)
    Td = double(T);
    Ad = double(A);
    [m, n] = size(Td);
    [M, N] = size(Ad);
    outH = M + m - 1;
    outW = N + n - 1;
    meanT = mean(Td);
    T0    = Td - meanT;
    denT  = sqrt(sum(T0 .^ 2));
    C = zeros(outH, outW);
    if denT == 0 then; return; end
    for i = 1:outH
        for j = 1:outW
            // Build zero-padded m*n window of A centred at (i,j)
            Apad = zeros(m, n);
            for ti = 1:m
                for tj = 1:n
                    ai = i - m + ti;    // row index in A
                    aj = j - n + tj;    // col index in A
                    if ai >= 1 & ai <= M & aj >= 1 & aj <= N then
                        Apad(ti, tj) = Ad(ai, aj);
                    end
                end
            end
            meanA = sum(Apad) / (m * n);
            A0    = Apad - meanA;
            denA  = sqrt(sum(A0 .^ 2));
            if denA > 0 then
                C(i, j) = sum(A0 .* T0) / (denA * denT);
            end
        end
    end
    C = min(max(C, -1), 1);
endfunction

// =============================================================================
// SECTION 1 - Output Range: all values in [-1, 1]
// =============================================================================
mprintf("\n--- SECTION 1: Output Range ---\n");

// TC-01  Random double inputs.
A = rand(10, 12);
T = rand(3, 4);
C = normxcorr2(T, A);
test_in_range("TC-01 rand inputs: min(C) in [-1,1]", min(C), -1, 1);
test_in_range("TC-02 rand inputs: max(C) in [-1,1]", max(C), -1, 1);

// TC-03  uint8 inputs.
A8 = uint8(round(rand(8,8) * 255));
T8 = uint8(round(rand(3,3) * 255));
C8 = normxcorr2(T8, A8);
test_in_range("TC-03 uint8 inputs: all in [-1,1]",   max(abs(C8)), 0, 1);

// =============================================================================
// SECTION 2 - Output Size
// =============================================================================
mprintf("\n--- SECTION 2: Output Size ---\n");

// TC-04  Size: rows(C) = rows(A)+rows(T)-1, cols(C) = cols(A)+cols(T)-1.
A = rand(6, 7);
T = rand(2, 3);
C = normxcorr2(T, A);
expH = 6 + 2 - 1;
expW = 7 + 3 - 1;
ok_size = isequal(size(C), [expH, expW]);
_printRow("TC-04 output size (6,7)+(2,3)-1 = (7,9)", ...
          msprintf("[%d,%d]", expH, expW), ...
          msprintf("[%d,%d]", size(C,1), size(C,2)), _s(ok_size));
if ~ok_size then
    error("FAIL: TC-04 wrong output size");
end

// TC-05  Template same size as image => output is (2M-1) x (2N-1).
A5 = rand(4, 5);
T5 = rand(4, 5);
C5 = normxcorr2(T5, A5);
ok5 = isequal(size(C5), [7, 9]);
_printRow("TC-05 T same size as A => (7,9)", "[7,9]", ...
          msprintf("[%d,%d]", size(C5,1), size(C5,2)), _s(ok5));
if ~ok5 then; error("FAIL: TC-05 wrong output size"); end

// TC-06  Scalar template => output same size as A.
C6 = normxcorr2([2], rand(5,6));
ok6 = isequal(size(C6), [5, 6]);
_printRow("TC-06 scalar T => output size = size(A)", "[5,6]", ...
          msprintf("[%d,%d]", size(C6,1), size(C6,2)), _s(ok6));
if ~ok6 then; error("FAIL: TC-06 wrong output size for scalar T"); end

// =============================================================================
// SECTION 3 - Degenerate Inputs
// =============================================================================
mprintf("\n--- SECTION 3: Degenerate Inputs ---\n");

// TC-07  Constant (zero-variance) template => C must be all zeros.
T_const = ones(3, 3) * 5;
A_rand  = rand(8, 8);
C7 = normxcorr2(T_const, A_rand);
test_almost("TC-07 const template => all zeros", max(abs(C7)), 0, 1e-10);

// TC-08  Constant image: interior windows (fully inside A) have zero variance => C=0.
//        Border positions in the full-correlation output include zero-padded partial
//        windows, which ARE non-constant (mix of 7s and 0s), so C may be nonzero
//        there.  We test only the interior sub-matrix where the full template fits.
//        Interior rows: m..M = 3..8, interior cols: n..N = 3..8  (1-based, in output).
T_rand  = rand(3, 3);
A_const = ones(8, 8) * 7;
[m8, n8] = size(T_rand);   // template size = (3,3)
[M8, N8] = size(A_const);  // image size   = (8,8)
C8b      = normxcorr2(T_rand, A_const);
C8b_interior = C8b(m8:M8, n8:N8);   // sub-matrix where full template fits inside A
test_almost("TC-08 const image interior => all zeros", max(abs(C8b_interior)), 0, 1e-6);

// TC-09  Scalar image (1x1 A): output is 1x1; value is 0 (window = single
//        pixel, mean = pixel, so zero-mean window => zero std => C=0).
C9 = normxcorr2([1], [5]);
test_almost("TC-09 scalar A, scalar T => 0",      C9(1,1), 0, 1e-10);

// =============================================================================
// SECTION 4 - Known Exact Values (FFT vs Brute Force)
// =============================================================================
mprintf("\n--- SECTION 4: FFT vs Brute-Force Ground Truth ---\n");

// TC-10  Small integer matrix: verify every output element matches.
A10 = [1 2 3 4 5;
       5 4 3 2 1;
       1 2 1 2 1;
       0 1 0 1 0;
       2 2 2 2 2];
T10 = [1 2 1;
       0 1 0;
       1 2 1];
C_fft   = normxcorr2(T10, A10);
C_brute = _normxcorr2_brute(T10, A10);
max_diff10 = max(abs(C_fft - C_brute));
test_almost("TC-10 5x5 A, 3x3 T: FFT==brute", max_diff10, 0, 1e-8);

// TC-11  Randomly generated 7x8 image with 2x3 template.
rand("seed", 42);    
A11 = rand(7, 8);
T11 = rand(2, 3);
C11_fft   = normxcorr2(T11, A11);
C11_brute = _normxcorr2_brute(T11, A11);
max_diff11 = max(abs(C11_fft - C11_brute));
test_almost("TC-11 rand 7x8 A, 2x3 T: FFT==brute", max_diff11, 0, 1e-7);

// TC-12  Template is a subimage of A => perfect match (C=1) expected at
//        the corresponding position.
A12 = [0 0 0 0 0;
       0 1 2 0 0;
       0 3 4 0 0;
       0 0 0 0 0];
T12 = [1 2; 3 4];
C12 = normxcorr2(T12, A12);
// The perfect-match position in full correlation is at (rows(T)+row_in_A-1,
// cols(T)+col_in_A-1) = (2+2-1, 2+2-1) = (3,3).
test_almost("TC-12 template = subimage of A => C(3,3)=1", C12(3,3), 1, 1e-8);

// TC-13  Negated template => negated image patch => C = -1 at match position.
T13 = -T12;
C13 = normxcorr2(T13, A12);
test_almost("TC-13 negated template => C(3,3)=-1", C13(3,3), -1, 1e-8);

// =============================================================================
// SECTION 5 - Symmetry and Properties
// =============================================================================
mprintf("\n--- SECTION 5: Symmetry and Properties ---\n");

// TC-14  Scaling A by a positive constant does not change C.
//        normxcorr2(T, k*A) == normxcorr2(T, A)  for k > 0.
A14 = rand(6, 6);
T14 = rand(3, 3);
C14a = normxcorr2(T14, A14);
C14b = normxcorr2(T14, 5 * A14);
test_almost("TC-14 scale A: C unchanged", max(abs(C14a - C14b)), 0, 1e-8);

// TC-15  Adding a constant to A does not change C.
C15 = normxcorr2(T14, A14 + 100);
[m14, n14] = size(T14);
[M14, N14] = size(A14);
C14a_interior = C14a(m14:M14, n14:N14);
C15_interior  = C15(m14:M14, n14:N14);
test_almost("TC-15 offset A: C unchanged", max(abs(C14a_interior - C15_interior)), 0, 1e-8);

// TC-16  Scaling T by a positive constant does not change C.
C16 = normxcorr2(7 * T14, A14);
test_almost("TC-16 scale T: C unchanged", max(abs(C14a - C16)), 0, 1e-8);

// =============================================================================
// SECTION 6 - Error Paths (all must THROW)
// =============================================================================
mprintf("\n--- SECTION 6: Error Paths ---\n");

// TC-17  Too few arguments.
test_throws("TC-17 1 arg          → throw", "normxcorr2(rand(3,3));");

// TC-18  Too many arguments.
test_throws("TC-18 3 args         → throw", "normxcorr2(rand(2,2),rand(4,4),1);");

// TC-19  Complex T.
test_throws("TC-19 complex T      → throw", "normxcorr2([1+%i 2;3 4], rand(5,5));");

// TC-20  Complex A.
test_throws("TC-20 complex A      → throw", "normxcorr2(rand(2,2), [1+%i 2;3 4]);");

// TC-21  3-D T.
test_throws("TC-21 3-D T          → throw", "normxcorr2(ones(2,2,2), rand(5,5));");

// TC-22  3-D A.
test_throws("TC-22 3-D A          → throw", "normxcorr2(rand(2,2), ones(5,5,3));");

// TC-23  Empty T.
test_throws("TC-23 empty T        → throw", "normxcorr2([], rand(5,5));");

// TC-24  Empty A.
test_throws("TC-24 empty A        → throw", "normxcorr2(rand(2,2), []);");

// =============================================================================
// SUMMARY
// =============================================================================
mprintf("\n");
mprintf("----------------------------------------------------------------\n");
mprintf("RESULT: normxcorr2 test suite PASSED  (24 tests)\n");
mprintf("----------------------------------------------------------------\n");

// =============================================================================
// END OF normxcorr2_test.sce
// =============================================================================
