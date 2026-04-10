// =============================================================================
// immse_test.sce  –  Comprehensive Test Suite for immse()
// =============================================================================
//
// Copyright (C) 2026 – FOSSEE Summer Fellowship Submission
// Released under the GNU General Public License v3 or later.
//
// Description:
//   Standalone test script for the immse() function implemented in immse.sci.
//   Covers correctness, edge cases, integer overflow safety, and error paths.
//   Designed to run in Scilab 6.1.0+ without any external toolbox dependency.
//
// Usage:
//   exec("immse_test.sce", -1)
//
//   All tests print a structured line:
//       TEST: <name> | expected=<val> | actual=<val> | PASS=<T/F>
//   The suite aborts immediately on the first failure with a descriptive error.
//   On full success it prints:
//       RESULT: immse test suite PASSED  (N tests)
//
// Scilab: 6.1.0+ (64-bit, Windows & Linux)
// =============================================================================

funcprot(0)
exec("immse.sci", -1);

// =============================================================================
// SECTION 0 – Minimal Portable Test Harness
// =============================================================================
// Scilab 6.1.0 does not have a built-in unit-test framework; the helpers
// below provide consistent, human-readable output and fast-fail behaviour.

// _s(x) – safely coerce any value to a single string scalar for mprintf.
function s = _s(x)
    if type(x) == 10 then
        if size(x, "*") == 1 then
            s = x;
        else
            s = x(1);
        end
    else
        sv = string(x);
        if size(sv, "*") == 1 then
            s = sv;
        else
            s = sv(1);
        end
    end
endfunction

// _fmt(v) – format a scalar or array for display in test output.
function s = _fmt(v)
    if size(v, "*") == 1 then
        s = msprintf("%0.12g", double(v));
    else
        sz = size(v);
        parts = string(sz);
        s = msprintf("[%s %s]", strcat(parts, "x"), typeof(v));
    end
endfunction

// _printRow – emit one structured test-result line.
function _printRow(name, expectedStr, actualStr, statusStr)
    mprintf("TEST: %-42s | expected=%-22s | actual=%-22s | PASS=%s\n", ...
            _s(name), _s(expectedStr), _s(actualStr), _s(statusStr));
endfunction

// test_equal – exact equality check (integers, booleans, exact doubles).
function test_equal(name, actual, expected)
    ok = isequal(actual, expected);
    _printRow(name, _fmt(expected), _fmt(actual), _s(ok));
    if ~ok then
        error(msprintf("FAIL: %s  |  expected=%s  got=%s", ...
                       _s(name), _fmt(expected), _fmt(actual)));
    end
end

// test_almost – floating-point proximity check with configurable tolerance.
function test_almost(name, actual, expected, tol)
    if argn(2) < 4 then tol = 1e-12; end
    err = abs(double(actual) - double(expected));
    ok  = (err <= tol);
    mprintf("TEST: %-42s | expected=%-22s | actual=%-22s | err=%0.3g  tol=%g  PASS=%s\n", ...
            _s(name), _s(_fmt(expected)), _s(_fmt(actual)), err, tol, _s(ok));
    if ~ok then
        error(msprintf("FAIL: %s  (err=%g > tol=%g)", _s(name), err, tol));
    end
endfunction

// test_throws – verify that a Scilab command string raises an error.
function test_throws(name, cmd)
    ok = %f;
    try
        execstr(cmd);
    catch
        ok = %t;
    end
    _printRow(name, "throw", "throw", _s(ok));
    if ~ok then
        error(msprintf("FAIL: %s  (no error was thrown)", _s(name)));
    end
endfunction

// Running test counter
global _N_PASS;
_N_PASS = 0;

// Wrapper so we can count tests automatically without instrumenting every call.
// (Not needed for simple scripts – we'll just tally at the end via a final count.)

// =============================================================================
// SECTION 1 – Basic Correctness Tests
// =============================================================================
mprintf("\n--- SECTION 1: Basic Correctness ---\n");

// TC-01  Identical arrays → MSE must be exactly 0.
//        This tests the zero-path and verifies no floating-point drift when
//        the inputs are the same object.
A = [1 2; 3 4];
test_equal("TC-01 identical 2x2 double  → 0", immse(A, A), 0);

// TC-02  Known small matrix with hand-computable MSE.
//        A=[1 2;3 4], B=[2 2;2 2]
//        D=[-1 0;1 2], D²=[1 0;1 4], mean=6/4=1.5
A = [1 2; 3 4];
B = [2 2; 2 2];
test_almost("TC-02 known 2x2 double  → 1.5", immse(A, B), 1.5, 1e-12);

// TC-03  All-zeros vs all-ones (flat images) – MSE should be 1.0.
test_almost("TC-03 zeros vs ones (2x2)  → 1.0", ...
            immse(zeros(2,2), ones(2,2)), 1.0, 1e-12);

// TC-04  Single-element scalar inputs.
test_almost("TC-04 scalar  immse(3, 7) → 16", immse(3, 7), 16.0, 1e-12);

// TC-05  Row-vector inputs.
//        A=[0 1 2 3], B=[1 1 1 1], D=[-1 0 1 2], D²=[1 0 1 4], mean=1.5
test_almost("TC-05 1x4 row-vector  → 1.5", ...
            immse([0 1 2 3], [1 1 1 1]), 1.5, 1e-12);

// TC-06  Column-vector inputs – same arithmetic, different shape.
test_almost("TC-06 4x1 column-vector  → 1.5", ...
            immse([0;1;2;3], [1;1;1;1]), 1.5, 1e-12);

// TC-07  Large uniform offset: A = zeros(3,3), B = 5*ones(3,3), MSE = 25.
test_almost("TC-07 uniform offset 5  → 25", ...
            immse(zeros(3,3), 5*ones(3,3)), 25.0, 1e-12);

// TC-08  Negative values – immse must handle signed doubles.
//        A=[-2 -1], B=[1 1], D=[-3 -2], D²=[9 4], mean=6.5
test_almost("TC-08 negative values  → 6.5", ...
            immse([-2 -1], [1 1]), 6.5, 1e-12);

// =============================================================================
// SECTION 2 – Integer-Type Safety (Overflow Prevention)
// =============================================================================
mprintf("\n--- SECTION 2: Integer-Type Safety ---\n");

// TC-09  uint8: subtraction WITHOUT double cast would overflow on many values.
//        A=[0 10;20 30], B=[0 0;0 0], D=[0 10;20 30], D²=[0 100;400 900]
//        mean = 1400/4 = 350
A8 = uint8([0 10; 20 30]);
B8 = uint8([0  0;  0  0]);
test_almost("TC-09 uint8 no-overflow  → 350", immse(A8, B8), 350, 1e-12);

// TC-10  uint8: values near the 255 boundary.
//        A=[200 250], B=[10 5], D=[190 245], D²=[36100 60025], mean=48062.5
A8b = uint8([200 250]);
B8b = uint8([10    5]);
test_almost("TC-10 uint8 near-max  → 48062.5", immse(A8b, B8b), 48062.5, 1e-12);

// TC-11  uint16: typical 16-bit medical image range.
//        A=[1000 2000], B=[0 0], D=[1000 2000], D²=[1e6 4e6], mean=2.5e6
A16 = uint16([1000 2000]);
B16 = uint16([0 0]);
test_almost("TC-11 uint16 range  → 2500000", immse(A16, B16), 2500000, 1e-12);

// TC-12  int16: signed integer type, A positive, B negative.
//        A=[100], B=[-100], D=[200], D²=[40000], mean=40000
test_almost("TC-12 int16 signed  → 40000", ...
            immse(int16([100]), int16([-100])), 40000, 1e-12);

// =============================================================================
// SECTION 3 – Multi-Dimensional Array Support
// =============================================================================
mprintf("\n--- SECTION 3: Multi-Dimensional Arrays ---\n");

// TC-13  3D array: ones(2,2,3) vs zeros(2,2,3).
//        All 12 differences = 1, so MSE = 1.
test_almost("TC-13 3D ones vs zeros  → 1.0", ...
            immse(ones(2,2,3), zeros(2,2,3)), 1.0, 1e-12);

// TC-14  3D RGB-like array with known per-pixel error.
//        Channel 1: diff=0, Channel 2: diff=1 everywhere, Channel 3: diff=2.
//        All channels 2x2 → 12 elements total.
//        D²: 4 zeros + 4 ones + 4 fours = 20, mean = 20/12 = 5/3
A3 = cat(3, zeros(2,2), ones(2,2), zeros(2,2));
B3 = cat(3, zeros(2,2), zeros(2,2), 2*ones(2,2));
test_almost("TC-14 3D mixed-channel  → 5/3", ...
            immse(A3, B3), 5/3, 1e-12);

// TC-15  Non-square 2D image (4x5).
//        All differences = 3, so MSE = 9.
test_almost("TC-15 4x5 non-square, diff=3  → 9", ...
            immse(3*ones(4,5), zeros(4,5)), 9.0, 1e-12);

// =============================================================================
// SECTION 4 – Symmetry and Commutativity
// =============================================================================
mprintf("\n--- SECTION 4: Symmetry / Commutativity ---\n");

// TC-16  immse(A,B) == immse(B,A)  (MSE is symmetric).
A = [1 3 5; 7 9 11];
B = [2 2 6; 6 10 10];
test_almost("TC-16 symmetry immse(A,B)==immse(B,A)", ...
            immse(A, B), immse(B, A), 0);

// TC-17  Constant offset commutes: immse(A, A+k) == immse(A+k, A).
A = rand(5, 5);
k = 0.3;
test_almost("TC-17 constant-offset symmetry (5x5 rand)", ...
            immse(A, A+k), immse(A+k, A), 0);

// =============================================================================
// SECTION 5 – Numerical Precision
// =============================================================================
mprintf("\n--- SECTION 5: Numerical Precision ---\n");

// TC-18  Very small differences (near floating-point epsilon).
//        D = eps for each element, MSE = eps²  (≈ 0 to machine precision).
e = 1e-14;
test_almost("TC-18 tiny diff (~eps)  → ~0", ...
            immse(zeros(1,4), e*ones(1,4)), e^2, 1e-30);

// TC-19  Very large values – check for no spurious overflow in double.
//        A=[1e15], B=[0], MSE = 1e30
test_almost("TC-19 very large values  → 1e30", ...
            immse([1e15], [0]), 1e30, 1e18);

// TC-20  Fractional / non-integer inputs.
//        A=[0.1 0.4], B=[0.3 0.1], D=[-0.2 0.3], D²=[0.04 0.09], mean=0.065
test_almost("TC-20 fractional inputs  → 0.065", ...
            immse([0.1 0.4], [0.3 0.1]), 0.065, 1e-15);

// =============================================================================
// SECTION 6 – Error-Path Tests  (all must THROW)
// =============================================================================
mprintf("\n--- SECTION 6: Error Paths ---\n");

// TC-21  Too few arguments.
test_throws("TC-21 too few args  (1 arg)  → throw", "immse([1 2 3]);");

// TC-22  Too many arguments.
test_throws("TC-22 too many args (3 args) → throw", "immse([1],[1],[1]);");

// TC-23  Size mismatch: 1x2 vs 1x3.
test_throws("TC-23 size mismatch 1x2 vs 1x3 → throw", "immse([1 2],[1 2 3]);");

// TC-24  Size mismatch: 2x2 vs 4x1 (same numel, different shape).
test_throws("TC-24 shape mismatch 2x2 vs 4x1 → throw", "immse(ones(2,2),ones(4,1));");

// TC-25  Empty arrays.
test_throws("TC-25 empty arrays [] → throw", "immse([],[]);");

// TC-26  Complex input rejected.
test_throws("TC-26 complex A → throw", "immse([1+2*%i 3],[1 3]);");

// TC-27  Complex B rejected.
test_throws("TC-27 complex B → throw", "immse([1 2],[1+%i 2]);");

// =============================================================================
// SUMMARY
// =============================================================================
mprintf("\n----------------------------------------------------------------\n");
mprintf("RESULT: immse test suite PASSED  (27 tests)\n");
mprintf("----------------------------------------------------------------\n");

// =============================================================================
// END OF immse_test.sce
// =============================================================================
