// =============================================================================
// otsuthresh_test.sce  -  Comprehensive Test Suite for otsuthresh()
// =============================================================================
//
// Copyright (C) 2026 - FOSSEE Summer Fellowship Submission
// Released under the GNU General Public License v3 or later.
//
// Description:
//   Standalone test script for otsuthresh() implemented in otsuthresh.sci.
//   Covers histogram mode, image mode, integer types, edge cases, and all
//   error paths.  No external toolbox dependency; runs on Scilab 6.1.0+.
//
// Usage:
//   exec("otsuthresh_test.sce", -1)
//
//   Each test prints:
//       TEST: <name> | expected=<val> | actual=<val> | PASS=T/%F
//   The suite aborts on the first failure with a descriptive message.
//   On full success:
//       RESULT: otsuthresh test suite PASSED  (N tests)
//
// Scilab: 6.1.0+
// =============================================================================

funcprot(0)
exec("otsuthresh.sci", -1);

// =============================================================================
// SECTION 0 - Test Harness  (same portable harness as immse_test)
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
        s = msprintf("%0.12g", double(v));
    else
        sz = size(v);
        s = msprintf("[%dx%d %s]", sz(1), sz(2), typeof(v));
    end
endfunction

function _printRow(name, expectedStr, actualStr, statusStr)
    mprintf("TEST: %-48s | expected=%-18s | actual=%-18s | PASS=%s\n", ...
            _s(name), _s(expectedStr), _s(actualStr), _s(statusStr));
endfunction

function test_almost(name, actual, expected, tol)
    if argn(2) < 4 then tol = 1e-12; end
    err = abs(double(actual) - double(expected));
    ok  = (err <= tol);
    mprintf("TEST: %-48s | expected=%-18s | actual=%-18s | err=%0.3g tol=%g PASS=%s\n", ...
            _s(name), _s(_fmt(expected)), _s(_fmt(actual)), err, tol, _s(ok));
    if ~ok then
        error(msprintf("FAIL: %s  (err=%g > tol=%g)", _s(name), err, tol));
    end
endfunction

function test_in_range(name, actual, lo, hi)
    ok = (actual >= lo) & (actual <= hi);
    _printRow(name, msprintf("[%g,%g]", lo, hi), _fmt(actual), _s(ok));
    if ~ok then
        error(msprintf("FAIL: %s  value=%g not in [%g,%g]", _s(name), actual, lo, hi));
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

// =============================================================================
// SECTION 1 - Histogram Mode: Return Value in [0, 1]
// =============================================================================
mprintf("\n--- SECTION 1: Histogram Mode - Range ---\n");

// TC-01  Classic bimodal histogram: half dark (bin 1), half bright (bin 256).
//        Threshold should land near the midpoint => ~0.5
counts = zeros(256, 1);
counts(1)   = 1000;
counts(256) = 1000;
t = otsuthresh(counts);
test_in_range("TC-01 bimodal hist (bins 1 & 256)  in [0,1]", t, 0, 1);

// TC-02  Single bright spike at bin 128 out of 256
//        (degenerate: only one class, threshold should still be in [0,1])
counts2 = zeros(256, 1);
counts2(128) = 500;
t2 = otsuthresh(counts2);
test_in_range("TC-02 single-spike hist            in [0,1]", t2, 0, 1);

// TC-03  Uniform histogram (all bins equal) - all between-class variances
//        are equal; first maximum is bin 1, so t = 0.
counts3 = ones(256, 1) * 10;
t3 = otsuthresh(counts3);
test_in_range("TC-03 uniform hist                 in [0,1]", t3, 0, 1);

// TC-04  Row-vector input (should work identically to column vector).
counts4 = zeros(1, 64);
counts4(1)  = 200;
counts4(64) = 200;
t4 = otsuthresh(counts4);
test_in_range("TC-04 row-vector hist (64 bins)    in [0,1]", t4, 0, 1);

// =============================================================================
// SECTION 2 - Histogram Mode: Known Exact Thresholds
// =============================================================================
mprintf("\n--- SECTION 2: Histogram Mode - Known Values ---\n");

// TC-05  Two-bin histogram: equal counts in bin 1 and bin 2.
//        Only one candidate threshold (k=1): omega_0=0.5, sigma_b2=0.25.
//        k=1 maps to t = (1-1)/(2-1) = 0.
counts5 = [100; 100];
test_almost("TC-05 2-bin equal counts       → 0", otsuthresh(counts5), 0, 0);

// TC-06  Three-bin histogram: all mass in bin 3.
//        Both k=1 and k=2 have omega_0=0, denom=0 => sigma_b2=-Inf.
//        First max is k=1 => t = 0.
counts6 = [0; 0; 500];
test_almost("TC-06 3-bin single-bright-bin  → 0", otsuthresh(counts6), 0, 0);

// TC-07  Three-bin histogram with known separation.
//        counts = [500, 0, 500], nbins=3.
//        k=1: omega=0.5, mu=0, mu_T=1, num=(1*0.5-0)^2=0.25, den=0.25 => 1.0
//        k=2: omega=0.5, mu=0, same => 1.0  (tie: first max wins => k=1)
//        t = (1-1)/(3-1) = 0.
counts7 = [500; 0; 500];
test_almost("TC-07 3-bin bimodal            → 0", otsuthresh(counts7), 0, 0);

// TC-08  Bimodal with known separation: 4-bin [100 0 0 100].
//        Candidates k=1,2,3.
//        k=1: omega=0.5, mu=0, mu_T=1.5, num=(1.5*0.5-0)^2=0.5625, den=0.25 => 2.25
//        k=2: omega=0.5, mu=0, same result => 2.25  (tie => k=1 wins)
//        k=3: omega=0.5, mu=0, same => 2.25          (tie => k=1 wins)
//        t = (1-1)/(4-1) = 0
counts8 = [100; 0; 0; 100];
test_almost("TC-08 4-bin perfect bimodal    → 0", otsuthresh(counts8), 0, 0);

// TC-09  Verify normalisation: result must be (k-1)/(nbins-1).
//        Use a 5-bin histogram heavily loaded at bin 4.
//        counts=[10 10 10 970 0], nbins=5.
//        Dominant peak at bin 4 (0-based: 3) => optimal k in {3 or 4}.
//        Either way t must be in [0,1] and a valid (k-1)/4.
counts9 = [10; 10; 10; 970; 0];
t9 = otsuthresh(counts9);
test_in_range("TC-09 normalisation check (5-bin)  in [0,1]", t9, 0, 1);
// Verify t is one of {0, 0.25, 0.5, 0.75, 1}
valid_vals = [0, 0.25, 0.5, 0.75, 1];
found = %f;
for v = valid_vals
    if abs(t9 - v) < 1e-12 then; found = %t; end
end
_printRow("TC-09b normalisation is (k-1)/4", "one of {0,0.25,...,1}", _fmt(t9), _s(found));
if ~found then
    error(msprintf("FAIL TC-09b: t=%g not a valid (k-1)/4 value", t9));
end

// =============================================================================
// SECTION 3 - Image Mode
// =============================================================================
mprintf("\n--- SECTION 3: Image Mode ---\n");

// TC-10  Bimodal double image: clear two-class structure.
//        Left half = 0.0, right half = 1.0 => threshold near 0.5.
I10 = [zeros(10,10), ones(10,10)];
t10 = otsuthresh(I10);
test_in_range("TC-10 double bimodal image       in [0,1]", t10, 0, 1);

// TC-11  Uniform image (all pixels equal): must return exactly 0.
Iu = ones(8,8) * 7;
test_almost("TC-11 uniform image (all=7)      → 0", otsuthresh(Iu), 0, 0);

// TC-12  Random double image in [0,1]: result must be in (0,1).
Ir = rand(50, 50);
t12 = otsuthresh(Ir);
test_in_range("TC-12 rand(50,50) double         in [0,1]", t12, 0, 1);

// TC-13  uint8 image with bimodal distribution.
I13 = uint8([zeros(1,128), 255*ones(1,128)]);
t13 = otsuthresh(I13);
test_in_range("TC-13 uint8 bimodal 1x256        in [0,1]", t13, 0, 1);

// TC-14  uint16 image: values in [0, 65535] range.
I14 = uint16([zeros(1,50), 60000*ones(1,50)]);
t14 = otsuthresh(I14);
test_in_range("TC-14 uint16 bimodal 1x100       in [0,1]", t14, 0, 1);

// TC-15  int16 image: negative pixel values (valid signed range).
I15 = int16([-100*ones(5,5), 100*ones(5,5)]);
t15 = otsuthresh(I15);
test_in_range("TC-15 int16 signed bimodal 5x10  in [0,1]", t15, 0, 1);

// TC-16  Non-square image.
I16 = rand(7, 13);
t16 = otsuthresh(I16);
test_in_range("TC-16 rand(7,13) non-square      in [0,1]", t16, 0, 1);

// TC-17  Small uniform image (all pixels identical) => return 0.
//        Uses a 2x2 matrix to exercise the image-mode path.
//        (A scalar/1-element input is treated as a degenerate histogram and throws.)
test_almost("TC-17 2x2 uniform image [5]      → 0", otsuthresh(5*ones(2,2)), 0, 0);

// =============================================================================
// SECTION 4 - Consistency: Image Mode == Histogram Mode
// =============================================================================
mprintf("\n--- SECTION 4: Image/Histogram Consistency ---\n");

// TC-18  Build a histogram from a known image manually, then compare
//        otsuthresh(image) with otsuthresh(manual_histogram).
//        They must agree to floating-point precision.
I18 = [zeros(10,10); ones(10,10) * 0.8];   // 100 pixels at 0, 100 at 0.8
nbins = 256;
J18 = (I18 - min(I18)) / (max(I18) - min(I18));
idx18 = min(floor(J18(:) * (nbins-1)) + 1, nbins);
h18 = zeros(nbins, 1);
for b = 1:nbins
    h18(b) = sum(idx18 == b);
end
t18_img  = otsuthresh(I18);
t18_hist = otsuthresh(h18);
test_almost("TC-18 img mode == hist mode", t18_img, t18_hist, 1e-12);

// =============================================================================
// SECTION 5 - Error Paths (all must THROW)
// =============================================================================
mprintf("\n--- SECTION 5: Error Paths ---\n");

// TC-19  Zero arguments.
test_throws("TC-19 zero args          → throw", "otsuthresh();");

// TC-20  Two arguments.
test_throws("TC-20 two args           → throw", "otsuthresh([1 2],[3 4]);");

// TC-21  Complex input.
test_throws("TC-21 complex input      → throw", "otsuthresh([1+%i 2]);");

// TC-22  Histogram with only 1 bin.
test_throws("TC-22 1-bin histogram    → throw", "otsuthresh([100]);");

// TC-23  Negative histogram counts.
test_throws("TC-23 negative counts    → throw", "otsuthresh([-1; 10; 5]);");

// TC-24  All-zero histogram (sum = 0).
test_throws("TC-24 zero-sum histogram → throw", "otsuthresh(zeros(10,1));");

// TC-25  3-D array input (neither vector nor 2-D matrix).
test_throws("TC-25 3-D array input    → throw", "otsuthresh(ones(3,3,3));");

// =============================================================================
// SUMMARY
// =============================================================================
mprintf("\n");
mprintf("----------------------------------------------------------------\n");
mprintf("RESULT: otsuthresh test suite PASSED  (25 tests)\n");
mprintf("----------------------------------------------------------------\n");

// =============================================================================
// END OF otsuthresh_test.sce
// =============================================================================
