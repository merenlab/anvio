#!/usr/bin/env python
# tests for the distribution of coverage (discov) calculation in utils.CoverageStats. every test builds a small
# coverage array for which all intermediate values can be computed by hand, and compares them to what
# CoverageStats reports. B (the number of aligned bases) is always the sum of the coverage array, and unless a
# test says otherwise, N (the number of reads) is set to match 100 bp reads, i.e. N = B / 100. every test that
# reaches C and discov checks them with both the linear and the geometric formula (alpha = 0.5)

import math
import numpy as np

import anvio.utils as utils
import anvio.terminal as terminal

from anvio.errors import ConfigError

run = terminal.Run(width=50)

num_tests = 0


def make_coverage(G, regions):
    """Returns a coverage array of length G, with each (start, stop, depth) region set to depth"""

    coverage = np.zeros(G, dtype=int)
    for start, stop, depth in regions:
        coverage[start:stop] = depth

    return coverage


def get_stats(coverage, N, w=None, pct=None, min_w=None, max_w=None, alpha=0.5, formula='geometric', min_reads=0):
    return utils.CoverageStats(coverage, skip_outliers=True, discov_window_length=w, discov_window_percentage=pct,
                               discov_min_window_len=min_w, discov_max_window_len=max_w, discov_alpha=alpha, discov_formula=formula,
                               discov_min_reads=min_reads,
                               num_reads=N, report_discov=True)


def check(label, observed, expected, tol=1e-9):
    if expected is None or observed is None:
        passed = observed is expected
        run.info(label, f"{observed} (expected {expected})", mc='green' if passed else 'red')
    elif isinstance(expected, int):
        passed = observed == expected
        run.info(label, f"{observed} (expected {expected})", mc='green' if passed else 'red')
    else:
        passed = math.isclose(observed, expected, rel_tol=0, abs_tol=tol)
        run.info(label, f"{observed:.6f} (expected {expected:.6f})", mc='green' if passed else 'red')

    assert passed, f"{label}: observed {observed}, expected {expected}"


def check_all(stats, expected, context=None):
    for attribute, label in [('num_windows', 'num windows'),
                             ('prop_win_covered', 'S'),
                             ('spread_expected', 'S_exp'),
                             ('spread_norm', 'S_norm'),
                             ('window_depth_evenness', 'E'),
                             ('detection', 'D'),
                             ('spread_evenness_score', 'C'),
                             ('discov', 'discov')]:
        if attribute in expected:
            check(f"{label} ({context})" if context else label, getattr(stats, attribute), expected[attribute])


def expected_C(S_norm, E, formula, alpha=0.5):
    if formula == 'linear':
        return alpha * S_norm + (1 - alpha) * E

    return S_norm ** alpha * E ** (1 - alpha)


def check_formulas(get_stats_for, S_norm, E, D, alpha=0.5, show_alpha=False, context=None):
    """Checks C and discov for both formulas. get_stats_for(formula) returns a CoverageStats object"""

    for formula in ['linear', 'geometric']:
        stats = get_stats_for(formula)
        C = expected_C(S_norm, E, formula, alpha)
        label = f"{formula}, alpha = {alpha}" if show_alpha or alpha != 0.5 else formula
        if context:
            label = f"{label}, {context}"
        check(f"C ({label})", stats.spread_evenness_score, C)
        check(f"discov ({label})", stats.discov, min(1.0, C + D ** 2))


def test(title, description):
    global num_tests
    num_tests += 1
    run.warning(description, header=f"TEST {num_tests}: {title}", lc='cyan')


#########################################

test("Uniform coverage everywhere",
     "G = 10,000, w = 1,000, depth 10 at every position, N = 1,000 reads (so B = 100,000). All 10 windows "
     "are covered (S = 1) at the same depth (CV = 0, so E = 1). S_exp = 1 - exp(-(100,000 + 1,000 x 1,000) / 10,000) "
     "= 1 - exp(-110), which is 1 for all practical purposes. C = 1 with both formulas (0.5 x 1 + 0.5 x 1 linear, "
     "1^0.5 x 1^0.5 geometric), and discov = min(1, C + D^2) = min(1, 2) = 1 (capped).")

uniform = make_coverage(10000, [(0, 10000, 10)])
check_all(get_stats(uniform, N=1000, w=1000), {'num_windows': 10,
                                               'prop_win_covered': 1.0,
                                               'spread_expected': 1 - math.exp(-110),
                                               'spread_norm': 1.0,
                                               'window_depth_evenness': 1.0,
                                               'detection': 1.0})
check_formulas(lambda formula: get_stats(uniform, N=1000, w=1000, formula=formula), S_norm=1.0, E=1.0, D=1.0)

#########################################

E_uneven = 1 / (1 + np.std([5, 10, 20, 50, 100]) / 37)
S_norm_half = 0.5 / (1 - math.exp(-203.5))

test("Half of the windows covered at uneven depths",
     "G = 10,000, w = 1,000, windows 1, 3, 5, 7 and 9 are fully covered at depths 5, 10, 20, 50 and 100, and the "
     "others have no coverage. B = 185,000 and N = 1,850. S = 5 / 10 = 0.5. S_exp = 1 - exp(-(185,000 + 1,850 x "
     "1,000) / 10,000) = 1 - exp(-203.5) ≈ 1, so S_norm = 0.5. E only considers the nonzero window means: mean = 37, "
     f"std = sqrt(1,236) = 35.157, CV = 0.9502, E = 1 / (1 + 0.9502) = {E_uneven:.6f}. D = 0.5. Linear: C = 0.5 x 0.5 + "
     f"0.5 x E = {0.25 + 0.5 * E_uneven:.6f}, discov = C + 0.5^2 = {0.25 + 0.5 * E_uneven + 0.25:.6f}. Geometric: C = "
     f"sqrt(0.5 x E) = {math.sqrt(0.5 * E_uneven):.6f}, discov = {math.sqrt(0.5 * E_uneven) + 0.25:.6f}.")

half_covered = make_coverage(10000, [(0, 1000, 5), (2000, 3000, 10), (4000, 5000, 20), (6000, 7000, 50), (8000, 9000, 100)])
check_all(get_stats(half_covered, N=1850, w=1000), {'num_windows': 10,
                                                    'prop_win_covered': 0.5,
                                                    'spread_expected': 1 - math.exp(-203.5),
                                                    'spread_norm': S_norm_half,
                                                    'window_depth_evenness': E_uneven,
                                                    'detection': 0.5})
check_formulas(lambda formula: get_stats(half_covered, N=1850, w=1000, formula=formula), S_norm=S_norm_half, E=E_uneven, D=0.5)

#########################################

S_exp_alpha = 1 - math.exp(-(185000 + 1850 * 1000) / 5000)

test("Changing alpha",
     "G = 5,000, w = 1,000, and each of the 5 windows has a 100 bp stretch covered at depths 50, "
     "100, 200, 500 and 1,000, so the window means are again 5, 10, 20, 50 and 100 (E = "
     f"{E_uneven:.6f}), but every window is covered (S = 1). B = 185,000, N = 1,850, S_exp = 1 - exp(-407) ≈ 1, so "
     "S_norm = 1. D = 0.1, so discov = C + 0.01 (capped at 1). As alpha goes from 0 to 1, C should move from E alone to "
     "S_norm alone, with both formulas: linear C = alpha + (1 - alpha) x E, and geometric C = E^(1 - alpha). With alpha "
     f"= 0.99, E has almost no influence: linear C = {0.99 + 0.01 * E_uneven:.6f}, geometric C = "
     f"{E_uneven ** 0.01:.6f}, and discov = 1 with both.")

all_windows_uneven = make_coverage(5000, [(0, 100, 50), (1000, 1100, 100), (2000, 2100, 200), (3000, 3100, 500), (4000, 4100, 1000)])
check_all(get_stats(all_windows_uneven, N=1850, w=1000), {'num_windows': 5,
                                                          'prop_win_covered': 1.0,
                                                          'spread_expected': S_exp_alpha,
                                                          'spread_norm': 1.0,
                                                          'window_depth_evenness': E_uneven,
                                                          'detection': 0.1})
for alpha in [0, 0.5, 0.9, 0.99, 1]:
    check_formulas(lambda formula: get_stats(all_windows_uneven, N=1850, w=1000, alpha=alpha, formula=formula), S_norm=1.0, E=E_uneven, D=0.1, alpha=alpha, show_alpha=True)

#########################################

S_exp_short = 1 - math.exp(-(2000 + 20 * 1000) / 100000)
S_norm_short = 0.1 / S_exp_short

test("Sparse short reads, where the expected spread matters",
     "G = 100,000, w = 1,000, and 20 reads of 100 bp (N = 20, B = 2,000). They are placed two per window in 10 "
     f"windows, so S = 10 / 100 = 0.1. S_exp = 1 - exp(-(2,000 + 20 x 1,000) / 100,000) = {S_exp_short:.6f}, "
     f"S_norm = 0.1 / S_exp = {S_norm_short:.6f}. All covered windows have the same mean (200 / 1,000 = 0.2), "
     f"so E = 1. D = 2,000 / 100,000 = 0.02. Linear: C = 0.5 x S_norm + 0.5 = {0.5 * S_norm_short + 0.5:.6f}, "
     f"discov = C + 0.02^2 = {0.5 * S_norm_short + 0.5 + 0.0004:.6f}. Geometric: C = sqrt(S_norm) = "
     f"{math.sqrt(S_norm_short):.6f}, discov = {math.sqrt(S_norm_short) + 0.0004:.6f}.")

regions = []
for window in range(10):
    regions.extend([(window * 1000, window * 1000 + 100, 1), (window * 1000 + 500, window * 1000 + 600, 1)])
sparse_short = make_coverage(100000, regions)
check_all(get_stats(sparse_short, N=20, w=1000), {'num_windows': 100,
                                                  'prop_win_covered': 0.1,
                                                  'spread_expected': S_exp_short,
                                                  'spread_norm': S_norm_short,
                                                  'window_depth_evenness': 1.0,
                                                  'detection': 0.02})
check_formulas(lambda formula: get_stats(sparse_short, N=20, w=1000, formula=formula), S_norm=S_norm_short, E=1.0, D=0.02)

#########################################

S_exp_long = 1 - math.exp(-(2000 + 2 * 1000) / 100000)
S_norm_long = 0.02 / S_exp_long

test("Same number of aligned bases, but fewer (longer) reads",
     "G = 100,000, w = 1,000, and 2 reads of 1,000 bp, each filling one window (N = 2, B = 2,000, same as the "
     f"previous test). S = 2 / 100 = 0.02. S_exp = 1 - exp(-(2,000 + 2 x 1,000) / 100,000) = {S_exp_long:.6f}, "
     f"which is much lower than with 20 short reads, so S_norm = 0.02 / S_exp = {S_norm_long:.6f}. Only 2 "
     "windows are covered, so E = 0 (see the test about fewer than 3 covered windows). D = 0.02. Linear: C = 0.5 x "
     f"S_norm = {0.5 * S_norm_long:.6f}, discov = {0.5 * S_norm_long + 0.0004:.6f}. Geometric: C = 0 because E = 0, "
     "so discov = D^2 = 0.0004.")

long_reads = make_coverage(100000, [(0, 1000, 1), (50000, 51000, 1)])
check_all(get_stats(long_reads, N=2, w=1000), {'prop_win_covered': 0.02,
                                               'spread_expected': S_exp_long,
                                               'spread_norm': S_norm_long,
                                               'window_depth_evenness': 0.0,
                                               'detection': 0.02})
check_formulas(lambda formula: get_stats(long_reads, N=2, w=1000, formula=formula), S_norm=S_norm_long, E=0.0, D=0.02)

#########################################

S_exp_pacbio = 1 - math.exp(-(50000 + 5 * 1000) / 200000)
S_exp_pacbio_as_short = 1 - math.exp(-(50000 + 500 * 1000) / 200000)
S_norm_pacbio_as_short = 0.25 / S_exp_pacbio_as_short

test("Long reads (PacBio, mean read length 10,000 bp)",
     "G = 200,000, w = 1,000, and 5 PacBio reads of 6,000, 8,000, 10,000, 12,000 and 14,000 bp (mean 10,000 bp) at "
     "depth 1, starting at window boundaries and not overlapping (N = 5, B = 50,000). They cover 6 + 8 + 10 + 12 + 14 "
     "= 50 windows, so S = 50 / 200 = 0.25. S_exp = 1 - exp(-(50,000 + 5 x 1,000) / 200,000) = "
     f"{S_exp_pacbio:.6f}, which is lower than S, so S_norm = 1. All covered windows have depth 1, so E = 1. "
     "D = 0.25. C = 1 with both formulas and discov = min(1, 1 + 0.25^2) = 1. Since each long read covers many "
     "neighboring windows, covering only a quarter of the windows is what we expect from 5 reads, and the score is not "
     "penalized. If the same coverage came from 100 bp reads (N = 500), S_exp = 1 - exp(-(50,000 + 500 x 1,000) / "
     f"200,000) = {S_exp_pacbio_as_short:.6f}, so S_norm = {S_norm_pacbio_as_short:.6f}. Linear: C = "
     f"{0.5 * S_norm_pacbio_as_short + 0.5:.6f}, discov = {0.5 * S_norm_pacbio_as_short + 0.5 + 0.0625:.6f}. "
     f"Geometric: C = sqrt(S_norm) = {math.sqrt(S_norm_pacbio_as_short):.6f}, discov = "
     f"{math.sqrt(S_norm_pacbio_as_short) + 0.0625:.6f}.")

regions, start = [], 0
for read_length in [6000, 8000, 10000, 12000, 14000]:
    regions.append((start, start + read_length, 1))
    start += read_length + 30000
pacbio = make_coverage(200000, regions)
run.info_single("Long reads (N = 5)", mc='cyan', nl_after=1)
check_all(get_stats(pacbio, N=5, w=1000), {'num_windows': 200,
                                           'prop_win_covered': 0.25,
                                           'spread_expected': S_exp_pacbio,
                                           'spread_norm': 1.0,
                                           'window_depth_evenness': 1.0,
                                           'detection': 0.25})
check_formulas(lambda formula: get_stats(pacbio, N=5, w=1000, formula=formula), S_norm=1.0, E=1.0, D=0.25)
run.info_single("Same coverage from 100 bp reads (N = 500)", mc='cyan', nl_before=1, nl_after=1)
check_all(get_stats(pacbio, N=500, w=1000), {'spread_expected': S_exp_pacbio_as_short,
                                             'spread_norm': S_norm_pacbio_as_short})
check_formulas(lambda formula: get_stats(pacbio, N=500, w=1000, formula=formula), S_norm=S_norm_pacbio_as_short, E=1.0, D=0.25)

#########################################

test("Observed spread higher than expected (S_norm capped at 1)",
     "G = 100,000, w = 1,000, and 20 reads of 100 bp, each in a different window (N = 20, B = 2,000). S = 0.2, which "
     f"is higher than S_exp = {S_exp_short:.6f}, so S_norm = min(1, 0.2 / S_exp) = 1. All covered windows have the "
     "same mean (0.1), so E = 1. D = 0.02. C = 1 with both formulas, and discov = min(1, 1 + 0.02^2) = 1.")

capped = make_coverage(100000, [(window * 1000, window * 1000 + 100, 1) for window in range(20)])
check_all(get_stats(capped, N=20, w=1000), {'prop_win_covered': 0.2,
                                            'spread_expected': S_exp_short,
                                            'spread_norm': 1.0,
                                            'window_depth_evenness': 1.0,
                                            'detection': 0.02})
check_formulas(lambda formula: get_stats(capped, N=20, w=1000, formula=formula), S_norm=1.0, E=1.0, D=0.02)

#########################################

test("No coverage at all",
     "G = 5,000, w = 1,000, no coverage and N = 0. Then S = 0, S_exp = 1 - exp(0) = 0 (and S_norm is set to 0 "
     "rather than dividing by 0), E = 0, D = 0, and C and discov are 0 with both formulas.")

no_coverage = make_coverage(5000, [])
check_all(get_stats(no_coverage, N=0, w=1000), {'num_windows': 5,
                                                'prop_win_covered': 0.0,
                                                'spread_expected': 0.0,
                                                'spread_norm': 0.0,
                                                'window_depth_evenness': 0.0,
                                                'detection': 0.0})
check_formulas(lambda formula: get_stats(no_coverage, N=0, w=1000, formula=formula), S_norm=0.0, E=0.0, D=0.0)

#########################################

test("Sequence shorter than the window length",
     "G = 500, w = 1,000, with depth 1 over the first 100 bases and N = 1 (B = 100). The whole sequence is a single "
     "window, and the window length used in S_exp is the length of that window (500), not the requested one: S_exp "
     f"= 1 - exp(-(100 + 1 x 500) / 500) = {1 - math.exp(-1.2):.6f} (it would be {1 - math.exp(-2.2):.6f} with "
     "w = 1,000). S = 1, S_norm = 1, E = 0 (only one window), D = 0.2. Linear: C = 0.5 x 1 + 0.5 x 0 = 0.5, discov = "
     "0.5 + 0.2^2 = 0.54. Geometric: C = 0 because E = 0, so discov = 0.04.")

short_sequence = make_coverage(500, [(0, 100, 1)])
check_all(get_stats(short_sequence, N=1, w=1000), {'num_windows': 1,
                                                   'prop_win_covered': 1.0,
                                                   'spread_expected': 1 - math.exp(-1.2),
                                                   'spread_norm': 1.0,
                                                   'window_depth_evenness': 0.0,
                                                   'detection': 0.2})
check_formulas(lambda formula: get_stats(short_sequence, N=1, w=1000, formula=formula), S_norm=1.0, E=0.0, D=0.2)

#########################################

test("Fewer than 3 covered windows",
     "G = 10,000, w = 1,000, N = 300 (S_exp = 1 - exp(-33) ≈ 1). With 2 windows fully covered (depths 10 and 20), E "
     "is set to 0. S_norm = 0.2 and D = 0.2, so the linear C = 0.1 (discov = 0.14) "
     "and the geometric C = 0 (discov = 0.04). With a third window covered at the same depth as the others (10, 10, "
     "10), CV = 0 and E = 1. S_norm = 0.3 and D = 0.3, so the linear C = 0.65 (discov = 0.74) and the geometric C = "
     f"sqrt(0.3) = {math.sqrt(0.3):.6f} (discov = {math.sqrt(0.3) + 0.09:.6f}).")

two_windows = make_coverage(10000, [(0, 1000, 10), (5000, 6000, 20)])
check("E (2 covered windows)", get_stats(two_windows, N=300, w=1000).window_depth_evenness, 0.0)
check_formulas(lambda formula: get_stats(two_windows, N=300, w=1000, formula=formula), S_norm=0.2 / (1 - math.exp(-33)), E=0.0, D=0.2)
three_windows = make_coverage(10000, [(0, 1000, 10), (4000, 5000, 10), (8000, 9000, 10)])
check("E (3 covered windows)", get_stats(three_windows, N=300, w=1000).window_depth_evenness, 1.0)
check_formulas(lambda formula: get_stats(three_windows, N=300, w=1000, formula=formula), S_norm=0.3 / (1 - math.exp(-33)), E=1.0, D=0.3)

#########################################

S_exp_final = 1 - math.exp(-(500 + 5 * 1000) / 10100)

test("Short final window",
     "w = 1,000. With G = 10,050, the final window is 50 bp (< 10% of w), so it is dropped and there are 10 windows. "
     "Since only that final window has coverage, S = 0 although D > 0, and C = 0 with both formulas (discov = D^2). "
     "With G = 10,100, the final window is 100 bp (exactly 10% of w), so it is kept and there are 11 windows, one of "
     f"which is covered (S = 1 / 11). S_exp = 1 - exp(-(500 + 5 x 1,000) / 10,100) = {S_exp_final:.6f}, so S_norm = "
     f"{1 / 11 / S_exp_final:.6f}, and E = 0 (one covered window). Linear: C = 0.5 x S_norm = "
     f"{0.5 / 11 / S_exp_final:.6f}. Geometric: C = 0. In both cases the coverage comes from 5 reads (N = 5) that "
     "overlap the end of the sequence.")

dropped_final = make_coverage(10050, [(10000, 10050, 5)])
check_all(get_stats(dropped_final, N=5, w=1000), {'num_windows': 10,
                                                  'prop_win_covered': 0.0,
                                                  'detection': 50 / 10050})
check_formulas(lambda formula: get_stats(dropped_final, N=5, w=1000, formula=formula), S_norm=0.0, E=0.0, D=50 / 10050)
kept_final = make_coverage(10100, [(10000, 10100, 5)])
check_all(get_stats(kept_final, N=5, w=1000), {'num_windows': 11,
                                               'prop_win_covered': 1 / 11,
                                               'spread_expected': S_exp_final})
check_formulas(lambda formula: get_stats(kept_final, N=5, w=1000, formula=formula), S_norm=1 / 11 / S_exp_final, E=0.0, D=100 / 10100)

#########################################

test("Window length as a percentage of the sequence length",
     "G = 4,000 with uniform coverage and a window length of 5%, which is 200 bp. Without a minimum window length, "
     "there are 4,000 / 200 = 20 windows. With a minimum of 300 bp, the windows are raised to 300 bp: 13 full windows "
     "plus a final window of 100 bp (>= 10% of 300), so 14 windows. In both cases S = E = D = 1, so C and discov are "
     "1 with both formulas.")

coverage = make_coverage(4000, [(0, 4000, 5)])
check("num windows (no minimum)", get_stats(coverage, N=200, pct=5, min_w=None).num_windows, 20)
check_formulas(lambda formula: get_stats(coverage, N=200, pct=5, min_w=None, formula=formula), S_norm=1.0, E=1.0, D=1.0)
check("num windows (minimum 300)", get_stats(coverage, N=200, pct=5, min_w=300).num_windows, 14)
check_formulas(lambda formula: get_stats(coverage, N=200, pct=5, min_w=300, formula=formula), S_norm=1.0, E=1.0, D=1.0)

#########################################

test("Maximum window length for percentage-based windows",
     "G = 4,000 with uniform coverage and a window length of 50%, which is 2,000 bp, so there are 2 windows. With a "
     "maximum of 300 bp, the windows are capped at 300 bp: 13 full windows plus a final window of 100 bp (>= 10% of "
     "300), so 14 windows. With a maximum larger than the percentage-based window length, nothing changes. In all "
     "cases S = E = D = 1, so C and discov are 1 with both formulas.")

coverage = make_coverage(4000, [(0, 4000, 5)])
check("num windows (no maximum)", get_stats(coverage, N=200, pct=50).num_windows, 2)
check("num windows (maximum 300)", get_stats(coverage, N=200, pct=50, max_w=300).num_windows, 14)
check_formulas(lambda formula: get_stats(coverage, N=200, pct=50, max_w=300, formula=formula), S_norm=1.0, E=1.0, D=1.0)
check("num windows (maximum 5000)", get_stats(coverage, N=200, pct=50, max_w=5000).num_windows, 2)
check("num windows (minimum 100, maximum 300)", get_stats(coverage, N=200, pct=50, min_w=100, max_w=300).num_windows, 14)

#########################################

test("Minimum number of mapped reads",
     "G = 10,000 with w = 1,000, and 4 reads of 100 bp each landing in separate windows (N = 4, B = 400). Without a "
     "read threshold, S = 4 / 10 = 0.4 and S_exp = 1 - e^(-(400 + 4,000) / 10,000) = 1 - e^(-0.44), so S_norm = "
     "min(1, 0.4 / S_exp), and E = 1 because all covered windows have the same mean coverage. With the threshold "
     "of 10 reads, S_norm and E are set to 0 (while S and S_exp are still reported), so C = 0 and discov = D^2 with "
     "both formulas. With the threshold at exactly N = 4, S_norm and E are computed as usual.")

coverage = make_coverage(10000, [(0, 100, 1), (2000, 2100, 1), (5000, 5100, 1), (8000, 8100, 1)])
S_exp_few = 1 - math.exp(-0.44)
S_norm_few = min(1, 0.4 / S_exp_few)
D_few = 400 / 10000

for min_reads, context, S_norm, E in [(0, 'no read threshold', S_norm_few, 1.0),
                                      (10, 'min reads = 10', 0.0, 0.0),
                                      (4, 'min reads = 4', S_norm_few, 1.0)]:
    check_all(get_stats(coverage, N=4, w=1000, min_reads=min_reads), {'num_windows': 10,
                                                                      'prop_win_covered': 0.4,
                                                                      'spread_expected': S_exp_few,
                                                                      'spread_norm': S_norm,
                                                                      'window_depth_evenness': E,
                                                                      'detection': D_few},
              context=context)
    check_formulas(lambda formula: get_stats(coverage, N=4, w=1000, min_reads=min_reads, formula=formula),
                   S_norm=S_norm, E=E, D=D_few, context=context)

#########################################

test("Discov attributes when discov is not requested, and missing N",
     "Without report_discov, all discov attributes should be None. With report_discov but without the number of "
     "reads, CoverageStats should raise a ConfigError, with either formula.")

stats = utils.CoverageStats(make_coverage(1000, [(0, 1000, 1)]), skip_outliers=True)
for attribute in ['num_windows', 'prop_win_covered', 'spread_expected', 'spread_norm', 'window_depth_evenness', 'spread_evenness_score', 'discov']:
    check(attribute, getattr(stats, attribute), None)

for formula in ['linear', 'geometric']:
    try:
        get_stats(make_coverage(1000, [(0, 1000, 1)]), N=None, w=100, formula=formula)
        raised = False
    except ConfigError:
        raised = True
    run.info(f'ConfigError raised without N ({formula})', raised, mc='green' if raised else 'red')
    assert raised, f"CoverageStats did not raise a ConfigError when N was missing ({formula} formula)"

#########################################

run.info_single(f"All {num_tests} discov tests passed.", mc='green', nl_before=1, nl_after=1)
