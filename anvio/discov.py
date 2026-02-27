# coding: utf-8
# pylint: disable=line-too-long
"""A class to compute distribution of coverage by evaluating various metrics for 'evenness' of breadth/depth."""

import numpy as np

import anvio
import anvio.utils as utils
import anvio.filesnpaths as filesnpaths

from anvio.terminal import Run
from anvio.errors import ConfigError

# use --debug flag to allow things to be written to terminal
run = Run()
run.verbose = False

class DisCov:
    """Computes and reports various metrics for coverage evenness."""

    def __init__(self, coverage, contig_name=None, sample_name=None):
        self.median: float = np.median(coverage)
        self.mean: float = np.mean(coverage)
        self.std: float = np.std(coverage)
        self.detection: float = np.sum(coverage > 0) / len(coverage)
        self.name = contig_name # optional parameter, only used for debug output
        self.sample = sample_name

        # establish output files for testing
        unfilt_output = "TEST_UNFILTERED.txt"
        filt_output = "TEST_FILTERED.txt"
        header = ["contig", "sample", "Non-specific Filter Removed Bases", 
                  "Midpoint Range", "SW Depth Evenness MAD (fine)",
                  "SW Depth Evenness MAD (medium)", "SW Depth Evenness MAD (coarse)",
                  "SW Depth Evenness CV (fine)", "SW Depth Evenness CV (medium)",
                  "SW Depth Evenness CV (coarse)", "SW Proportion Covered (fine)",
                  "SW Proportion Covered (medium)", "SW Proportion Covered (coarse)",
                  "Window-Scaling Variance", "Dispersion of Counts (num bins = 50)",
                  "Dispersion of Counts (num bins = 30)", "Dispersion of Counts (num bins = 10)",
                  "Shannon Entropy Evenness"]
        for outfile in [unfilt_output, filt_output]:
            if not filesnpaths.is_file_exists(outfile, dont_raise=True):
                with open(outfile, 'w') as f:
                    f.write("\t".join(header) + "\n")

        # compute all sub-metrics on regular coverage array
        self.compute_all(coverage, unfilt_output, filter_nonspecific_mapping=False)
        # compute all sub-metrics on filtered coverage array
        self.compute_all(coverage, filt_output, filter_nonspecific_mapping=True)

    def compute_all(self, cov_array, output_file, filter_nonspecific_mapping: bool = False):
        """Compute all sub-metrics for the contig and print to the specified output file"""

        # skip any contigs without any mapping for which all values would be 0
        if self.detection == 0:
            # append all metrics to file
            output_list = [self.name, self.sample,
                        f"0", # filter removed 0 bases
                        f"0", # midpoint range is 0
                        "NA\tNA\tNA", # sliding window evenness MAD, NA because we can't divide by median of 0
                        "NA\tNA\tNA", # sliding window evenness CV, NA because we can't divide by mean of 0
                        "0\t0\t0", # sliding window proportion covered
                        "NA", # window-scaling variance, NA because we can't take log of 0 variance
                        "NA", # distribution of counts in bins, NA because we can't divide by mean of 0
                        "NA" # shannon entropy evenness, NA because we have no data to compute entropy on
                      ]
            with open(output_file, 'a') as f:
                f.write("\t".join(output_list) + "\n")
            return

        run.warning(f"Anvi'o is now attempting to compute how coverage is distributed across a contig. The "
                    f"terminal output below relates to this calculation.", header="Computing Distribution of Coverage", 
                    lc='green', overwrite_verbose=anvio.DEBUG)
        run.info(f"Contig name", self.name if self.name else 'Unknown', overwrite_verbose=anvio.DEBUG)
        run.info(f"Sample", self.sample if self.sample else 'Unknown', overwrite_verbose=anvio.DEBUG)
        run.info("Filtering out non-specific read recruitment", filter_nonspecific_mapping, overwrite_verbose=anvio.DEBUG)
        run.info("Detection of contig", self.detection, overwrite_verbose=anvio.DEBUG)

        # this will be a list of tuples like (start_position, stop_position, mean_coverage)
        # if mean_coverage is 0, then the region is a gap region. Otherwise, it is a covered region
        regions = self.get_list_of_coverage_and_gap_regions(cov_array, min_gap_length=3)

        # filter out regions with unusually-high coverage that might be non-specific read recruitment
        # these regions and their surrounding gaps get replaced with one longer gap
        ns_filter_removed_something = 0
        filtered_detection = self.detection
        if filter_nonspecific_mapping:
            regions, new_cov_array = self.filter_nonspecific_regions(cov_array, regions, mod_z_score_threshold_external=4, mod_z_score_threshold_internal=3)
            filtered_detection = np.sum(new_cov_array > 0) / len(new_cov_array)
            run.info("Detection of contig (post-filter)", filtered_detection, overwrite_verbose=anvio.DEBUG)
            if not np.array_equal(new_cov_array, cov_array):
                ns_filter_removed_something = len(cov_array) - len(new_cov_array)
                cov_array = new_cov_array

        # compute all metrics
        ## region metrics
        num_coverage_regions = len([r for r in regions if r[2] > 0])
        num_gaps = len([r for r in regions if r[2] == 0])
        gap_gini = self.gap_evenness_gini(regions, min_num_gaps_for_gini=10)
        mp_range, mp_evenness = self.midpoint_spread_and_evenness(regions)
        ## window metrics
        window_scales = {
            'fine': max(100, len(cov_array) // 100),      # ~1% of contig
            'medium': max(500, len(cov_array) // 20),     # ~5% of contig  
            'coarse': max(1000, len(cov_array) // 10)     # ~10% of contig
        }
        sliding_window_metrics = self.sliding_window_evenness(cov_array, window_scales)
        sw_mad_vals = [sliding_window_metrics[scale]['Depth_Evenness_MAD'] for scale in ['fine','medium','coarse']]
        sw_mad = [f"{m:.04}" if m else "NA" for m in sw_mad_vals ]
        sliding_window_evenness_mad = "\t".join(sw_mad)
        sw_cv_vals = [sliding_window_metrics[scale]['Depth_Evenness_CV'] for scale in ['fine','medium','coarse']]
        sw_cv = [f"{m:.04}" if m else "NA" for m in sw_cv_vals ]
        sliding_window_evenness_cv = "\t".join(sw_cv)
        sw_prop_vals = [sliding_window_metrics[scale]['Proportion_Covered'] for scale in ['fine','medium','coarse']]
        sw_prop = [f"{m:.04}" if m else "NA" for m in sw_prop_vals ]
        sliding_window_proportion_covered = "\t".join(sw_prop)
        window_scaling_variance = self.compute_window_scaling_variance(cov_array)
        disp_50 = self.binned_count_dispersion(cov_array, num_windows=50) # fine-scale clustering, small windows
        disp_30 = self.binned_count_dispersion(cov_array, num_windows=30)
        disp_10 = self.binned_count_dispersion(cov_array, num_windows=10) # coarse-scale clustering, large windows
        disp_all = [f"{m:.04}" if m else "NA" for m in [disp_50, disp_30, disp_10]]
        disp_counts = "\t".join(disp_all)

        ## whole-contig metrics
        shannon = self.Shannon_entropy_evenness(cov_array)

        # append all metrics to file
        output_list = [self.name, self.sample,
                        f"{ns_filter_removed_something}",
                        f"{mp_range:.4}" if mp_range else "NA",
                        sliding_window_evenness_mad,
                        sliding_window_evenness_cv,
                        sliding_window_proportion_covered,
                        f"{window_scaling_variance:.4}" if window_scaling_variance else "NA",
                        disp_counts,
                        f"{shannon:.4}" if shannon else "NA"
                      ]
        with open(output_file, 'a') as f:
            f.write("\t".join(output_list) + "\n")


    ##### METRIC FUNCTIONS #####
    def gap_evenness_gini(self, region_list, min_num_gaps_for_gini=10):
        """Evenness of gap lengths. E = 1 - Gini(gap lengths)."""
        gaplens = []
        for start, stop, meancov in region_list:
            if meancov == 0:
                gaplens.append(stop - start)
        ngaps = len(gaplens)

        if ngaps >= min_num_gaps_for_gini:
            G = self.compute_gini_coeff(gaplens)
            return 1 - G
        else:
            return None
    
    def sliding_window_evenness(self, coverage, window_len_dict):
        """Evenness of coverage depth and proportion of covered windows at multiple window lengths.
        Given a dictionary of labeled window lengths, returns a dictionary of values for each window length:
            Depth_Evenness_MAD : 1 / (1 + (MAD / median))
            Depth_Evenness_CV : 1 / (1 + CV); where CV = std-dev / mean
            Proportion_Covered : # windows with nonzero mean coverage / # windows

            Note that the 1 / (1 + x) is a strategy to get an evenness score between 0 and 1.
            Note also that we don't worry about dividing by zero because we are already working with nonzero means.
        """

        results = {scale_name : {} for scale_name in window_len_dict}
        for scale_name, wlen in window_len_dict.items():
            results[scale_name]['Depth_Evenness_MAD'] = None
            results[scale_name]['Depth_Evenness_CV'] = None
            results[scale_name]['Proportion_Covered'] = None
            if wlen > len(coverage):
                continue

            windows = self.get_sliding_window_regions(coverage, wlen)
            nonzero_window_means = [x for x in windows if x[2] > 0]
            median_depth, depth_MAD = self.compute_median_and_MAD(nonzero_window_means) # will be None if <3 windows
            cv_of_means = self.compute_CV(nonzero_window_means) # will be None if mean of means is None

            if median_depth:
                mad_dispersion = depth_MAD / median_depth
                results[scale_name]['Depth_Evenness_MAD'] = 1 / (1 + mad_dispersion)
            if cv_of_means:
                results[scale_name]['Depth_Evenness_CV'] = 1 / (1 + cv_of_means)
            if windows:
                results[scale_name]['Proportion_Covered'] = len(nonzero_window_means) / len(windows)
        return results

    def compute_window_scaling_variance(self, coverage, max_window_fraction=0.10, num_window_sizes=20):
        """Computes the 'window-scaling variance' slope from an array of coverage values.

        Divides the coverage array into non-overlapping windows of logarithmically-spaced
        sizes ranging from max(100, 1% of array length) to max_window_fraction * array length.
        For each window size, computes the variance of the per-window mean coverages.
        Returns the slope of log(variance) vs log(window_size), fit by linear regression. The closer 
        to -1 this slope is, the more stationary the coverage distribution. The closer to 0, the more 
        sparse it is.

        Parameters
        ==========
        coverage : array-like
            Array of read recruitment coverage values.
        max_window_fraction : float
            Upper bound on window size as a fraction of the array length. Default 0.10
            (i.e., 10%), which yields ~10 windows at the largest scale. Going above ~15%
            risks too few windows for reliable variance estimates.
        num_window_sizes : int
            Number of logarithmically-spaced window sizes to evaluate. Default 20.

        Returns
        =======
        float
            The scaling slope (i.e., the exponent relating variance to window size).
        """

        L = len(coverage)
        min_w = max(100, int(0.01 * L))
        max_w = int(max_window_fraction * L)

        if min_w >= max_w:
            raise ConfigError(f"compute_window_scaling_variance() could not construct a valid window size range: "
                                f"min window size ({min_w}) >= max window size ({max_w}). Your coverage array "
                                f"(length {L}) may be too short, or max_window_fraction ({max_window_fraction}) too small.")

        # Logarithmically spaced window sizes, deduplicated after rounding to integers
        window_sizes = np.unique(
            np.round(np.geomspace(min_w, max_w, num_window_sizes)).astype(int)
        )

        log_w_vals = []
        log_var_vals = []

        for w in window_sizes:
            windows = self.get_sliding_window_regions(coverage, w)
            means = [region[2] for region in windows]

            if len(means) < 2:
                # Too few windows to compute meaningful variance; skip this size
                continue

            var = np.var(means, ddof=1) #ddof means we use unbiased sample variance

            if var <= 0:
                # Zero or negative variance (e.g., perfectly flat coverage); skip to avoid log(0)
                continue

            log_w_vals.append(np.log(w))
            log_var_vals.append(np.log(var))

        if len(log_w_vals) < 2:
            run.warning(f"compute_window_scaling_variance() could not collect enough valid (window_size, variance) "
                        f"pairs to fit a slope (got {len(log_w_vals)}). This may indicate that coverage is "
                        f"extremely flat or that the array is too short.", overwrite_verbose=anvio.DEBUG)
            return None

        slope, _intercept = np.polyfit(log_w_vals, log_var_vals, 1)
        return slope
        
    def midpoint_spread_and_evenness(self, regions):
        """Returns the range and evenness based on the midpoints of coverage regions.
        
        Given a list of region tuples, determines the midpoint of each region, normalizes by 
        contig length (so midpoints become fractions), and computes:
            spread = fraction of contig between first and last coverage region
            evenness = 1 / (1 + CV of midpoint distances)
        """

        starts = np.sort(np.array([r[0] for r in regions if r[2] > 0]))
        stops = np.sort(np.array([r[1] for r in regions if r[2] > 0]))
        midpoints = (starts + stops) / 2
        contig_length = regions[-1][1] # ending index of the last region
        normalized_midpoints = midpoints / contig_length

        midpoint_range = normalized_midpoints.max() - normalized_midpoints.min()
        inter_region_distances = np.diff(normalized_midpoints)
        distance_CV = self.compute_CV(inter_region_distances)
        midpoint_evenness = None
        if distance_CV is not None and len(normalized_midpoints) > 2:
            midpoint_evenness = 1 / (1 + distance_CV)

        return midpoint_range, midpoint_evenness

    def binned_count_dispersion(self, coverage, num_windows=30, min_window_len=100):
        """Returns index of dispersion D = variance(count)/mean(count), where count is the number of 
        bases with coverage in each equal-sized window.

        D should be low (moderate dispersion) for present populations and high (overdispersion) for absent populations.
        """

        window_len = len(coverage) // num_windows
        if window_len < min_window_len:
            return None
        windows = self.get_sliding_window_regions(coverage, window_length=window_len)
        remainder = len(coverage) % num_windows
        if remainder != 0 and remainder < 0.9*window_len:
            windows = windows[:-1] # throw out the last (smaller) region so it doesn't throw off the variance too much

        counts = []
        for w_start, w_stop, w_mean in windows:
            w_bases = coverage[w_start:w_stop]
            w_covered = w_bases[w_bases > 0]
            counts.append(len(w_covered))

        if np.mean(counts) == 0:
            return None
        
        return np.var(counts) / np.mean(counts)

    def Shannon_entropy_evenness(self, coverage):
        """Computes the Shannon entropy of the coverage distribution, normalized to [0, 1].

        Coverage values are treated as a probability distribution by dividing by their sum.
        Entropy is then normalized by log(N), where N is the length of the coverage array,
        which is the maximum possible entropy for a distribution of that length.

        Interpretation:
        ~0.9–1.0: very even coverage, reads are distributed across the genome with little positional bias
        ~0.7–0.9: moderate unevenness, perhaps some peaks or troughs but broadly distributed
        below ~0.5–0.6: substantial coverage concentration, likely indicating strong mapping bias, repetitive regions, or highly variable depth

        Parameters
        ==========
        coverage : array-like
            Array of read recruitment coverage values. Must not sum to zero.

        Returns
        =======
        float
            Evenness score in [0, 1]. Values near 1 indicate coverage is spread
            uniformly across positions; values near 0 indicate coverage is concentrated
            at very few positions.
        """

        coverage = np.asarray(coverage, dtype=float)
        total = coverage.sum()

        if total <= 0:
            raise ConfigError("compute_coverage_evenness() received a coverage array with zero total coverage. "
                              "Shannon entropy is undefined when there are no reads to distribute.")

        N = len(coverage)
        p = coverage / total

        # Compute entropy, ignoring zero-coverage positions (0 * log(0) := 0 by convention)
        nonzero = p > 0
        entropy = -np.sum(p[nonzero] * np.log(p[nonzero]))

        return entropy / np.log(N)


    ##### HELPER FUNCTIONS #####
    def get_list_of_coverage_and_gap_regions(self, coverage, min_gap_length=0):
        """Given an array of coverage values, divides the sequence into gap regions and regions of nonzero coverage.
        
        Each region is described as a tuple of (start_position, stop_position, mean_coverage). If the mean_coverage 
        is 0, then the region is a gap region. Otherwise, it is a covered region.

        Note that start and stop positions follow Python indexing rules to enable slicing of the coverage array: 
        indices start from 0, and the last index of the region is stop_position - 1

        If a gap is too small, the two coverage regions around it are combined into one.
        """
        regions = []

        current_start = 0
        current_stop = 0
        current_state = None
        # create the list of regions
        for i,p in enumerate(coverage):
            current_state = 'covered' if p > 0 else 'gap'
                
            if i == len(coverage) - 1: # EXIT CASE: end of the input sequence
                region_data = (current_start, i+1, np.mean(coverage[current_start:i+1]))
                regions.append(region_data)
                break

            next_state = 'covered' if coverage[i+1] > 0 else 'gap'
            current_stop = i+1
            if next_state != current_state: # STATE TRANSITION: store the previous region and reset for the next one
                region_data = (current_start, current_stop, np.mean(coverage[current_start:current_stop]))
                regions.append(region_data)
                current_start = i+1

        if min_gap_length > 0: # optionally combine coverage regions around small gaps
            small_gap_indices = []
            for i,r in enumerate(regions):
                if i == 0 or i == len(regions) - 1:
                    continue # we skip last and first regions because we can't combine around them
                if r[2] == 0 and (r[1] - r[0] < min_gap_length):
                    small_gap_indices.append(i)
            if small_gap_indices:
                new_regions = []
                skip_next = False
                for i,r in enumerate(regions):
                    if skip_next:
                        skip_next = False
                        continue
                    if i in small_gap_indices:
                        # edit the last entry in the new regions to be the combined coverage region around this small gap
                        prev_cov_start = regions[i-1][0]
                        next_cov_stop = regions[i+1][1]
                        overall_mean = np.mean(coverage[prev_cov_start:next_cov_stop])
                        new_regions[-1] = (prev_cov_start, next_cov_stop, overall_mean)
                        skip_next = True
                    else:
                        new_regions.append(r)

                run.warning(f"Found {len(small_gap_indices)} small gaps (< {min_gap_length} bp) and combined "
                            f"their surrounding coverage regions. Old region list: {regions}. New region "
                            f"list: {new_regions}", overwrite_verbose=anvio.DEBUG)
                regions = new_regions

        return regions

    def get_sliding_window_regions(self, coverage, window_length):
        """Given an array of coverage values, divides it into non-overlapping windows of the requested length.

        Each region is described as a tuple of (start_position, stop_position, mean_coverage), with start and stop 
        positions following Python indexing rules to enable slicing.
        """
        if window_length > len(coverage):
            raise ConfigError(f"get_sliding_window_regions() function was requested to make windows of length {window_length}, "
                              f"but the input coverage array is smaller than that (length {len(coverage)}).")
        windows = []
        current_start = 0
        current_stop = current_start + window_length
        while current_stop < len(coverage):
            region_data = (current_start, current_stop, np.mean(coverage[current_start:current_stop]))
            windows.append(region_data)
            current_start = current_stop
            current_stop = current_start + window_length
            # EDGE CASE: final region is incomplete window
            if current_stop >= len(coverage):
                final_region = (current_start, len(coverage), np.mean(coverage[current_start:]))
                windows.append(final_region)

        return windows
    
    def compute_CV(self, array):
        """Returns coefficient of variation from an input array. Returns None if mean is 0."""

        mean = np.mean(array)
        std = np.std(array)
        if mean == 0:
            return None
        return std / mean

    def compute_gini_coeff(self, array):
        """Computes the Gini coefficient from an input array.

        The equation we use comes from https://en.wikipedia.org/wiki/Gini_coefficient#Alternative_expressions:
        for a list of m values sorted such that g_1 <= g_2 <= .... <= g_m, 
        Gini = (2 x Σᵢ₌₁ᵐ (i x g₍ᵢ₎)) / (m x Σᵢ₌₁ᵐ g₍ᵢ₎) - (m+1)/m
        """

        array = np.sort(array)
        m = array.size
        index = np.arange(1, m + 1)
        
        gini = (2 * np.sum(index * array)) / (m * np.sum(array)) - (m+1)/m
        return gini

    def compute_median_and_MAD(self, array):
        """Returns the median and the median absolute deviation of the input array. If there are less than 3 data points, returns None."""

        if len(array) < 3:
            return None, None
        median = np.median(array)
        mad = np.median(np.abs(array - median))
        return median, mad


    ##### FILTER FUNCTIONS #####
    def filter_nonspecific_regions(self, coverage, region_tuples, mod_z_score_threshold_external, mod_z_score_threshold_internal, min_regions=5, min_detection_for_internal_filter=0.5):
        """Given a set of coverage regions, this function finds and removes outlier coverage values.
        
        If there are many regions, it first selects which of them to investigate further by looking at their mean coverage values. This 
        is the 'external filter': subsetting regions of coverage with suspiciously high mean coverage relative to the other regions.
        Otherwise, if the detection is high enough (ie, there are a few long regions), it examines every single region. 
        
        Selected regions of interest undergo the 'internal filter': identifying outliers in an array of per-nucleotide coverages.
        These outlier regions are simply removed from the original input array of coverages before recomputing the set of coverage 
        vs gap regions in the sequence. Thus, the internal filter returns a new set of coverage regions based on a 'shorter' input
        sequence without outlier bases.

        Parameters
        ==========
        coverage: List of int
            The coverage array of the entire input sequence
        region_tuples : List of tuples
            A list of alternating gap regions and regions of nonzero coverage. Each tuple contains 
            (start_position, stop_position, mean_coverage) for a region. This function works only on
            the regions of nonzero coverage.
        mod_z_score_threshold (internal and external): float
            A value will be flagged as an outlier if its modified z-score is greater than this number.
        min_regions : int
            The minimum number of coverage regions needed for the 'external' filter. If we have fewer
            than this many regions, we instead choose to filter all regions or none depending on the 
            overall detection value.
        min_detection_for_internal_filter : float
            When there are not enough regions for the 'external' filter, we use this detection threshold
            to determine if the regions are likely to be long enough for the 'internal' filter.
        """
        filtered_regions, new_cov_array = None, None
        nonzero_mean_coverages = np.array([m for (x,y,m) in region_tuples if m > 0])
        run.info("EXTERNAL FILTER: number of nonzero coverage regions", len(nonzero_mean_coverages), overwrite_verbose=anvio.DEBUG)
        run.info("EXTERNAL FILTER: minimum number to detect outlier regions", min_regions, overwrite_verbose=anvio.DEBUG)
        run.info("EXTERNAL FILTER will be run", len(nonzero_mean_coverages) >= min_regions, overwrite_verbose=anvio.DEBUG)
        if len(nonzero_mean_coverages) >= min_regions:
            outlier_mean_coverages = utils.get_list_of_outliers(nonzero_mean_coverages, threshold=mod_z_score_threshold_external, only_positive_outliers=True)
            min_outlier_mean_coverage = None # we need this min so we can identify the outlier regions in the original list
            for i, cov in enumerate(nonzero_mean_coverages):
                if outlier_mean_coverages[i]:
                    if not min_outlier_mean_coverage or min_outlier_mean_coverage > cov:
                        min_outlier_mean_coverage = cov 
            run.info("EXTERNAL FILTER: smallest mean coverage of outlier region", min_outlier_mean_coverage, overwrite_verbose=anvio.DEBUG)
            # get the indices of regions with outlier mean coverage in region_tuples
            if not min_outlier_mean_coverage: # there are no outlier regions to filter out
                filtered_regions = region_tuples
                new_cov_array = coverage
                run.info("EXTERNAL FILTER: num outlier regions found", 0, overwrite_verbose=anvio.DEBUG)
            else:
                outlier_region_indices = [i for i,r_tuple in enumerate(region_tuples) if r_tuple[2] >= min_outlier_mean_coverage]
                run.info("EXTERNAL FILTER: num outlier regions found", len(outlier_region_indices), overwrite_verbose=anvio.DEBUG)
                filtered_regions, new_cov_array = self.filter_nonspecific_coverage_internal(coverage, region_tuples, mod_z_score_threshold_internal, 
                                        regions_to_check_indices = outlier_region_indices, drop_entire_region_if_no_internal_outliers=True)
        elif self.detection >= min_detection_for_internal_filter:
            filtered_regions, new_cov_array = self.filter_nonspecific_coverage_internal(coverage, region_tuples, mod_z_score_threshold_internal)
        else:
            run.warning(f"CoverageStats class speaking. While attempting to filter out regions of non-specific "
                        f"read recruitment, we realized that there is not enough data to reliably detect outlier "
                        f"coverage values. Specifically, there were only {len(nonzero_mean_coverages)} regions of "
                        f"nonzero coverage spanning only {self.detection*100:.2f}% of the input sequence. Rather than "
                        f"make any bad calls, we skipped the filtering entirely.",
                        header="FILTERING OF NON-SPECIFIC MAPPING FAILED", overwrite_verbose=anvio.DEBUG)
            filtered_regions = region_tuples
            new_cov_array = coverage

        return filtered_regions, new_cov_array


    def filter_nonspecific_coverage_internal(self, coverage, region_tuples, mod_z_score_threshold, regions_to_check_indices = None, 
                                             drop_entire_region_if_no_internal_outliers=False, min_region_length=20):
        """This function goes through each region of nonzero coverage and removes nucleotides with suspiciously high coverage.

        It creates a subset of the coverage array with the bases of outlier coverage values removed, and returns a new 
        region list created that subset of coverage values. It also returns the subset coverage array.
        
        Parameters
        ==========
        coverage: List of int
            The coverage array of the entire input sequence
        region_tuples : List of tuples
            A list of alternating gap regions and regions of nonzero coverage. Each tuple contains 
            (start_position, stop_position, mean_coverage) for a region. This function works only on
            the regions of nonzero coverage.
        mod_z_score_threshold : float
            A coverage value will be flagged as an outlier if its modified z-score is greater than this number.
        regions_to_check_indices : List of int
            Indices (in region_tuples) of the specific regions to check. If None, check all regions.
        drop_entire_region_if_no_internal_outliers : bool
            If True, it means the parent function flagged these regions as having suspiciously high coverage 
            (ie, due to the external filter) and we should remove the entire region (because if the entire
            region consists entirely of outliers, the internal check won't find them).
        min_region_length : int
            Only detect and remove outliers from regions of nonzero coverage that are at least this long.
        """
        run.info("INTERNAL FILTER: number of regions to examine", len(regions_to_check_indices) if regions_to_check_indices else len(region_tuples), overwrite_verbose=anvio.DEBUG)
        cov_indices_to_remove = []
        for i, r_tuple in enumerate(region_tuples):
            start = r_tuple[0]
            stop = r_tuple[1]
            if (not regions_to_check_indices) or (i in regions_to_check_indices):
                # skip gaps and regions that are too small
                if r_tuple[2] == 0:
                    run.info(f"\tREGION {i}", f"gap region (SKIPPING)", overwrite_verbose=anvio.DEBUG)
                    continue
                elif (stop - start) < min_region_length:
                    run.info(f"\tREGION {i}", f"too small to check (<{min_region_length} bp)", overwrite_verbose=anvio.DEBUG)
                    continue
                region_cov_array = coverage[start:stop]
                outlier_cov_in_region = np.where(utils.get_list_of_outliers(region_cov_array, threshold=mod_z_score_threshold, only_positive_outliers=True) == True)[0]
                # if we didn't find any outliers, the entire region could be outliers. We remove the whole thing when requested.
                if outlier_cov_in_region.shape[0] == 0 and drop_entire_region_if_no_internal_outliers:
                    # we use the region length here because we will add the start position to every element in the array later
                    outlier_cov_in_region = np.arange(0, len(region_cov_array))
                    run.info(f"\tREGION {i}", f"no internal outliers; will remove entire region", overwrite_verbose=anvio.DEBUG)

                if outlier_cov_in_region.shape[0] > 0:
                    outlier_cov_in_original_array = outlier_cov_in_region + start # we add the region start position to each index
                    cov_indices_to_remove.extend(list(outlier_cov_in_original_array))
                    run.info(f"\tREGION {i}", f"removing {outlier_cov_in_region.shape[0]} bases with outlier coverage", overwrite_verbose=anvio.DEBUG)
                else:
                    run.info(f"\tREGION {i}", f"no outliers found", overwrite_verbose=anvio.DEBUG)

        # remove all bases with outlier coverage
        new_coverage = np.delete(coverage, cov_indices_to_remove)
        new_regions = self.get_list_of_coverage_and_gap_regions(new_coverage, min_gap_length=3)
        run.info("INTERNAL FILTER: total bases removed", len(cov_indices_to_remove), overwrite_verbose=anvio.DEBUG)

        return new_regions, new_coverage
