# coding: utf-8
# pylint: disable=line-too-long
"""A class to compute distribution of coverage by evaluating various metrics for 'evenness' of breadth/depth."""

import numpy as np

import anvio
import anvio.utils as utils
import anvio.filesnpaths as filesnpaths

from anvio.terminal import Run

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
        header = ["contig", "sample", "gap_evenness_gini"]
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
        #FIXME: use sliding windows to detect regions
        regions = self.get_list_of_coverage_and_gap_regions(cov_array)

        # filter out regions with unusually-high coverage that might be non-specific read recruitment
        # these regions and their surrounding gaps get replaced with one longer gap
        filtered_detection = self.detection
        if filter_nonspecific_mapping:
            #FIXME: two thresholds for external vs internal filter
            regions, new_cov_array = self.filter_nonspecific_regions(cov_array, regions, mod_z_score_threshold=3)
            filtered_detection = np.sum(new_cov_array > 0) / len(new_cov_array)
            run.info("Detection of contig (post-filter)", filtered_detection, overwrite_verbose=anvio.DEBUG)
        
        # compute all metrics
        gap_gini = self.gap_evenness_gini(regions, min_num_gaps_for_gini=10)

        # append all metrics to file
        output_list = [self.name, self.sample, 
                        f"{gap_gini:.4}\t" if gap_gini else "NA"
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
        

    ##### HELPER FUNCTIONS #####
    def get_list_of_coverage_and_gap_regions(self, coverage):
        """Given an array of coverage values, divides the sequence into gap regions and regions of nonzero coverage.
        
        Each region is described as a tuple of (start_position, stop_position, mean_coverage). If the mean_coverage 
        is 0, then the region is a gap region. Otherwise, it is a covered region.

        Note that start and stop positions follow Python indexing rules to enable slicing of the coverage array: 
        indices start from 0, and the last index of the region is stop_position - 1
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
    def filter_nonspecific_regions(self, coverage, region_tuples, mod_z_score_threshold, min_regions=5, min_detection_for_internal_filter=0.5):
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
        mod_z_score_threshold : float
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
            outlier_mean_coverages = utils.get_list_of_outliers(nonzero_mean_coverages, threshold=mod_z_score_threshold, only_positive_outliers=True)
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
                filtered_regions, new_cov_array = self.filter_nonspecific_coverage_internal(coverage, region_tuples, mod_z_score_threshold, 
                                        regions_to_check_indices = outlier_region_indices, drop_entire_region_if_no_internal_outliers=True)
        elif self.detection >= min_detection_for_internal_filter:
            filtered_regions, new_cov_array = self.filter_nonspecific_coverage_internal(coverage, region_tuples, mod_z_score_threshold)
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
        new_regions = self.get_list_of_coverage_and_gap_regions(new_coverage)
        run.info("INTERNAL FILTER: total bases removed", len(cov_indices_to_remove), overwrite_verbose=anvio.DEBUG)

        return new_regions, new_coverage
