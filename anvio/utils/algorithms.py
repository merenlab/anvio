import gzip

import numpy as np
import pandas as pd
import itertools as it

from numba import jit

from anvio.errors import ConfigError

# for full output
pd.options.display.max_columns=100
pd.options.display.max_rows=100

def convert_numpy_array_to_binary_blob(array, compress=True):
    if compress:
        return gzip.compress(memoryview(array), compresslevel=1)
    else:
        return memoryview(array)



def convert_binary_blob_to_numpy_array(blob, dtype, decompress=True):
    if decompress:
        return np.frombuffer(gzip.decompress(blob), dtype=dtype)
    else:
        return np.frombuffer(blob, dtype=dtype)



@jit(nopython=True)
def add_to_2D_numeric_array(x, y, a, count=1):
    """just-in-time compiled function

    Parameters
    ==========
    x : array
        array of row indices
    y : array
        array of corresponding y indices
    count : int, 1
        How much to add to each coordinate

    Examples
    ========

    Make a 5x20000 array (a) and define 95 coordinate positions to update (i and p)

    >>> a = np.zeros((5, 20000))
    >>> i = np.random.choice(range(5), size=95, replace=True)
    >>> p = np.random.choice(range(100), size=95, replace=False) + 1000

    For comparison, define the slow method

    >>> def add_to_2D_numeric_array_slow(x, y, a, count=1):
    >>>     for idx, pos in zip(x, y):
    >>>         a[idx, pos] += count
    >>>     return a

    Compare the speeds

    >>> %timeit add_to_2D_numeric_array_slow(i, p, a)
    74.5 µs ± 4.42 µs per loop (mean ± std. dev. of 7 runs, 10000 loops each)
    >>> %timeit _add_to_2D_numeric_array(i, p, a)
    798 ns ± 12.7 ns per loop (mean ± std. dev. of 7 runs, 1000000 loops each)
    """
    for idx, pos in zip(x, y):
        a[idx, pos] += count

    return a



def apply_and_concat(df, fields, func, column_names, func_args=tuple([])):
    """ This function has been taken from https://tinyurl.com/y9ylqy4l
        and has been modified for speed considerations using this blog post:
        https://tinyurl.com/ya4e5tz3. Its utility is to append multiple columns to an existing
        dataframe row by row. This is usually a bad idea because usually operations can be
        vectorized. However when they cannot, looping through each row becomes a necessary evil.

        df: pandas DataFrame object
            An existing dataframe to loop through append columns to.
        fields: list
            A list of columns in the existing dataframe used to calculate the new columns
        func: function
            A function that takes as its first argument a row of `df` (i.e. a pd.Series
            object) and potential additional positional arguments `func_args`. It should return a
            tuple of values with with the same length as `column_names`.
        func_args: tuple
            A tuple of arguments passed to `func` besides the assumed first argument (a pd.Series
            object). For example, is `def func(row, a)`, then `func_args = (a,)`. If func_args is an
            empty tuple, `func` should take no other args.
        column_names: list
            A list of column headers for the newly appended columns
    """
    d = {column_name: [] for column_name in column_names}
    for _, row in df[fields].iterrows():
        out_values = func(row, *func_args)
        for ind, column_name in enumerate(column_names):
            d[column_name].append(out_values[ind])

    df2 = pd.DataFrame(d, index=df.index)
    return pd.concat((df, df2), axis=1, sort=True)



def get_indices_for_outlier_values(c):
    is_outlier = get_list_of_outliers(c)
    return set([p for p in range(0, c.size) if is_outlier[p]])



def get_list_of_outliers(values, threshold=None, zeros_are_outliers=False, median=None):
    """Return boolean array of whether values are outliers (True means outlier)

    Modified from Joe Kington's (https://stackoverflow.com/users/325565/joe-kington)
    implementation computing absolute deviation around the median.

    Parameters
    ==========
    values : array-like
        A num_observations by num_dimensions array of observations.

    threshold : number, None
        The modified z-score to use as a thresholdold. Observations with
        a modified z-score (based on the median absolute deviation) greater
        than this value will be classified as outliers.

    median : array-like, None
        Pass median of values if you already calculated it to save time.

    Returns
    =======
    mask : numpy array (dtype=bool)
        A num_observations-length boolean array. True means outlier

    Examples
    ========

    Create an array with 5 manually created outliers:

    >>> import numpy as np
    >>> import anvio.utils as utils
    >>> array = 10*np.ones(30) + np.random.rand(30)
    >>> array[5] = -10
    >>> array[9] = -10
    >>> array[12] = -10
    >>> array[15] = -10
    >>> array[23] = -10
    >>> mask = utils.get_list_of_outliers(array, threshold=5)
    >>> print(mask)
    [False False False False False  True False False False  True False False
     True False False  True False False False False False False False  True
     False False False False False False]

    As can be seen, mask returns a numpy array of True/False values, where True corresponds to
    outlier values.

    >>> print(outlier_indices)
    >>> outlier_indices = np.where(mask == True)[0]
    [ 5  9 12 15 23]

    The `True` values indeed occur at the indices where the values hand been manually changed to
    represent outliers.

    References
    ==========
    Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
    Handle Outliers", The ASQC Basic References in Quality Control:
    Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.

    http://www.sciencedirect.com/science/article/pii/S0022103113000668
    """

    if threshold is None:
        threshold = 1.5

    if len(values.shape) == 1:
        values = values[:, None]

    if not median: median = np.median(values, axis=0)

    diff = np.sum((values - median) ** 2, axis=-1)
    diff = np.sqrt(diff)
    median_absolute_deviation = np.median(diff)

    if not median_absolute_deviation:
       if values[0] == 0:
            # A vector of all zeros is considered "all outliers"
            return np.array([True] * values.size)
       else:
            # A vector of uniform non-zero values is "all non-outliers"
            # This could be important for silly cases (like in megahit) in which there is a maximum value for coverage
            return np.array([False] * values.size)

    modified_z_score = 0.6745 * diff / median_absolute_deviation
    non_outliers = modified_z_score > threshold

    if not zeros_are_outliers:
        return non_outliers
    else:
        zero_positions = [x for x in range(len(values)) if values[x] == 0]
        for i in zero_positions:
            non_outliers[i] = True
        return non_outliers



def get_stretches_for_numbers_list(numbers_list, discard_singletons=False):
    """Takes a array of numbers, and turns reports back stretches

    For example, for a `numbers_list` that looks like this:

        >>> [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 15, 16, 17, 18, 99]

    This function will return the following:

        >>> [(2, 12), (15, 18), (99, 99)]

    If the user chooses to `discard_singletons`, then for the same input list,
    they will get the following:

        >>> [(2, 12), (15, 18), (99, 99)]

    The output of this function can be the input of `merge_stretches` function.
    """

    if not len(numbers_list):
        return []

    numbers_list = sorted(numbers_list)

    stretches = []
    groups = it.groupby(numbers_list, key=lambda item, c=it.count():item-next(c))

    for k, g in groups:
        stretch = list(g)

        if discard_singletons and len(stretch) == 1:
            continue

        stretches.append((stretch[0], stretch[-1]), )

    return stretches



def merge_stretches(stretches, min_distance_between_independent_stretches=0):
    """A function to merge stretches of indices in an array.

    It takes an array, `stretches`, that looks like this:

        >>> [(3, 9), (14, 27), (32, 36), (38, 42)]

    And returns an array like this, if `min_distance_between_independent_stretches`, say, 3:

        >>> [(3, 9), (14, 27), (32, 42)]

    If all you have is an array of numbers rather than stretches to be merbed, see
    the function `get_stretches_for_numbers_list`.
    """

    if not len(stretches):
        return None

    if not isinstance(stretches[0], tuple):
        raise ConfigError("This function expect a list of tuples :/")

    STRETCHES_FINAL = stretches
    STRETCHES_CURRENT = stretches

    while 1:
        stretches_to_merge = []

        CURRENT = 0
        START, END = 0, 1
        while 1:
            if not len(STRETCHES_CURRENT):
                break

            NEXT = CURRENT + 1

            if NEXT == len(STRETCHES_CURRENT):
                stretches_to_merge.append([STRETCHES_CURRENT[CURRENT]])
                break

            while 1:
                if NEXT > len(STRETCHES_CURRENT):
                    break

                if STRETCHES_CURRENT[NEXT][START] - STRETCHES_CURRENT[CURRENT][END] < min_distance_between_independent_stretches:
                    NEXT = NEXT + 1

                    if NEXT == len(STRETCHES_CURRENT):
                        break
                else:
                    break

            if NEXT > len(STRETCHES_CURRENT):
                break
            elif NEXT - CURRENT == 1:
                stretches_to_merge.append([STRETCHES_CURRENT[CURRENT]])
                CURRENT += 1
            else:
                stretches_to_merge.append(STRETCHES_CURRENT[CURRENT:NEXT])
                CURRENT = NEXT

        # here the array `stretches_to_merge` contains all the lists of
        # stretches that need to be merged.
        STRETCHES_MERGED = [(s[0][0], s[-1][1]) for s in stretches_to_merge]

        if STRETCHES_FINAL == STRETCHES_MERGED:
            break
        else:
            STRETCHES_FINAL = STRETCHES_MERGED
            STRETCHES_CURRENT = STRETCHES_MERGED

    return STRETCHES_FINAL



@jit(nopython=True)
def get_constant_value_blocks(array, value):
    """Generator that returns blocks of consecutive numbers

    Parameters
    ==========
    array : array
        a numerical numpy array. If a list is passed, this function is very slow

    value : number
        The number you want to get constant blocks for.

    Examples
    ========

    >>> a = np.array([47, 47, 47, 49, 50, 47, 47, 99])
    >>> for i in get_constant_value_blocks(a, 47): print(i)
    (0, 3)
    (5, 7)
    """
    ans = []

    matching = False
    for i in range(len(array)):
        if array[i] == value:
            if not matching:
                start = i
                matching = True
        else:
            if matching:
                matching = False
                ans.append((start, i))

    if matching:
        ans.append((start, i + 1))

    return ans



@jit(nopython=True)
def find_value_index(x, val, reverse_search=False):
    """Returns first instance of indices where a value is found

    Created this because unlike np.where, this stops after the first instance is found. If you only
    want the first instance, this algorithm is therefore preferrable for array sizes < 1000 (see
    examples)

    Parameters
    ==========
    x : 1D array

    val : number
        return index of x where the value == val.

    reverse_search : bool, False
        Search in reverse order

    Examples
    ========
    >>> import numpy as np
    >>> import anvio.utils as utils
    >>> x = np.arange(1000)
    >>> %timeit utils.find_value_index(x, 999, rev=True)
    574 ns ± 15.8 ns per loop (mean ± std. dev. of 7 runs, 1000000 loops each)
    >>> %timeit utils.find_value_index(x, 999)
    2.21 µs ± 36.7 ns per loop (mean ± std. dev. of 7 runs, 100000 loops each)
    >>> %timeit np.where(x == 999)[0][0]
    2.91 µs ± 563 ns per loop (mean ± std. dev. of 7 runs, 100000 loops each)
    """

    for i in range(len(x)) if not reverse_search else range(len(x)-1, -1, -1):
        if x[i] == val:
            return i



def get_N50(contig_lengths):
    h, S = sum(contig_lengths) / 2.0, 0

    for l in sorted(contig_lengths, reverse=True):
        S += l
        if h < S:
            return l



def split_by_delim_not_within_parens(d, delims, return_delims=False):
    """Takes a string, and splits it on the given delimiter(s) as long as the delimeter is not within parentheses.

    This function exists because regular expressions don't handle nested parentheses very well.

    The function can also be used to determine if the parentheses in the string are unbalanced (it will return False
    instead of the list of splits in this situation)

    PARAMETERS
    ==========
    d : str
        string to split
    delims : str or list of str
        a single delimiter, or a list of delimiters, to split on
    return_delims : boolean
        if this is true then the list of delimiters found between each split is also returned

    RETURNS
    =======
    If parentheses are unbalanced in the string, this function returns False. Otherwise:
    splits : list
        strings that were split from d
    delim_list : list
        delimiters that were found between each split (only returned if return_delims is True)
    """

    parens_level = 0
    last_split_index = 0
    splits = []
    delim_list = []
    for i in range(len(d)):
        # only split if not within parentheses
        if d[i] in delims and parens_level == 0:
            splits.append(d[last_split_index:i])
            delim_list.append(d[i])
            last_split_index = i + 1 # we add 1 here to skip the space
        elif d[i] == "(":
            parens_level += 1
        elif d[i] == ")":
            parens_level -= 1

        # if parentheses become unbalanced, return False to indicate this
        if parens_level < 0:
            return False
    splits.append(d[last_split_index:len(d)])

    if return_delims:
        return splits, delim_list
    return splits



def get_split_start_stops(contig_length, split_length, gene_start_stops=None):
    """Wrapper function for get_split_start_stops_with_gene_calls and get_split_start_stops_without_gene_calls"""
    if gene_start_stops:
        return get_split_start_stops_with_gene_calls(contig_length, split_length, gene_start_stops)
    else:
        return get_split_start_stops_without_gene_calls(contig_length, split_length)



def get_split_start_stops_with_gene_calls(contig_length, split_length, gene_start_stops):
    """Here we have a contig of `contig_length`, and a desired split length of `split_length`. also
       we know where genes start and stop in this contigs. we would like to split this contig into
       smaller pieces, i.e. sizes of `splits_length`, but in such a way that that splits do not
       break contigs in the middle of a gene."""

    # if the contig is too short, return it back.
    if contig_length < 2 * split_length:
        return [(0, contig_length)]

    coding_positions_in_contig = []

    # Pretend the beginning and end are coding (even if they aren't) so that we prevent very short pieces.
    for position in it.chain(range(int(split_length / 2)), range(contig_length - int(split_length / 2), contig_length)):
        coding_positions_in_contig.append(position)

    # Track positions that code for genes.
    for gene_unique_id, start, stop in gene_start_stops:
        start = start - 5
        stop = stop + 5

        for position in range(start, stop):
            coding_positions_in_contig.append(position)

    non_coding_positions_in_contig = set(range(contig_length)) - set(coding_positions_in_contig)

    # what would be our break points in an ideal world? compute an initial list of break
    # points based on the length of the contig and desired split size:
    optimal_number_of_splits = int(contig_length / split_length)

    optimal_split_length = int(contig_length / optimal_number_of_splits)
    optimal_break_points = list(range(optimal_split_length, contig_length - optimal_split_length + 1, optimal_split_length))

    # now we will identify the very bad break points that we can't find a way to split around
    bad_break_points = set([])
    for i in range(0, len(optimal_break_points)):
        break_point = optimal_break_points[i]

        if break_point not in non_coding_positions_in_contig:
            # optimal break point hits a gene. we shall search towards both directions
            # to find a better break point:
            new_break_point = None
            for s in range(0, int(split_length / 2)):
                if break_point + s in non_coding_positions_in_contig:
                    new_break_point = break_point + s
                    break
                if break_point - s in non_coding_positions_in_contig:
                    new_break_point = break_point - s
                    break

            if not new_break_point:
                # nope. we failed. this is a bad bad break point.
                bad_break_points.add(break_point)
            else:
                # we are satisfied with the new one we found for now. it may be a shitty one,
                # but we will learn about that later. for now, let's replace the previous
                # optimal break pont with this one.
                optimal_break_points[i] = new_break_point

    # remove all the bad breakpoints from our 'optimal' break points:
    optimal_break_points = [p for p in optimal_break_points if p not in bad_break_points]

    if not len(optimal_break_points):
        # we have nothing left to work with after removal of the crappy break points. we will
        # keep this bad boy the way it is.
        return [(0, contig_length)]

    # create start/stop positions from these break points
    chunks = list(zip([0] + optimal_break_points[:-1], optimal_break_points)) + [(optimal_break_points[-1], contig_length)]

    return chunks



def get_split_start_stops_without_gene_calls(contig_length, split_length):
    """Returns split start stop locations for a given contig length."""
    num_chunks = int(contig_length / split_length)

    if num_chunks < 2:
        return [(0, contig_length)]

    chunks = []
    for i in range(0, num_chunks):
        chunks.append((i * split_length, (i + 1) * split_length),)
    chunks.append(((i + 1) * split_length, contig_length),)

    if (chunks[-1][1] - chunks[-1][0]) < (split_length / 2):
        # last chunk is too small :/ merge it to the previous one.
        last_tuple = (chunks[-2][0], contig_length)
        chunks.pop()
        chunks.pop()
        chunks.append(last_tuple)

    return chunks


def human_readable_file_size(nbytes):
    suffixes = ['B', 'KB', 'MB', 'GB', 'TB', 'PB']
    if nbytes == 0: return '0 B'
    i = 0
    while nbytes >= 1024 and i < len(suffixes)-1:
        nbytes /= 1024.
        i += 1
    f = ('%.2f' % nbytes).rstrip('0').rstrip('.')
    return '%s %s' % (f, suffixes[i])
