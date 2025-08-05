import os
import yaml
import tarfile
import hashlib

import numpy as np

import anvio
import anvio.constants as constants
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.terminal import Run, Progress, pluralize
from anvio.dbinfo import is_profile_db_and_contigs_db_compatible

from anvio.utils.validation import is_ascii_only
from anvio.utils.validation import check_sample_id


def gzip_compress_file(input_file_path, output_file_path=None, keep_original=False):
    filesnpaths.is_file_exists(input_file_path)

    if not output_file_path:
        output_file_path = input_file_path + '.gz'

    filesnpaths.is_output_file_writable(output_file_path)

    import gzip
    with open(input_file_path, 'rb') as f_in, gzip.open(output_file_path, 'wb') as f_out:
        f_out.writelines(f_in)

    if not keep_original:
        os.remove(input_file_path)



def gzip_decompress_file(input_file_path, output_file_path=None, keep_original=True):
    filesnpaths.is_file_exists(input_file_path)

    if not input_file_path.endswith('.gz'):
        raise ConfigError("gzip_decompress_file function is upset because your input file ('%s') does not "
                          "end with a '.gz' extension :(")

    if not output_file_path:
        output_file_path = input_file_path[:-3]

    filesnpaths.is_output_file_writable(output_file_path)

    import gzip
    with gzip.open(input_file_path, 'rb') as f_in, open(output_file_path, 'wb') as f_out:
        f_out.writelines(f_in)

    if not keep_original:
        os.remove(input_file_path)

    return output_file_path


def tar_extract_file(input_file_path, output_file_path=None, keep_original=True):
    filesnpaths.is_file_tar_file(input_file_path)

    if not output_file_path:
        raise ConfigError("The tar_extract_file function is displeased because an output file path has not been specified. "
                          "If you are seeing this message, you are probably a developer, so go fix your code please, and "
                          "everyone will be happy then.")

    tf = tarfile.open(input_file_path)
    tf.extractall(path = output_file_path)

    if not keep_original:
        os.remove(input_file_path)



def store_array_as_TAB_delimited_file(a, output_path, header, exclude_columns=[]):
    filesnpaths.is_output_file_writable(output_path)

    num_fields = len(a[0])

    if len(header) != num_fields:
        raise ConfigError("store array: header length (%d) differs from data (%d)..." % (len(header), num_fields))

    for col in exclude_columns:
        if not col in header:
            raise ConfigError("store array: column %s is not in the header array...")

    exclude_indices = set([header.index(c) for c in exclude_columns])

    header = [header[i] for i in range(0, len(header)) if i not in exclude_indices]

    f = open(output_path, 'w')
    f.write('%s\n' % '\t'.join(header))

    for row in a:
        f.write('\t'.join([str(row[i]) for i in range(0, num_fields) if i not in exclude_indices]) + '\n')

    f.close()
    return output_path



def store_dataframe_as_TAB_delimited_file(d, output_path, columns=None, include_index=False, index_label="index", naughty_characters=[-np.inf, np.inf], rep_str=""):
    """ Stores a pandas DataFrame as a tab-delimited file.

    Parameters
    ==========
    d: pandas DataFrame
        DataFrame you want to save.
    output_path: string
        Output_path for the file. Checks if file is writable.
    columns: list, pandas.Index, tuple (default = d.columns)
        Columns in DataFrame to write. Default is all, in the order they appear.
    include_index: Boolean (default = False)
        Should the index be included as the first column? Default is no.
    index_label: String (default = "index")
        If include_index is True, this is the header for the index.
    naughty_characters: list (default = [np.inf, -np.inf])
        A list of elements that are replaced with rep_str. Note that all np.nan's (aka NaN's) are also replaced with
        rep_str.
    rep_str: String (default = "")
        The string that elements belonging to naughty_characters are replaced by.

    Returns
    =======
    output_path
    """

    filesnpaths.is_output_file_writable(output_path)

    if not columns:
        columns = d.columns

    d.replace(naughty_characters, np.nan, inplace=True)

    d.to_csv(output_path, sep="\t", columns=columns, index=include_index, index_label=index_label, na_rep=rep_str)
    return output_path



def store_dict_as_TAB_delimited_file(d, output_path, headers=None, file_obj=None, key_header=None, keys_order=None,
                                     header_item_conversion_dict=None, do_not_close_file_obj=False, do_not_write_key_column=False,
                                     none_value=''):
    """Store a dictionary of dictionaries as a TAB-delimited file.

    Parameters
    ==========
    d: dictionary
        A dictionary of dictionaries where each first order key represents a row,
        and each key in the subdictionary represents a column.
    output_path: string
        Output path for the TAB delmited file.path
    headers: list
        Headers of the subdictionary to include (by default include all)
        these are the columns that will be included in the output file (this
        doesn't include the first column which is the keys of the major dictionary)
    file_obj: file_object
        A file object ot write (instead of the output file path)
    key_header: string
        The header for the first column ('key' if None)
    keys_order: list
        The order in which to write the rows (if None, first order keys will be sorted to get the row order)
    header_item_conversion_dict: dictionary
        To replace the column names at the time of writing.
    do_not_close_file_obj: boolean
        If True, file object will not be closed after writing the dictionary to the file
    do_not_write_key_column: boolean
        If True, the first column (keys of the dictionary) will not be written to the file. For use in
        instances when the key is meaningless or arbitrary.
    none_value : string
        What value to write for entries that are None. Default is empty string ('').

    Returns
    =======
    output_path
    """

    if not file_obj:
        filesnpaths.is_output_file_writable(output_path)

    if not file_obj:
        f = open(output_path, 'w')
    else:
        f = file_obj

    key_header = key_header if key_header else 'key'
    if not headers:
        headers = [key_header] + sorted(list(d.values())[0].keys())

    # write header after converting column names (if necessary)
    if header_item_conversion_dict:
        missing_headers = [h for h in headers[1:] if h not in header_item_conversion_dict]
        if len(missing_headers):
            raise ConfigError("Your header item conversion dict is missing keys for one or "
                              "more headers :/ Here is a list of those that do not have any "
                              "entry in the dictionary you sent: '%s'." % (', '.join(missing_headers)))
        if do_not_write_key_column:
            header_text = '\t'.join([header_item_conversion_dict[h] for h in headers[1:]])
        else:
            header_text = '\t'.join([headers[0]] + [header_item_conversion_dict[h] for h in headers[1:]])
    else:
        if do_not_write_key_column:
            header_text = '\t'.join(headers[1:])
        else:
            header_text = '\t'.join(headers)

    if anvio.AS_MARKDOWN:
        tab = '\t'
        f.write(f"|{header_text.replace(tab, '|')}|\n")
        f.write(f"|{':--|' + '|'.join([':--:'] * (len(headers[1:])))}|\n")
    else:
        f.write(f"{header_text}\n")

    if not keys_order:
        keys_order = sorted(d.keys())
    else:
        missing_keys = [k for k in keys_order if k not in d]
        if len(missing_keys):
            if anvio.DEBUG:
                if len(missing_keys) > 10:
                    raise ConfigError("Some keys (n=%d) are not in your dictionary :/ Here is the first ten "
                                      " of them: %s" % (len(missing_keys), missing_keys[:10].__str__()))
                else:
                    raise ConfigError("Some keys are not in your dictionary :/ Here they are: %s" % missing_keys.__str__())
            else:
                raise ConfigError("Some keys are not in your dictionary :/ Use `--debug` to see where this "
                                  "error is coming from the codebase with a list of example keys that are "
                                  "missing.")

    for k in keys_order:
        if do_not_write_key_column:
            line = []
        else:
            line = [str(k)]
        for header in headers[1:]:
            try:
                val = d[k][header]
            except KeyError:
                raise ConfigError("Header ('%s') is not found in the dict :/" % (header))
            except TypeError:
                raise ConfigError("Your dictionary is not properly formatted to be exported "
                                   "as a TAB-delimited file :/ You ask for '%s', but it is not "
                                   "even a key in the dictionary" % (header))

            line.append(str(val) if not isinstance(val, type(None)) else none_value)

        if anvio.AS_MARKDOWN:
            f.write(f"|{'|'.join(map(str, line))}|\n")
        else:
            f.write('%s\n' % '\t'.join(line))

    if not do_not_close_file_obj:
        f.close()

    return output_path



def get_column_data_from_TAB_delim_file(input_file_path, column_indices=[], expected_number_of_fields=None, separator='\t'):
    """Returns a dictionary where keys are the column indices, and items are the list of entries
    found in that that column"""
    filesnpaths.is_file_exists(input_file_path)
    filesnpaths.is_file_tab_delimited(input_file_path, expected_number_of_fields=expected_number_of_fields)

    d = {}

    for index in column_indices:
        d[index] = []

    with open(input_file_path, "r") as input_file:
        for line in input_file.readlines():
            fields = line.strip('\n').split(separator)

            for index in column_indices:
                try:
                    d[index].append(fields[index])
                except:
                    raise ConfigError("get_column_data_from_TAB_delim_file is speaking: The file you sent "
                                       "does not have data for the column index %d. Something is wrong :/" % (index))

    return d



def get_columns_of_TAB_delim_file(file_path, include_first_column=False):
    filesnpaths.is_file_exists(file_path)

    if include_first_column:
        return open(file_path, 'r').readline().strip('\n').split('\t')
    else:
        return open(file_path, 'r').readline().strip('\n').split('\t')[1:]



def is_all_columns_present_in_TAB_delim_file(columns, file_path, including_first_column=False):
    columns_in_file = get_columns_of_TAB_delim_file(file_path, include_first_column=including_first_column)
    return False if len([False for c in columns if c not in columns_in_file]) else True



def transpose_tab_delimited_file(input_file_path, output_file_path, remove_after=False):
    filesnpaths.is_file_tab_delimited(input_file_path)
    filesnpaths.is_output_file_writable(output_file_path)

    file_content = [line.strip('\n').split('\t') for line in open(input_file_path, 'r').readlines()]

    output_file = open(output_file_path, 'w')
    for entry in zip(*file_content):
        output_file.write('\t'.join(entry) + '\n')
    output_file.close()

    if remove_after:
        os.remove(input_file_path)

    return output_file_path



def get_vectors_from_TAB_delim_matrix(file_path, cols_to_return=None, rows_to_return=[], transpose=False, pad_with_zeros=False, run=Run()):
    filesnpaths.is_file_exists(file_path)
    filesnpaths.is_file_tab_delimited(file_path)

    run.warning("Anvi'o is recovering your data from your TAB-delimited file, and it is"
                "instructed to pad your input vectors with 0 values probably because you "
                "used the flag `--pad-input-with-zeros` somewhere. Just so you know.")

    if cols_to_return and pad_with_zeros:
        raise ConfigError("Dear developer, you can't use `cols_to_return` and `pad_with_zeros` at the same "
                          "time with this function. The `pad_with_zeros` header variable in this function "
                          "is a mystery at this point. But the only way to be able to use it requires one "
                          "to not use `cols_to_return`. More mystery .. but essentially this is a necessity "
                          "because we have to update fields_of_interest value if pad_with_zeros is true, so "
                          "anvi'o clustering step DOES NOT IGNORE THE LAST SAMPLE IN THE MATRIX BECUASE PAD "
                          "WITH ZEROS SHIFT EVERYTHING, and we can't do it blindly if the programmer requests "
                          "only specific columnts to be returned with `cols_to_return` :/")

    if transpose:
        transposed_file_path = filesnpaths.get_temp_file_path()
        transpose_tab_delimited_file(file_path, transposed_file_path)
        file_path = transposed_file_path

    rows_to_return = set(rows_to_return)
    vectors = []
    id_to_sample_dict = {}
    sample_to_id_dict = {}

    input_matrix = open(file_path, 'r')
    columns = input_matrix.readline().strip('\n').split('\t')[1:]

    fields_of_interest = []
    if cols_to_return:
        fields_of_interest = [columns.index(col) for col in cols_to_return]
    else:
        fields_of_interest = [f for f in range(0, len(columns)) if constants.IS_ESSENTIAL_FIELD(columns[f])]

    # update columns:
    columns = [columns[i] for i in fields_of_interest]

    if not len(columns):
        raise ConfigError("Only a subset (%d) of fields were requested by the caller, but none of them was found "
                           "in the matrix (%s) :/" % (len(cols_to_return), file_path))

    id_counter = 0
    for line in input_matrix.readlines():
        row_name = line.strip().split('\t')[0]

        if rows_to_return and row_name not in rows_to_return:
            continue

        id_to_sample_dict[id_counter] = row_name
        fields = line.strip('\n').split('\t')[1:]

        # because stupid stuff. see warning above.
        if pad_with_zeros:
            fields = [0] + fields + [0]
            fields_of_interest = list(range(0, len(fields)))

        try:
            if fields_of_interest:
                vector = [float(fields[i]) if fields[i] != '' else None for i in fields_of_interest]
            else:
                # the code will literally never enter here:
                vector = [float(f) if f != '' else None for f in fields]
        except ValueError:
            raise ConfigError("Matrix should contain only numerical values.")

        vectors.append(vector)

        id_counter += 1

    input_matrix.close()

    if transpose:
        # remove clutter
        os.remove(file_path)

    sample_to_id_dict = dict([(v, k) for k, v in id_to_sample_dict.items()])

    return id_to_sample_dict, sample_to_id_dict, columns, vectors



def get_TAB_delimited_file_as_dictionary(file_path, expected_fields=None, dict_to_append=None, column_names=None,\
                                        column_mapping=None, indexing_field=0, separator='\t', no_header=False,\
                                        ascii_only=False, only_expected_fields=False, assign_none_for_missing=False,\
                                        none_value=None, empty_header_columns_are_OK=False, return_failed_lines=False,
                                        ignore_duplicated_keys=False):
    """Takes a file path, returns a dictionary.

       - If `return_failed_lines` is True, it the function will not throw an exception, but instead
         return a list of `failed_lines` along with a dictionary of final results.
    """

    if expected_fields and (not isinstance(expected_fields, list) and not isinstance(expected_fields, set)):
        raise ConfigError("'expected_fields' variable must be a list (or a set).")

    if only_expected_fields and not expected_fields:
        raise ConfigError("'only_expected_fields' variable guarantees that there are no more fields present "
                           "in the input file but the ones requested with 'expected_fields' variable. If you "
                           "need to use this flag, you must also be explicit about what fields you expect to "
                           "find in the file.")

    filesnpaths.is_file_plain_text(file_path)
    filesnpaths.is_file_tab_delimited(file_path, separator=separator)

    failed_lines = []
    column_mapping_for_line_failed = None

    f = open(file_path, 'r')

    # learn the number of fields and reset the file:
    num_fields = len(f.readline().strip('\n').split(separator))
    f.seek(0)

    # if there is no file header, make up a columns list:
    if no_header and not column_names:
        column_names = ['column_%05d' % i for i in range(0, num_fields)]

    if column_names:
        columns = column_names

        if num_fields != len(columns):
            raise  ConfigError("Number of column names declared (%d) differs from the number of columns "
                                "found (%d) in the matrix ('%s') :/" % (len(columns), num_fields, file_path))

        # now we set the column names. if the file had its header, we must discard
        # the first line. so here we go:
        if not no_header:
            f.readline()
    else:
        columns = f.readline().strip('\n').split(separator)

    if not empty_header_columns_are_OK and min(map(len, columns)) == 0:
        raise ConfigError("At least one of the column headers in your tab delimited file '%s' "
                          "is empty." % file_path)

    if expected_fields:
        for field in expected_fields:
            if field not in columns:
                raise ConfigError("The file '%s' does not contain the right type of header. It was expected "
                                   "to have these: '%s', however it had these: '%s'" % (file_path,
                                                                                        ', '.join(expected_fields),
                                                                                        ', '.join(columns[1:])))

    if only_expected_fields:
        for field in columns:
            if field not in expected_fields:
                raise ConfigError("There are more fields in the file '%s' than the expected fields :/ "
                                  "Anvi'o is telling you about this because get_TAB_delimited_file_as_dictionary "
                                  "function is called with `only_expected_fields` flag turned on." % (file_path))

    d = {}
    line_counter = 0

    for line in f.readlines():
        if ascii_only:
            if not is_ascii_only(line):
                raise ConfigError("The input file conitans non-ascii characters at line number %d. Those lines "
                                   "either should be removed, or edited." % (line_counter + 2))

        line_fields = [f if f else None for f in line.strip('\n').split(separator)]

        if line_fields and line_fields[0] == None:
            raise ConfigError("The line number %d in '%s' has no data in its first column, and this doesn't "
                              "seem right at all :/" % (line_counter + 1, file_path))

        if column_mapping:
            column_mapping_for_line_failed = False
            updated_line_fields = []
            for i in range(0, len(line_fields)):
                try:
                    if line_fields[i] == None and column_mapping[i] in [float, int]:
                        updated_line_fields.append(column_mapping[i](0))
                    else:
                        updated_line_fields.append(column_mapping[i](line_fields[i]))
                except NameError:
                    if return_failed_lines:
                        failed_lines.append(line_counter + 1)
                        column_mapping_for_line_failed = True
                    else:
                        raise ConfigError("Mapping function '%s' did not work on value '%s'. These functions can be native "
                                          "Python functions, such as 'str', 'int', or 'float', or anonymous functions "
                                          "defined using lambda notation." % (column_mapping[i], line_fields[i]))
                except TypeError:
                    if return_failed_lines:
                        failed_lines.append(line_counter + 1)
                        column_mapping_for_line_failed = True
                    else:
                        raise ConfigError("Mapping function '%s' does not seem to be a proper Python function :/" % column_mapping[i])
                except ValueError:
                    if return_failed_lines:
                        failed_lines.append(line_counter + 1)
                        column_mapping_for_line_failed = True
                    else:
                        raise ConfigError("Mapping function '%s' did not like the value '%s' in column number %d "
                                          "of the input matrix '%s' :/" % (column_mapping[i], line_fields[i], i + 1, file_path))

            line_fields = updated_line_fields

        if column_mapping_for_line_failed:
            continue

        if indexing_field == -1:
            entry_name = 'line__%09d__' % line_counter
        else:
            entry_name = line_fields[indexing_field]

        if entry_name in d and not ignore_duplicated_keys:
            raise ConfigError("The entry name %s appears more than once in the TAB-delimited file '%s'. There may be more "
                              "instances of duplicated entries, but anvi'o stopped in the first instance. If this is expected "
                              "and if you are a programmer, you can turn on the flag `ignore_duplicated_keys`, which will "
                              "ensure that duplicated entries are ignored, and only the last one is represented in the "
                              "resulting dictionary. If this is expected behavior of the input file and yet you are a user "
                              "then feel free to contact us and we will happily take a look at your case and perhaps update "
                              "the anvi'o program. But if you have gotten this error while working with HMMs, do not contact "
                              "us since helping you in that case is beyond us (see https://github.com/merenlab/anvio/issues/1206 "
                              "for details))." % (entry_name, file_path))

        d[entry_name] = {}

        for i in range(0, len(columns)):
            if i == indexing_field:
                continue
            d[entry_name][columns[i]] = line_fields[i]

        line_counter += 1

    # we have the dict, but we will not return it the way it is if its supposed to be appended to an
    # already existing dictionary.
    if dict_to_append:
        # we don't want to through keys in d each time we want to add stuff to 'dict_to_append', so we keep keys we
        # find in the first item in the dict in another variable. this is potentially very dangerous if not every
        # item in 'd' has identical set of keys.
        keys = list(d.values())[0].keys()

        for entry in dict_to_append:
            if entry not in d:
                # so dict to append is missing a key that is in the dict to be appended. if the user did not
                # ask us to add None for these entries via none_for_missing, we are going to make a noise,
                # otherwise we will tolerate it.
                if not assign_none_for_missing:
                    raise ConfigError("Appending entries to the already existing dictionary from file '%s' failed "
                                       "as the entry %s does not appear to be in the file." % (file_path, entry))
                else:
                    for key in keys:
                        dict_to_append[entry][key] = none_value
            else:
                for key in keys:
                    dict_to_append[entry][key] = d[entry][key]

        return dict_to_append

    # this is here for backward compatibility.
    failed_lines = list(set(failed_lines)) # confirming we are not printing multiple instances of the same line
    if return_failed_lines:
        return d, failed_lines

    return d



def get_BLAST_tabular_output_as_dict(tabular_output_path, target_id_parser_func=None, query_id_parser_func=None):
    """Takes a BLAST output, returns a dict where each query appears only once!!

    If there are multiple hits for a given query, the one with lower e-value.
    remains in the dict.

    Notes
    =====
    - Works only for the default "-outfmt 6"
    """

    results_dict = {}

    for line in open(tabular_output_path):
        fields = line.strip().split('\t')
        query_id = fields[0] if not query_id_parser_func else query_id_parser_func(fields[0])
        target_id = fields[1] if not target_id_parser_func else target_id_parser_func(fields[1])
        e_value = float(fields[10])

        if query_id in results_dict:
            if e_value > results_dict[query_id]['evalue']:
                continue

        results_dict[query_id] = {'hit': target_id, 'evalue': e_value}

    return results_dict



def get_bams_and_profiles_txt_as_data(file_path, no_profile_and_bam_column_is_ok=False):
    """bams-and-profiles.txt is an anvi'o artifact with four columns.

    This function will sanity check one, process it, and return data.

    Updates to this function may require changes in the artifact description at
    anvio/docs/artifacts/bams-and-profiles-txt.md

    Parameters
    ==========
    no_profile_and_bam_column_is_ok : bool
        In some specific instances (e.g., as specific as wanting to compute
        inversion activities in new metagenomes for pre-computed inversion sites),
        downstream analyses may not require profile-db files or BAM files. In such,
        cases the programmer can use this parameter to relax the sanity checks
    """

    COLUMN_DATA = lambda x: get_column_data_from_TAB_delim_file(file_path, [columns_found.index(x)])[columns_found.index(x)][1:]

    if not filesnpaths.is_file_tab_delimited(file_path, dont_raise=True):
        raise ConfigError(f"The bams and profiles txt file must be a TAB-delimited flat text file :/ "
                          f"The file you have at '{file_path}' is nothing of that sorts.")

    if no_profile_and_bam_column_is_ok:
        expected_columns = ['name', 'contigs_db_path']
    else:
        expected_columns = ['name', 'contigs_db_path', 'profile_db_path', 'bam_file_path']

    columns_found = get_columns_of_TAB_delim_file(file_path, include_first_column=True)

    if not set(expected_columns).issubset(set(columns_found)):
        raise ConfigError(f"A bams and profiles txt file is supposed to have at least the following "
                          f"{len(expected_columns)} columns: \"{', '.join(expected_columns)}\".")

    has_profile_db_column = 'profile_db_path' in columns_found
    has_bam_file_column = 'bam_file_path' in columns_found

    names = COLUMN_DATA('name')
    if len(set(names)) != len(names):
        raise ConfigError("Every name listed in the `names` column in a bams and profiles txt must be unique :/ "
                          "You have some redundant names in yours.")

    contigs_db_paths = COLUMN_DATA('contigs_db_path')
    if len(set(contigs_db_paths)) != 1:
        raise ConfigError("All single profiles in bams and profiles file must be associated with the same "
                          "contigs database. Meaning, you have to use the same contigs database path for "
                          "every entry. Confusing? Yes. Still a rule? Yes.")

    if not has_profile_db_column:
        pass
    else:
        profile_db_paths = COLUMN_DATA('profile_db_path')
        if len(set(profile_db_paths)) != len(profile_db_paths):
            raise ConfigError("You listed the same profile database more than once in your bams and profiles txt file :/")

        bam_file_paths = COLUMN_DATA('bam_file_path')
        if len(set(bam_file_paths)) != len(bam_file_paths):
            raise ConfigError("You listed the same BAM file more than once in your bams and profiles txt file :/")

    contigs_db_path = contigs_db_paths[0]
    profiles_and_bams = get_TAB_delimited_file_as_dictionary(file_path)
    for sample_name in profiles_and_bams:
        profiles_and_bams[sample_name].pop('contigs_db_path')
        if has_bam_file_column:
            filesnpaths.is_file_bam_file(profiles_and_bams[sample_name]['bam_file_path'])

        if has_profile_db_column:
            is_profile_db_and_contigs_db_compatible(profiles_and_bams[sample_name]['profile_db_path'], contigs_db_path)

    # this file can optionally contain `r1` and `r2` for short reads
    for raw_reads in ['r1', 'r2']:
        if raw_reads in columns_found:
            file_paths = COLUMN_DATA(raw_reads)
            if '' in file_paths:
                raise ConfigError("If you are using r1/r2 columns in your `bams-and-profiles-txt` file, then you "
                                  "must have a valid file path for every single sample. In your current file there "
                                  "are some blank ones. Sorry.")
            missing_files = [f for f in file_paths if not os.path.exists(f)]
            if len(missing_files):
                raise ConfigError(f"Anvi'o could not find some of the {raw_reads.upper()} files listed in your "
                                  f"`bams-and-profiles-txt` at {file_path} where they were supposed to be: "
                                  f"{missing_files}.")

    return contigs_db_path, profiles_and_bams



def get_samples_txt_file_as_dict(file_path, run=Run(), progress=Progress()):
    "Samples txt file is a commonly-used anvi'o artifact to describe FASTQ file paths for input samples"

    filesnpaths.is_file_tab_delimited(file_path)

    columns_found = get_columns_of_TAB_delim_file(file_path, include_first_column=True)

    if columns_found[0] == 'sample':
        expected_columns = ['sample', 'r1', 'r2']
    elif columns_found[0] == 'name':
        expected_columns = ['name', 'r1', 'r2']
    else:
        raise ConfigError("The first column of any samples-txt must be either `sample` or `name` :/")

    possible_columns = expected_columns + ['group']

    extra_columns = set(columns_found).difference(set(possible_columns))

    if not set(expected_columns).issubset(set(columns_found)):
        raise ConfigError(f"A samples txt file is supposed to have at least the columns {', '.join(expected_columns)}.")

    if len(extra_columns):
        run.warning(f"Your samples txt file contains {pluralize('extra column', len(extra_columns))}: "
                    f"{', '.join(extra_columns)} compared to what is expected of a `samples-txt` file, "
                    f"which is absolutely fine. You're reading this message becasue anvi'o wanted to "
                    f"make sure you know that it knows that it is the case. Classic anvi'o virtue "
                    f"signaling.", lc="yellow")

    samples_txt = get_TAB_delimited_file_as_dictionary(file_path)

    samples_with_missing_files = []
    samples_with_identical_r1_r2_files = []
    for sample_name in samples_txt:
        check_sample_id(sample_name)

        r1_sample_paths = samples_txt[sample_name]['r1'].split(',')
        r2_sample_paths = samples_txt[sample_name]['r2'].split(',')

        if len(r1_sample_paths) != len(r2_sample_paths):
            raise ConfigError(f"Uh oh. The sample {sample_name} has a different number of R1 ({len(r1_sample_paths)}) "
                              f"and R2 ({len(r2_sample_paths)}) paths. Anvi'o expects these to be the same, so please "
                              f"fix this in your samples-txt file.")

        for path in r1_sample_paths + r2_sample_paths:
            if not os.path.exists(path):
                samples_with_missing_files.append(sample_name)

        for i in range(len(r1_sample_paths)):
            if r1_sample_paths[i] == r2_sample_paths[i]:
                samples_with_identical_r1_r2_files.append(sample_name)

    if len(samples_with_missing_files):
        raise ConfigError(f"Bad news. Your samples txt contains {pluralize('sample', len(samples_with_missing_files))} "
                          f"({', '.join(samples_with_missing_files)}) with missing files (by which we mean that the "
                          f"r1/r2 paths are there, but the files they point to are not).")

    if len(samples_with_identical_r1_r2_files):
        raise ConfigError(f"Interesting. Your samples txt contains {pluralize('sample', len(samples_with_missing_files))} "
                          f"({', '.join(samples_with_identical_r1_r2_files)}) where r1 and r2 file paths are identical. Not OK.")

    return samples_txt



def get_primers_txt_file_as_dict(file_path, run=Run(), progress=Progress()):
    """Primers-txt is an anvi'o artifact for primer sequencs."""

    filesnpaths.is_file_tab_delimited(file_path)

    columns_found = get_columns_of_TAB_delim_file(file_path, include_first_column=True)

    if 'name' not in columns_found:
        progress.reset()
        raise ConfigError("A primers-txt file should have a column that is called `name` for the primer name.")

    if 'primer_sequence' not in columns_found:
        progress.reset()
        raise ConfigError("A primers-txt file should have a column that is called `primer_sequence` for the primer sequence.")

    if len(columns_found) < 2:
        progress.reset()
        raise ConfigError("A primers-txt file should have at least two columns - one for primer names, and one for primer sequences.")

    item_column = columns_found[0]
    if item_column != 'name':
        progress.reset()
        raise ConfigError("The first column in your primers-txt file does not seem to be `name`. Anvi'o expects the first "
                          "column to have sequence names.")

    primers_txt = get_TAB_delimited_file_as_dictionary(file_path)

    return primers_txt



def get_groups_txt_file_as_dict(file_path, run=Run(), progress=Progress(), include_missing_samples_is_true=False):
    """Groups-txt is an anvi'o artifact associating items with groups. This function extracts this file into a set of dictionaries.

    Note that it only extracts the first column of the file (which will contain the 'item' or 'sample' information and can have any
    header - let's call these the items) and the 'group' column of the file. Then it will return the following:

    Returns
    =======
    item_to_group_dict : dict
        Dictionary in which keys are items and values are groups
    group_to_item_dict : dict
        Dictionary in which keys are groups and values are lists of items in that group
    include_missing_samples_is_true : Boolean
        set this to True if samples not in this file will be included in the downstream analysis as
        an 'UNGROUPED' group, just to let this function know that it should not crash if less
        than 2 group names are in the groups txt file.
    """

    filesnpaths.is_file_tab_delimited(file_path)

    columns_found = get_columns_of_TAB_delim_file(file_path, include_first_column=True)

    if 'group' not in columns_found:
        raise ConfigError("A groups-txt file should have a column that is called `group`.")

    if len(columns_found) < 2:
        raise ConfigError("A groups-txt file should have at least two columns - one for item names, and one for group names.")

    item_column = columns_found[0]
    if item_column == 'group':
        raise ConfigError("The first column in your groups-txt file appears to be called 'group'. Sadly, anvi'o rather rigidly "
                          "expects the first column to have item names, not group names, so you will have to re-format it. Sorry "
                          "for any inconvenience.")

    groups_txt = get_TAB_delimited_file_as_dictionary(file_path)

    group_to_item_dict = {}
    item_to_group_dict = {}
    for item in groups_txt:
        group_name = groups_txt[item]['group']

        if item in item_to_group_dict:
            raise ConfigError(f"Uh oh. The item {item} occurs more than once in your groups-txt file. This could explode things "
                              f"downstream, so we will stop you right there. Please remove all duplicate items from this file. :)")
        item_to_group_dict[item] = group_name

        if not group_name in group_to_item_dict:
            group_to_item_dict[group_name] = []
        group_to_item_dict[group_name].append(item)

    num_groups = len(group_to_item_dict.keys())
    if num_groups < 2:
        if include_missing_samples_is_true and num_groups == 1:
            run.warning("There is only one group in your groups-txt file, but we have been told that samples not in this file "
                             "will be included in a group called 'UNGROUPED', so that means you have 2 groups in total. Everything is "
                             "fine, as far as we know, but if you look at this and think 'Wait, this is very much NOT fine', well, then "
                             "the power is in your hands to fix it.")
        else:
            raise ConfigError("We notice that there is only one group in your groups-txt file. In the current applications that require "
                          "a groups-txt, we expect to have at least two groups, so we think this is an error. If the context you are "
                          "working in should allow for only one group in this file, please feel free to let us know.")

    return item_to_group_dict, group_to_item_dict



def get_yaml_as_dict(file_path):
    """YAML parser"""

    filesnpaths.is_file_plain_text(file_path)

    try:
        return yaml.load(open(file_path), Loader=yaml.FullLoader)
    except Exception as e:
        raise ConfigError(f"Anvi'o run into some trouble when trying to parse the file at "
                          f"{file_path} as a YAML file. It is likely that it is not a properly "
                          f"formatted YAML file and it needs editing, but here is the error "
                          f"message in case it clarifies things: '{e}'.")



def get_chunk(stream, separator, read_size=4096):
    """Read from a file chunk by chunk based on a separator substring

    This utility of this function is to avoid reading in the entire contents of a file all at once.
    Instead, you can read in a chunk, process it, then read in the next chunk, and repeat this until
    the EOF.

    Parameters
    ==========
    stream : _io.TextIOWrapper
        A file handle, e.g. stream = open('<path_to_file>', 'r')

    separator : str
        Each value returned will be the string from the last `separator` to the next `separator`

    read_size : int, 4096
        How big should each read size be? Bigger means faster reading, but higher memory usage. This
        has no effect on what is returned, but can greatly influence speed. Default is 4MB.

    References
    ==========
    https://stackoverflow.com/questions/47927039/reading-a-file-until-a-specific-character-in-python
    """

    contents_buffer = ''
    while True:
        chunk = stream.read(read_size)
        if not chunk:
            yield contents_buffer
            break

        contents_buffer += chunk
        while True:
            try:
                part, contents_buffer = contents_buffer.split(separator, 1)
            except ValueError:
                break
            else:
                yield part



def concatenate_files(dest_file, file_list, remove_concatenated_files=False):
    if not dest_file:
        raise ConfigError("Destination cannot be empty.")

    filesnpaths.is_output_file_writable(dest_file)

    if not len(file_list):
        raise ConfigError("File list cannot be empty.")

    for f in file_list:
        filesnpaths.is_file_exists(f)

    dest_file_obj = open(dest_file, 'w')
    for chunk_path in file_list:
        for line in open(chunk_path, 'r'):
            dest_file_obj.write(line)

    dest_file_obj.close()

    if remove_concatenated_files:
        for f in file_list:
            os.remove(f)

    return dest_file



def get_file_md5(file_path):
    hash_md5 = hashlib.md5()

    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)

    return hash_md5.hexdigest()
