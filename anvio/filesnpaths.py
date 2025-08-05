# -*- coding: utf-8
# pylint: disable=line-too-long
"""File/Path operations"""

import os
import json
import time
import shutil
import tarfile
import tempfile

import anvio
import anvio.fastalib as u
import anvio.constants as constants

from anvio.terminal import Run
from anvio.terminal import Progress
from anvio.terminal import SuppressAllOutput
from anvio.errors import FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


allowed_chars = constants.allowed_chars.replace('.', '').replace('-', '')
is_bad_column_name = lambda col: len([char for char in col if char not in allowed_chars])


def is_gene_clusters_txt(file_path):
    is_file_tab_delimited(file_path)

    header_proper = ['genome_name', 'gene_caller_id', 'gene_cluster_name']

    with open(file_path, 'r') as input_file:
        header = input_file.readline().strip().split('\t')

    if sorted(header_proper) != sorted(header):
        raise FilesNPathsError(f"The file '{file_path}' does not seem to contain the column headers anvi'o "
                               f"expects to find in a proper gene-clusters-txt file :/ Please check the artifact "
                               f"documentation online.")

    return True


def is_proper_newick(newick_data, dont_raise=False, names_with_only_digits_ok=False):

    with SuppressAllOutput():
        from ete3 import Tree
    try:
        tree = Tree(newick_data, format=1)

        seen = set([])
        duplicates = set([])
        for leaf in tree.get_leaves():
            name = leaf.name
            if name in seen:
                duplicates.add(name)
            seen.add(name)

        if len(duplicates):
            raise Exception("Your newick tree contains duplicate leaves, here is a list of them: %s" % ", ".join(duplicates))

    except Exception as e:
        if dont_raise:
            return False
        else:
            raise FilesNPathsError("Your tree doesn't seem to be properly formatted. Here is what ETE had "
                                   "to say about this: '%s'. Pity :/" % e)

    names_with_only_digits = [n.name for n in tree.get_leaves() if n.name.isdigit()]
    if len(names_with_only_digits) and not names_with_only_digits_ok:
        raise FilesNPathsError("Your tree contains names that are composed of only digits (like this one: '%s'). Sadly, anvi'o "
                               "is not happy with such names in newick trees or clustering dendrograms :( Anvi'o developers "
                               "apologize for the inconvenience." % (names_with_only_digits[0]))

    return True


def is_proper_external_gene_calls_file(file_path):
    is_file_tab_delimited(file_path)

    headers_proper = ['gene_callers_id', 'contig', 'start', 'stop', 'direction', 'partial', 'call_type', 'source', 'version', 'aa_sequence']
    call_types_allowed = set(list(constants.gene_call_types.values()))

    with open(file_path, 'r') as input_file:
        headers = input_file.readline().strip().split('\t')

        if len(headers) == 10:
            missing_headers = [h for h in headers_proper if h not in headers]
            has_aa_sequences = True
        elif len(headers) == 9:
            missing_headers = [h for h in headers_proper[:-1] if h not in headers]
            has_aa_sequences = False
        else:
            raise FilesNPathsError("Your external gene calls file does not contain the right number of columns :/ Here is how "
                                   "your header line should look like (the `aa_sequence` is optional): '%s'." % ', '.join(headers_proper))

        if len(missing_headers):
            raise FilesNPathsError("The headers in your external gene calls file looks wrong :/ Here is how "
                                   "your header line should look like (the `aa_sequence` is optional): '%s'." % ', '.join(headers_proper))

        for i,h in enumerate(headers_proper):
            if not has_aa_sequences and i == 9:
                break
            if headers[i] != h:
                raise FilesNPathsError(f"The headers in your external gene calls file are out of order, so we can't associate each line's fields " 
                                       f"to the right data type. Please re-order the columns to match this order: \"{', '.join(headers_proper)}. "
                                       f"Anvi'o is sorry to make you jump through these hoops, but promises that it is the best way for more "
                                       f"efficient processing of your data.")

        while 1:
            line = input_file.readline()
            if not line:
                break

            fields = line.strip('\n').split('\t')

            try:
                start, stop = int(fields[2]), int(fields[3])
            except ValueError:
                raise FilesNPathsError("All start/stop positions in an external gene calls file must contain integer values (duh). "
                                       "Guess whose file has gene calls with start/stop positions of nope?")

            if start < 0:
                raise FilesNPathsError("At least one gene call in your external genes file ('%s') contains a start position "
                                       "smaller than 0. Anvi'o could extend your contigs with imaginary nucleotides "
                                       "to make things work. Admittedly it would have been much more fun to do that "
                                       "instead of asking you to go back and correct your external gene calls file. "
                                       "But we are burdened to act as adults here :(" % (fields[0]))

            if start >= stop:
                raise FilesNPathsError("At least one gene call in your external genes calls file ('%s') has a stop "
                                       "position that is not larger than the start position. No, says anvi'o. "
                                       "If you need to reverse your genes, the way to do it is to use the `direction`"
                                       "column as it is instructed on our web resources." % (fields[0]))
            try:
                call_type = int(fields[6])
            except ValueError:
                raise FilesNPathsError("Values in the call_type column must be integers :/ Please see "
                                       "http://anvio.org/help/main/artifacts/external-gene-calls/")

            if call_type not in call_types_allowed:
                raise FilesNPathsError("Each call type in an external gene calls file must have a value of either "
                                       "of these: '%s'." % (', '.join([str(e) for e in sorted(list(call_types_allowed))])))

            if call_type is not constants.gene_call_types["CODING"] and has_aa_sequences and len(fields[9].strip()) > 0:
                raise FilesNPathsError("At least one gene call in your external gene calls file ('%s') has amino acid "
                                       "sequence listed despite the fact that it is not marked as 'coding' (1) in `call_type` "
                                       "column. Not OK." % fields[0])

    return True


def is_file_exists(file_path, dont_raise=False):
    if not file_path:
        raise FilesNPathsError("No input file is declared...")
    if not os.path.exists(os.path.abspath(file_path)):
        if dont_raise:
            return False
        else:
            Progress().reset()
            raise FilesNPathsError("No such file: '%s' :/" % file_path)
    return True


def is_output_file_writable(file_path, ok_if_exists=True):
    if not file_path:
        raise FilesNPathsError("No output file is declared...")
    if os.path.isdir(file_path):
        raise FilesNPathsError(f"The path you have provided for your output file ('{os.path.abspath(file_path)}') "
                               f"already is used by a directory :/")
    if not os.access(os.path.dirname(os.path.abspath(file_path)), os.W_OK):
        raise FilesNPathsError(f"It seems you are not autorhized to create an output file at '{file_path}' (lol).")
    if os.path.exists(file_path) and not os.access(file_path, os.W_OK):
        raise FilesNPathsError(f"You do not have write access to the file at '{file_path}' :/")
    if os.path.exists(file_path) and not ok_if_exists:
        if anvio.FORCE_OVERWRITE:
            try:
                os.remove(file_path)
            except Exception as e:
                raise FilesNPathsError(f"As per your instructions, anvi'o was trying to delete the file at '{file_path}', "
                                       f"yet it failed (typical anvi'o?). Here is the error message from another programmer "
                                       f"in the matrix: {e}.")
        else:
            raise FilesNPathsError(f"The output file '{file_path}' already exists. Generally speaking anvi'o tries to "
                                   f"avoid overwriting stuff. But you can always use the flag `--force-overwrite` "
                                   f"to instruct anvi'o to delete the existing file first.")

    return True


def is_output_dir_writable(dir_path):
    if not dir_path:
        raise FilesNPathsError("No output directory path is declared...")
    if not os.path.isdir(dir_path):
        raise FilesNPathsError("'%s' is not a directory..." % dir_path)
    if not os.access(os.path.abspath(dir_path), os.W_OK):
        raise FilesNPathsError("You do not have permission to generate files in '%s'" % dir_path)
    return True


def is_dir_empty(dir_path):
    if not dir_path:
        raise FilesNPathsError("is_dir_empty: No directory path is declared...")
    if not os.path.isdir(dir_path):
        raise FilesNPathsError("is_dir_empty: '%s' is not a directory..." % dir_path)

    return False if len(os.listdir(dir_path)) else True


def is_file_empty(file_path):
    return os.stat(file_path).st_size == 0


def is_file_tab_delimited(file_path, separator='\t', expected_number_of_fields=None, dont_raise=False):
    is_file_exists(file_path)
    f = open(file_path, 'r')

    try:
        while True:
            line = f.readline().strip(' ')
            if line.startswith('#'):
                continue
            else:
                break
    except UnicodeDecodeError:
        if dont_raise:
            return False
        else:
            raise FilesNPathsError(f"The probability of the file at '{file_path}' to be a TAB-delimited file "
                                   f"is zero (unless, of course the locale is not properly set for your shell).")

    if len(line.split(separator)) == 1 and expected_number_of_fields != 1:
        if dont_raise:
            return False
        else:
            raise FilesNPathsError(f"You (or some code on your behalf) asked anvi'o if the file at '{file_path}' "
                                   f"was a TAB-delimited file. Anvi'o took the very first line in it that did not "
                                   f"start with the character '#' (as in commented out lines), and found zero TAB "
                                   f"in it. This is not how we make TAB-delimited files :(")

    f.seek(0)
    num_fields_set = set([len(line.split(separator)) for line in f.readlines()])
    if len(num_fields_set) != 1:
        if dont_raise:
            return False
        else:
            raise FilesNPathsError("Not all lines in the file '%s' have equal number of fields..." % file_path)

    if expected_number_of_fields:
        num_fields_in_file = list(num_fields_set)[0]
        if num_fields_in_file != expected_number_of_fields:
            if dont_raise:
                return False
            else:
                raise FilesNPathsError("The expected number of columns for '%s' is %d. Yet, it has %d "
                                       "of them :/" % (file_path, expected_number_of_fields, num_fields_in_file))

    f.close()
    return True


def is_file_json_formatted(file_path):
    is_file_exists(file_path)

    try:
        json.load(open(file_path, 'r'))
    except ValueError as e:
        raise FilesNPathsError("File '%s' does not seem to be a properly formatted JSON "
                           "file ('%s', cries the library)." % (file_path, e))

    return True


def is_file_fasta_formatted(file_path):
    is_file_exists(file_path)

    try:
        f = u.SequenceSource(file_path)
    except u.FastaLibError as e:
        raise FilesNPathsError("Someone is not happy with your FASTA file '%s' (this is "
                           "what the lib says: '%s'." % (file_path, e))

    f.close()

    return True


def is_file_plain_text(file_path, dont_raise=False):
    is_file_exists(file_path)

    try:
        open(os.path.abspath(file_path), 'r').read(512)
    except IsADirectoryError:
        if dont_raise:
            return False
        else:
            raise FilesNPathsError("Someone want's to make sure %s is a plain text file, however, it is actually a "
                                   "directory :(" % file_path)
    except UnicodeDecodeError:
        if dont_raise:
            return False
        else:
            raise FilesNPathsError("The file at '%s' does not seem to be plain a text file :/" % file_path)

    return True


def is_file_tar_file(file_path, dont_raise=False):
    is_file_exists(file_path)

    if tarfile.is_tarfile(file_path):
        return True
    else:
        if dont_raise:
            return False
        else:
            raise FilesNPathsError("The file at '%s' does not seem to be a tarfile." % file_path)


def is_file_bam_file(file_path, dont_raise=False, ok_if_not_indexed=False):
    """Checks if a BAM file is a proper BAM file, AND if it is intexed"""

    is_file_exists(file_path)

    from pysam import AlignmentFile

    try:
        bam_file = AlignmentFile(file_path, "rb")
    except Exception as e:
        if dont_raise:
            return False
        else:
            raise FilesNPathsError(f"The BAM file you have there upsets samtools very much: '{e}'.")

    if not ok_if_not_indexed:
        try:
            bam_file.mapped
        except ValueError:
            if dont_raise:
                return False
            else:
                raise FilesNPathsError(f"The BAM file at '{file_path}' does not seem to be indexed (when a BAM file) "
                                       f"is indexed, you usually find a file with the same name that ends with '.bam.bai' "
                                       f"extention in the same directory). You can do it via `samtools`, or using the "
                                       f"anvi'o program 'anvi-init-bam'.")

    return True


def is_program_exists(program):
    """adapted from http://stackoverflow.com/a/377028"""
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return True
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return True

    raise FilesNPathsError("'%s' is not found" % program)


def get_temp_directory_path(just_the_path=False):
    temp_directory_path = tempfile.mkdtemp()

    if just_the_path:
        shutil.rmtree(temp_directory_path)

    return temp_directory_path


def get_temp_file_path(prefix=None, just_the_path=True):
    f = tempfile.NamedTemporaryFile(delete=False, prefix=prefix)
    temp_file_name = f.name
    f.close()

    if just_the_path:
        os.remove(temp_file_name)

    return temp_file_name


def get_num_lines_in_file(file_path):
    if os.stat(file_path).st_size == 0:
        return 0

    def blocks(files, size=65536):
        while True:
            b = files.read(size)
            if not b: break
            yield b

    with open(file_path, "r") as f:
        return sum(bl.count("\n") for bl in blocks(f))


def check_output_directory(output_directory, ok_if_exists=False):
    if not output_directory:
        raise FilesNPathsError("Sorry. You must declare an output directory path.")

    output_directory = os.path.abspath(output_directory)

    if os.path.exists(output_directory) and not ok_if_exists:
        if anvio.FORCE_OVERWRITE:
            try:
                shutil.rmtree(output_directory)
            except Exception as e:
                raise FilesNPathsError(f"As per your instructions, anvi'o was trying to delete the directory at '{output_directory}', "
                                       f"yet it failed (typical anvi'o?). Here is the error message from another programmer in "
                                       f"the matrix: {e}.")
        else:
            raise FilesNPathsError(f"The output directory '{output_directory}' already exists (and anvi'o does not like "
                                   f"overwriting stuff (except when it does (typical anvi'o))). But you can always use "
                                   f"the flag `--force-overwrite` to assert your dominance. In which case anvi'o would "
                                   f"first remove the existing output directory (a flag that deserves extreme caution "
                                   f"for obvious reasons).")

    return output_directory


def gen_output_directory(output_directory, progress=Progress(verbose=False), run=Run(), delete_if_exists=False, dont_warn=False):
    if not output_directory:
        raise FilesNPathsError("Someone called `gen_output_directory` function without an output "
                               "directory name :( An embarrassing moment for everyone involved.")

    if os.path.exists(output_directory) and delete_if_exists and not is_dir_empty(output_directory):
        try:
            if not dont_warn:
                run.warning('The existing directory "%s" is about to be removed... (You have '
                            '20 seconds to press CTRL + C). [filesnpaths::gen_output_directory]' % output_directory,
                             header = '!!! READ THIS NOW !!!')
                time.sleep(20)
            shutil.rmtree(output_directory)
        except:
            progress.end()
            raise FilesNPathsError("I was instructed to remove this directory, but I failed: '%s' :/" % output_directory)

    if not os.path.exists(output_directory):
        try:
            os.makedirs(output_directory)
        except:
            progress.end()
            raise FilesNPathsError("Output directory does not exist (attempt to create one failed as well): '%s'" % \
                                                            (output_directory))
    if not os.access(output_directory, os.W_OK):
        progress.end()
        raise FilesNPathsError("You do not have write permission for the output directory: '%s'" % output_directory)

    return output_directory


def get_name_from_file_path(file_path, postfix_separator="."):
    """Return a decent name for a given file at file at `file_path`"""

    splits = os.path.basename(os.path.abspath(file_path)).split(postfix_separator)

    return splits[0] if len(splits) in [1, 2] else '.'.join(splits[:-1])
