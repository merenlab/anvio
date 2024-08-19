# coding utf-8
"""An interface for Vmatch"""

import os
import time
import queue
import shutil
import pandas as pd

# multiprocess is a fork of multiprocessing that uses the dill serializer instead of pickle
# using the multiprocessing module directly results in a pickling error in Python 3.10 which
# goes like this:
#
#   >>> AttributeError: Can't pickle local object 'SOMEFUNCTION.<locals>.<lambda>' multiprocessing
#
import multiprocess as multiprocessing

from collections import defaultdict

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Samuel Miller"
__email__ = "samuelmiller10@gmail.com"


pp = terminal.pretty_print


class Vmatch(object):
    """Parent class for all Vmatch usage"""

    QUERY_CHUNK_SIZE_DEFAULT = 100000

    def __init__(self, args, run=None, progress=None, skip_sanity_check=False):
        self.index_program_name = 'mkvtree'
        self.search_program_name = 'vmatch'

        self.tested_versions = ('2.3.0', )
        self.supported_match_modes = ('exact_query_substring', 'query_substring_with_mismatches', 'query_substring_with_indels')

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.match_mode = A('match_mode')
        self.fasta_db_path = A('fasta_db_file')
        self.fasta_query_path = A('fasta_query_file')

        # Multithreading is not internal to Vmatch. The query file is chunked and different
        # instances are run. Defining the number of sequences in each chunk prevents unexpectedly
        # massive output files from being generated -- the output from each chunk is processed and
        # deleted before the next chunk runs.
        self.num_threads = A('num_threads') or 1
        # To manage disk space, scale the default chunk size with the number of threads using a rule
        # of thumb. Intermediate files for chunks can build up dramatically though they are
        # regularly deleted.
        self.query_chunk_size = A('query_chunk_size') or self.QUERY_CHUNK_SIZE_DEFAULT // self.num_threads

        self.max_hamming_dist = A('max_hamming_dist')
        self.max_edit_dist = A('max_edit_dist') # Levenshtein distance
        # Reject indels within a certain distance of the ends of the alignment by providing left and
        # right buffers.
        self.edit_left_buffer = A('edit_left_buffer')
        self.edit_right_buffer = A('edit_right_buffer')
        self.min_ident = A('min_ident') or 100

        # With this option, sequence alignments are included in the output. They are wrapped to the
        # given length on each line.
        self.align_output_length = A('align_output_length')

        self.quiet = A('quiet') or True
        self.temp_dir = A('temp_dir') or filesnpaths.get_temp_directory_path()
        self.keep_temp_dir = A('keep_temp_dir') or False
        self.log_path = A('log_path') or os.path.join(self.temp_dir, '00_log.txt')

        self.run = run or terminal.Run(verbose=(not self.quiet))
        self.progress = progress or terminal.Progress(verbose=(not self.quiet))

        if not skip_sanity_check:
            self.sanity_check()


    def sanity_check(self):
        """Do basic checks before proceeding."""
        self.check_programs()

        for name, variable in (("the match mode", self.match_mode),
                               ("the query sequence FASTA file path", self.fasta_db_path),
                               ("the target sequence FASTA file path", self.fasta_db_path),
                               ("the minimum identity", self.min_ident)):
            if not variable:
                raise ConfigError(f"A proper instance of the anvi'o Vmatch driver must have {name} variable set.")

        if self.match_mode not in self.supported_match_modes:
            raise ConfigError(f"The supported Vmatch match modes are: {', '.join(self.supported_match_modes)}. "
                              f"No mode named {self.match_mode} is allowed.")

        if self.match_mode == 'query_substring_with_mismatches':
            if self.max_hamming_dist is None:
                raise ConfigError("The match mode, 'query_substring_with_mismatches', requires the maximum Hamming distance argument.")
            if self.max_hamming_dist % 1 != 0 or self.max_hamming_dist <= 1:
                raise ConfigError("The maximum allowed Hamming distance must be a positive integer.")
        else:
            if self.max_hamming_dist is not None:
                raise ConfigError(f"The match mode, '{self.match_mode}', is incompatible with the maximum Hamming distance argument.")

        if self.match_mode == 'query_substring_with_indels':
            if self.max_edit_dist is None:
                raise ConfigError("The match mode, 'query_substring_with_indels', requires the maximum edit (Levenshtein) distance argument.")
            if self.max_edit_dist % 1 != 0 or self.max_edit_dist <= 1:
                raise ConfigError("The maximum allowed edit (Levenshtein) distance must be a positive integer.")
        else:
            if self.max_edit_dist is not None:
                raise ConfigError(f"The match mode, '{self.match_mode}', is incompatible with the maximum edit (Levenshtein) distance argument.")
            if self.edit_left_buffer is not None:
                raise ConfigError(f"The match mode, '{self.match_mode}', is incompatible with the edit left buffer argument.")
            if self.edit_right_buffer is not None:
                raise ConfigError(f"The match mode, '{self.match_mode}', is incompatible with the edit right buffer argument.")

        if self.match_mode == 'query_substring_with_mismatches' or self.match_mode == 'query_substring_with_indels':
            if not self.align_output_length:
                raise ConfigError(f"The match mode, '{self.match_mode}', requires the argument, 'align_output_length'.")
        else:
            if self.align_output_length is not None:
                raise ConfigError(f"The match mode, '{self.match_mode}', is incompatible with the alignment output length argument.")

        filesnpaths.is_file_exists(self.fasta_db_path)
        self.index_path = os.path.join(self.temp_dir, 'index')

        if self.min_ident % 1 != 0:
            raise ConfigError("The minimum identity of an alignment must be an integer between 1 and 100. "
                              f"The provided value was {self.min_ident}.")


    def check_programs(self):
        self.installed_index_program_version = None
        self.installed_search_program_version = None

        for program_name, attr_name in zip((self.index_program_name, self.search_program_name),
                                           ('installed_index_program_version', 'installed_search_program_version')):
            utils.is_program_exists(program_name)

            output, ret_code = utils.get_command_output_from_shell(f'{program_name} -version')
            try:
                version_found = output.split(b'\n')[0].split(b' ')[2].decode("utf-8")
                if not self.quiet:
                    self.run.info(f"{program_name} version found", version_found, mc='green')
            except:
                version_found = 'Unknown'
                self.run.warning(f"Anvi'o failed to learn the version of {program_name} installed on this system :/")

            if version_found not in self.tested_versions:
                self.run.warning(f"The version of {program_name} installed on your system ('{version_found}') "
                                 f"has been not tested with the anvi'o {program_name} driver. "
                                 "Anvi'o will continue to try to run everything as if this didn't happen. "
                                 "If you see this warning but everything works fine, "
                                 "let us know so we can include this version number into the list of 'tested' version numbers. "
                                 "If you see an unexpected error, please consider installing one of these versions of Vmatch "
                                 "(and again please let us know anyway so we can address it for later): "
                                 f"'{', '.join(list(self.tested_versions))}'")

            setattr(self, attr_name, version_found)

        if self.installed_index_program_version != self.installed_search_program_version:
            self.run.warning("For some hopefully good reason, the Vmatch commands, "
                             f"'{self.index_program_name}' and '{self.search_program_name}', are not the same on your system. "
                             f"The version of {self.index_program_name} that was found is {self.installed_index_program_version}. "
                             f"The version of {self.search_program_name} that was found is {self.installed_search_program_version}. "
                             "If you see an unexpected error, please consider installing one of these versions of Vmatch: "
                             f"'{', '.join(list(self.tested_versions))}'")


    def make_index(self):
        self.progress.new("Vmatch vmktree")
        self.progress.update(f"Creating a persistent index from {self.fasta_db_path}")
        if self.match_mode == 'exact_query_substring':
            # `-pl` means that sequence suffixes are placed into buckets by their prefix length,
            # here automatically determined.
            # `-tis -suf -sti1 -lcp -bck` are necessary index tables.
            # `-v` produces verbose commented output.
            command = ['mkvtree',
                       '-db', f'{self.fasta_db_path}',
                       '-dna',
                       '-indexname', f'{self.index_path}',
                       '-pl',
                       '-tis', '-suf', '-sti1', '-lcp', '-bck',
                       '-v']
        else:
            # `-ois` produces an index table of original sequences and is needed to output sequence alignments.
            # `-skp` produces an index table used in finding whole queries with mismatches.
            command = ['mkvtree',
                       '-db', f'{self.fasta_db_path}',
                       '-dna',
                       '-indexname', f'{self.index_path}',
                       '-pl',
                       '-tis', '-suf', '-sti1', '-lcp', '-bck', '-ois', '-skp',
                       '-v']
        ret_val = utils.run_command(command, self.log_path, remove_log_file_if_exists=False)

        if int(ret_val):
            raise ConfigError("vmktree returned with non-zero exit code, there may be some errors. "
                              f"Please check the log file for details: {self.log_path}")
        self.progress.end()


    def search_queries(self, make_index=True):
        if make_index:
            self.make_index()

        pid = "Vmatch"
        self.progress.new(pid)
        self.progress.update(f"Setting up search")

        with open(self.fasta_query_path) as query_file:
            for line_num, line in enumerate(query_file, 1):
                pass
        total_lines = line_num
        pp_total_lines = pp(total_lines)
        query_chunk_size = self.query_chunk_size * 2 # convert sequences to FASTA lines
        chunk_start = 0
        chunk_stop = query_chunk_size
        unprocessed_chunk_num = 1
        unprocessed_chunk_dict = defaultdict(dict)
        running_line_count = 0

        match_mode = self.match_mode
        if chunk_stop < total_lines:
            query_chunked = True
        else:
            query_chunked = False
            # All query sequences fit in one chunk, so use the original query file as vmatch input.
            query_chunk_path = self.fasta_query_path
            unprocessed_chunk_dict[unprocessed_chunk_num]['query_path'] = query_chunk_path
            output_chunk_path = os.path.join(self.temp_dir, 'output.txt')
            unprocessed_chunk_dict[unprocessed_chunk_num]['output_path'] = output_chunk_path
            # Change the chunk stop to the number of lines in the query file.
            chunk_stop = total_lines

            if match_mode == 'exact_query_substring':
                command = ['vmatch',
                           '-q', f'{query_chunk_path}',
                           '-complete',
                           '-identity', '100',
                           '-nodist', '-noevalue', '-noidentity', '-noscore', # do not include alignment distance, E-value, identity and score in the output, but include the alignment distance
                           '-showdesc', '0', # query and target sequence names included in the output
                           '-v', # verbose commented output
                           f'{self.index_path}']
            elif match_mode == 'query_substring_with_mismatches':
                command = ['vmatch',
                           '-q', f'{query_chunk_path}',
                           '-complete',
                           '-h', f'{self.max_hamming_dist}',
                           '-identity', f'{self.min_ident}',
                           '-nodist', '-noevalue', '-noidentity', '-noscore',
                           '-showdesc', '0',
                           '-s', f'{self.align_output_length}',
                           '-v',
                           f'{self.index_path}']
            elif match_mode == 'query_substring_with_indels':
                command = ['vmatch',
                           '-q', f'{query_chunk_path}',
                           '-complete',
                           '-e', f'{self.max_edit_dist}',
                           '-identity', f'{self.min_ident}',
                           '-nodist', '-noevalue', '-noidentity', '-noscore',
                           '-showdesc', '0',
                           '-s', f'{self.align_output_length}',
                           '-v',
                           f'{self.index_path}']
            output_chunk_file = open(output_chunk_path, 'w', encoding='utf-8')
            unprocessed_chunk_dict[unprocessed_chunk_num]['output_file'] = output_chunk_file
            subprocess = utils.start_command(command, self.log_path, stdout=output_chunk_file, remove_log_file_if_exists=False)
            unprocessed_chunk_dict[unprocessed_chunk_num]['subprocess'] = subprocess

        # Chunk output parsing can lag behind vmatch, so distribute parsing jobs to multiple
        # processes. At any time during the search, the allocated cores will either be running
        # vmatch or parsing output.
        manager = multiprocessing.Manager()
        input_queue = manager.Queue()
        output_queue = manager.Queue()
        # Rather than passing parsed DataFrames through a pipe for each chunk, append them to a full
        # output table file. Use a lock to prevent multiple processes from writing at once.
        full_output_path = os.path.join(self.temp_dir, 'parsed_output.tsv')
        lock = multiprocessing.Lock()
        parsing_processes = [multiprocessing.Process(target=parsing_worker,
                                                     args=(input_queue,
                                                           output_queue,
                                                           self.align_output_length,
                                                           self.match_mode,
                                                           full_output_path,
                                                           lock,
                                                           self.edit_left_buffer,
                                                           self.edit_right_buffer))
                             for _ in range(self.num_threads)]
        for p in parsing_processes:
            p.start()
        num_unparsed_chunks = 0

        query_file = open(self.fasta_query_path)
        output_dfs = []
        if unprocessed_chunk_dict:
            self.progress.update_pid(pid)
            self.progress.update(f"Queuing queries 1-{pp_total_lines}/{pp_total_lines}")
        while True:
            # Stay in the loop while there are chunks to be processed. Fill up as many available
            # cores as possible with vmatch processes.
            if len(unprocessed_chunk_dict) + num_unparsed_chunks < self.num_threads:
                if chunk_start < total_lines:
                    # Write the next chunk of query sequences to a new file.
                    query_chunk_path = os.path.join(self.temp_dir, 'query_chunk_' + str(unprocessed_chunk_num) + '.fa')
                    unprocessed_chunk_dict[unprocessed_chunk_num]['query_path'] = query_chunk_path
                    with open(query_chunk_path, 'w') as query_chunk_file:
                        while running_line_count < chunk_stop:
                            query_chunk_file.write(next(query_file))
                            running_line_count += 1

                    output_chunk_path = os.path.join(self.temp_dir, 'output_' + str(unprocessed_chunk_num) + '.txt')
                    unprocessed_chunk_dict[unprocessed_chunk_num]['output_path'] = output_chunk_path

                    if match_mode == 'exact_query_substring':
                        command = ['vmatch',
                                   '-q', f'{query_chunk_path}',
                                   '-complete',
                                   '-identity', '100',
                                   '-nodist', '-noevalue', '-noidentity', '-noscore',
                                   '-showdesc', '0',
                                   '-v',
                                   f'{self.index_path}']
                    elif match_mode == 'query_substring_with_mismatches':
                        command = ['vmatch',
                                   '-q', f'{query_chunk_path}',
                                   '-complete',
                                   '-h', f'{self.max_hamming_dist}',
                                   '-identity', f'{self.min_ident}',
                                   '-nodist', '-noevalue', '-noidentity', '-noscore',
                                   '-showdesc', '0',
                                   '-s', f'{self.align_output_length}',
                                   '-v',
                                   f'{self.index_path}']
                    elif match_mode == 'query_substring_with_indels':
                        command = ['vmatch',
                                   '-q', f'{query_chunk_path}',
                                   '-complete',
                                   '-e', f'{self.max_edit_dist}',
                                   '-identity', f'{self.min_ident}',
                                   '-nodist', '-noevalue', '-noidentity', '-noscore',
                                   '-showdesc', '0',
                                   '-s', f'{self.align_output_length}',
                                   '-v',
                                   f'{self.index_path}']
                    output_chunk_file = open(output_chunk_path, 'w', encoding='utf-8')
                    unprocessed_chunk_dict[unprocessed_chunk_num]['output_file'] = output_chunk_file
                    subprocess = utils.start_command(command, self.log_path, stdout=output_chunk_file, remove_log_file_if_exists=False)
                    unprocessed_chunk_dict[unprocessed_chunk_num]['subprocess'] = subprocess

                    self.progress.update_pid(pid)
                    self.progress.update(f"Queuing queries {pp(chunk_start + 1)}-{pp(chunk_stop)}/{pp_total_lines}")

                    chunk_start = chunk_stop if chunk_stop < total_lines else total_lines
                    chunk_stop = chunk_start + query_chunk_size if chunk_start + query_chunk_size <= total_lines else total_lines
                    unprocessed_chunk_num += 1

            # Check if any chunks have finished.
            ret_val = None
            for chunk_num, chunk_subdict in unprocessed_chunk_dict.items():
                ret_val = chunk_subdict['subprocess'].poll()
                if ret_val == 0:
                    break
            else:
                # Slow down the loop.
                time.sleep(0.1)

            if ret_val == 0:
                # Parse the output of a finished chunk.
                chunk_subdict['output_file'].close()
                input_queue.put(chunk_subdict['output_path'])
                num_unparsed_chunks += 1

                if query_chunked:
                    # If chunking occurred, remove the chunked query file.
                    os.remove(chunk_subdict['query_path'])
                unprocessed_chunk_dict.pop(chunk_num)

            if num_unparsed_chunks:
                try:
                    output_queue.get(False)
                    # Make room for another vmatch job.
                    num_unparsed_chunks -= 1
                except queue.Empty:
                    pass

            if (chunk_start >= total_lines) and (not unprocessed_chunk_dict) and (num_unparsed_chunks == 0):
                # All queries have been processed.
                break
        query_file.close()

        for p in parsing_processes:
            p.terminate()
            p.join()

        # Load the concatenated chunked output.
        self.progress.update_pid(pid)
        self.progress.update("Finalizing matches")
        if self.match_mode == 'exact_query_substring':
            output_df = pd.read_csv(full_output_path, sep='\t', header=None, names=['query_name',
                                                                                    'target_name',
                                                                                    'query_start_in_target',
                                                                                    'query_length'])
        elif self.match_mode == 'query_substring_with_mismatches':
            output_df = pd.read_csv(full_output_path, sep='\t', header=None, names=['query_name',
                                                                                    'target_name',
                                                                                    'query_start_in_target',
                                                                                    'mismatch_positions'])
        elif self.match_mode == 'query_substring_with_indels':
            output_df = pd.read_csv(full_output_path,
                                    sep='\t',
                                    header=None,
                                    names=['query_name',
                                           'target_name',
                                           'query_start_in_target',
                                           'query_align_insert_starts',
                                           'target_align_insert_starts',
                                           'insert_lengths',
                                           'query_align_del_starts',
                                           'target_align_del_starts',
                                           'del_lengths'],
                                    dtype={'query_align_insert_starts': str,
                                           'target_align_insert_starts': str,
                                           'insert_lengths': str,
                                           'query_align_del_starts': str,
                                           'target_align_del_starts': str,
                                           'del_lengths': str}).fillna('')
        self.progress.end()

        if anvio.DEBUG:
            self.run.warning(f"The temp directory, '{self.temp_dir}', is kept. Don't forget to clean it up later!", header="Debug")
        elif not self.keep_temp_dir:
            if not self.quiet:
                self.run.info_single("Cleaning up the temp directory (use `--debug` to keep it for testing purposes)", nl_before=1, nl_after=1)
            shutil.rmtree(self.temp_dir)

        return output_df


def parsing_worker(input_queue, output_queue, align_output_length, match_mode, full_output_path, lock, edit_left_buffer=None, edit_right_buffer=None):
    """This worker is called from `Vmatch.search_queries` to parse chunked vmatch output. It is
    located outside `Vmatch` to allow multiprocessing."""
    while True:
        output_path = input_queue.get()
        if match_mode == 'exact_query_substring':
            output_df = parse_exact_query_substrings(output_path)
        elif match_mode == 'query_substring_with_mismatches':
            output_df = parse_query_substrings_with_mismatches(output_path, align_output_length)
        elif match_mode == 'query_substring_with_indels':
            if edit_left_buffer is None:
                edit_left_buffer = 0
            if edit_right_buffer is None:
                edit_right_buffer = 0
            output_df = parse_query_substrings_with_indels(output_path, align_output_length, edit_left_buffer, edit_right_buffer)

        with lock:
            # Sending the DataFrame through the multiprocessing pipe is much slower than writing
            # to file.
            output_df.to_csv(full_output_path, mode='a', sep='\t', index=False, header=False)

        os.remove(output_path)
        output_queue.put(True)


def parse_exact_query_substrings(output_path):
    """Parse queries that are exact substrings of one or more targets. All this method does to the
    output is remove commented lines and extraneous columns."""
    output_df = pd.read_csv(output_path,
                            delim_whitespace=True, # there is a variable number of spaces between fields, and there are spaces at the beginning of each line
                            comment='#', # remove comment lines at the beginning of the file
                            header=None,
                            usecols=[0, 1, 2, 4, 5, 6], # disregard the fourth column, indicating a direct match, "D", or palindromic match, "P", as the latter does not occur here
                            names=['target_align_length', 'target_name', 'query_start_in_target', 'query_length', 'query_name', 'align_start_in_query'])

    return output_df[['query_name', 'target_name', 'query_start_in_target', 'query_length']]


def parse_query_substrings_with_mismatches(output_path, align_output_length):
    """Parse queries that are fully contained by targets and that have mismatches but not gaps.
    Exact matches without any mismatches are not retained! Finding the positions of mismatches
    requires parsing actual sequence alignments in addition to rows of summary data."""
    align_records = []
    output_file = open(output_path)
    for line in output_file:
        if line[0] != '#':
            # Ignore comment lines at the beginning of the file.
            break
    # The first non-comment line contains summary data for the first alignment.
    align_record = line.rstrip().split()

    # The following is True when the alignment lines are being parsed.
    parsing_align = True
    # The first line of the alignment is (part of) the target sequence.
    prev_line_was_target = False
    # Alignments longer than the determined width are split between lines with a blank line between
    # sections. A second blank line signifies the end of the alignment record.
    prev_line_was_blank_after_align = False
    align_index = 0
    mismatch_positions = []
    for line in output_file:
        if parsing_align:
            if line[0] == 'S':
                prev_line_was_blank_after_align = False
                prev_line_was_target = True
            elif line[0] == 'Q':
                if prev_line_was_target:
                    prev_line_was_target = False
                    # If the query line comes immediately after the target line, then this section
                    # of the alignment has no mismatches.
                    align_index += align_output_length
            elif prev_line_was_target:
                # A line between the target and query is under consideration. Its existence means
                # there is a mismatch in this section of the alignment.
                for section_index, char in enumerate(line[7: 7 + align_output_length]):
                    if char == '!':
                        mismatch_positions.append(str(align_index + section_index))
                prev_line_was_target = False
                align_index += align_output_length
            elif prev_line_was_blank_after_align:
                # A second blank line has now been encountered, indicating the end of a record.
                if mismatch_positions:
                    # Retain the query name, target name, query start in target, and mismatch positions in the query.
                    align_records.append((align_record[5],
                                          align_record[1],
                                          align_record[2],
                                          ','.join(mismatch_positions)))
                    mismatch_positions = []
                parsing_align = False
                prev_line_was_blank_after_align = False
                align_index = 0
                continue
            else:
                prev_line_was_blank_after_align = True
                continue
        else:
            align_record = line.rstrip().split()
            parsing_align = True
    output_file.close()

    return pd.DataFrame(align_records, columns=['query_name',
                                                'target_name',
                                                'query_start_in_target',
                                                'mismatch_positions'])


def parse_query_substrings_with_indels(output_path, align_output_length, edit_left_buffer, edit_right_buffer):
    """Parse queries that are fully contained by targets and that have indels but not mismatches.
    Exact matches without any indels are not retained! Finding the positions of indels requires
    parsing actual sequence alignments in addition to rows of summary data."""
    align_records = []
    output_file = open(output_path)
    for line in output_file:
        if line[0] != '#':
            # Ignore comment lines at the beginning of the file.
            break
    # The first non-comment line contains summary data for the first alignment.
    align_record = line.rstrip().split()

    # The following is True when the alignment lines are being parsed.
    parsing_align = True
    # The first line of the alignment is (part of) the target sequence.
    prev_line_was_target = False
    # Alignments longer than the determined width are split between lines with a blank line between
    # sections. A second blank line signifies the end of the alignment record.
    prev_line_was_blank_after_align = False
    target_align = ''
    query_align = ''
    align_index = 0
    edit_positions = []
    for line in output_file:
        if parsing_align:
            if line[0] == 'S':
                target_align += line[7: 7 + align_output_length]
                prev_line_was_blank_after_align = False
                prev_line_was_target = True
            elif line[0] == 'Q':
                query_align += line[7: 7 + align_output_length]
                if prev_line_was_target:
                    prev_line_was_target = False
                    # If the query line comes immediately after the target line, then this section
                    # of the alignment has no mismatches or indels.
                    align_index += align_output_length
            elif prev_line_was_target:
                # A line between the target and query is under consideration. Its existence means
                # there is a mismatch or indel in this section of the alignment.
                for section_index, char in enumerate(line[7: 7 + align_output_length]):
                    if char == '!':
                        edit_positions.append(align_index + section_index)
                prev_line_was_target = False
                align_index += align_output_length
            elif prev_line_was_blank_after_align:
                # A second blank line has now been encountered, indicating the end of a record.
                if edit_positions:
                    align_records.extend(get_alignments_with_indels(align_record,
                                                                    query_align.rstrip(' '),
                                                                    target_align.rstrip(' '),
                                                                    edit_positions,
                                                                    edit_left_buffer,
                                                                    edit_right_buffer))

                parsing_align = False
                prev_line_was_blank_after_align = False
                align_index = 0
                query_align = ''
                target_align = ''
                edit_positions = []
                continue
            else:
                prev_line_was_blank_after_align = True
                continue
        else:
            align_record = line.rstrip().split()
            parsing_align = True
    output_file.close()

    return pd.DataFrame(align_records, columns=['query_name',
                                                'target_name',
                                                'query_start_in_target',
                                                'query_align_insert_starts',
                                                'target_align_insert_starts',
                                                'insert_lengths',
                                                'query_align_del_starts',
                                                'target_align_del_starts',
                                                'del_lengths'])


def get_alignments_with_indels(align_record, query_align, target_align, edit_positions, edit_left_buffer, edit_right_buffer):
    """Determine the positions of indels in the alignment. Reject alignments with any mismatches.
    Also reject alignments that have indels within the left and right buffer distances of the
    respective ends of the alignment."""
    # Find coordinates and lengths of indels in both query and target seqs.
    # For example:
    #     0123456  789012345
    # Q: -ACGTTAC--GTACGACTA
    #    0123 4567890123  4
    # S: TACG-TACGCGTACG--T-
    # Del positions in aligned query, relative to the nt before the gap, which can result in a value of -1.
    # [-1, 6]
    # Del positions in aligned target:
    # [0, 7]
    # Del lengths:
    # [1, 2]
    # Insertion positions in aligned query:
    # [3, 12, 15]
    # Insertion positions in aligned target, relative to the nt before the gap:
    # [3, 13, 14]
    # Insertion lengths:
    # [1, 2, 1]

    # There may be equivalent alignments that are not reported by Vmatch but which are recorded
    # here. Vmatch reports alignments with the gap in the first possible positions. Exclamation
    # marks show the transposable gaps from the same example as before:
    #        !!
    # Q: -ACGTTAC--GTACGACTA
    # S: TACGT-ACGCGTACG--T-

    #            !!!
    # Q: -ACGTTACG--TACGACTA
    # S: TACG-TACGCGTACG--T-

    #        !!  !!!
    # Q: -ACGTTACG--TACGACTA
    # S: TACGT-ACGCGTACG--T-

    align_records = []
    indel_configs = [{'query_align_insert_starts': [],
                      'target_align_insert_starts': [],
                      'insert_lengths': [],
                      'query_align_del_starts': [],
                      'target_align_del_starts': [],
                      'del_lengths': []}]
    insert_pos_in_align = -1
    insert_length = 0
    del_pos_in_align = -1
    del_length = 0
    # Count the gaps recorded through the query and target.
    running_query_gap_count = 0
    running_target_gap_count = 0
    edit_right_bound = len(query_align) - edit_right_buffer
    for edit_pos in edit_positions:
        if not edit_left_buffer <= edit_pos < edit_right_bound:
            # Do not allow indels within the left and right buffer distances of the
            # respective ends of the alignment.
            break

        if target_align[edit_pos] == '-': # insertion
            if del_length:
                # The insertion is the next edit after a del. The del is recorded here.
                query_align_del_start = del_pos_in_align - running_query_gap_count + 1
                target_align_del_start = del_pos_in_align - running_target_gap_count - del_length + 1
                del_extensions = []
                for indel_config_dict in indel_configs:
                    if indel_config_dict['query_align_del_starts']:
                        if query_align_del_start == indel_config_dict['query_align_del_starts'][-1]:
                            # The del being recorded is actually adjacent to and part of the last
                            # recorded del. This can happen when the last recorded del was an
                            # "alternative" configuration. Therefore, extend the last recorded del.
                            indel_config_dict['del_lengths'][-1] = indel_config_dict['del_lengths'][-1] + del_length
                            del_extensions.append(True)
                            continue
                    indel_config_dict['query_align_del_starts'].append(query_align_del_start)
                    indel_config_dict['target_align_del_starts'].append(target_align_del_start)
                    indel_config_dict['del_lengths'].append(del_length)
                    del_extensions.append(False)

                # Record any alternative configurations of the del.
                pos_in_align = del_pos_in_align
                alt_configs = []
                while True:
                    # Loop through each subsequent nt in the alignment.
                    next_pos_in_align = pos_in_align + 1
                    try:
                        next_query_value = query_align[next_pos_in_align]
                    except IndexError:
                        # The del is at the end of the alignment.
                        break

                    if next_pos_in_align >= edit_right_bound:
                        # A transposed del would violate the right bound.
                        break

                    if next_query_value != target_align[next_pos_in_align - del_length]:
                        # The next nt in the query does not match the nt in the target to which it
                        # would be paired, so the del cannot be shifted without creating a mismatch.
                        break

                    # Make a new indel config with the present del shifted to the right by ≥1
                    # positions.
                    pos_in_align = next_pos_in_align
                    indel_shift = pos_in_align - del_pos_in_align
                    for basis_indel_config_dict, del_extension in zip(indel_configs, del_extensions):
                        if del_extension:
                            # The shifted del is no longer an extension of the previous del, so
                            # split it off.
                            alt_configs.append({'query_align_insert_starts': basis_indel_config_dict['query_align_insert_starts'].copy(),
                                                'target_align_insert_starts': basis_indel_config_dict['target_align_insert_starts'].copy(),
                                                'insert_lengths': basis_indel_config_dict['insert_lengths'].copy(),
                                                'query_align_del_starts': basis_indel_config_dict['query_align_del_starts'] + [query_align_del_start + indel_shift],
                                                'target_align_del_starts': basis_indel_config_dict['target_align_del_starts'] + [target_align_del_start + indel_shift],
                                                'del_lengths': basis_indel_config_dict['del_lengths'][: -1] + [basis_indel_config_dict['del_lengths'][-1] - del_length, del_length]})
                        else:
                            alt_configs.append({'query_align_insert_starts': basis_indel_config_dict['query_align_insert_starts'].copy(),
                                                'target_align_insert_starts': basis_indel_config_dict['target_align_insert_starts'].copy(),
                                                'insert_lengths': basis_indel_config_dict['insert_lengths'].copy(),
                                                'query_align_del_starts': basis_indel_config_dict['query_align_del_starts'][: -1] + [basis_indel_config_dict['query_align_del_starts'][-1] + indel_shift],
                                                'target_align_del_starts': basis_indel_config_dict['target_align_del_starts'][: -1] + [basis_indel_config_dict['target_align_del_starts'][-1] + indel_shift],
                                                'del_lengths': basis_indel_config_dict['del_lengths'].copy()})
                indel_configs.extend(alt_configs)

                del_pos_in_align = -1
                del_length = 0
            elif insert_length:
                # A distinct insertion rather than a del or mismatch is the next edit after an
                # insertion. The prior insertion is recorded here.
                if edit_pos - insert_pos_in_align > 1:
                    query_align_insert_start = insert_pos_in_align - running_query_gap_count - insert_length + 1
                    target_align_insert_start = insert_pos_in_align - running_target_gap_count + 1
                    insert_extensions = []
                    for indel_config_dict in indel_configs:
                        if indel_config_dict['query_align_insert_starts']:
                            if query_align_insert_start == indel_config_dict['query_align_insert_starts'][-1] + indel_config_dict['insert_lengths'][-1]:
                                # The insertion being recorded is actually adjacent to and part of the
                                # last recorded insertion. This can happen when the last recorded
                                # insertion was an "alternative" configuration. Therefore, extend the
                                # last recorded insertion.
                                indel_config_dict['insert_lengths'][-1] = indel_config_dict['insert_lengths'][-1] + insert_length
                                insert_extensions.append(True)
                                continue
                        indel_config_dict['query_align_insert_starts'].append(query_align_insert_start)
                        indel_config_dict['target_align_insert_starts'].append(target_align_insert_start)
                        indel_config_dict['insert_lengths'].append(insert_length)
                        insert_extensions.append(False)

                    # Record any alternative configurations of the insertion.
                    pos_in_align = insert_pos_in_align
                    alt_configs = []
                    while True:
                        # Loop through each subsequent nt in the alignment.
                        next_pos_in_align = pos_in_align + 1
                        try:
                            next_target_value = target_align[next_pos_in_align]
                        except IndexError:
                            # The insertion is at the end of the alignment.
                            break

                        if next_pos_in_align >= edit_right_bound:
                            # A transposed insertion would violate the right bound.
                            break

                        if next_target_value != query_align[next_pos_in_align - insert_length]:
                            # The next nt in the target does not match the nt in the query to which
                            # it would be paired, so the insertion cannot be shifted without
                            # creating a mismatch.
                            break

                        # Make a new indel config with the present insertion shifted to the right by
                        # ≥1 positions.
                        pos_in_align = next_pos_in_align
                        indel_shift = pos_in_align - insert_pos_in_align
                        for basis_indel_config_dict, insert_extension in zip(indel_configs, insert_extensions):
                            if insert_extension:
                                # The shifted insertion is no longer an extension of the previous
                                # insertion, so split it off.
                                alt_configs.append({'query_align_insert_starts': basis_indel_config_dict['query_align_insert_starts'] + [query_align_insert_start + indel_shift],
                                                    'target_align_insert_starts': basis_indel_config_dict['target_align_insert_starts'] + [target_align_insert_start + indel_shift],
                                                    'insert_lengths': basis_indel_config_dict['insert_lengths'][: -1] + [basis_indel_config_dict['insert_lengths'][-1] - insert_length, insert_length],
                                                    'query_align_del_starts': basis_indel_config_dict['query_align_del_starts'].copy(),
                                                    'target_align_del_starts': basis_indel_config_dict['target_align_del_starts'].copy(),
                                                    'del_lengths': basis_indel_config_dict['del_lengths'].copy()})
                            else:
                                alt_configs.append({'query_align_insert_starts': basis_indel_config_dict['query_align_insert_starts'][: -1] + [basis_indel_config_dict['query_align_insert_starts'][-1] + indel_shift],
                                                    'target_align_insert_starts': basis_indel_config_dict['target_align_insert_starts'][: -1] + [basis_indel_config_dict['target_align_insert_starts'][-1] + indel_shift],
                                                    'insert_lengths': basis_indel_config_dict['insert_lengths'].copy(),
                                                    'query_align_del_starts': basis_indel_config_dict['query_align_del_starts'].copy(),
                                                    'target_align_del_starts': basis_indel_config_dict['target_align_del_starts'].copy(),
                                                    'del_lengths': basis_indel_config_dict['del_lengths'].copy()})
                    indel_configs.extend(alt_configs)

                    insert_length = 0
            running_target_gap_count += 1
            insert_pos_in_align = edit_pos
            insert_length += 1
        elif query_align[edit_pos] == '-': # del
            if insert_length:
                # The del is the next edit after an insertion. The insertion is recorded here.
                query_align_insert_start = insert_pos_in_align - running_query_gap_count - insert_length + 1
                target_align_insert_start = insert_pos_in_align - running_target_gap_count + 1
                insert_extensions = []
                for indel_config_dict in indel_configs:
                    if indel_config_dict['query_align_insert_starts']:
                        if query_align_insert_start == indel_config_dict['query_align_insert_starts'][-1] + indel_config_dict['insert_lengths'][-1]:
                            # The insertion being recorded is actually adjacent to and part of the
                            # last recorded insertion. This can happen when the last recorded
                            # insertion was an "alternative" configuration. Therefore, extend the
                            # last recorded insertion.
                            indel_config_dict['insert_lengths'][-1] = indel_config_dict['insert_lengths'][-1] + insert_length
                            insert_extensions.append(True)
                            continue
                    indel_config_dict['query_align_insert_starts'].append(query_align_insert_start)
                    indel_config_dict['target_align_insert_starts'].append(target_align_insert_start)
                    indel_config_dict['insert_lengths'].append(insert_length)
                    insert_extensions.append(False)

                # Record any alternative configurations of the insertion.
                pos_in_align = insert_pos_in_align
                alt_configs = []
                while True:
                    # Loop through each subsequent nt in the alignment.
                    next_pos_in_align = pos_in_align + 1
                    try:
                        next_target_value = target_align[next_pos_in_align]
                    except IndexError:
                        # The insertion is at the end of the alignment.
                        break

                    if next_pos_in_align >= edit_right_bound:
                        # A transposed insertion would violate the right bound.
                        break

                    if next_target_value != query_align[next_pos_in_align - insert_length]:
                        # The next nt in the target does not match the nt in the query to which it
                        # would be paired, so the insertion cannot be shifted without creating a
                        # mismatch.
                        break

                    # Make a new indel config with the present insertion shifted
                    # to the right by ≥1 positions.
                    pos_in_align = next_pos_in_align
                    indel_shift = pos_in_align - insert_pos_in_align
                    for basis_indel_config_dict, insert_extension in zip(indel_configs, insert_extensions):
                        if insert_extension:
                            alt_configs.append({'query_align_insert_starts': basis_indel_config_dict['query_align_insert_starts'] + [query_align_insert_start + indel_shift],
                                                'target_align_insert_starts': basis_indel_config_dict['target_align_insert_starts'] + [target_align_insert_start + indel_shift],
                                                'insert_lengths': basis_indel_config_dict['insert_lengths'][: -1] + [basis_indel_config_dict['insert_lengths'][-1] - insert_length, insert_length],
                                                'query_align_del_starts': basis_indel_config_dict['query_align_del_starts'].copy(),
                                                'target_align_del_starts': basis_indel_config_dict['target_align_del_starts'].copy(),
                                                'del_lengths': basis_indel_config_dict['del_lengths'].copy()})
                        else:
                            alt_configs.append({'query_align_insert_starts': basis_indel_config_dict['query_align_insert_starts'][: -1] + [basis_indel_config_dict['query_align_insert_starts'][-1] + indel_shift],
                                                'target_align_insert_starts': basis_indel_config_dict['target_align_insert_starts'][: -1] + [basis_indel_config_dict['target_align_insert_starts'][-1] + indel_shift],
                                                'insert_lengths': basis_indel_config_dict['insert_lengths'].copy(),
                                                'query_align_del_starts': basis_indel_config_dict['query_align_del_starts'].copy(),
                                                'target_align_del_starts': basis_indel_config_dict['target_align_del_starts'].copy(),
                                                'del_lengths': basis_indel_config_dict['del_lengths'].copy()})
                indel_configs.extend(alt_configs)

                insert_pos_in_align = -1
                insert_length = 0
            elif del_length:
                # A distinct del rather than an insertion or mismatch is the next edit after a del.
                # The prior del is recorded here.
                if edit_pos - del_pos_in_align > 1:
                    query_align_del_start = del_pos_in_align - running_query_gap_count + 1
                    target_align_del_start = del_pos_in_align - running_target_gap_count - del_length + 1
                    del_extensions = []
                    for indel_config_dict in indel_configs:
                        if indel_config_dict['query_align_del_starts']:
                            if query_align_del_start == indel_config_dict['query_align_del_starts'][-1]:
                                # The del being recorded is actually adjacent to and part of the last
                                # recorded del. This can happen when the last recorded del was an
                                # "alternative" configuration. Therefore, extend the last recorded del.
                                indel_config_dict['del_lengths'][-1] = indel_config_dict['del_lengths'][-1] + del_length
                                del_extensions.append(True)
                                continue
                        indel_config_dict['query_align_del_starts'].append(query_align_del_start)
                        indel_config_dict['target_align_del_starts'].append(target_align_del_start)
                        indel_config_dict['del_lengths'].append(del_length)
                        del_extensions.append(False)

                    # Record any alternative configurations of the del.
                    pos_in_align = del_pos_in_align
                    alt_configs = []
                    while True:
                        # Loop through each subsequent nt in the alignment.
                        next_pos_in_align = pos_in_align + 1
                        try:
                            next_query_value = query_align[next_pos_in_align]
                        except IndexError:
                            # The del is at the end of the alignment.
                            break

                        if next_pos_in_align >= edit_right_bound:
                            # A transposed del would violate the right bound.
                            break

                        if next_query_value != target_align[next_pos_in_align - del_length]:
                            # The next nt in the query does not match the nt in the target to which
                            # it would be paired, so the del cannot be shifted without creating a
                            # mismatch.
                            break

                        # Make a new indel config with the present del shifted to the right by ≥1
                        # positions.
                        pos_in_align = next_pos_in_align
                        indel_shift = pos_in_align - del_pos_in_align
                        for basis_indel_config_dict, del_extension in zip(indel_configs, del_extensions):
                            if del_extension:
                                alt_configs.append({'query_align_insert_starts': basis_indel_config_dict['query_align_insert_starts'].copy(),
                                                    'target_align_insert_starts': basis_indel_config_dict['target_align_insert_starts'].copy(),
                                                    'insert_lengths': basis_indel_config_dict['insert_lengths'].copy(),
                                                    'query_align_del_starts': basis_indel_config_dict['query_align_del_starts'] + [query_align_del_start + indel_shift],
                                                    'target_align_del_starts': basis_indel_config_dict['target_align_del_starts'] + [target_align_del_start + indel_shift],
                                                    'del_lengths': basis_indel_config_dict['del_lengths'][: -1] + [basis_indel_config_dict['del_lengths'][-1] - del_length, del_length]})
                            else:
                                alt_configs.append({'query_align_insert_starts': basis_indel_config_dict['query_align_insert_starts'].copy(),
                                                    'target_align_insert_starts': basis_indel_config_dict['target_align_insert_starts'].copy(),
                                                    'insert_lengths': basis_indel_config_dict['insert_lengths'].copy(),
                                                    'query_align_del_starts': basis_indel_config_dict['query_align_del_starts'][: -1] + [basis_indel_config_dict['query_align_del_starts'][-1] + indel_shift],
                                                    'target_align_del_starts': basis_indel_config_dict['target_align_del_starts'][: -1] + [basis_indel_config_dict['target_align_del_starts'][-1] + indel_shift],
                                                    'del_lengths': basis_indel_config_dict['del_lengths'].copy()})
                    indel_configs.extend(alt_configs)

                    del_length = 0
            running_query_gap_count += 1
            del_pos_in_align = edit_pos
            del_length += 1
        else:
            # The edit is a mismatch.
            break
    else:
        # Record the last edit in the alignment.
        if insert_length:
            # The last edit was an insertion.
            query_align_insert_start = insert_pos_in_align - running_query_gap_count - insert_length + 1
            target_align_insert_start = insert_pos_in_align - running_target_gap_count + 1
            insert_extensions = []
            for indel_config_dict in indel_configs:
                if indel_config_dict['query_align_insert_starts']:
                    if query_align_insert_start == indel_config_dict['query_align_insert_starts'][-1] + indel_config_dict['insert_lengths'][-1]:
                        # The insertion being recorded is actually adjacent to and part of the
                        # last recorded insertion. This can happen when the last recorded
                        # insertion was an "alternative" configuration. Therefore, extend the
                        # last recorded insertion.
                        indel_config_dict['insert_lengths'][-1] = indel_config_dict['insert_lengths'][-1] + insert_length
                        insert_extensions.append(True)
                        continue
                indel_config_dict['query_align_insert_starts'].append(query_align_insert_start)
                indel_config_dict['target_align_insert_starts'].append(target_align_insert_start)
                indel_config_dict['insert_lengths'].append(insert_length)
                insert_extensions.append(False)

            # Record any alternative configurations of the insertion.
            pos_in_align = insert_pos_in_align
            alt_configs = []
            while True:
                # Loop through each subsequent nt in the alignment.
                next_pos_in_align = pos_in_align + 1
                try:
                    next_target_value = target_align[next_pos_in_align]
                except IndexError:
                    # The insertion is at the end of the alignment.
                    break

                if next_pos_in_align >= edit_right_bound:
                    # A transposed insertion would violate the right bound.
                    break

                if next_target_value != query_align[next_pos_in_align - insert_length]:
                    # The next nt in the target does not match the nt in the query to which it
                    # would be paired, so the insertion cannot be shifted without creating a
                    # mismatch.
                    break

                # Make a new indel config with the present insertion shifted to the right by ≥1
                # positions.
                pos_in_align = next_pos_in_align
                indel_shift = pos_in_align - insert_pos_in_align
                for basis_indel_config_dict, insert_extension in zip(indel_configs, insert_extensions):
                    if insert_extension:
                        alt_configs.append({'query_align_insert_starts': basis_indel_config_dict['query_align_insert_starts'] + [query_align_insert_start + indel_shift],
                                            'target_align_insert_starts': basis_indel_config_dict['target_align_insert_starts'] + [target_align_insert_start + indel_shift],
                                            'insert_lengths': basis_indel_config_dict['insert_lengths'][: -1] + [basis_indel_config_dict['insert_lengths'][-1] - insert_length, insert_length],
                                            'query_align_del_starts': basis_indel_config_dict['query_align_del_starts'].copy(),
                                            'target_align_del_starts': basis_indel_config_dict['target_align_del_starts'].copy(),
                                            'del_lengths': basis_indel_config_dict['del_lengths'].copy()})
                    else:
                        alt_configs.append({'query_align_insert_starts': basis_indel_config_dict['query_align_insert_starts'][: -1] + [basis_indel_config_dict['query_align_insert_starts'][-1] + indel_shift],
                                            'target_align_insert_starts': basis_indel_config_dict['target_align_insert_starts'][: -1] + [basis_indel_config_dict['target_align_insert_starts'][-1] + indel_shift],
                                            'insert_lengths': basis_indel_config_dict['insert_lengths'].copy(),
                                            'query_align_del_starts': basis_indel_config_dict['query_align_del_starts'].copy(),
                                            'target_align_del_starts': basis_indel_config_dict['target_align_del_starts'].copy(),
                                            'del_lengths': basis_indel_config_dict['del_lengths'].copy()})
            indel_configs.extend(alt_configs)
        elif del_length:
            # The last edit was a deletion.
            query_align_del_start = del_pos_in_align - running_query_gap_count + 1
            target_align_del_start = del_pos_in_align - running_target_gap_count - del_length + 1
            del_extensions = []
            for indel_config_dict in indel_configs:
                if indel_config_dict['query_align_del_starts']:
                    if query_align_del_start == indel_config_dict['query_align_del_starts'][-1]:
                        # The del being recorded is actually adjacent to and part of the last
                        # recorded del. This can happen when the last recorded del was an
                        # "alternative" configuration. Therefore, extend the last recorded del.
                        indel_config_dict['del_lengths'][-1] = indel_config_dict['del_lengths'][-1] + del_length
                        del_extensions.append(True)
                        continue
                indel_config_dict['query_align_del_starts'].append(query_align_del_start)
                indel_config_dict['target_align_del_starts'].append(target_align_del_start)
                indel_config_dict['del_lengths'].append(del_length)
                del_extensions.append(False)

            # Record any alternative configurations of the del.
            pos_in_align = del_pos_in_align
            alt_configs = []
            while True:
                # Loop through each subsequent nt in the alignment.
                next_pos_in_align = pos_in_align + 1
                try:
                    next_query_value = query_align[next_pos_in_align]
                except IndexError:
                    # The del is at the end of the alignment.
                    break

                if next_pos_in_align >= edit_right_bound:
                    # A transposed del would violate the right bound.
                    break

                if next_query_value != target_align[next_pos_in_align - del_length]:
                    # The next nt in the query does not match the nt in the target to which it
                    # would be paired, so the del cannot be shifted without creating a mismatch.
                    break

                # Make a new indel config with the present del shifted to the right by ≥1 positions.
                pos_in_align = next_pos_in_align
                indel_shift = pos_in_align - del_pos_in_align
                for basis_indel_config_dict, del_extension in zip(indel_configs, del_extensions):
                    if del_extension:
                        alt_configs.append({'query_align_insert_starts': basis_indel_config_dict['query_align_insert_starts'].copy(),
                                            'target_align_insert_starts': basis_indel_config_dict['target_align_insert_starts'].copy(),
                                            'insert_lengths': basis_indel_config_dict['insert_lengths'].copy(),
                                            'query_align_del_starts': basis_indel_config_dict['query_align_del_starts'] + [query_align_del_start + indel_shift],
                                            'target_align_del_starts': basis_indel_config_dict['target_align_del_starts'] + [target_align_del_start + indel_shift],
                                            'del_lengths': basis_indel_config_dict['del_lengths'][: -1] + [basis_indel_config_dict['del_lengths'][-1] - del_length, del_length]})
                    else:
                        alt_configs.append({'query_align_insert_starts': basis_indel_config_dict['query_align_insert_starts'].copy(),
                                            'target_align_insert_starts': basis_indel_config_dict['target_align_insert_starts'].copy(),
                                            'insert_lengths': basis_indel_config_dict['insert_lengths'].copy(),
                                            'query_align_del_starts': basis_indel_config_dict['query_align_del_starts'][: -1] + [basis_indel_config_dict['query_align_del_starts'][-1] + indel_shift],
                                            'target_align_del_starts': basis_indel_config_dict['target_align_del_starts'][: -1] + [basis_indel_config_dict['target_align_del_starts'][-1] + indel_shift],
                                            'del_lengths': basis_indel_config_dict['del_lengths'].copy()})
            indel_configs.extend(alt_configs)

        if indel_configs[0]['insert_lengths'] or indel_configs[0]['del_lengths']:
            # Retain the query name, target name, and query start in target in
            # addition to edit information.
            query_name = align_record[5]
            target_name = align_record[1]
            query_start = align_record[2]
            for indel_config_dict in indel_configs:
                align_records.append((query_name,
                                      target_name,
                                      query_start,
                                      ','.join(map(str, indel_config_dict['query_align_insert_starts'])),
                                      ','.join(map(str, indel_config_dict['target_align_insert_starts'])),
                                      ','.join(map(str, indel_config_dict['insert_lengths'])),
                                      ','.join(map(str, indel_config_dict['query_align_del_starts'])),
                                      ','.join(map(str, indel_config_dict['target_align_del_starts'])),
                                      ','.join(map(str, indel_config_dict['del_lengths']))))
    return align_records
