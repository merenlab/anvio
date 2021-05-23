# coding utf-8
"""An interface for Vmatch"""

import os
import time
import queue
import shutil
import pandas as pd
import multiprocessing as mp

from collections import defaultdict

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2020, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Samuel Miller"
__email__ = "samuelmiller10@gmail.com"


class Vmatch(object):
    """Parent class for all Vmatch usage"""

    def __init__(self, args, run=None, progress=None, skip_sanity_check=False):
        self.index_program_name = 'mkvtree'
        self.search_program_name = 'vmatch'

        self.tested_versions = ('2.3.0', )
        self.supported_match_modes = ('exact_query_substring', 'query_substring_with_mismatches')

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
        self.query_chunk_size = A('query_chunk_size') or 100000 // self.num_threads

        self.max_hamming_dist = A('max_hamming_dist')
        self.min_ident = A('min_ident') or 100

        # With this option, sequence alignments are included in the output. They are wrapped to the
        # given length on each line.
        self.align_output_length = A('align_output_length')

        self.quiet = A('quiet')
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

        if self.match_mode == 'exact_query_substring':
            if self.max_hamming_dist is not None:
                raise ConfigError("The match mode, 'exact_query_substring', is incompatible with the maximum Hamming distance argument.")

            if self.align_output_length is not None:
                raise ConfigError("The match mode, 'exact_query_substring', is incompatible with the alignment output length argument.")
        else:
            if self.max_hamming_dist is None:
                raise ConfigError("The match mode, 'query_substring_with_mismatches', requires the maximum Hamming distance argument.")

            if self.max_hamming_dist % 1 != 0 or self.max_hamming_dist <= 1:
                raise ConfigError("The maximum allowed Hamming distance must be a positive integer.")

            if not self.align_output_length:
                raise ConfigError("The match mode, 'query_substring_with_mismatches', requires the argument, 'align_output_length'.")

        filesnpaths.is_file_exists(self.fasta_db_path)
        self.index_path = os.path.join(self.temp_dir, 'index')

        if self.min_ident % 1 != 0:
            raise ConfigError("The minimum identity of an alignment must be an integer between 1 and 100. "
                              f"The provided value was {self.min_ident}.")


    def check_programs(self, quiet=False):
        self.installed_index_program_version = None
        self.installed_search_program_version = None

        for program_name, attr_name in zip((self.index_program_name, self.search_program_name),
                                           ('installed_index_program_version', 'installed_search_program_version')):
            utils.is_program_exists(program_name)

            output, ret_code = utils.get_command_output_from_shell(f'{program_name} -version')
            try:
                version_found = output.split(b'\n')[0].split(b' ')[2].decode("utf-8")
                if not quiet:
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
            command = ('mkvtree '
                       f'-db {self.fasta_db_path} '
                       '-dna '
                       f'-indexname {self.index_path} '
                       '-pl '
                       f'-tis -suf -sti1 -lcp -bck '
                       '-v'
                       ).split(' ')
        else:
            # `-ois` produces an index table of original sequences and is needed to output sequence alignments.
            # `-skp` produces an index table used in finding whole queries with mismatches.
            command = ('mkvtree '
                       f'-db {self.fasta_db_path} '
                       '-dna '
                       f'-indexname {self.index_path} '
                       '-pl '
                       f'-tis -suf -sti1 -lcp -bck -ois -skp '
                       '-v'
                       ).split(' ')
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

        query_chunk_size = self.query_chunk_size * 2 # convert sequences to FASTA lines
        next_chunk_start = query_chunk_size
        unprocessed_chunk_num = 1
        unprocessed_chunk_dict = defaultdict(dict)
        running_line_count = 0
        with open(self.fasta_query_path) as query_file:
            for line_num, line in enumerate(query_file, 1):
                pass
        total_lines = line_num

        match_mode = self.match_mode
        if next_chunk_start >= total_lines:
            # All query sequences fit in one chunk, so use the original query file as vmatch input.
            query_chunk_path = self.fasta_query_path
            unprocessed_chunk_dict[unprocessed_chunk_num]['query_path'] = query_chunk_path
            output_chunk_path = os.path.join(self.temp_dir, 'output.txt')
            unprocessed_chunk_dict[unprocessed_chunk_num]['output_path'] = output_chunk_path

            if match_mode == 'exact_query_substring':
                command = ('vmatch '
                           f'-q {query_chunk_path} '
                           '-complete '
                           '-identity 100 '
                           '-nodist -noevalue -noidentity -noscore ' # do not include alignment distance, E-value, identity and score in the output, but include the alignment distance
                           '-showdesc 0 ' # query and target sequence names included in the output
                           '-v ' # verbose commented output
                           f'{self.index_path}').split(' ')
            elif match_mode == 'query_substring_with_mismatches':
                command = ('vmatch '
                           f'-q {query_chunk_path} '
                           '-complete '
                           f'-h {self.max_hamming_dist} '
                           f'-identity {self.min_ident} '
                           '-nodist -noevalue -noidentity -noscore '
                           '-showdesc 0 '
                           f'-s {self.align_output_length} '
                           '-v '
                           f'{self.index_path}').split(' ')
            output_chunk_file = open(output_chunk_path, 'w', encoding='utf-8')
            unprocessed_chunk_dict[unprocessed_chunk_num]['output_file'] = output_chunk_file
            subprocess = utils.start_command(command, self.log_path, stdout=output_chunk_file, remove_log_file_if_exists=False)
            unprocessed_chunk_dict[unprocessed_chunk_num]['subprocess'] = subprocess

            # Ensure that `next_chunk_start` is greater than the total number of lines in the query
            # file as a way of distinguishing this case from the chunking case.
            next_chunk_start += 1

        # Chunk output parsing can lag behind vmatch, so distribute parsing jobs to multiple
        # processes. At any time during the search, the allocated cores will either be running
        # vmatch or parsing output.
        manager = mp.Manager()
        input_queue = manager.Queue()
        output_queue = manager.Queue()
        # Rather than passing parsed DataFrames through a pipe for each chunk, append them to a full
        # output table file. Use a lock to prevent multiple processes from writing at once.
        full_output_path = os.path.join(self.temp_dir, 'parsed_output.tsv')
        lock = mp.Lock()
        parsing_processes = [mp.Process(target=parsing_worker,
                                        args=(input_queue, output_queue, self.align_output_length, self.match_mode, full_output_path, lock))
                             for _ in range(self.num_threads)]
        for p in parsing_processes:
            p.start()
        num_unparsed_chunks = 0

        query_file = open(self.fasta_query_path)
        output_dfs = []
        if unprocessed_chunk_dict:
            self.progress.update_pid(pid)
            self.progress.update(f"Queuing queries 1-{total_lines}/{total_lines}")
        while True:
            # Stay in the loop while there are chunks to be processed. Fill up as many available
            # cores as possible with vmatch processes.
            if len(unprocessed_chunk_dict) + num_unparsed_chunks < self.num_threads:
                if next_chunk_start < total_lines:
                    # Write the next chunk of query sequences to a new file.
                    query_chunk_path = os.path.join(self.temp_dir, 'query_chunk_' + str(unprocessed_chunk_num) + '.fa')
                    unprocessed_chunk_dict[unprocessed_chunk_num]['query_path'] = query_chunk_path
                    with open(query_chunk_path, 'w') as query_chunk_file:
                        while running_line_count < next_chunk_start:
                            query_chunk_file.write(next(query_file))
                            running_line_count += 1

                    output_chunk_path = os.path.join(self.temp_dir, 'output_' + str(unprocessed_chunk_num) + '.txt')
                    unprocessed_chunk_dict[unprocessed_chunk_num]['output_path'] = output_chunk_path

                    if match_mode == 'exact_query_substring':
                        command = ('vmatch '
                                   f'-q {query_chunk_path} '
                                   '-complete '
                                   '-identity 100'
                                   '-nodist -noevalue -noidentity -noscore '
                                   '-showdesc 0 '
                                   '-v '
                                   f'{self.index_path}').split(' ')
                    elif match_mode == 'query_substring_with_mismatches':
                        command = ('vmatch '
                                   f'-q {query_chunk_path} '
                                   '-complete '
                                   f'-h {self.max_hamming_dist} '
                                   f'-identity {self.min_ident} '
                                   '-nodist -noevalue -noidentity -noscore '
                                   '-showdesc 0 '
                                   f'-s {self.align_output_length} '
                                   '-v '
                                   f'{self.index_path}').split(' ')
                    output_chunk_file = open(output_chunk_path, 'w', encoding='utf-8')
                    unprocessed_chunk_dict[unprocessed_chunk_num]['output_file'] = output_chunk_file
                    subprocess = utils.start_command(command, self.log_path, stdout=output_chunk_file, remove_log_file_if_exists=False)
                    unprocessed_chunk_dict[unprocessed_chunk_num]['subprocess'] = subprocess

                    self.progress.update_pid(pid)
                    self.progress.update(f"Queuing queries {next_chunk_start - query_chunk_size + 1}-{next_chunk_start}/{total_lines}")
                    # When query sequences are chunked, the maximum value of `next_chunk_start` is
                    # `total_lines`. In the absence of chunking, `next_chunk_start` exceeds
                    # `total_lines`.
                    next_chunk_start = next_chunk_start + query_chunk_size if next_chunk_start + query_chunk_size < total_lines else total_lines
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

                # Clean up.
                if next_chunk_start <= total_lines:
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

            if (next_chunk_start >= total_lines) and (not unprocessed_chunk_dict) and (num_unparsed_chunks == 0):
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
            output_df = pd.read_csv(full_output_path, sep='\t', header=None, names=['query_name', 'target_name', 'query_start_in_target', 'query_length'])
        elif self.match_mode == 'query_substring_with_mismatches':
            output_df = pd.read_csv(full_output_path, sep='\t', header=None, names=['query_name', 'target_name', 'query_start_in_target', 'mismatch_positions'])
        self.progress.end()

        if anvio.DEBUG:
            self.run.warning(f"The temp directory, '{self.temp_dir}', is kept. Don't forget to clean it up later!", header="Debug")
        elif not self.keep_temp_dir:
            self.run.info_single("Cleaning up the temp directory (use `--debug` to keep it for testing purposes)", nl_before=1, nl_after=1)
            shutil.rmtree(self.temp_dir)

        return output_df


def parsing_worker(input_queue, output_queue, align_output_length, match_mode, full_output_path, lock):
    """This worker is called from `Vmatch.search_queries` to parse chunked vmatch output. It is
    located outside `Vmatch` to allow multiprocessing."""
    while True:
        output_path = input_queue.get()
        if match_mode == 'exact_query_substring':
            output_df = parse_exact_query_substrings(output_path)
        elif match_mode == 'query_substring_with_mismatches':
            output_df = parse_query_substrings_with_mismatches(output_path, align_output_length)

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
    """Parse queries that are fully contained by targets and that contain mismatches but not gaps.
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
    # The first line of the alignment is (part of) the subject sequence.
    prev_line_was_subject = False
    # Alignments longer than the determined width occur are split up with a blank line between
    # sections. A second blank line signifies the end of the alignment record.
    prev_line_was_blank_after_align = False
    align_index = 0
    mismatch_positions = []
    for line in output_file:
        if parsing_align:
            if line[0] == 'S':
                prev_line_was_blank_after_align = False
                prev_line_was_subject = True
            elif line[0] == 'Q':
                if prev_line_was_subject:
                    prev_line_was_subject = False
                    # If the query line comes immediately after the subject line, then this section
                    # of the alignment has no mismatches.
                    align_index += align_output_length
            elif prev_line_was_subject:
                # A line between the subject and query is under consideration. Its existence means
                # there is a mismatch in this section of the alignment.
                for section_index, char in enumerate(line[7: 7 + align_output_length]):
                    if char == '!':
                        mismatch_positions.append(str(align_index + section_index))
                prev_line_was_subject = False
                align_index += align_output_length
            elif prev_line_was_blank_after_align:
                # A second blank line has now been encountered, indicating the end of a record.
                if mismatch_positions:
                    # Retain the query name, target name, query start in target, and mismatch positions in the query.
                    align_records.append((align_record[5], align_record[1], align_record[2], ','.join(mismatch_positions)))
                parsing_align = False
                prev_line_was_blank_after_align = False
                align_index = 0
                mismatch_positions = []
                continue
            else:
                prev_line_was_blank_after_align = True
                continue
        else:
            align_record = line.rstrip().split()
            parsing_align = True
    output_file.close()

    return pd.DataFrame(align_records, columns=['query_name', 'target_name', 'query_start_in_target', 'mismatch_positions'])
