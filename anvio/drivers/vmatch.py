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
        # Vmatch also can run in all-against-all search mode, but this is not supported yet by the driver.
        self.supported_search_modes = ('query', )
        # In the match mode, "exact_query_substring", the query names should end in
        # "-<query length>", e.g., "c_000000000001-75".
        self.supported_match_modes = ('exact_query_substring', )

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.search_mode = A('search_mode')
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

        self.min_align_length = A('min_align_length')
        self.align_seed_length = A('align_seed_length')
        # When extending a seed in both directions, allowing both mismatches and indels, this is the
        # maximum drop in alignment score allowed before extension stops. Mismatches are scored as
        # -1 and indels as -2.
        self.exdrop = A('exdrop')
        self.min_ident = A('min_ident')

        self.progress_interval = A('progress_interval')
        self.quiet = A('quiet')

        self.temp_dir = filesnpaths.get_temp_directory_path()
        self.log_path = os.path.join(self.temp_dir, '00_log.txt')

        self.run = run or terminal.Run(verbose=(not self.quiet))
        self.progress = progress or terminal.Progress(verbose=(not self.quiet))

        if not skip_sanity_check:
            self.sanity_check()


    def sanity_check(self):
        """Do basic checks before proceeding."""
        self.check_programs()

        for name, variable in (("the search mode", self.search_mode),
                               ("the match mode", self.match_mode),
                               ("the target sequence FASTA file path", self.fasta_db_path),
                               ("the minimum alignment length", self.min_align_length),
                               ("the alignment seed length", self.align_seed_length),
                               ("the exdrop", self.exdrop),
                               ("the minimum identity", self.min_ident)):
            if not variable:
                raise ConfigError(f"A proper instance of the anvi'o Vmatch driver must have {name} variable set.")

        if self.search_mode not in self.supported_search_modes:
            raise ConfigError(f"The supported Vmatch search modes are: {', '.join(self.supported_search_modes)}. "
                              f"No mode named {self.search_mode} is allowed.")

        if self.match_mode not in self.supported_match_modes:
            raise ConfigError(f"The supported Vmatch match modes are: {', '.join(self.supported_match_modes)}. "
                              f"No mode named {self.match_mode} is allowed.")

        filesnpaths.is_file_exists(self.fasta_db_path)
        self.index_path = os.path.join(self.temp_dir, 'index')

        if self.min_align_length < self.align_seed_length:
            raise ConfigError("The seed from which the sequence alignment is found cannot be longer than the minimum alignment length. "
                              f"The provided seed length was {self.align_seed_length} and the minimum alignment length was {self.min_align_length}.")

        if self.exdrop < 1:
            raise ConfigError("Exdrop must have an integer value of at least 1. "
                              f"The provided value was {self.exdrop}.")

        if self.min_ident % 1 != 0:
            raise ConfigError("The minimum identity of an alignment must be an integer between 1 and 100. "
                              f"The provided value was {self.min_ident}.")

        # Check variables specific to the "query" search mode.
        if self.search_mode == 'query':
            for name, variable in (("the query sequence FASTA file path", self.fasta_db_path), ):
                if not variable:
                    raise ConfigError(f"A proper instance of the anvi'o Vmatch driver in query search mode must have {name} variable set.")


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

        command = ['mkvtree',
                   '-db', self.fasta_db_path,
                   '-dna',
                   '-indexname', self.index_path,
                   '-pl', # sequence suffixes are placed into buckets by their prefix length, here automatically determined
                   '-tis', '-suf', '-sti1', '-lcp', '-bck', # index tables, excluding some optional tables -- all tables can be produced with "-allout" instead
                   '-v']

        ret_val = utils.run_command(command, self.log_path)

        if int(ret_val):
            raise ConfigError("vmktree returned with non-zero exit code, there may be some errors. "
                              f"Please check the log file for details: {self.log_path}")
        self.progress.end()


    def search_queries(self, make_index=True):
        if make_index:
            self.make_index()

        self.progress.new("Vmatch")
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

        if next_chunk_start >= total_lines:
            # All query sequences fit in one chunk, so use the original query file as vmatch input.
            query_chunk_path = self.fasta_query_path
            unprocessed_chunk_dict[unprocessed_chunk_num]['query_path'] = query_chunk_path
            output_chunk_path = 'output.txt'
            unprocessed_chunk_dict[unprocessed_chunk_num]['output_path'] = output_chunk_path

            command = ['vmatch',
                       '-q', query_chunk_path,
                       '-l', self.min_align_length,
                       '-seedlength', self.align_seed_length,
                       '-exdrop', self.exdrop,
                       '-identity', self.min_ident,
                       '-noevalue', '-noidentity', '-noscore', # do not include alignment E-value, identity and score in the output, but include the alignment distance
                       '-showdesc', 0, # query and target sequence names included in the output
                       '-v', # verbose commented output
                       self.index_path]
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
                                        args=(input_queue, output_queue, self.match_mode, full_output_path, lock))
                             for _ in range(self.num_threads)]
        for p in parsing_processes:
            p.start()
        num_unparsed_chunks = 0

        query_file = open(self.fasta_query_path)
        output_dfs = []
        if unprocessed_chunk_dict:
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

                    command = ['vmatch',
                               '-q', query_chunk_path,
                               '-l', self.min_align_length,
                               '-seedlength', self.align_seed_length,
                               '-exdrop', self.exdrop,
                               '-identity', self.min_ident,
                               '-noevalue', '-noidentity', '-noscore', # do not include alignment E-value, identity and score in the output, but include the alignment distance
                               '-showdesc', 0, # query and target sequence names included in the output
                               '-v', # verbose commented output
                               self.index_path]
                    output_chunk_file = open(output_chunk_path, 'w', encoding='utf-8')
                    unprocessed_chunk_dict[unprocessed_chunk_num]['output_file'] = output_chunk_file
                    subprocess = utils.start_command(command, self.log_path, stdout=output_chunk_file, remove_log_file_if_exists=False)
                    unprocessed_chunk_dict[unprocessed_chunk_num]['subprocess'] = subprocess

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

        # Load the concatenated chunked output and clip query lengths from query names.
        self.progress.update("Finalizing matches")
        output_df = pd.read_csv(full_output_path, sep='\t', header=None, names=['query_name', 'target_name'])
        output_df.loc[:, 'query_name'] = output_df['query_name'].apply(lambda name: '_'.join(name.split('_')[: -1]))

        if anvio.DEBUG:
            self.run.warning(f"The temp directory, '{self.temp_dir}', is kept. Don't forget to clean it up later!", header="Debug")
        else:
            self.run.info_single("Cleaning up the temp directory (use `--debug` to keep it for testing purposes)", nl_before=1, nl_after=1)
            shutil.rmtree(self.temp_dir)

        self.progress.end()
        return output_df


def parsing_worker(input_queue, output_queue, match_mode, full_output_path, lock):
    """This worker is called from `Vmatch.search_queries` to parse chunked vmatch output. It is
    located outside `Vmatch` to allow multiprocessing."""
    while True:
        output_path = input_queue.get()
        output_df = pd.read_csv(output_path,
                                delim_whitespace=True, # there is a variable number of spaces between fields, and there are spaces at the beginning of each line
                                comment='#', # remove comment lines at the beginning of the file
                                header=None,
                                usecols=[0, 1, 2, 4, 5, 6, 7], # disregard the fourth column, indicating a direct match, "D", or palindromic match, "P", as the latter does not occur here
                                names=['target_align_length', 'target_name', 'target_start', 'query_align_length', 'query_name', 'query_start', 'distance'])
        os.remove(output_path)
        if match_mode == 'exact_query_substring':
            output_df = parse_exact_query_substrings(output_df)
        with lock:
            output_df.to_csv(full_output_path, mode='a', sep='\t', index=False, header=False)
        output_queue.put(True)


def parse_exact_query_substrings(output_df):
    """Parse queries that are exact substrings of one or more targets. The closest option that
    Vmatch provides, `-complete`, only retains queries that are the same as targets."""
    t = time.time()
    output_df.loc[:, 'query_length'] = output_df['query_name'].apply(lambda name: int(name.split('-')[-1]))
    output_df.loc[:, 'unalign_length'] = output_df['query_length'] - output_df['query_align_length']

    # A full match is indicated by an alignment spanning the entire query with an edit distance
    # of zero.
    output_df = output_df[(output_df['unalign_length'] == 0) & (output_df['distance'] == 0)].copy()

    # Only the query and target names are needed.
    return output_df[['query_name', 'target_name']]
