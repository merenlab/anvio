# coding: utf-8
"""Interface to HMMer."""

import os
import io
import glob
import shutil

# multiprocess is a fork of multiprocessing that uses the dill serializer instead of pickle
# using the multiprocessing module directly results in a pickling error in Python 3.10 which
# goes like this:
#
#   >>> AttributeError: Can't pickle local object 'SOMEFUNCTION.<locals>.<lambda>' multiprocessing
#
import multiprocess as multiprocessing

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print
P = terminal.pluralize


class HMMer:
    def __init__(self, target_files_dict, num_threads_to_use=1, program_to_use='hmmscan', progress=progress, run=run):
        """A class to streamline HMM runs.

        Notes
        =====
        - HMMer user guide: http://eddylab.org/software/hmmer/Userguide.pdf
        """

        self.num_threads_to_use = num_threads_to_use
        self.program_to_use = program_to_use
        self.progress = progress
        self.run = run

        self.tmp_dirs = []
        self.target_files_dict = {}
        self.number_of_sequences = {}

        acceptable_programs = ["hmmscan", "hmmsearch"]
        if self.program_to_use not in acceptable_programs:
            raise ConfigError("HMMer class here. You are attempting to use the program %s to run HMMs, but we don't recognize it. The currently "
                              "supported programs are: %s" % (self.program_to_use, ", ".join(acceptable_programs)))

        for source in target_files_dict:
            tmp_dir = filesnpaths.get_temp_directory_path()
            self.tmp_dirs.append(tmp_dir)
            
            # create splitted fasta files inside tmp directory
            split_fastas, number_of_sequences = utils.split_fasta(target_files_dict[source],
                                                               parts=self.num_threads_to_use,
                                                               output_dir=tmp_dir,
                                                               return_number_of_sequences=True)
            
            self.target_files_dict[source] = split_fastas
            self.number_of_sequences[source] = number_of_sequences

    def verify_hmmpress_output(self, hmm_path):
        """This function verifies that the HMM profiles located at hmm_path have been successfully hmmpressed.

        What this means is that every .hmm profile in the directory has an associated .h3f, .h3i, .h3m, and
        .h3p file.

        Parameters
        ==========
        hmm_path: string
            the path at which the HMM profiles are located
        """

        for file_path in glob.glob(os.path.join(hmm_path, '*.hmm')):
            base_path = file_path[:-3]
            expected_extensions = ['h3f', 'h3i', 'h3m', 'h3p']
            for ext in expected_extensions:
                if not os.path.exists(base_path + ext):
                    raise ConfigError("It appears that hmmpress was not properly run on the hmm profiles at %s. The "
                                      "file %s does not exist. It is likely that you will have to set up your profiles "
                                      "again by running a program such as `anvi-setup-pfams` or `anvi-setup-kegg-data`. "
                                      "We are very sorry about this." % (hmm_path, base_path + ext))


    def run_hmmer(self, source, alphabet, context, kind, domain, num_genes_in_model, hmm, ref, noise_cutoff_terms,
                  desired_output='table', hmmer_output_dir=None):
        """Run the program

        Parameters
        ==========
        source : str
            A name for your HMM effort.

        alphabet : str
            Which alphabet are you using? Choose from {'AA', 'DNA', 'RNA'}

        context : str
            This will determine how your output is processed. FIXME Documentation is lacking. Choose
            from {'GENE', 'CONTIG', 'DOMAIN'}.

        kind : str
            Used for user stdout info. Don't by afraid to pass None

        domain : str
            Used for user stdout info. Don't by afraid to pass None

        num_genes_in_model : int
            Used for user stdout info. Don't by afraid to pass None

        hmm : str
            Path to the input .hmm file

        ref : int
            Used for user stdout info. Don't by afraid to pass None

        noise_cutoff_terms : str
            Filter out hits with built-in flags. e.g. '--cut_ga'

        desired_output : str OR list, 'table'
            HMMER programs have a couple of outputs. For the standard output (specified by the hmmer
            program flag `-o`), pass 'standard'. For the regular tabular output (specified by the hmmer
            program flag `--tblout`), pass 'table'. For the domain tabular output (specified by the hmmer
            program flag `--domtblout`), pass 'domtable'. If you want to use multiple, pass a tuple like
            ('standard', 'table')

        hmmer_output_dir : str
            The path at which to store the HMMER output files, if desired. After all HMMER workers are
            done and their partial output files have been combined into one (for each type), those combined
            output files will be moved to this location.
        """

        target = ':'.join([alphabet, context])

        if target not in self.target_files_dict:
            raise ConfigError("You have an unknown target :/ Target, which defines an alphabet and context "
                               "to clarify whether the HMM search is supposed to be done using alphabets DNA, "
                               "RNA, or AA sequences, and contexts of GENEs or CONTIGs. Yours is %s, and it "
                               "doesn't work for anvi'o." % target)

        if not self.target_files_dict[target]:
            raise ConfigError("HMMer class does not know about Sequences file for the target %s :/" % target)

        if isinstance(desired_output, str):
            desired_output = (desired_output, )

        for output in desired_output:
            if output not in ['standard', 'table', 'domtable']:
                raise ConfigError("HMMer.run_hmmer :: Unknown desired_output, '%s'" % output)

        if hmmer_output_dir:
            if not os.path.exists(hmmer_output_dir):
                filesnpaths.gen_output_directory(hmmer_output_dir)
            else:
                filesnpaths.is_output_dir_writable(hmmer_output_dir)
                for output in desired_output:
                    file_path = os.path.join(hmmer_output_dir, f"hmm.{output}")
                    if filesnpaths.is_file_exists(file_path, dont_raise=True):
                        raise ConfigError(f"The file {file_path} already exists, and anvi'o does not like to "
                                          "to overwrite things. Please either remove the file or rename your "
                                          "desired output.")

        number_of_sequences = self.number_of_sequences[target]
        
        self.run.warning('', header='HMM Profiling for %s' % source, lc='green')
        self.run.info('Reference', ref if ref else 'unknown')
        self.run.info('Kind', kind if kind else 'unknown')
        self.run.info('Alphabet', alphabet)
        self.run.info('Context', context)
        self.run.info('Domain', domain if domain else 'N/A')
        self.run.info('HMM model path', hmm)
        self.run.info('Number of genes in HMM model', num_genes_in_model or 'unknown')
        self.run.info('Number of sequences in database', number_of_sequences or 'unknown')
        self.run.info('Noise cutoff term(s)', noise_cutoff_terms)
        self.run.info('Number of CPUs will be used for search', self.num_threads_to_use)
        if alphabet in ['DNA', 'RNA']:
            self.run.info('HMMer program used for search', 'nhmmscan')
            if 'domtable' in desired_output:
                raise ConfigError("Oh, dear. Someone (probably a programmer) has requested domain table output from "
                                  f"the run_hmmer() function when the alphabet is {alphabet}. Sadly, this will not "
                                  "work because that alphabet requires the use of `nhmmscan`, which does not have "
                                  "the --domtblout parameter.")
        else:
            self.run.info('HMMer program used for search', self.program_to_use)

        tmp_dir = os.path.dirname(self.target_files_dict[target][0])
        self.run.info('Temporary work dir', tmp_dir)

        # check if all hmmpress files are in the HMM directory
        self.verify_hmmpress_output(hmm)

        workers = []
        manager = multiprocessing.Manager() # this dude holds the shared objects that will be modified by workers
        ret_value_queue = manager.Queue(maxsize=self.num_threads_to_use)
        output_queue = manager.Queue()

        # Holds buffer and write lock for each output
        merged_files_dict = {}
        for output in desired_output:
            merged_files_dict[output] = {'buffer': io.StringIO(), 'lock': manager.Lock()}

        num_parts = len(self.target_files_dict[target])
        cores_per_process = 1
        original_num_threads_requested = None
        if num_parts < self.num_threads_to_use:
            cores_per_process = self.num_threads_to_use // num_parts

            self.run.warning(f"You requested {P('core', self.num_threads_to_use)} but there were only {P('sequence', num_parts)} "
                             f"in the FASTA file for the target '{target}'. Anvi'o will use {P('process', num_parts, sfp='es')} "
                             f"with {P('core', cores_per_process)} instead. And that's that.")

            # if we need to change the number of threads for a SINGLE run, then we need to keep
            # in mind and set the originally reqeusted number of threads. not doing that leads
            # to an extremely tricky bug that is described here thanks to help from Daan Speth:
            # https://github.com/merenlab/anvio/issues/1748
            original_num_threads_requested = self.num_threads_to_use
            self.num_threads_to_use = num_parts

        if alphabet in ['DNA', 'RNA'] and self.program_to_use == 'hmmsearch':
            self.run.warning("You requested to use the program `%s`, but because you are working with %s sequences Anvi'o will use `nhmmscan` instead. "
                             "We hope that is alright." % (self.program_to_use, alphabet))

        thread_num = 0
        for partial_input_file in self.target_files_dict[target]:
            log_file = partial_input_file + '_log'
            output_file = partial_input_file + '_output'
            table_file = partial_input_file + '_table'
            if 'domtable' in desired_output:
                domtable_file = partial_input_file + '_domtable'
            else:
                domtable_file = None

            self.run.info('Log file for thread %s' % thread_num, log_file)
            thread_num += 1

            if noise_cutoff_terms:
                if 'domtable' in desired_output:
                    cmd_line = ['nhmmscan' if alphabet in ['DNA', 'RNA'] else self.program_to_use,
                                '-o', output_file, *noise_cutoff_terms.split(),
                                '--cpu', cores_per_process,
                                '--tblout', table_file,
                                '--domtblout', domtable_file,
                                '-Z', number_of_sequences,
                                '--domZ', number_of_sequences,
                                hmm, partial_input_file]
                else:
                    cmd_line = ['nhmmscan' if alphabet in ['DNA', 'RNA'] else self.program_to_use,
                                '-o', output_file, *noise_cutoff_terms.split(),
                                '--cpu', cores_per_process,
                                '--tblout', table_file,
                                '-Z', number_of_sequences,
                                '--domZ', number_of_sequences,
                                hmm, partial_input_file]
            else: # if we didn't pass any noise cutoff terms, here we don't include them in the command line
                if 'domtable' in desired_output:
                    cmd_line = ['nhmmscan' if alphabet in ['DNA', 'RNA'] else self.program_to_use,
                                '-o', output_file,
                                '--cpu', cores_per_process,
                                '--tblout', table_file,
                                '--domtblout', domtable_file,
                                '-Z', number_of_sequences,
                                '--domZ', number_of_sequences,
                                hmm, partial_input_file]
                else:
                    cmd_line = ['nhmmscan' if alphabet in ['DNA', 'RNA'] else self.program_to_use,
                                '-o', output_file,
                                '--cpu', cores_per_process,
                                '--tblout', table_file,
                                '-Z', number_of_sequences,
                                '--domZ', number_of_sequences,
                                hmm, partial_input_file]

            t = multiprocessing.Process(target=self.hmmer_worker, args=(partial_input_file,
                                                       cmd_line,
                                                       table_file,
                                                       output_file,
                                                       desired_output,
                                                       log_file,
                                                       output_queue,
                                                       ret_value_queue,
                                                       domtable_file))
            t.start()
            workers.append(t)

        self.progress.new('Processing')
        self.progress.update(f'Running {self.program_to_use} in {P("thread", self.num_threads_to_use)}...')

        finished_workers = 0
        while finished_workers < self.num_threads_to_use:
            try:
                ret_value = ret_value_queue.get()

                if isinstance(ret_value, Exception):
                    # If thread returns an exception, we raise it and kill the main thread.
                    raise ret_value

                finished_workers += 1
                if ret_value == 0:
                    if anvio.DEBUG:
                        self.run.info_single(f"{finished_workers} out of {self.num_threads_to_use} have finished")
                else:
                    raise ConfigError("An HMMER worker thread came back with an unexpected return value of {ret_value}. "
                                      "Something is probably wrong, so you should contact a developer for help.")

                # if worker finished successfully we can take its individual output file(s) and append them to the main file(s)
                output_dict = output_queue.get()
                for file_type, file in output_dict.items():
                    main_file_buffer = merged_files_dict[file_type]['buffer']
                    main_file_lock = merged_files_dict[file_type]['lock']
                    worker_file = file
                    if file_type == 'table':
                        append_function = self.append_to_main_table_file
                    elif file_type == 'standard':
                        append_function = self.append_to_main_standard_file
                    elif file_type == 'domtable':
                        append_function = self.append_to_main_table_file

                    append_function(main_file_buffer, worker_file, main_file_lock)

            except KeyboardInterrupt:
                self.run.info_single("HMMER driver received SIGINT, terminating all threads...", nl_before=2)
                break

            except Exception as worker_error:
                # An exception was thrown in one of the threads so we kill all of them
                self.progress.end()
                self.run.warning("An exception was thrown in one of the worker threads (see output below for details).")
                for worker in workers:
                    worker.terminate()
                raise worker_error

        for worker in workers:
            worker.terminate()

        self.progress.end()

        if original_num_threads_requested:
            self.num_threads_to_use = original_num_threads_requested
            self.run.info_single(f'Done with {source} 🎊 (and num threads requested is set back to {self.num_threads_to_use}).', level=0, nl_before=1, nl_after=1, mc="cyan")
        else:
            self.run.info_single(f'Done with {source} 🎊', level=0, nl_before=1, nl_after=1, mc="cyan")

        output_file_paths = []
        for output in desired_output:
            if hmmer_output_dir:
                output_file_path = os.path.join(hmmer_output_dir, f"hmm.{output}")
            else:
                output_file_path = os.path.join(tmp_dir, f"hmm.{output}")

            with open(output_file_path, 'w') as out:
                merged_files_dict[output]['buffer'].seek(0)
                out.write(merged_files_dict[output]['buffer'].read())

            if output == 'table' or output == 'domtable':
                num_raw_hits = filesnpaths.get_num_lines_in_file(output_file_path)
                self.run.info(f'Number of raw hits in {output} file', num_raw_hits, progress=self.progress)
                output_file_path = output_file_path if num_raw_hits else None

            output_file_paths.append(output_file_path)

        # Return output path as string if desired_output is len 1. Else return tuple of output paths
        output = output_file_paths[0] if len(output_file_paths) == 1 else tuple(output_file_paths)

        return output


    def hmmer_worker(self, partial_input_file, cmd_line, table_output_file, standard_output_file, desired_output, log_file,
                     output_queue, ret_value_queue, domtable_output_file=None):

        try:
            # First we run the command
            utils.run_command(cmd_line, log_file)

            if not os.path.exists(table_output_file) or not os.path.exists(standard_output_file) or \
                                 (domtable_output_file and not os.path.exists(domtable_output_file)):
                self.progress.end()
                raise ConfigError("Something went wrong with %s and it failed to generate the expected output :/ Fortunately "
                                  "we have this log file which should clarify the problem: '%s'. Please do not forget to include this "
                                  "file in your question if you were to seek help from the community." % (self.program_to_use, log_file))

            # Then we send the results back to the main thread to be appended to the main files
            output_dict = {}
            for output in desired_output:
                if output == 'table':
                    output_dict['table'] = table_output_file
                elif output == 'standard':
                    output_dict['standard'] = standard_output_file
                elif output == 'domtable':
                    output_dict['domtable'] = domtable_output_file
            output_queue.put(output_dict)

            # return value of 0 to indicate success
            ret_value_queue.put(0)

        except Exception as e:
            # This thread encountered an error. We send the error back to the main thread which
            # will terminate the job.
            ret_value_queue.put(e)


    def append_to_main_standard_file(self, merged_file_buffer, standard_output_file, buffer_write_lock):
        """Append standard output to the main file.

        Notes
        =====
        - The resulting file may not be universally parseable because there will be as many headers
          as there are threads, whereas in a proper output file there is only one header. A header
          looks like this, for example (note it ends in a \n character):

          >>> # hmmsearch :: search profile(s) against a sequence database
          >>> # HMMER 3.2.1 (June 2018); http://hmmer.org/
          >>> # Copyright (C) 2018 Howard Hughes Medical Institute.
          >>> # Freely distributed under the BSD open source license.
          >>> # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          >>> # query HMM file:                  /Users/evan/Software/anvio/anvio/data/misc/Interacdome/Pfam-A.hmm
          >>> # target sequence database:        /var/folders/58/mpjnklbs5ql_y2rsgn0cwwnh0000gn/T/tmpsvyamen6/AA_gene_sequences.fa.3
          >>> # output directed to file:         /var/folders/58/mpjnklbs5ql_y2rsgn0cwwnh0000gn/T/tmpsvyamen6/AA_gene_sequences.fa.3_output
          >>> # per-dom hits tabular output:     /var/folders/58/mpjnklbs5ql_y2rsgn0cwwnh0000gn/T/tmpsvyamen6/AA_gene_sequences.fa.3_table
          >>> # model-specific thresholding:     GA cutoffs
          >>> # number of worker threads:        1
          >>> # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          >>>

          Additionally, there will be as many "[ok]"'s as there are threads, whereas in a proper
          output file there is only one that marks the EOF.
        """

        with open(standard_output_file, 'r') as f:
            with buffer_write_lock:
                merged_file_buffer.write(f.read())


    def append_to_main_table_file(self, merged_file_buffer, table_output_file, buffer_write_lock, lines_per_chunk=1000000):
        """Append table output to the main file.

        Lines starting with '#' (i.e., header lines) are ignored.
        """

        detected_non_ascii = False
        lines_with_non_ascii = []
        output_lines = []

        def process_batch(line_batch, base_line_number):
            nonlocal detected_non_ascii
            local_output = []
            for i, line_bytes in enumerate(line_batch):
                line_number_actual = base_line_number + i
                try:
                    line = line_bytes.decode('ascii')
                except UnicodeDecodeError:
                    decoded = line_bytes.decode('ascii', 'ignore')
                    if len(decoded) != len(line_bytes):
                        detected_non_ascii = True
                        lines_with_non_ascii.append(line_number_actual)
                    line = decoded
                if not line.startswith('#'):
                    local_output.append(line.rstrip())
            return local_output

        with open(table_output_file, 'rb') as hmm_hits_file:
            batch = []
            base_line_number = 1

            for line_bytes in hmm_hits_file:
                batch.append(line_bytes)
                if len(batch) >= lines_per_chunk:
                    output_lines.extend(process_batch(batch, base_line_number))
                    base_line_number += len(batch)
                    batch = []

            if batch:
                filtered_lines = process_batch(batch, base_line_number)
                if filtered_lines:
                    filtered_text = "\n".join(filtered_lines) + "\n"
                    with buffer_write_lock:
                        merged_file_buffer.write(filtered_text)

        # Log warning if non-ASCII characters were detected
        if detected_non_ascii:
            self.run.warning(
                "Just a heads-up, Anvi'o HMMer parser detected non-ASCII characters while processing "
                f"the file '{table_output_file}' and cleared them. Here are the line numbers with non-ASCII characters: "
                f"{', '.join(map(str, lines_with_non_ascii))}. You may want to check those lines with a command like "
                "\"awk 'NR==<line number>' <file path> | cat -vte\"."
            )



    def clean_tmp_dirs(self):
        for tmp_dir in self.tmp_dirs:
            shutil.rmtree(tmp_dir)
