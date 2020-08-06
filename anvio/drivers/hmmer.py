# coding: utf-8
"""Interface to HMMer."""

import os
import io
import gzip
import shutil
import glob
from threading import Thread, Lock


import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


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

        acceptable_programs = ["hmmscan", "hmmsearch"]
        if self.program_to_use not in acceptable_programs:
            raise ConfigError("HMMer class here. You are attempting to use the program %s to run HMMs, but we don't recognize it. The currently "
                              "supported programs are: %s" % (self.program_to_use, ", ".join(acceptable_programs)))

        for source in target_files_dict:
            tmp_dir = filesnpaths.get_temp_directory_path()
            self.tmp_dirs.append(tmp_dir)

            part_file_name = os.path.join(tmp_dir, os.path.basename(target_files_dict[source]))

            # create splitted fasta files inside tmp directory
            self.target_files_dict[source] = utils.split_fasta(target_files_dict[source],
                                                               parts=self.num_threads_to_use,
                                                               prefix=part_file_name)

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
                                      "again by running a program such as `anvi-setup-pfams` or `anvi-setup-kegg-kofams`. "
                                      "We are very sorry about this." % (hmm_path, base_path + ext))


    def run_hmmer(self, source, alphabet, context, kind, domain, num_genes_in_model, hmm, ref, noise_cutoff_terms,
                  desired_output='table', out_fmt='--tblout'):
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
            program flag `-o`), pass 'standard'. For the tabular output (specified by the hmmer
            program flag `--tblout` or `--domtblout`), pass 'table'. If you want to use both, pass
            ('standard', 'table')

        out_fmt : str, '--tblout'
            HMMer programs have different table output formats. For example, choose from --tblout or
            --domtblout.
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
            if output not in ['standard', 'table']:
                raise ConfigError("HMMer.run_hmmer :: Unknown desired_output, '%s'" % output)

        if out_fmt not in ['--tblout', '--domtblout']:
            raise ConfigError("HMMer.run_hmmer :: Unknown out_fmt, '%s'" % out_fmt)

        self.run.warning('', header='HMM Profiling for %s' % source, lc='green')
        self.run.info('Reference', ref if ref else 'unknown')
        self.run.info('Kind', kind if kind else 'unknown')
        self.run.info('Alphabet', alphabet)
        self.run.info('Context', context)
        self.run.info('Domain', domain if domain else 'N/A')
        self.run.info('HMM model path', hmm)
        self.run.info('Number of genes in HMM model', num_genes_in_model or 'unknown')
        self.run.info('Noise cutoff term(s)', noise_cutoff_terms)
        self.run.info('Number of CPUs will be used for search', self.num_threads_to_use)
        if alphabet in ['DNA', 'RNA']:
            self.run.info('HMMer program used for search', 'nhmmscan')
        else:
            self.run.info('HMMer program used for search', self.program_to_use)

        tmp_dir = os.path.dirname(self.target_files_dict[target][0])
        self.run.info('Temporary work dir', tmp_dir)

        # check if all hmmpress files are in the HMM directory
        self.verify_hmmpress_output(hmm)

        workers = []

        # Holds buffer and write lock for each output
        merged_files_dict = {}
        for output in desired_output:
            merged_files_dict[output] = {'buffer': io.StringIO(), 'lock': Lock()}

        num_parts = len(self.target_files_dict[target])
        cores_per_process = 1
        if num_parts < self.num_threads_to_use:
            cores_per_process = self.num_threads_to_use // num_parts

            self.run.warning("You requested %s cores but there were only %s entries in the fasta for the target '%s'. "
                             "Anvi'o will use %s process with %s cores each instead. I hope thats okay for you. " %
                             (str(self.num_threads_to_use), str(num_parts), target, str(num_parts), cores_per_process))

        if alphabet in ['DNA', 'RNA'] and self.program_to_use == 'hmmsearch':
            self.run.warning("You requested to use the program `%s`, but because you are working with %s sequences Anvi'o will use `nhmmscan` instead. "
                             "We hope that is alright." % (self.program_to_use, alphabet))

        thread_num = 0
        for partial_input_file in self.target_files_dict[target]:
            log_file = partial_input_file + '_log'
            output_file = partial_input_file + '_output'
            table_file = partial_input_file + '_table'

            self.run.info('Log file for thread %s' % thread_num, log_file)
            thread_num += 1

            if noise_cutoff_terms:
                cmd_line = ['nhmmscan' if alphabet in ['DNA', 'RNA'] else self.program_to_use,
                            '-o', output_file, *noise_cutoff_terms.split(),
                            '--cpu', cores_per_process,
                            out_fmt, table_file,
                            hmm, partial_input_file]
            else: # if we didn't pass any noise cutoff terms, here we don't include them in the command line
                cmd_line = ['nhmmscan' if alphabet in ['DNA', 'RNA'] else self.program_to_use,
                            '-o', output_file,
                            '--cpu', cores_per_process,
                            out_fmt, table_file,
                            hmm, partial_input_file]

            t = Thread(target=self.hmmer_worker, args=(partial_input_file,
                                                       cmd_line,
                                                       table_file,
                                                       output_file,
                                                       desired_output,
                                                       log_file,
                                                       merged_files_dict))
            t.start()
            workers.append(t)

        self.progress.new('Processing')
        self.progress.update('Running %s in %d threads...' % (self.program_to_use, self.num_threads_to_use))

        # Wait for all workers to finish.
        for worker in workers:
            worker.join()

        self.progress.end()

        output_file_paths = []
        for output in desired_output:
            output_file_path = os.path.join(tmp_dir, f"hmm.{output}")

            with open(output_file_path, 'w') as out:
                merged_files_dict[output]['buffer'].seek(0)
                out.write(merged_files_dict[output]['buffer'].read())

            if output == 'table':
                num_raw_hits = filesnpaths.get_num_lines_in_file(output_file_path)
                self.run.info('Number of raw hits', num_raw_hits)
                output_file_path = output_file_path if num_raw_hits else None

            output_file_paths.append(output_file_path)

        self.progress.end()

        # Return output path as string if desired_output is len 1. Else return tuple of output paths
        output = output_file_paths[0] if len(output_file_paths) == 1 else tuple(output_file_paths)
        return output


    def hmmer_worker(self, partial_input_file, cmd_line, table_output_file, standard_output_file, desired_output, log_file,
                     merged_files_dict):

        # First we run the command
        utils.run_command(cmd_line, log_file)

        if not os.path.exists(table_output_file) or not os.path.exists(standard_output_file):
            self.progress.end()
            raise ConfigError("Something went wrong with %s and it failed to generate the expected output :/ Fortunately "
                              "we have this log file which should clarify the problem: '%s'. Please do not forget to include this "
                              "file in your question if you were to seek help from the community." % (self.program_to_use, log_file))

        # Then we append the results to the main file(s)
        for output in desired_output:
            main_file_buffer = merged_files_dict[output]['buffer']
            main_file_lock = merged_files_dict[output]['lock']

            if output == 'table':
                worker_file = table_output_file
                append_function = self.append_to_main_table_file
            elif output == 'standard':
                worker_file = standard_output_file
                append_function = self.append_to_main_standard_file

            append_function(main_file_buffer, worker_file, main_file_lock)


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


    def append_to_main_table_file(self, merged_file_buffer, table_output_file, buffer_write_lock):
        """Append table output to the main file.

        FIXME In addition to appending, this function also pre-processes the data, which should not
        be done here. That qualifies as hmmer output parsing, and should be in
        anvio/parsers/hmmer.py.
        """

        detected_non_ascii = False
        lines_with_non_ascii = []
        clip_description_index = None
        clip_index_found = False

        with open(table_output_file, 'rb') as hmm_hits_file:
            line_counter = 0
            for line_bytes in hmm_hits_file:
                line_counter += 1
                line = line_bytes.decode('ascii', 'ignore')

                if not len(line) == len(line_bytes):
                    lines_with_non_ascii.append(line_counter)
                    detected_non_ascii = True

                if line.startswith('#'):
                    if not clip_index_found and line.find('description') != -1:
                        # This parser removes the description column from the data
                        clip_description_index = line.find('description')
                        clip_index_found = True

                    continue

                with buffer_write_lock:
                    merged_file_buffer.write('\t'.join(line[:clip_description_index].split()) + '\n')

        if detected_non_ascii:
            self.run.warning("Just a heads-up, Anvi'o HMMer parser detected non-ascii characters while processing "
                             "the file '%s' and cleared them. Here are the line numbers with non-ascii characters: %s. "
                             "You may want to check those lines with a command like \"awk 'NR==<line number>' <file path> | cat -vte\"." %
                                                 (table_output_file, ", ".join(map(str, lines_with_non_ascii))))


    def clean_tmp_dirs(self):
        for tmp_dir in self.tmp_dirs:
            shutil.rmtree(tmp_dir)
