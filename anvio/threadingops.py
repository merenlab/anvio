# coding: utf-8
"""A namespace for various threading operations.

Exports
=======

- State(UserDict): Light wrapper around UserDict
- AnviThread(Thread): Light wrapper around Thread
- ThreadedCommandRunner(abc.ABC): Abstract class/interface for multithreaded command runners
"""

# todo check args of everything!

import os
import abc
import anvio

from anvio import fastalib
from anvio import filesnpaths

from anvio.errors import ConfigError, CommandError
from collections import UserDict
from threading import Thread
from anvio.utils.commandline import run_command
from anvio.utils.fasta import split_fasta
from anvio.utils.files import concatenate_files


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Ryan Moore"
__email__ = "moorer@udel.edu"


class State(UserDict):
    """Light wrapper around UserDict.

    Raises ValueError if the key isn't a str.

    Useful for typing.
    """

    def __setitem__(self, key, value):
        if isinstance(key, str):
            super().__setitem__(key, value)
        else:
            raise ValueError(f'key must be str.  Got {type(key).__name__}.')


# Todo: add type hints to this.
class AnviThread(Thread):
    """Like regular Thread, except that it tracks the return value of the `target` function call.

    You can then use this return value so the caller can decide what to do.  For example, if you're running a bash
    command, you may want to retry it if the return value is > 0, (i.e., indicating failure).
    """

    def __init__(self, group=None, target=None, name=None, args=(), kwargs=None, *, daemon=None):

        super().__init__(group, target, name, args, kwargs, daemon=daemon)

        # Instance variable to track return value of target function.
        self.target_return_value = None

    def run(self):
        """Method representing the thread's activity.

        You may override this method in a subclass. The standard run() method
        invokes the callable object passed to the object's constructor as the
        target argument, if any, with sequential and keyword arguments taken
        from the args and kwargs arguments, respectively.

        """
        try:
            if self._target:
                # Save the output of the target function in an instance variable for later use.
                self.target_return_value = self._target(*self._args, **self._kwargs)
        finally:
            # Avoid a refcycle if the thread is running a function with
            # an argument that has a member that points to the thread.
            del self._target, self._args, self._kwargs


# todo fix all the self.logger checks for None


class ThreadedCommandRunner(abc.ABC):
    """An abstract class for making mutlithreaded command runners.

    Overview
    ========

    It abstracts the process of making embarrassingly parallel (
    https://en.wikipedia.org/wiki/Embarrassingly_parallel) types of command runners.  E.g., splitting an input file,
    running a command on all the splits, and then collating the results.

    This class provides default implementations of certain methods, but feel free to override them to meet your needs!

    Some methods are abstract and must be given a concrete implementation in all subclasses.  Abstract methods that you
    must implement include:
        - _split_input_file: method for splitting the input file
        - _make_commands: method for making the actual commands to run
        - _collate_output_files: method for collating any output files

    See `ThreadedProdigalRunner` for an example of a concrete implementation of this interface.

    See the `run` method for how work gets done.  One thing to note is that any method that is called by the `run`
    method should return `State`.  This can be empty (but shouldn't be None).  This provides a convenient way to get
    any state out of the runner object into the caller.  You could then use this info for logging, debugging,
    returning info from parsing output files, etc.

    See individual method docstrings for guidelines on how to implement abstract methods or override existing
    concrete implementations.

    Things to keep in mind:
        - This is an abstact class (with some default implementations), so it's not meant to be used directly.  You
          should subclass it for each particular command you want to run!  See
          `anvio.drivers.prodigal.MultithreadedProdigalRunner` for an example.
        - For now, it only handles the case where a single input file needs to be spilt.
        - If you only want a single thread/split, it will go through the whole process of "splitting" the input file
          into a single split, running the command in a thread, and collating the single output.  It sounds wasteful,
          but in my tests, it generally doesn't make that much difference in the total runtime, as you would only be
          using this class for commands that take a long time to run anyway.  (Or else why bother splitting them up?)

    Parameters
    ==========

    See __init__ definition for expected types!

    - `input_file_path`: Path to file to split and run the commands on.
    - `collated_output_file_paths`: Mapping of output file name (key) to output file path (value).  The key is for
       referencing in the code, and won't be reflected in the filesystem.  Some software you want to run will output
       more than one important file, so The mapping allows you to specify the location of multiple outfiles.
       E.g., {'peptide_path': '/path/to/peptides.faa', 'gff_path': '/path/to/file.gff'}.
    - `log_file_path`: Path to the logfile (e.g., '/path/to/logs.txt')
    - `number_of_splits`: How many splits do you want?  Keep it reasonable: a `Thread` will be spawned for each split.
    - `logger`: a `anvio.terminal.Logger`.  Let's the commands access anvio's nice logging.
    """

    def __init__(self,
                 input_file_path,
                 collated_output_file_paths,
                 log_file_path,
                 number_of_splits,
                 logger,
                 ):
        """In addition to the named parameters, the following instance variables are also initialized to empty values:

        - threads
        - commands
        - input_file_splits
        - output_file_split_paths
        """
        self.input_file_path = input_file_path
        self.collated_output_file_paths = collated_output_file_paths
        self.number_of_splits = number_of_splits
        self.log_file_path = log_file_path
        self.logger = logger

        self.threads = []
        self.commands = []
        self.input_file_splits = []
        self.output_file_split_paths = {}

    # Main entrypoint for this class
    def run(self):
        """Run the multithreaded command wrapper, and return the `State`.

        Each of the intermediate function calls should return `State`.  Some of the steps may return empty `State`
        mappings, this is fine!  Note that latter steps may overwrite keys from earlier steps.

        State is returned to make it easier for the calling context to access any important state generated in any step
        of the pipeline.
        """
        # Here comes what looks like a clever hack, but is really just python dictionary merging.
        # See https://www.python.org/dev/peps/pep-0448/.  We use the {} literal instead of dict() in case one of the
        # later functions needs to override the value in one of the earlier.
        state = State({**self._split_input_file(),
                       **self._set_output_file_split_paths(),
                       **self._make_commands(),
                       **self._run_commands(),
                       **self._collate_output_files(),
                       **self._remove_input_file_splits(),
                       **self._remove_output_file_splits()})

        return state

    # Feel free to redefine this one as well!  You'll need to handle these errors because the Thread won't bubble errors
    # to the calling context.
    #
    # Needs its own special logfile because each thread wants to write to it.
    # the collection will most likely hold int, str, float, that sort of thing
    @staticmethod
    def _command_runner(command, log_file_path):
        """Run `command`, writing any logs to `log_file_path`.

        If the command returns a zero exit code, _command_runner returns 0, otherwise it returns `CommandError`.

        Feel free to override this function.  However, if you do, I would suggest that you return CommandError rather
        than raising the error.  This is because _command_runner will be run in its own Thread.  Any Exception that
        is raised in a Thread (or AnviThread) will not bubble up into the calling context.  What that means is that
        if you get an error in your job that's running in the thread, it will just crash the thread but the calling
        context will happily go on with it's work.  This is often NOT what you want.  So to handle errors in the
        calling context, you need to return something like CommandError, and then handle that.  See the
        implementation for _run_commands for an example of how to do this properly.
        """
        try:
            return_value = run_command(command, log_file_path)
        except ConfigError as e:
            # run_command can raise ConfigError.  So pass the message from that to a CommandError to keep it
            # consistent.
            return CommandError(e.e)

        if return_value < 0 or return_value > 0:
            # Technically, run_command will return ConfigError if the return code was < 0, but just do this
            # sanity check here as well to be sure.
            return CommandError(f"Failed to run '{command}'.  Exit code: {return_value}")
        else:
            return return_value

    # this method should set self.input_file_splits
    @abc.abstractmethod
    def _split_input_file(self):
        """Takes `input_file_path` pointing to the file you want to split up.

        It splits this file into a certain `number_of_splits`, and writes the splits to disk.

        Returns `State` with any info needed by the caller.

        IMPORTANT:  Unless you're planning on overriding all the methods, when implementing `_split_input_file` you
        must set self.input_file_splits to track the paths of the splits you wrote to disk.
        """
        pass

    def _set_output_file_split_paths(self):
        """Sets the output_file_split_paths instance variable.

        Keys of `self.output_file_split_paths` will match keys in `self.collated_output_file_paths`.

        Output file split paths are based on the collated output file paths set in the constructor
        (`collated_output_file_paths`).  For example, if you passed in
        `collated_output_file_paths={'gff_path': 'silly.txt'}`, then you would get something like
        `self.collated_output_file_paths['gff_path'] = [...paths to all the splits...]`.

        If you override this function, be sure to set self.output_file_split_paths, and have the keys match those in
        self.collated_output_file_paths, unless you want to override everything.
        """
        for which, path in self.collated_output_file_paths.items():
            paths = ['.'.join([path, f'split_{i}']) for i in range(self.number_of_splits)]
            self.output_file_split_paths[which] = paths

        return State(output_file_split_paths=self.output_file_split_paths)

    @abc.abstractmethod
    def _make_commands(self):
        """Build commands to be run and set them to `self.commands`.

        Should return State, e.g., to return the commands in `self.commands` to the calling context.
        """
        pass

    def _run_commands(self):
        """Run commands built with `_make_commands` each in its own `AnviThread`.

        Uses `run_command` to run all the commands, each one in its own `AnviThread`.

        Will block until all threads are finished.

        Checks for CommandErrors that occurred in individual threads.

        Returns `State` with key `threads` and values being the actual `AnviThreads` created.

        If you override this function, be sure to set `self.threads` with the AnviThreads you created to run the jobs.
        """
        intermediate_log_file_paths = []
        for index, command in enumerate(self.commands):
            # Each thread will get its own logfile so we don't have to worry about locking on the main logfile.
            this_log_file_path = '.'.join([self.log_file_path, f'thread-{index}'])
            intermediate_log_file_paths.append(this_log_file_path)

            # Spawn a new thread and start it.
            thread = AnviThread(target=self._command_runner, args=(command, this_log_file_path))
            thread.start()
            self.threads.append(thread)

        # Wait for all jobs to finish.
        for thread in self.threads:
            thread.join()

        self._check_threads_for_errors()

        # Todo pretty sure this will overwrite logfile.  Is this behavior what we really want?
        anvio.concatenate_files(self.log_file_path, intermediate_log_file_paths, remove_concatenated_files=True)

        return State(threads=self.threads)

    @abc.abstractmethod
    def _collate_output_files(self):
        """Collate all the intermediate output files from each of the threads.

        The result will be that all output files specified in `self.collated_output_file_paths` will contain the results
        from each individual thread.

        Should return State with any important info required from the output files.
        """
        pass

    def _remove_input_file_splits(self):
        """Remove the input file splits specified in `self.input_file_splits`."""
        self._remove_files(self.input_file_splits)

        return State()

    def _remove_output_file_splits(self):
        """Remove all the intermediate output file splits."""
        for _, paths in self.output_file_split_paths.items():
            self._remove_files(paths)

        return State()

    # Helper functions
    def _check_threads_for_errors(self):
        """Check threads for any exceptions that may have occurred during their time running commands.

        Any errors will be raised.
        """
        # Check threads for exceptions
        for thread in self.threads:
            if isinstance(thread.target_return_value, CommandError):
                if self.logger:
                    self.logger.progress.end()
                raise thread.target_return_value

    @staticmethod
    def _remove_files(paths) -> None:
        """Helper function to remove all the given `paths`."""
        for path in paths:
            if filesnpaths.is_file_exists(path, dont_raise=True):
                os.remove(path)


class ThreadedProdigalRunner(ThreadedCommandRunner):
    def __init__(self, args):
        """See ThreadedCommandRunner for more info.

        `args` should have the `__dict__` attribute, e.g., `argparse.Namespace`.

        Parameters not shared with Superclass
        =====================================

        - `installed_version`: The version string of Prodigal
        - `parser`:  The Callable used to parse Prodigal peptide file deflines.  (See signature for expected types.)
        - `translation_table`: Either a string or int specifying the specific translation table, or None.  If None, then
           Prodigal's `-p meta` option will be used instead.
        """
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None

        required_args = ['input_file_path', 'collated_output_file_paths', 'number_of_splits', 'log_file_path',
                         'installed_version', 'parser', 'prodigal_single_mode']

        # Check that the required arguments are present.
        for arg in required_args:
            if A(arg) is None:
                raise KeyError(f"'{arg}' is a required argument")

        super().__init__(input_file_path=A('input_file_path'),
                         collated_output_file_paths=A('collated_output_file_paths'),
                         number_of_splits=A('number_of_splits'),
                         log_file_path=A('log_file_path'),
                         logger=A('logger'))

        self.prodigal_single_mode = A('prodigal_single_mode')
        self.installed_version = A('installed_version')
        self.parser = A('parser')

        # if it's an int we need to cast to str
        translation_table = A('translation_table')
        if isinstance(translation_table, int):
            translation_table = str(translation_table)

        self.translation_table = translation_table


    # Implement the abstract methods
    #
    #
    def _split_input_file(self):
        """Split input fasta into the correct number of splits.

        Returns `State` with the paths to each of the fasta splits.

        Raises `ValueError` if `self.number_of_splits < 0.

        See Superclass `_split_input_file` for more info.
        """
        # Todo: probably should move this check into the constructor.
        if self.number_of_splits <= 0:
            raise ValueError(f'number_of_splits must be > 0.  Got {self.number_of_splits}')

        # Todo are there errors to catch here?
        self.input_file_splits = split_fasta(self.input_file_path, parts=self.number_of_splits, shuffle=True)

        return State(input_file_splits=self.input_file_splits)


    def _make_commands(self):
        """Make commands and store them in `self.commands`.

        Reads input file split paths from `self.input_file_splits`.  Reads output file split paths from
        `self.output_file_split_paths`.  Ensures that they conform.

        Returns commands in `State`.

        See Superclass `_make_commands` for more info.
        """
        # Todo: This should probably be done in the constructor of the abstract superclass.
        # Sanity check.
        if len(self.output_file_split_paths) == 0:
            raise ValueError('Did you forget to initialize `set_output_file_split_paths`?')

        # Check that the proper keyword arguments are given.
        for key in ['peptide_path', 'gff_path']:
            if key not in self.output_file_split_paths:
                raise ValueError(f"self.output_file_split_paths should have the key '{key}', but it is missing.")

            if not len(self.input_file_splits) == len(self.output_file_split_paths[key]):
                raise ConfigError("The number input files do not match to the number of files expected :/")

        for i in range(len(self.input_file_splits)):
            command = ['prodigal',
                       '-m', # <- this one is to treat runs of N as masked sequences and not to build genes across them.
                             #    see a longer discussion here: https://github.com/merenlab/anvio/issues/1641
                       '-i', self.input_file_splits[i],
                       '-o', self.output_file_split_paths['gff_path'][i],
                       '-a', self.output_file_split_paths['peptide_path'][i]]

            # Use either custom translation table or 'meta' mode.
            if self.translation_table:
                # Translation tables and meta mode cannot be used together.
                command.extend(['-g', self.translation_table])
                if self.logger:
                    self.logger.run.warning(
                        "Prodigal translation table is set to '%s' (whatever you did has worked so far, but "
                        "keep an eye for errors from prodigal in case it doesn't like your translation table "
                        "parameter). This means we will not use prodigal in the metagenomics mode, due to this"
                        "issue: https://github.com/hyattpd/Prodigal/issues/19. If that issue is closed, and "
                        "you are reading this message, then please contact an anvi'o developer."
                        % str(self.translation_table))
            else:
                if self.prodigal_single_mode:
                    # the user explicitly requested to not use the `-p meta` flag to run
                    # prodigal (the default procedure is single, so prodigal will fall back
                    # to 'single' mode in this case)
                    pass
                else:
                    # Use 'meta' mode if no translation tables are given.
                    command.extend(['-p', 'meta'])

            self.commands.append(command)

        return State(commands=self.commands)

    def _collate_output_files(self):
        """Collates output files specified in `self.output_file_split_paths` to those specified in `self.collated_output_file_paths`.

        Returns `gene_calls_dict` and `amino_acid_sequences_dict` in `State`.

        See Superclass `_collate_output_files` for more info"""

        # todo check this in the constructor
        # Check that the proper keyword arguments are given.
        for key in ['peptide_path', 'gff_path']:
            if key not in self.output_file_split_paths:
                raise ValueError(f"kwargs should have the key '{key}', but it is missing.")

        peptide_paths = self.output_file_split_paths['peptide_path']
        gff_paths = self.output_file_split_paths['gff_path']

        gene_calls_dict, amino_acid_sequences_dict = self._process_peptide_files(peptide_paths)

        self._process_gff_files(gff_paths)

        return State({'gene_calls_dict': gene_calls_dict,
                      'amino_acid_sequences_dict': amino_acid_sequences_dict})

    # Helper methods
    #
    #
    def _process_peptide_files(self, peptide_paths):
        """Checks that `peptide_paths` files actually exist, then combines them."""
        # For renaming fasta headers
        hit_id = 0

        # Set up data storage.
        # Todo: These are keyed with ints, so idk why they aren't lists.  Maybe used as a dict later in the code?
        gene_calls_dict = {}  # each entry must contain {'contig', 'start', stop, 'direction', 'partial'} items.
        amino_acid_sequences_dict = {}

        for peptide_path in peptide_paths:
            if not os.path.exists(peptide_path):
                if self.logger:
                    self.logger.progress.end()

                raise ConfigError("Something went wrong with prodigal, and it failed to generate the "
                                  "expected output ('%s') :/ Fortunately, this log file should tell you what "
                                  "might be the problem: '%s'. Please do not forget to include this "
                                  "file if you were to ask for help." % (peptide_path, self.log_file_path))

            # Some splits may not actually have gene calls.  Skip them.
            if filesnpaths.is_file_empty(peptide_path):
                continue

            # If we get here, the fasta file will not be empty.
            fasta = fastalib.SequenceSource(peptide_path)

            while next(fasta):
                gene_calls_dict[hit_id] = self.parser(fasta.id)
                amino_acid_sequences_dict[hit_id] = fasta.seq.replace('*', '')
                hit_id += 1

            fasta.close()

            # todo i think this is removed elsewhere
            # Remove the split peptide file.
            os.remove(peptide_path)

        # If no genes were predicted across all output files, warn the user.
        if len(amino_acid_sequences_dict) == 0:
            if self.logger:
                self.logger.run.info('Result',
                                     f'Prodigal ({self.installed_version}) has identified no genes :/',
                                     nl_after=1,
                                     mc="red")
        else:  # Write out the final gene file
            assert 'peptide_path' in self.collated_output_file_paths
            with open(self.collated_output_file_paths['peptide_path'], 'w') as f:
                for hit_id, sequence in amino_acid_sequences_dict.items():
                    f.write(f">{hit_id}\n{sequence}\n")

        return gene_calls_dict, amino_acid_sequences_dict

    # Check the gff files and append them into the single outfile.
    def _process_gff_files(self, gff_paths):
        """Check the gff files and combine them.

        Note:  Currently genes_in_contigs (the gff file) is not used downstream.  This is important as each of the
        threads will start renumbering sequences from 1 in each subfile.  If the gff file itself is to be used, you MUST
        deal with this issue first!  See: https://github.com/merenlab/anvio/pull/1437#discussion_r440514864.
        """
        assert self.collated_output_file_paths['gff_path']

        # This function also checks if the files exist.
        concatenate_files(self.collated_output_file_paths['gff_path'], gff_paths, remove_concatenated_files=True)
