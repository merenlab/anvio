# coding: utf-8
"""
Interface to MODELLER (https://salilab.org/modeller/).
"""

import os
import tempfile

import anvio
import anvio.utils as utils
import anvio.fastalib as u
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError


__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2016, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print

class MODELLER:

    def __init__(self, target_fasta_path, directory=None, run=run, progress=progress):
        self.run = run
        self.progress = progress

        # check that MODELLER exists
        utils.is_program_exists('mod9.19') # FIXME
        self.mod = "mod9.19"

        # All MODELLER scripts are housed in self.script_folder
        self.scripts_folder = "/Users/evan/Software/anvio/anvio/data/misc/MODELLER/scripts" # FIXME
        if utils.filesnpaths.is_dir_empty(self.scripts_folder):
            raise ConfigError("Anvi'o houses all its MODELLER scripts in {}, but your directory \
                               contains no scripts. You should change that.")

        # does target_fasta_path point to a fasta file?
        self.target_fasta_path = target_fasta_path
        utils.filesnpaths.is_file_fasta_formatted(self.target_fasta_path)
        run.info("target fasta", self.target_fasta_path)

        # make sure target_fasta is valid
        target_fasta = u.SequenceSource(target_fasta_path, lazy_init=False)
        if target_fasta.total_seq != 1:
            raise ConfigError("MODELLER::The input FASTA file must have exactly one sequence.\
                               You provided one with {}.".format(target_fasta.total_seq))

        # get gene_id while target_fasta is opened
        while next(target_fasta):
            self.gene_id = target_fasta.id
        target_fasta.close()

        # self.directory is defined because there are a lot of extraneous output and log
        # files output by MODELLER so we cd into directory, do our business, and then cd back
        # into the starting directory.
        self.directory = directory
        self.start_dir = os.getcwd()
        os.chdir(self.directory)


    def validate_and_return_script_path(self, script_name):
        """
        All MODELLER scripts are housed in anvio/data/misc/MODELLER/scripts/
        """
        script_path = os.path.join(self.scripts_folder, script_name)
        try:
            utils.filesnpaths.is_file_exists(script_path)
        except:
            raise ConfigError("MODELLER :: The script {} is not in {}".format(script_name, self.scripts_folder))

        return script_path


    def fasta_to_pir(self):
        """
        MODELLER uses their own .pir format for search and alignment instead of .fasta. This script
        does the conversion. An example is found at https://salilab.org/modeller/tutorial/basic.html
        """

        script_name = "fasta_to_pir.py"
        script_path = self.validate_and_return_script_path(script_name)
        #self.target_pir = os.path.join(self.directory, "{}.pir".format(self.gene_id))
        self.target_pir = "{}.pir".format(self.gene_id)

        command = [self.mod,
                   script_path,
                   self.target_fasta_path,
                   self.target_pir]
        command = " ".join(command)

        self.run_command(command, script_path)
        run.info("target pir", self.target_pir)


    def run_command(self, command, script_path):
        os.system(command)
        self.move_log_to_working_directory(script_path)


    def move_log_to_working_directory(self, script_path):
        """ MODELLER outputs a log file in the directory of the script. Move it to self.directory"""
        log_name = os.path.splitext(os.path.basename(script_path))[0] + ".log"
        source = os.path.join(os.path.dirname(script_path), log_name)
        destination = log_name
        os.rename(source, destination)

    def search(self):
        pass


    def close(self):
        """
        Change directories back to self.start_dir. Could also have options selectively bring back
        log files and PDB files to self.start_dir if self.directory != self.start_dir.
        """
        # make back to starting directory
        os.chdir(self.start_dir)




