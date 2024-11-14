# coding: utf-8
"""Interface to PyANI."""

import os

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Özcan Esen"
__email__ = "ozcanesen@gmail.com"


class PyANI:
    def __init__(
        self,
        args={},
        run=terminal.Run(),
        progress=terminal.Progress(),
        program_name="average_nucleotide_identity.py",
    ):
        self.run = run
        self.progress = progress
        self.program_name = program_name

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.num_threads = A("num_threads") or 1
        self.method = A("method") or "ANIb"
        self.log_file_path = os.path.abspath(
            A("log_file") or filesnpaths.get_temp_file_path()
        )
        self.quiet = A("quiet")

        self.check_programs()

        self.run.warning(
            "Anvi'o will use 'PyANI' by Pritchard et al. (DOI: 10.1039/C5AY02550H) to compute ANI. If you publish your findings, \
                            please do not forget to properly credit their work.",
            lc="green",
            header="CITATION",
        )

        self.run.info("[PyANI] Num threads to use", self.num_threads)
        self.run.info("[PyANI] Alignment method", self.method)
        self.run.info("[PyANI] Log file path", self.log_file_path, nl_after=1)

    def check_programs(self):
        utils.is_program_exists(self.program_name)

        if self.method == "ANIb":
            utils.is_program_exists("blastn")
        elif self.method == "ANIblastall":
            utils.is_program_exists("blastall")
        elif self.method == "ANIm":
            utils.is_program_exists("nucmer")

    def run_command(self, input_path):
        # backup the old working directory before changing the directory
        old_wd = os.getcwd()
        os.chdir(input_path)

        full_command = [
            self.program_name,
            "--outdir",
            "output",
            "--indir",
            input_path,
            "--method",
            self.method,
            "--workers",
            self.num_threads,
        ]

        self.progress.new("PyANI")
        self.progress.update("Running ...")
        exit_code = utils.run_command(full_command, self.log_file_path)
        self.progress.end()

        if int(exit_code):
            raise ConfigError(
                "PyANI returned with non-zero exit code, there may be some errors. "
                "please check the log file for details."
            )

        output_matrix_names = [
            "alignment_coverage",
            "alignment_lengths",
            "hadamard",
            "percentage_identity",
            "similarity_errors",
            "correlations",
        ]

        full_matrix_path = lambda name: os.path.join(
            input_path, "output", self.method + "_" + name + ".tab"
        )

        matrices = {}
        for matrix_name in output_matrix_names:
            output_matrix_path = full_matrix_path(matrix_name)
            if os.path.exists(output_matrix_path):
                matrices[matrix_name] = utils.get_TAB_delimited_file_as_dictionary(
                    output_matrix_path, empty_header_columns_are_OK=True
                )

        if not len(matrices):
            raise ConfigError(
                "None of the output matrices pyANI was supposed to generate was found in the "
                "output directory :( You may find some clues in the log file?"
            )

        if not self.quiet:
            self.run.info_single(
                "Output matrices for the following items are stored in the output "
                "directory: %s <success kid meme.png>."
                % (", ".join(["'%s'" % m.replace("_", " ") for m in matrices])),
                nl_before=1,
                mc="green",
            )

        # restore old working directory
        os.chdir(old_wd)

        return matrices
