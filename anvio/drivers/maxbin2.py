# coding: utf-8
"""Interface to MaxBin2."""
import os
import glob

import anvio
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.errors import ConfigError

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Özcan Esen"
__email__ = "ozcanesen@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


class MaxBin2:
    arguments = {
        "min_contig_length": (
            ["--min-contig-length"],
            {
                "metavar": "INT",
                "required": False,
                "help": "Minimum contig length. Default: 1000.",
            },
        ),
        "max_iteration": (
            ["--max-iteration"],
            {
                "metavar": "INT",
                "required": False,
                "help": "Maximum Expectation-Maximization algorithm iteration number. Default 50.",
            },
        ),
        "prob_threshold": (
            ["--prob-threshold"],
            {
                "metavar": "FLOAT",
                "required": False,
                "help": "Probability threshold for EM final classification. Default 0.9.",
            },
        ),
        "markerset": (
            ["--merkerset"],
            {
                "metavar": "INT",
                "required": False,
                "help": "Minimum contig length. Default: 1000",
            },
        ),
    }

    citation = "Yu-Wei Wu, Blake A. Simmons, Steven W. Singer, MaxBin 2.0: \
                an automated binning algorithm to recover genomes from multiple \
                metagenomic datasets, Bioinformatics, Volume 32, Issue 4, 15 February 2016,\
                Pages 605–607, https://doi.org/10.1093/bioinformatics/btv638"

    cluster_type = "contig"

    def __init__(self, run=run, progress=progress):
        self.run = run
        self.progress = progress
        self.program_name = "run_MaxBin.pl"

        utils.is_program_exists(self.program_name)

    def cluster(self, input_files, args, work_dir, threads=1, log_file_path=None):
        J = lambda p: os.path.join(work_dir, p)

        output_file_prefix = J("MAXBIN_")

        if not log_file_path:
            log_file_path = J("logs.txt")

        cmd_line = [
            self.program_name,
            "-contig",
            input_files.contigs_fasta,
            "-abund",
            input_files.contig_coverages,
            "-out",
            output_file_prefix,
            "-thread",
            str(threads),
            *utils.serialize_args(args, single_dash=True, use_underscore=True),
        ]

        self.progress.new(self.program_name)
        self.progress.update("Running using %d threads..." % threads)
        utils.run_command(cmd_line, log_file_path)
        self.progress.end()

        output_file_paths = glob.glob(J(output_file_prefix + "*.fasta"))
        if not len(output_file_paths):
            raise ConfigError(
                "Some critical output files are missing. Please take a look at the "
                "log file: %s" % (log_file_path)
            )

        clusters = {}
        bin_count = 0

        for bin_file in output_file_paths:
            bin_count += 1
            with open(bin_file, "r") as f:
                bin_name = os.path.basename(bin_file).replace(".fasta", "")
                bin_name = bin_name.replace(".", "_")

                clusters[bin_name] = []

                for line in f.readlines():
                    if line.startswith(">"):
                        clusters[bin_name].append(line[1:].strip())

        return clusters
