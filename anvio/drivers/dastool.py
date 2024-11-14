# -*- coding: utf-8
# pylint: disable=line-too-long
"""
"""
import os

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.ccollections as ccollections

from anvio.errors import ConfigError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Özcan Esen"
__email__ = "ozcanesen@gmail.com"


run = terminal.Run()
progress = terminal.Progress()


class DAS_Tool:
    arguments = {
        "source-collections": (
            ["-S", "--source-collections"],
            {
                "metavar": "COLLECTION_LIST",
                "required": True,
                "help": "Comma-separated list of collections, case sensitive.",
            },
        ),
        "search-engine": (
            ["--search-engine"],
            {
                "metavar": "PROGRAM",
                "required": False,
                "default": "usearch",
                "help": "Engine used for single copy gene identification [blast/diamond/usearch].\
                              (default: usearch)",
            },
        ),
        "score-threshold": (
            ["--score-threshold"],
            {
                "metavar": "FLOAT",
                "required": False,
                "help": "Score threshold until selection algorithm will keep selecting bins [0..1].\
                              (default: 0.5)",
            },
        ),
        "duplicate-penalty": (
            ["--duplicate-penalty"],
            {
                "metavar": "FLOAT",
                "required": False,
                "help": "Penalty for duplicate single copy genes per bin (weight b).\
                              Only change if you know what you're doing. [0..3]\
                              (default: 0.6)",
            },
        ),
        "megabin-penalty": (
            ["--megabin-penalty"],
            {
                "metavar": "FLOAT",
                "required": False,
                "help": "Penalty for megabins (weight c). Only change if you know what you're doing. [0..3]\
                              (default: 0.5)",
            },
        ),
        "db-directory": (
            ["--db-directory"],
            {
                "metavar": "PATH",
                "required": False,
                "default": None,
                "help": "Directory of single copy gene database. (default: install_dir/db)",
            },
        ),
    }

    citation = "Christian M. K. Sieber, Alexander J. Probst, Allison Sharrar, Brian C. Thomas, \
                Matthias Hess, Susannah G. Tringe & Jillian F. Banfield (2018). \
                Recovery of genomes from metagenomes via a dereplication, aggregation \
                and scoring strategy. Nature Microbiology. https://doi.org/10.1038/s41564-018-0171-1."

    cluster_type = "split"

    def __init__(self, run=run, progress=progress):
        self.run = run
        self.progress = progress

        self.program_name = "DAS_Tool"

        utils.is_program_exists(self.program_name)

    def cluster(self, input_files, args, work_dir, threads=1, log_file_path=None):
        J = lambda p: os.path.join(work_dir, p)

        cwd_backup = os.getcwd()
        os.chdir(work_dir)

        if not log_file_path:
            log_file_path = J("logs.txt")

        c = ccollections.Collections(r=run, p=progress)
        c.populate_collections_dict(input_files.profile_db)

        source_collections = set(map(str.strip, args.source_collections.split(",")))

        missing_collections = source_collections - set(c.collections_dict.keys())

        if len(missing_collections):
            raise ConfigError(
                "Some of the collections you wanted are missing in the database. "
                "Here is the list of missing collections: %s"
                % (", ".join(missing_collections))
            )

        c_names = []
        c_files = []

        for collection_name in source_collections:
            prefix = J(collection_name)

            c_names.append(collection_name)
            c_files.append(prefix + ".txt")

            c.export_collection(
                collection_name, output_file_prefix=prefix, include_unbinned=False
            )

        cmd_line = [
            self.program_name,
            "-c",
            input_files.splits_fasta,
            "-i",
            ",".join(c_files),
            "-l",
            ",".join(c_names),
            "-o",
            J("OUTPUT"),
            "--threads",
            str(threads),
            *utils.serialize_args(
                args, use_underscore=True, skip_keys=["source_collections"]
            ),
        ]

        self.progress.new(self.program_name)
        self.progress.update("Running using %d threads..." % threads)
        utils.run_command(cmd_line, log_file_path)
        self.progress.end()

        output_file_name = "OUTPUT_DASTool_scaffolds2bin.txt"
        output_file_path = J(output_file_name)
        if not os.path.exists(output_file_path):
            # if this output file is missing, we may find the other one
            # perhaps
            output_file_name = "OUTPUT_DASTool_contig2bin.tsv"
            output_file_path = J(output_file_name)
            if not os.path.exists(output_file_path):
                raise ConfigError(
                    "One of the critical output files is missing ('%s'). Please take a look at the "
                    "log file: %s" % (output_file_name, log_file_path)
                )

        clusters = {}
        with open(output_file_path, "r") as f:
            lines = f.readlines()

            for entry in lines:
                contig, bin_name = map(str.strip, entry.split())

                pretty_bin_name = "Bin_" + bin_name.replace(".", "_")

                if pretty_bin_name not in clusters:
                    clusters[pretty_bin_name] = []

                clusters[pretty_bin_name].append(contig)

        # restore cwd
        os.chdir(cwd_backup)

        return clusters
