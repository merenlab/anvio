# coding: utf-8
"""Interface to MCL."""

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
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


class MCL:
    def __init__(self, mcl_input_file_path, run=run, progress=progress, num_threads=1):
        self.run = run
        self.progress = progress

        self.inflation = 2.0
        self.mcl_input_file_path = mcl_input_file_path
        self.num_threads = num_threads

        utils.is_program_exists("mcl")

        # if the programmer wishes to store the clusters output file in a particular
        # location, they will have to explicitly set the path after gettinga an
        # instance of this class
        self.clusters_file_path = filesnpaths.get_temp_file_path()

        if not self.run.log_file_path:
            self.run.log_file_path = filesnpaths.get_temp_file_path()

    def check_output(self, expected_output, process="diamond"):
        if not os.path.exists(expected_output):
            self.progress.end()
            raise ConfigError(
                "Pfft. Something probably went wrong with MCL's '%s' since one of the expected output files are missing. "
                "Please check the log file here: '%s."
                % (process, self.run.log_file_path)
            )

    def get_clusters_dict(self, name_prefix="GC"):
        self.cluster()

        clusters_dict = {}

        line_no = 1
        for line in open(self.clusters_file_path).readlines():
            cluster_name = f"{name_prefix}_{line_no:08d}"
            clusters_dict[cluster_name] = line.strip().split("\t")

            line_no += 1

        self.run.info("Number of MCL clusters", "%s" % pp(len(clusters_dict)))

        return clusters_dict

    def cluster(self):
        self.run.warning(None, header="MCL", lc="green")
        self.run.info("MCL inflation", self.inflation)

        self.progress.new("MCL")
        self.progress.update("clustering (using %d thread(s)) ..." % self.num_threads)

        cmd_line = [
            "mcl",
            self.mcl_input_file_path,
            "--abc",
            "-I",
            self.inflation,
            "-o",
            self.clusters_file_path,
            "-te",
            self.num_threads,
        ]

        self.run.info("mcl cmd", " ".join([str(x) for x in cmd_line]), quiet=True)

        utils.run_command(cmd_line, self.run.log_file_path)

        self.progress.end()

        self.check_output(self.clusters_file_path, "makedb")

        self.run.info("MCL output", self.clusters_file_path)
