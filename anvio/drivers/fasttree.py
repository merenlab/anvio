# coding: utf-8
"""Interface to FastTree."""

from subprocess import Popen, PIPE

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths


__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Özcan Esen"
__email__ = "ozcanesen@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


class FastTree:
    def __init__(self, run=run):
        self.run = run
        self.progress = progress
        self.command = ['FastTree']

        utils.is_program_exists('FastTree')

    def run_command(self, input_file_path, output_file_path):
        input_file = open(input_file_path, 'rb')

        fasttree = Popen(self.command, stdout=PIPE, stdin=PIPE, stderr=PIPE)
        output = fasttree.communicate(input=input_file.read())
        input_file.close()

        output_stdout = output[0].decode().rstrip()
        output_stderr = output[1].decode().splitlines()

        run.info("Version", output_stderr[0])
        warning = ""
        for line in output_stderr[1:]:
            if len(warning) > 0 or line.startswith("WARNING! "):
                warning += line + "\n"
                if line == "":
                    run.warning(warning)
                    warning = ""
            elif line.startswith("      "):
                pass
            elif 'seconds' in line:
                pass
            else:
                line = line.split(":")
                if len(line) == 2:
                    run.info(line[0], line[1].strip())
                else:
                    run.info("Info", ":".join(line))

        if filesnpaths.is_proper_newick(output_stdout):
            output_file = open(output_file_path, 'w')
            output_file.write(output_stdout + '\n')
            output_file.close()

            run.info('FastTree output newick file', output_file_path, mc='green', nl_before=1, nl_after=1)
