"""Interface to mmseqs2."""

import anvio

import anvio.terminal as terminal

__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Samuel Miller"
__email__ = "samuelmiller10@gmail.com"

run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print

class MMseqs2:
    def __init__(self, target_files_dict, num_threads_to_use=1, progress=progress, run=run):
        """A class to streamline MMseqs2 clustering."""
        self.num_threads_to_use = num_threads_to_use
        self.progress = progress
        self.run = run
        
        self.target_files_dict = {}
        
        for source in target_files_dict:
            pass