import sys

import anvio.terminal as terminal

from anvio.drivers.fasttree import FastTree
from anvio.drivers.muscle import Muscle
from anvio.drivers.famsa import FAMSA
from anvio.drivers.concoct import CONCOCT
from anvio.drivers.metabat2 import MetaBAT2
from anvio.drivers.maxbin2 import MaxBin2
from anvio.drivers.dastool import DAS_Tool
from anvio.drivers.binsanity import BinSanity
from anvio.errors import ConfigError

run = terminal.Run()
progress = terminal.Progress()

driver_modules = {}

phylogeny_default = "fasttree"
driver_modules["phylogeny"] = {"default": FastTree, "fasttree": FastTree}

driver_modules["binning"] = {
    "concoct": CONCOCT,
    "metabat2": MetaBAT2,
    "maxbin2": MaxBin2,
    "dastool": DAS_Tool,
    "binsanity": BinSanity,
}


class Aligners:
    def __init__(self, run=run, progress=progress):
        self.run = run
        self.progress = progress

        self.default = "muscle"

        self.aligners = {"default": Muscle, "muscle": Muscle, "famsa": FAMSA}

    def list(self):
        available_options = [o for o in self.aligners.keys() if o != "default"]

        self.run.warning(
            "The default anvi'o driver for multiple sequence alignment is '%s'. Available drivers "
            "are:" % (self.default),
            header="Multiple sequence alignment",
            lc="yellow",
        )
        for option in available_options:
            self.run.info_single(
                option,
                nl_after=1 if available_options[-1] == option else 0,
                mc="yellow",
            )

    def check(self, aligner):
        if aligner not in self.aligners:
            raise ConfigError(
                "Sorry, anvi'o knows nothing of the aligner '%s'. Nice try though :/"
                % aligner
            )

    def select(self, aligner=None, quiet=False):
        if not aligner:
            aligner = self.default

        self.check(aligner)

        _aligner = self.aligners[aligner]

        if not quiet:
            self.run.warning(
                "The workflow you are using will likely use '%s' by %s (%s) to align your sequences. "
                "If you publish your findings, please do not forget to properly credit this tool."
                % (aligner, _aligner().citation, _aligner().web),
                lc="yellow",
                header="CITATION",
            )

        return self.aligners[aligner]


if "--list-aligners" in sys.argv:
    aligners = Aligners()
    aligners.list()
    sys.exit()
