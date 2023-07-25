# -*- coding: utf-8
# pylint: disable=line-too-long
"""Generate a metabolic reaction network from gene annotations."""

import os

from anvio.errors import ConfigError
from anvio import __version__ as VERSION


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2023, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = VERSION
__maintainer__ = "Samuel Miller"
__email__ = "samuelmiller10@gmail.com"
__status__ = "Development"


class Gene:
    """Representation of a gene in the metabolic network."""
    def __init__(self) -> None:
        self.gcid: int = None
        self.kos: List[KO] = []

class KO:
    """Representation of a KEGG Ortholog in the network."""
    def __init__(self) -> None:
        self.id: str = None
        self.name: str = None
        self.reactions: List[ModelSEEDReaction] = []
        # Record how KOs were associated with ModelSEED reactions, e.g., KEGG REACTION ID alias or EC number.
        self.modelseed_associations: List[str] = []
class ModelSEEDReaction:
    """Representation of a reaction in the network, with properties given by the ModelSEED Biochemistry database."""
    def __init__(self) -> None:
        self.modelseed_id: str = None
        self.modelseed_name: str = None
        self.kegg_id_aliases: Tuple[str] = []
        self.compounds: Tuple[ModelSEEDCompound] = []
        self.compartments: List[str] = []
        self.reversibility: bool = None

class ModelSEEDCompound:
    """Representation of a chemical in the network, with properties given by the ModelSEED Biochemistry database."""
    def __init__(self) -> None:
        self.modelseed_id: str = None
        self.modelseed_name: str = None
        self.kegg_id_aliases: List[str] = []
        self.charge: int = None
        self.formula: str = None
class Constructor:
    """
    Construct a metabolic reaction network within an anvi'o database.

    This currently uses KO annotations and the ModelSEED Biochemistry database.
    """
    def __init__(
        self,
        kegg_dir: str,
        modelseed_dir: str,
        progress: terminal.Progress = terminal.Progress()
    ) -> None:
        self.kegg_dir = kegg_dir
        self.modelseed_dir = modelseed_dir
        self.progress = progress
