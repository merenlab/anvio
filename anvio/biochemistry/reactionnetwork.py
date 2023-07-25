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

class KEGGDatabase:
    """The KEGG KO and REACTION databases set up by anvi'o."""
    default_dir = os.path.join(os.path.dirname(ANVIO_PATH), 'data/MISC/PROTEIN_DATA/kegg')

    def __init__(self) -> None:
        # The KO and reaction tables are derived from the downloaded definition files. They
        # facilitate the lookup of KO IDs, names, EC numbers, and KEGG reactions.
        self.ko_table: pd.DataFrame = None
        self.reaction_table: pd.DataFrame = None

    def load(self, db_dir: str = None) -> None:
        """Load KO and reaction tables from the data directory."""
        if db_dir:
            if not os.path.isdir(db_dir):
                raise ConfigError(f"The provided KEGG database directory, '{db_dir}', was not recognized as a directory.")
        else:
            db_dir = KEGGDatabase.default_dir
        ko_data_path = os.path.join(db_dir, 'ko_data.tsv')
        if not os.path.isfile(ko_data_path):
            raise ConfigError(f"The KO data table, 'ko_data.tsv', was not found in the database directory, '{db_dir}'.")
        reaction_data_path = os.path.join(db_dir, 'reaction_data.tsv')
        if not os.path.isfile(reaction_data_path):
            raise ConfigError(f"The KEGG REACTION data table, 'reaction_data.tsv', was not found in the database directory, '{db_dir}'.")

        self.ko_table = pd.read_csv(ko_data_path, sep='\t', header=0, index_col=0, low_memory=False)
        self.reaction_table = pd.read_csv(reaction_data_path, sep='\t', header=0, index_col=0, low_memory=False)

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

        self.kegg_db = KEGGDatabase()
        self.kegg_db.load(self.kegg_dir)
