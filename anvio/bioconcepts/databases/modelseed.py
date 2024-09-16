import anvio.bioconcepts.biology as biology
import anvio.bioconcepts.databases as databases

class Database:
    def __init__(self):
        self.compounds: dict[str, Compound] = {}
        self.reactions: dict[str, Reaction] = {}

class Compound(databases.Molecule):
    def __init__(self):
        self._molecules: tuple[biology.Molecule] = None

        self.id: str = None
        self.abbreviation: str = None
        self.name: str = None
        self.mass: float = None
        self.source: str = None
        self.is_core: int = None
        self.is_obsolete: int = None
        self.linked_compound: tuple[str] = None
        self.is_cofactor: int = None
        self.deltag: float = None
        self.deltagerr: float = None
        self.pka: tuple[tuple[int, int, float]] = None
        self.pkb: tuple[tuple[int, int, float]] = None
        self.abstract_compound: int = None
        self.comprised_of: int = None
        self.aliases: dict[str, tuple[str]] = None
        self.notes: str = None

    @property
    def molecules(self) -> tuple[biology.Molecule]:
        return self._molecules

class Reaction(databases.Reaction):
    def __init__(self):
        self._reactions: tuple[biology.Reaction] = None

        self.id: str = None
        self.abbreviation: str = None
        self.name: str = None
        self.code: str = None
        self.stoichiometry: str = None
        self.is_transport: int = None
        self.equation: str = None
        self.definition: str = None
        self.reversibility: str = None
        self.direction: str = None
        self.pathways: dict[str, tuple[str]] = None
        self.aliases: dict[str, tuple[str]] = None
        self.ec_numbers: tuple[str] = None
        self.deltag: float = None
        self.deltagerr: float = None
        self.compound_ids: tuple[str] = None
        self.status: str = None
        self.is_obsolete: int = None
        self.linked_reaction: tuple[str] = None
        self.notes: str = None
        self.source: str = None

    @property
    def reactions(self) -> tuple[biology.Reaction]:
        return self._reactions
