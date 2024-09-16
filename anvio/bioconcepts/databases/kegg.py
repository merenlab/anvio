import anvio.bioconcepts.biology as biology
import anvio.bioconcepts.databases as databases

class Database:
    def __init__(self):
        self.reactions: dict[str, Reaction] = {}
        self.orthologs: dict[str, Ortholog] = {}

class Reaction(databases.Reaction):
    def __init__(self):
        self._reactions: tuple[biology.Reaction] = None

        self.id: str = None
        self.name: str = None
        self.definition: str = None
        self.equation: str = None
        self.remark: str = None
        self.comment: str = None
        self.rclass: dict[str, tuple[str]] = None
        self.enzyme: tuple[str] = None
        self.pathway: dict[str, str] = None
        self.module: dict[str, str] = None
        self.brite: dict[str, str] = None
        self.dblinks: dict[str, tuple[str]] = None
        self.orthology: dict[str, str] = None
        self.reference: str = None
        self.authors: str = None
        self.title: str = None
        self.journal: str = None
        self.sequence: str = None

    @property
    def reactions(self) -> tuple[biology.Reaction]:
        return self._reactions

class Ortholog(databases.Ortholog):
    def __init__(self):
        self._genes: tuple[biology.Gene] = None

        self.id: str = None
        self.symbol: tuple[str] = None
        self.name: str = None
        self.pathway: dict[str, str] = None
        self.module: dict[str, str] = None
        self.reaction: tuple[str] = None
        self.network: dict[str, str] = None
        self.element: dict[str, str] = None
        self.disease: dict[str, str] = None
        self.brite: dict[str, str] = None
        self.dblinks: dict[str, tuple[str]] = None
        self.genes: dict[str, tuple[str]] = None
        self.reference: str = None
        self.authors: str = None
        self.title: str = None
        self.journal: str = None
        self.sequence: str = None

    @property
    def genes(self) -> tuple[biology.Gene]:
        return self._genes
