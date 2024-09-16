from abc import ABC, abstractmethod

import anvio.bioconcepts.biology as biology

class Entry(ABC):
    pass

class Molecule(Entry):
    @property
    @abstractmethod
    def molecules(self) -> tuple[biology.Molecule]:
        pass

class Reaction(Entry):
    @property
    @abstractmethod
    def reactions(self) -> tuple[biology.Reaction]:
        pass

class Ortholog(Entry):
    @property
    @abstractmethod
    def genes(self) -> tuple[biology.Gene]:
        pass
