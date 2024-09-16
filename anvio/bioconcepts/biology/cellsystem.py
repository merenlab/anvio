from abc import ABC, abstractmethod

import anvio.bioconcepts.biology as biology

from anvio.bioconcepts.biology import Infomer, DNA, RNA, GeneProduct

class CellSystem(biology.System):
    """Cellular systems involve the flow of genetic information from DNA."""
    def __init__(self):
        super().__init__()

        self.premRNAs: dict[str, PremRNA] = {}
        self.mRNAs: dict[str, mRNA] = {}
        self.polypeptides: dict[str, Polypeptide] = {}
        self.proteins: dict[str, Protein] = {}

class Replicon(DNA):
    def __init__(self):
        biology.Replicon.__init__()
        self.genes: tuple[str] = None

class Gene(Infomer, DNA, biology.Gene):
    def __init__(self):
        biology.Gene.__init__()
        self.premRNAs: tuple[str] = None

    @property
    def products(self) -> tuple[str]:
        return self.premRNAs

class PremRNA(Infomer, GeneProduct, RNA):
    def __init__(self):
        self._sequence = None
        self.genes: tuple[str] = None
        self.mRNAs: tuple[str] = None

    @property
    def sequence(self):
        return self._sequence

    @property
    def products(self) -> tuple[str]:
        return self.mRNAs

    @property
    def precursors(self) -> tuple[str]:
        return self.genes

class mRNA(Infomer, GeneProduct, RNA):
    def __init__(self):
        self._sequence: str = None
        self.premRNAs: tuple[str] = None
        self.positions: dict[str, tuple[tuple[int]]] = None
        self.polypeptides: tuple[str] = None

    @property
    def sequence(self):
        return self._sequence

    @property
    def products(self) -> tuple[str]:
        return self.polypeptides

    @property
    def precursors(self) -> tuple[str]:
        return self.premRNAs

class Polypeptide(Infomer, GeneProduct, biology.Polypeptide):
    def __init__(self):
        self._sequence: str = None
        self.mRNAs: tuple[str] = None
        self.positions: dict[str, tuple[tuple[int]]] = None
        self.proteins: tuple[str] = None

    @property
    def sequence(self):
        return self._sequence

    @property
    def products(self) -> tuple[str]:
        return self.proteins

    @property
    def precursors(self) -> tuple[str]:
        return self.mRNAs

class Protein(Infomer, GeneProduct, biology.Polypeptide):
    def __init__(self):
        self._sequence: str = None
        self.polypeptides: tuple[str] = None
        self.preproteins: tuple[str] = None
        self.positions: tuple[tuple[int]] = None
        self.processed_proteins: tuple[str] = None

    @property
    def sequence(self):
        return self._sequence

    @property
    def products(self) -> tuple[str]:
        return self.processed_proteins

    @property
    def precursors(self) -> tuple[str]:
        return self.polypeptides + self.preproteins
