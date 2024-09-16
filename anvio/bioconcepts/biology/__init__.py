from __future__ import annotations
from abc import ABC, abstractmethod

class System:
    def __init__(self):
        self.molecules: dict[str, Molecule] = {}
        self.reactions: dict[str, Reaction] = {}
        self.enzymes: dict[str, Enzyme] = {}

        self.replicons: dict[str, Replicon] = {}
        self.genes: dict[str, Gene] = {}

class Molecule:
    def __init__(self):
        self.charge: int = None
        self.formula: str = None
        self.smiles: str = None
        self.inchikey: str = None

class Reaction:
    def __init__(self):
        self.molecules: tuple[str] = None
        self.coefficients: tuple[int] = None

class Enzyme:
    def __init__(self):
        self.reactions: tuple[str] = None

class DNA(ABC):
    pass

class RNA(ABC):
    pass

class Replicon:
    def __init__(self):
        pass

class Infomer(ABC):
    @property
    @abstractmethod
    def sequence(self) -> str:
        pass

    @property
    @abstractmethod
    def products(self) -> tuple[str]:
        pass

class Gene(Infomer):
    def __init__(self):
        self._sequence: str = None
        self.replicon: str = None
        self.positions: tuple[tuple[int]] = None

    def sequence(self):
        return self._sequence

    def products(self):
        pass

class GeneProduct(ABC):
    @property
    @abstractmethod
    def precursors(self, system: System):
        pass

class Polypeptide(ABC):
    pass
