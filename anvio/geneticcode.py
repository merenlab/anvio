#!/usr/bin/env python
# -*- coding: utf-8
"""
Genetic code objects are especially useful for getting information on synonymous encodings for codon
usage analyses.
"""

import anvio.terminal as terminal
import anvio.constants as constants

from anvio.errors import ConfigError

class GeneticCode:
    """
    A genetic code relating codons to amino acids and termination encoding.

    Attributes
    ==========
    codon_amino_acid : dict[str, str]
        Keys are codons. Values are three-letter codes for encoded amino acids, plus "STP" for stop
        codons.

    amino_acid_codons : dict[str, list[str]]
        Keys are three-letter codes for encoded amino acids, plus "STP" for stop codons. Values are
        codons.

    nonstop_codon_amino_acid : dict[str, str]
        Keys are codons. Values are three-letter codes for encoded amino acids, excluding stop
        codons.

    nonstop_amino_acid_codons : dict[str, list[str]]
        Keys are three-letter codes for encoded amino acids, excluding stop codons. Values are
        codons.

    synonymous_codon_amino_acid : dict[str, str]
        Keys are synonymous codons. Values are three-letter codes for encoded degenerate amino
        acids, excluding stop codons.

    synonymous_amino_acid_codons : dict[str, list[str]]
        Keys are three-letter codes for encoded degenerate amino acids, excluding stop codons.
        Values are synonymous codons.

    run : anvio.terminal.Run, anvio.terminal.Run()
        Prints run information to the terminal.
    """
    standard_code: dict[str, str] = {codon: aa for codon, aa in constants.codon_to_AA.items()}

    def __init__(self, code: dict[str, str] = None, run: terminal.Run = terminal.Run()):
        """
        Parameters
        ==========
        code : dict[str, str]
            Keys are codons. Values are three-letter codes for encoded amino acids, which can
            include "STP" for stop codons. If None, the standard genetic code is used.

        run : anvio.terminal.Run, anvio.terminal.Run()
            Prints run information to the terminal.
        """
        self.run = run
        self.set(code)

    def set(self, code: dict[str, str] = None):
        """
        Set the genetic code.

        Parameters
        ==========
        code : dict[str, str], None
            Keys are codons. Values are three-letter codes for encoded amino acids, which can
            include "STP" for stop codons. If None, the standard genetic code is used.
        """
        if code is None:
            self.codon_amino_acid = self.standard_code
        else:
            self.codon_amino_acid = code
            self.check()
            if self.codon_amino_acid != self.standard_code:
                self.run.info_single("A nonstandard genetic code was set.")

        self.amino_acid_codons: dict[str, list[str]] = {}
        for codon, aa in self.codon_amino_acid.items():
            try:
                self.amino_acid_codons[aa].append(codon)
            except KeyError:
                self.amino_acid_codons[aa] = [codon]

        self.nonstop_codon_amino_acid: dict[str, str] = {}
        self.nonstop_amino_acid_codons: dict[str, list[str]] = {}
        for codon, aa in self.codon_amino_acid.items():
            if aa == 'STP':
                continue
            self.nonstop_codon_amino_acid[codon] = aa
            try:
                self.nonstop_amino_acid_codons[aa].append(codon)
            except KeyError:
                self.nonstop_amino_acid_codons[aa] = [codon]

        self.synonymous_codon_amino_acid: dict[str, str] = {}
        self.synonymous_amino_acid_codons: dict[str, list[str]] = {}
        for aa, codons in self.nonstop_amino_acid_codons.items():
            if len(codons) == 1:
                continue
            for codon in codons:
                self.synonymous_codon_amino_acid[codon] = aa
            self.synonymous_amino_acid_codons[aa] = codons.copy()

    def check(self):
        """
        Check that recognized codons and three-letter amino acid codes are used in the dictionary
        mapping codon to amino acid, raising an exception if untrue.
        """
        unrecognized_codons = []
        unrecognized_amino_acids = []
        for codon, amino_acid in self.codon_amino_acid.items():
            if codon not in constants.codons:
                unrecognized_codons.append(codon)
            if amino_acid not in constants.amino_acids:
                unrecognized_amino_acids.append(amino_acid)

        if unrecognized_codons:
            unrecognized_codon_message = (
                "The following codons in the genetic code are not recognized: "
                f"{', '.join(unrecognized_codons)}."
            )
        else:
            unrecognized_codon_message = ""

        if unrecognized_amino_acids:
            unrecognized_amino_acid_message = (
                "The following amino acids in the genetic code are not recognized: "
                f"{', '.join(unrecognized_amino_acids)}. These should be three-letter codes and "
                "\"STP\" for stop codons."
            )
            if unrecognized_codon_message:
                unrecognized_amino_acid_message = " " + unrecognized_amino_acid_message
        else:
            unrecognized_amino_acid_message = ""

        if unrecognized_codon_message or unrecognized_amino_acid_message:
            raise ConfigError(f"{unrecognized_codon_message}{unrecognized_amino_acid_message}")

    def modify(self, encodings: dict[str, str] = None, encodings_txt: str = None):
        """
        Modify the genetic code.

        Parameters
        ==========
        encodings : dict[str, str]
            Changes to the genetic code. Keys are codons. Values are three-letter codes for encoded
            amino acids, which can include "STP" for stop codons.

        encodings_txt : str
            Path to a tab-delimited file with information like that passed to the 'encodings'
            option. The file should have no header, codons in the first column, and amino acids in
            the second column.
        """
        if (encodings is None) + (encodings_txt is None) != 1:
            raise ConfigError("Provide either 'encodings' or 'encodings_txt'.")

        new_code = self.codon_amino_acid.copy()

        if encodings is not None:
            new_code.update(encodings)
            self.set(new_code)
            return

        filesnpaths.is_file_tab_delimited(encodings_txt)
        encodings_df = pd.read_csv(encodings_txt, sep='\t', header=None)
        new_code.update(dict(zip(encodings_df.iloc[:, 0], encodings_df.iloc[:, 1])))
        self.set(new_code)
