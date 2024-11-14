#!/usr/bin/env python
# -*- coding: utf-8

import anvio
import anvio.terminal as terminal
import anvio.constants as constants

from anvio.parsers.base import Parser


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


class KrakenUniq(Parser):
    def __init__(
        self,
        input_files,
        taxonomy_table_structure,
        run=terminal.Run(),
        progress=terminal.Progress(),
    ):
        self.run = run
        self.progress = progress

        if type(input_files) != type(list()):
            input_files = [input_files]

        files_expected = {"kraken_output": input_files[0]}

        files_structure = {
            "kraken_output": {
                "col_names": ["taxonomy", "count"],
                "col_mapping": [str, int],
                "separator": "\t",
                "indexing_field": -1,
            },
        }

        Parser.__init__(
            self, "KrakenUniq", input_files, files_expected, files_structure
        )

    def process(self):
        """The file to be parsed looks like this:

        #Sample ID	krakenhll_files/CP_R01_CDI_P_01_PRE.krakenhll
        d__Viruses	63623
        d__Viroids	0
        d__Bacteria	6535697
        d__Archaea	665
        d__Eukaryota	0
        d__Eukaryota|k__Viridiplantae	0
        d__Eukaryota|k__Fungi	0
        d__Eukaryota|k__Metazoa	0
        d__Bacteria|p__Proteobacteria	2601002
        d__Bacteria|p__Candidatus_Acetothermia	0
        d__Bacteria|p__Candidatus_Aminicenantes	0
        d__Bacteria|p__Candidatus_Atribacteria	0
        d__Bacteria|p__Candidatus_Saccharibacteria	4
        d__Bacteria|p__Candidatus_Poribacteria	0
        d__Bacteria|p__candidate_division_WWE3	0
        d__Bacteria|p__candidate_division_ZB3	0
        d__Bacteria|p__candidate_division_NC10	0
        d__Bacteria|p__Candidatus_Calescamantes	0
        d__Bacteria|p__Candidatus_Aerophobetes	0

        """

        taxonomy_dict = {}
        taxonomic_descendants = {}
        for level in constants.levels_of_taxonomy:
            taxonomy_dict[level] = {}

        kraken_output = self.dicts["kraken_output"]

        self.run.info("Total num hits found", len(kraken_output))

        # Kraken file contains crazy number of hits with 0 count.
        entries_with_zero_hits = [
            e for e in kraken_output if not kraken_output[e]["count"]
        ]
        for entry in entries_with_zero_hits:
            kraken_output.pop(entry)

        self.run.info("Total num non-zero-count hits", len(kraken_output))

        self.progress.new("Bleep kraken stuff bloop")
        self.progress.update("Processing the input data ...")

        for entry in kraken_output.values():
            tax_string_list = [
                f for f in entry["taxonomy"].split("|") if not f.startswith("k__")
            ]
            count = entry["count"]

            num_levels_tax_entry = len(tax_string_list)

            # don't include the synthetic construct
            if (
                num_levels_tax_entry == 1
                and tax_string_list[0] == "s__synthetic_construct"
            ):
                continue

            tax_level_described = constants.levels_of_taxonomy[num_levels_tax_entry - 1]
            taxon_names = [t[3:] for t in tax_string_list]
            taxon = taxon_names[-1]
            ancestor = taxon_names[-2] if num_levels_tax_entry > 1 else None

            taxonomy_dict[tax_level_described][taxon] = count

            if ancestor:
                if ancestor not in taxonomic_descendants:
                    taxonomic_descendants[ancestor] = set([])

                taxonomic_descendants[ancestor].add(taxon)

        self.progress.end()

        # some taxa will only resolve to a higher taxonomic level. for instance, you may have X number of
        # phylum Bacteroidetes hits, but all classes that belong to this phylum may sum up to `X - d`. Here
        # we add that `d` as 'Unknown_phlum_Bacteroidetes' in the next taxonomic level.
        for level in constants.levels_of_taxonomy[:-1]:
            next_level = constants.levels_of_taxonomy[
                constants.levels_of_taxonomy.index(level) + 1
            ]

            for taxon in taxonomy_dict[level]:
                if taxon.startswith("Unknown"):
                    name_for_the_unknown_in_descendant = taxon
                else:
                    name_for_the_unknown_in_descendant = "Unknown_%s_%s" % (
                        level,
                        taxon,
                    )

                if taxon not in taxonomic_descendants:
                    taxonomic_descendants[taxon] = set(
                        [name_for_the_unknown_in_descendant]
                    )
                    taxonomy_dict[next_level][name_for_the_unknown_in_descendant] = (
                        taxonomy_dict[level][taxon]
                    )
                    continue

                total_counts_descendants_describe = 0
                for descendant in taxonomic_descendants[taxon]:
                    if descendant in taxonomy_dict[next_level]:
                        total_counts_descendants_describe += taxonomy_dict[next_level][
                            descendant
                        ]

                total_counts_described_by_ancestor = taxonomy_dict[level][taxon]
                unaccounted_for = (
                    total_counts_described_by_ancestor
                    - total_counts_descendants_describe
                )

                if unaccounted_for:
                    taxonomy_dict[next_level][
                        name_for_the_unknown_in_descendant
                    ] = unaccounted_for
                    taxonomic_descendants[taxon].add(name_for_the_unknown_in_descendant)

        return taxonomy_dict
