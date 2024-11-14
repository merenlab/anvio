#!/usr/bin/env python
# -*- coding: utf-8

import os
import random

import anvio
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.parsers.base import Parser
from anvio.parsers.base import TaxonomyHelper


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


class Kaiju(Parser):
    def __init__(
        self,
        input_files,
        taxonomy_table_structure,
        run=terminal.Run(),
        progress=terminal.Progress(),
        skip_fix_input=False,
    ):
        self.run = run
        self.progress = progress
        self.just_do_it = False

        if type(input_files) != type(list()):
            input_files = [input_files]

        if not skip_fix_input:
            input_files[0] = self.fix_input_file(input_files[0])

        files_expected = {"kaiju_output": input_files[0]}

        files_structure = {
            "kaiju_output": {
                "col_names": [
                    "_",
                    "gene_callers_id",
                    "_",
                    "_",
                    "_",
                    "_",
                    "_",
                    "taxonomy",
                ],
                "col_mapping": [str, int, str, int, str, str, str, str],
                "separator": "\t",
                "indexing_field": -1,
            },
        }

        Parser.__init__(self, "Kaiju", input_files, files_expected, files_structure)

        if not skip_fix_input:
            os.remove(input_files[0])

    def fix_input_file(self, input_file_path):
        """This is sadly necessary because Kaiju output contains either three,
        or eight TAB-delimited columns ... not very parser friendly"""

        filesnpaths.is_file_exists(input_file_path)

        self.progress.new("Fixing the broken kaiju output")
        self.progress.update("...")

        corrected_temp_file_path = filesnpaths.get_temp_file_path()
        corrected_temp_file = open(corrected_temp_file_path, "w")
        input_file = open(input_file_path, "r")

        num_correct_lines = 0
        for line in input_file.readlines():
            if len(line.split("\t")) == 8:
                corrected_temp_file.write(line)
                num_correct_lines += 1

        corrected_temp_file.close()

        self.progress.end()

        if not num_correct_lines:
            os.remove(corrected_temp_file_path)
            raise ConfigError(
                "Something must have been wrong with you input file. Not a single line in it "
                "matched to what the kaiju parsers expects (a proper file should have eight "
                "TAB-delimited columns)."
            )

        return corrected_temp_file_path

    def process(self):
        """The file to be parsed looks like this:

                C	1	340016	248	340016,	ASF00353.1,	KPDGLIIDIEGLENVQLGKGGELQPLDLHDIYEQTGVFYYRSKNPEGG,	NA; NA; NA; NA; NA; uncultured virus;
                C	2	745014	427	745014,	WP_009471064.1,	GISENIRAISVIDRYLEHPRVYIAYSRGEPKYYMGSADLMTRNIDYRVEVLCPVHDPKAQKTLQDVLDQQWNDNVKARVIDASQ,	Proteobacteria; Cellvibrionales; Gammaproteobacteria; Halieaceae; NA; gamma proteobacterium HIMB55;
                C	3	745014	564	745014,	WP_009471063.1,	MSLIASMARGMDESNLKMNAAGNIQSENTSGVAERRMRMPPIPDDFFSLITHEQRIALNQKQQFGWFVKFVRRPLFQPIEVVLSNPEGSEFLLLETDGITRPFFNVRTDDLR,	Proteobacteria; Cellvibrionales; Gammaproteobacteria; Halieaceae; NA; gamma proteobacterium HIMB55;
                C	4	745014	991	745014,	WP_009471062.1,	SNPRWEQLLRYRYIELIALWEGRLTTRQLCETFGIGRQQANKDLTSYRRGLTRGDLVYDAVAKYYSPSEDFAPTLTQGLASEYLQMAAQQSDVQQILGDLPVASANVEVIAAPLREVPASLLRPIIRAMAESRRIDVDYVSLNNPDREGRIIVPHTLVWTGYRWHVRAWCEKNLDFRDFVLSRFRGDADLMD,FSNPRWEQLLRYRYIELIALWEGRLTTRQLCETFGIGRQQANKDLTSYRRGLTRGDLVYDAVAKYYSPSEDFAPTLTQGLASEYLQMAAQQSDVQQILGDLPVASANVEVIAAPLREVPASLLRPIIRAMAESRRIDVDYVSLNNPDREGRIIVPHTLVWTGYRWHVRAWCEKNLDFRDFVLSRFRGDADLMD,	Proteobacteria; Cellvibrionales; Gammaproteobacteria; Halieaceae; NA; gamma proteobacterium HIMB55;
                C	5	745014	737	745014,	WP_009471061.1,	GLGGLEALRLFTDVLAPACISTDHPKFLAFVPAAPTEAATLFDLIVGASSICGTSWLESAGATYAENQALQWIADLAGFGPEAGGTFVSGGTAGNLSALVAARHKWRRGNESRDALRGLVISSKGAHASIKQATYVMDVDLLEVGGD,Proteobacteria; Cellvibrionales; Gammaproteobacteria; Halieaceae; NA; gamma proteobacterium HIMB55;
                C	6	745014	462	745014,	WP_009471060.1,	SFRECLSNPLFVGYQLSGTFSFCGVFVYISTVAFFLRDVFDVSTEFFGVVFAMTAAGFIVGSLSSSRLVLKWGADRTLRRGAFICALSTTSAL,	Proteobacteria; Cellvibrionales; Gammaproteobacteria; Halieaceae; NA; gamma proteobacterium HIMB55;
                C	7	745014	583	745014,	WP_009471059.1,	MSTNRSYVSATLTADENKAAIEAHLHEILERSLTPMEPGQAKVYMEHTAVRMAEEAGAGVTTFQMVEVKHANTAYMIRLAVLTNGSAIGLDLMDMENGQFFIPEVCPVIPLETPTVN,	Proteobacteria; Cellvibrionales; Gammaproteobacteria; Halieaceae; NA; gamma proteobacterium HIMB55;
                U	8	0
                C	11	1898104	296	1898104,	PTL97260.1,	YFEPWVKGGNSIIRAIHYPPITTDPGDSVRAGQHEDINLITLLMGASAEGLEVLNKQG,	Bacteroidetes; NA; NA; NA; NA; Bacteroidetes bacterium;
                U	13	0

        Where, according to https://github.com/bioinformatics-centre/kaiju, each column corresponds to,

                1. either C or U, indicating whether the read is classified or unclassified.
                2. name of the read
                3. NCBI taxon identifier of the assigned taxon
                4. the length or score of the best match used for classification
                5. the taxon identifiers of all database sequences with the best match
                6. the accession numbers of all database sequences with the best match
                7. matching fragment sequence(s)
        """

        if not self.just_do_it:
            raise ConfigError(
                "Anvi'o assumes you used this exact parameter during your kaiju run: "
                "'-r superkingdom,phylum,class,order,family,genus,species'. If you "
                "haven't, you will run into trouble later. If you are positive that "
                "you did include that parameter to your run, re-run this program with "
                "`--just-do-it` flag."
            )

        # THIS IS IMPORTANT.
        levels_of_taxonomy = constants.levels_of_taxonomy

        taxonomy_dict = {}

        kaiju_output = self.dicts["kaiju_output"]

        self.run.info("Total num hits found", len(kaiju_output))

        self.progress.new("Bleep kaiju stuff bloop")
        self.progress.update("Processing the input data ...")

        for entry in kaiju_output.values():
            tax_string_list = [e.strip() for e in entry["taxonomy"].split(";")]
            gene_callers_id = entry["gene_callers_id"]

            last_known = ""
            for i in range(0, len(tax_string_list)):
                if tax_string_list[i] == "NA":
                    if last_known:
                        tax_string_list[i] = "Unknown_" + last_known
                    else:
                        tax_string_list[i] = "Unknown"
                else:
                    last_known = tax_string_list[i].replace(" ", "_")

            taxonomy_dict[gene_callers_id] = {}
            for i in range(0, len(levels_of_taxonomy)):
                level = levels_of_taxonomy[i]
                taxonomy_dict[gene_callers_id][level] = tax_string_list[i]

        self.progress.end()

        random_phylum_names = set(
            [
                taxonomy_dict[e]["t_phylum"]
                for e in random.sample(list(taxonomy_dict.keys()), 20)
            ]
        )

        self.run.warning(
            "Good news: anvi'o finished parsing kaiju taxonomy output. Bad news: it has no idea whether "
            "it did well or not. Because the user can ask kaiju to report certain taxonomic levels, but "
            "can't ask anvi'o to utilize that information. So anvi'o always assumes you started from "
            "the domain-level, and followed the conventional levels of taxonomy. Here is your question. We "
            "randomly picked some phylum names from your input taoxnomy as anvi'o parsed them. Here they are: "
            "'%s'. Do they look like phylum names to you? If they don't, you are in very big trouble :( The "
            "best way to get yourself out of trouble is to immediately press CTRL-C, turn your computer off, "
            "and move permanently to Cuba since you are done here (if you are already in Cuba, please let "
            "us know for more instructions)." % ", ".join(random_phylum_names),
            header="Good news, bad news, and a question",
            lc="yellow",
        )

        return TaxonomyHelper(taxonomy_dict).get_genes_taxonomy_and_taxon_names_dicts()
