# pylint: disable=line-too-long
"""Classes to work with functional annotations stored in anvi'o contigs databases"""

import sys
import argparse
import textwrap

import anvio
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import anvio.ccollections as ccollections

from anvio.errors import ConfigError
from anvio.genomedescriptions import GenomeDescriptions


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print
P = terminal.pluralize


class SequencesForFunctionHits(object):
    """Get sequences for genes with hits in a single functional annotation source.

    This class resolves every supported input mode (a single contigs-db, bins in a collection,
    external genomes, internal genomes) into a list of 'units' where each unit is a tuple of
    `(name, contigs_db_path, split_names_of_interest)`, and runs a single code path over them.
    The last item is `None` when all genes in the contigs-db belong to the unit. It is similar
    to SequencesForHMMHits in anvio/hmmops.py in spirit, but it is much simpler as this class
    does not offer any concatenation or alignment of sequences on purpose: genes annotated with
    the same function can share extremely low sequence identity, and at the time of writing this
    it made no sense to make their sequences to be used for phylogenomics as if they were homologs.

    If that changes, and someone wishes to make use of function hits for phylogenomics, instead of
    extending this class, they should consider creating a base class with the features of
    SequencesForHMMHits and then create two subclasses for HMM hits and function hits. Otherwise
    it will be a mess to maintain and extend this in the future as there will be a crap ton of
    redundancy.
    """

    # variables users can use in their `--defline-format` strings as this differs a bit from
    # other use cases.
    DEFLINE_VARIABLES = ['gene_callers_id', 'genome_name', 'source', 'accession', 'function',
                         'contig_name', 'start', 'stop', 'length', 'direction']

    # the default defline. the first token is unique across the entire output (it look ugly but
    # but it be comprehensive)
    DEFAULT_DEFLINE_FORMAT = "{genome_name}_{gene_callers_id} source:{source}|accession:{accession}|function:{function}"

    # how multiple accessions or functions of a single gene are joined in deflines. the never ending
    # torture of exclamation marks as a separator in anvi'o :*
    MULTI_VALUE_SEPARATOR = '!!!'


    def __init__(self, args, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.contigs_db_path = A('contigs_db')
        self.profile_db_path = A('profile_db')
        self.collection_name = A('collection_name')
        self.bin_id = A('bin_id')
        self.bin_ids_file_path = A('bin_ids_file')
        self.external_genomes_path = A('external_genomes')
        self.internal_genomes_path = A('internal_genomes')

        self.annotation_source = A('annotation_source')
        self.accession_ids_input = A('accession_ids')
        self.output_file_path = A('output_file')
        self.report_aa_sequences = A('get_aa_sequences')
        self.wrap = None if A('no_wrap') else 120
        self.defline_format = A('defline_format') or self.DEFAULT_DEFLINE_FORMAT

        # will be filled during self.sanity_check(), so make sure to wire your code
        # through `self.process()` if you use this class in the future.
        self.accession_ids = []

        self.check_input_combinations()


    @classmethod
    def list_defline_variables(cls, run=run):
        """Show the variables that can be used in defline format strings, and the default format"""

        run.warning("Here are the variables you can use to provide a user-defined defline template: ")
        for key in cls.DEFLINE_VARIABLES:
            run.info_single("{%s}" % key)

        run.info_single("Remember, by default, anvi'o will use the following template to format the deflines of "
                        "the FASTA file it produces:", level=0, nl_before=1, nl_after=1, mc='red')
        run.info_single(cls.DEFAULT_DEFLINE_FORMAT, level=0, nl_after=1, cut_after=0)


    def check_input_combinations(self):
        """Make sure the input files make sense together. Bad combinations are errors, not warnings."""

        genomes_files_used = self.external_genomes_path or self.internal_genomes_path

        if not (genomes_files_used or self.contigs_db_path or self.profile_db_path):
            raise ConfigError("You gotta give this program some input files :/ Please provide either a contigs-db (with or "
                              "without a profile-db and a collection), or an external and/or internal genomes file.")

        if genomes_files_used and (self.contigs_db_path or self.profile_db_path):
            raise ConfigError("You can either work with external/internal genomes files, or with a contigs-db (and optionally "
                              "a profile-db). Not both :/")

        if genomes_files_used and (self.collection_name or self.bin_id or self.bin_ids_file_path):
            raise ConfigError("Collection and bin parameters are only relevant when you are working with a contigs-db and a "
                              "profile-db. For internal genomes, collections and bins are already described in the internal "
                              "genomes file.")

        if self.profile_db_path and not self.contigs_db_path:
            raise ConfigError("A profile-db is not useful without the contigs-db it belongs to. Please provide one with "
                              "`--contigs-db`.")

        if self.profile_db_path and not self.collection_name:
            raise ConfigError("If you provide a profile-db, you must also provide a collection name with `--collection-name`. "
                              "Otherwise anvi'o has no idea what you want from the profile-db.")

        if self.collection_name and not self.profile_db_path:
            raise ConfigError("You can't use a collection name without a profile-db :/ Please provide one with `--profile-db`.")

        if (self.bin_id or self.bin_ids_file_path) and not self.collection_name:
            raise ConfigError("Bin parameters are only meaningful together with a collection name (and a profile-db).")

        if self.bin_id and self.bin_ids_file_path:
            raise ConfigError("Either use a file to list all the bin ids (`--bin-ids-file`), or declare a single bin "
                              "(`--bin-id`) you would like to focus. Not both :/")


    def sanity_check(self):
        """Checks that are only relevant when we are about to do actual work (and not just listing things)"""

        if not self.annotation_source:
            raise ConfigError("You must declare a single functional annotation source with `--annotation-source`. You can "
                              "use the flag `--list-annotation-sources` to see which ones are available (or simply use "
                              "`anvi-db-info` on your contigs-db).")

        if ',' in self.annotation_source:
            raise ConfigError(f"This program wants exactly one functional annotation source, but your input "
                              f"('{self.annotation_source}') looks like a list of them :/")

        if not self.output_file_path:
            raise ConfigError("Please don't make anvi'o make up output file names for you by providing an output "
                              ("file path with with the parameter `--output-file`.")

        filesnpaths.is_output_file_writable(self.output_file_path)

        # make sure the defline format is composed of known variables
        utils.get_f_string_evaluated_by_dict(self.defline_format, dict([(v, None) for v in self.DEFLINE_VARIABLES]))

        self.accession_ids = self.get_accession_ids_from_args()


    def get_accession_ids_from_args(self):
        """Returns a list of unique accession ids from a comma-separated string or a file with one id per line"""

        if not self.accession_ids_input:
            return []

        if filesnpaths.is_file_exists(self.accession_ids_input, dont_raise=True):
            accession_ids = [a.strip() for a in open(self.accession_ids_input).readlines()]
        else:
            accession_ids = [a.strip() for a in self.accession_ids_input.split(',')]

        # remove empty ones and redundancy, but keep the order
        accession_ids = list(dict.fromkeys([a for a in accession_ids if a]))

        if not len(accession_ids):
            raise ConfigError("Anvi'o could not find any accession ids in your `--accession-ids` input :/")

        return accession_ids


    def get_units(self):
        """Resolve all input modes into a list of `(name, contigs_db_path, split_names_of_interest)` tuples"""

        units = []

        if self.external_genomes_path or self.internal_genomes_path:
            genome_descriptions = GenomeDescriptions(self.args, run=self.run, progress=self.progress)
            genome_descriptions.load_genomes_descriptions(skip_functions=True, init=False)

            for genome_name in genome_descriptions.genomes:
                entry = genome_descriptions.genomes[genome_name]

                if genome_name in genome_descriptions.internal_genome_names:
                    split_names = set(genome_descriptions.get_split_names_of_interest_for_internal_genome(entry))
                else:
                    split_names = None

                units.append((genome_name, entry['contigs_db_path'], split_names))
        elif self.profile_db_path:
            utils.is_profile_db_and_contigs_db_compatible(self.profile_db_path, self.contigs_db_path)

            splits_dict = ccollections.GetSplitNamesInBins(self.args).get_dict()

            for bin_name in splits_dict:
                units.append((bin_name, self.contigs_db_path, set(splits_dict[bin_name])))
        else:
            contigs_db = dbops.ContigsDatabase(self.contigs_db_path, run=self.run, progress=self.progress)
            project_name = contigs_db.meta['project_name_str']
            contigs_db.disconnect()

            units.append((project_name, self.contigs_db_path, None))

        return units


    def list_annotation_sources(self):
        """List functional annotation sources available in the input databases"""

        if self.external_genomes_path or self.internal_genomes_path:
            genome_descriptions = GenomeDescriptions(self.args, run=self.run, progress=self.progress)
            genome_descriptions.load_genomes_descriptions(init=False)
            genome_descriptions.list_function_sources()
        else:
            dbops.ContigsDatabase(self.contigs_db_path, run=self.run, progress=self.progress).list_function_sources()


    def check_annotation_source_in_contigs_dbs(self, units):
        """Make sure the annotation source is in every contigs-db we will work with"""

        unit_names_missing_the_source = []

        for unit_name, contigs_db_path, _ in units:
            contigs_db = dbops.ContigsDatabase(contigs_db_path, run=self.run, progress=self.progress)
            sources_in_db = contigs_db.meta['gene_function_sources'] or []
            contigs_db.disconnect()

            if self.annotation_source not in sources_in_db:
                unit_names_missing_the_source.append(unit_name)

        if len(unit_names_missing_the_source):
            raise ConfigError(f"The functional annotation source '{self.annotation_source}' is missing from the contigs-db "
                              f"of {P('input', len(unit_names_missing_the_source))}: {', '.join(unit_names_missing_the_source)}. "
                              f"You can use the flag `--list-annotation-sources` to see what is available.")


    def get_hits_for_contigs_db(self, contigs_db, units):
        """Returns hits in the annotation source for each unit that shares a single contigs-db.

        Parameters
        ==========
        contigs_db : dbops.ContigsSuperclass
            An instance initialized for the contigs-db all `units` share.
        units : list
            A list of `(name, contigs_db_path, split_names_of_interest)` tuples.

        Returns
        =======
        hits : dict
            A dictionary of `{unit_name: {gene_callers_id: [(accession, function), ...]}}`. Each gene
            appears once per unit with a non-redundant list of its accession/function pairs.
        """

        # we want ALL hits from this source for each gene, and not only the best one, so we set this
        # global for the duration of the call and restore it right after.
        original_value = anvio.RETURN_ALL_FUNCTIONS_FROM_SOURCE_FOR_EACH_GENE
        anvio.RETURN_ALL_FUNCTIONS_FROM_SOURCE_FOR_EACH_GENE = True
        try:
            contigs_db.init_functions(requested_sources=[self.annotation_source])
        finally:
            anvio.RETURN_ALL_FUNCTIONS_FROM_SOURCE_FOR_EACH_GENE = original_value

        # collapse multiple rows and `!!!`-joined values into a non-redundant list of
        # (accession, function) pairs for each gene
        hits_per_gene = {}
        for gene_callers_id, sources_dict in contigs_db.gene_function_calls_dict.items():
            entry = sources_dict.get(self.annotation_source)
            if not entry:
                continue

            accessions = str(entry[0] or '').split('!!!')
            functions = str(entry[1] or '').split('!!!')

            # these two lists should always be of the same length, but let's not crash if they are not
            if len(accessions) != len(functions):
                accessions = accessions + [''] * (len(functions) - len(accessions))
                functions = functions + [''] * (len(accessions) - len(functions))

            pairs = list(dict.fromkeys([(a.strip(), f.strip()) for a, f in zip(accessions, functions)]))

            if self.accession_ids and not any([a in self.accession_ids for a, _ in pairs]):
                continue

            hits_per_gene[gene_callers_id] = pairs

        # assign genes to units
        hits = {}
        for unit_name, _, split_names in units:
            if split_names is not None:
                gene_to_split = contigs_db.gene_callers_id_to_split_name_dict
                hits[unit_name] = dict([(g, p) for g, p in hits_per_gene.items() if gene_to_split.get(g) in split_names])
            else:
                hits[unit_name] = dict(hits_per_gene)

        return hits


    def process(self):
        """Do everything and write the FASTA file"""

        self.sanity_check()

        units = self.get_units()

        self.check_annotation_source_in_contigs_dbs(units)

        # group units by contigs-db so each contigs-db is read only once, no matter
        # how many bins or genomes it describes
        units_per_contigs_db = {}
        for unit in units:
            units_per_contigs_db.setdefault(unit[1], []).append(unit)

        # here we will collect the final entries for the FASTA file, in the order of units
        entries = dict([(unit_name, []) for unit_name, _, _ in units])
        accession_ids_found = set([])
        genes_without_sequences = []

        for contigs_db_path, units_of_contigs_db in units_per_contigs_db.items():
            contigs_db = dbops.ContigsSuperclass(argparse.Namespace(contigs_db=contigs_db_path), r=self.run, p=self.progress)

            hits = self.get_hits_for_contigs_db(contigs_db, units_of_contigs_db)

            gene_caller_ids = sorted(set([g for unit_hits in hits.values() for g in unit_hits]))

            # nothing to get sequences for. NOTE: we must not call `get_sequences_for_gene_callers_ids`
            # with an empty list, since it would then return every gene in the contigs-db.
            if not len(gene_caller_ids):
                continue

            _, sequences_dict = contigs_db.get_sequences_for_gene_callers_ids(gene_caller_ids_list=gene_caller_ids,
                                                                              report_aa_sequences=self.report_aa_sequences,
                                                                              simple_headers=True)

            for unit_name, _, _ in units_of_contigs_db:
                for gene_callers_id in sorted(hits[unit_name]):
                    pairs = hits[unit_name][gene_callers_id]
                    gene_call = sequences_dict[gene_callers_id]

                    accession_ids_found.update([a for a, _ in pairs if a in self.accession_ids])

                    sequence = gene_call['aa_sequence'] if self.report_aa_sequences else gene_call['sequence']
                    if not sequence:
                        genes_without_sequences.append(f"{unit_name}:{gene_callers_id}")
                        continue

                    defline_data_dict = {'gene_callers_id': gene_callers_id,
                                         'genome_name': unit_name,
                                         'source': self.annotation_source,
                                         'accession': self.MULTI_VALUE_SEPARATOR.join([a for a, _ in pairs]),
                                         'function': self.MULTI_VALUE_SEPARATOR.join([f for _, f in pairs]),
                                         'contig_name': gene_call['contig'],
                                         'start': gene_call['start'],
                                         'stop': gene_call['stop'],
                                         'length': gene_call['length'],
                                         'direction': gene_call['direction']}

                    header = utils.get_f_string_evaluated_by_dict(self.defline_format, defline_data_dict)
                    entries[unit_name].append((header, sequence))

        # did we find every accession id the user asked for in at least one input?
        if self.accession_ids:
            missing_accession_ids = [a for a in self.accession_ids if a not in accession_ids_found]
            if len(missing_accession_ids):
                raise ConfigError(f"Some of the accession ids you requested did not match to any gene in any of your "
                                  f"inputs for the annotation source '{self.annotation_source}' o_O Anvi'o does not want to "
                                  f"silently ignore them, so here they are: {', '.join(missing_accession_ids)}.")

        if self.report_aa_sequences and len(genes_without_sequences):
            if not sum([len(e) for e in entries.values()]):
                raise ConfigError("You requested amino acid sequences with `--get-aa-sequences`, but none of the genes with "
                                  "hits have amino acid sequences associated with them :/ Which is weird and should never "
                                  "happen under normal circumstances. But it is difficult the code here to figure out what "
                                  "might be the reason for that. As the user, you can get the DNA sequences of these hits by "
                                  "dropping that flag.")
            else:
                self.run.warning(f"{P('gene', len(genes_without_sequences))} with hits had no amino acid sequences, and will "
                                 f"not be in the output FASTA file. Here they are (as genome_name:gene_callers_id): "
                                 f"{', '.join(genes_without_sequences)}.")

        num_entries = sum([len(e) for e in entries.values()])
        if not num_entries:
            raise ConfigError(f"None of the genes in your input(s) had a hit in the annotation source '{self.annotation_source}'"
                              f"{' that matched to your accession ids' if self.accession_ids else ''}. There is nothing to "
                              f"report. #SAD.")

        # make sure deflines are unique across the entire output
        all_headers = [header for unit_name in entries for header, _ in entries[unit_name]]
        if len(all_headers) != len(set(all_headers)):
            seen, duplicates = set([]), []
            for header in all_headers:
                if header in seen and header not in duplicates:
                    duplicates.append(header)
                seen.add(header)

            raise ConfigError(f"Your defline format ('{self.defline_format}') results in {P('defline', len(duplicates))} that "
                              f"occur more than once in the output.. Please include variables in your format that make each "
                              f"entry unique (such as `{{genome_name}}` and `{{gene_callers_id}}`). Here are some of the "
                              f"duplicates: {', '.join(duplicates[:5])}.")

        # write the output
        with open(self.output_file_path, 'w') as output:
            for unit_name in entries:
                for header, sequence in entries[unit_name]:
                    if self.wrap:
                        sequence = textwrap.fill(sequence, self.wrap, break_on_hyphens=False)

                    output.write(f">{header}\n{sequence}\n")

        # report
        self.run.info('Annotation source', self.annotation_source, nl_before=1)
        self.run.info('Accession ids requested', len(self.accession_ids) if self.accession_ids else 'None (all hits are reported)')
        self.run.info('Sequence type', 'AA sequences' if self.report_aa_sequences else 'DNA sequences', mc='green')
        for unit_name in entries:
            self.run.info(f"Genes reported for '{unit_name}'", len(entries[unit_name]))
        self.run.info('Output FASTA', self.output_file_path, nl_before=1)
