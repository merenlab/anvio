# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Common functions and classes for SCG/TRNA taxonomy.
"""

import os
import gzip
import pandas as pd
import multiprocessing

from collections import OrderedDict

import anvio
import anvio.tables as t
import anvio.utils as utils
import anvio.hmmops as hmmops
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.drivers.blast import BLAST
from anvio.dbops import ContigsSuperclass
from anvio.drivers.diamond import Diamond
from anvio.tables.scgtaxonomy import TableForSCGTaxonomy
from anvio.tables.trnataxonomy import TableForTRNATaxonomy

__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


run_quiet = terminal.Run(log_file_path=None, verbose=False)
progress_quiet = terminal.Progress(verbose=False)
pp = terminal.pretty_print

proper_foci = ['trnas', 'scgs']


class XXXTaxonomyExampleBaseClass(object):
    def __init__(self, args):
        """An example base class to fill in common arguments for SCG/TRNA Taxonomy classes."""

        if not hasattr(self, 'focus'):
            raise ConfigError("You are lost :/ The class that initializes `XXXTaxonomyArgs` must have a `focus` "
                              "self variable. Focus can be one of these: %s." % ', '.join(proper_foci))

        if self.focus not in proper_foci:
            raise ConfigError("Unknown focus %s :/" % (self.focus))

        self.scgs_focus = self.focus == 'scgs'
        self.trna_focus = self.focus == 'trnas'

        if self.scgs_focus:
            pass

        if self.trna_focus:
            pass


class PopulateContigsDatabaseWithXXXTaxonomy(object):
    def __init__(self, args):
        if not hasattr(self, 'focus'):
            raise ConfigError("You are lost :/ The class that initializes `PopulateContigsDatabaseWithXXXTaxonomy` must have a `focus` "
                              "self variable. Focus can be one of these: %s." % ', '.join(proper_foci))

        if self.focus not in proper_foci:
            raise ConfigError("Unknown focus %s :/" % (self.focus))

        self.scgs_focus = self.focus == 'scgs'
        self.trna_focus = self.focus == 'trnas'

        self._DELIV = "SCG taxonomy" if self.scgs_focus else "tRNA taxonomy"
        self._ITEMS = "SCGs" if self.scgs_focus else 'anticodons'
        self._DATA = "single-copy core gene sequences" if self.scgs_focus else "tRNA gene sequences"
        self._SETIP = "anvi-setup-scg-taxonomy" if self.scgs_focus else "anvi-setup-trna-taxonomy"

        if self.scgs_focus:
            pass

        if self.trna_focus:
            pass

        self.taxonomy_dict = OrderedDict()

        self.mutex = multiprocessing.Lock()


    def get_sequences_dict_from_contigs_db(self):
        """Returns a dictionary of all HMM hits per SCG of interest"""

        contigs_db = ContigsSuperclass(self.args, r=run_quiet, p=progress_quiet)
        splits_dict = {contigs_db.a_meta['project_name']: list(contigs_db.splits_basic_info.keys())}

        if self.scgs_focus:
            hmm_source = self.ctx.hmm_source_for_scg_taxonomy
            default_items_for_taxonomy = self.ctx.default_scgs_for_taxonomy
            return_amino_acid_sequences = True
        elif self.trna_focus:
            hmm_source = self.ctx.hmm_source_for_trna_genes
            default_items_for_taxonomy = self.ctx.default_anticodons_for_taxonomy
            return_amino_acid_sequences = False
        else:
            raise ConfigError("Anvi'o is lost.")

        s = hmmops.SequencesForHMMHits(self.args.contigs_db, sources=hmm_source, run=run_quiet, progress=progress_quiet)
        hmm_sequences_dict = s.get_sequences_dict_for_hmm_hits_in_splits(splits_dict, return_amino_acid_sequences=return_amino_acid_sequences)

        # a trick for trnas specifically since tRNA gene names contain two features and formed as `AMINOACID_ANTICODON`.
        # here we replace those gene names with `ANTICODON` only
        if self.trna_focus:
            for entry in hmm_sequences_dict:
                hmm_sequences_dict[entry]['gene_name'] = hmm_sequences_dict[entry]['gene_name'].split('_')[1]

        hmm_sequences_dict = utils.get_filtered_dict(hmm_sequences_dict, 'gene_name', set(default_items_for_taxonomy))

        if not len(hmm_sequences_dict):
            return None

        self.progress.reset()
        self.run.info(f'Num relevant {self._ITEMS} in contigs db', '%s' % (pp(len(hmm_sequences_dict))))

        item_sequences_dict = {}
        for entry_id in hmm_sequences_dict:
            entry = hmm_sequences_dict[entry_id]

            item_name = entry['gene_name']
            if item_name in item_sequences_dict:
                item_sequences_dict[item_name][entry_id] = entry
            else:
                item_sequences_dict[item_name] = {entry_id: entry}

        return item_sequences_dict


    def populate_contigs_database(self):
        """Populates SCG/tRNA taxonomy tables in a contigs database"""

        # get an instnce for the tables for taxonomy early on.
        if self.scgs_focus:
            anvio_taxonomy_table_name = t.scg_taxonomy_table_name
            self.tables_for_taxonomy = TableForSCGTaxonomy(self.contigs_db_path, self.run, self.progress)
            database_version = self.ctx.scg_taxonomy_database_version
        elif self.trna_focus:
            anvio_taxonomy_table_name = t.trna_taxonomy_table_name
            self.tables_for_taxonomy = TableForTRNATaxonomy(self.contigs_db_path, self.run, self.progress)
            database_version = self.ctx.trna_taxonomy_database_version
        else:
            anvio_taxonomy_table_name = None
            self.tables_for_taxonomy = None
            database_version = None

        # get the dictionary that shows all hits for each SCG of interest
        self.progress.new('Contigs bleep bloop')
        self.progress.update(f'Recovering the {self._ITEMS} dictionary')
        item_sequences_dict = self.get_sequences_dict_from_contigs_db()
        self.progress.end()

        if not item_sequences_dict:
            self.run.warning(f"This contigs database contains no {self._DATA} that are used by the "
                              "anvi'o taxonomy headquarters in Lausanne. Somewhat disappointing but totally OK.")

            # even if there are no SCGs to use for taxonomy later, we did attempt ot populate the
            # contigs database, so we shall note that in the self table to make sure the error from
            # `anvi-estimate-genome-taxonomy` is not "you seem to have not run taxonomy".
            self.tables_for_taxonomy.update_db_self_table_values(taxonomy_was_run=True, database_version=database_version)

            # return empty handed like a goose in the job market in 2020
            return None

        log_file_path = filesnpaths.get_temp_file_path()

        self.run.info('Taxonomy', self.ctx.accession_to_taxonomy_file_path)
        self.run.info('Database reference', self.ctx.search_databases_dir_path)
        self.run.info(f'Number of {self._ITEMS}', len(item_sequences_dict))

        self.run.warning('', header='Parameters for search', lc='green')
        self.run.info('Max number of target sequences', self.max_target_seqs)
        self.run.info('Max e-value to report alignments', self.evalue)
        self.run.info('Min percent identity to report alignments', self.min_pct_id)
        self.run.info('Num aligment tasks running in parallel', self.num_parallel_processes)
        self.run.info('Num CPUs per aligment task', self.num_threads)
        self.run.info('Log file path', log_file_path)

        self.tables_for_taxonomy.delete_contents_of_table(anvio_taxonomy_table_name, warning=False)
        self.tables_for_taxonomy.update_db_self_table_values(taxonomy_was_run=False, database_version=None)

        total_num_processes = len(item_sequences_dict)

        self.progress.new('Performing search', progress_total_items=total_num_processes)
        self.progress.update('Initializing %d process...' % int(self.num_parallel_processes))

        manager = multiprocessing.Manager()
        input_queue = manager.Queue()
        output_queue = manager.Queue()
        error_queue = manager.Queue()

        search_output = []

        for item_name in item_sequences_dict:
            sequence = ""
            for entry in item_sequences_dict[item_name].values():
                if 'sequence' not in entry or 'gene_name' not in entry:
                    raise ConfigError("The `get_filtered_dict` function got a parameter that "
                                      "does not look like the way we expected it. This function "
                                      "expects a dictionary that contains keys `gene_name` and `sequence`.")

                sequence = sequence + ">" + str(entry['gene_callers_id']) + "\n" + entry['sequence'] + "\n"
                entry['hits'] = []

            input_queue.put([item_name, sequence])

        workers = []
        for i in range(0, int(self.num_parallel_processes)):
            worker = multiprocessing.Process(target=self.blast_search_worker, args=(input_queue, output_queue, error_queue, log_file_path))

            workers.append(worker)
            worker.start()

        num_finished_processes = 0
        while num_finished_processes < total_num_processes:
            # check error
            error_text = error_queue.get()
            if error_text:
                self.progress.reset()

                for worker in workers:
                    worker.terminate()

                if 'incompatible' in error_text:
                    raise ConfigError(f"Your current databases are incompatible with the diamond version you have on your computer. "
                                       "Please run the command `{self._SETUP} --redo-databases` and come back.")
                else:
                    if self.scgs_focus:
                        raise ConfigError("Bad news. The database search operation failed somewhere :( It is very hard for anvi'o "
                                          "to know what happened, but the MOST LIKELY reason is that you have a diamond version "
                                          "installed on your system that is incompatible with anvi'o :/ The best course of action for that "
                                          "is to make sure running `diamond --version` on your terminal returns `0.9.14`. If not, "
                                          "try to upgrade/downgrade your diamond to match this version. If you are in a conda environmnet "
                                          "you can try running `conda install diamond=0.9.14`. Please feel free to contact us if the problem "
                                          "persists. We apologize for the inconvenience.")
                    else:
                        raise ConfigError("Bad news. The database search operation failed somewhere :( It is very hard for anvi'o "
                                          "to know what happened, but the MOST LIKELY reason is that you have a BLAST version "
                                          "installed on your system that is incompatible with anvi'o :/ If you see `2.6.0` or higher "
                                          "version numbers when you type `blastn -version` in your terminal, it may be wortwhile to "
                                          "get in touch with anvi'o developers.")

            try:
                search_output += output_queue.get()

                if self.write_buffer_size > 0 and len(search_output) % self.write_buffer_size == 0:
                    self.tables_for_taxonomy.add(search_output)
                    search_output = []

                num_finished_processes += 1

                self.progress.increment(increment_to=num_finished_processes)
                self.progress.update(f"%s of %s {self._ITEMS} are finished in %s processes with %s threads." \
                                        % (num_finished_processes, total_num_processes, int(self.num_parallel_processes), self.num_threads))

            except KeyboardInterrupt:
                print("Anvi'o profiler recieved SIGINT, terminating all processes...")
                break

        for worker in workers:
            worker.terminate()

        # finally the remaining hits are written to the database, and we are done
        self.tables_for_taxonomy.add(search_output)

        # time to update the self table:
        self.tables_for_taxonomy.update_db_self_table_values(taxonomy_was_run=True, database_version=database_version)

        self.progress.end()


    def show_hits_gene_callers_id(self, gene_callers_id, item_name, hits):
        self.progress.reset()
        self.run.warning(None, header='Hits for gene caller id %s' % gene_callers_id, lc="green")

        if len(hits):
            header = ['%id', 'bitscore', 'accession', 'taxonomy']
            table = []

            self.run.info_single("For '%s'" % item_name, nl_before=1, nl_after=1)

            for hit in hits:
                table.append([str(hit['percent_identity']), str(hit['bitscore']), hit['accession'], ' / '.join([hit[l] if hit[l] else '' for l in self.ctx.levels_of_taxonomy])])

            anvio.TABULATE(table, header)
        else:
            self.run.info_single("No hits :/")


    def update_dict_with_taxonomy(self, d, mode=None):
        """Takes a dictionary that includes a key `accession` and populates the dictionary with taxonomy"""

        if not mode:
            if not 'accession' in d:
                raise ConfigError("`add_taxonomy_to_dict` is speaking: the dictionary sent here does not have a member "
                                  "with key `accession`.")

            if d['accession'] in self.ctx.accession_to_taxonomy_dict:
                d.update(self.ctx.accession_to_taxonomy_dict[d['accession']])
            else:
                d.update(self.ctx.accession_to_taxonomy_dict['unknown_accession'])

        elif mode == 'list_of_dicts':
            if len([entry for entry in d if 'accession' not in entry]):
                raise ConfigError("`add_taxonomy_to_dict` is speaking: you have a bad formatted data here :/")

            for entry in d:
                print(self.taxonomy_dict[entry['accession']])

        else:
            raise ConfigError("An unknown mode (%s) is set to `add_taxonomy_to_dict` :/" % (mode))

        return d


    def blast_search_worker(self, input_queue, output_queue, error_queue, log_file_path):
        """BLAST each SCG identified in the contigs database against the corresopinding
           target local database of GTDB seqeunces
        """

        while True:
            item_name, fasta_formatted_sequence = input_queue.get(True)

            if self.scgs_focus:
                target_database_path = self.ctx.SCGs[item_name]['db']

                diamond = Diamond(target_database_path, run=run_quiet, progress=progress_quiet)
                diamond.max_target_seqs = self.max_target_seqs
                diamond.evalue = self.evalue
                diamond.min_pct_id = self.min_pct_id
                diamond.num_threads = self.num_threads
                diamond.run.log_file_path = log_file_path

                raw_search_output = diamond.blastp_stdin_multi(fasta_formatted_sequence)
            elif self.trna_focus:
                # FIXME FIXME FIXME FIXME FIXME
                # the following line must be deleted. but we can't delete it unless we first fix the
                # anvi-scan-trna's code. curently it uses gene names as aminoacid_codon, when it should
                # use aminoacid_anticodon.
                item_name = utils.rev_comp(item_name)

                target_database_path = self.ctx.anticodons[item_name]['db']

                blast = BLAST(None, target_database_path, run=run_quiet, progress=progress_quiet)
                blast.search_program = 'blastn'
                blast.max_target_seqs = self.max_target_seqs
                blast.evalue = self.evalue
                blast.min_pct_id = self.min_pct_id
                blast.num_threads = self.num_threads
                blast.run.log_file_path = log_file_path

                raw_search_output = blast.blast_stdin(fasta_formatted_sequence)
            else:
                raise ConfigError("You must be joking, Mr. Feynman.")

            hits_per_gene = {}
            genes_estimation_output=[]

            for blastp_hit in raw_search_output.split('\n'):
                if len(blastp_hit) and not blastp_hit.startswith('Query'):
                    fields = blastp_hit.split('\t')

                    try:
                        gene_callers_id = int(fields[0])
                        error_queue.put(None)
                    except:
                        error_queue.put(raw_search_output)

                    hit = dict(zip(['accession', 'percent_identity', 'bitscore'], [fields[1], float(fields[2]), float(fields[11])]))
                    hit = self.update_dict_with_taxonomy(hit)

                    if gene_callers_id not in hits_per_gene:
                        hits_per_gene[gene_callers_id] = {}

                    if item_name not in hits_per_gene[gene_callers_id]:
                        hits_per_gene[gene_callers_id][item_name] = []

                    hits_per_gene[gene_callers_id][item_name].append(hit)
                else:
                    error_queue.put(None)

            for gene_callers_id, raw_hits in hits_per_gene.items():
                if len(raw_hits.keys()) > 1:
                    self.run.warning("As crazy as it sounds, the gene callers id `%d` seems to have hit more than one SCG o_O Anvi'o will only use "
                                     "one of them almost absolutely randomly. Here are the SCGs the gene sequence matches: '%s'" % [s for s in raw_hits.keys()])

                item_name = list(raw_hits.keys())[0]
                raw_hits = raw_hits[item_name]

                onsensus_hit = self.get_consensus_hit(raw_hits)
                onsensus_hit['accession'] = 'CONSENSUS'

                if anvio.DEBUG:
                    # avoid race conditions when priting this information when `--debug` is true:
                    with self.mutex:
                        self.progress.reset()
                        self.show_hits_gene_callers_id(gene_callers_id, item_name, raw_hits + [onsensus_hit])

                genes_estimation_output.append([gene_callers_id, item_name, [onsensus_hit]])

            output_queue.put(genes_estimation_output)


    def get_consensus_hit(self, raw_hits):
        pd.set_option('mode.chained_assignment', None)

        df = pd.DataFrame.from_records(raw_hits)

        # remove hits that are null at the phylum level if there are still hits
        # in the df that are not null:
        not_null_hits = df[df.t_phylum.notnull()]
        if len(not_null_hits):
            df = not_null_hits

        # find the max percent identity score in the df
        max_percent_identity = max(df['percent_identity'])

        # subset the data frame to those with percent identity that match to `max_percent_identity`
        df_max_identity = df.loc[df.percent_identity == max_percent_identity]

        # if some of the competing names have null species deignations, remove them from consideration
        if len(df_max_identity.t_species.unique()) > 1:
            df_max_identity = df_max_identity[df_max_identity.t_species.notnull()]

        # find the taxonomic level where the number of unique taxon names is one
        for taxonomic_level in self.ctx.levels_of_taxonomy[::-1]:
            if len(df_max_identity[taxonomic_level].unique()) == 1:
                break

        # take one of the hits from `df_max_identity`, and assign None to all taxonomic levels
        # beyond `taxonomic_level`, which, after the loop above shows the proper level of
        # assignment for this set
        final_hit = df_max_identity.head(1)
        for taxonomic_level_to_nullify in self.ctx.levels_of_taxonomy[self.ctx.levels_of_taxonomy.index(taxonomic_level) + 1:]:
            final_hit.at[0, taxonomic_level_to_nullify] = None

        # FIXME: final hit is still not what we can trust. next, we should find out whether the percent identity
        # for the level of taxonomy at `taxonomic_level` is higher than the minimum percent identity for all sequences
        # considered that are affiliated with final_hit[taxonomic_level]

        # turn it into a Python dict before returning
        final_hit_dict = final_hit.to_dict('records')[0]

        return final_hit_dict


def get_accession_to_taxonomy_dict(accession_to_taxonomy_file_path, levels_of_taxonomy, progress, run):
    if not os.path.exists(accession_to_taxonomy_file_path):
        return None

    letter_to_level = dict([(l.split('_')[1][0], l) for l in levels_of_taxonomy])

    progress.new("Reading the accession to taxonomy file")
    progress.update('...')

    accession_to_taxonomy_dict = {}
    with gzip.open(accession_to_taxonomy_file_path, 'rb') as taxonomy_file:
        for line in taxonomy_file.readlines():
            line = line.decode('utf-8')

            if line.startswith('#'):
                continue

            accession, taxonomy_text = line.strip('\n').split('\t')
            # taxonomy_text kinda looks like these:
            #
            #    d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Alcaligenes;s__Alcaligenes faecalis_C
            #    d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Enterococcaceae;g__Enterococcus_B;s__Enterococcus_B faecalis
            #    d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Moraxellaceae;g__Acinetobacter;s__Acinetobacter sp1
            #    d__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Bacillaceae_G;g__Bacillus_A;s__Bacillus_A cereus_AU
            #    d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Tissierellales;f__Helcococcaceae;g__Finegoldia;s__Finegoldia magna_H

            d = {}
            for letter, taxon in [e.split('__', 1) for e in taxonomy_text.split(';')]:
                if letter in letter_to_level:
                    # NOTE: This is VERY important. Here we are basically removing subclades GTDB defines for
                    # simplicity. We may have to change this behavior later. So basically, Enterococcus_B will
                    # become Enterococcus
                    if '_' in taxon:
                        if letter != 's':
                            d[letter_to_level[letter]] = '_'.join(taxon.split('_')[:-1])
                        else:
                            # special treatment for species level taxonomy string.
                            # the genus is copied for the species level taxonomy, such as this one, 'Bacillus_A cereus', or
                            # species itmay have a subclade, such as this one, 'Corynebacterium aurimucosum_C', so we
                            # neeed to make sure the subclades are removed from all words in the species level
                            # taxonomy string.
                            d[letter_to_level[letter]] = ' '.join(['_'.join(word.split('_')[:-1]) if '_' in word else word for word in taxon.split(' ')])
                    else:
                        d[letter_to_level[letter]] = taxon
                else:
                    run.warning("Some weird letter found in '%s' :(" % taxonomy_text)

            accession_to_taxonomy_dict[accession] = d

    # let's add one more accession for all those missing accessions
    accession_to_taxonomy_dict['unknown_accession'] = dict([(taxon, None) for taxon in levels_of_taxonomy])

    progress.end()

    return accession_to_taxonomy_dict
