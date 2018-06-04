# -*- coding: utf-8
# pylint: disable=line-too-long

import os
import hashlib

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.drivers.hmmer import HMMer
from anvio.tables.tableops import Table
from anvio.parsers import parser_modules
from anvio.dbops import ContigsSuperclass
from anvio.tables.genecalls import TablesForGeneCalls


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


class TablesForHMMHits(Table):
    def __init__(self, db_path, num_threads_to_use=1, run=run, progress=progress):
        self.num_threads_to_use = num_threads_to_use
        self.db_path = db_path

        utils.is_contigs_db(self.db_path)

        Table.__init__(self, self.db_path, anvio.__contigs__version__, run, progress)

        if not self.genes_are_called:
            raise ConfigError("It seems the contigs database '%s' was created with '--skip-gene-calling' flag.\
                                Nothing to do here :/" % (self.db_path))

        self.init_gene_calls_dict()

        if not len(self.gene_calls_dict):
            raise ConfigError("Tables that should contain gene calls are empty. Which probably means the gene\
                                caller reported no genes for your contigs.")

        self.set_next_available_id(t.hmm_hits_table_name)
        self.set_next_available_id(t.hmm_hits_splits_table_name)


    def populate_search_tables(self, sources={}):
        # if we end up generating a temporary file for amino acid sequences:
        if not len(sources):
            import anvio.data.hmm
            sources = anvio.data.hmm.sources

        if not sources:
            return

        target_files_dict = {}

        tmp_directory_path = filesnpaths.get_temp_directory_path()

        # here we will go through targets and populate target_files_dict based on what we find among them.
        targets = set([s['target'] for s in list(sources.values())])
        for target in targets:

            alphabet, context = utils.anvio_hmm_target_term_to_alphabet_and_context(target)

            self.run.info('Target found', '%s:%s' % (alphabet, context))

            class Args: pass
            args = Args()
            args.contigs_db = self.db_path
            contigs_db = ContigsSuperclass(args)

            if context == 'GENE':
                target_files_dict['%s:GENE' % alphabet] = os.path.join(tmp_directory_path, '%s_gene_sequences.fa' % alphabet)
                contigs_db.gen_FASTA_file_of_sequences_for_gene_caller_ids(output_file_path=target_files_dict['%s:GENE' % alphabet],
                                                                           simple_headers=True,
                                                                           rna_alphabet=True if alphabet=='RNA' else False,
                                                                           report_aa_sequences=True if alphabet=='AA' else False)
            elif context == 'CONTIG':
                if alphabet == 'AA':
                    raise ConfigError("You are somewhere you shouldn't be. You came here because you thought it would be OK\
                                       to ask for AA sequences in the CONTIG context. The answer to that is 'no, thanks'. If\
                                       you think this is dumb, please let us know.")
                else:
                    target_files_dict['%s:CONTIG' % alphabet] = os.path.join(tmp_directory_path, '%s_contig_sequences.fa' % alphabet)
                    utils.export_sequences_from_contigs_db(self.db_path,
                                                           target_files_dict['%s:CONTIG' % alphabet],
                                                           rna_alphabet=True if alphabet=='RNA' else False)

        commander = HMMer(target_files_dict, num_threads_to_use=self.num_threads_to_use)

        for source in sources:
            alphabet, context = utils.anvio_hmm_target_term_to_alphabet_and_context(sources[source]['target'])

            kind_of_search = sources[source]['kind']
            domain = sources[source]['domain']
            all_genes_searched_against = sources[source]['genes']
            hmm_model = sources[source]['model']
            reference = sources[source]['ref']
            noise_cutoff_terms = sources[source]['noise_cutoff_terms']

            hmm_scan_hits_txt = commander.run_hmmscan(source,
                                                      alphabet,
                                                      context,
                                                      kind_of_search,
                                                      domain,
                                                      len(all_genes_searched_against),
                                                      hmm_model,
                                                      reference,
                                                      noise_cutoff_terms)

            if not hmm_scan_hits_txt:
                search_results_dict = {}
            else:
                parser = parser_modules['search']['hmmscan'](hmm_scan_hits_txt, alphabet=alphabet, context=context)
                search_results_dict = parser.get_search_results()

            if not len(search_results_dict):
                run.info_single("The HMM source '%s' returned 0 hits. SAD (but it's stil OK)." % source, nl_before=1)


            if context == 'CONTIG':
                # we are in trouble here. because our search results dictionary contains no gene calls, but contig
                # names that contain our hits. on the other hand, the rest of the code outside of this if statement
                # expects a `search_results_dict` with gene callers id in it. so there are two things we need to do
                # to do. one is to come up with some new gene calls and add them to the contigs database. so things
                # will go smoothly downstream. two, we will need to update our `search_results_dict` so it looks
                # like a a dictionary the rest of the code expects with `gene_callers_id` fields. both of these
                # steps are going to be taken care of in the following function. magic.

                self.run.warning("Alright! You just called an HMM profile that runs on contigs. Because it is not\
                                 working with anvi'o gene calls directly, the resulting hits will need to be added\
                                 as 'new gene calls' into the contigs database. This is a new feature, and if it\
                                 starts screwing things up for you please let us know. Other than that you're pretty\
                                 much golden. Carry on.",
                                 header="Psst. Your fancy HMM profile '%s' speaking" % source,
                                 lc="green")

                num_hits_before = len(search_results_dict)
                search_results_dict = utils.get_pruned_HMM_hits_dict(search_results_dict)
                num_hits_after = len(search_results_dict)

                if num_hits_before != num_hits_after:
                    self.run.info('Pruned', '%d out of %d hits were removed due to redundancy' % (num_hits_before - num_hits_after, num_hits_before))

                search_results_dict = self.add_new_gene_calls_to_contigs_db_and_update_serach_results_dict(kind_of_search, search_results_dict)

            self.append(source, reference, kind_of_search, domain, all_genes_searched_against, search_results_dict)

        # FIXME: I have no clue why importing the anvio module is necessary at this point,
        #        but without this, mini test fails becasue "`anvio.DEBUG` is being used
        #        before initialization". nonsense.
        import anvio
        if not anvio.DEBUG:
            commander.clean_tmp_dirs()
            for v in list(target_files_dict.values()):
                os.remove(v)


    def add_new_gene_calls_to_contigs_db_and_update_serach_results_dict(self, source, search_results_dict):
        """Add new gene calls to the contigs database and update the HMM `search_results_dict`.

           When we are looking for HMM hits in the context of CONTIGS, our hits do not
           related to the gene calls we already have in a given contigs database. One
           slution is to add additional gene calls for a given set of HMM hits to keep
           them in the database."""

        # we will first learn the next available id in the gene callers table
        database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))
        next_id = database.get_max_value_in_column('genes_in_contigs', 'gene_callers_id') + 1
        database.disconnect()

        additional_gene_calls = {}
        for e in search_results_dict.values():
            start = e['start']
            stop = e['stop']

            if stop > start:
                direction = 'f'
            else:
                direction = 'r'
                stop, start = start, stop

            partial = 0 if ((stop - start) % 3 == 0) else 1

            # add a new gene call in to the dictionary
            additional_gene_calls[next_id] = {'contig': e['contig_name'],
                                              'start': start,
                                              'stop': stop,
                                              'direction': direction,
                                              'partial': partial,
                                              'source': source,
                                              'version': 'unknown'
                                            }

            # update the search results dictionary with gene callers id:
            e['gene_callers_id'] = next_id

            # update the next available gene callers id:
            next_id += 1

        if not len(additional_gene_calls):
            return search_results_dict

        # update the contigs db with the gene calls in `additional_gene_calls` dict.
        gene_calls_table = TablesForGeneCalls(self.db_path, run=terminal.Run(verbose=False))
        gene_calls_table.use_external_gene_calls_to_populate_genes_in_contigs_table(input_file_path=None,
                                                                                    gene_calls_dict=additional_gene_calls,
                                                                                    ignore_internal_stop_codons=True)
        gene_calls_table.populate_genes_in_splits_tables(gene_calls_dict=additional_gene_calls)

        # refresh the gene calls dict
        self.init_gene_calls_dict()

        self.run.info('Gene calls added to db', '%d (from source "%s")' % (len(additional_gene_calls), source))

        return search_results_dict


    def remove_source(self, source):
        self.delete_entries_for_key('source', source, [t.hmm_hits_info_table_name, t.hmm_hits_table_name, t.hmm_hits_splits_table_name])


    def append(self, source, reference, kind_of_search, domain, all_genes, search_results_dict):
        # we want to define unique identifiers for each gene first. this information will be used to track genes that will
        # break into multiple pieces due to arbitrary split boundaries. while doing that, we will add the 'source' info
        # into the dictionary, so it perfectly matches to the table structure

        for entry_id in search_results_dict:
            hit = search_results_dict[entry_id]

            gene_call = self.gene_calls_dict[hit['gene_callers_id']]

            hit['gene_unique_identifier'] = hashlib.sha224('_'.join([gene_call['contig'], hit['gene_name'], str(gene_call['start']), str(gene_call['stop'])]).encode('utf-8')).hexdigest()
            hit['source'] = source

        self.remove_source(source)

        database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))

        # push information about this search result into serach_info table.
        db_entries = [source, reference, kind_of_search, domain, ', '.join(all_genes)]
        database._exec('''INSERT INTO %s VALUES (?,?,?,?,?)''' % t.hmm_hits_info_table_name, db_entries)

        # if our search results were empty, we can return from here.
        if not len(search_results_dict):
            database.disconnect()
            return

        # then populate serach_data table for each contig.
        db_entries = []
        for hit in list(search_results_dict.values()):
            entry_id = self.next_id(t.hmm_hits_table_name)
            db_entries.append(tuple([entry_id] + [hit[h] for h in t.hmm_hits_table_structure[1:]]))
            # tiny hack here: for each hit, we are generating a unique id (`entry_id`), and feeding that information
            #                 back into the dictionary to pass it to processing of splits, so each split-level
            #                 entry knows who is their parent.
            hit['hmm_hit_entry_id'] = entry_id

        database._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?,?,?)''' % t.hmm_hits_table_name, db_entries)

        db_entries = self.process_splits(search_results_dict)
        database._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?)''' % t.hmm_hits_splits_table_name, db_entries)

        database.disconnect()


    def process_splits(self, search_results_dict):
        hits_per_contig = {}
        for hit in list(search_results_dict.values()):
            contig_name = self.gene_calls_dict[hit['gene_callers_id']]['contig']

            if contig_name in hits_per_contig:
                hits_per_contig[contig_name].append(hit)
            else:
                hits_per_contig[contig_name] = [hit]

        db_entries_for_splits = []

        for contig in self.contigs_info:
            if contig not in hits_per_contig:
                # no hits for this contig. pity!
                continue

            for split_name in self.contig_name_to_splits[contig]:
                split_start = self.splits_info[split_name]['start']
                split_stop = self.splits_info[split_name]['end']

                # FIXME: this really needs some explanation.
                for hit in hits_per_contig[contig]:
                    hit_start = self.gene_calls_dict[hit['gene_callers_id']]['start']
                    hit_stop = self.gene_calls_dict[hit['gene_callers_id']]['stop']

                    if hit_stop > split_start and hit_start < split_stop:
                        gene_length = hit_stop - hit_start
                        # if only a part of the gene is in the split:
                        start_in_split = (split_start if hit_start < split_start else hit_start) - split_start
                        stop_in_split = (split_stop if hit_stop > split_stop else hit_stop) - split_start
                        percentage_in_split = (stop_in_split - start_in_split) * 100.0 / gene_length

                        db_entry = tuple([self.next_id(t.hmm_hits_splits_table_name), hit['hmm_hit_entry_id'], split_name, percentage_in_split, hit['source']])
                        db_entries_for_splits.append(db_entry)

        return db_entries_for_splits


