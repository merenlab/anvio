# -*- coding: utf-8
# pylint: disable=line-too-long

import os
import hashlib
import gzip
import shutil

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.utils as utils
import anvio.hmmops as hmmops
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.filesnpaths as filesnpaths

from anvio.drivers.hmmer import HMMer
from anvio.tables.tableops import Table
from anvio.parsers import parser_modules
from anvio.dbops import ContigsSuperclass
from anvio.errors import ConfigError, StupidHMMError
from anvio.tables.genecalls import TablesForGeneCalls
from anvio.tables.genefunctions import TableForGeneFunctions


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
    def __init__(self, db_path, num_threads_to_use=1, run=run, progress=progress, initializing_for_deletion=False, just_do_it=False,
                 hmm_program_to_use='hmmscan', hmmer_output_directory=None, get_domain_table_output=False, add_to_functions_table=False):
        self.num_threads_to_use = num_threads_to_use
        self.db_path = db_path
        self.just_do_it = just_do_it
        self.hmm_program = hmm_program_to_use or 'hmmscan'
        self.hmmer_output_dir = hmmer_output_directory
        self.hmmer_desired_output = ('table', 'domtable') if get_domain_table_output else 'table'
        self.add_to_functions_table = add_to_functions_table

        utils.is_contigs_db(self.db_path)
        filesnpaths.is_program_exists(self.hmm_program)

        self.contigs_db_hash = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path)).get_meta_value('contigs_db_hash')

        Table.__init__(self, self.db_path, anvio.__contigs__version__, run, progress)

        self.init_gene_calls_dict()

        if not len(self.gene_calls_dict):
            if self.genes_are_called:
                self.run.warning("Tables in this contigs database that should contain gene calls are empty despite the fact that "
                                 "you didn't skip the gene calling step while generating this contigs database. This probably means "
                                 "that the gene caller did not find any genes among contigs. This is OK for now. But might explode "
                                 "later. If it does explode and you decide to let us know about that problem, please remember to mention "
                                 "this warning. By the way, this warning probably has been seen by like only 2 people on the planet. Who "
                                 "works with contigs with no gene calls? A better implementation of anvi'o will unite researchers who "
                                 "study weird stuff.")
            else:
                self.run.warning("It seems you have skipped gene calling step while generating your contigs database, and you have no "
                                 "genes calls in tables that should contain gene calls. Anvi'o will let you go with this since some HMM "
                                 "sources only operate on DNA sequences, and at this point it doesn't know which HMMs you wish to run. "
                                 "If the lack of genes causes a problem, you will get another error message later probably :/")

        if not initializing_for_deletion:
            self.set_next_available_id(t.hmm_hits_table_name)


    def check_sources(self, sources):

        if self.add_to_functions_table: # check that source is not already in gene_functions table
            gene_function_sources_in_db = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path)).get_meta_value('gene_function_sources')
            sources_in_db = set(gene_function_sources_in_db.split(',') if gene_function_sources_in_db else [])
            sources_need_to_be_removed = set(sources.keys()).intersection(sources_in_db)

            if len(sources_need_to_be_removed):
                source_string = ', '.join(sources_need_to_be_removed)
                raise ConfigError("Some of the HMM sources are already in the gene functions table in the database and anvi'o "
                                  "doesn't want to overwrite them. If YOU want to overwrite them, however, (because you do you, "
                                  "friend) you can do that by "
                                  "running `anvi-delete-functions` first, and then re-running this program. Here are the sources "
                                  f"that you would need to delete: {source_string}")
        else: # default checks for hmm_hits table
            sources_in_db = list(hmmops.SequencesForHMMHits(self.db_path).hmm_hits_info.keys())

            if 'Ribosomal_RNAs' in sources_in_db and len([s for s in sources if s.startswith('Ribosomal_RNA_')]):
                raise ConfigError(f"Here is one more additional step we need to you take care of before we can go forward: Your contigs database "
                                  f"already contains HMMs from an older `Ribosomal_RNAs` model anvi'o no longer uses AND you are about to run "
                                  f"its newer models that do the same thing (but better). Since Ribosomal RNA models add new gene calls to the "
                                  f"database, running newer models without first cleaning up the old ones will result in duplication of gene calls "
                                  f"as examplified here: https://github.com/merenlab/anvio/issues/1598. Anvi'o could've removed the `Ribosomal_RNAs` "
                                  f"model for you automatically, but the wisdom tells us that the person who passes the sentence should swing the "
                                  f"sword. Here it is for your grace: \"anvi-delete-hmms -c {self.db_path} --hmm-source Ribosomal_RNAs\".")

            sources_need_to_be_removed = set(sources.keys()).intersection(sources_in_db)

            if len(sources_need_to_be_removed):
                if self.just_do_it:
                    for source_name in sources_need_to_be_removed:
                        self.remove_source(source_name)
                else:
                    raise ConfigError("Some of the HMM sources you wish to run on this database are already in the database and anvi'o "
                                      "refuses to overwrite them without your explicit input. You can either use `anvi-delete-hmms` "
                                      "to remove them first, or run this program with `--just-do-it` flag so anvi'o would remove all "
                                      "for you. Here are the list of HMM sources that need to be removed: '%s'." % (', '.join(sources_need_to_be_removed)))


    def hmmpress_sources(self, sources, tmp_dir):
        """This function runs hmmpress on the hmm profiles.

        It returns the locations of each hmmpressed file path in a dictionary keyed by the source.
        """
        hmmpressed_file_paths = {}
        for source in sources:
            model_file = sources[source]['model']
            hmm_file_path = os.path.join(tmp_dir, source + '.hmm')
            hmm_file = open(hmm_file_path, 'wb')
            hmm_file.write(gzip.open(model_file, 'rb').read())
            hmm_file.close()

            log_file_path = log_file_path = os.path.join(tmp_dir, 'hmmpress.log')
            cmd_line = ['hmmpress', hmm_file_path]
            ret_val = utils.run_command(cmd_line, log_file_path)

            hmmpressed_file_paths[source] = hmm_file_path

            if ret_val:
                raise ConfigError("Sadly, anvi'o failed while attempting to compress the HMM model for source %s. You can check out the log file (%s) for "
                                  "more detailed information on why this happened." % (source, log_file_path))
        return hmmpressed_file_paths


    def populate_search_tables(self, sources={}):
        # make sure the output file is OK to write.
        filesnpaths.is_output_file_writable(self.db_path, ok_if_exists=True)

        # if we end up generating a temporary file for amino acid sequences:
        if not len(sources):
            import anvio.data.hmm
            sources = anvio.data.hmm.sources

        if not sources:
            return

        self.check_sources(sources)

        target_files_dict = {}

        tmp_directory_path = filesnpaths.get_temp_directory_path()

        hmmpressed_files = self.hmmpress_sources(sources, tmp_directory_path)

        self.run.info("Contigs DB", self.db_path)
        self.run.info("HMM sources", ', '.join(sources.keys()))

        # here we will go through targets and populate target_files_dict based on what we find among them.
        targets = set([s['target'] for s in list(sources.values())])
        have_hmm_sources_with_non_RNA_contig_context = False
        for target in targets:
            alphabet, context = utils.anvio_hmm_target_term_to_alphabet_and_context(target)

            if not self.genes_are_called and context != "CONTIG":
                raise ConfigError(f"You are in trouble. The gene calling was skipped for this contigs database, yet anvi'o asked to run an "
                                  f"HMM profile that wishes to operate on {context} context using the {alphabet} alphabet. It is not OK. You "
                                  f"still could run HMM profiles that does not require gene calls to be present (such as the HMM profile that "
                                  f"identifies Ribosomal RNAs in contigs, but for that you would have to explicitly ask for it by using the "
                                  f"additional parameter '--installed-hmm-profile PROFILE_NAME_HERE').")

            self.run.info('Alphabet/context target found', f"{alphabet}:{context}")

            if context == 'CONTIG' and alphabet != 'RNA':
                have_hmm_sources_with_non_RNA_contig_context =True

            class Args: pass
            args = Args()
            args.contigs_db = self.db_path
            contigs_db = ContigsSuperclass(args, r=terminal.Run(verbose=False))

            if context == 'GENE':
                target_file_path = os.path.join(tmp_directory_path, f'{alphabet}_gene_sequences.fa')

                contigs_db.get_sequences_for_gene_callers_ids(output_file_path=target_file_path,
                                                              simple_headers=True,
                                                              rna_alphabet=True if alphabet=='RNA' else False,
                                                              report_aa_sequences=True if alphabet=='AA' else False)

                target_files_dict[f'{alphabet}:GENE'] = target_file_path
            elif context == 'CONTIG':
                if alphabet == 'AA':
                    raise ConfigError("You are somewhere you shouldn't be. You came here because you thought it would be OK "
                                      "to ask for AA sequences in the CONTIG context. The answer to that is 'no, thanks'. If "
                                      "you think this is dumb, please let us know.")
                else:
                    target_file_path = os.path.join(tmp_directory_path, f'{alphabet}_contig_sequences.fa')

                    utils.export_sequences_from_contigs_db(self.db_path,
                                                           target_file_path,
                                                           rna_alphabet=True if alphabet=='RNA' else False)

                    target_files_dict[f'{alphabet}:CONTIG'] = target_file_path

        # now we know our sequences
        self.run.info('Target sequences determined',
                      '; '.join([f"{pp(utils.get_num_sequences_in_fasta(file_path))} sequences for {target}" \
                                    for target, file_path in target_files_dict.items()]))

        if have_hmm_sources_with_non_RNA_contig_context:
            # in that case, we should remind people what's up.
            self.run.warning("The HMM profiles that are about to be run includes at least one HMM profile that runs on "
                             "contigs and not genes. Thus, this HMM operation will not be working with gene calls anvi'o "
                             "already knows about. Which means, the resulting hits will need to be added as 'new gene calls' "
                             "into the contigs database. So far so good. But because we are in the realm of contigs rather "
                             "than genes, the resulting HMM hits will unlikely correspond to open reading frames that are "
                             "supposed to be translated (such as ribosomal RNAs). While anvi'o adds new gene calls to your "
                             "contigs database for these hits, it will NOT report amino acid sequences for the "
                             "new gene calls that will emerge from these HMMs, expecting you to judge whether this will "
                             "influence your pangenomic analyses or other things you thought you would be doing with the "
                             "result of this HMM search downstream. If you do not feel like being the judge of anything today "
                             "you can move on yet remember to remember this if things look somewhat weird later on.",
                             header="THE MORE YOU KNOW 🌈", lc="green")

        commander = HMMer(target_files_dict, num_threads_to_use=self.num_threads_to_use, program_to_use=self.hmm_program)

        for source in sources:
            alphabet, context = utils.anvio_hmm_target_term_to_alphabet_and_context(sources[source]['target'])

            if alphabet in ['DNA', 'RNA'] and 'domtable' in self.hmmer_desired_output:
                raise ConfigError(f"Domain table output was requested (probably with the --get-domtable-output flag, "
                                  f"does that look familiar?) but unfortunately this option is incompatible with the "
                                  f"current source of HMM profiles, {source}, because this source uses a nucleotide "
                                  f"alphabet.")

            kind_of_search = sources[source]['kind']
            domain = sources[source]['domain']
            all_genes_searched_against = sources[source]['genes']
            hmm_model = hmmpressed_files[source]
            reference = sources[source]['ref']
            noise_cutoff_terms = sources[source]['noise_cutoff_terms']

            hmmer_output = commander.run_hmmer(source,
                                              alphabet,
                                              context,
                                              kind_of_search,
                                              domain,
                                              len(all_genes_searched_against),
                                              hmm_model,
                                              reference,
                                              noise_cutoff_terms,
                                              desired_output=self.hmmer_desired_output,
                                              hmmer_output_dir=self.hmmer_output_dir)

            if self.hmmer_output_dir:
                self.run.info("HMMER output directory", self.hmmer_output_dir)

            if not isinstance(hmmer_output, tuple):
                hmm_scan_hits_txt = hmmer_output
            else:
                hmm_scan_hits_txt,domain_hits_txt = hmmer_output
                self.run.info("Domain table output", domain_hits_txt)

            if not hmm_scan_hits_txt:
                search_results_dict = {}
            else:
                try:
                    parser = parser_modules['search']['hmmer_table_output'](hmm_scan_hits_txt, alphabet=alphabet, context=context, program=self.hmm_program)
                except StupidHMMError as e:
                    raise ConfigError(f"Unfortunately something went wrong while anvi'o was trying to parse some HMM output for your data. "
                                      f"This error is typically due to contig names that are long and variable in length, which that "
                                      f"confuses HMMER and so it generates output tables that are simply unparseable. Anvi'o does its best, "
                                      f"but occasionally fails, which leads to this error. If you are curious why is this happening, you can take a "
                                      f"look at this issue where this issue is described: https://github.com/merenlab/anvio/issues/1564. "
                                      f"Solution to this is relatively easy: use `anvi-script-reformat-fasta` with `--simplify-names` flag "
                                      f"BEFORE generating your contigs database as we advice you to. Sorry you came all this way just to "
                                      f"find out about this :/ Here is the origial error message anvi'o produced from the code beneath: {e}.")

                search_results_dict = parser.get_search_results()

            if not len(search_results_dict):
                run.info_single("The HMM source '%s' returned 0 hits. SAD (but it's OK)." % source, nl_before=1)

            if context == 'CONTIG':
                # we are in trouble here. because our search results dictionary contains no gene calls, but contig
                # names contain our hits. on the other hand, the rest of the code outside of this if statement
                # expects a `search_results_dict` with gene caller ids in it. so there are two things we need to do.
                # one is to come up with some new gene calls and add them to the contigs database. so things
                # will go smoothly downstream. two, we will need to update our `search_results_dict` so it looks
                # like a a dictionary the rest of the code expects with `gene_callers_id` fields. both of these
                # steps are going to be taken care of in the following function. magic.
                num_hits_before = len(search_results_dict)
                search_results_dict = utils.get_pruned_HMM_hits_dict(search_results_dict)
                num_hits_after = len(search_results_dict)

                if num_hits_before != num_hits_after:
                    self.run.info('Pruned', '%d out of %d hits were removed due to redundancy' % (num_hits_before - num_hits_after, num_hits_before))

                search_results_dict = self.add_new_gene_calls_to_contigs_db_and_update_serach_results_dict(kind_of_search,
                                                                                                           search_results_dict,
                                                                                                           skip_amino_acid_sequences=True)

            if self.add_to_functions_table: # add to gene_functions table (upon request)
                self.append_to_gene_functions_table(source, search_results_dict)
            else:                           # add to hmm_hits table (default)
                self.append_to_hmm_hits_table(source, reference, kind_of_search, domain, all_genes_searched_against, search_results_dict)


        # FIXME: I have no clue why importing the anvio module is necessary at this point,
        #        but without this, mini test fails becasue "`anvio.DEBUG` is being used
        #        before initialization". nonsense.
        import anvio
        if not anvio.DEBUG:
            commander.clean_tmp_dirs()
            for v in list(target_files_dict.values()):
                os.remove(v)

            shutil.rmtree(tmp_directory_path)


    def add_new_gene_calls_to_contigs_db_and_update_serach_results_dict(self, source, search_results_dict, skip_amino_acid_sequences=False):
        """Add new gene calls to the contigs database and update the HMM `search_results_dict`.

           When we are looking for HMM hits in the context of CONTIGS, our hits do not
           related to the gene calls we already have in a given contigs database. One
           slution is to add additional gene calls for a given set of HMM hits to keep
           them in the database."""

        if not len(search_results_dict):
            return search_results_dict

        # we will first learn the next available id in the gene callers table
        database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))
        next_id = database.get_max_value_in_column('genes_in_contigs', 'gene_callers_id', value_if_empty=0) + 1
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
                                              'call_type': constants.gene_call_types['NONCODING'] if skip_amino_acid_sequences else constants.gene_call_types['CODING'],
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
                                                                                    ignore_internal_stop_codons=True,
                                                                                    skip_amino_acid_sequences=skip_amino_acid_sequences)
        gene_calls_table.populate_genes_in_splits_tables(gene_calls_dict=additional_gene_calls)

        # refresh the gene calls dict
        self.init_gene_calls_dict()

        self.run.info('Gene calls added to db', '%d (from source "%s")' % (len(additional_gene_calls), source))

        return search_results_dict


    def remove_source(self, source):
        """Remove an HMM source from the database."""

        tables_with_source = [
            t.hmm_hits_info_table_name,
            t.hmm_hits_table_name,
            t.hmm_hits_splits_table_name,
            t.genes_in_contigs_table_name,
            t.gene_function_calls_table_name,
        ]

        tables_with_gene_callers_id = [
            t.gene_amino_acid_sequences_table_name,
            t.genes_taxonomy_table_name,
            t.genes_in_splits_table_name
        ]

        # delete entries from tables with 'source' column
        self.delete_entries_for_key('source', source, tables_with_source)

        # collect gene caller ids that were added to the db via the HMM source
        gene_caller_ids_to_remove = set(key for key, val in self.gene_calls_dict.items() if val['source'] == source)

        # if there are any, remove them from tables with 'gene_callers_id' column
        if len(gene_caller_ids_to_remove):
            database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))

            CLAUSE = "gene_callers_id in (%s)" % (','.join([str(x) for x in gene_caller_ids_to_remove]))
            for table in tables_with_gene_callers_id:
                database.remove_some_rows_from_table(table, CLAUSE)

            database.disconnect()

            run.warning("%d gene caller ids that were added via the HMM source have been removed from \"%s\"" \
                        % (len(gene_caller_ids_to_remove), ', '.join(tables_with_gene_callers_id)))


    def append_to_gene_functions_table(self, source, search_results_dict):
        """Append custom HMM hits to the gene functions table in the contigs database."""

        # get an instance of gene functions table
        gene_function_calls_table = TableForGeneFunctions(self.db_path, self.run, self.progress)

        # first we massage the hmm_hits dictionary to match expected input for the gene_functions table
        if search_results_dict:
            for entry_id in search_results_dict:
                hit = search_results_dict[entry_id]
                hit['source'] = source
                hit['accession'] = hit['gene_hmm_id']
                hit['function'] = hit['gene_name']

            gene_function_calls_table.create(search_results_dict)
        else:
            self.run.warning("There are no hits to add to the database. Returning empty handed, "
                             f"but still adding {source} as a functional source.")
            gene_function_calls_table.add_empty_sources_to_functional_sources({source})


    def append_to_hmm_hits_table(self, source, reference, kind_of_search, domain, all_genes, search_results_dict):
        """Append a new HMM source in the contigs database."""

        # just to make 100% sure.
        if source in list(hmmops.SequencesForHMMHits(self.db_path).hmm_hits_info.keys()):
            raise ConfigError("The source '%s' you're trying to append is already in the database :( "
                              "You should have never been able to come here in the code unless you "
                              "have passed the `check_sources` sanity check. Very good but not "
                              "good really. Bad. Bad you." % source)

        # we want to define unique identifiers for each gene first. this information will be used to track genes that will
        # break into multiple pieces due to arbitrary split boundaries. while doing that, we will add the 'source' info
        # into the dictionary, so it perfectly matches to the table structure
        for entry_id in search_results_dict:
            hit = search_results_dict[entry_id]

            gene_call = self.gene_calls_dict[hit['gene_callers_id']]

            hit['gene_unique_identifier'] = hashlib.sha224('_'.join([str(self.contigs_db_hash),
                                                                     gene_call['contig'],
                                                                     hit['gene_name'],
                                                                     str(gene_call['start']),
                                                                     str(gene_call['stop'])]).encode('utf-8')).hexdigest()
            hit['source'] = source

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
        database._exec_many('''INSERT INTO %s VALUES (?,?,?,?)''' % t.hmm_hits_splits_table_name, db_entries)

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

                        db_entry = tuple([hit['hmm_hit_entry_id'], split_name, percentage_in_split, hit['source']])
                        db_entries_for_splits.append(db_entry)

        return db_entries_for_splits
