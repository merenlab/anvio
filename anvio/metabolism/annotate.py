import os
import pickle
import shutil

import multiprocess as multiprocessing

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.drivers.hmmer import HMMer
from anvio.parsers import parser_modules
from anvio.dbops import ContigsDatabase, ContigsSuperclass
from anvio.tables.genefunctions import TableForGeneFunctions

from anvio.metabolism.context import KeggContext
from anvio.metabolism.modulesdb import ModulesDatabase



class AnnotationWorker(multiprocessing.Process):
    """Worker process that annotates one FASTA partition end-to-end.

    Runs HMM search, bitscore filtering, optional bitscore heuristic, and module/BRITE
    lookups entirely within the worker process. Read-only shared data (KO dicts and
    module lookups) is loaded from a shared memory segment at worker startup rather than
    being pickled into each process separately. Results are returned via a queue as
    compact annotation dicts — raw hits are never sent back to the main process.

    Parameters
    ==========
    Most of the arguments provided to the `anvi-run-kegg-kofams` CLI plus
    partition_fasta : str
        Path to the FASTA file containing a subset of gene sequences
    partition_index : int
        The numerical ID of the worker
    log_file_path : str
        Path to the log file where this worker will report progress. If log_bitscores
        is True, the output file for the bitscore log is derived from this parameter
    result_queue : multiprocess.Manager.Queue
        The queue where this worker will send its results (or exceptions)
    shared_data_shm_name : str
        The name of the multiprocessing.shared_memory.SharedMemory object where all
        workers can access KEGG information, such as the dictionaries of KOfam thresholds
    """

    def __init__(self, partition_fasta, partition_index, hmm_program, kofam_hmm_path,
                 stray_ko_hmm_path, include_stray_kos, keep_all_hits, skip_bitscore_heuristic,
                 no_hmmer_prefiltering, log_bitscores, skip_brite_hierarchies,
                 bitscore_heuristic_e_value, bitscore_heuristic_bitscore_fraction,
                 log_file_path, result_queue, shared_data_shm_name, shared_data_size):
        super().__init__()
        self.partition_fasta = partition_fasta
        self.partition_index = partition_index
        self.hmm_program = hmm_program
        self.kofam_hmm_path = kofam_hmm_path
        self.stray_ko_hmm_path = stray_ko_hmm_path
        self.include_stray_kos = include_stray_kos
        self.keep_all_hits = keep_all_hits
        self.skip_bitscore_heuristic = skip_bitscore_heuristic
        self.no_hmmer_prefiltering = no_hmmer_prefiltering
        self.log_bitscores = log_bitscores
        self.skip_brite_hierarchies = skip_brite_hierarchies
        self.bitscore_heuristic_e_value = bitscore_heuristic_e_value
        self.bitscore_heuristic_bitscore_fraction = bitscore_heuristic_bitscore_fraction
        self.log_file_path = log_file_path
        self.result_queue = result_queue
        self.shared_data_shm_name = shared_data_shm_name
        self.shared_data_size = shared_data_size


    def _load_shared_data(self):
        """Attach to the shared memory segment and unpickle the KO and module dicts.

        Workers close (but do not unlink) the segment — only the main process unlinks.
        """
        from multiprocessing.shared_memory import SharedMemory
        from multiprocessing import resource_tracker

        shm = SharedMemory(name=self.shared_data_shm_name, create=False)
        resource_tracker.unregister('/' + shm.name, 'shared_memory')
        bundle = pickle.loads(bytes(shm.buf[:self.shared_data_size]))
        shm.close()

        self.ko_dict = bundle['ko_dict']
        self.stray_ko_dict = bundle['stray_ko_dict']
        self.ko_to_modules = bundle['ko_to_modules']
        self.ko_to_module_names = bundle['ko_to_module_names']
        self.ko_to_module_classes = bundle['ko_to_module_classes']
        self.ko_to_brite = bundle['ko_to_brite']
        self.total_seqs = bundle['total_seqs']


    def _get_kofam_annotation(self, knum):
        """Returns the function definition string for a given KO accession."""

        if knum in self.ko_dict:
            return self.ko_dict[knum]['definition']
        if self.stray_ko_dict and knum in self.stray_ko_dict:
            return self.stray_ko_dict[knum]['definition']
        return f"Unknown function with KO num {knum}"


    def _add_module_entries(self, knum, gcid, counter, kegg_module_names_dict, kegg_module_classes_dict):
        """Append KEGG_Module and KEGG_Class annotation entries for a KO to the running dicts.

        If the KO belongs to no modules, nothing is added. Both the module name and class
        strings use `!!!` as the inter-module separator, consistent with how other
        multi-value annotations are stored in anvi'o.

        Parameters
        ==========
        knum : str
            KO accession (e.g. 'K00001').
        gcid : int
            gene_callers_id for the annotated gene.
        counter : int
            Integer key under which to store the new entries.
        kegg_module_names_dict : dict
            Running KEGG_Module annotation dict; updated in-place.
        kegg_module_classes_dict : dict
            Running KEGG_Class annotation dict; updated in-place.
        """

        mods = self.ko_to_modules.get(knum)
        if not mods:
            return
        mod_annotation = "!!!".join(mods)
        names = self.ko_to_module_names.get(knum, {})
        classes = self.ko_to_module_classes.get(knum, [])
        mod_name_annotation = "!!!".join(names.get(m, '') for m in mods)
        mod_class_annotation = "!!!".join(classes)
        kegg_module_names_dict[counter] = {
            'gene_callers_id': gcid,
            'source': 'KEGG_Module',
            'accession': mod_annotation,
            'function': mod_name_annotation,
            'e_value': None,
        }
        kegg_module_classes_dict[counter] = {
            'gene_callers_id': gcid,
            'source': 'KEGG_Class',
            'accession': mod_annotation,
            'function': mod_class_annotation,
            'e_value': None,
        }


    def _build_brite_entry(self, knum, gcid):
        """Build a KEGG_BRITE annotation dict entry for a KO, or return None if no BRITE data exists for it.

        Parameters
        ==========
        knum : str
            KO accession (e.g. 'K00001').
        gcid : int
            gene_callers_id for the annotated gene.

        Returns
        =======
        entry : dict or None
            A gene-functions-table-style dict with 'source' set to 'KEGG_BRITE', or None
            if the KO has no entries in the BRITE table.
        """

        brite_rows = self.ko_to_brite.get(knum)
        if not brite_rows:
            return None
        hierarchy_accession = ""
        categorizations = ""
        for row in brite_rows[:-1]:
            hierarchy_accession += f"{row['hierarchy_accession']}!!!"
            categorizations += f"{row['hierarchy_name']}>>>{row['categorization']}!!!"
        last_row = brite_rows[-1]
        hierarchy_accession += last_row['hierarchy_accession']
        categorizations += f"{last_row['hierarchy_name']}>>>{last_row['categorization']}"
        return {
            'gene_callers_id': gcid,
            'source': 'KEGG_BRITE',
            'accession': hierarchy_accession,
            'function': categorizations,
            'e_value': None,
        }


    def _parse_hits(self, hits_dict, hits_label, next_key,
                    functions_dict, kegg_module_names_dict, kegg_module_classes_dict,
                    kegg_brite_categorizations_dict, gcids_to_hits_dict, gcids_to_functions_dict):
        """Apply bitscore threshold filtering to a set of HMM hits and populate annotation dicts.

        Each hit whose bitscore meets the KOfam-defined threshold is stored in
        functions_dict along with any matching KEGG_Module, KEGG_Class, and KEGG_BRITE
        entries. Hits that fail the threshold are counted but not stored (unless
        self.keep_all_hits is True). All hits and gene-to-hit mappings are recorded in
        the tracking dicts so the bitscore heuristic can inspect them afterward.

        Parameters
        ==========
        hits_dict : dict
            Parsed HMM output keyed by hit index; values are dicts with at least
            'gene_name' (KO accession), 'gene_callers_id', 'bit_score',
            'domain_bit_score', and 'e_value'.
        hits_label : str
            Label for this hit set (e.g. 'KOfam' or 'No-threshold KO').
        next_key : int
            Starting integer key for new entries in the annotation dicts.
        functions_dict : dict
            Running KOfam annotation dict; updated in-place.
        kegg_module_names_dict : dict
            Running KEGG_Module annotation dict; updated in-place.
        kegg_module_classes_dict : dict
            Running KEGG_Class annotation dict; updated in-place.
        kegg_brite_categorizations_dict : dict
            Running KEGG_BRITE annotation dict; updated in-place.
        gcids_to_hits_dict : dict
            Maps gene_callers_id → {label → [hit_key, ...]}; updated in-place for
            later use by the bitscore heuristic.
        gcids_to_functions_dict : dict
            Maps gene_callers_id → [annotation key, ...] for genes that passed
            filtering; updated in-place.

        Returns
        =======
        counter : int
            The next available integer key after all entries added in this call.
        num_removed : int
            Number of hits that failed the bitscore threshold.
        """

        counter = next_key
        num_removed = 0

        for hit_key, hmm_hit in hits_dict.items():
            knum = hmm_hit['gene_name']
            gcid = hmm_hit['gene_callers_id']
            keep = False

            if gcid not in gcids_to_hits_dict:
                gcids_to_hits_dict[gcid] = {hits_label: [hit_key]}
            else:
                if hits_label not in gcids_to_hits_dict[gcid]:
                    gcids_to_hits_dict[gcid][hits_label] = [hit_key]
                else:
                    gcids_to_hits_dict[gcid][hits_label].append(hit_key)

            if knum in self.ko_dict:
                score_type = self.ko_dict[knum]['score_type']
                threshold = float(self.ko_dict[knum]['threshold'])
                if score_type == 'domain':
                    keep = hmm_hit['domain_bit_score'] >= threshold
                elif score_type == 'full':
                    keep = hmm_hit['bit_score'] >= threshold
                else:
                    raise ConfigError(f"The KO noise cutoff dictionary for {knum} has a strange score type "
                                      f"unknown to anvi'o: {score_type}")
            elif self.stray_ko_dict and knum in self.stray_ko_dict:
                score_type = self.stray_ko_dict[knum]['score_type']
                threshold = float(self.stray_ko_dict[knum]['threshold'])
                if score_type == 'domain':
                    keep = hmm_hit['domain_bit_score'] >= threshold
                elif score_type == 'full':
                    keep = hmm_hit['bit_score'] >= threshold
                else:
                    raise ConfigError(f"The KO noise cutoff dictionary for nt-KO {knum} has a strange score type "
                                      f"unknown to anvi'o: {score_type}")
            else:
                raise ConfigError(f"KO {knum} from {hits_label} HMM hits is not in the noise cutoff dictionary. "
                                  f"Anvi'o does not know how to distinguish strong hits from weak ones for this KO. "
                                  f"Please contact a developer.")

            if keep or self.keep_all_hits:
                functions_dict[counter] = {
                    'gene_callers_id': gcid,
                    'source': 'KOfam',
                    'accession': knum,
                    'function': self._get_kofam_annotation(knum),
                    'e_value': hmm_hit['e_value'],
                }
                if gcid not in gcids_to_functions_dict:
                    gcids_to_functions_dict[gcid] = [counter]
                else:
                    gcids_to_functions_dict[gcid].append(counter)
                self._add_module_entries(knum, gcid, counter, kegg_module_names_dict, kegg_module_classes_dict)
                if not self.skip_brite_hierarchies:
                    brite_entry = self._build_brite_entry(knum, gcid)
                    if brite_entry:
                        kegg_brite_categorizations_dict[counter] = brite_entry
                counter += 1
            else:
                num_removed += 1

        return counter, num_removed


    def _apply_heuristic(self, gcids_list, super_hits_dict, next_key,
                         functions_dict, kegg_module_names_dict, kegg_module_classes_dict,
                         kegg_brite_categorizations_dict, gcids_to_hits_dict, gcids_to_functions_dict):
        """Annotate unannotated genes using a relaxed bitscore heuristic.

        For each gene that received no annotation from threshold filtering, considers
        all of its HMM hits. If every hit that passes the relaxed criteria (e-value <=
        self.bitscore_heuristic_e_value and bitscore > self.bitscore_heuristic_bitscore_fraction
        * KO threshold) points to the same KO, that gene is annotated with that KO.
        This rescues likely-valid annotations that narrowly missed the strict threshold.

        Parameters
        ==========
        gcids_list : list
            All gene_callers_ids in this FASTA partition.
        super_hits_dict : dict
            Maps label → raw hits dict as returned by the HMM parser.
        next_key : int
            Starting integer key for new entries in the annotation dicts.
        functions_dict : dict
            Running KOfam annotation dict; updated in-place.
        kegg_module_names_dict : dict
            Running KEGG_Module annotation dict; updated in-place.
        kegg_module_classes_dict : dict
            Running KEGG_Class annotation dict; updated in-place.
        kegg_brite_categorizations_dict : dict
            Running KEGG_BRITE annotation dict; updated in-place.
        gcids_to_hits_dict : dict
            Maps gene_callers_id → {label → [hit_key, ...]} from _parse_hits.
        gcids_to_functions_dict : dict
            Maps gene_callers_id → [annotation key, ...]; genes already present
            here were annotated by threshold filtering and are skipped.

        Returns
        =======
        next_key : int
            The next available integer key after all entries added in this call.
        num_added : int
            Total number of annotations added by the heuristic.
        num_stray_added : int
            Subset of added annotations that came from nt-KO (no-threshold) hits.
        """

        num_added = 0
        num_stray_added = 0

        for gcid in gcids_list:
            if gcid in gcids_to_functions_dict:
                continue
            if gcid not in gcids_to_hits_dict:
                continue

            decent_hit_kos = set()
            best_e_value = 100
            best_hit_key = None
            best_hit_label = None

            for hit_label in gcids_to_hits_dict[gcid]:
                for hit_key in gcids_to_hits_dict[gcid][hit_label]:
                    knum = super_hits_dict[hit_label][hit_key]['gene_name']

                    if knum in self.ko_dict:
                        ko_threshold = float(self.ko_dict[knum]['threshold'])
                        ko_score_type = self.ko_dict[knum]['score_type']
                    elif self.stray_ko_dict and knum in self.stray_ko_dict:
                        ko_threshold = float(self.stray_ko_dict[knum]['threshold'])
                        ko_score_type = self.stray_ko_dict[knum]['score_type']
                    else:
                        raise ConfigError(f"_apply_heuristic() cannot find bit score threshold for {knum}.")

                    if ko_score_type == 'domain':
                        hit_bitscore = super_hits_dict[hit_label][hit_key]['domain_bit_score']
                        hit_eval = super_hits_dict[hit_label][hit_key]['domain_e_value']
                    else:
                        hit_bitscore = super_hits_dict[hit_label][hit_key]['bit_score']
                        hit_eval = super_hits_dict[hit_label][hit_key]['e_value']

                    if hit_eval <= self.bitscore_heuristic_e_value and \
                            hit_bitscore > (self.bitscore_heuristic_bitscore_fraction * ko_threshold):
                        decent_hit_kos.add(knum)
                        if hit_eval <= best_e_value:
                            best_e_value = hit_eval
                            best_hit_key = hit_key
                            best_hit_label = hit_label

            if len(decent_hit_kos) == 1:
                best_knum = super_hits_dict[best_hit_label][best_hit_key]['gene_name']
                functions_dict[next_key] = {
                    'gene_callers_id': gcid,
                    'source': 'KOfam',
                    'accession': best_knum,
                    'function': self._get_kofam_annotation(best_knum),
                    'e_value': super_hits_dict[best_hit_label][best_hit_key]['e_value'],
                }
                gcids_to_functions_dict[gcid] = [next_key]
                self._add_module_entries(best_knum, gcid, next_key, kegg_module_names_dict, kegg_module_classes_dict)
                if not self.skip_brite_hierarchies:
                    brite_entry = self._build_brite_entry(best_knum, gcid)
                    if brite_entry:
                        kegg_brite_categorizations_dict[next_key] = brite_entry
                if best_hit_label == "No-threshold KO":
                    num_stray_added += 1
                next_key += 1
                num_added += 1

        return next_key, num_added, num_stray_added


    def run(self):
        """Execute the full annotation pipeline for this worker's FASTA partition.

        Loads shared KO/module data from shared memory, runs KOfam (and optionally
        nt-KO) HMM searches, applies bitscore threshold filtering, optionally applies
        the bitscore heuristic, and puts compact result dicts on the result queue.
        On any unhandled exception, logs the failure and puts an error dict on the
        queue instead. Progress at each stage is written to self.log_file_path.
        """

        hmmer = None
        ohmmer = None
        worker_run = terminal.Run(log_file_path=self.log_file_path, verbose=False)
        try:
            self._load_shared_data()

            gcids_in_partition = []
            with open(self.partition_fasta) as fasta_f:
                for line in fasta_f:
                    if line.startswith('>'):
                        gcids_in_partition.append(int(line.strip()[1:]))
            num_seqs = len(gcids_in_partition)
            worker_run.log(f"Worker {self.partition_index} started: {num_seqs} sequences in partition")

            silent_run = terminal.Run(verbose=False)
            silent_progress = terminal.Progress(verbose=False)
            noise_cutoff_terms = "-T -20 --domT -20" if self.no_hmmer_prefiltering else None
            target_files_dict = {'AA:GENE': self.partition_fasta}

            worker_run.log(f"Worker {self.partition_index}: running KOfam HMMs")
            z_param = {'AA:GENE': self.total_seqs}
            hmmer = HMMer(target_files_dict, num_threads_to_use=1, program_to_use=self.hmm_program,
                          run=silent_run, progress=silent_progress,
                          total_number_of_sequences_per_target=z_param)
            hmm_hits_file = hmmer.run_hmmer('KOfam', 'AA', 'GENE', None, None,
                                            len(self.ko_dict), self.kofam_hmm_path, None, noise_cutoff_terms)
            worker_run.log(f"Worker {self.partition_index}: KOfam HMM search complete")

            stray_hits_file = None
            if self.include_stray_kos and self.stray_ko_hmm_path and self.stray_ko_dict:
                worker_run.log(f"Worker {self.partition_index}: running nt-KO HMMs")
                ohmmer = HMMer(target_files_dict, num_threads_to_use=1, program_to_use=self.hmm_program,
                               run=silent_run, progress=silent_progress,
                               total_number_of_sequences_per_target=z_param)
                stray_hits_file = ohmmer.run_hmmer('Stray KOs', 'AA', 'GENE', None, None,
                                                   len(self.stray_ko_dict), self.stray_ko_hmm_path, None, noise_cutoff_terms)
                worker_run.log(f"Worker {self.partition_index}: nt-KO HMM search complete")

            functions_dict = {}
            kegg_module_names_dict = {}
            kegg_module_classes_dict = {}
            kegg_brite_categorizations_dict = {}
            gcids_to_hits_dict = {}
            gcids_to_functions_dict = {}
            super_hits_dict = {}
            next_key = 0
            num_weak_removed = 0

            if hmm_hits_file:
                parser = parser_modules['search']['hmmer_table_output'](
                    hmm_hits_file, alphabet='AA', context='GENE', program=self.hmm_program)
                search_results_dict = parser.get_search_results()
                worker_run.log(f"Worker {self.partition_index}: parsing {len(search_results_dict)} KOfam hits")
                super_hits_dict['KOfam'] = search_results_dict
                next_key, removed = self._parse_hits(
                    search_results_dict, 'KOfam', next_key,
                    functions_dict, kegg_module_names_dict, kegg_module_classes_dict,
                    kegg_brite_categorizations_dict, gcids_to_hits_dict, gcids_to_functions_dict)
                num_weak_removed += removed
                worker_run.log(f"Worker {self.partition_index}: KOfam hits parsed — {removed} weak hits removed, {next_key} annotations added")

            if stray_hits_file:
                oparser = parser_modules['search']['hmmer_table_output'](
                    stray_hits_file, alphabet='AA', context='GENE', program=self.hmm_program)
                stray_results_dict = oparser.get_search_results()
                worker_run.log(f"Worker {self.partition_index}: parsing {len(stray_results_dict)} nt-KO hits")
                super_hits_dict['No-threshold KO'] = stray_results_dict
                prev_key = next_key
                next_key, removed = self._parse_hits(
                    stray_results_dict, 'No-threshold KO', next_key,
                    functions_dict, kegg_module_names_dict, kegg_module_classes_dict,
                    kegg_brite_categorizations_dict, gcids_to_hits_dict, gcids_to_functions_dict)
                num_weak_removed += removed
                worker_run.log(f"Worker {self.partition_index}: nt-KO hits parsed — {removed} weak hits removed, {next_key - prev_key} annotations added")

            num_heuristic_added = 0
            num_stray_heuristic_added = 0
            if not self.keep_all_hits and not self.skip_bitscore_heuristic and super_hits_dict:
                worker_run.log(f"Worker {self.partition_index}: applying bitscore heuristic to {num_seqs} genes")
                next_key, num_heuristic_added, num_stray_heuristic_added = self._apply_heuristic(
                    gcids_in_partition, super_hits_dict, next_key,
                    functions_dict, kegg_module_names_dict, kegg_module_classes_dict,
                    kegg_brite_categorizations_dict, gcids_to_hits_dict, gcids_to_functions_dict)
                worker_run.log(f"Worker {self.partition_index}: bitscore heuristic complete — "
                               f"{num_heuristic_added} annotations added ({num_stray_heuristic_added} from nt-KOs)")

            bitscore_log_file = None
            if self.log_bitscores and 'KOfam' in super_hits_dict:
                bitscore_log_file = self.log_file_path + '_bitscores.txt'
                bitscore_data = super_hits_dict['KOfam']
                for entry in bitscore_data.values():
                    entry['gene_hmm_id'] = entry['gene_name']
                    entry['gene_name'] = self._get_kofam_annotation(entry['gene_hmm_id'])
                utils.store_dict_as_TAB_delimited_file(bitscore_data, bitscore_log_file, do_not_write_key_column=True)

            worker_run.log(f"Worker {self.partition_index} complete: {len(functions_dict)} total annotations, "
                           f"{num_weak_removed} weak hits removed, {num_heuristic_added} added by heuristic")

            self.result_queue.put({
                'success': True,
                'functions_dict': functions_dict,
                'kegg_module_names_dict': kegg_module_names_dict,
                'kegg_module_classes_dict': kegg_module_classes_dict,
                'kegg_brite_categorizations_dict': kegg_brite_categorizations_dict,
                'bitscore_log_file': bitscore_log_file,
                'num_weak_removed': num_weak_removed,
                'num_heuristic_added': num_heuristic_added,
                'num_stray_heuristic_added': num_stray_heuristic_added,
            })

        except Exception as e:
            worker_run.log(f"Worker {self.partition_index} FAILED: {e}")
            self.result_queue.put({'success': False, 'error': e})

        finally:
            if not anvio.DEBUG:
                if hmmer:
                    hmmer.clean_tmp_dirs()
                if ohmmer:
                    ohmmer.clean_tmp_dirs()


class RunKOfams(KeggContext):
    """Class for running `hmmscan` against the KOfam database and adding the resulting hits to contigs DB for later metabolism prediction.

    Parameters
    ==========
    args: Namespace object
        All the arguments supplied by user to anvi-run-kegg-kofams
    """

    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.contigs_db_path = A('contigs_db')
        self.num_threads = A('num_threads') or 1
        self.hmm_program = A('hmmer_program') or 'hmmsearch'
        self.include_stray_kos = True if A('include_nt_KOs') else False
        self.keep_all_hits = True if A('keep_all_hits') else False
        self.log_bitscores = True if A('log_bitscores') else False
        self.skip_bitscore_heuristic = True if A('skip_bitscore_heuristic') else False
        self.no_hmmer_prefiltering = True if A('no_hmmer_prefiltering') else False
        self.bitscore_heuristic_e_value = A('heuristic_e_value')
        self.bitscore_heuristic_bitscore_fraction = A('heuristic_bitscore_fraction')
        self.skip_brite_hierarchies = A('skip_brite_hierarchies')
        self.ko_dict = None # should be set up by setup_ko_dict()
        self.stray_ko_dict = None # should be set up by setup_stray_ko_dict(), if possible

        # init the base class
        KeggContext.__init__(self, self.args)

        filesnpaths.is_program_exists(self.hmm_program)

        # verify that Kofam HMM profiles have been set up
        if not os.path.exists(self.kofam_hmm_file_path):
            raise ConfigError(f"Anvi'o is unable to find any KEGG files around :/ It is likely you need to first run the program "
                              f"`anvi-setup-kegg-data` to set things up. If you already have run it, but instructed anvi'o to "
                              f"store the output to a specific directory, then instead of running `anvi-setup-kegg-data` again, "
                              f"you simply need to specify the location of the KEGG data using the flag `--kegg-data-dir`. Just for "
                              f"your information, anvi'o was looking for the KEGG data here: {self.kegg_data_dir}")

        utils.is_contigs_db(self.contigs_db_path)
        filesnpaths.is_output_file_writable(self.contigs_db_path)

        # reminder to be a good citizen
        self.run.warning("Anvi'o will annotate your database with the KEGG KOfam database, as described in "
                         "Aramaki et al (doi:10.1093/bioinformatics/btz859) When you publish your findings, "
                         "please do not forget to properly credit this work.", lc='green', header="CITATION")

        self.setup_ko_dict() # read the ko_list file into self.ko_dict
        self.run.info("nt-KOs will be annotated", self.include_stray_kos)
        if self.include_stray_kos:
            self.setup_stray_ko_dict()
            self.run.warning("Please note! Because you used the flag `--include-nt-KOs`, anvi'o will annotate "
                             "your genes with KO models that do not come with a bit score threshold defined by KEGG. "
                             "We have generated new models and estimated rather conservative thresholds for them ourselves. To learn "
                             "how we did that, please read the documentation page for `anvi-setup-kegg-data`: "
                             "https://anvio.org/help/main/programs/anvi-setup-kegg-data/#what-are-nt-kos-and-what-happens-when-i-include-them")

        # load existing kegg modules db, if one exists
        if os.path.exists(self.kegg_modules_db_path):
            self.kegg_modules_db = ModulesDatabase(self.kegg_modules_db_path, module_data_directory=self.kegg_module_data_dir,
                                                   brite_data_directory=self.brite_data_dir, skip_brite_hierarchies=self.skip_brite_hierarchies,
                                                   args=self.args)

            if not self.skip_brite_hierarchies and not self.kegg_modules_db.db.get_meta_value('is_brite_setup'):
                self.run.warning("The KEGG Modules database does not contain BRITE hierarchy data, "
                                 "which could very well be useful to you. BRITE is guaranteed to be set up "
                                 "when downloading the latest version of KEGG with `anvi-setup-kegg-data`.")
        else:
            self.run.warning("No modules database was found in the KEGG data directory you specified. This is fine, but "
                             "you will not get functional annotations related to KEGG MODULES or BRITE hierarchies in your "
                             "contigs database. If you want to include these annotations later, you will have to rerun this "
                             "program with a data directory including a modules database (which you can obtain by running "
                             "`anvi-setup-kegg-data` again with the right mode(s).")
            self.kegg_modules_db = None


    def check_hash_in_contigs_db(self):
        """Checks the contigs DB self table to make sure it was not already annotated"""

        A = lambda x: self.args.__dict__[x] if x in self.args.__dict__ else None
        self.contigs_db_path = A('contigs_db')

        contigs_db = ContigsDatabase(self.contigs_db_path)
        current_module_hash_in_contigs_db = contigs_db.db.get_meta_value('modules_db_hash', return_none_if_not_in_table=True)

        if current_module_hash_in_contigs_db and not self.just_do_it:
            contigs_db.disconnect()
            raise ConfigError("The contigs database (%s) has already been annotated with KOfam hits. If you really want to "
                              "overwrite these annotations with new ones, please re-run the command with the flag --just-do-it. "
                              "For those who need this information, the Modules DB used to annotate this contigs database previously "
                              "had the following hash: %s" % (self.contigs_db_path, current_module_hash_in_contigs_db))


    def set_hash_in_contigs_db(self):
        """Modifies the contigs DB self table to indicate which MODULES.db has been used to annotate it."""

        A = lambda x: self.args.__dict__[x] if x in self.args.__dict__ else None
        self.contigs_db_path = A('contigs_db')

        hash_to_add = "only_KOfams_were_annotated"
        if self.kegg_modules_db:
            hash_to_add = self.kegg_modules_db.db.get_meta_value('hash')

        contigs_db = ContigsDatabase(self.contigs_db_path)
        contigs_db.db.set_meta_value('modules_db_hash', hash_to_add)
        contigs_db.disconnect()


    def store_annotations_in_db(self):
        """Takes the dictionary of function annotations (already parsed, if necessary) and puts them in the DB.

        Should be called after the function that parses the HMM hits and creates self.functions_dict :) which is
        parse_kofam_hits()
        """

        # get an instance of gene functions table
        gene_function_calls_table = TableForGeneFunctions(self.contigs_db_path, self.run, self.progress)

        if self.functions_dict:
            gene_function_calls_table.create(self.functions_dict)
            if self.kegg_module_names_dict:
                gene_function_calls_table.create(self.kegg_module_names_dict)
            if self.kegg_module_classes_dict:
                gene_function_calls_table.create(self.kegg_module_classes_dict)
            if self.kegg_brite_categorizations_dict:
                gene_function_calls_table.create(self.kegg_brite_categorizations_dict)
        else:
            self.run.warning("There are no KOfam hits to add to the database. Returning empty handed, "
                             "but still adding KOfam as a functional source.")
            gene_function_calls_table.add_empty_sources_to_functional_sources({'KOfam'})


    def _prefetch_modules_data(self):
        """Build module and BRITE lookup dicts via bulk DB queries for use by annotation workers."""

        ko_to_modules = {}
        ko_to_module_names = {}
        ko_to_module_classes = {}
        ko_to_brite = {}

        if not self.kegg_modules_db:
            return ko_to_modules, ko_to_module_names, ko_to_module_classes, ko_to_brite

        mdb = self.kegg_modules_db

        # Build knum → [module, ...] from all ORTHOLOGY rows in one query
        orth_rows = mdb.db.get_some_rows_from_table_as_dict(
            mdb.module_table_name, "data_name = 'ORTHOLOGY'",
            row_num_as_key=True, error_if_no_data=False)
        for row in orth_rows.values():
            knum = row['data_value']
            mnum = row['module']
            if knum not in ko_to_modules:
                ko_to_modules[knum] = []
            if mnum not in ko_to_modules[knum]:
                ko_to_modules[knum].append(mnum)

        # Bulk fetch module names and class strings
        name_rows = mdb.db.get_some_rows_from_table_as_dict(
            mdb.module_table_name, "data_name = 'NAME'",
            row_num_as_key=True, error_if_no_data=False)
        module_names = {row['module']: row['data_value'] for row in name_rows.values()}

        class_rows = mdb.db.get_some_rows_from_table_as_dict(
            mdb.module_table_name, "data_name = 'CLASS'",
            row_num_as_key=True, error_if_no_data=False)
        module_classes = {row['module']: row['data_value'] for row in class_rows.values()}

        for knum, mnums in ko_to_modules.items():
            ko_to_module_names[knum] = {mnum: module_names.get(mnum, '') for mnum in mnums}
            ko_to_module_classes[knum] = [module_classes.get(mnum, '') for mnum in mnums]

        # Bulk fetch all BRITE rows and group by ortholog accession
        if not self.skip_brite_hierarchies:
            brite_data = mdb.db.get_table_as_dict(
                mdb.brite_table_name, row_num_as_key=True, error_if_no_data=False)
            for row in brite_data.values():
                knum = row['ortholog_accession']
                if knum not in ko_to_brite:
                    ko_to_brite[knum] = []
                ko_to_brite[knum].append(row)

        return ko_to_modules, ko_to_module_names, ko_to_module_classes, ko_to_brite


    def process_kofam_hmms(self):
        """Driver for running HMMs against the KOfam database and processing hits into the contigs DB."""

        from multiprocessing.shared_memory import SharedMemory

        tmp_directory_path = filesnpaths.get_temp_directory_path()
        contigs_db = ContigsSuperclass(self.args)

        self.check_hash_in_contigs_db()

        # Export all AA sequences to a single FASTA
        aa_fasta_path = os.path.join(tmp_directory_path, 'AA_gene_sequences.fa')
        contigs_db.get_sequences_for_gene_callers_ids(
            output_file_path=aa_fasta_path,
            simple_headers=True,
            report_aa_sequences=True)

        # Split the FASTA once; each worker gets one partition
        split_dir = filesnpaths.get_temp_directory_path()
        partition_files, total_seqs = utils.split_fasta(
            aa_fasta_path,
            parts=self.num_threads,
            output_dir=split_dir,
            return_number_of_sequences=True)
        num_partitions = len(partition_files)

        self.run.info("Number of amino acid sequences", total_seqs)
        self.run.info("Number of partitions for parallel processing", num_partitions)

        # Pre-fetch all module and BRITE lookup data via bulk DB queries
        ko_to_modules, ko_to_module_names, ko_to_module_classes, ko_to_brite = self._prefetch_modules_data()

        if not self.keep_all_hits and not self.skip_bitscore_heuristic:
            self.run.warning(
                "Anvi'o will re-visit genes without KOfam annotations to see if potentially valid "
                "functional annotations were missed. These genes will be annotated with a KO only if "
                f"all KOfam hits to this gene with e-value <= {self.bitscore_heuristic_e_value} and bitscore > "
                f"({self.bitscore_heuristic_bitscore_fraction} * KEGG threshold) are hits to the same KO. Just "
                "so you know what is going on here. If this sounds like A Very Bad Idea to you, then please "
                "feel free to turn off this behavior with the flag --skip-bitscore-heuristic or to change "
                "the e-value/bitscore parameters (see the help page for more info).")

        # Pack all read-only shared data into a single pickle blob in shared memory.
        # Workers attach and deserialize once at startup — no per-argument pickling overhead.
        shared_bundle = {
            'ko_dict': self.ko_dict,
            'stray_ko_dict': self.stray_ko_dict,
            'ko_to_modules': ko_to_modules,
            'ko_to_module_names': ko_to_module_names,
            'ko_to_module_classes': ko_to_module_classes,
            'ko_to_brite': ko_to_brite,
            'total_seqs': total_seqs,
        }
        bundle_bytes = pickle.dumps(shared_bundle)
        shared_data_shm = SharedMemory(create=True, size=max(len(bundle_bytes), 1))
        shared_data_shm.buf[:len(bundle_bytes)] = bundle_bytes
        shared_data_shm_name = shared_data_shm.name
        shared_data_size = len(bundle_bytes)
        del bundle_bytes  # free local copy; data now lives only in shared memory

        # Spawn one annotation worker per partition
        manager = multiprocessing.Manager()
        result_queue = manager.Queue()

        workers = []
        for i, partition_fasta in enumerate(partition_files):
            log_file = os.path.join(tmp_directory_path, f'worker_{i}.log')
            worker = AnnotationWorker(
                partition_fasta=partition_fasta,
                partition_index=i,
                hmm_program=self.hmm_program,
                kofam_hmm_path=self.kofam_hmm_file_path,
                stray_ko_hmm_path=self.stray_ko_hmm_file_path if self.include_stray_kos else None,
                include_stray_kos=self.include_stray_kos,
                keep_all_hits=self.keep_all_hits,
                skip_bitscore_heuristic=self.skip_bitscore_heuristic,
                no_hmmer_prefiltering=self.no_hmmer_prefiltering,
                log_bitscores=self.log_bitscores,
                skip_brite_hierarchies=self.skip_brite_hierarchies,
                bitscore_heuristic_e_value=self.bitscore_heuristic_e_value,
                bitscore_heuristic_bitscore_fraction=self.bitscore_heuristic_bitscore_fraction,
                log_file_path=log_file,
                result_queue=result_queue,
                shared_data_shm_name=shared_data_shm_name,
                shared_data_size=shared_data_size)
            worker.start()
            workers.append(worker)

        # Collect results as workers finish
        self.progress.new('Annotating with KOfam HMMs', progress_total_items=num_partitions)
        self.progress.update(f'Waiting for {num_partitions} annotation worker(s)...')

        all_results = []
        bitscore_log_files = []
        total_weak_removed = 0
        total_heuristic_added = 0
        total_stray_heuristic_added = 0

        finished = 0
        try:
            while finished < num_partitions:
                result = result_queue.get()
                finished += 1

                if not result['success']:
                    raise result['error']

                self.progress.update(f'Finished {finished} of {num_partitions} partition(s)')
                self.progress.increment(increment_to=finished)
                all_results.append(result)

                if result.get('bitscore_log_file'):
                    bitscore_log_files.append(result['bitscore_log_file'])
                total_weak_removed += result.get('num_weak_removed', 0)
                total_heuristic_added += result.get('num_heuristic_added', 0)
                total_stray_heuristic_added += result.get('num_stray_heuristic_added', 0)

        finally:
            # always end progress and terminate workers, even on error
            self.progress.end()
            for w in workers:
                w.terminate()
            # release and unlink shared memory segment
            try:
                shared_data_shm.close()
                shared_data_shm.unlink()
            except Exception as shm_err:
                self.run.warning(f"Failed to clean up shared memory segment '{shared_data_shm_name}': {shm_err}")

        # Merge annotation dicts across partitions with non-overlapping integer keys
        self.functions_dict = {}
        self.kegg_module_names_dict = {}
        self.kegg_module_classes_dict = {}
        self.kegg_brite_categorizations_dict = {}

        key_offset = 0
        for result in all_results:
            fd = result['functions_dict']
            for k, v in fd.items():
                self.functions_dict[k + key_offset] = v
            for k, v in result['kegg_module_names_dict'].items():
                self.kegg_module_names_dict[k + key_offset] = v
            for k, v in result['kegg_module_classes_dict'].items():
                self.kegg_module_classes_dict[k + key_offset] = v
            for k, v in result['kegg_brite_categorizations_dict'].items():
                self.kegg_brite_categorizations_dict[k + key_offset] = v
            # worker keys are 0..N-1, so len(fd) advances the offset past this partition
            key_offset += len(fd)

        # Report aggregate stats
        self.run.info("Total number of annotations", len(self.functions_dict))
        self.run.info("Total weak hits removed", total_weak_removed)
        if not self.keep_all_hits and not self.skip_bitscore_heuristic:
            self.run.info("Annotations added via bitscore heuristic", total_heuristic_added)
            if self.include_stray_kos:
                self.run.info("... of which are from nt-KOs", total_stray_heuristic_added)

        if not self.functions_dict:
            self.run.info_single("The HMM search returned no hits :/ So there is nothing to add to the contigs "
                                 "database. But now anvi'o will add KOfam as a functional source with no hits, "
                                 "clean the temporary directories and gracefully quit.", nl_before=1, nl_after=1)

        # Store annotations in DB
        self.store_annotations_in_db()

        # Write bitscore log — each worker wrote a partition file; concatenate them here
        if self.log_bitscores:
            self.bitscore_log_file = os.path.splitext(os.path.basename(self.contigs_db_path))[0] + "_bitscores.txt"
            if bitscore_log_files:
                with open(self.bitscore_log_file, 'w') as out_f:
                    for i, log_path in enumerate(bitscore_log_files):
                        with open(log_path, 'r') as in_f:
                            if i == 0:
                                out_f.write(in_f.read())
                            else:
                                lines = in_f.readlines()
                                if len(lines) > 1:
                                    out_f.writelines(lines[1:])
                self.run.info("Bit score information file", self.bitscore_log_file)
            else:
                self.run.info_single("No KOfam hits found; no bitscore log written.")

        # Mark contigs DB with hash of modules DB content for version tracking
        self.set_hash_in_contigs_db()

        if anvio.DEBUG:
            self.run.warning(f"The temp directories '{tmp_directory_path}' and '{split_dir}' are kept. "
                            "Please don't forget to clean those up later.", header="Debug")
        else:
            self.run.info_single("Cleaning up the temp directory (you can use `--debug` if you would "
                                "like to keep it for testing purposes)", nl_before=1, nl_after=1)
            shutil.rmtree(tmp_directory_path)
            shutil.rmtree(split_dir)
