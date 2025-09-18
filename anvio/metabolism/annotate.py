import os
import shutil

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
        self.num_threads = A('num_threads')
        self.hmm_program = A('hmmer_program') or 'hmmsearch'
        self.include_stray_kos = True if A('include_stray_KOs') else False
        self.keep_all_hits = True if A('keep_all_hits') else False
        self.log_bitscores = True if A('log_bitscores') else False
        self.skip_bitscore_heuristic = True if A('skip_bitscore_heuristic') else False
        self.no_hmmer_prefiltering = True if A('no_hmmer_prefiltering') else False
        self.bitscore_heuristic_e_value = A('heuristic_e_value')
        self.bitscore_heuristic_bitscore_fraction = A('heuristic_bitscore_fraction')
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
        self.run.info("Stray KOs will be annotated", self.include_stray_kos)
        if self.include_stray_kos:
            self.setup_stray_ko_dict()
            self.run.warning("Please note! Because you used the flag `--include-stray-KOs`, anvi'o will annotate "
                             "your genes with KO models that do not come with a bit score threshold defined by KEGG. "
                             "We have generated new models and estimated rather conservative thresholds for them ourselves. To learn "
                             "how we did that, please read the documentation page for `anvi-setup-kegg-data`: "
                             "https://anvio.org/help/main/programs/anvi-setup-kegg-data/#what-are-stray-kos-and-what-happens-when-i-include-them")

        # load existing kegg modules db, if one exists
        if os.path.exists(self.kegg_modules_db_path):
            self.kegg_modules_db = ModulesDatabase(self.kegg_modules_db_path, module_data_directory=self.kegg_module_data_dir,
                                                   args=self.args)

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


    def get_annotation_from_ko_dict(self, knum, ok_if_missing_from_dict=False):
        """Returns the functional annotation of the provided KO number.

        Parameters
        ==========
        knum : str
            The KO number for which to get an annotation for
        ok_if_missing_from_dict : bool
            If false, not finding the KO will raise an error. If true, the function will quietly return an "Unknown" annotation string for the missing KO

        Returns
        =======
        annotation : str
        """

        if not self.ko_dict:
            raise ConfigError("Oops! The ko_list file has not been properly loaded, so get_annotation_from_ko_dict() is "
                              "extremely displeased and unable to function properly. Please refrain from calling this "
                              "function until after setup_ko_dict() has been called.")
        if self.include_stray_kos and not self.stray_ko_dict:
            raise ConfigError("Oops! The bit score thresholds for stray KOs have not been properly loaded, so "
                              "get_annotation_from_ko_dict() is unable to work properly. If you plan to use "
                              "--include-stray-KOs, then make sure you run the setup_stray_ko_dict() function before "
                              "calling this one.")

        ret_value = None
        if knum in self.ko_dict:
            ret_value = self.ko_dict[knum]['definition']
        elif self.include_stray_kos and knum in self.stray_ko_dict:
            ret_value = self.stray_ko_dict[knum]['definition']
        else:
            if ok_if_missing_from_dict:
                return "Unknown function with KO num %s" % knum
            else:
                raise ConfigError("It seems %s found a KO number that does not exist "
                                  "in the KOfam ko_list file: %s" % (self.hmm_program, knum))

        return ret_value


    def parse_kofam_hits(self, hits_dict, hits_label = "KOfam", next_key=None):
        """This function applies bitscore thresholding (if requested) to establish the self.functions_dict
        which can then be used to store annotations in the contigs DB.

        If self.keep_all_hits is True, all hits will be added to the self.functions_dict regardless of bitscore
        threshold.

        Note that the input hits_dict contains bitscores, but the self.functions_dict does not (because the DB
        tables do not have a column for it, at least at the time of writing this).

        PARAMETERS
        ===========
        hits_dict : dictionary
            The output from the hmmsearch parser, which should contain all hits (ie, weak hits not yet removed)
        hits_label : str
            A label for the set of hits we are working on, used to keep sets separate from each other and to enable
            us to match gene caller ids to hits in different dictionaries later
        next_key : int
            The next integer key that is available for adding functions to self.functions_dict. If None is provided,
            the keys will start at 0.

        RETURNS
        ========
        counter : int
            The number of functions added to self.functions_dict. Useful for downstream functions that want to
            add to this dictionary, since it is the next available integer key.
        """

        total_num_hits = len(hits_dict.values())
        starting_annotations_in_dict = len(self.functions_dict.keys())
        self.progress.new(f"Parsing {hits_label} hits", progress_total_items=total_num_hits)
        counter = 0
        if next_key:
            counter = next_key
        num_hits_removed = 0
        cur_num_hit = 0
        for hit_key,hmm_hit in hits_dict.items():
            cur_num_hit += 1
            knum = hmm_hit['gene_name']
            gcid = hmm_hit['gene_callers_id']
            keep = False

            if cur_num_hit % 1000 == 0:
                self.progress.update("Removing weak hits [%d of %d KOs]" % (cur_num_hit, total_num_hits))
                self.progress.increment(increment_to=cur_num_hit)

            # later, we will need to quickly access the hits for each gene call. So we map gcids to the keys in the raw hits dictionary
            if gcid not in self.gcids_to_hits_dict:
                self.gcids_to_hits_dict[gcid] = {hits_label : [hit_key]}
            else:
                if hits_label not in self.gcids_to_hits_dict[gcid]:
                    self.gcids_to_hits_dict[gcid][hits_label] = [hit_key]
                else:
                    self.gcids_to_hits_dict[gcid][hits_label].append(hit_key)

            if (knum not in self.ko_dict and (self.stray_ko_dict is not None and knum not in self.stray_ko_dict)) or \
                (knum not in self.ko_dict and self.stray_ko_dict is None):
                self.progress.reset()
                raise ConfigError(f"Something went wrong while parsing the {hits_label} HMM hits. It seems that KO "
                                  f"{knum} is not in the noise cutoff dictionary for KOs. That means we do "
                                  "not know how to distinguish strong hits from weak ones for this KO. "
                                  "Anvi'o will fail now :( Please contact a developer about this error to "
                                  "get this mess fixed. ")
            # if hit is above the bitscore threshold, we will keep it
            if knum in self.ko_dict:
                if self.ko_dict[knum]['score_type'] == 'domain':
                    if hmm_hit['domain_bit_score'] >= float(self.ko_dict[knum]['threshold']):
                        keep = True
                elif self.ko_dict[knum]['score_type'] == 'full':
                    if hmm_hit['bit_score'] >= float(self.ko_dict[knum]['threshold']):
                        keep = True
                else:
                    self.progress.reset()
                    raise ConfigError(f"The KO noise cutoff dictionary for {knum} has a strange score type which "
                                    f"is unknown to anvi'o: {self.ko_dict[knum]['score_type']}")
            elif knum in self.stray_ko_dict:
                if self.stray_ko_dict[knum]['score_type'] == 'domain':
                    if hmm_hit['domain_bit_score'] >= float(self.stray_ko_dict[knum]['threshold']):
                        keep = True
                elif self.stray_ko_dict[knum]['score_type'] == 'full':
                    if hmm_hit['bit_score'] >= float(self.stray_ko_dict[knum]['threshold']):
                        keep = True
                else:
                    self.progress.reset()
                    raise ConfigError(f"The KO noise cutoff dictionary for the stray KO {knum} has a strange score type which "
                                    f"is unknown to anvi'o: {self.stray_ko_dict[knum]['score_type']}")
            else:
                raise ConfigError(f"We cannot find KO {knum} in either self.ko_dict or in self.stray_ko_dict. This is likely a "
                                  f"problem for the developers.")

            if keep or self.keep_all_hits:
                self.functions_dict[counter] = {
                    'gene_callers_id': gcid,
                    'source': 'KOfam',
                    'accession': knum,
                    'function': self.get_annotation_from_ko_dict(knum, ok_if_missing_from_dict=True),
                    'e_value': hmm_hit['e_value'],
                }
                # later, we will need to know if a particular gene call has hits or not. So here we are just saving for each
                # gene caller id the keys for its corresponding hits in the function dictionary.
                if gcid not in self.gcids_to_functions_dict:
                    self.gcids_to_functions_dict[gcid] = [counter]
                else:
                    self.gcids_to_functions_dict[gcid].append(counter)

                # add associated KEGG module information to database
                mods = None
                if self.kegg_modules_db:
                    mods = self.kegg_modules_db.get_modules_for_knum(knum)
                    names = self.kegg_modules_db.get_module_names_for_knum(knum)
                    classes = self.kegg_modules_db.get_module_classes_for_knum_as_list(knum)

                if mods:
                    mod_annotation = "!!!".join(mods)
                    mod_class_annotation = "!!!".join(classes) # why do we split by '!!!'? Because that is how it is done in COGs. So so sorry. :'(
                    mod_name_annotation = ""

                    for mod in mods:
                        if mod_name_annotation:
                            mod_name_annotation += "!!!" + names[mod]
                        else:
                            mod_name_annotation = names[mod]

                    self.kegg_module_names_dict[counter] = {
                        'gene_callers_id': gcid,
                        'source': 'KEGG_Module',
                        'accession': mod_annotation,
                        'function': mod_name_annotation,
                        'e_value': None,
                    }
                    self.kegg_module_classes_dict[counter] = {
                        'gene_callers_id': gcid,
                        'source': 'KEGG_Class',
                        'accession': mod_annotation,
                        'function': mod_class_annotation,
                        'e_value': None,
                    }

                counter += 1
            else:
                num_hits_removed += 1

        self.progress.end()
        ending_annotations_in_dict = len(self.functions_dict.keys())
        self.run.info(f"Number of weak hits removed by {hits_label} parser", num_hits_removed)
        self.run.info(f"Number of annotations added for {hits_label}", ending_annotations_in_dict - starting_annotations_in_dict)

        return counter


    def update_dict_for_genes_with_missing_annotations(self, gcids_list, super_hits_dict, next_key):
        """This function adds functional annotations for genes with missing hits to the dictionary.

        The reason this is necessary is that the bitscore thresholds can be too stringent, causing
        us to miss legitimate annotations. To find these annotations, we adopt the following heuristic:
            For every gene without a KOfam annotation, we examine all the hits with an e-value below X
            and a bitscore above Y% of the threshold. If those hits are all to a unique KO profile,
            then we annotate the gene call with that KO.

            X is self.bitscore_heuristic_e_value, Y is self.bitscore_heuristic_bitscore_fraction

        For reasons that are hopefully obvious, this function must be called after parse_kofam_hits(),
        which establishes the self.functions_dict attribute.

        PARAMETERS
        ===========
        gcids_list : list
            The list of gene caller ids in the contigs database. We will use this to figure out which
            genes have no annotations
        super_hits_dict : dictionary
            A two-level dictionary in which keys are the labels for each set of hits and values are the dictionary output
            from the hmmsearch parser, which should contain all hits from the set (ie, weak hits not yet removed)
        next_key : int
            The next integer key that is available for adding functions to self.functions_dict
        """

        self.run.warning("Anvi'o will now re-visit genes without KOfam annotations to see if potentially valid "
                         "functional annotations were missed. These genes will be annotated with a KO only if "
                         f"all KOfam hits to this gene with e-value <= {self.bitscore_heuristic_e_value} and bitscore > "
                         f"({self.bitscore_heuristic_bitscore_fraction} * KEGG threshold) are hits to the same KO. Just "
                         "so you know what is going on here. If this sounds like A Very Bad Idea to you, then please "
                         "feel free to turn off this behavior with the flag --skip-bitscore-heuristic or to change "
                         "the e-value/bitscore parameters (see the help page for more info).")

        num_annotations_added = 0
        num_stray_KOs_added = 0
        total_num_genes = len(gcids_list)
        self.progress.new("Relaxing bitscore threshold", progress_total_items=total_num_genes)

        # for each gene call, check for annotation in self.functions_dict
        current_gene_num = 0
        for gcid in gcids_list:
            current_gene_num += 1
            if current_gene_num % 1000 == 0:
                self.progress.update("Adding back decent hits [%d of %d gene calls]" % (current_gene_num, total_num_genes))
                self.progress.increment(increment_to=current_gene_num)

            if gcid not in self.gcids_to_functions_dict:
                decent_hit_kos = set()
                best_e_value = 100 # just an arbitrary positive value that will be larger than any evalue
                best_hit_key = None
                best_hit_label = None

                # if no annotation, get all hits for gene caller id from the hits dictionaries
                if gcid in self.gcids_to_hits_dict:
                    for hit_label in self.gcids_to_hits_dict[gcid]:
                        for hit_key in self.gcids_to_hits_dict[gcid][hit_label]:
                            knum = super_hits_dict[hit_label][hit_key]['gene_name']
                            ko_threshold = None
                            ko_score_type = None
                            if knum in self.ko_dict:
                                ko_threshold = float(self.ko_dict[knum]['threshold'])
                                ko_score_type = self.ko_dict[knum]['score_type']
                            elif self.include_stray_kos and knum in self.stray_ko_dict:
                                ko_threshold = float(self.stray_ko_dict[knum]['threshold'])
                                ko_score_type = self.stray_ko_dict[knum]['score_type']
                            else:
                                raise ConfigError(f"panik. the function update_dict_for_genes_with_missing_annotations() "
                                                  f"cannot find the bit score threshold for {knum}.")

                            # get set of hits that fit specified heuristic parameters
                            if ko_score_type == 'domain':
                                hit_bitscore = super_hits_dict[hit_label][hit_key]['domain_bit_score']
                                hit_eval = super_hits_dict[hit_label][hit_key]['domain_e_value']
                            elif ko_score_type == 'full':
                                hit_bitscore = super_hits_dict[hit_label][hit_key]['bit_score']
                                hit_eval = super_hits_dict[hit_label][hit_key]['e_value']
                            if hit_eval <= self.bitscore_heuristic_e_value and hit_bitscore > (self.bitscore_heuristic_bitscore_fraction * ko_threshold):
                                decent_hit_kos.add(knum)
                                # keep track of hit with lowest e-value we've seen so far
                                if hit_eval <= best_e_value:
                                    best_e_value = hit_eval
                                    best_hit_key = hit_key
                                    best_hit_label = hit_label

                    # if unique KO, add annotation with best e-value to self.functions_dict
                    if len(decent_hit_kos) == 1:
                        best_knum = super_hits_dict[best_hit_label][best_hit_key]['gene_name']
                        ## TODO: WE NEED A GENERIC FUNCTION FOR THIS SINCE IT IS SAME AS ABOVE
                        self.functions_dict[next_key] = {
                            'gene_callers_id': gcid,
                            'source': 'KOfam',
                            'accession': best_knum,
                            'function': self.get_annotation_from_ko_dict(best_knum, ok_if_missing_from_dict=True),
                            'e_value': super_hits_dict[best_hit_label][best_hit_key]['e_value'],
                        }
                        # we may never access this downstream but let's add to it to be consistent
                        self.gcids_to_functions_dict[gcid] = [next_key]

                        # track how many of stray KOs are added back
                        if best_hit_label == "Stray KO":
                            num_stray_KOs_added += 1

                        # add associated KEGG module information to database
                        mods = None
                        if self.kegg_modules_db:
                            mods = self.kegg_modules_db.get_modules_for_knum(best_knum)
                            names = self.kegg_modules_db.get_module_names_for_knum(best_knum)
                            classes = self.kegg_modules_db.get_module_classes_for_knum_as_list(best_knum)

                        if mods:
                            mod_annotation = "!!!".join(mods)
                            mod_class_annotation = "!!!".join(classes) # why do we split by '!!!'? Because that is how it is done in COGs. So so sorry. :'(
                            mod_name_annotation = ""

                            for mod in mods:
                                if mod_name_annotation:
                                    mod_name_annotation += "!!!" + names[mod]
                                else:
                                    mod_name_annotation = names[mod]

                            self.kegg_module_names_dict[next_key] = {
                                'gene_callers_id': gcid,
                                'source': 'KEGG_Module',
                                'accession': mod_annotation,
                                'function': mod_name_annotation,
                                'e_value': None,
                            }
                            self.kegg_module_classes_dict[next_key] = {
                                'gene_callers_id': gcid,
                                'source': 'KEGG_Class',
                                'accession': mod_annotation,
                                'function': mod_class_annotation,
                                'e_value': None,
                            }

                        next_key += 1
                        num_annotations_added += 1

        self.progress.end()
        self.run.info("Number of decent hits added back after relaxing bitscore threshold", num_annotations_added)
        if self.include_stray_kos:
            self.run.info("... of these, number of regular KOs is", num_annotations_added - num_stray_KOs_added)
            self.run.info("... of these, number of stray KOs is", num_stray_KOs_added)
        self.run.info("Total number of hits in annotation dictionary after adding these back", len(self.functions_dict.keys()))


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


    def process_kofam_hmms(self):
        """This is a driver function for running HMMs against the KOfam database and processing the hits into the provided contigs DB."""

        tmp_directory_path = filesnpaths.get_temp_directory_path()
        contigs_db = ContigsSuperclass(self.args) # initialize contigs db
        # we will need the gene caller ids later
        all_gcids_in_contigs_db = contigs_db.genes_in_contigs_dict.keys()

        # safety check for previous annotations so that people don't overwrite those if they don't want to
        self.check_hash_in_contigs_db()

        # get AA sequences as FASTA
        target_files_dict = {'AA:GENE': os.path.join(tmp_directory_path, 'AA_gene_sequences.fa')}
        contigs_db.get_sequences_for_gene_callers_ids(output_file_path=target_files_dict['AA:GENE'],
                                                      simple_headers=True,
                                                      report_aa_sequences=True)

        # turn off HMMER's default reporting thresholds if requested. We use an extremely low bitscore threshold instead.
        noise_cutoff_terms = None
        if self.no_hmmer_prefiltering:
            noise_cutoff_terms = "-T -20 --domT -20"
        # run hmmscan
        hmmer = HMMer(target_files_dict, num_threads_to_use=self.num_threads, program_to_use=self.hmm_program)
        hmm_hits_file = hmmer.run_hmmer('KOfam', 'AA', 'GENE', None, None, len(self.ko_dict), self.kofam_hmm_file_path, None, noise_cutoff_terms)

        has_stray_hits = False
        stray_hits_file = None
        if self.include_stray_kos:
            ohmmer = HMMer(target_files_dict, num_threads_to_use=self.num_threads, program_to_use=self.hmm_program)
            stray_hits_file = ohmmer.run_hmmer('Stray KOs', 'AA', 'GENE', None, None, len(self.stray_ko_dict), self.stray_ko_hmm_file_path, None, noise_cutoff_terms)
            has_stray_hits = True if stray_hits_file else False

        if not hmm_hits_file and not has_stray_hits:
            self.run.info_single("The HMM search returned no hits :/ So there is nothing to add to the contigs database. But "
                             "now anvi'o will add KOfam as a functional source with no hits, clean the temporary directories "
                             "and gracefully quit.", nl_before=1, nl_after=1)
            if not anvio.DEBUG:
                shutil.rmtree(tmp_directory_path)
                hmmer.clean_tmp_dirs()
            else:
                self.run.warning("Because you ran this script with the --debug flag, anvi'o will not clean up the temporary "
                                 "directories located at %s and %s. Please be responsible for cleaning up this directory yourself "
                                 "after you are finished debugging :)" % (tmp_directory_path, ', '.join(hmmer.tmp_dirs)), header="Debug")
            gene_function_calls_table = TableForGeneFunctions(self.contigs_db_path, self.run, self.progress)
            gene_function_calls_table.add_empty_sources_to_functional_sources({'KOfam'})
            return

        # set up some attributes that we'll need later
        self.functions_dict = {}
        self.kegg_module_names_dict = {}
        self.kegg_module_classes_dict = {}
        self.kegg_brite_categorizations_dict = {}
        self.gcids_to_hits_dict = {}
        self.gcids_to_functions_dict = {}
        super_hits_dict = {} # will store the hits from each set of HMMs

        # parse hmmscan output
        self.run.warning('', header='HMM hit parsing for KOfams', lc='green')
        self.run.info("HMM output table", hmm_hits_file)
        parser = parser_modules['search']['hmmer_table_output'](hmm_hits_file, alphabet='AA', context='GENE', program=self.hmm_program)
        search_results_dict = parser.get_search_results()
        next_key_in_functions_dict = self.parse_kofam_hits(search_results_dict)
        super_hits_dict["KOfam"] = search_results_dict
        self.run.info("Current number of annotations in functions dictionary", len(self.functions_dict))

        if has_stray_hits:
            self.run.warning('', header='HMM hit parsing for Stray KOs', lc='green')
            self.run.info("HMM output table", stray_hits_file)
            oparser = parser_modules['search']['hmmer_table_output'](stray_hits_file, alphabet='AA', context='GENE', program=self.hmm_program)
            stray_search_results = oparser.get_search_results()
            next_key_in_functions_dict = self.parse_kofam_hits(stray_search_results, hits_label = "Stray KO", next_key=next_key_in_functions_dict)
            super_hits_dict["Stray KO"] = stray_search_results
            self.run.info("Current number of annotations in functions dictionary", len(self.functions_dict))

        if not self.skip_bitscore_heuristic:
            self.update_dict_for_genes_with_missing_annotations(all_gcids_in_contigs_db, super_hits_dict, next_key=next_key_in_functions_dict)

        # add functions and KEGG modules info to database
        self.store_annotations_in_db()

        # If requested, store bit scores of each hit in file
        if self.log_bitscores:
            self.bitscore_log_file = os.path.splitext(os.path.basename(self.contigs_db_path))[0] + "_bitscores.txt"
            # we have to change some things in the raw search results because the KOfam models don't store the gene functions
            # and the HMM ID is actually stored in a model's gene name element
            for key, entry in search_results_dict.items():
                entry['gene_hmm_id'] = entry['gene_name']
                entry['gene_name'] = self.get_annotation_from_ko_dict(entry['gene_hmm_id'], ok_if_missing_from_dict=True)
            anvio.utils.store_dict_as_TAB_delimited_file(search_results_dict, self.bitscore_log_file, do_not_write_key_column=True)
            self.run.info("Bit score information file: ", self.bitscore_log_file)

        # mark contigs db with hash of modules.db content for version tracking
        self.set_hash_in_contigs_db()

        if anvio.DEBUG:
            self.run.warning("The temp directories, '%s' and '%s' are kept. Please don't forget to clean those up "
                        "later" % (tmp_directory_path, ', '.join(hmmer.tmp_dirs)), header="Debug")
        else:
            self.run.info_single("Cleaning up the temp directory (you can use `--debug` if you would "
                            "like to keep it for testing purposes)", nl_before=1, nl_after=1)
            shutil.rmtree(tmp_directory_path)
            hmmer.clean_tmp_dirs()


