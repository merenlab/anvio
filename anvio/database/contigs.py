
class ContigsDatabase(object):
    """To create an empty contigs database and/or access one."""
    def __init__(self, db_path, run=run, progress=progress, quiet=True, create_new=False, meta_values=None):
        self.db = None
        self.db_path = db_path

        self.run = run
        self.progress = progress
        self.create_new = create_new
        self.quiet = quiet

        self.meta = {}
        self.meta_for_new_db = meta_values

        if self.create_new and self.meta_for_new_db:
            self.create()
        else:
            raise ConfigError("You want to create a new contigs database but did not provide any 'meta_values'.")
        
        self.init_database()
        self.init_tables()


    def init_database(self):
        if os.path.exists(self.db_path):
            is_contigs_db(self.db_path)
            self.db = db.DB(self.db_path, anvio.__contigs__version__)
            meta_table = self.db.get_table_as_dict('self')
            self.meta = dict([(k, meta_table[k]['value']) for k in meta_table])

            for key in ['split_length', 'kmer_size', 'total_length', 'num_splits', 'num_contigs', 'genes_are_called', 'splits_consider_gene_calls']:
                try:
                    self.meta[key] = int(self.meta[key])
                except:
                    pass

            self.meta['gene_function_sources'] = [s.strip() for s in self.meta['gene_function_sources'].split(',')] if self.meta['gene_function_sources'] else None

            if 'creation_date' not in self.meta:
                raise ConfigError("The contigs database ('%s') seems to be corrupted :/ This happens if the process that\
                                    that generates the database ends prematurely. Most probably, you will need to generate\
                                    the contigs database from scratch. Sorry!" % (self.db_path))
            else:
                self.a_meta['creation_date'] = utils.get_time_to_date(self.a_meta['creation_date'])


            if not self.create_new:
                # Actually we can print these information even though create_new flag is on, 
                # but since newly created database is empty it will not give any useful information
                self.run.info('Contigs database', 'An existing database, %s, has been initiated.' % self.db_path, quiet=self.quiet)
                self.run.info('Number of contigs', self.meta['num_contigs'], quiet=self.quiet)
                self.run.info('Number of splits', self.meta['num_splits'], quiet=self.quiet)
                self.run.info('Total number of nucleotides', self.meta['total_length'], quiet=self.quiet)
                self.run.info('Split length', self.meta['split_length'], quiet=self.quiet)
        else:
            # WHY we dont raise an exception here?
            self.db = None


    def init_tables(self):
        self.splits_basic_info = {}
        self.splits_taxonomy_dict = {}
        self.split_sequences = {}
        self.contigs_basic_info = {}
        self.contig_sequences = {}

        self.genes_in_contigs_dict = {}
        self.gene_lengths = {}
        self.contig_name_to_genes = {}
        self.genes_in_splits = {} # keys of this dict are NOT gene caller ids. they are ids for each entry.
        self.genes_in_splits_summary_dict = {}
        self.genes_in_splits_summary_headers = []
        self.split_name_to_genes_in_splits_entry_ids = {} # for fast access to all self.genes_in_splits entries for a given split
        self.gene_callers_id_to_split_name_dict = {} # for fast access to a split name that contains a given gene callers id

        self.auxiliary_contigs_data_available = False
        self.auxiliary_contigs_data_path = None
        self.nt_positions_info = None

        self.gene_function_call_sources = []
        self.gene_function_calls_dict = {}
        self.gene_function_calls_initiated = False

        self.hmm_sources_info = {}
        self.hmm_searches_dict = {}   # <--- upon initiation, this dict only keeps hmm hits for non-singlecopy
        self.hmm_searches_header = [] #      gene searches... single-copy gene info is accessed through completeness.py

        self.singlecopy_gene_hmm_sources = set([])
        self.non_singlecopy_gene_hmm_sources = set([])

        self.progress.new('Loading the contigs DB')
        self.progress.update('Reading contigs basic info')
        self.contigs_basic_info = self.db.get_table_as_dict(t.contigs_info_table_name, string_the_key=True)

        self.progress.update('Reading splits basic info')
        self.splits_basic_info = self.db.get_table_as_dict(t.splits_info_table_name)

        self.progress.update('Reading genes in contigs table')
        self.genes_in_contigs_dict = self.db.get_table_as_dict(t.genes_in_contigs_table_name)

        self.progress.update('Populating gene lengths dict')
        self.gene_lengths = dict([(g, (self.genes_in_contigs_dict[g]['stop'] - self.genes_in_contigs_dict[g]['start'])) for g in self.genes_in_contigs_dict])

        self.progress.update('Populating contig name to gene IDs dict')
        for contig_name in self.contigs_basic_info:
            self.contig_name_to_genes[contig_name] = set([])
        for gene_unique_id in self.genes_in_contigs_dict:
            e = self.genes_in_contigs_dict[gene_unique_id]
            self.contig_name_to_genes[e['contig']].add((gene_unique_id, e['start'], e['stop']), )

        self.progress.update('Reading genes in splits table')
        self.genes_in_splits = self.db.get_table_as_dict(t.genes_in_splits_table_name)

        self.progress.update('Reading genes in splits summary table')
        self.genes_in_splits_summary_dict = self.db.get_table_as_dict(t.genes_in_splits_summary_table_name)
        self.genes_in_splits_summary_headers = self.db.get_table_structure(t.genes_in_splits_summary_table_name)

        self.progress.update('Identifying HMM searches for single-copy genes and others')
        self.hmm_sources_info = self.db.get_table_as_dict(t.hmm_hits_info_table_name)
        for hmm_source in self.hmm_sources_info:
            self.hmm_sources_info[hmm_source]['genes'] = sorted([g.strip() for g in self.hmm_sources_info[hmm_source]['genes'].split(',')])

        self.singlecopy_gene_hmm_sources = set([s for s in list(self.hmm_sources_info.keys()) if self.hmm_sources_info[s]['search_type'] == 'singlecopy'])
        self.non_singlecopy_gene_hmm_sources = set([s for s in list(self.hmm_sources_info.keys()) if self.hmm_sources_info[s]['search_type'] != 'singlecopy'])

        self.progress.update('Generating "split name" to "gene entry ids" mapping dict')
        for entry_id in self.genes_in_splits:
            split_name = self.genes_in_splits[entry_id]['split']
            if split_name in self.split_name_to_genes_in_splits_entry_ids:
                self.split_name_to_genes_in_splits_entry_ids[split_name].add(entry_id)
            else:
                self.split_name_to_genes_in_splits_entry_ids[split_name] = set([entry_id])

        for split_name in self.splits_basic_info:
            if split_name not in self.split_name_to_genes_in_splits_entry_ids:
                self.split_name_to_genes_in_splits_entry_ids[split_name] = set([])

        self.progress.update('Generating "gene caller id" to "split name" mapping dict')
        for entry in list(self.genes_in_splits.values()):
            self.gene_callers_id_to_split_name_dict[entry['gene_callers_id']] = entry['split']

        contigs_db.disconnect()

        self.progress.update('Accessing the auxiliary data file')
        auxiliary_contigs_data_path = get_auxiliary_data_path_for_contigs_db(self.contigs_db_path)
        if os.path.exists(auxiliary_contigs_data_path):
            self.auxiliary_contigs_data_available = True
            self.auxiliary_contigs_data_path = auxiliary_contigs_data_path
            self.nt_positions_info = auxiliarydataops.AuxiliaryDataForNtPositions(auxiliary_contigs_data_path, self.a_meta['contigs_db_hash'])
            self.progress.end()
        else:
            self.progress.end()
            self.run.warning("Auxiliary contigs data ('%s') is not available. Some operations related to\
                              variability analyses will not be available." % auxiliary_contigs_data_path)

        if self.auxiliary_contigs_data_available:
            self.run.info('Auxiliary Data', 'Found: %s (v. %s)' % (auxiliary_contigs_data_path, anvio.__auxiliary_data_version__))

        self.run.info('Contigs DB', 'Initialized: %s (v. %s)' % (self.contigs_db_path, anvio.__contigs__version__))


    def get_date(self):
        return time.time()


    def get_hash(self):
        return '%08x' % random.randrange(16**8)


    def create(self):
        """Creates an empty contigs database on disk, and sets `self.db` to access to it.

        At some point self.db.disconnect() must be called to complete the creation of the new db."""

        is_db_ok_to_create(self.db_path, 'contigs')

        self.db = db.DB(self.db_path, anvio.__contigs__version__, new_database=True)

        # creating empty default tables
        self.db.create_table(t.hmm_hits_table_name, t.hmm_hits_table_structure, t.hmm_hits_table_types)
        self.db.create_table(t.hmm_hits_info_table_name, t.hmm_hits_info_table_structure, t.hmm_hits_info_table_types)
        self.db.create_table(t.hmm_hits_splits_table_name, t.hmm_hits_splits_table_structure, t.hmm_hits_splits_table_types)
        self.db.create_table(t.collections_info_table_name, t.collections_info_table_structure, t.collections_info_table_types)
        self.db.create_table(t.collections_bins_info_table_name, t.collections_bins_info_table_structure, t.collections_bins_info_table_types)
        self.db.create_table(t.collections_contigs_table_name, t.collections_contigs_table_structure, t.collections_contigs_table_types)
        self.db.create_table(t.collections_splits_table_name, t.collections_splits_table_structure, t.collections_splits_table_types)
        self.db.create_table(t.genes_in_contigs_table_name, t.genes_in_contigs_table_structure, t.genes_in_contigs_table_types)
        self.db.create_table(t.genes_in_splits_table_name, t.genes_in_splits_table_structure, t.genes_in_splits_table_types)
        self.db.create_table(t.splits_taxonomy_table_name, t.splits_taxonomy_table_structure, t.splits_taxonomy_table_types)
        self.db.create_table(t.taxon_names_table_name, t.taxon_names_table_structure, t.taxon_names_table_types)
        self.db.create_table(t.genes_taxonomy_table_name, t.genes_taxonomy_table_structure, t.genes_taxonomy_table_types)
        self.db.create_table(t.contig_sequences_table_name, t.contig_sequences_table_structure, t.contig_sequences_table_types)
        self.db.create_table(t.gene_function_calls_table_name, t.gene_function_calls_table_structure, t.gene_function_calls_table_types)
        self.db.create_table(t.gene_protein_sequences_table_name, t.gene_protein_sequences_table_structure, t.gene_protein_sequences_table_types)
        self.db.create_table(t.genes_in_splits_summary_table_name, t.genes_in_splits_summary_table_structure, t.genes_in_splits_summary_table_types)
        self.db.create_table(t.splits_info_table_name, t.splits_info_table_structure, t.splits_info_table_types)
        self.db.create_table(t.contigs_info_table_name, t.contigs_info_table_structure, t.contigs_info_table_types)

        # know thyself
        self.db.set_meta_value('db_type', 'contigs')
        self.db.set_meta_value('split_length', self.meta_for_new_db['split_length'])
        self.db.set_meta_value('project_name', self.meta_for_new_db['project_name'])
        self.db.set_meta_value('description', self.meta_for_new_db['description'])
        self.db.set_meta_value('kmer_size', self.meta_for_new_db['kmer_size'])


    def disconnect(self):
        self.db.disconnect()