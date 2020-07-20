#!/usr/bin/env python
# -*- coding: utf-8
"""Setup and utilize Interacdome.

The Interacdome is from Mona Singh's group at Princeton. In brief, they attribute empirical
residue-by-residue binding frequencies for Pfam families. Only families with members that have
crystallized structures with ligand co-complexes are used, since otherwise there are no binding
frequencies to attribute. The assumption is that if a residue is close to a bound ligand, it gets a
non-zero weight for being important for binding.

References
==========
- https://interacdome.princeton.edu/
- Systematic domain-based aggregation of protein structures highlights DNA-, RNA- and other
  ligand-binding positions.  Shilpa Nadimapalli Kobren and Mona Singh. Nucleic Acids Research (2019)
  47: 582-593.
"""

import anvio
import anvio.pfam as pfam
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.tables as tables
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.filesnpaths as filesnpaths

from anvio.pfam import Pfam
from anvio.errors import ConfigError
from anvio.parsers import parser_modules
from anvio.drivers.hmmer import HMMer

import os
import shutil
import argparse
import numpy as np
import pandas as pd


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Evan Kiefl"
__email__ = "kiefl.evan@gmail.com"


class InteracdomeSuper(Pfam):
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):

        self.run = run
        self.progress = progress
        self.run.warning("", header='INITIALIZATION', lc='green')

        A = lambda x, t: t(args.__dict__[x]) if x in args.__dict__ else None
        null = lambda x: x
        self.interacdome_data_dir = A('interacdome_data_dir', null) or constants.default_interacdome_data_path
        self.information_content_cutoff = A('information_content_cutoff', null) or 4

        self.hmm_filepath = os.path.join(self.interacdome_data_dir, 'Pfam-A.hmm')

        # Init the InteracDome table
        self.interacdome_table = InteracdomeTableData(kind='representable', interacdome_data_dir=self.interacdome_data_dir)
        self.interacdome_table.load()

        # Init the Pfam baseclass
        args.hmmer_program = 'hmmsearch' # Force use of hmmsearch
        args.pfam_data_dir = self.interacdome_data_dir
        Pfam.__init__(self, args, run=self.run, progress=self.progress)

        # Init contigs database
        args = argparse.Namespace(contigs_db=self.contigs_db_path)
        self.contigs_db = dbops.ContigsSuperclass(args)

        # Init the HMM profile
        self.hmms = pfam.HMMProfile(self.hmm_filepath)

        self.run.warning("Anvi'o will use 'InteracDome' by Kobren and Singh (DOI: 10.1093/nar/gky1224) to attribute binding frequencies. "
                         "If you publish your findings, please do not forget to properly credit their work.", lc='green', header="CITATION")

        # This dictionary is populated and stored as entries in `amino_acid_additional_data`
        self.bind_freq = {}


    def is_database_exists(self):
        """Checks pfam database and interacdome table data exist. Overwrites Pfam.is_database_exists"""

        try:
            Pfam.is_database_exists(self)
        except ConfigError:
            raise ConfigError("It seems you do not have the associated Pfam data required to use "
                              "Interacdome, please run 'anvi-setup-interacdome' to download it. Then "
                              "run this command again.")


    def process(self):
        """Runs Interacdome."""

        tmp_directory_path = filesnpaths.get_temp_directory_path()

        # export AA sequences for genes
        target_files_dict = {'AA:DOMAIN': os.path.join(tmp_directory_path, 'AA_gene_sequences.fa')}
        self.contigs_db.gen_FASTA_file_of_sequences_for_gene_caller_ids(
            output_file_path=target_files_dict['AA:DOMAIN'],
            simple_headers=True,
            rna_alphabet=False,
            report_aa_sequences=True,
        )

        # run hmmer
        hmmer = HMMer(target_files_dict, num_threads_to_use=self.num_threads, program_to_use=self.hmm_program)
        hmm_hits_file = hmmer.run_hmmer(
            source='Interacdome',
            alphabet='AA',
            context='DOMAIN',
            kind=None,
            domain=None,
            num_genes_in_model=len(self.function_catalog),
            hmm=self.hmm_filepath,
            ref=None,
            noise_cutoff_terms='--cut_ga',
            desired_output='standard',
            out_fmt='--domtblout',
        )

        self.hmm_out = parser_modules['search']['hmmer_std_output'](hmm_hits_file, context='interacdome')

        self.run.info('num total domain hits', self.hmm_out.dom_hits.shape[0])
        self.run.info('num unique genes', self.hmm_out.dom_hits['corresponding_gene_call'].unique().shape[0])
        self.run.info('num unique HMMs', self.hmm_out.dom_hits['pfam_id'].unique().shape[0])

        if self.hmm_out.dom_hits.shape[0] == 0:
            self.run.info_single("The HMM search returned no hits :/ So there is nothing to do. Anvi'o "
                                 "will now clean the temporary directories and gracefully quit.",
                                 nl_before=1, nl_after=1)
            shutil.rmtree(tmp_directory_path)
            hmmer.clean_tmp_dirs()
            return

        self.filter_hits()
        self.attribute_binding_frequencies()

        if self.bind_freq.empty:
            self.run.warning("There are 0 HMM hits, so there is nothing to do :( Binding frequencies were not "
                             "added to your database", header="Oh no...")
        else:
            self.store()

        if anvio.DEBUG:
            self.run.warning("The temp directories, '%s' and '%s' are kept. Please don't forget to "
                             "clean those up later" % (tmp_directory_path, ', '.join(hmmer.tmp_dirs)),
                             header="Debug")
        else:
            self.run.info_single("Cleaning up the temp directory (you can use `--debug` if you would "
                                 "like to keep it for testing purposes)", nl_before=1, nl_after=1)
            shutil.rmtree(tmp_directory_path)
            hmmer.clean_tmp_dirs()


    def filter_hits(self):
        """Filter out hits that do not pass criteria

        Removes hits from self.hmm_out.dom_hits and self.self.hmm_out.ali_info
        """

        num_hits = self.hmm_out.dom_hits.shape[0]
        self.progress.new('Filtering HMM hits', progress_total_items=num_hits)

        # lists to remove hits from self.hmm_out.dom_hits
        # each holds 3 length tuples (gene_callers_id, pfam_id, domain_id)
        no_binding_freqs = []
        conserved_match_states_disagree = []
        indices_to_keep = []

        for i, row in self.hmm_out.dom_hits.iterrows():
            pfam_id, gene_callers_id, dom_id = row['pfam_id'], row['corresponding_gene_call'], row['domain']

            if i % 100 == 0:
                self.progress.increment(i)
                self.progress.update('%d/%d processed' % (i, num_hits))

            if pfam_id not in self.interacdome_table.bind_freqs:
                # pfam hit to gene, but has no binding_frequency data
                no_binding_freqs.append((gene_callers_id, pfam_id, dom_id))
                indices_to_keep.append(False)
                continue

            if not self.conserved_positions_match(row):
                # gene sequence did not match the match state consensus values in conserved regions
                conserved_match_states_disagree.append((gene_callers_id, pfam_id, dom_id))
                indices_to_keep.append(False)
                continue

            indices_to_keep.append(True)

        self.progress.increment(len(self.hmm_out.ali_info))
        self.progress.reset()

        self.run.info("filtered hits (no Interacdome data)", len(no_binding_freqs))
        self.run.info("filtered hits (disagrees with conserved match states)", len(conserved_match_states_disagree))

        hits_to_delete = no_binding_freqs + conserved_match_states_disagree
        if not len(hits_to_delete):
            self.progress.end()
            return

        self.progress.update("Removing above hits from data structures")

        self.hmm_out.dom_hits = self.hmm_out.dom_hits.loc[indices_to_keep, :]
        for gene_callers_id, pfam_id, dom_id in hits_to_delete:
            del self.hmm_out.ali_info[gene_callers_id][(pfam_id, dom_id)]

        self.progress.end()


    def conserved_positions_match(self, dom_hit):
        """Returns True if gene sequence matches the consensus of conserved HMM match states

        Conserved is defined as having information content >= self.information_content_cutoff

        Parameters
        ==========
        dom_hit : pandas Series
            A row from self.hmm_out.dom_hits
        """

        pfam_id, gene_callers_id, dom_id = dom_hit['pfam_id'], dom_hit['corresponding_gene_call'], dom_hit['domain']

        # Match states with information content > self.information_content_cutoff (conserved)
        match_states = self.hmms.data[pfam_id]['MATCH_STATES']

        conserved = match_states[match_states['IC'] >= self.information_content_cutoff]
        conserved = conserved[(conserved.index >= dom_hit['hmm_start']) & (conserved.index <= dom_hit['hmm_stop'])]

        if conserved.empty:
            # There are no conserved positions
            return True

        ali = self.hmm_out.ali_info[gene_callers_id][(pfam_id, dom_id)]

        conserved_pos_has_partner = [cons in ali['hmm_positions'].values for cons in conserved.index]
        if not all(conserved_pos_has_partner):
            # Some of the conserved positions had a gap in the sequence, and therefore do not match
            return False

        ali_conserved = ali.loc[ali['hmm_positions'].isin(conserved.index), :]
        conserved_pos_matches = ali_conserved['hmm'].values == ali_conserved['seq'].values
        if not all(conserved_pos_matches):
            # Some of the conserved positions are different in the sequence to the consensus
            return False

        return True


    def attribute_binding_frequencies(self):
        """Populate self.bind_freq dict

        Notes
        =====
        - Must have self.hmm_out (see self.process for example)
        """

        self.bind_freq = {
            'codon_order_in_gene': [],
            'gene_callers_id': [],
            'data_key': [],
            'data_value': [],
        }

        self.run.warning("", header="InteracDome Results", lc='green')

        self.progress.new('Matching binding frequencies to residues', progress_total_items=len(self.hmm_out.ali_info))

        for i, gene_callers_id in enumerate(self.hmm_out.ali_info):

            if i % 100 == 0:
                self.progress.increment(i)
                self.progress.update('%d/%d done' % (i, len(self.hmm_out.ali_info)))

            for pfam_id, dom_id in self.hmm_out.ali_info[gene_callers_id]:

                if pfam_id not in self.interacdome_table.bind_freqs:
                    # pfam hit to gene, but has no binding_frequency data
                    continue

                # The alignment of the sequence with the pfam
                ali = self.hmm_out.ali_info[gene_callers_id][(pfam_id, dom_id)]

                for ligand, freqs in self.interacdome_table.bind_freqs[pfam_id].items():
                    # ali['hmm_positions'] is the 0-indexed positions (match states) in the HMM profile that
                    # aligned to the gene sequence. Hence, the fancy-indexing below subsets freqs to
                    # this subset of match states.
                    attributed_freqs = freqs[ali['hmm_positions'].values]

                    # attributed_freqs and ali['seq_positions'] are now in 1 to 1 correspondence, i.e.
                    # dict(zip(ali['seq_positions'], attributed_freqs)) is a dictionary where keys are
                    # codon_orders and attributed freqs are binding frequencies.

                    # Many of these attributed binding frequenices are 0. We do not bother reporting
                    # these, so we filter them out
                    codon_orders = ali['seq_positions'].values[attributed_freqs > 0]
                    attributed_freqs = attributed_freqs[attributed_freqs > 0]

                    # Populate binding frequencies for this ligand/gene combo
                    self.bind_freq['gene_callers_id'].extend([gene_callers_id] * len(attributed_freqs))
                    self.bind_freq['codon_order_in_gene'].extend(codon_orders)
                    self.bind_freq['data_key'].extend([ligand] * len(attributed_freqs))
                    self.bind_freq['data_value'].extend(list(attributed_freqs))

        self.progress.increment(len(self.hmm_out.ali_info))

        self.progress.update('Casting data as dataframe')
        self.bind_freq = pd.DataFrame(self.bind_freq)

        self.progress.update('Averaging multi-hit positions')
        self.bind_freq = self.bind_freq.groupby(['gene_callers_id', 'data_key', 'codon_order_in_gene']).mean().reset_index()
        self.bind_freq['data_type'] = 'float'
        self.bind_freq['data_group'] = 'Interacdome'

        self.progress.end()

        self.run.info("Num genes with at least 1 binding score", self.bind_freq['gene_callers_id'].nunique())
        self.run.info("Num positions with >0.00 binding scores", (self.bind_freq['data_value'] >= 0).sum())
        self.run.info("Num positions with >0.25 binding scores", (self.bind_freq['data_value'] >= 0.25).sum())
        self.run.info("Num positions with >0.50 binding scores", (self.bind_freq['data_value'] >= 0.50).sum())
        self.run.info("Num positions with >0.75 binding scores", (self.bind_freq['data_value'] >= 0.75).sum())
        self.run.info("Num positions with =1.00 binding score", (self.bind_freq['data_value'] == 1.00).sum())


    def store(self):
        self.progress.new("Storing binding frequencies in contigs database")
        self.progress.update("...")

        contigs_db = dbops.ContigsDatabase(self.contigs_db_path, run=self.run, progress=self.progress)

        # This is horrible. But we hijacked the miscdata framework to create
        # amino_acid_addtional_data table in the contigs DB, so we must live with the consequences
        self.bind_freq['item_name'] = self.bind_freq['gene_callers_id'].astype(str) + ':' + self.bind_freq['codon_order_in_gene'].astype(str)

        contigs_db.db.insert_rows_from_dataframe(
            tables.amino_acid_additional_data_table_name,
            self.bind_freq[tables.amino_acid_additional_data_table_structure]
        )

        contigs_db.disconnect()
        self.progress.end()


class InteracdomeTableData(object):
    """Manages query and access of the interacdome tabular data sets"""

    def __init__(self, kind='representable', interacdome_data_dir=None):
        self.kind = kind
        self.interacdome_data_dir = interacdome_data_dir

        if not self.interacdome_data_dir:
            self.interacdome_data_dir = constants.default_interacdome_data_path

        self.files = {
            'representable': os.path.join(self.interacdome_data_dir, 'representable_interactions.txt'),
            'confident': os.path.join(self.interacdome_data_dir, 'confident_interactions.txt'),
        }

        if kind not in self.files:
            raise ConfigError("Unknown kind '%s' of Interacdome data. Known kinds: %s" \
                              % (kind, ','.join(list(self.files.keys()))))

        self.filepath = self.files[self.kind]
        filesnpaths.is_file_exists(self.filepath)

        self.df = None
        self.bind_freqs = {}


    def load(self):
        # Load the table as a df
        df = pd.read_csv(os.path.join(self.interacdome_data_dir, self.filepath), sep='\t', comment='#')
        df.index = df['pfam_id'].str.split('_', n=1, expand=True)[0]
        self.df = df

        # Populate bind_freqs arrays
        for pfam_id, row in self.df.iterrows():
            if pfam_id not in self.bind_freqs:
                self.bind_freqs[pfam_id] = {}

            self.bind_freqs[pfam_id][row['ligand_type']] = np.array(row['binding_frequencies'].split(',')).astype(float)


class InteracdomeSetup(object):
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        """Setup a Pfam database for anvi'o

        Parameters
        ==========
        args : argparse.Namespace
            See `bin/anvi-setup-interacdome` for available arguments
        """

        self.run = run
        self.progress = progress
        self.interacdome_data_dir = args.interacdome_data_dir

        self.pfam_setup = None

        self.interacdome_files = {
            'representable_interactions.txt': 'https://interacdome.princeton.edu/session/344807c477180230032a2a3807ba47c3/download/downloadBP?w=',
            'confident_interactions.txt': 'https://interacdome.princeton.edu/session/344807c477180230032a2a3807ba47c3/download/downloadConfidentBP?w=',
        }

        if self.interacdome_data_dir and args.reset:
            raise ConfigError("You are attempting to run Interacdome setup on a non-default data directory (%s) using the --reset flag. "
                              "To avoid automatically deleting a directory that may be important to you, anvi'o refuses to reset "
                              "directories that have been specified with --interacdome-data-dir. If you really want to get rid of this "
                              "directory and regenerate it with Interacdome data inside, then please remove the directory yourself using "
                              "a command like `rm -r %s`. We are sorry to make you go through this extra trouble, but it really is "
                              "the safest way to handle things." % (self.interacdome_data_dir, self.interacdome_data_dir))

        if not self.interacdome_data_dir:
            self.interacdome_data_dir = constants.default_interacdome_data_path

        self.run.warning('', header='Setting up Interacdome', lc='yellow')
        self.run.info('Data directory', self.interacdome_data_dir)
        self.run.info('Reset contents', args.reset)

        filesnpaths.is_output_dir_writable(os.path.dirname(os.path.abspath(self.interacdome_data_dir)))

        if not args.reset and not anvio.DEBUG:
            self.is_database_exists()

        filesnpaths.gen_output_directory(self.interacdome_data_dir, delete_if_exists=args.reset)


    def is_database_exists(self):
        """Raise ConfigError if database exists

        Currently, this primitively decides if the database exists by looking at whether Pfam-A.hmm
        or Pfam-A.hmm.gz exists.
        """

        if (os.path.exists(os.path.join(self.interacdome_data_dir, 'Pfam-A.hmm') or 
            os.path.exists(os.path.join(self.interacdome_data_dir, 'Pfam-A.hmm.gz')))):
            raise ConfigError("It seems you already have the Interacdome data downloaded in '%s', please "
                              "use --reset flag if you want to re-download it." % self.interacdome_data_dir)


    def setup(self):
        """The main method of this class. Sets up the Interacdome data directory for usage"""

        self.run.warning('', header='Downloading Interacdome tables', lc='yellow')
        self.download_interacdome_files()

        self.run.warning('', header='Downloading associated Pfam HMM profiles', lc='yellow')
        self.download_pfam()

        self.run.warning('', header='Filtering Pfam HMM profiles', lc='yellow')
        self.filter_pfam()


    def download_interacdome_files(self):
        """Download the confident and representable non-redundant Interacdome datasets

        These datasets can be found at the interacdome webpage: https://interacdome.princeton.edu/
        """

        for path, url in self.interacdome_files.items():
            utils.download_file(
                url,
                os.path.join(self.interacdome_data_dir, path),
                check_certificate=False,
                progress=self.progress,
                run=self.run
            )


    def download_pfam(self):
        """Setup the pfam data subset used by interacdome

        Currently, interacdome only works for pfam version 31.0, so that is the version downloaded here.
        """

        pfam_args = argparse.Namespace(
            pfam_data_dir=self.interacdome_data_dir,
            pfam_version='31.0',
            reset=False,
        )

        self.pfam_setup = pfam.PfamSetup(pfam_args)
        self.pfam_setup.get_remote_version()
        self.pfam_setup.download(hmmpress_files=False)


    def load_interacdome(self, kind='representable'):
        """Loads the representable interacdome dataset as pandas df"""

        data = InteracdomeTableData(kind=kind, interacdome_data_dir=self.interacdome_data_dir)
        return data.get_as_dataframe()


    def get_interacdome_pfam_accessions(self):
        """Get the representable interacdome Pfam accessions"""

        return set(self.load_interacdome(kind='representable').index.tolist())


    def filter_pfam(self):
        """Filter Pfam data according to whether the ACC is in the Interacdome dataset"""

        interacdome_pfam_accessions = self.get_interacdome_pfam_accessions()

        hmm_profiles = pfam.HMMProfile(os.path.join(self.interacdome_data_dir, 'Pfam-A.hmm'))
        hmm_profiles.filter(by='ACC', subset=interacdome_pfam_accessions)
        hmm_profiles.write(filepath=None) # overwrites

        # hmmpresses the new .hmm
        self.pfam_setup.hmmpress_files()

        # We also filter out the Pfam-A.clans.tsv file, since it is used as a catalog
        clans_file = os.path.join(self.interacdome_data_dir, 'Pfam-A.clans.tsv')
        clans = pd.read_csv(clans_file, sep='\t', header=None)
        clans[clans[0].isin(interacdome_pfam_accessions)].to_csv(clans_file, sep='\t', index=False, header=False)
