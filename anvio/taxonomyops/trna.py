#!/usr/bin/env python3
# -*- coding: utf-8
"""
Classes to use local databases for tRNA taxonomy to affiliate tRNA seqeunces in anvi'o
databases with taxon names.
"""

import os
import shutil
import hashlib

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.drivers.blast import BLAST
from anvio.dbops import ContigsDatabase

from anvio.taxonomyops import AccessionIdToTaxonomy
from anvio.taxonomyops import PopulateContigsDatabaseWithTaxonomy

__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


run_quiet = terminal.Run(log_file_path=None, verbose=False)
progress_quiet = terminal.Progress(verbose=False)
pp = terminal.pretty_print

HASH = lambda d: str(hashlib.sha224(''.join([str(d[level]) for level in constants.levels_of_taxonomy]).encode('utf-8')).hexdigest()[0:8])

class TRNATaxonomyContext(AccessionIdToTaxonomy):
    """The purpose of this base class is ot define file paths and constants for trna taxonomy ops."""

    def __init__(self, trna_taxonomy_data_dir=None, scgs_taxonomy_remote_database_url=None, run=terminal.Run(), progress=terminal.Progress()):
        self.run = run
        self.progress = progress

        # know thyself.
        self.focus = "trnas"

        # hard-coded GTDB variables. poor design, but I don't think we are going do need an
        # alternative to GTDB.
        self.target_database = "GTDB"
        self.target_database_URL = "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/"

        # some variables from anvi'o constants
        self.hmm_source_for_trna_genes = constants.default_hmm_source_for_trna_genes
        self.default_trna_taxonomy_data_dir = constants.default_trna_taxonomy_data_dir
        self.default_anticodons_for_taxonomy = constants.default_anticodons_for_taxonomy
        self.levels_of_taxonomy = constants.levels_of_taxonomy

        # these are all the user accessible paths. defaults will serve well for all applications,
        # but these can be used for debugging.
        self.tRNA_taxonomy_data_dir = (os.path.abspath(trna_taxonomy_data_dir) if trna_taxonomy_data_dir else None) or (os.path.join(self.default_trna_taxonomy_data_dir, self.target_database))
        self.accession_to_taxonomy_file_path = os.path.join(self.tRNA_taxonomy_data_dir, 'ACCESSION_TO_TAXONOMY.txt.gz')
        self.database_version_file_path = os.path.join(self.tRNA_taxonomy_data_dir, 'VERSION')
        self.search_databases_dir_path = os.path.join(self.tRNA_taxonomy_data_dir, 'ANTICODON_SEARCH_DATABASES')
        self.target_database_URL = scgs_taxonomy_remote_database_url or self.target_database_URL

        # some dictionaries for convenience. we set them up here, but the proper place to sanity check
        # them may be somewhere else. for instance, when this class is inheritded by SetupLocalTRNATaxonomyData
        # the paths will not point to an actual file, but when it is inherited by PopulateContigsDatabaseWithTRNATaxonomy,
        # they better point to actual files.
        self.anticodons = dict([(anticodon, {'db': os.path.join(self.search_databases_dir_path, anticodon)}) for anticodon in self.default_anticodons_for_taxonomy])

        self.accession_to_taxonomy_dict = {}

        # set version for ctx, so we know what version of the databases are on disk
        if os.path.exists(self.database_version_file_path):
            self.trna_taxonomy_database_version = open(self.database_version_file_path).readline().strip()
        else:
            self.trna_taxonomy_database_version = None

        # populate `self.accession_to_taxonomy_dict`
        AccessionIdToTaxonomy.__init__(self)


# here we create an instance for the module. the idea is to overwrite it if
# it is necessary to overwrite some of the defaults
ctx = TRNATaxonomyContext()


class SanityCheck(object):
    def __init__(self):
        if self.skip_sanity_check:
            self.run.warning("We are skipping all sanity checks :( Dangerous stuff is happening.")
        else:
            self.sanity_check()


    def sanity_check(self):
        if sorted(constants.default_anticodons_for_taxonomy) != sorted(self.ctx.default_anticodons_for_taxonomy):
            raise ConfigError("Oh no. The anticodons designated to be used for all tRNA taxonomy tasks in the constants.py "
                              "are not the same names described in locally known HMMs to remote FASTA files "
                              "conversion table definedd in SetupLocalTRNATaxonomyData module. If this makes zero "
                              "sense to you please ask a developer.")

        if not self.ctx.tRNA_taxonomy_data_dir:
            raise ConfigError("`SetupLocalTRNATaxonomyData` class is upset because it was inherited without "
                              "a directory for tRNA taxonomy data to be stored :( This variable can't be None.")

        if self.user_taxonomic_level and self.user_taxonomic_level not in constants.levels_of_taxonomy:
            raise ConfigError("The taxonomic level %s is not a level anvi'o knows about. Here is the list of "
                              "taxonomic levels anvi'o recognizes: %s" % (', '.join(constants.levels_of_taxonomy)))

        # sanity checks specific to classes start below
        if self.__class__.__name__ in ['SetupLocalTRNATaxonomyData']:
            pass

        if self.__class__.__name__ in ['SetupLocalTRNATaxonomyData', 'PopulateContigsDatabaseWithTRNATaxonomy']:
            if self.user_taxonomic_level:
                raise ConfigError("There is no need to set a taxonomic level while working with the class SetupLocalTRNATaxonomyData "
                                  "or PopulateContigsDatabaseWithTRNATaxonomy. Something fishy is going on :/")

        if self.__class__.__name__ in ['PopulateContigsDatabaseWithTRNATaxonomy', 'TRNATaxonomyEstimatorSingle', 'TRNATaxonomyEstimatorMulti']:
            if not os.path.exists(self.ctx.tRNA_taxonomy_data_dir):
                raise ConfigError("Anvi'o could not find the data directory for the tRNA taxonomy setup. If you have "
                                  "a non-default location for your tRNA taxonomy databases, please use the parameter "
                                  "`--trna-taxonomy-data-dir` parameter). Anvi'o tried to find your files here: '%s'" % (self.ctx.tRNA_taxonomy_data_dir))

            if not os.path.exists(self.ctx.accession_to_taxonomy_file_path) or not os.path.exists(self.ctx.database_version_file_path):
                raise ConfigError("While your tRNA taxonomy data dir seems to be in place, it is missing at least one critical "
                                  "file. This is not someting you can add or remove as this file is distributed with anvi'o "
                                  "releases :( Please get in touch with a developer, or fix it if you are one.")

            if not os.path.exists(self.ctx.accession_to_taxonomy_file_path):
                raise ConfigError("While your tRNA taxonomy data dir seems to be in place, it is missing at least one critical "
                                  "file (in this case, the file to resolve accession IDs to taxon names). This is not someting "
                                  "you can add or remove as this file is distributed with anvi'o releases :( Please get in touch "
                                  "with a developer, or fix it if you are one.")

            ###########################################################
            # PopulateContigsDatabaseWithTRNATaxonomy
            ###########################################################
            if self.__class__.__name__ in ['PopulateContigsDatabaseWithTRNATaxonomy']:
                for prefix in ['.nhr', '.nin', '.nsq']:
                    missing_anticodon_databases = [anticodon for anticodon in self.ctx.anticodons if not os.path.exists(self.ctx.anticodons[anticodon]['db'] + '.nhr')]
                    if len(missing_anticodon_databases):
                        raise ConfigError("OK. It is very likley that if you run `anvi-setup-scg-taxonomy` first you will be golden. "
                                          "Because even though anvi'o found the directory for taxonomy headquarters, "
                                          "your setup seems to be missing %d of %d databases required for everything to work "
                                          "with the current genes configuration of this class (sources say this is a record, FYI)." % \
                                                    (len(missing_anticodon_databases), len(self.ctx.anticodons)))

class TRNATaxonomyArgs(object):
    def __init__(self, args, format_args_for_single_estimator=False):
        """A base class to fill in common arguments for tRNA Taxonomy classes.

        The purpose of this class is to reduce the complexity of setting member variables
        for various classes in this module so we can get away with a single multi-talented
        sanity check base class without any complaints regarding missing member variables.

        Parameters
        ==========
        format_args_for_single_estimator: bool
            This is a special case where an args instance is generated to be passed to the
            single estimator from within multi estimator. More specifically, the multi estimator
            class is nothing but one that iterates through all metagenomes
            given to it using the single estimator class. So it needs to create instances of
            single estimators, and collect results at upstream. The problem is, if a single
            estimtor is initiated with the args of a multi estimator, the sanity check will
            go haywire. This flag nullifies most common offenders.
        """

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.output_file_path = A('output_file')
        self.output_file_prefix = A('output_file_prefix')
        self.just_do_it = A('just_do_it')
        self.simplify_taxonomy_information = A('simplify_taxonomy_information')
        self.metagenome_mode = True if A('metagenome_mode') else False
        self.anticodon_for_metagenome_mode = A('anticodon_for_metagenome_mode')
        self.compute_anticodon_coverages = A('compute_anticodon_coverages')
        self.report_anticodon_frequencies_path = A('report_anticodon_frequencies')
        self.metagenomes = A('metagenomes')
        self.user_taxonomic_level = A('taxonomic_level')
        self.matrix_format = A('matrix_format')
        self.raw_output = A('raw_output')

        if format_args_for_single_estimator:
            # so you're here to get an args instance to fool a single estimator class.
            # very cute. we shall make that happen.
            self.metagenomes = None
            self.output_file_path = None
            self.output_file_prefix = None
            self.matrix_format = None
            self.raw_output = None

        self.skip_sanity_check = A('skip_sanity_check')


class TRNATaxonomyEstimatorMulti(TRNATaxonomyArgs, SanityCheck):
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress(), skip_init=False):
        """Iterate through metagenome descriptions using TRNATaxonomyEstimatorSingle"""

        self.args = args
        self.run = run
        self.progress = progress

        # update your self args
        TRNATaxonomyArgs.__init__(self, self.args)

        # set your context
        self.ctx = ctx

        # intiate sanity check
        SanityCheck.__init__(self)


class TRNATaxonomyEstimatorSingle(TRNATaxonomyArgs, SanityCheck, TaxonomyEstimatorSingle):
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress(), skip_init=False):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.contigs_db_path = A('contigs_db')
        self.profile_db_path = A('profile_db')
        self.collection_name = A('collection_name')
        self.update_profile_db_with_taxonomy = A('update_profile_db_with_taxonomy')
        self.bin_id = A('bin_id')

        TRNATaxonomyArgs.__init__(self, self.args)

        self.ctx = ctx

        SanityCheck.__init__(self)

        TaxonomyEstimatorSingle.__init__(self, skip_init=skip_init)


class SetupLocalTRNATaxonomyData(TRNATaxonomyArgs, SanityCheck):
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        # update your self args
        TRNATaxonomyArgs.__init__(self, self.args)

        # user accessible variables
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.redo_databases = A("redo_databases")
        self.num_threads = A('num_threads')

        self.ctx = ctx

        SanityCheck.__init__(self)


    def setup(self):
        """This function downloads all GTDB files necessary to setup the SCG databases anvi'o will rely upon.

           In addition to downloading the original files, the setup will make sure everything, including the
           DIAMOND search databases are in place.
        """

        self.run.warning("Please remember that the data anvi'o uses for SCG taxonomy is a "
                         "courtesy of The Genome Taxonomy Database (GTDB), an initiative to establish a "
                         "standardised microbial taxonomy based on genome phylogeny, primarly funded by "
                         "tax payers in Australia. Please don't forget to cite the original work, "
                         "doi:10.1038/nbt.4229 by Parks et al to explicitly mention the source of databases "
                         "anvi'o relies upon to estimate genome level taxonomy. If you are not sure how it "
                         "should look like in your methods sections, anvi'o developers will be happy to "
                         "help you if you can't find any published example to get inspiration.", lc = 'yellow')

        # let's try to figure out what do we want to do here. download files or just create/re-create
        # databases.
        an_anticodon_fasta, an_anticodon_database = list(self.ctx.anticodons.values())[0]['db'] + '.gz', list(self.ctx.anticodons.values())[0]['db'] + '.nhr'

        if os.path.exists(an_anticodon_fasta) and os.path.exists(an_anticodon_database) and self.redo_databases:
            self.run.warning("Anvi'o is removing all the previous databases so it can regenerate them from their "
                             "ashes.")

            db_paths = [v['db'] for v in self.ctx.anticodons.values()]
            for db_path in db_paths:
                os.remove(db_path) if os.path.exists(db_path) else None

        elif os.path.exists(an_anticodon_fasta) and os.path.exists(an_anticodon_database) and not self.redo_databases:
            raise ConfigError("It seems you have both the FASTA files and the search databases for anvi'o SCG taxonomy "
                              "in place. If you want to regenerate the databases for some reason, please use the flag "
                              "`--redo-databases`")

        elif os.path.exists(an_anticodon_fasta) and not os.path.exists(an_anticodon_database):
            self.run.warning("Anvi'o found your FASTA files in place, but not the databases. Now it will generate all "
                             "the search databases using the existing FASTA files.")

        self.create_search_databases()

        if not anvio.DEBUG:
            self.clean_up()


    def download_and_format_files(self):
        """There is no downloading files for this module.

        Essentially this is a manual operation. Please see the README.md file.
        """

        pass


    def create_search_databases(self):
        """Creates all the search databases"""

        self.progress.new("Creating search databases")
        self.progress.update("Removing any database that still exists in the output directory...")
        for prefix in ['.nhr', '.nin', '.nsq']:
            [os.remove(database_path) for database_path in [s['db'] + prefix for s in self.ctx.anticodons.values()] if os.path.exists(database_path)]

        # compresssing and decompressing FASTA files changes their hash and make them look like
        # modified in git. to avoid that, we will do the database generation in a temporary directory.
        temp_dir = filesnpaths.get_temp_directory_path()

        self.progress.update("Copying FASTA files to %s ..." % (temp_dir))
        # the following line basically returns a dictionary that shows the new path
        # of the FASTA file under temp_dir for a given anticodon .. apologies for the
        # incomprehensible list comprehension
        new_paths = dict([(os.path.basename(fasta_path), shutil.copy((fasta_path + '.gz'), os.path.join(temp_dir, os.path.basename(fasta_path) + '.gz'))) for fasta_path in [s['db'] for s in self.ctx.anticodons.values()]])

        missing_FASTA_files = [anticodon for anticodon in self.ctx.anticodons if not os.path.exists(new_paths[anticodon])]
        if len(missing_FASTA_files):
            raise ConfigError("Weird news :( Anvi'o is missing some FASTA files that were supposed to be somewhere. Since this "
                              "can't be your fault, it is not easy to advice what could be the solution to this. If you are not "
                              "an anvi'o programmer working on this problem this very moment, please get in touch with one.")

        self.progress.update("Decompressing FASTA files in %s" % (temp_dir))
        new_paths = dict([(anticodon, utils.gzip_decompress_file(new_paths[anticodon], keep_original=False)) for anticodon in new_paths])

        for anticodon in self.ctx.anticodons:
            self.progress.update("Working on %s in %d threads" % (anticodon, self.num_threads))

            FASTA_file_path_for_anticodon = new_paths[anticodon]

            # create a BLAST search database for `FASTA_file_path_for_anticodon`
            blast = BLAST(query_fasta=FASTA_file_path_for_anticodon, run=run_quiet, progress=progress_quiet, num_threads=self.num_threads)
            blast.log_file_path = os.path.join(os.path.dirname(FASTA_file_path_for_anticodon), '%s.log' % anticodon)
            blast.makedb(dbtype='nucl')

            for prefix in ['.nhr', '.nin', '.nsq']:
                if not os.path.exists(FASTA_file_path_for_anticodon + prefix):
                    raise ConfigError("Something went wrong and BLAST did not create the database file it was supposed to "
                                      "for %s :(" % anticodon)
                else:
                    shutil.move(FASTA_file_path_for_anticodon + prefix, os.path.dirname(self.ctx.anticodons[anticodon]['db']))

        shutil.rmtree(temp_dir)

        self.progress.end()
        self.run.info_single("Every FASTA is now turned into a fancy search database. It means you are now allowed to run "
                             "`anvi-run-trna-taxonomy` on anvi'o contigs databases. This workflow is very new, and there are "
                             "caveats to it just like every other computational approach you use to make sense of complex 'omics "
                             "data. To better understand those caveats you should read our online documentation a bit. If you see "
                             "things that concerns you, please let anvi'o developers know. They love bad news. If you get good "
                             "results from this workflow, thank to those who contributed to the GTDB.", nl_after=1, mc="green")


    def clean_up(self):
        pass


class PopulateContigsDatabaseWithTRNATaxonomy(TRNATaxonomyArgs, SanityCheck, PopulateContigsDatabaseWithTaxonomy):
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        # update your self args
        TRNATaxonomyArgs.__init__(self, self.args)

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.write_buffer_size = int(A('write_buffer_size') if A('write_buffer_size') is not None else 1000)
        self.contigs_db_path = A('contigs_db')
        self.num_parallel_processes = int(A('num_parallel_processes')) if A('num_parallel_processes') else 1
        self.num_threads = int(A('num_threads')) if A('num_threads') else 1

        self.ctx = ctx

        self.max_target_seqs = 20
        self.evalue = float(A('e_value')) if A('e_value') else 1e-05
        self.min_pct_id = float(A('min_percent_identity')) if A('min_percent_identity') else 90

        SanityCheck.__init__(self)

        self.focus = 'trnas'
        PopulateContigsDatabaseWithTaxonomy.__init__(self, self.args)
