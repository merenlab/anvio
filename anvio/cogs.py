# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Making sense of COGs.
"""

import os
import gzip
import glob
import shutil
import hashlib

import anvio
import anvio.utils as utils
import anvio.dbops as dbops
import anvio.dictio as dictio
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.drivers.blast import BLAST
from anvio.drivers.diamond import Diamond
from anvio.tables.genefunctions import TableForGeneFunctions

# just to make sure things don't break too far when they do:
COG_DATA_VERSION='2'


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print
P = terminal.pluralize

J = lambda x, y: os.path.join(x, y)

# if you add a new database version, please do not forget to add its reference here
COG_REFERENCES = {'COG20': {
                            'ref_text': 'Galperin et al. 2021',
                            'doi_link': 'https://doi.org/10.1093/nar/gkaa1018',
                           },
                  'COG14': {
                            'ref_text': 'Galperin et al. 2015',
                            'doi_link': 'https://doi.org/10.1093/nar/gku1223',
                           },
                  'arCOG14': {
                            'ref_text': 'Makarova, Wolf, and Koonin 2015',
                            'doi_link': 'https://doi.org/10.3390%2Flife5010818',
                             },
}


class Args():
    pass


class COGs:
    """A class to run COGs"""
    def __init__(self, args=Args(), run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.num_threads = A('num_threads')
        self.contigs_db_path = A('contigs_db')
        self.search_with = A('search_with') or 'diamond'
        self.temp_dir_path = A('temporary_dir_path')

        self.log_file_path = None

        self.default_search_method = 'diamond'
        self.search_methods_factory = {'diamond': self.search_with_diamond,
                                       'blastp': self.search_with_ncbi_blast}
        self.available_search_methods = [p for p in self.search_methods_factory.keys() if utils.is_program_exists(p, dont_raise=True)]

        if not len(self.available_search_methods):
            raise ConfigError("None of the search methods this class could use, which include '%s', seem to be "
                              "available on your system :/" % (', '.join(list(self.search_methods_factory.keys()))))

        if self.default_search_method not in self.available_search_methods:
            self.default_search_method = self.available_search_methods[0]

        if len(args.__dict__):
            self.initialize(args)


    def initialize(self, args):
        self.hits = None # the search function will take care of this one.

        # get an instance of the setup class
        self.COG_setup = COGsSetup(args, run=self.run, progress=self.progress)

        self.COG_data_dir = self.COG_setup.COG_data_dir
        self.COG_base_dir = self.COG_setup.COG_base_dir
        self.COG_version = self.COG_setup.COG_version
        self.available_db_search_program_targets = self.COG_setup.get_formatted_db_paths()
        self.essential_files = self.COG_setup.get_essential_file_paths()

        # this is an additional check after the 2020 COG release.
        if os.path.exists(os.path.join(self.COG_base_dir, 'COG.txt')):
            # This means that this installation still has the old COG data. Let's try to solve it:
            self.move_old_COG_data_to_its_new_location()

        # Check whether there is data for the version requested
        data_available_for_cog_versions = [os.path.basename(d) for d in glob.glob(os.path.join(self.COG_base_dir, '*COG*'))]
        if not len(data_available_for_cog_versions):
            raise ConfigError("You don't seem to have any COG data setup in your COG data directory. Please first run "
                              "`anvi-setup-ncbi-cogs` to take care of that.")
        elif self.COG_version not in data_available_for_cog_versions:
            if self.COG_version == "COG20" and len(data_available_for_cog_versions) == 1 and data_available_for_cog_versions[0] == "COG14":
                raise ConfigError("OK. There is a problem, but it is easy to solve: On this system you have the data for the 2014 "
                                  "release of the NCBI's COGs. But the NCBI recently made a new release, and `anvi-run-ncbi-cogs` "
                                  "sets the COG data version to that release. If you want to run the older data, you can rerun your "
                                  "command by adding the parameter `--cog-version COG14`. If you want to take this opportunity to "
                                  "setup the new COG data (which you totally should), please go ahead and run `anvi-setup-ncbi-cogs`.")
            else:
                raise ConfigError(f"You requested to run the NCBI COGs using the version '{self.COG_version}', but you don't have "
                                  f"the data to be able to run that version of COGs setup on this system yet. Here is the "
                                  f"{P('version', len(data_available_for_cog_versions), pfs='only')} available "
                                  f"for COGs on this system: {', '.join(data_available_for_cog_versions)}. You can "
                                  f"either specify which version you wish to use through the parameter `--cog-version`, or you can "
                                  f"run the following command to setup the data files for the version you wish to use: "
                                  f"`anvi-setup-ncbi-cogs --cog-version {self.COG_version}`, and re-run the same command that brought "
                                  f"you here.")

        # citation, if possible
        if self.COG_version in COG_REFERENCES:
            ref_str = f"{COG_REFERENCES[self.COG_version]['ref_text']} ({COG_REFERENCES[self.COG_version]['doi_link']})"
            self.run.warning(f"Anvi'o will annotate your database with the {self.COG_version} version of NCBI COGs. "
                             f"We recommend citing the following reference for this database when you publish your findings : "
                             f"{ref_str}", lc='green', header="CITATION")


    def move_old_COG_data_to_its_new_location(self):
        try:
            filesnpaths.is_output_dir_writable(self.COG_base_dir)
        except:
            raise ConfigError(f"Please read this carefully: The NCBI has made a new release of COGs. To make room for that "
                              f"while maintaining the old COG data from 2014 version, anvi'o needs to move some files around. "
                              f"While anvi'o can do it automatically, your user does not seem to have permission to do that. "
                              f"One alternative is to ask your system administrator to run this program on your behalf. It will "
                              f"solve everything. OR you can ask them to do exactly these steps: (1) go to the directory "
                              f"'{self.COG_base_dir}', (2) create a new directory called `COG14`, and (3) move everything in "
                              f"'{self.COG_base_dir}' (WHICH INCLUDES the files: CATEGORIES.txt, COG.txt, DB_BLAST/ "
                              f"DB_DIAMOND/, MISSING_COG_IDs.cPickle, PID-TO-CID.cPickle, and RAW_DATA_FROM_NCBI/ as well as the "
                              f"hidden file .VERSION) into the new `COG14` directory. Then you will be golden.")

        # we have the write permission, so let's do this.
        tmp_dir = filesnpaths.get_temp_directory_path(just_the_path=True)
        self.run.warning(f"This is a bit important: The NCBI has made a new release of COGs. To make room for that "
                         f"while maintaining the old COG data from 2014 version, anvi'o needs to move some files around. "
                         f"It seems you have the necessary permissions to write into anvi'o misc data directory, so anvi'o "
                         f"will now attempt to do it automatically by first moving things to a temporary directory "
                         f"('{tmp_dir}') and then moving them back into their new target location. If you have not been "
                         f"having an exceptionally bad day, this should go smoothly. But if you see an error below, anvi'o is "
                         f"very sorry for breaking itself on your system :( In which case please find us on our Discord channel "
                         f"and we will try to help you to sort things out.")
        self.progress.new("Moving files around")
        shutil.move(self.COG_base_dir, tmp_dir)
        os.makedirs(self.COG_base_dir)
        shutil.move(tmp_dir, os.path.join(self.COG_base_dir, 'COG14'))

        self.run.info_single("Congratulations! Anvi'o managed to migrate your old data into its new location without breaking "
                             "things. We are all very proud here but let's never do this again.", mc='green', nl_after=1)


    def process(self, aa_sequences_file_path=None):
        if self.search_with not in self.available_search_methods:
            raise ConfigError("Let us start by making it clear that we probably like '%s' as much as you do, but it doesn't "
                              "seem to be available on your system OR recognized by the COGs class since anvi'o couldn't "
                              "find it among the available search methods. You probably need to try something else :/" \
                                                                                                    % self.search_with)

        if self.search_with not in self.available_db_search_program_targets:
            raise ConfigError("Anvi'o understands that you want to use '%s' to search for COGs, however, there is no "
                              "database formatted under the COGs data directory for that program :/ You may need to "
                              "re-run the COGs setup (anvi-setup-ncbi-cogs), UNLESS, you set up your COG data directory "
                              "somewhere else than what anvi'o attempts to use at the moment ('%s'). If that is the case, "
                              "this may be the best time to point the right directory using the --cog-data-dir parameter, "
                              "or the environmental variable 'ANVIO_COG_DATA_DIR'." % (self.search_with, self.COG_data_dir))

        if not aa_sequences_file_path and not self.contigs_db_path:
            raise ConfigError("You either need to provide an anvi'o contigs database path, or a FASTA file for AA "
                              "sequences")

        if aa_sequences_file_path and self.contigs_db_path:
            raise ConfigError("You can't provide both an AA sequences file and a contigs database. Choose one!")

        if self.contigs_db_path:
            utils.is_contigs_db(self.contigs_db_path)

        if not self.temp_dir_path:
            self.temp_dir_path = filesnpaths.get_temp_directory_path()
            self.remove_temp_dir_path = True
        else:
            filesnpaths.is_file_exists(self.temp_dir_path)
            filesnpaths.is_output_dir_writable(self.temp_dir_path)

            self.run.warning("Because you set the temporary directory path by hand, anvi'o will not remove its content "
                             "when it is done. But she certainly hopes that you will clean those files later.")

            self.remove_temp_dir_path = False

        self.run.info('COG data directory', self.COG_data_dir)
        self.run.info('Searching with', self.search_with)
        self.run.info('Directory to store temporary files', self.temp_dir_path)
        self.run.info('Directory will be removed after the run', self.remove_temp_dir_path)


        if not aa_sequences_file_path:
            aa_sequences_file_path = J(self.temp_dir_path, 'aa_sequences.fa')
            dbops.ContigsSuperclass(self.args, r=terminal.Run(verbose=False)).get_sequences_for_gene_callers_ids(output_file_path=aa_sequences_file_path,
                                                                                  report_aa_sequences=True,
                                                                                  simple_headers=True)

        # do the search
        search_results_tabular = self.search_methods_factory[self.search_with](aa_sequences_file_path)

        # convert the output to a hits dict
        if self.COG_version == 'COG14' or self.COG_version == 'arCOG14':
            self.hits = utils.get_BLAST_tabular_output_as_dict(search_results_tabular, target_id_parser_func=lambda x: x.split('|')[1])
        elif self.COG_version == 'COG20':
            self.hits = utils.get_BLAST_tabular_output_as_dict(search_results_tabular)
        else:
            raise ConfigError("You need to edit all the if/else statements with COG version checks to ensure proper "
                              "parsing of a new generation of COG files.")

        # store hits into the contigs database
        self.store_hits_into_contigs_db()

        if self.remove_temp_dir_path:
            shutil.rmtree(self.temp_dir_path)


    def store_hits_into_contigs_db(self):
        if not self.hits:
            self.run.warning("COGs class has no hits to process. Returning empty handed, but still adding COGs as "
                             "functional sources.")
            gene_function_calls_table = TableForGeneFunctions(self.contigs_db_path, self.run, self.progress)

            if self.COG_version == 'COG14':
                gene_function_calls_table.add_empty_sources_to_functional_sources({'COG14_FUNCTION', 'COG14_CATEGORY'})
            elif self.COG_version == 'COG20':
                gene_function_calls_table.add_empty_sources_to_functional_sources({'COG20_FUNCTION', 'COG20_CATEGORY', 'COG20_PATHWAY'})
            else:
                raise ConfigError("You need to edit all the if/else statements with COG version checks to ensure proper "
                                  "parsing of a new generation of COG files.")
            return

        cogs_data = COGsData(self.args)
        cogs_data.init_p_id_to_cog_id_dict()

        functions_dict = {}
        self.__entry_id = 0


        def add_entry(gene_callers_id, source, accession, function, e_value):
            functions_dict[self.__entry_id] = {'gene_callers_id': int(gene_callers_id),
                                               'source': source,
                                               'accession': accession,
                                               'function': function,
                                               'e_value': float(e_value)}
            self.__entry_id += 1

        # let's keep track of hits that match to missing COGs
        hits_for_missing_cogs = 0
        hits_for_missing_ncbi_protein_ids = 0

        missing_cogs_found = set([])
        missing_ncbi_protein_ids_found = set([])

        for gene_callers_id in self.hits:
            ncbi_protein_id = self.hits[gene_callers_id]['hit']

            in_proteins_FASTA_not_in_cogs_CSV = []
            if ncbi_protein_id not in cogs_data.p_id_to_cog_id:
                in_proteins_FASTA_not_in_cogs_CSV.append((ncbi_protein_id, gene_callers_id),)
                missing_ncbi_protein_ids_found.add(ncbi_protein_id)
                hits_for_missing_ncbi_protein_ids += 1
                continue

            COG_ids = cogs_data.p_id_to_cog_id[ncbi_protein_id]

            annotations = []
            categories = []
            pathways = []
            category_descriptions = []
            for COG_id in COG_ids:
                # is missing?
                if COG_id in cogs_data.missing_cogs:
                    missing_cogs_found.add(COG_id)
                    hits_for_missing_cogs += 1
                    continue

                # resolve categories
                for category in cogs_data.cogs[COG_id]['categories']:
                    categories.append(category)
                    category_descriptions.append(cogs_data.categories[category])

                # append annotation
                annotations.append(cogs_data.cogs[COG_id]['annotation'])

                if self.COG_version == 'COG14' or self.COG_version == 'arCOG14':
                    pass
                elif self.COG_version == 'COG20':
                    if cogs_data.cogs[COG_id]['pathway']:
                        pathways.append(cogs_data.cogs[COG_id]['pathway'])
                else:
                    raise ConfigError("You need to edit all the if/else statements with COG version checks to ensure proper "
                                      "parsing of a new generation of COG files.")

            # all these shitty heuristics... If there are multiple COG ids or categories, separate them from each other by '!!!' so parsing
            # them later is possible. Am I embarrassed? Yes. Is there a better way of doing this efficiently? Absolutely. What time is it?
            # 9pm. Where am I? In the lab. Is it OK for me to let this slip away if it means for me to go home sooner? Yes, probably. Am I
            # gonna remember this crap in the code for the next two months at random times in the shower and feel bad about myself? Fuck yes.
            if self.COG_version == 'COG14' or self.COG_version == 'arCOG14':
                add_entry(gene_callers_id, f'{self.COG_version}_FUNCTION', '!!!'.join(COG_ids), '!!!'.join(annotations), self.hits[gene_callers_id]['evalue'])
                add_entry(gene_callers_id, f'{self.COG_version}_CATEGORY', '!!!'.join(categories), '!!!'.join(category_descriptions), 0.0)
            elif self.COG_version == 'COG20':
                add_entry(gene_callers_id, f'{self.COG_version}_FUNCTION', '!!!'.join(COG_ids), '!!!'.join(annotations), self.hits[gene_callers_id]['evalue'])
                add_entry(gene_callers_id, f'{self.COG_version}_CATEGORY', '!!!'.join(categories), '!!!'.join(category_descriptions), 0.0)
                if len(pathways):
                    add_entry(gene_callers_id, f'{self.COG_version}_PATHWAY', '!!!'.join(COG_ids), '!!!'.join(pathways), self.hits[gene_callers_id]['evalue'])
            else:
                raise ConfigError("You need to edit all the if/else statements with COG version checks to ensure proper "
                                  "parsing of a new generation of COG files.")

        # store hits in contigs db.
        gene_function_calls_table = TableForGeneFunctions(self.contigs_db_path, self.run, self.progress)
        gene_function_calls_table.create(functions_dict)

        if len(missing_cogs_found):
            self.run.warning('Although your COGs are successfully added to the database, there were some COG IDs your genes hit '
                             'were among the ones that were not described in the raw data. Here is the list of %d COG IDs that '
                             'were hit %d times: %s.' % (len(missing_cogs_found), hits_for_missing_cogs, ', '.join(missing_cogs_found)))

        if len(missing_ncbi_protein_ids_found):
            self.run.warning("Well. Your COGs were successfully added to the database, but there were some garbage anvi'o brushed "
                             "off under the rug. There were %d genes in your database that hit %d protein IDs in NCBIs COGs database, "
                             "but since NCBI did not release what COGs they correspond to in the database they made available (that "
                             "helps us to resolve protein IDs to COG ids), we could not annotate those genes with functions. Anvi'o "
                             "apologizes on behalf of all computer scientists for half-done stuff we often force biologists to deal "
                             "with. If you want to do some Googling, these were the offending protein IDs: '%s'." % \
                                        (hits_for_missing_ncbi_protein_ids, len(missing_ncbi_protein_ids_found), ', '.join([str(s) for s in missing_ncbi_protein_ids_found])))

        if len(in_proteins_FASTA_not_in_cogs_CSV):
            # so some of the hits represented in the FASTA file from the NCBI were not put in the
            # CSV file from NCBI to associate them with COGs
            report_output_file_path = filesnpaths.get_temp_file_path()
            report_output = open(report_output_file_path, 'w')
            report_output.write('anvio_gene_callers_id\tNCBI_protein_id\n')

            for protein_id, gene_callers_id in in_proteins_FASTA_not_in_cogs_CSV:
                report_output.write('%s\t%s\n' % (gene_callers_id, protein_id))

            report_output.close()

            self.run.warning("This is important. %s hits for your genes that appeared in the proteins FASTA file from the NCBI had protein "
                             "IDs that were not described in the CSV file from the NCBI that associates each protein ID with a COG function. "
                             "That's OK if you don't care. But if you would like to take a look, anvi'o stored a report "
                             "file for you at %s" \
                        % (len(in_proteins_FASTA_not_in_cogs_CSV), report_output_file_path))


    def search_with_diamond(self, aa_sequences_file_path):
        diamond = Diamond(aa_sequences_file_path, run=self.run, progress=self.progress, num_threads=self.num_threads)

        diamond.target_fasta = self.available_db_search_program_targets['diamond']
        self.run.log_file_path = self.log_file_path or J(self.temp_dir_path, 'log.txt')
        diamond.search_output_path = J(self.temp_dir_path, 'diamond-search-results')
        diamond.tabular_output_path = J(self.temp_dir_path, 'diamond-search-results.txt')

        diamond.max_target_seqs = 1

        diamond.blastp()
        diamond.view()

        return diamond.tabular_output_path


    def search_with_ncbi_blast(self, aa_sequences_file_path):
        blast = BLAST(aa_sequences_file_path, run=self.run, progress=self.progress, num_threads=self.num_threads)

        blast.target_fasta = self.available_db_search_program_targets['blastp']
        self.run.log_file_path = self.log_file_path or J(self.temp_dir_path, 'log.txt')
        blast.search_output_path = J(self.temp_dir_path, 'blast-search-results.txt')
        blast.max_target_seqs = 1

        blast.blast()

        return blast.search_output_path


class COGsData:
    """A class to make sense of COG ids and categories"""
    def __init__(self, args=Args(), cog_data_dir = None, run=run, progress=progress, panic_on_failure_to_init=False):
        self.run = run
        self.progress = progress

        if cog_data_dir:
            args.cog_data_dir = cog_data_dir

        self.setup = COGsSetup(args, run=self.run, progress=self.progress)
        self.essential_files = self.setup.get_essential_file_paths()
        self.COG_version = self.setup.COG_version

        self.cogs = None
        self.p_id_to_cog_id = None
        self.missing_cogs = None
        self.categories = None
        self.initialized = False

        if self.essential_files:
            self.init()
        elif panic_on_failure_to_init:
            raise ConfigError("It seems you don't have your COG data set up on this system. Whatever you were "
                              "trying to do is not going to continue being done :( Did you setup your COGs? If "
                              "not, you can take a look at the program `anvi-setup-ncbi-cogs`. Maybe you did "
                              "setup into another directory than the default destination? If that is the case "
                              "maybe you can use the `--cog-data-dir` parameter if it is applicable? No? None "
                              "of these work? Well. Anvi'o hates it as much as you do when things come to this.")


    def init(self):
        self.progress.new('Initializing COGs Data')
        self.progress.update('Reading COG functions ...')

        if self.COG_version == 'COG14' or self.COG_version == 'arCOG14':
            self.cogs = utils.get_TAB_delimited_file_as_dictionary(self.essential_files['COG.txt'], no_header=True, column_names=['COG', 'categories', 'annotation'])
        elif self.COG_version == 'COG20':
            self.cogs = utils.get_TAB_delimited_file_as_dictionary(self.essential_files['COG.txt'], no_header=True, column_names=['COG', 'categories', 'annotation', 'pathway'])
        else:
            raise ConfigError("You need to edit all the if/else statements with COG version checks to ensure proper "
                              "parsing of a new generation of COG files.")


        self.progress.update('Reading COG categories ...')
        self.categories = utils.get_TAB_delimited_file_as_dictionary(self.essential_files['CATEGORIES.txt'], no_header=True, column_names=['category', 'description'])

        self.progress.update('Reading missing COG IDs ...')
        self.missing_cogs = dictio.read_serialized_object(self.essential_files['MISSING_COG_IDs.cPickle'])

        self.progress.end()

        for cog in self.cogs:
            self.cogs[cog]['categories'] = [c.strip() for c in self.cogs[cog]['categories'].split(',')]

        for cat in self.categories:
            self.categories[cat] = self.categories[cat]['description']

        self.initialized = True


    def init_p_id_to_cog_id_dict(self):
        self.progress.new('Initializing COGs Data')
        self.progress.update('Reading NCBI Protein ID to COG id converter ...')

        self.p_id_to_cog_id = dictio.read_serialized_object(self.essential_files['PID-TO-CID.cPickle'])

        self.progress.end()


class COGsSetup:
    """A class to download and setup the COG data from NCBI."""
    def __init__(self, args=Args(), cog_data_dir=None, run=run, progress=progress):
        self.run = run
        self.progress = progress

        # description of COG primary files per version
        self.cog_files = {'COG14':
                             {'cog2003-2014.csv': {
                                  'url': 'ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/cog2003-2014.csv',
                                  'func': self.format_p_id_to_cog_id_cPickle,
                                  'type': 'essential',
                                  'formatted_file_name': 'PID-TO-CID.cPickle'},
                              'cognames2003-2014.tab': {
                                  'url': 'ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/cognames2003-2014.tab',
                                  'func': self.format_cog_names,
                                  'type': 'essential',
                                  'formatted_file_name': 'COG.txt'},
                              'fun2003-2014.tab': {
                                  'url': 'ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/fun2003-2014.tab',
                                  'func': self.format_categories,
                                  'type': 'essential',
                                  'formatted_file_name': 'CATEGORIES.txt'},
                              'prot2003-2014.fa.gz': {
                                  'url': 'ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/prot2003-2014.fa.gz',
                                  'func': self.format_protein_db,
                                  'type': 'database',
                                  'formatted_file_name': 'IGNORE_THIS_AND_SEE_THE_FUNCTION'}
                             },
                        'arCOG14':
                             {'ar14.arCOG.csv': {
                                  'url': 'ftp://ftp.ncbi.nih.gov/pub/wolf/COGs/arCOG/ar14.arCOG.csv',
                                  'func': self.format_p_id_to_cog_id_cPickle,
                                  'type': 'essential',
                                  'formatted_file_name': 'PID-TO-CID.cPickle'},
                              'ar14.arCOGdef.tab': {
                                  'url': 'ftp://ftp.ncbi.nih.gov/pub/wolf/COGs/arCOG/ar14.arCOGdef.tab',
                                  'func': self.format_cog_names,
                                  'type': 'essential',
                                  'formatted_file_name': 'COG.txt'},
                              'funclass.tab': {
                                  'url': 'ftp://ftp.ncbi.nih.gov/pub/wolf/COGs/arCOG/funclass.tab',
                                  'func': self.format_categories,
                                  'type': 'essential',
                                  'formatted_file_name': 'CATEGORIES.txt'},
                              'ar14.fa.gz': {
                                  'url': 'ftp://ftp.ncbi.nih.gov/pub/wolf/COGs/arCOG/ar14.fa.gz',
                                  'func': self.format_protein_db,
                                  'type': 'database',
                                  'formatted_file_name': 'IGNORE_THIS_AND_SEE_THE_FUNCTION'},
                             },
                        'COG20':
                             {'cog-20.cog.csv': {
                                  'url': 'ftp://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.cog.csv',
                                  'func': self.format_p_id_to_cog_id_cPickle,
                                  'type': 'essential',
                                  'formatted_file_name': 'PID-TO-CID.cPickle'},
                              'cog-20.def.tab': {
                                  'url': 'ftp://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.def.tab',
                                  'func': self.format_cog_names,
                                  'type': 'essential',
                                  'formatted_file_name': 'COG.txt'},
                              'fun-20.tab': {
                                  'url': 'ftp://ftp.ncbi.nih.gov/pub/COG/COG2020/data/fun-20.tab',
                                  'func': self.format_categories,
                                  'type': 'essential',
                                  'formatted_file_name': 'CATEGORIES.txt'},
                              'checksum.md5.txt': {  # No func as it is called by the setup_raw_data function
                                   'url': 'ftp://ftp.ncbi.nih.gov//pub/COG/COG2020/data/checksums.md5.txt',
                                   'type': 'non-essential',
                                   'formatted_file_name': 'CHECKSUMS.txt'},
                              'cog-20.fa.gz': {
                                  'url': 'ftp://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.fa.gz',
                                  'func': self.format_protein_db,
                                  'type': 'database',
                                  'formatted_file_name': 'IGNORE_THIS_AND_SEE_THE_FUNCTION'},
                             },
                        }

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.num_threads = A('num_threads') or 1
        self.just_do_it = A('just_do_it')
        self.reset = A('reset')
        self.COG_data_source = 'unknown'
        self.COG_version = A('cog_version') or 'COG20'

        if self.COG_version not in self.cog_files:
            raise ConfigError(f"The COG versions known to anvi'o do not include '{self.COG_version}' :/ This is "
                              f"what we know of: {', '.join(self.cog_files.keys())}. This is one of those things "
                              f"that should have never happened. We salute you.")

        if cog_data_dir:
            self.COG_base_dir = cog_data_dir
            self.COG_data_source = 'The function call.'
        elif A('cog_data_dir'):
            self.COG_base_dir = A('cog_data_dir')
            self.COG_data_source = 'The command line parameter.'
        elif 'ANVIO_COG_DATA_DIR' in os.environ:
            self.COG_base_dir = os.environ['ANVIO_COG_DATA_DIR']
            self.COG_data_source = 'The environmental variable.'
        else:
            self.COG_base_dir = J(os.path.dirname(anvio.__file__), 'data/misc/COG')
            self.COG_data_source = "The anvi'o default."

        self.COG_base_dir = os.path.abspath(os.path.expanduser(self.COG_base_dir))
        self.COG_data_dir = os.path.join(self.COG_base_dir, self.COG_version)

        self.run.info('COG version', self.COG_version, mc='green')
        self.run.info('COG data source', self.COG_data_source)
        self.run.info('COG base directory', self.COG_base_dir)

        self.COG_data_dir_version = J(self.COG_data_dir, '.VERSION')
        self.raw_NCBI_files_dir = J(self.COG_data_dir, 'RAW_DATA_FROM_NCBI')

        self.files = self.cog_files[self.COG_version]

        self.cogs_found_in_proteins_fasta = set([])
        self.cogs_found_in_cog_names_file = set([])

        # citation, if possible
        if self.COG_version in COG_REFERENCES:
            ref_str = f"{COG_REFERENCES[self.COG_version]['ref_text']} ({COG_REFERENCES[self.COG_version]['doi_link']})"
            self.run.warning(f"Anvi'o will set up the {self.COG_version} version of NCBI COGs. "
                             f"We recommend citing the following reference for this database when you publish your findings : "
                             f"{ref_str}", lc='green', header="CITATION")


    def get_formatted_db_paths(self):
        formatted_db_paths = {}

        diamond_db_path = J(self.COG_data_dir, 'DB_DIAMOND')
        if os.path.exists(diamond_db_path):
            formatted_db_paths['diamond'] = J(diamond_db_path, 'COG')

        blast_db_path = J(self.COG_data_dir, 'DB_BLAST')
        if os.path.exists(blast_db_path):
            formatted_db_paths['blastp'] = J(blast_db_path, 'COG/COG.fa')

        return formatted_db_paths


    def get_essential_file_paths(self):
        if not os.path.exists(self.COG_data_dir):
            # the COG_data_dir is not there
            return None

        essential_files = {}
        for v in list(self.files.values()):
            if v['type'] == 'essential':
                essential_files[v['formatted_file_name']] = J(self.COG_data_dir, v['formatted_file_name'])

        # add the missing COG IDs file into the list:
        essential_files['MISSING_COG_IDs.cPickle'] = J(self.COG_data_dir, 'MISSING_COG_IDs.cPickle')

        for file_name in essential_files:
            if not os.path.exists(essential_files[file_name]):
                raise ConfigError("At least one essential formatted file that is necesary for COG operations is not where it should "
                                   "be ('%s'). You should run COG setup, with the flag `--reset` if necessary, to make sure things "
                                   "are in order." % essential_files[file_name])

        return essential_files


    def create(self):
        if not os.path.exists(self.COG_data_dir):
            try:
                os.makedirs(self.COG_data_dir)
                open(self.COG_data_dir_version, 'w').write(COG_DATA_VERSION)
            except Exception as e:
                raise ConfigError(f"So the COG data directory is not there, and anvi'o wants to create one. But it didn't "
                                  f"go that well. It could be due to permissions (which may require you to run this with sudo "
                                  f"or may need to ask your sys admin to do it for you since this is a one time operation), or "
                                  f"it could be due to something totally irrelevant. Here is the error message: {e}.")

        filesnpaths.is_output_dir_writable(self.COG_data_dir)

        if self.reset:
            run.warning('This program will remove everything in the COG data directory, then download and reformat '
                        'everything from scratch.')

            # OK. reset the crap out of it.
            shutil.rmtree(self.COG_data_dir)
            os.mkdir(self.COG_data_dir)
            open(self.COG_data_dir_version, 'w').write(COG_DATA_VERSION)
        else:
            run.warning("This program will first check whether you have all the raw files, and then will attempt to "
                        "regenerate everything that is necessary from them.")

        if not os.path.exists(self.COG_data_dir_version) or open(self.COG_data_dir_version).read().strip() != COG_DATA_VERSION:
            raise ConfigError("The version of your COG data directory is different than what anvi'o hoping to see. "
                               "It seems you need to (re)run anvi'o script to download and format COG data from NCBI.")

        # get raw files
        self.get_raw_data()

        # format raw files
        self.setup_raw_data()

        # identify missing COGs
        self.generate_missing_cog_ids_file()


    def generate_missing_cog_ids_file(self):
        missing_cog_ids = self.cogs_found_in_proteins_fasta.difference(self.cogs_found_in_cog_names_file)

        if len(missing_cog_ids):
            self.run.warning("%d of %d COG IDs that appear in the list of orthology domains file (which links protein IDs "
                             "to COG names), are missing from the COG names file (which links COG IDs to function names and "
                             "categories). Because clearly even the files that are distributed together should not be expected to "
                             "be fully compatible. Anvi'o thanks everyone for their contributions." % \
                                                        (len(missing_cog_ids), len(self.cogs_found_in_proteins_fasta)))

        dictio.write_serialized_object(missing_cog_ids, J(self.COG_data_dir, 'MISSING_COG_IDs.cPickle'))


    def format_p_id_to_cog_id_cPickle(self, input_file_path, output_file_path):
        num_lines_in_file = filesnpaths.get_num_lines_in_file(input_file_path)

        def raise_error(line_num, line_content, fields, e):
            raise ConfigError(f"Bad news :( While parsing a COG input file, anvi'o encountered an error (which said: [{e}]) "
                              f"while processing the line {line_counter} in your file. Where the fields in that file looked "
                              f"looked like this: {fields}. Sadly, this has been a long-standing and very annoying issue that "
                              f"anvi'o developers were unable to reproduce. But we recently learned that the issue is likely due "
                              f"to your internet speed (https://github.com/merenlab/anvio/issues/1738). Slower connections lead "
                              f"to broken connections with the NCBI servers, and leave you with an unfinished file :/ The only "
                              f"working solution so far is to try again with a faster internet connection.")

        progress.new('Formatting protein ids to COG ids file', progress_total_items=num_lines_in_file)

        p_id_to_cog_id = {}
        p_id_without_cog_id = set([])

        line_counter = 0
        for line in open(input_file_path, 'rU').readlines():
            line_counter += 1

            if line_counter % 500 == 0:
                self.progress.increment(line_counter)
                progress.update(f"{line_counter * 100 / num_lines_in_file:.2f}%")

            fields = line.strip('\n').split(',')

            # the arCOG14 release comes with a broken file. Starting on line 388147, all following lines
            # include only the first 6 fields. What we need is in the 7th field (the COG ID), so we are
            # basically forced to ignore all COGs after that. Sucks.
            if self.COG_version == 'arCOG14' and line_counter >= 388147:
                if line_counter == 388147: # once we hit this line, print the warning
                    num_lines_to_end = num_lines_in_file - line_counter + 1
                    self.run.warning(f"There is a problem with the {input_file_path} file downloaded from NBCI. "
                                    f"Basically, starting from line {line_counter}, the arCOG ID number is not provided, which "
                                    f"means that we cannot match those protein sequences to their COG IDs. The only solution we "
                                    f"have at the moment is to skip the {num_lines_to_end} protein IDs that are affected by this "
                                    f"issue. Sorry.")

                # let's save the protein IDs that don't have an associated arCOG ID
                p_id_without_cog_id.add(fields[2].replace('.', '_'))

                continue

            # `p_id` should look just like the FASTA ids, and its location has changed between
            # 2014 release and 2020 release.
            if self.COG_version == 'COG14':
                try:
                    p_id = fields[0]
                    COG = fields[6]
                except Exception as e:
                    raise_error(line_counter, line, fields, e)
            elif self.COG_version == 'COG20' or self.COG_version == 'arCOG14':
                try:
                    p_id = fields[2].replace('.', '_')
                    COG = fields[6]
                except Exception as e:
                    raise_error(line_counter, line, fields, e)
            else:
                raise ConfigError("You need to edit all the if/else statements with COG version checks to ensure proper "
                                  "parsing of a new generation of COG files.")

            self.cogs_found_in_proteins_fasta.add(COG)

            if p_id in p_id_to_cog_id:
                if COG not in p_id_to_cog_id[p_id]:
                    p_id_to_cog_id[p_id].append(COG)
            else:
                p_id_to_cog_id[p_id] = [COG]

        progress.update("Serializing the data dictionary for future use (a.k.a, very pro stuff).")
        dictio.write_serialized_object(p_id_to_cog_id, output_file_path)

        progress.end()

        if p_id_without_cog_id:
                self.run.warning(f"There were {len(p_id_without_cog_id)} protein IDs without an associated COG ID. "
                                 f"This may cause issues later, so please keep this warning in mind. Here are a few examples "
                                 f"of the affected protein IDs: {', '.join(list(p_id_without_cog_id)[:5])}")


    def format_cog_names(self, input_file_path, output_file_path):
        progress.new('Formatting COG names file')
        progress.update('...')
        output = open(output_file_path, 'w')

        try:
            lines = open(input_file_path).readlines()
        except UnicodeDecodeError:
            lines = open(input_file_path, encoding='ISO-8859-1').readlines()

        for line in lines:
            if line.startswith('#'):
                continue

            if self.COG_version == 'COG14':
                # example line from 2014:
                #
                # COG0059 EH      Ketol-acid reductoisomerase
                COG, function, name = line.strip('\n').split('\t')
                name = ''.join([i if ord(i) < 128 else '' for i in name])
                output.write('\t'.join([COG, ', '.join(list(function)), name]) + '\n')
            elif self.COG_version == 'arCOG14':
                # example line from arCOG 2014:
                #
                # arCOG00024	C	GlpK	Glycerol kinase	COG00554	pfam00370,pfam02782	cd07786	TIGR01311
                COG, category, gene, function, superclust, pfam_profiles, cdd_profiles, tigrfam_profiles = line.strip('\n').split('\t')

                function = ''.join([i if ord(i) < 128 else '' for i in function])
                function = function if not gene else f"{function} ({gene})"

                output.write('\t'.join([COG, ', '.join(list(category)), function]) + '\n')

            elif self.COG_version == 'COG20':
                # example line from 2020:
                #
                # COG0059	EH	Ketol-acid reductoisomerase	IlvC	Isoleucine, leucine, valine biosynthesis		1NP3
                COG, category, function, nn, pathway, pubmed_id, PDB_id = line.strip('\n').split('\t')

                function = ''.join([i if ord(i) < 128 else '' for i in function])
                function = function if not nn else f"{function} ({nn})"
                function = function if not PDB_id else f"{function} (PDB:{PDB_id})"
                function = function if not pubmed_id else f"{function} (PUBMED:{pubmed_id})"

                output.write('\t'.join([COG, ', '.join(list(category)), function, pathway]) + '\n')

            else:
                raise ConfigError("You need to edit all the if/else statements with COG version checks to ensure proper "
                                  "parsing of a new generation of COG files.")

            self.cogs_found_in_cog_names_file.add(COG)
            progress.end()


    def format_categories(self, input_file_path, output_file_path):
        progress.new('Formatting COG categories file')
        progress.update('...')

        output = open(output_file_path, 'w')
        for line in open(input_file_path, 'rU').readlines():
            if line.startswith('#'):
                continue

            if self.COG_version == 'COG14':
                category, description = line.strip('\n').split('\t')
            elif self.COG_version == 'COG20' or self.COG_version == 'arCOG14':
                category, _, description = line.strip('\n').split('\t')
            else:
                raise ConfigError("You need to edit all the if/else statements with COG version checks to ensure proper "
                                  "parsing of a new generation of COG files.")

            # get rid of non-ascii chars:
            description = ''.join([i if ord(i) < 128 else '' for i in description])

            output.write('\t'.join([category, description]) + '\n')

        progress.end()


    def format_protein_db(self, input_file_path, output_file_path):
        progress.new('Formatting raw files')
        progress.update('Decompressing protein sequences')

        # poor man's uncompress
        temp_fasta_path = filesnpaths.get_temp_file_path()
        try:
            with open(temp_fasta_path, 'wb') as f_out, gzip.open(input_file_path, 'rb') as f_in:
                f_out.write(f_in.read())
        except Exception as e:
            progress.end()
            raise ConfigError(f"Something went wrong while decompressing the downloaded file :/ It is likely that "
                              f"the download failed and only part of the file was downloaded. If you would like to "
                              f"try again, please run the setup command with the flag `--reset`. Here is what the "
                              f"downstream library said: '{e}'.")

        progress.end()

        if utils.is_program_exists('diamond', dont_raise=True):
            output_dir = J(self.COG_data_dir, 'DB_DIAMOND')
            if os.path.exists(output_dir):
                shutil.rmtree(output_dir)

            os.mkdir(output_dir)

            output_db_path = J(output_dir, 'COG')
            log_file_path = J(output_dir, 'log.txt')

            self.run.info('Diamond log', log_file_path)

            diamond = Diamond(temp_fasta_path)
            diamond.num_threads = self.num_threads
            diamond.run.log_file_path = log_file_path
            diamond.makedb(output_db_path)
        else:
            self.run.warning("DIAMOND does not seem to be installed on this system, so anvi'o is not going to "
                             "generate a search database for it. Remember this when/if things go South.")

        if utils.is_program_exists('makeblastdb', dont_raise=True) and utils.is_program_exists('blastp', dont_raise=True):
            output_dir = J(self.COG_data_dir, 'DB_BLAST')
            if os.path.exists(output_dir):
                shutil.rmtree(output_dir)

            os.mkdir(output_dir)

            output_db_path = J(output_dir, 'COG')
            log_file_path = J(output_dir, 'log.txt')

            self.run.info('BLAST log', log_file_path)

            blast = BLAST(temp_fasta_path)
            blast.run.log_file_path = log_file_path
            blast.num_threads = self.num_threads
            blast.makedb(os.path.join(output_db_path, 'COG.fa'))
        else:
            self.run.warning("BLAST tools do not seem to be installed on this system, so anvi'o is not going to "
                             "generate a search database for them to be used. Keep this in mind for later.")

        os.remove(temp_fasta_path)


    def get_raw_data(self):
        if not os.path.exists(self.raw_NCBI_files_dir):
            os.mkdir(self.raw_NCBI_files_dir)
            open(self.COG_data_dir_version, 'w').write(COG_DATA_VERSION)

        for file_name in self.files:
            if not 'url' in self.files[file_name]:
                continue

            file_path = J(self.raw_NCBI_files_dir, file_name)
            if not os.path.exists(file_path):
                utils.download_file(self.files[file_name]['url'], file_path, progress=progress, run=run)


    def setup_raw_data(self):
        # Check hash of downloaded files before any setup
        self.check_raw_data_hash_and_existence(None, None)

        for file_name in self.files:
            file_path = J(self.raw_NCBI_files_dir, file_name)

            if not 'func' in self.files[file_name]:
                continue

            self.files[file_name]['func'](file_path, J(self.COG_data_dir, self.files[file_name]['formatted_file_name']))


    def check_raw_data_hash_and_existence(self, input_file_path, output_file_path):
        """Checks the cheksum of each downloaded file to ensure succesful download."""
        progress.new('Checking checksums and file existence')

        # Checksum file either provided by NCBI or us
        if self.COG_version == 'COG20':
            input_file_path = J(self.raw_NCBI_files_dir, "checksum.md5.txt")

        elif self.COG_version == 'COG14' or self.COG_version == 'arCOG14':
            # Get check_.md5.txt file from anvio/misc
            input_file_path = J(os.path.dirname(anvio.__file__), 'data/misc/checksum.md5.txt')

        else:
            self.run.warning(f"Anvio does not know how to check the checksums of the COG version `{self.COG_version}`."
                             f"Thus, anvi'o cannot check file integrity after download. This is not a big "
                             f"danger since failed downloads often lead to catastrophic errors instead of "
                             f"quiet omission of such problems :)")

        # Get a dictionnary of checksums, the file is formatted as "checksum filename" per line
        checksums = {}
        for line in open(input_file_path, 'rU').readlines():
            stripped = line.strip('\n').split(' ')
            file_name = stripped[-1].strip('*')
            checksums[file_name] = stripped[0]


        # For each file, check existence and check checksum
        for file_name in self.files:
            file_path = J(self.raw_NCBI_files_dir, file_name)

            # Check file exists
            if not os.path.exists(file_path):
                raise ConfigError("Something is wrong :/ Raw files are not in place...")

            # Check file present in checksum
            if file_name not in checksums.keys() and file_name != "checksum.md5.txt":
                self.run.warning(f"The file name `{file_name}` is not present in the checksum file. You should be able to "
                                 f"continue despite this, but this is unexpected.")

            # Check checksum
            if self.COG_version == 'COG20' and file_name != "checksum.md5.txt":
                if not hashlib.md5(open(file_path, "rb").read()).hexdigest() == checksums[file_name]:
                    raise ConfigError(f"Something went wrong with your download :/ The checksum we calculated for `{file_name}` "
                                      f"anvi'o just finished downloading does not match to the checksum provided by the NCBI. "
                                      f"This is most likely due to an interrupted download, as the NCBI servers often prematurely "
                                      f"end data transfers. Please try running the same command again with the `--reset` flag.")

        progress.end()
