# coding: utf-8
"""Interface to eggnog-mapper (https://github.com/jhcepas/eggnog-mapper)."""

import os
import sys
import time

import anvio
import anvio.tables as t
import anvio.utils as utils
import anvio.dbops as dbops
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import anvio.cogs as cogs

from anvio.errors import ConfigError
from anvio.tables.genefunctions import TableForGeneFunctions


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

class EggNOGMapper:
    """An interface class between eggnog-mapper and anvi'o.

       The default client of this class is `anvi-script-run-eggnog-mapper`. It may have changed already, but
       if it didn't, you should be able to run it this way to run eggnog-mapper on stuff in a contigs database,
       and store results in it:

            $ anvi-script-run-eggnog-mapper -c CONGITS_DB --num-threads NUM_THREADS

       If you are interested in some ad hoc testing with existing annotation files, you can do stuff like this after
       running `anvi-self-test --suite pangenomics` once (so you have the contigs database to annotate):

            $ ipython
            >>> from anvio.drivers.emapper import EggNOGMapper
            >>> e = EggNOGMapper('tests/sandbox/test-output/pan_test/01.db', num_threads=1)
            >>> e.populate_annotations_dict('tests/sandbox/mock_data_for_pangenomics/emapper_hits/aa_sequences_01.emapper.annotations')
            >>> e.store_annotations_in_db()

        Happy? Good. Not happy? Write to me about it!
       """

    def __init__(self, args, database='bact', executable = 'emapper.py', usemem=True, use_version=None, progress=progress, run=run):
        self.executable = executable
        self.progress = progress
        self.run = run

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.contigs_db_path = A('contigs_db')
        self.num_threads = A('num_threads')
        self.usemem = usemem


        self.COGs_data = cogs.COGsData(args)

        if not self.COGs_data.initialized:
            raise ConfigError("It seems you don't have your COG data set up on this system. Unfortunately EggNOGmapper class\
                                depends on it, so this is the end of the road for you. If you set up your COG directory to\
                                a specific path, you can use `--cog-data-dir` parameter to show anvi'o where it is. If you\
                                never set up one, then maybe it is time for you to take a look at the program\
                                `anvi-setup-ncbi-cogs`.")

        try:
            self.num_threads = int(self.num_threads) if self.num_threads else None
        except Exception as e:
            raise ConfigError("Someone didn't like the number of threads, and said [%s]. Shame on you :/" % e)

        if database not in ['euk', 'bact', 'arch']:
            raise ConfigError("Wrong database (%s). eggnog-mapper knows only about euk, bact, or arch db types..." % (database))
        else:
            self.database = database

        if self.contigs_db_path:
            utils.is_contigs_db(self.contigs_db_path)

        self.parser = None
        self.entry_id = 0
        self.installed_version = None
        self.aa_sequences_file_name = 'aa_sequences.fa'
        self.log_file_path = 'log.txt'
        self.output_file_prefix = 'project'
        self.annotations_file_name = self.output_file_prefix + '.emapper.annotations'
        self.annotations_dict = {}

        # this is a shitty workaround to make sure integers used as gene caller ids by anvi'o will not
        # cause any issues downstream (because they did in the past when silly programs started treating
        # them as numerical data and then converted them to float, and then storing them as 1.0, 2.0, etc).
        self.gene_caller_id_prefix = 'g'

        self.available_parsers = {'0.12.6': self.__parser_1}

        self.check_version(use_version)

        if not self.num_threads:
            try:
                run.warning("You have not set the number of threads, and the default is whatever the default is for eggnog-mapper. You\
                             may really want to change that since if you have a large number of genes to annotate, this may take a very\
                             long time. If you don't want to see this message again, just set the number of threads you want eggnog-mapper\
                             to use explicitly. You can press CTRL + C to cancel this run, or simply do nothing since your operation\
                             will contine in probably like 2 seconds or less ... depending how fast you read.")
                time.sleep(25)
            except KeyboardInterrupt:
                sys.exit()


    def check_version(self, use_version=None):
        """checks the installed version of eggnog-mapper, sets the parser"""

        utils.is_program_exists(self.executable)

        if not use_version:
            output, ret_code = utils.get_command_output_from_shell('%s --version' % self.executable)
            version_found = output.split('\n')[0].split('-')[1]
        else:
            version_found = use_version

        if version_found not in self.available_parsers:
            if use_version:
                raise ConfigError("Anvi'o does not know about the version you requested. Here are the ones available: %s" % \
                                                        (', '.join(list(self.available_parsers.keys()))))
            else:
                raise ConfigError("Bad news :( This version of anvi'o does not have a parser for the eggnog-mapper installed\
                                    on your system. This is the version you have on your system (if this looks totally alien\
                                    to you it may indicate another problem, in which case consider writing to anvi'o developers):\
                                    %s. For your reference, these are the versions anvi'o knows what to do with: %s" % \
                                                        (version_found, ', '.join(list(self.available_parsers.keys()))))

        self.installed_version = version_found
        self.parser = self.available_parsers[version_found]


    def add_entry(self, gene_callers_id, source, accession, function, e_value):
        self.annotations_dict[self.entry_id] = {'gene_callers_id': int(gene_callers_id),
                                                'source': source,
                                                'accession': accession,
                                                'function': function,
                                                'e_value': float(e_value)}
        self.entry_id += 1


    def __parser_1(self, defline):
        """parses this:

            0          1                    2                    3                   4                   5                     6                          7                    8                                                    9                    10      11
            query_name seed_eggNOG_ortholog seed_ortholog_evalue seed_ortholog_score predicted_gene_name GO_terms              KEGG_pathways              Annotation_tax_scope OGs                                                  bestOG|evalue|score  COG_cat eggNOG_annot
            g0         469008.B21_00002     0.0                  1825.0              THRA                GO:0003674,GO:0003824 map00061,map00540,map01100 bactNOG[38]          05C7C@bactNOG,0QIDI@gproNOG,16PBE@proNOG,COG0774@NOG 05C7C|6.1e-130|437.2       M involved in the biosynthesis of lipid A

        """

        fields = defline.strip('\n').split('\t')

        if len(fields) != 12:
            raise ConfigError("The parser for eggnog-mapper version %s does not know how to deal with this annotation fiel because the\
                                number of fields in the file (%d) is not matching to what is expected (%s)." % (self.installed_version, len(fields), 12))

        if not fields[0].startswith(self.gene_caller_id_prefix):
            raise ConfigError("Gene caller ids found in this annotation file does not start with the expected prefix. Anvi'o can not trust\
                                this file :(. This is the prefix that anvi'o expected: %s" % self.gene_caller_id_prefix)

        try:
            gene_callers_id = int(fields[0][len(self.gene_caller_id_prefix):])
        except:
            raise ConfigError("At least one gene caller id in this annotation file (%s) does not seem to look like what anvi'o is uses for\
                                gene calls. Hint: what should remain after removing gene caller id prefix (%s) should be an integer value." %\
                                                (fields[0], self.gene_caller_id_prefix))

        if fields[11] and fields[11] != 'NA' and not fields[11].startswith('Protein of unknown function'):
            self.add_entry(gene_callers_id, 'EGGNOG_%s' % self.database.upper(), fields[1], fields[11], fields[2])

        if fields[8]:
            COG_ids=[og[:-4] for og in fields[8].split(',') if og.endswith('@NOG') and og.startswith('COG')]

            if COG_ids:
                annotations = '; '.join([self.COGs_data.cogs[COG_id]['annotation'] for COG_id in COG_ids if COG_id in self.COGs_data.cogs])
                self.add_entry(gene_callers_id, 'COG_FUNCTION', ', '.join(COG_ids), annotations, 0.0)

        if fields[10]:
            self.add_entry(gene_callers_id, 'COG_CATEGORY', '', fields[10], 0.0)

        if fields[5]:
            self.add_entry(gene_callers_id, 'GO_TERMS', '', ', '.join(fields[5].split(',')), 0.0)

        if fields[6]:
            self.add_entry(gene_callers_id, 'KEGG_PATHWAYS', '', ', '.join(fields[6].split(',')), 0.0)


    def store_annotations_in_db(self, drop_previous_annotations=False):
        if not self.contigs_db_path:
            raise ConfigError("EggNOGMapper::store_annotations_in_db() is speaking: you can't really call this function if you inherited\
                                this class without a contigs database path :/ Well, yes, OK, you can call it, but it wouldn't work. Happy?")

        if not len(self.annotations_dict):
            raise ConfigError('Annotations dictionary is empty :/ There is nothing to add to the database.')

        gene_functions_table = TableForGeneFunctions(self.contigs_db_path)
        gene_functions_table.create(self.annotations_dict, drop_previous_annotations_first=drop_previous_annotations)


    def populate_annotations_dict(self, annotations_file_path):
        filesnpaths.is_file_exists(annotations_file_path)

        num_entries_processed = 0
        self.progress.new('Parisng the annotations file')
        for line in open(annotations_file_path, 'rU').readlines():
            if line.startswith('#') or line == '\n':
                continue

            self.parser(line)
            num_entries_processed += 1

            if num_entries_processed % 100 == 0:
                self.progress.update('%d ...' % num_entries_processed)

        self.progress.end()


    def process(self, output_dir, drop_previous_annotations=False):
        """Takes an anvi'o contigs database, and does its magic.

        Which involves exporting amino acid sequences for gene calls, running emapper.py on them,\
        parsing the output, and storing the results in the contigs database.
        """

        if not self.contigs_db_path:
            raise ConfigError("EggNOGMapper::process() is speaking: you can't really call this function if you inherited\
                                this class without a contigs database path :/ What are you doing?")

        filesnpaths.is_output_dir_writable(output_dir)

        contigs_db = dbops.ContigsDatabase(self.contigs_db_path)
        if not contigs_db.meta['genes_are_called']:
            raise ConfigError("It seems genes were not called for this contigs database (%s). This is a\
                                total no-no since we will need them to get amino acid seqeunces for functional\
                                annotationd :/" % self.contigs_db_path)

        aa_sequences_list = contigs_db.db.get_table_as_list_of_tuples(t.gene_amino_acid_sequences_table_name)
        num_aa_sequences = len(aa_sequences_list)
        contigs_db.disconnect()

        # change the current work directory
        work_dir = os.getcwd()
        os.chdir(output_dir)

        self.run.info('Work directory for temporary files', output_dir)
        self.run.info('Num threads to use', self.num_threads)
        self.run.info('Target database', self.database, mc='red')
        self.run.info('Use memomory', self.usemem)
        self.run.info('Genes found', num_aa_sequences, mc='green')
        self.run.info('AA sequences', self.aa_sequences_file_name)

        self.progress.new('Processing')
        self.progress.update('Storing gene sequences ...')

        aa_sequences_fp = open(self.aa_sequences_file_name, 'w')
        for gene_callers_id, aa_sequence in aa_sequences_list:
            aa_sequences_fp.write('>%s%d\n%s\n' % (self.gene_caller_id_prefix, gene_callers_id, aa_sequence))
        aa_sequences_fp.close()
        del aa_sequences_list

        cmd_line = [self.executable, '-i', self.aa_sequences_file_name, '--output', self.output_file_prefix]

        # num threads
        cmd_line.extend(['--cpu', self.num_threads]) if self.num_threads else None

        # usemem
        cmd_line.extend(['--usemem']) if self.usemem else None

        # database
        cmd_line.extend(['--database', self.database])

        self.progress.update('Running eggnog-mapper on %d sequences. This may take a while ...' % num_aa_sequences)
        utils.run_command(cmd_line, self.log_file_path)

        if not os.path.exists(self.annotations_file_name):
            self.progress.end()
            raise ConfigError("Something went wrong with eggnog-mapper :( The annotations file is not where it is supposed to be.\
                                If you are lucky, this log file will have enough output information for you to make sense of\
                                what went wrong: '%s'. Due to this error, the output directory will be kept as is, and you\
                                will have to remove it manually. Sorry about the inconvenience! Anvi'o developers know how much\
                                it sucks when things just don't work." % os.path.join(output_dir, self.log_file_path))

        self.progress.end()

        # we are done, and the annotations file is there.
        self.populate_annotations_dict(os.path.join(output_dir, self.annotations_file_name))
        os.chdir(work_dir)

        # alright. store annotations into the database
        self.store_annotations_in_db(drop_previous_annotations=drop_previous_annotations)

