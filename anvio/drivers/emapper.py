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

from anvio.errors import ConfigError
from anvio.tables.genefunctions import TableForGeneFunctions


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
        self.args = args
        self.executable = executable

        self.progress = progress
        self.run = run

        A = lambda x: args.__dict__[x] if args and x in args.__dict__ else None
        self.contigs_db_path = A('contigs_db')
        self.num_threads = A('num_threads')
        self.annotation = A('annotation')
        self.use_version = use_version
        self.usemem = usemem

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
        self.version_to_use = None
        self.aa_sequences_file_name = 'aa_sequences.fa'
        self.log_file_path = 'log.txt'
        self.output_file_prefix = 'project'
        self.annotations_file_name = self.output_file_prefix + '.emapper.annotations'
        self.annotations_dict = {}

        # this is a shitty workaround to make sure integers used as gene caller ids by anvi'o will not
        # cause any issues downstream (because they did in the past when silly programs started treating
        # them as numerical data and then converted them to float, and then storing them as 1.0, 2.0, etc).
        self.gene_caller_id_prefix = 'g'

        self.available_parsers = {'2.0.0': self.__parser_3,
                                  '2.0.1': self.__parser_3,
                                  '2.0.5': self.__parser_4,
                                  '2.1.4': self.__parser_5,
                                  '2.1.6': self.__parser_5,
                                  '2.1.8': self.__parser_5,
                                  '2.1.9': self.__parser_5}


    def check_version(self):
        """checks the installed version of eggnog-mapper, sets the parser"""

        if not self.args:
            return

        if self.annotation and not self.use_version:
            raise ConfigError("You must provide a version number to use if you have your own annotations.")
        elif not self.annotation and self.use_version:
            raise ConfigError("If you are not providing any annotations, you must let anvi'o figure out what "
                              "version of emapper to use.")


        if self.annotation:
            version_to_use = self.use_version
            pass
        else:
            utils.is_program_exists(self.executable)
            output, ret_code = utils.get_command_output_from_shell('%s --version' % self.executable)
            version_to_use = output.split('\n')[0].split('-')[1]

        if version_to_use not in self.available_parsers:
            if self.annotation:
                raise ConfigError("Anvi'o does not know about the version you requested. Here are the ones available: %s" % \
                                                        (', '.join(list(self.available_parsers.keys()))))
            else:
                raise ConfigError("Bad news :( This version of anvi'o does not have a parser for the eggnog-mapper installed "
                                   "on your system. This is the version you have on your system (if this looks totally alien "
                                   "to you it may indicate another problem, in which case consider writing to anvi'o developers): "
                                   "%s. For your reference, these are the versions anvi'o knows what to do with: %s" % \
                                                        (version_to_use, ', '.join(list(self.available_parsers.keys()))))

        self.version_to_use = version_to_use
        self.parser = self.available_parsers[version_to_use]


    def add_entry(self, gene_callers_id, source, accession, function, e_value):
        self.annotations_dict[self.entry_id] = {'gene_callers_id': int(gene_callers_id),
                                                'source': source,
                                                'accession': accession,
                                                'function': function,
                                                'e_value': float(e_value)}
        self.entry_id += 1


    def check_prefix_and_get_gene_callers_id(self, fields):
        if not fields[0].startswith(self.gene_caller_id_prefix):
            raise ConfigError("Gene caller ids found in this annotation file does not start with the expected prefix. This is a historical "
                              "glitch that is not quite easy to address programmatically, so anvi'o asks you to add the expected prefix as the "
                              "first character of every gene call in your annotations file. This is the prefix what you need to add manually "
                              "to the very beginning of every line (anvi'o developers are very sorry for this step): '%s'." % self.gene_caller_id_prefix)

        try:
            gene_callers_id = int(fields[0][len(self.gene_caller_id_prefix):])
        except:
            raise ConfigError("At least one gene caller id in this annotation file (%s) does not look like how anvi'o likes its gene calls. "
                              "Hint: what should remain after removing gene caller id prefix (%s) should be an integer value." %\
                                                (fields[0], self.gene_caller_id_prefix))

        return gene_callers_id


    def __parser_3(self, defline):
        """parses this:
        0           1                           2                     3                    4                          5                       6                    7          8        9                 10           11             12           13      14                                15    16             17               18                                                         19              20                      21
        query_name  seed_eggNOG_ortholog        seed_ortholog_evalue  seed_ortholog_score  best_tax_level             Preferred_name          GOs                  EC         KEGG_ko  KEGG_Pathway      KEGG_Module  KEGG_Reaction  KEGG_rclass  BRITE   KEGG_TC                           CAZy  BiGG_Reaction  taxonomic scope  eggNOG OGs                                                 best eggNOG OG  COG Functional cat.     eggNOG free text desc.
        g51         562982.HMPREF0432_01129     1.4e-63               248.8                Bacillales incertae sedis  rpsK                    GO:00028,GO:00462    X          ko:K048  ko03010,map03010  M00178,M0017 X              X            br01610 ko00000,ko00001,ko00002,ko03011                        Bacteria         1V3IK@1239,3WESA@539002,4HH2T@91061,COG0100@1,COG0100@2    NA|NA|NA        J                       Located on the platform of the 30S subunit, it bridges several disparate RNA helices of the 16S rRNA. Forms part of the Shine-Dalgarno cleft in the 70S ribosome
        """

        fields = defline.strip('\n').split('\t')

        if len(fields) != 22:
            raise ConfigError("The parser for eggnog-mapper version %s does not know how to deal with this annotation fiel because the "
                               "number of fields in the file (%d) is not matching to what is expected (%s)." % (self.version_to_use, len(fields), 22))

        gene_callers_id = self.check_prefix_and_get_gene_callers_id(fields)

        if (fields[21] and fields[21] != 'NA' and not fields[21].startswith('Protein of unknown function')) or fields[5]:
            if fields[5]:
                self.add_entry(gene_callers_id, 'EGGNOG_%s' % self.database.upper(), fields[1], "%s :: %s" % (fields[5], fields[21]), fields[2])
            else:
                self.add_entry(gene_callers_id, 'EGGNOG_%s' % self.database.upper(), fields[1], fields[21], fields[2])

        if fields[4]:
            self.add_entry(gene_callers_id, 'EGGNOG_BEST_TAX', fields[1], ', '.join(fields[4].split(',')), fields[2])

        if fields[5]:
            self.add_entry(gene_callers_id, 'EGGNOG_GENE_FUNCTION_NAME', fields[1], ', '.join(fields[5].split(',')), fields[2])

        if fields[6]:
            self.add_entry(gene_callers_id, 'EGGNOG_GO_TERMS', fields[1], ', '.join(fields[6].split(',')), fields[2])

        if fields[7]:
            self.add_entry(gene_callers_id, 'EGGNOG_EC_NUMBER', fields[1], ', '.join(fields[7].split(',')), fields[2])

        if fields[8]:
            self.add_entry(gene_callers_id, 'EGGNOG_KEGG_KO', fields[1], ', '.join(fields[8].split(',')), fields[2])

        if fields[9]:
            self.add_entry(gene_callers_id, 'EGGNOG_KEGG_PATHWAYS', fields[1], ', '.join(fields[9].split(',')), fields[2])

        if fields[10]:
            self.add_entry(gene_callers_id, 'EGGNOG_KEGG_MODULE', fields[1], ', '.join(fields[10].split(',')), fields[2])

        if fields[13]:
            self.add_entry(gene_callers_id, 'EGGNOG_BRITE', fields[1], ', '.join(fields[13].split(',')), fields[2])

        if fields[14]:
            self.add_entry(gene_callers_id, 'EGGNOG_KEGG_TC', fields[1], ', '.join(fields[14].split(',')), fields[2])

        if fields[15]:
            self.add_entry(gene_callers_id, 'EGGNOG_CAZy', fields[1], ', '.join(fields[15].split(',')), fields[2])

        if fields[16]:
            self.add_entry(gene_callers_id, 'EGGNOG_BiGG_REACTIONS', fields[1], ', '.join(fields[16].split(',')), fields[2])


    def __parser_4(self, defline):
        """parses this:

           0             1                         2                       3                     4                                                                                             5                                6              7                                  8                                 9              10                                 11                12                             13          14           15                                  16             17              18                 19                                                 20         21      22                                                         23
           query_name    seed_eggNOG_ortholog      seed_ortholog_evalue    seed_ortholog_score   eggNOG OGs                                                                                    narr_og_name                     narr_og_cat    narr_og_desc                       best_og_name                      best_og_cat    best_og_desc                       Preferred_name    GOs                            EC          KEGG_ko      KEGG_Pathway                        KEGG_Module    KEGG_Reaction   KEGG_rclass        BRITE                                              KEGG_TC    CAZy    BiGG_Reaction                                              PFAMs
           g84988        573.JG24_24155            1.1e-178                632.5                 COG0463@1|root,COG0463@2|Bacteria,1MWE5@1224|Proteobacteria,1RPCE@1236|Gammaproteobacteria    1RPCE@1236|Gammaproteobacteria   M              Catalyzes the transfer of (...)    1RPCE@1236|Gammaproteobacteria    M              Catalyzes the transfer of (...)    arnC              GO:0003674,GO:0003824,(...)    2.4.2.53    ko:K10012    ko00520,ko01503,map00520,map01503   M00721,M00761  R07661          RC00005,RC02954    ko00000,ko00001,ko00002,ko01000,ko01005,ko02000    4.D.2.1.8  GT2     iAPECO1_1312.APECO1_4307,iE2348C_1286.E2348C_2398,(...)    Glyco_tranf_2_2,Glyco_tranf_2_3,Glycos_transf_2
        """

        fields = defline.strip('\n').split('\t')

        if len(fields) != 24:
            raise ConfigError("The parser for eggnog-mapper version %s does not know how to deal with this annotation fiel because the "
                               "number of fields in the file (%d) is not matching to what is expected (%s)." % (self.version_to_use, len(fields), 24))

        gene_callers_id = self.check_prefix_and_get_gene_callers_id(fields)

        #
        F = lambda x: fields[x] and fields[x] != '-'

        if (F(21) and not fields[21].startswith('Protein of unknown function')) or F(10):
            if F(11):
                self.add_entry(gene_callers_id, 'EGGNOG_%s' % self.database.upper(), fields[1], "%s :: %s" % (fields[11], fields[10]), fields[2])
            else:
                self.add_entry(gene_callers_id, 'EGGNOG_%s' % self.database.upper(), fields[1], fields[10], fields[2])

        if F(4):
            try:
                # because there is crap like this where the delimiter (',') is used in text fields ....
                #
                #    "4QF6A@10239|Viruses,4QV1N@35237|dsDNA viruses, no RNA stage,4QRT3@28883|Caudovirales"
                #
                tax = ' / '.join([x.split('|')[1] for x in fields[4].split(',')])
                self.add_entry(gene_callers_id, 'EGGNOG_BEST_TAX', fields[1], tax, fields[2])
            except:
                pass

        if F(11):
            self.add_entry(gene_callers_id, 'EGGNOG_GENE_FUNCTION_NAME', fields[1], ', '.join(fields[11].split(',')), fields[2])

        if F(12):
            self.add_entry(gene_callers_id, 'EGGNOG_GO_TERMS', fields[1], ', '.join(fields[12].split(',')), fields[2])

        if F(13):
            self.add_entry(gene_callers_id, 'EGGNOG_EC_NUMBER', fields[1], ', '.join(fields[13].split(',')), fields[2])

        if F(14):
            self.add_entry(gene_callers_id, 'EGGNOG_KEGG_KO', fields[1], ', '.join(fields[14].split(',')), fields[2])

        if F(15):
            self.add_entry(gene_callers_id, 'EGGNOG_KEGG_PATHWAYS', fields[1], ', '.join(fields[15].split(',')), fields[2])

        if F(16):
            self.add_entry(gene_callers_id, 'EGGNOG_KEGG_MODULE', fields[1], ', '.join(fields[16].split(',')), fields[2])

        if F(19):
            self.add_entry(gene_callers_id, 'EGGNOG_BRITE', fields[1], ', '.join(fields[19].split(',')), fields[2])

        if F(20):
            self.add_entry(gene_callers_id, 'EGGNOG_KEGG_TC', fields[1], ', '.join(fields[20].split(',')), fields[2])

        if F(21):
            self.add_entry(gene_callers_id, 'EGGNOG_CAZy', fields[1], ', '.join(fields[21].split(',')), fields[2])

        if F(22):
            self.add_entry(gene_callers_id, 'EGGNOG_BiGG_REACTIONS', fields[1], ', '.join(fields[22].split(',')), fields[2])


    def __parser_5(self, defline):
        """parses this:

           0             1                         2                       3                     4                                                                                             5                                6              7                                  8                                 9              10                                 11                    12                         13           14               15                16                               17         18      19               20
           query         seed_ortholog             evalue                  score                 eggNOG_OGs                                                                                    max_annot_lvl                    COG_category   Description                        Preferred_name                    GOs            EC                                 KEGG_ko               KEGG_Pathway               KEGG_Module  KEGG_Reaction    KEGG_rclass       BRITE                            KEGG_TC    CAZy    BiGG_Reaction    PFAMs
           g2183         858619.CVAR_2845          2.4e-178                632.1                 COG1021@1|root,COG1021@2|Bacteria,2I2IR@201174|Actinobacteria,22MBM@1653|Corynebacteriaceae   201174|Actinobacteria            Q	           23-dihydroxybenzoate-AMP ligase    entE                              -              2.7.7.58,6.3.2.14                  ko:K02363,ko:K12238   ko01053,ko01110,(...)      -            R07644           RC00162,RC03046   ko00000,ko00001,ko01000,ko01008  -          -	   -                AMP-binding,AMP-binding_C

        """

        fields = defline.strip('\n').split('\t')

        if len(fields) != 21:
            raise ConfigError("The parser for eggnog-mapper version %s does not know how to deal with this annotation fiel because the "
                               "number of fields in the file (%d) is not matching to what is expected (%s)." % (self.version_to_use, len(fields), 21))

        gene_callers_id = self.check_prefix_and_get_gene_callers_id(fields)

        #
        F = lambda x: fields[x] and fields[x] != '-'
        A = lambda x, y: self.add_entry(gene_callers_id, y, fields[1], ', '.join(fields[x].split(',')), fields[2]) if (fields[x] and fields[x] != '-') else None

        if (F(7) and not fields[7].startswith('Protein of unknown function')):
            if F(8):
                self.add_entry(gene_callers_id, 'EGGNOG_%s' % self.database.upper(), fields[1], "%s :: %s" % (fields[8], fields[7]), fields[2])
            else:
                A(7, 'EGGNOG_%s' % self.database.upper())

        if F(4):
            try:
                # because there is crap like this where the delimiter (',') is used in text fields ....
                #
                #    "4QF6A@10239|Viruses,4QV1N@35237|dsDNA viruses, no RNA stage,4QRT3@28883|Caudovirales"
                #
                tax = ' / '.join([x.split('|')[1] for x in fields[4].split(',')])
                self.add_entry(gene_callers_id, 'EGGNOG_BEST_TAX', fields[1], tax, fields[2])
            except:
                pass

        A(8, 'EGGNOG_GENE_FUNCTION_NAME')
        A(9, 'EGGNOG_GO_TERMS')
        A(6, 'EGGNOG_COG_CATEGORY')
        A(10, 'EGGNOG_EC_NUMBER')
        A(11, 'EGGNOG_KEGG_KO')
        A(12, 'EGGNOG_KEGG_PATHWAYS')
        A(13, 'EGGNOG_KEGG_MODULE')
        A(13, 'EGGNOG_KEGG_REACTION')
        A(15, 'EGGNOG_KEGG_RCLASS')
        A(16, 'EGGNOG_BRITE')
        A(17, 'EGGNOG_KEGG_TC')
        A(18, 'EGGNOG_CAZy')
        A(19, 'EGGNOG_BiGG_REACTIONS')
        A(20, 'EGGNOG_PFAMs')


    def store_annotations_in_db(self, drop_previous_annotations=False):
        if not self.contigs_db_path:
            raise ConfigError("EggNOGMapper::store_annotations_in_db() is speaking: you can't really call this function if you inherited "
                               "this class without a contigs database path :/ Well, yes, OK, you can call it, but it wouldn't work. Happy?")

        if not len(self.annotations_dict):
            raise ConfigError('Annotations dictionary is empty :/ There is nothing to add to the database.')

        gene_functions_table = TableForGeneFunctions(self.contigs_db_path)
        gene_functions_table.create(self.annotations_dict, drop_previous_annotations_first=drop_previous_annotations)


    def populate_annotations_dict(self, annotations_file_path):
        filesnpaths.is_file_exists(annotations_file_path)

        self.check_version()

        num_entries_processed = 0
        self.progress.new('Parsing the annotations file')
        for line in open(annotations_file_path, 'r').readlines():
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

        # check number of threads first.
        if not self.num_threads:
            try:
                run.warning("You have not set the number of threads, and the default is whatever the default is for eggnog-mapper. You "
                            "may really want to change that since if you have a large number of genes to annotate, this may take a very "
                            "long time. If you don't want to see this message again, just set the number of threads you want eggnog-mapper "
                            "to use explicitly. You can press CTRL + C to cancel this run, or simply do nothing since your operation "
                            "will contine in probably like 10 seconds or less ... depending how fast you read.")
                time.sleep(25)
            except KeyboardInterrupt:
                sys.exit()

        if not self.contigs_db_path:
            raise ConfigError("EggNOGMapper::process() is speaking: you can't really call this function if you inherited "
                               "this class without a contigs database path :/ What are you doing?")

        filesnpaths.is_output_dir_writable(output_dir)

        contigs_db = dbops.ContigsDatabase(self.contigs_db_path)
        if not contigs_db.meta['genes_are_called']:
            raise ConfigError("It seems genes were not called for this contigs database (%s). This is a "
                               "total no-no since we will need them to get amino acid seqeunces for functional "
                               "annotationd :/" % self.contigs_db_path)

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
            raise ConfigError("Something went wrong with eggnog-mapper :( The annotations file is not where it is supposed to be. "
                               "If you are lucky, this log file will have enough output information for you to make sense of "
                               "what went wrong: '%s'. Due to this error, the output directory will be kept as is, and you "
                               "will have to remove it manually. Sorry about the inconvenience! Anvi'o developers know how much "
                               "it sucks when things just don't work." % os.path.join(output_dir, self.log_file_path))

        self.progress.end()

        # we are done, and the annotations file is there.
        self.populate_annotations_dict(os.path.join(output_dir, self.annotations_file_name))
        os.chdir(work_dir)

        # alright. store annotations into the database
        self.store_annotations_in_db(drop_previous_annotations=drop_previous_annotations)

