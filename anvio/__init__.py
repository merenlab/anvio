# -*- coding: utf-8
# pylint: disable=line-too-long

"""Lots of under-the-rug, operational garbage in here. Run. Run away.."""

import sys
import json
import copy
import platform

# yes, this library is imported but never used, but don't remove it
# unless you want to explode `bottle`:
import pkg_resources

anvio_version = '6'
anvio_codename = 'esther'

DEBUG = '--debug' in sys.argv
FORCE = '--force' in sys.argv

def P(d, dont_exit=False):
    """Poor man's debug output printer during debugging."""

    print(json.dumps(d, indent=2))

    if not dont_exit:
        sys.exit()


# Make sure the Python environment hasn't changed since the installation (happens more often than you'd think
# on systems working with multiple Python installations that are managed through modules):
try:
    if sys.version_info.major != 3 or sys.version_info.minor < 5:
        sys.stderr.write("Sad face :( Your active Python version is %s, but anvi'o only works with Python version 3.5.0 or later.\n" % (platform.python_version()))
        sys.exit(-1)
except Exception:
    sys.stderr.write("(anvi'o failed to learn about your Python version, but it will pretend as if nothing happened)\n\n")


def get_args(parser):
    """A helper function to parse args anvi'o way.

       This function allows us to make sure some ad hoc parameters such as `--debug`
       can be used with any anvi'o program spontaneously even if they are not explicitly
       defined as an accepted argument, yet flags (or parameters) anvi'o does not expect
       to see can still be sorted out.
    """

    allowed_ad_hoc_flags = ['--version', '--debug', '--force']

    args, unknown = parser.parse_known_args()

    # if there are any args in the unknown that we do not expect to find
    # we we will make argparse complain about those.
    if len([f for f in unknown if f not in allowed_ad_hoc_flags]):
        for f in allowed_ad_hoc_flags:
            parser.add_argument(f, action='store_true')
        parser.parse_args()

    return args


import anvio.tables as tables
import anvio.constants as constants


# a comprehensive arguments dictionary that provides easy access from various programs that interface anvi'o modules:
D = {
    'profile-db': (
            ['-p', '--profile-db'],
            {'metavar': "PROFILE_DB",
             'required': True,
             'help': "Anvi'o profile database"}
                ),
    'pan-db': (
            ['-p', '--pan-db'],
            {'metavar': "PAN_DB",
             'required': True,
             'help': "Anvi'o pan database"}
                ),
    'pan-or-profile-db': (
            ['-p', '--pan-or-profile-db'],
            {'metavar': "PAN_OR_PROFILE_DB",
             'required': True,
             'help': "Anvi'o pan or profile database (and even genes database in appropriate contexts)."}
                ),
    'genomes-storage': (
            ['-g', '--genomes-storage'],
            {'metavar': "GENOMES_STORAGE",
             'required': False,
             'help': "Anvi'o genomes storage file"}
                ),
    'structure-db': (
            ['-s', '--structure-db'],
            {'metavar': "STRUCTURE_DB",
             'required': True,
             'help': "Anvi'o structure database."}
                ),
    'only-if-structure': (
            ['--only-if-structure'],
            {'default': False,
             'action': 'store_true',
             'help': "If provided, your genes of interest will be further subset to only include\
                      genes with structures in your structure database, and therefore must be supplied in\
                      conjunction with a structure database, i.e. `-s <your_structure_database>`. If you did\
                      not specify genes of interest, ALL genes will be subset to those that have\
                      structures."}
                ),
    'genomes-names': (
            ['-G', '--genome-names'],
            {'metavar': "GENOME_NAMES",
             'required': False,
             'help': "Genome names to 'focus'. You can use this parameter to limit the genomes included in your analysis.\
                      You can provide these names as a comma-separated list of names, or you can put them in a file,\
                      where you have a single genome name in each line, and provide the file path."}
                ),
    'blank-profile': (
            ['--blank-profile'],
            {'default': False,
             'action': 'store_true',
             'help': "If you only have contig sequences, but no mapping data (i.e., you found a genome and would like to\
                      take a look from it), this flag will become very hand. After creating a contigs database for your\
                      contigs, you can create a blank anvi'o profile database to use anvi'o interactive\
                      interface with that contigs database without any mapping data."}
                ),
    'contigs-db': (
            ['-c', '--contigs-db'],
            {'metavar': 'CONTIGS_DB',
             'required': True,
             'help': "Anvi'o contigs database generated by 'anvi-gen-contigs'"}
                ),
    'runinfo': (
            ['-r', '--runinfo'],
            {'metavar': 'RUNINFO_PATH',
             'required': True,
             'help': "Anvi'o runinfo file path."}
                ),
    'description': (
            ['--description'],
            {'metavar': 'TEXT_FILE',
             'required': False,
             'help': "A plain text file that contains some description about the project. You can use Markdwon syntax.\
                      The description text will be rendered and shown in all relevant interfaces, including the\
                      anvi'o interactive interface, or anvi'o summary outputs."}
                ),

    'additional-view': (
            ['-V', '--additional-view'],
            {'metavar': 'ADDITIONAL_VIEW',
             'help': "A TAB-delimited file for an additional view to be used in the interface. This\
                      file should contain all split names, and values for each of them in all\
                      samples. Each column in this file must correspond to a sample name. Content\
                      of this file will be called 'user_view', which will be available as a new item\
                      in the 'views' combo box in the interface"}
                ),
    'fasta-file': (
            ['-f', '--fasta-file'],
            {'metavar': 'FASTA',
             'help': "A FASTA-formatted input file"}
                ),
    'fasta-text-file': (
            ['-f', '--fasta-text-file'],
            {'metavar': 'FASTA_TEXT_FILE',
            'dest': 'fasta_text_file',
            'help': "A two-column TAB-delimited file that lists multiple FASTA files to import\
                     for analysis. If using for `anvi-dereplicate-genomes` or `anvi-compute-distance`,\
                     each FASTA is assumed to be a genome. The first item in the header line\
                     should read 'name', and the second item should read 'fasta_path'. Each line\
                     in the field should describe a single entry, where the first column is\
                     the name of the FASTA file or corresponding sequence, and the second column\
                     is the path to the FASTA file itself."}
                ),
    'layers-information-file': (
            ['-D', '--layers-information-file'],
            {'metavar': 'FILE',
             'help': "A TAB-delimited file with information about layers in your dataset. Each row in this\
                      file must correspond to a sample name. Each column must contain a unique attribute.\
                      Please refer to the documentation to learn more about the structure and purpose of\
                      this file."}
                ),
    'layers-order-file': (
            ['-R', '--layers-order-file'],
            {'metavar': 'FILE',
             'help': "A TAB-delimited file with three columns: 'attribute', 'basic', 'newick'. For each attribute,\
                      the order of samples must be defined either in the 'basic' form or via a 'newick'-formatted\
                      tree structure that describes the organization of each sample. Anvi'o will look for a\
                      comma-separated list of sample names for the 'basic' form. Please refer to the online docs\
                      for more info. Also you shouldn't hesitate to try to find the right file format until you get\
                      it working. There are stringent checks on this file, and you will not break anything while trying!."}
                ),
    'split-length': (
            ['-L', '--split-length'],
            {'metavar': 'INT',
             'default': 20000,
             'type': int,
             'help': "Anvi'o splits very long contigs into smaller pieces, without actually splitting them for real. These\
                      'virtual' splits improve the efficacy of the visualization step, and changing the split size gives\
                      freedom to the user to adjust the resolution of their display when necessary. The default value is\
                      (%(default)d). If you are planning to use your contigs database for metagenomic binning, we advise you\
                      to not go below 10,000 (since the lower the split size is, the more items to show in the display, and\
                      decrasing the split size does not really help much to binning). But if you are thinking about using this\
                      parameter for ad hoc investigations other than binning, you should ignore our advice, and set the split\
                      size as low as you want. If you do not want your contigs to be split, you can set the split size to '0'\
                      or any other negative integer (lots of unnecessary freedom here, enjoy!)."}
                ),
    'kmer-size': (
            ['-K', '--kmer-size'],
            {'metavar': 'INT',
             'default': 4,
             'type': int,
             'help': "K-mer size for k-mer frequency calculations. The default k-mer size for composition-based\
                      analyses is 4, historically. Although tetra-nucleotide frequencies seem to offer the\
                      the sweet spot of sensitivity, information density, and manageable number of dimensions\
                      for clustering approaches, you are welcome to experiment (but maybe you should leave\
                      it as is for your first set of analyses)."}
                ),
    'prodigal-translation-table': (
            ['--prodigal-translation-table'],
            {'metavar': 'INT',
             'default': None,
             'help': "This is a parameter to pass to the Prodigal for a specific translation table. This parameter\
                      corresponds to the parameter `-g` in Prodigal, the default value of which is 11 (so if you do\
                      not set anything, it will be set to 11 in Prodigal runtime. Please refer to the Prodigal\
                      documentation to determine what is the right translation table for you if you think you need\
                      it.)"}
                ),

    'skip-gene-calling': (
            ['--skip-gene-calling'],
            {'default': False,
             'action': 'store_true',
             'help': "By default, generating an anvi'o contigs database includes the identification of open reading\
                      frames in contigs by running a bacterial gene caller. Declaring this flag will by-pass that\
                      process. If you prefer, you can later import your own gene calling results into the database."}
                ),
    'remove-partial-hits': (
            ['--remove-partial-hits'],
            {'default': False,
             'action': 'store_true',
             'help': "By default anvi'o will return hits even if they are partial. Declaring this flag will make\
                      anvi'o filter all hits that are partial. Partial hits are hits in which you asked for n1\
                      genes before and n2 genes after the gene that matched the search criteria but the search\
                      hits the end of the contig before finding the number of genes that you asked."}
            ),
    'never-reverse-complement': (
            ['--never-reverse-complement'],
            {'default': False,
             'action': 'store_true',
             'help': "By default, if a gene that is found by the search criteria is reverse in it's direction,\
                      then the sequence of the entire locus is reversed before it is saved to the output.\
                      If you wish to prevent this behavior then use the flag --never-reverse-complement.",}
             ),
    'zeros-are-outliers': (
            ['--zeros-are-outliers'],
            {'default': False,
             'action': 'store_true',
             'help': "If you want all zero coverage positions to be treated like outliers\
                      then use this flag. The reason to treat zero coverage as outliers\
                      is because when mapping reads to a reference we could get many zero\
                      positions due to accessory genes. These positions then skew the average\
                      values that we compute."}
            ),
    'outliers-threshold': (
            ['--outliers-threshold'],
            {'default': 1.5,
             'type': float,
             'metavar': 'NUM',
             'help': "Threshold to use for the outlier detection. The default value is '%(default).1f'.\
                      Absolute deviation around the median is used. To read more about the method please\
                      refer to: 'How to Detect and Handle Outliers' by Boris Iglewicz and David Hoaglin \
                      (doi:10.1016/j.jesp.2013.03.013)."}
            ),
    'external-gene-calls': (
            ['--external-gene-calls'],
            {'metavar': 'GENE-CALLS',
             'help': "A TAB-delimited file to utilize external gene calls. The file must have these columns: 'gene_callers_id'\
                      (a unique integer number for each gene call, start from 1), 'contig' (the contig name the gene call is found),\
                      'start' (start position, integer), 'stop' (stop position, integer), 'direction' (the direction of the gene open reading\
                      frame; can be 'f' or 'r'), 'partial' (whether it is a complete gene call, or a partial one; must be 1 for partial\
                      calls, and 0 for complete calls), 'source' (the gene caller), and 'version' (the version of the gene caller, i.e.,\
                      v2.6.7 or v1.0). An example file can be found via the URL https://goo.gl/TqCWT2"}
                ),
    'external-genomes': (
            ['-e', '--external-genomes'],
            {'metavar': 'FILE_PATH',
             'help': "A two-column TAB-delimited flat text file that lists anvi'o contigs databases. The first item\
                      in the header line should read 'name', and the second should read 'contigs_db_path'. Each line in the\
                      file should describe a single entry, where the first column is the name of the genome (or MAG), and\
                      the second column is the anvi'o contigs database generated for this genome."}
                ),
    'internal-genomes': (
            ['-i', '--internal-genomes'],
            {'metavar': 'FILE_PATH',
             'help': "A five-column TAB-delimited flat text file. The header line must contain thse columns: 'name', 'bin_id',\
                      'collection_id', 'profile_db_path', 'contigs_db_path'. Each line should list a single entry, where 'name'\
                      can be any name to describe the anvi'o bin identified as 'bin_id' that is stored in a collection."}
                ),
    'gene-caller': (
            ['--gene-caller'],
            {'metavar': 'GENE-CALLER',
             'help': "The gene caller to utilize. Anvi'o supports multiple gene callers, and some operations (including this one)\
                      requires an explicit mentioning of which one to use. The default is '%s', but it will not be enough if you\
                      if you were a rebel and have used `--external-gene-callers` or something." % constants.default_gene_caller}
                ),

    'ignore-internal-stop-codons': (
            ['--ignore-internal-stop-codons'],
            {'default': False,
             'action': 'store_true',
             'help': "This is only relevant when you have an external gene calls file. If anvi'o figures out that your custom gene calls\
                      result in amino acid sequences with stop codons in the middle, it will complain about it. You can use this flag\
                      to tell anvi'o to don't check for internal stop codons, EVEN THOUGH IT MEANS THERE IS MOST LIKELY SOMETHING\
                      WRONG WITH YOUR EXTERNAL GENE CALLS FILE. Anvi'o will understand that sometimes we don't want to care, and will\
                      not judge you. Instead, it will replace every stop codon residue in the amino acid sequence with an 'X' character.\
                      Please let us know if you used this and things failed, so we can tell you that you shouldn't have really used it\
                      if you didn't like failures at the first place (smiley)."}
                ),
    'get-samples-stats-only': (
            ['--get-samples-stats-only'],
            {'default': False,
             'action': 'store_true',
             'help': "If you only wish to get statistics regarding the occurrence of bins in samples, then use this flag. \
                      Especially when dealing with many samples or large genomes, gene stats could be a long time to compute. \
                      By using this flag you could save a lot of computation time."}
                ),
    'gen-figures': (
            ['--gen-figures'],
            {'default': False,
             'action': 'store_true',
             'help': "For those of you who wish to dig deeper, a collection of figures could be created to allow\
                      you to get insight into how the classification was generated. This is especially useful to\
                      identify cases in which you shouldn't trust the classification (for example due to a large\
                      number of outliers). NOTICE: if you ask anvi'o to generate these figures then it will\
                      significantly extend the execution time. To learn about which figures are created and what\
                      they mean, contact your nearest anvi'o developer, because currently it is a well-hidden secret."}
                ),
    'skip-SNV-profiling': (
            ['--skip-SNV-profiling'],
            {'default': False,
             'action': 'store_true',
             'help': "By default, anvi'o characterizes single-nucleotide variation in each sample. The use of this flag\
                      will instruct profiler to skip that step. Please remember that parameters and flags must be\
                      identical between different profiles using the same contigs database for them to merge properly."}
                ),
    'return-AA-frequencies-instead': (
            ['--return-AA-frequencies-instead'],
            {'default': False,
             'action': 'store_true',
             'help': "By default, anvi'o will return codon frequencies (as the name suggests), but you can ask for amino\
                      acid frequencies instead, simply because you always need more data and more stuff. You're lucky\
                      this time, but is there an end to this? Will you ever be satisfied with what you have?\
                      Anvi'o needs answers."}
                ),
    'profile-SCVs': (
            ['--profile-SCVs'],
            {'default': False,
             'action': 'store_true',
             'help': "Anvi'o can perform accurate characterization of codon frequencies in genes during profiling. While having\
                      codon frequencies opens doors to powerful evolutionary insights in downstream analyses, due to its\
                      computational complexity, this feature comes 'off' by default. Using this flag you can rise against the\
                      authority, as you always should, and make anvi'o profile codons."}
                ),
    'ignore-orphans': (
            ['--ignore-orphans'],
            {'default': False,
             'action': 'store_true',
             'help': "Ignore orphan reads (paired reads that are not in a proper pair). The default is to include orphans."}
                ),
    'max-coverage-depth': (
            ['-m', '--max-coverage-depth'],
            {'default': 8000,
             'metavar': 'INT',
             'type': int,
             'help': "Max depth of coverage to consider when reading from the BAM file. It means, nucleotide positions with\
                      coverages that exceed this value will have a flat coverage that is equal to this value. The default\
                      is %(default)d."}
                ),
    'drop-previous-annotations': (
            ['--drop-previous-annotations'],
            {'default': False,
             'action': 'store_true',
             'help': "Use this flag if you want anvi'o to remove ALL previous functional annotations for your genes,\
                      and then import the new data. The default behavior will add any annotation source into the db\
                      incrementally unless there are already annotations from this source. In which case, it will first\
                      remove previous annotations for that source only (i.e., if source X is both in the db and in the\
                      incoming annotations data, it will replace the content of source X in the db)."}
                ),
    'skip-mindful-splitting': (
            ['--skip-mindful-splitting'],
            {'default': False,
             'action': 'store_true',
             'help': "By default, anvi'o attempts to prevent soft-splitting large contigs by cutting proper gene calls\
                      to make sure a single gene is not broken into multiple splits. This requires a careful\
                      examination of where genes start and end, and to find best locations to split contigs with respect\
                      to this information. So, when the user asks for a split size of, say, 1,000, it serves as a\
                      mere suggestion. When this flag is used, anvi'o does what the user wants and creates splits at\
                      desired lengths (although some functionality may become unavailable for the projects that rely on\
                      a contigs database that is initiated this way)."}
                ),
    'contigs-fasta': (
            ['-f', '--contigs-fasta'],
            {'metavar': 'FASTA',
             'required': True,
             'help': "The FASTA file that contains reference sequences you mapped your samples against. This\
                      could be a reference genome, or contigs from your assembler. Contig names in this file\
                      must match to those in other input files. If there is a problem anvi'o will gracefully\
                      complain about it."}
                ),
    'view-data': (
            ['-d', '--view-data'],
            {'metavar': 'VIEW_DATA',
             'help': "A TAB-delimited file for view data"}
                ),
    'tree': (
            ['-t', '--tree'],
            {'metavar': 'NEWICK',
             'help': "NEWICK formatted tree structure"}
                ),
    'items-order': (
            ['--items-order'],
            {'metavar': 'FLAT_FILE',
             'help': "A flat file that contains the order of items you wish the display using the interactive interface. You\
                      may want to use this if you have a specific order of items in your mind, and do not want to display a\
                      tree in the middle (or simply you don't have one). The file format is simple: each line should have an\
                      item name, and there should be no header."}
                ),
    'additional-layers': (
            ['-A', '--additional-layers'],
            {'metavar': 'ADDITIONAL_LAYERS',
             'help': "A TAB-delimited file for additional layers for splits. The first column of this file\
                      must be split names, and the remaining columns should be unique attributes.\
                      The file does not need to contain all split names, or values for each split in\
                      every column. Anvi'o will try to deal with missing data nicely. Each column in this\
                      file will be visualized as a new layer in the tree."}
                ),
    'target-data-group': (
            ['-D', '--target-data-group'],
            {'metavar': 'NAME',
             'default': None,
             'help': "Data group to focus. Anvi'o misc data tables support associating a set of data keys\
                      with a data group. If you have no idea what this is, then probably you don't need it,\
                      and anvi'o will take care of you. Note: this flag is IRRELEVANT if you are working with\
                      additioanl order data tables."}
                ),
    'target-data-table': (
            ['-t', '--target-data-table'],
            {'metavar': 'NAME',
             'help': "The target table is the table you are interested in accessing. Currently it can be 'items','layers', or\
                      'layer_orders'. Please see most up-to-date online documentation for more information."}
                ),
    'view': (
            ['--view'],
            {'metavar': 'NAME',
             'help': "Start the interface with a pre-selected view. To see a list of available views,\
                      use --show-views flag."}
                ),
    'category-variable': (
            ['--category-variable'],
            {'default': None,
             'metavar': 'CATEGORY',
             'help': "The additional layers data variable name that divides layers into multiple categories."}
                ),
    'exclude-ungrouped': (
            ['--exclude-ungrouped'],
            {'default': False,
             'action': 'store_true',
             'help': "Use this flag if you want anvi'o to ignore genomes with no value set for the catergory variable \
                      (which you specified using --category-variable). By default all variables with no value will be \
                      considered as a single group when preforming the statistical analysis."}
                ),
    'functional-occurrence-table-output': (
            ['-F', '--functional-occurrence-table-output'],
            {'metavar': 'FILE',
             'default': None,
             'type': str,
             'help': "Saves the occurrence frequency information for functions in genomes in a TAB-delimited format.\
                      A file name must be provided. To learn more about how the functional occurrence is computed, please\
                      refer to the tutorial."}
                ),
    'table': (
            ['--table'],
            {'metavar': 'TABLE_NAME',
             'help': "Table name to export."}
                ),
    'fields': (
            ['-f', '--fields'],
            {'metavar': 'FIELD(S)',
             'help': "Fields to report. USe --list-tables parameter with a table name to see available\
                      fields  You can list fields using this notation: --fields 'field_1, field_2, ... field_N'."}
                ),
    'list': (
            ['-l', '--list'],
            {'default': False,
             'action': 'store_true',
             'help': "Gives a list of tables in a database and quits. If a table is already declared\
                      this time it lists all the fields in a given table, in case you would to export\
                      only a specific list of fields from the table using --fields parameter."}
                ),
    'title': (
            ['--title'],
            {'metavar': 'NAME',
             'help': "Title for the interface. If you are working with a RUNINFO dict, the title\
                      will be determined based on information stored in that file. Regardless,\
                      you can override that value using this parameter."}
                ),
    'split-hmm-layers': (
            ['--split-hmm-layers'],
            {'default': False,
             'action': 'store_true',
             'help': "When declared, this flag tells the interface to split every gene found in HMM\
                      searches that were performed against non-singlecopy gene HMM profiles into\
                      their own layer. Please see the documentation for details."}
                ),
    'taxonomic-level': (
            ['--taxonomic-level'],
            {'default': 't_genus',
             'type': str,
             'choices': tables.taxon_names_table_structure[1:],
             'help': "The taxonomic level to use. The default is '%(default)s'. Only relevant if the\
                      anvi'o ontigs database contains taxonomic annotations."}
                ),
    'taxonomy-file': (
            ['-t', '--taxonomy-file'],
            {'default': None,
             'type': str,
             'help': "Path to The taxonomy file format tsv containe:\
              ID\td__domaine;p__phylum;[..];s__genus species"}
                ),
    'metagenome-mode': (
            ['-m', '--metagenome-mode'],
            {'default': False,
             'action': 'store_true',
             'help': "Treat a given contigs database as a metagenome rather than treating it as a single genome."}
                ),
    'scg-name-for-metagenome-mode': (
            ['-S','--scg-name-for-metagenome-mode'],
            {'default': None,
             'type': str,
             'metavar': 'SCG_NAME',
             'help': "When running in metagenome mode, anvi'o automatically chooses the most frequent single-copy\
                      core gene to estimate the taxonomic composition within a contigs database. If you have a\
                      different preference you can use this parameter to communicate that."}
                ),
    'compute-scg-coverages': (
            ['--compute-scg-coverages'],
            {'default': False,
             'action': 'store_true',
             'help': "When this flag is declared, anvi'o will go back to the profile database to learn coverage\
                      statistics of single-copy core genes for which we have taxonomy information."}
                ),
    'update-profile-db-with-taxonomy': (
            ['--update-profile-db-with-taxonomy'],
            {'default': False,
             'action': 'store_true',
             'help': "When anvi'o knows all both taxonomic affiliations and coverages across samples for single-copy\
                      core genes, it can, in theory add this information to the profile database. With this flag you\
                      can isntruct anvi'o to do that and find information on taxonomy in the `layers` tab of your\
                      interactive interface."}
                ),
    'taxonomy-database': (
            ['-r', '--taxonomy-database'],
            {'default': None,
             'type': str,
             'metavar': 'PATH',
             'help': "Path to the directory that contains the BLAST databases for single-copy core\
                      genes. You will almost never need to use this parameter unless you are\
                      trying something very fancy. Anvi'o will know where its database files are."}
                ),
    'scgs-taxonomy-data-dir': (
            ['--scgs-taxonomy-data-dir'],
            {'default': None,
             'type': str,
             'metavar': 'PATH',
             'help': "The directory for SCGs data to be stored (or read from, depending on the context).\
                      If you leave it as is without specifying anything, anvi'o will set up everything in\
                      (or try to read things from) a pre-defined default directory. The advantage of using\
                      the default directory at the time of set up is that every user of anvi'o on a computer\
                      system will be using a single data directory, but then you may need to run the setup\
                      program with superuser privileges. If you don't have superuser privileges, then you can\
                      use this parameter to tell anvi'o the location you wish to use to setup your databases.\
                      If you are using a program (such as `anvi-run-scg-taxonomy` or `anvi-estimate-genome-taxonomy`)\
                      you will have to use this parameter to tell those programs where your data are."}
                ),
    'scgs-taxonomy-remote-database-url': (
            ['--scgs-taxonomy-remote-database-url'],
            {'default': None,
             'type': str,
             'metavar': 'URL',
             'help': "Anvi'o will always try to download the latest release, but if there is a problem with\
                      the latest release, feel free to run setup using a different URL. Just to note, anvi'o\
                      will expect to find the following files in the URL provided here: 'VERSION', \
                      'ar122_msa_individual_genes.tar.gz', 'ar122_taxonomy.tsv', 'bac120_msa_individual_genes.tar.gz', \
                      and 'bac120_taxonomy.tsv'. If everything fails, you can give this URL, which is supposed to work\
                      if teh server in which these databases are maintained is still online:\
                      https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/"}
                ),
    'reset': (
            ['--reset'],
            {'default': False,
             'action': 'store_true',
             'help': "Remove all the previously stored files and start over. If something is feels wrong\
                      for some reason and if you believe re-downloading files and setting them up could\
                      address the issue, this is the flag that will tell anvi'o to act like a real comptuer\
                      scientist challenged with a computational problem."}
                ),
    'cog-data-dir': (
            ['--cog-data-dir'],
            {'default': None,
             'type': str,
             'help': "The directory path for your COG setup. Anvi'o will try to use the default path\
                      if you do not specify anything."}
                ),
    'pfam-data-dir': (
            ['--pfam-data-dir'],
            {'default': None,
             'type': str,
             'help': "The directory path for your Pfam setup. Anvi'o will try to use the default path\
                      if you do not specify anything."}
                ),
    'hide-outlier-SNVs': (
            ['--hide-outlier-SNVs'],
            {'default': False,
             'action': 'store_true',
             'help': "During profiling, anvi'o marks positions of single-nucleotide variations (SNVs)\
                      that originate from places in contigs where coverage values are a bit 'sketchy'.\
                      If you would like to avoid SNVs in those positions of splits in applicable projects\
                      you can use this flag, and the interface would hide SNVs that are marked as 'outlier'\
                      (although it is clearly the best to see everything, no one will judge you if you end\
                      up using this flag) (plus, there may or may not be some historical data on this here: \
                      https://github.com/meren/anvio/issues/309)."}
                ),
    'hmm-source': (
            ['--hmm-source'],
            {'metavar': 'SOURCE NAME',
             'default': None,
             'help': "Use a specific HMM source. You can use '--list-hmm-sources' flag to see\
                      a list of available resources. The default is '%(default)s'."}
                ),
    'hmm-sources': (
            ['--hmm-sources'],
            {'metavar': 'SOURCE NAME',
             'help': "Get sequences for a specific list of HMM sources. You can list one or more\
                      sources by separating them from each other with a comma character (i.e., \
                      '--hmm-sources source_1,source_2,source_3'). If you would like to see a list\
                      of available sources in the contigs database, run this program with\
                      '--list-hmm-sources' flag."}
                ),
    'list-hmm-sources': (
            ['-l', '--list-hmm-sources'],
            {'default': False,
             'action': 'store_true',
             'help': "List available HMM sources in the contigs database and quit."}
                ),
    'annotation-source': (
            ['--annotation-source'],
            {'metavar': 'SOURCE NAME',
             'default': None,
             'help': "Get functional annotations for a specific annotation source. You can use the flag\
                      '--list-annotation-sources' to learn about what sources are available."}
                ),
    'annotation-sources': (
            ['--annotation-sources'],
            {'metavar': 'SOURCE NAME[S]',
             'default': None,
             'help': "Get functional annotations for a specific list of annotation sources. You\
                      can specify one or more sources by separating them from each other with a comma\
                      character (i.e., '--annotation-sources source_1,source_2,source_3'). The default\
                      behavior is to return everything"}
                ),
    'list-annotation-sources': (
            ['-l', '--list-annotation-sources'],
            {'default': False,
             'action': 'store_true',
             'help': "List available functional annotation sources."}
                ),
    'include-gc-identity-as-function': (
            ['--include-gc-identity-as-function'],
            {'default': False,
             'action': 'store_true',
             'help': "This is an option that asks anvi'o to treat gene cluster names as functions. By\
                      doing so, you are in fact creating an opportunity to study functional enrichment\
                      statistics for each gene cluster independently. For instance, multiple gene\
                      clusters may have the same COG function. But if you wish to use the same enrichment\
                      analysis in your pangenome without collapsing multiple gene clusters into a single\
                      function name, you can use this flag, and ask for 'IDENTITY' as the functional\
                      annotation source."}
                ),
    'gene-names': (
            ['--gene-names'],
            {'metavar': 'HMM HIT NAME',
             'help': "Get sequences only for a specific gene name. Each name should be separated from\
                      each other by a comma character. For instance, if you want to get back only RecA\
                      and Ribosomal_L27, you can type '--gene-names RecA,Ribosomal_L27', and you will\
                      get any and every hit that matches these names in any source. If you would like\
                      to see a list of available gene names, you can use '--list-available-gene-names'\
                      flag."}
                ),
    'get-aa-sequences': (
            ['--get-aa-sequences'],
            {'default': False,
             'action': 'store_true',
             'help': "Store amino acid sequences instead."}
                ),
    'return-best-hit': (
            ['--return-best-hit'],
            {'default': False,
             'action': 'store_true',
             'help': "A bin may contain more than one hit for a gene name in a given HMM source. For instance, there may\
                      be multiple RecA hits in a genome bin from Campbell et al.. Using this flag, will go through all of\
                      the gene names that appear multiple times, and remove all but the one with the lowest e-value. Good\
                      for whenever you really need to get only a single copy of single-copy core genes from a genome bin."}
                ),
    'max-num-genes-missing-from-bin': (
            ['--max-num-genes-missing-from-bin'],
            {'default': None,
             'metavar': 'INTEGER',
             'help': "This filter removes bins (or genomes) from your analysis. If you have a list of gene names, you can\
                      use this parameter to omit any bin (or external genome) that is missing more than a number of genes\
                      you desire. For instance, if you have 100 genome bins, and you are interested in working with 5\
                      ribosomal proteins, you can use '--max-num-genes-missing-from-bin 4' to remove the bins that\
                      are missing more than 4 of those 5 genes. This is especially useful for phylogenomic analyses.\
                      Parameter 0 will remove any bin that is missing any of the genes."}
                ),
    'min-num-bins-gene-occurs': (
            ['--min-num-bins-gene-occurs'],
            {'default': None,
             'metavar': 'INTEGER',
             'help': "This filter removes genes from your analysis. Let's assume you have 100 bins to get sequences for HMM\
                      hits. If you want to work only with genes among all the hits that occur in at least X number of bins,\
                      and discard the rest of them, you can use this flag. If you say '--min-num-bins-gene-occurs 90', each\
                      gene in the analysis will be required at least to appear in 90 genomes. If a gene occurs in less than\
                      that number of genomes, it simply will not be reported. This is especially useful for phylogenomic\
                      analyses, where you may want to only focus on genes that are prevalent across the set of genomes\
                      you wish to analyze."}
                ),
    'max-num-gene-clusters-missing-from-genome': (
            ['--max-num-gene-clusters-missing-from-genome'],
            {'default': 0,
             'metavar': 'INTEGER',
             'help': "This filter will remove genomes from your report. If you have a list of gene cluster names, you can\
                      use this parameter to omit any genome from your report if it is missing more than a number of genes\
                      you desire. For instance, if you have 100 genomes in your pan genome, and you are interested in\
                      working only with genomes that have all 5 specific gene clusters of your choice, you can use \
                      '--max-num-gene-clusters-missing-from-genome 4' to remove remove the bins that\
                      are missing more than 4 of those 5 genes. This is especially useful for phylogenomic analyses.\
                      Parameter 0 will remove any genome that is missing any of the genes."}
                ),
    'min-num-genomes-gene-cluster-occurs': (
            ['--min-num-genomes-gene-cluster-occurs'],
            {'default': 0,
             'metavar': 'INTEGER',
             'help': "This filter will remove gene clusters from your report. Let's assume you have 100 genomes in your pan\
                      genome analysis. You can use this parameter if you want to work only with gene clusters that occur in\
                      at least X number of genomes. If you say '--min-num-genomes-gene-cluster-occurs 90', each\
                      gene cluster in the analysis will be required at least to appear in 90 genomes. If a gene occurs in\
                      less than that number of genomes, it simply will not be reported. This is especially useful for\
                      phylogenomic analyses, where you may want to only focus on gene clusters that are prevalent across\
                      the set of genomes you wish to analyze."}
                ),
    'max-num-genomes-gene-cluster-occurs': (
            ['--max-num-genomes-gene-cluster-occurs'],
            {'default': sys.maxsize,
             'metavar': 'INTEGER',
             'help': "This filter will remove gene clusters from your report. Let's assume you have 100 genomes in your pan\
                      genome analysis. You can use this parameter if you want to work only with gene clusters that occur in\
                      at most X number of genomes. If you say '--max-num-genomes-gene-cluster-occurs 1', you will get gene\
                      clusters that are singletons. Combining this parameter with --min-num-genomes-gene-cluster-occurs can\
                      give you a very precise way to filter your gene clusters."}
                ),
    'min-num-genes-from-each-genome': (
            ['--min-num-genes-from-each-genome'],
            {'default': 0,
             'metavar': 'INTEGER',
             'help': "This filter will remove gene clusters from your report. If you say '--min-num-genes-from-each-genome 2',\
                      this filter will remove every gene cluster, to which every genome in your analysis contributed less than\
                      2 genes. This can be useful to find out gene clusters with many genes from many genomes (such as conserved\
                      multi-copy genes within a clade)."}
                ),
    'max-num-genes-from-each-genome': (
            ['--max-num-genes-from-each-genome'],
            {'default': sys.maxsize,
             'metavar': 'INTEGER',
             'help': "This filter will remove gene clusters from your report. If you say '--max-num-genes-from-each-genome 1',\
                      every gene cluster that has more than one gene from any genome that contributes to it will be removed\
                      from your analysis. This could be useful to remove gene clusters with paralogs from your report for\
                      appropriate phylogenomic analyses. For instance, using '--max-num-genes-from-each-genome 1' and \
                      'min-num-genomes-gene-cluster-occurs X' where X is the total number of your genomes, would give you the\
                      single-copy gene clusters in your pan genome."}
                ),
    'min-functional-homogeneity-index': (
            ['--min-functional-homogeneity-index'],
            {'default': -1,
             'metavar': 'FLOAT',
             'type': float,
             'help': "This filter will remove gene clusters from your report. If you say '--min-functional-homogeneity-index 0.3', \
                      every gene cluster with a functional homogeneity index less than 0.3 will be removed from your analysis. This \
                      can be useful if you only want to look at gene clusters that are highly conserved in resulting function"}
                ),
    'max-functional-homogeneity-index': (
            ['--max-functional-homogeneity-index'],
            {'default': 1,
             'metavar': 'FLOAT',
             'type': float,
             'help': "This filter will remove gene clusters from your report. If you say '--max-functional-homogeneity-index 0.5', \
                      every gene cluster with a functional homogeneity index greater than 0.5 will be removed from your analysis. This \
                      can be useful if you only want to look at gene clusters that don't seem to be functionally conserved"}
                ),
    'min-geometric-homogeneity-index': (
            ['--min-geometric-homogeneity-index'],
            {'default': -1,
             'metavar': 'FLOAT',
             'type': float,
             'help': "This filter will remove gene clusters from your report. If you say '--min-geometric-homogeneity-index 0.3', \
                      every gene cluster with a geometric homogeneity index less than 0.3 will be removed from your analysis. This \
                      can be useful if you only want to look at gene clusters that are highly conserved in geometric configuration"}
                ),
    'max-geometric-homogeneity-index': (
            ['--max-geometric-homogeneity-index'],
            {'default': 1,
             'metavar': 'FLOAT',
             'type': float,
             'help': "This filter will remove gene clusters from your report. If you say '--max-geometric-homogeneity-index 0.5', \
                      every gene cluster with a geometric homogeneity index greater than 0.5 will be removed from your analysis. This \
                      can be useful if you only want to look at gene clusters that have many not be as conserved as others"}
                ),
    'min-combined-homogeneity-index': (
            ['--min-combined-homogeneity-index'],
            {'default': -1,
             'metavar': 'FLOAT',
             'type': float,
             'help': "This filter will remove gene clusters from your report. If you say '--min-combined-homogeneity-index 0.3', \
                      every gene cluster with a combined homogeneity index less than 0.3 will be removed from your analysis. This \
                      can be useful if you only want to look at gene clusters that are highly conserved overall"}
                ),
    'max-combined-homogeneity-index': (
            ['--max-combined-homogeneity-index'],
            {'default': 1,
             'metavar': 'FLOAT',
             'type': float,
             'help': "This filter will remove gene clusters from your report. If you say '--max-combined-homogeneity-index 0.5', \
                      every gene cluster with a combined homogeneity index greater than 0.5 will be removed from your analysis. This \
                      can be useful if you only want to look at gene clusters that have many not be as conserved overall as others"}
                ),
    'add-into-items-additional-data-table': (
            ['--add-into-items-additional-data-table'],
            {'default': None,
             'metavar': 'NAME',
             'help': "If you use any of the filters, and would like to add the resulting item names into the items additional\
                      data table of your database, you can use this parameter. You will need to give a name for these results to\
                      be saved. If the given name is already in the items additional data table, its contents will be replaced\
                      with the new one. Then you can run anvi-interactive or anvi-display-pan to 'see' the results of your filters."}
                ),
    'concatenate-genes': (
            ['--concatenate-genes'],
            {'default': False,
             'action': 'store_true',
             'help': "Concatenate output genes in the same order to create a multi-gene alignment output that is suitable\
                      for phylogenomic analyses."}
                ),
    'separator': (
            ['--separator'],
            {'metavar': 'STRING',
             'default': None,
             'type': str,
             'help': "Characters to separate things (the default is whatever is most suitable)."}
                ),
    'align-with': (
            ['--align-with'],
            {'metavar': 'ALIGNER',
             'default': None,
             'type': str,
             'help': "The multiple sequence alignment program to use when multiple sequence alignment is necessary. To see\
                      all available options, use the flag `--list-aligners`."}
                ),
    'list-aligners': (
            ['--list-aligners'],
            {'default': False,
             'action': 'store_true',
             'help': "Show available software for multiple sequence alignment."}
                ),
    'concatenate-gene-clusters': (
            ['--concatenate-gene-clusters'],
            {'default': False,
             'action': 'store_true',
             'help': "Concatenate output gene clusters in the same order to create a multi-gene alignment output that is suitable\
                      for phylogenomic analyses."}
                ),
    'report-DNA-sequences': (
            ['--report-DNA-sequences'],
            {'default': False,
             'action': 'store_true',
             'help': "By default, this program reports amino acid sequences. You can change that behavior and as for DNA\
                      sequences instead using this flag."}
                ),
    'skip-multiple-gene-calls': (
            ['--skip-multiple-gene-calls'],
            {'default': False,
             'action': 'store_true',
             'help': "When generating concatenated output skip gene clusters contain multiple gene calls."}
                ),
    'list-available-gene-names': (
            ['-L', '--list-available-gene-names'],
            {'default': False,
             'action': 'store_true',
             'help': "List available gene names in HMM sources selection and quit."}
                ),
    'search-terms': (
            ['--search-terms'],
            {'metavar': 'SEARCH_TERMS',
             'help': "Search terms. Multiple of them can be declared separated by a delimiter (the default is a comma)."}
                ),
    'sensitive': (
            ['--sensitive'],
            {'default': False,
             'action': 'store_true',
             'help': "DIAMOND sensitivity. With this flag you can instruct DIAMOND to be 'sensitive', rather than 'fast'\
                      during the search. It is likely the search will take remarkably longer. But, hey, if you are doing\
                      it for your final analysis, maybe it should take longer and be more accurate. This flag is only\
                      relevant if you are running DIAMOND."}
                ),
    'gene-caller-ids': (
            ['--gene-caller-ids'],
            {'metavar': 'GENE_CALLER_IDS',
             'type': str,
             'help': "Gene caller ids. Multiple of them can be declared separated by a delimiter (the default is a comma).\
                      In anvi-gen-variability-profile, if you declare nothing you will get all genes matching your other\
                      filtering criteria. In other programs, you may get everything, nothing, or an error. It really depends\
                      on the situation. Fortunately, mistakes are cheap, so it's worth a try."}
                ),
    'flank-mode': (
            ['--flank-mode'],
            {'action': 'store_true',
             'help': "If in --flank-mode, anvi-export-locus will extract a locus based on the coordinates \
                     of flanking genes. You MUST provide 2 flanking genes in the form of TWO \
                     --search-term, --gene-caller-ids, or --hmm-sources. The --flank-mode option is  \
                     appropriate for extracting loci of variable gene number lengths, but are consistantly \
                     located between the same flanking genes in the genome(s) of interest."}
              ),
    'num-genes': (
            ['-n','--num-genes'],
            {'metavar': 'NUM_GENES',
             'type': str,
             'help': "Required for DEFAULT mode. For each match (to the function, or HMM that was searched) a sequence which includes \
                      a block of genes will be saved. The block could include either genes only in the forward direction of the gene (defined \
                      according to the direction of transcription of the gene) or reverse or both. \
                      If you wish to get both direction use a comma (no spaces) to define the block \
                      For example, \"-n 4,5\" will give you four genes before and five genes after. \
                      Whereas, \"-n 5\" will give you five genes after (in addition to the gene that matched). \
                      To get only genes preceeding the match use \"-n 5,0\". \
                      If the number of genes requested exceeds the length of the contig, then the output \
                      will include the sequence until the end of the contig."}
              ),
    'gene-mode': (
            ['--gene-mode'],
            {'default': False,
             'action': 'store_true',
             'help': "Initiate the interactive interface in \"gene mode\". In this mode, the items are genes (instead of\
                      splits of contigs). The following views are available: detection (the detection value of each gene\
                      in each sample). The mean_coverage (the mean coverage of genes). The non_outlier_mean_coverage\
                      (the mean coverage of the non-outlier nucleotide positions of each gene in each sample (median absolute\
                      deviation is used to remove outliers per gene per sample)). The non_outlier_coverage_std view (standard deviation\
                      of the coverage of non-outlier positions of genes in samples). You can also choose to order items\
                      and layers according to each one of the aforementioned views. In addition, all layer ordering\
                      that are available in the regular mode (i.e. the full mode where you have contigs/splits) are also\
                      available in \"gene mode\", so that, for example, you can choose to order the layers according to \"detection\", and that\
                      would be the order according to the detection values of splits, whereas if you choose \"genes_detections\"\
                      then the order of layers would be according to the detection values of genes. Inspection and sequence\
                      functionality are available (through the right-click menu), except now sequences are of the specific gene.\
                      Inspection has now two options available: \"Inspect Context\", which brings you to the inspection page of the split\
                      to which the gene belongs where the inspected gene will be highlighted in yellow in the bottom, and \"Inspect Gene\",\
                      which opens the inspection page only for the gene and 100 nts around each side of it (the purpose of this option\
                      is to make the inspection page load faster if you only want to look at the nucleotide coverage of a specific gene).\
                      NOTICE: You can't store states or collections in \"gene mode\". However, you still can make fake selections, and create\
                      fake bins for your viewing convenience only (smiley). Search options are available, and you can even search for functions\
                      if you have them in your contigs database. ANOTHER NOTICE: loading this mode might take a while if your bin\
                      has many genes, and your profile database has many samples, this is because the gene coverages stats are\
                      computed in an ad-hoc manner when you load this mode, we know this is not ideal and we plan to improve that\
                      (along with other things). If you have suggestions/complaints regarding this mode please comment on this\
                      github issue: https://goo.gl/yHhRei. Please refer to the online tutorial for more information."}
                ),
    'gene-caller-id': (
            ['--gene-caller-id'],
            {'metavar': 'GENE_CALLER_ID',
             'type': int,
             'help': "A single gene id."}
                ),
    'target-version': (
            ['-t', '--target-version'],
            {'metavar': 'VERSION',
             'type': int,
             'help': "Anvi'o will stop upgrading your database when it reaches to this version. "}
                ),
    'delimiter': (
            ['--delimiter'],
            {'metavar': 'CHAR',
             'default': ',',
             'help': "The delimiter to parse multiple input terms. The default is '%(default)s'."}
                ),
    'wrap': (
            ['--wrap'],
            {'metavar': 'WRAP',
             'default': 120,
             'type': int,
             'help': "When to wrap sequences when storing them in a FASTA file. The default is\
                      '%(default)d'. A value of '0' would be equivalent to 'do not wrap'."}
                ),
    'no-wrap': (
            ['--no-wrap'],
            {'default': False,
             'action': 'store_true',
             'help': "Do not be wrap sequences nicely in the output file."}
                ),
    'leeway': (
            ['--leeway'],
            {'metavar': 'LEEWAY_NTs',
             'default': 100,
             'type': int,
             'help': "The minimum number of nucleotides for a given short read mapping into\
                      the gene context for it to be reported. You must consider the length of\
                      your short reads, as well as the length of the gene you are targeting.\
                      The default is %(default)d nts."}
                ),
    'split-R1-and-R2': (
            ['-Q', '--split-R1-and-R2'],
            {'default': False,
             'action': 'store_true',
             'help': "When declared, this program outputs 3 FASTA files for paired-end reads: one\
                      for R1, one for R2, and one for unpaired reads."}
                ),
    'gzip-output': (
            ['-X', '--gzip-output'],
            {'default': False,
             'action': 'store_true',
             'help': "When declared, output file(s) will be gzip compressed and the extension `.gz` will be added."}
                ),
    'list-contigs': (
            ['--list-contigs'],
            {'default': False,
             'action': 'store_true',
             'help': "When declared, the program will list contigs in the BAM file and exit gracefully\
                      without any further analysis."}
                ),
    'list-splits': (
            ['--list-splits'],
            {'default': False,
             'action': 'store_true',
             'help': "When declared, the program will list split names in the profile database and quite"}
                ),

    'list-collections': (
            ['--list-collections'],
            {'default': False,
             'action': 'store_true',
             'help': "Show available collections and exit."}
                ),
    'list-bins': (
            ['--list-bins'],
            {'default': False,
             'action': 'store_true',
             'help': "List available bins in a collection and exit."}
                ),
    'list-states': (
            ['--list-states'],
            {'default': False,
             'action': 'store_true',
             'help': "Show available states and exit."}
                ),
    'show-views': (
            ['--show-views'],
            {'default': False,
             'action': 'store_true',
             'help': "When declared, the program will show a list of available views, and exit."}
                ),
    'list-completeness-sources': (
            ['--list-completeness-sources'],
            {'default': False,
             'action': 'store_true',
             'help': "Show available sources and exit."}
                ),
    'completeness-source': (
            ['--completeness-source'],
            {'metavar': 'NAME',
             'help': "Single-copy gene source to use to estimate completeness."}
                ),
    'split-name': (
            ['--split-name'],
            {'metavar': 'SPLIT_NAME',
             'help': "Split name."}
                ),
    'program': (
            ['--program'],
            {'metavar': 'PROGRAM_NAME',
             'help': "Program name.",
             'required': False,
             'default': 'default'}
                ),
    'splits-of-interest': (
            ['--splits-of-interest'],
            {'metavar': 'FILE',
             'help': "A file with split names. There should be only one column in the file, and each line\
                      should correspond to a unique split name."}
                ),
    'contigs-of-interest': (
            ['--contigs-of-interest'],
            {'metavar': 'FILE',
             'help': "It is possible to analyze only a group of contigs from a given BAM file. If you provide\
                      a text file, in which every contig of interest is listed line by line, the profiler would\
                      engine only on those contigs in the BAM file and ignore the rest. This can be used for\
                      debugging purposes, or to engine on a particular group of contigs that were identified as\
                      relevant during the interactive analysis."}
                ),
    'samples-of-interest': (
            ['--samples-of-interest'],
            {'metavar': 'FILE',
             'help': "A file with samples names. There should be only one column in the file, and each line\
                      should correspond to a unique sample name (without a column header)."}
                ),
    'genes-of-interest': (
            ['--genes-of-interest'],
            {'metavar': 'FILE',
             'help': "A file with anvi'o gene caller IDs. There should be only one column in the file, and each line\
                      should correspond to a unique gene caller id (without a column header)."}
                ),
    'gene-cluster-id': (
            ['--gene-cluster-id'],
            {'metavar': 'GENE_CLUSTER_ID',
             'help': "Gene cluster ID you are interested in."}
                ),
    'gene-cluster-ids-file': (
            ['--gene-cluster-ids-file'],
            {'metavar': 'FILE_PATH',
             'help': "Text file for gene clusters (each line should contain be a unique gene cluster id)."}
                ),
    'bin-id': (
            ['-b', '--bin-id'],
            {'metavar': 'BIN_NAME',
             'help': "Bin name you are interested in."}
                ),
    'bin-names-list': (
            ['-b', '--bin-names-list'],
            {'metavar': 'BIN NAMES',
             'help': "Comma-separated list of bin names."}
                ),
    'new-bin-name': (
            ['-B', '--new-bin-name'],
            {'metavar': 'BIN NAME',
             'help': "The new bin name."}
                ),
    'bin-ids-file': (
            ['-B', '--bin-ids-file'],
            {'metavar': 'FILE_PATH',
             'help': "Text file for bins (each line should be a unique bin id)."}
                ),
    'collection-name': (
            ['-C', '--collection-name'],
            {'metavar': 'COLLECTION_NAME',
             'help': "Collection name."}
                ),
    'num-positions-from-each-split': (
            ['--num-positions-from-each-split'],
            {'metavar': 'INT',
             'default': 0,
             'type': int,
             'help': "Each split may have one or more variable positions. By default, anvi'o will report every SNV\
                      position found in a given split. This parameter will help you to define a cutoff for the maximum\
                      number of SNVs to be reported from a split (if the number of SNVs is more than the number you\
                      declare using this parameter, the positions will be randomly subsampled)."}
             ),
    'min-scatter': (
            ['-m', '--min-scatter'],
            {'metavar': 'INT',
             'default': 0,
             'type': int,
             'help': "This one is tricky. If you have N samples in your dataset, a given variable position x in one\
                      of your splits can split your N samples into `t` groups based on the identity of the\
                      variation they harbor at position x. For instance, `t` would have been 1, if all samples had the same\
                      type of variation at position x (which would not be very interesting, because in this case\
                      position x would have zero contribution to a deeper understanding of how these samples differ\
                      based on variability. When `t` > 1, it would mean that identities at position x across samples\
                      do differ. But how much scattering occurs based on position x when t > 1? If t=2, how many\
                      samples ended in each group? Obviously, even distribution of samples across groups may tell\
                      us something different than uneven distribution of samples across groups. So, this parameter\
                      filters out any x if 'the number of samples in the second largest group' (=scatter) is less\
                      than -m. Here is an example: let's assume you have 7 samples. While 5 of those have AG, 2\
                      of them have TC at position x. This would mean scatter of x is 2. If you set -m to 2, this\
                      position would not be reported in your output matrix. The default value for -m is\
                      %(default)d, which means every `x` found in the database and survived previous filtering\
                      criteria will be reported. Naturally, -m cannot be more than half of the number of samples.\
                      Please refer to the user documentation if this is confusing."}
                ),
    'min-ratio-of-competings-nts': (
            ['-r', '--min-ratio-of-competings-nts'],
            {'metavar': 'RATIO',
             'default': 0,
             'type': float,
             'help': "Minimum ratio of the competing nucleotides at a given position. Default is %(default)d."}
                ),
    'max-num-unique-positions': (
            ['-n', '--max-num-unique-positions'],
            {'metavar': 'NUM_POSITIONS',
             'default': 0,
             'type': int,
             'help': "Maximum number of unique positions to be used in the network. This may be one way to avoid extremely\
                      large network descriptions that would defeat the purpose of a quick visualization. If there are more\
                      unique positions in the variability profile, the program will randomly select a subset of them to match\
                      the `max-num-unique-positions`. The default is %(default)d, which means all positions should be reported.\
                      Remember that the number of nodes in the network will also depend on the number of samples described in\
                      the variability profile."}
                ),
    'num-threads': (
            ['-T', '--num-threads'],
            {'metavar': 'NUM_THREADS',
             'default': 1,
             'type': int,
             'help': "Maximum number of threads to use for multithreading whenever possible. Very conservatively, the default\
                      is 1. It is a good idea to not exceed the number of CPUs / cores on your system. Plus, please\
                      be careful with this option if you are running your commands on a SGE --if you are clusterizing your runs,\
                      and asking for multiple threads to use, you may deplete your resources very fast."}
                ),
    'num-parallel-processes': (
            ['-P', '--num-parallel-processes'],
            {'metavar': 'NUM_PROCESSES',
             'default': 1,
             'type': int,
             'help': "Maximum number of processes to run in parallel. Please note that this is different than number of threads. If you\
                      ask for 4 parallel processes, and 5 threads, anvi'o will run four processes in parallel and assign 5 threads\
                      to each. For resource allocation you must multiply the number of processes and threads."}
                ),
    'variability-profile': (
            ['-V', '--variability-profile'],
            {'metavar': 'VARIABILITY_TABLE',
             'type': str,
             'required': False,
             'help': "FIXME"}
                ),
    'min-coverage-in-each-sample': (
            ['--min-coverage-in-each-sample'],
            {'metavar': 'INT',
             'default': 0,
             'type': int,
             'help': "Minimum coverage of a given variable nucleotide position in all samples. If a nucleotide position\
                      is covered less than this value even in one sample, it will be removed from the analysis. Default\
                      is %(default)d."}
                ),
    'min-departure-from-reference': (
            ['-r', '--min-departure-from-reference'],
            {'metavar': 'FLOAT',
             'default': 0,
             'type': float,
             'help': "Takes a value between 0 and 1, where 1 is maximum divergence from the reference. Default is %(default)f.\
                      The reference here observation that corresponds to a given position in the mapped context."}
                ),
    'max-departure-from-reference': (
            ['-z', '--max-departure-from-reference'],
            {'metavar': 'FLOAT',
             'default': 1,
             'type': float,
             'help': "Similar to '--min-departure-from-reference', but defines an upper limit for divergence. The\
                      default is %(default)f."}
                ),
    'min-departure-from-consensus': (
            ['-j', '--min-departure-from-consensus'],
            {'metavar': 'FLOAT',
             'default': 0,
             'type': float,
             'help': "Takes a value between 0 and 1, where 1 is maximum divergence from the consensus for a given position. The\
                      default is %(default)f. The consensus is the most frequent observation at a given positon."}
                ),
    'max-departure-from-consensus': (
            ['-a', '--max-departure-from-consensus'],
            {'metavar': 'FLOAT',
             'default': 1,
             'type': float,
             'help': "Similar to '--min-departure-from-consensus', but defines an upper limit for divergence. The\
                      default is %(default)f."}
                ),
    'min-occurrence-of-variable-positions': (
            ['-x', '--min-occurrence'],
            {'metavar': 'NUM_SAMPLES',
             'default': 1,
             'type': int,
             'help': "Minimum number of samples a nucleotide position should be reported as variable. Default is %(default)d.\
                      If you set it to 2, for instance, each eligible variable position will be expected to appear in at least\
                      two samples, which will reduce the impact of stochastic, or unintelligible variable positions."}
                ),
    'quince-mode': (
            ['--quince-mode'],
            {'default': False,
             'action': 'store_true',
             'help': "The default behavior is to report base frequencies of nucleotide positions only if there\
                      is any variation reported during profiling (which by default uses some heuristics to minimize\
                      the impact of error-driven variation). So, if there are 10 samples, and a given position has been\
                      reported as a variable site during profiling in only one of those samples, there will be no\
                      information will be stored in the database for the remaining 9. When this flag is\
                      used, we go back to each sample, and report base frequencies for each sample at this position\
                      even if they do not vary. It will take considerably longer to report when this flag is on, and the use\
                      of it will increase the file size dramatically, however it is inevitable for some statistical approaches\
                      (as well as for some beautiful visualizations)."}
                ),
    'include-contig-names': (
            ['--include-contig-names'],
            {'default': False,
             'action': 'store_true',
             'help': "Use this flag if you would like contig names for each variable position to be included in the\
                      output file as a column. By default, we do not include contig names since they can practically\
                      double the output file size without any actual benefit in most cases."}
                ),
    'include-split-names': (
            ['--include-split-names'],
            {'default': False,
             'action': 'store_true',
             'help': "Use this flag if you would like split names for each variable position to be included in the\
                      output file as a column."}
                ),
    'engine': (
            ['--engine'],
            {'default': 'NT',
             'metavar': 'ENGINE',
             'type': str,
             'help': "Varaibility engine. The default is '%(default)s'."}
                ),
    'driver': (
            ['--driver'],
            {'metavar': 'DRIVER',
             'type': str,
             'required': True,
             'help': "Automatic binning drivers. Available options '%(choices)s'."}
                ),
    'skip-synonymity': (
            ['--skip-synonymity'],
            {'default': False,
             'action': 'store_true',
             'help': "Computing synonymity can be an expensive operation for large data sets. Provide this flag to skip\
                      computing synonymity. It only makes sense to provide this flag when using --engine CDN."}
                ),
    'transpose': (
            ['--transpose'],
            {'default': False,
             'action': 'store_true',
             'help': "Transpose the input matrix file before clustering."}
                ),
    'skip-check-names': (
            ['--skip-check-names'],
            {'default': False,
             'action': 'store_true',
             'help': "For debugging purposes. You should never really need it."}
                ),
    'experimental-org-input-dir': (
            ['-i', '--input-directory'],
            {'metavar': 'DIR_PATH',
             'type': str,
             'help': "Input directory where the input files addressed from the configuration\
                      file can be found (i.e., the profile database, if PROFILE.db::TABLE\
                      notation is used in the configuration file)."}
                ),
    'clustering-name': (
            ['-N', '--name'],
            {'metavar': 'NAME',
             'type': str,
             'help': "The name to use when storing the resulting clustering in the database.\
                      This name will appear in the interactive interface and other relevant\
                      interfaces. Please consider using a short and descriptive single-word\
                      (if you do not do that you will make anvi'o complain)."}
                ),
    'distance': (
            ['--distance'],
            {'metavar': 'DISTANCE_METRIC',
             'type': str,
             'default': constants.distance_metric_default,
             'help': "The distance metric for the hierarchical clustering. The default distance\
                      metric is '%(default)s'. You can find the full list of distance metrics\
                      either by making a mistake (such as entering a non-existent distance metric\
                      and making anvi'o upset), or by taking a look at the help menu of the\
                      hierarchy.distance.pdist function in the scipy.cluster module."}
                ),
    'linkage': (
            ['--linkage'],
            {'metavar': 'LINKAGE_METHOD',
             'type': str,
             'default': constants.linkage_method_default,
             'help': "The linkage method for the hierarchical clustering. The default linkage\
                      method is '%(default)s', because that is the best one. It really is. We talked\
                      to a lot of people and they were all like 'this is the best one available' and\
                      it is just all out there. Honestly it is so good that we will build a wall around it\
                      and make other linkage methods pay for it. But if you want to see a full\
                      list of available ones you can check the hierarcy.linkage function in\
                      the scipy.cluster module. Up to you really. But then you can't use %(default)s\
                      anymore, and you would have to leave anvi'o right now."}
                ),
    'input-dir': (
            ['-i', '--input-dir'],
            {'metavar': 'DIR_PATH',
             'type': str,
             'help': "Directory path for input files"}
                ),
    'output-dir': (
            ['-o', '--output-dir'],
            {'metavar': 'DIR_PATH',
             'type': str,
             'help': "Directory path for output files"}
                ),
    'output-file': (
            ['-o', '--output-file'],
            {'metavar': 'FILE_PATH',
             'type': str,
             'help': "File path to store results."}
                ),
    'log-file': (
            ['--log-file'],
            {'metavar': 'FILE_PATH',
             'default': None,
             'type': str,
             'help': "File path to store debug/output messages."}
                ),
    'output-db-path': (
            ['-o', '--output-db-path'],
            {'metavar': 'DB_FILE_PATH',
             'type': str,
             'help': "Output file path for the new database."}
                ),
    'temporary-dir-path': (
            ['--temporary-dir-path'],
            {'metavar': 'PATH',
             'type': str,
             'help': "If you don't provide anything here, this program will come up with a temporary\
                      directory path by itself to store intermediate files, and clean it later. If you\
                      want to have full control over this, you can use this flag to define one.."}
                ),
    'output-file-prefix': (
            ['-O', '--output-file-prefix'],
            {'metavar': 'FILENAME_PREFIX',
             'type': str,
             'help': "A prefix to be used while naming the output files (no file type\
                      extensions please; just a prefix)."}
                ),
    'dry-run': (
            ['--dry-run'],
            {'default': False,
             'action': 'store_true',
             'help': "Don't do anything real. Test everything, and stop right before wherever the developer\
                      said 'well, this is enough testing', and decided to print out results."}
                ),
    'verbose': (
            ['--verbose'],
            {'default': False,
             'action': 'store_true',
             'help': "Be verbose, print more messages whenever possible."}
                ),
    'concise': (
            ['--concise'],
            {'default': False,
             'action': 'store_true',
             'help': "Don't be verbose, print less messages whenever possible."}
                ),
    'just-do-it': (
            ['--just-do-it'],
            {'default': False,
             'action': 'store_true',
             'help': "Don't bother me with questions or warnings, just do it."}
                ),
    'ip-address': (
            ['-I', '--ip-address'],
            {'metavar': 'IP_ADDR',
             'type': str,
             'default': '0.0.0.0',
             'help': "IP address for the HTTP server. The default ip address (%(default)s) should\
                      work just fine for most."}
                ),
   'browser-path': (
            ['--browser-path'],
            {'metavar': 'PATH',
             'type': str,
             'default': None,
             'help': "By default, anvi'o will use your default browser to launch the interactive interface. If you\
                      would like to use something else than your system default, you can provide a full path for an\
                      alternative browser using this parameter, and hope for the best. For instance we are using\
                      this parameter to call Google's experimental browser, Canary, which performs better with\
                      demanding visualizations."}
                ),
   'api-url': (
            ['--api-url'],
            {'metavar': 'API_URL',
             'type': str,
             'default': 'https://anvi-server.org',
             'help': "Anvi'server url"}
                ),
    'port-number': (
            ['-P', '--port-number'],
            {'metavar': 'INT',
             'default': None,
             'type': int,
             'help': "Port number to use for anvi'o services. If nothing is declared, anvi'o will try to find\
                      a suitable port number, starting from the default port number, %d." % constants.default_port_number}
                ),
    'user': (
            ['--user'],
            {'metavar': 'USERNAME',
             'default': None,
             'type': str,
             'help': "The user for an anvi'server."}
                ),
    'user-server-shutdown': (
            ['--user-server-shutdown'],
            {'default': False,
             'action': 'store_true',
             'help': "Allow users to shutdown an anvi'server via web interface."}
                ),
    'read-only': (
            ['--read-only'],
            {'default': False,
             'action': 'store_true',
             'help': "When the interactive interface is started with this flag, all 'database write'\
                      operations will be disabled."}
                ),
    'server-only': (
            ['--server-only'],
            {'default': False,
             'action': 'store_true',
             'help': "The default behavior is to start the local server, and fire up a browser that\
                      connects to the server. If you have other plans, and want to start the server\
                      without calling the browser, this is the flag you need."}
                ),
    'password-protected': (
            ['--password-protected'],
            {'default': False,
             'action': 'store_true',
             'help': "If this flag is set, command line tool will ask you to enter a password and interactive \
                      interface will be only accessible after entering same password. This option is recommended \
                      for shared machines like clusters or shared networks where computers are not isolated."}
                ),
    'store-in-db': (
            ['--store-in-db'],
            {'default': False,
             'action': 'store_true',
             'help': "Store analysis results into the database directly."}
                ),
    'skip-store-in-db': (
            ['--skip-store-in-db'],
            {'default': False,
             'action': 'store_true',
             'help': "By default, analysis results are stored in the profile database. The use of\
                      this flag will let you skip that"}
                ),
    'min-e-value': (
            ['-e', '--min-e-value'],
            {'metavar': 'E-VALUE',
             'default': 1e-15,
             'type': float,
             'help': "Minimum significance score of an HMM find to be considered as a valid hit.\
                      Default is %(default)g."}
                ),
    'min-percent-identity': (
            ['--min-percent-identity'],
            {'metavar': 'PERCENT_IDENTITY',
             'default': 80.0,
             'type': float,
             'help': "Minimum percent identity. The default is %(default)g."}
                ),
    'min-full-percent-identity': (
            ['--min-full-percent-identity'],
            {'metavar': 'FULL_PERCENT_IDENTITY',
             'default': 20.0,
             'type': float,
             'help': "In some cases you may get high raw ANI estimates (percent identity scores)\
                      between two genomes that have little to do with each other simply because only\
                      a small fraction of their content may be aligned. This can be partly\
                      alleviated by considering the *full* percent identity, which includes in its\
                      calculation regions that did not align. For example, if the alignment is a\
                      whopping 97 percent identity but only 8 percent of the genome aligned, the *full*\
                      percent identity is 0.970 * 0.080 = 0.078 OR 7.8 percent. *full* percent\
                      identity is always included in the report, but you can also use it as a filter\
                      for other metrics, such as percent identity. This filter will set all ANI\
                      measures between two genomes to 0 if the *full* percent identity is less than\
                      you deem trustable. When you set a value, anvi'o will go through the ANI\
                      results, and set all ANI measures between two genomes to 0 if the *full*\
                      percent identity *between either of them* is less than the parameter described\
                      here. The default is %(default)g."}
                ),
    'use-full-percent-identity': (
            ['--use-full-percent-identity'],
            {'action': 'store_true',
             'help': "Usually, percent identity is calculated only over aligned regions, and this\
                      is what is used as a distance metric by default. But with this flag,\
                      you can instead use the *full* percent identity as the distance metric. It is the\
                      same as percent identity, except that regions that did not align are included\
                      in the calculation. This means *full* percent identity will always be less than or\
                      equal to percent identity. How is it calculated? Well if P is the percentage identity\
                      calculated in aligned regions, L is the length of the genome, and A is the fraction\
                      of the genome that aligned to a compared genome, the full percent identity is\
                      P * (A/L). In other words, it is the percent identity multiplied by the alignment\
                      coverage. For example, if the alignment is a whopping 97 percent identity but\
                      only 8 percent of the genome aligned, the *full* percent identity is 0.970 * 0.080\
                      = 0.078, which is just 7.8 percent."}
                ),
    'min-alignment-fraction': (
            ['--min-alignment-fraction'],
            {'default': 0.0,
             'metavar': 'NUM',
             'type': float,
             'help': "In some cases you may get high raw ANI estimates\
                      (percent identity scores) between two genomes that have little to do with each other\
                      simply because only a small fraction of their content may be aligned. This filter will\
                      set all ANI scores between two genomes to 0 if the alignment fraction is less than you\
                      deem trustable. When you set a value, anvi'o will go through the ANI results, and set\
                      percent identity scores between two genomes to 0 if the alignment fraction *between either\
                      of them* is less than the parameter described here. The default is %(default)g."}
                ),
    'significant-alignment-length': (
            ['--significant-alignment-length'],
            {'default': None,
             'metavar': 'INT',
             'type': int,
             'help': "So --min-alignment-fraction\
                      discards any hit that is coming from alignments that represent shorter fractions of genomes,\
                      but what if you still don't want to miss an alignment that is longer than an X number of\
                      nucleotides regardless of what fraction of the genome it represents? Well, this parameter is\
                      to recover things that may be lost due to --min-alignment-fraction parameter. Let's say,\
                      if you set --min-alignment-fraction to '0.05', and this parameter to '5000', anvi'o will keep\
                      hits from alignments that are longer than 5000 nts, EVEN IF THEY REPRESENT less than 5 percent of\
                      a given genome pair. Basically if --min-alignment-fraction is your shield to protect yourself\
                      from incoming garbage, --significant-alignment-length is your chopstick to pick out those that\
                      may be interesting, and you are a true warrior here."}
                ),
    'bins-info': (
            ['--bins-info'],
            {'metavar': 'BINS_INFO',
             'help': "Additional information for bins. The file must contain three TAB-delimited columns,\
                      where the first one must be a unique bin name, the second should be a 'source', and the\
                      last one should be a 7 character HTML color code (i.e., '#424242'). Source column must\
                      contain information about the origin of the bin. If these bins are automatically\
                      identified by a program like CONCOCT, this column could contain the program name and\
                      version. The source information will be associated with the bin in various interfaces\
                      so in a sense it is not *that* critical what it says there, but on the other hand it is,\
                      becuse we should also think about people who may end up having to work with what we put\
                      together later."}
                ),
    'bins': (
            ['--bins'],
            {'metavar': 'BINS_DATA',
             'help': "Tab-delimited file, first column contains tree leaves (gene clusters, splits, contigs etc.) \
                      and second column contains which Bin they belong."
            }
      ),
    'contigs-mode': (
            ['--contigs-mode'],
            {'default': False,
             'action': 'store_true',
             'help': "Use this flag if your binning was done on contigs instead of splits. Please refer\
                      to the documentation for help."}
                ),
    'sample-name': (
            ['-S', '--sample-name'],
            {'metavar': 'NAME',
             'help': "It is important to set a sample name (using only ASCII letters and digits\
                      and without spaces) that is unique (considering all others). If you do not\
                      provide one, anvi'o will try to make up one for you based on other information,\
                      although, you should never let the software to decide these things)."}
                ),
    'project-name': (
            ['-n', '--project-name'],
            {'metavar': 'PROJECT_NAME',
             'help': "Name of the project. Please choose a short but descriptive name (so anvi'o can use\
                      it whenever she needs to name an output file, or add a new table in a database, or name\
                      her first born)."}
                ),
    'skip-hierarchical-clustering': (
            ['--skip-hierarchical-clustering'],
            {'default': False,
             'action': 'store_true',
             'help': "If you are not planning to use the interactive interface (or if you have other\
                      means to add a tree of contigs in the database) you may skip the step where\
                      hierarchical clustering of your items are preformed based on default clustering\
                      recipes matching to your database type."}
                ),
    'skip-variability-tables': (
            ['--skip-variability-tables'],
            {'default': False,
             'action': 'store_true',
             'help': "Processing variability tables in profile database might take a very long time. With\
                      this flag you will be asking anvi'o to skip them."}
                ),
    'enforce-hierarchical-clustering': (
            ['--enforce-hierarchical-clustering'],
            {'default': False,
             'action': 'store_true',
             'help': "If you have more than 25,000 splits in your merged profile, anvi-merge will automatically\
                      skip the hierarchical clustering of splits (by setting --skip-hierarchical-clustering flag\
                      on). This is due to the fact that computational time required for hierarchical clustering\
                      increases exponentially with the number of items being clustered. Based on our experience\
                      we decided that 25,000 splits is about the maximum we should try. However, this is not a\
                      theoretical limit, and you can overwrite this heuristic by using this flag, which would\
                      tell anvi'o to attempt to cluster splits regardless."}
                ),
    'compress-auxiliary-data': (
            ['--compress-auxiliary-data'],
            {'default': False,
             'action': 'store_true',
             'help': "When declared, the auxiliary data file in the resulting output will be compressed. This\
                      saves space, but it takes long. Also, if you are planning to compress the entire\
                      later using GZIP, it is even useless to do. But you are the boss!"}
                ),
    'cluster-contigs': (
            ['--cluster-contigs'],
            {'default': False,
             'action': 'store_true',
             'help': "Single profiles are rarely used for genome binning or visualization, and since\
                      clustering step increases the profiling runtime for no good reason, the default\
                      behavior is to not cluster contigs for individual runs. However, if you are\
                      planning to do binning on one sample, you must use this flag to tell anvio to\
                      run cluster configurations for single runs on your sample."}
                ),
    'num-clusters-requested': (
            ['--num-clusters-requested'],
            {'metavar': 'INT',
             'default': 400,
             'type': int,
             'help': "How many clusters do you request? Default is %(default)d."}
             ),
    'overwrite-output-destinations': (
            ['-W', '--overwrite-output-destinations'],
            {'default': False,
             'action': 'store_true',
             'help': "Overwrite if the output files and/or directories exist."}
                ),
    'delete-if-exists': (
            ['--delete-if-exists'],
            {'default': False,
             'action': 'store_true',
             'help': "Be bold (at your own risk), and delete if exists."}
                ),
    'report-variability-full': (
            ['--report-variability-full'],
            {'default': False,
             'action': 'store_true',
             'help': "One of the things anvi-profile does is to store information about variable\
                      nucleotide positions. Usually it does not report every variable position, since\
                      not every variable position is genuine variation. Say, if you have 1,000 coverage,\
                      and all nucleotides at that position are Ts and only one of them is a C, the\
                      confidence of that C being a real variation is quite low. anvio has a simple\
                      algorithm in place to reduce the impact of noise. However, using this flag\
                      you can disable it and ask profiler to report every single variation (which\
                      may result in very large output files and millions of reports, but you are the\
                      boss). Do not forget to take a look at '--min-coverage-for-variability' parameter"}
                ),
    'report-extended-deflines': (
            ['--report-extended-deflines'],
            {'default': False,
             'action': 'store_true',
             'help': "When declared, the deflines in the resulting FASTA file will contain more information."}
                ),
    'manual-mode': (
            ['--manual-mode'],
            {'default': False,
             'action': 'store_true',
             'help': "Using this flag, you can run the interactive interface in an ad hoc manner using\
                      input files you curated instead of standard output files generated by an anvi'o\
                      run. In the manual mode you will be asked to provide a profile database. In this\
                      mode a profile database is only used to store 'state' of the interactive interface\
                      so you can reload your visual settings when you re-analyze the same files again. If\
                      the profile database you provide does not exist, anvi'o will create an empty one for\
                      you."}
                ),
    'hmm-profile-dir': (
            ['-H', '--hmm-profile-dir'],
            {'metavar': 'HMM PROFILE PATH',
             'help': "You can use this parameter you can specify a directory path that contain an HMM profile.\
                      This way you can run HMM profiles that are not included in anvi'o. See the online\
                      to find out about the specifics of this directory structure ."}
                ),
    'installed-hmm-profile': (
            ['-I', '--installed-hmm-profile'],
            {'metavar': 'HMM PROFILE NAME'}
                ),
    'min-contig-length': (
            ['-M', '--min-contig-length'],
            {'metavar': 'INT',
             'default': 1000,
             'type': int,
             'help': "Minimum length of contigs in a BAM file to analyze. The minimum length should be long enough\
                      for tetra-nucleotide frequency analysis to be meaningful. There is no way to define a golden\
                      number of minimum length that would be applicable to genomes found in all environments, but we\
                      chose the default to be %(default)d, and have been happy with it. You are welcome to experiment,\
                      but we advise to never go below 1,000. You also should remember that the lower you go, the more\
                      time it will take to analyze all contigs. You can use --list-contigs parameter to have an idea how\
                      many contigs would be discarded for a given M."}
                ),
    'max-contig-length': (
            ['--max-contig-length'],
            {'metavar': 'INT',
             'default': 0,
             'type': int,
             'help': "Just like the minimum contig length parameter, but to set a maximum. Basically this will remove\
                      any contig longer than a certain value. Why would anyone need this? Who knows. But if you ever\
                      do, it is here."}
                ),
    'min-mean-coverage': (
            ['-X', '--min-mean-coverage'],
            {'metavar': 'INT',
             'default': 0,
             'type': int,
             'help': "Minimum mean coverage for contigs to be kept in the analysis. The default value is %(default)d,\
                      which is for your best interest if you are going to profile multiple BAM files which are then\
                      going to be merged for a cross-sectional or time series analysis. Do not change it if you are not\
                      sure this is what you want to do."}
                ),
    'min-coverage-for-variability': (
            ['-V', '--min-coverage-for-variability'],
            {'metavar': 'INT',
             'default': 10,
             'type': int,
             'help': "Minimum coverage of a nucleotide position to be subjected to SNV profiling. By default, anvio will\
                      not attempt to make sense of variation in a given nucleotide position if it is covered less than\
                      %(default)dX. You can change that minimum using this parameter."}
                ),
    'contigs-and-positions': (
            ['--contigs-and-positions'],
            {'metavar': 'CONTIGS_AND_POS',
             'required': True,
             'help': "This is the file where you list the contigs, and nucleotide positions you are interested in. This\
                      is supposed to be a TAB-delimited file with two columns. In each line, the first column should be\
                      the contig name, and the second column should be the comma-separated list of integers for nucleotide\
                      positions."}
                ),
    'state-autoload': (
            ['--state-autoload'],
            {'metavar': 'NAME',
             'help': "Automatically load previous saved state and draw tree. To see a list of available states,\
                      use --show-states flag."}
                ),
    'load-full-state': (
            ['--load-full-state'],
            {'required': False,
             'action': 'store_true',
             'help': "Often the minimum and maximum values defined for the an entire profile database that contains\
                      all contigs do not scale well when you wish to work with a single bin in the refine mode. For\
                      this reason, the default behavior of anvi-refine is to ignore min/max values set in the default\
                      state. This flag is your way of telling anvi'o to not do that, and load the state stored in the\
                      profile database as is."}
                ),
    'state': (
            ['-s', '--state'],
            {'metavar': 'STATE',
             'help': "State file, you can export states from database using anvi-export-state program"}
                ),
    'collection-autoload': (
            ['--collection-autoload'],
            {'metavar': 'NAME',
             'help': "Automatically load a collection and draw tree. To see a list of available collections,\
                      use --list-collections flag."}
                ),
    'full-report': (
            ['--full-report'],
            {'metavar': 'FILE_NAME',
             'default': None,
             'help': "Optional output file with a fuller description of findings."}
                ),
    'include-sequences': (
            ['--include-sequences'],
            {'default': False,
             'action': 'store_true',
             'help': "Include sequences in the report."}
                ),
    'show-states': (
            ['--show-states'],
            {'default': False,
             'action': 'store_true',
             'help': "When declared the program will print all available states and exit."}
                ),
    'skip-init-functions': (
            ['--skip-init-functions'],
            {'default': False,
             'action': 'store_true',
             'help': "When declared, function calls for genes will not be initialized (therefore will be missing from all\
                      relevant interfaces or output files). The use of this flag may reduce the memory fingerprint and\
                      processing time for large datasets."}
                ),
    'init-gene-coverages': (
            ['--init-gene-coverages'],
            {'default': False,
             'action': 'store_true',
             'help': "Initialize gene coverage and detection data. This is a very computationally expensive step, but it is\
                      necessary when you need gene level coverage data. The reason this is very computationally expensive\
                      is because anvi'o computes gene coverages by going back to actual coverage values of each gene to\
                      average them, instead of using contig average coverage values, for extreme accuracy."}
                ),
    'skip-auto-ordering': (
            ['--skip-auto-ordering'],
            {'default': False,
             'action': 'store_true',
             'help': "When declared, the attempt to include automatically generated orders of items based on additional data\
                      is skipped. In case those buggers cause issues with your data, and you still want to see your stuff and\
                      deal with the other issue maybe later."}
                ),
    'quick-summary': (
            ['--quick-summary'],
            {'default': False,
             'action': 'store_true',
             'help': "When declared the summary output will be generated as quickly as possible, with minimum amount\
                      of essential information about bins."}
                ),
    'only-complete-links': (
            ['--only-complete-links'],
            {'default': False,
             'action': 'store_true',
             'help': "When declared, only reads that cover all positions will be reported. It is necessary to use this\
                      flag if you want to perform oligotyping-like analyses on matching reads."}
                ),
    'users-data-dir': (
            ['-U', '--users-data-dir'],
            {'metavar': 'USERS_DATA_DIR',
             'type': str,
             'help': "Input directory where the user database is read and stored by the server. A new database will be\
                      created if no directory is found."}
                ),
    'smtp-config-file': (
            ['-E', '--smtp-config-file'],
            {'metavar': 'SMTP_CONFIG_INI',
             'type': str,
             'help': "The configuration file for SMTP server to send e-mails. The input file should be formatted as an INI\
                      file that starts with the header '[SMTP]', and should describe values of each of these variables in\
                      the following lines: 'from_address' (the e-mail address that should appear in the 'From' section of\
                      e-mails sent by the server), 'server_address' (the address of the SMTP server to connect), 'server_port'\
                      (the port number to connect), 'init_tls' (whether to initialize TLS protocol), 'username' (the username\
                      for the server to login, if necessary), 'password' (the password associated with the username for login,\
                      if the password is not blank)."}
                ),
    'validate-users-automatically': (
            ['--validate-users-automatically'],
            {'default': True,
             'action': 'store_true',
             'help': "If this is true, users will not receive a link via email to confirm their account but instead be validated\
                      automatically if there is no smtp configuration."}
                ),
    'queue-size': (
            ['--queue-size'],
            {'default': 0,
             'metavar': 'INT',
             'required': False,
             'help': "The queue size for worker threads to store data to communicate to the main thread. The default is set by the\
                      class based on the number of threads. If you have *any* hesitation about whether you know what you are doing,\
                      you should not change this value."}
                ),
    'write-buffer-size': (
            ['--write-buffer-size'],
            {'default': 500,
             'metavar': 'INT',
             'required': False,
             'help': "How many items should be kept in memory before they are written do the disk. The default is %(default)d.\
                      The larger the buffer size, the less frequent the program will access to the disk, yet the more memory\
                      will be consumed since the processed items will be cleared off the memory only after they are written\
                      to the disk. The default buffer size will likely work for most cases, but if you have very large\
                      contigs, you may need to decrease this value. Please keep an eye on the memory usage output to make sure\
                      the memory use never exceeds the size of the physical memory."}
                ),
    'export-gff3': (
        ['--export-gff3'],
        {
            'default': False,
            'action': 'store_true',
            'help': "If this is true, the output file will be in GFF3 format."
        }
    ),
    'export-svg': (
            ['--export-svg'],
            {'type': str,
             'metavar': 'FILE_PATH',
             'required': False,
             'help': "The SVG output file path."}
                ),
    'tab-delimited': (
            ['--tab-delimited'],
            {'default': False,
             'required': False,
             'action': 'store_true',
             'help': "Use the TAB-delimited format for the output file."}
                ),
    'splits-mode': (
            ['--splits-mode'],
            {'default': False,
             'action': 'store_true',
             'help': "Specify this flag if you would like to output coverages of individual 'splits', rather than their 'parent'\
                      contig coverages."}
                ),
    'report-as-text': (
            ['--report-as-text'],
            {'default': False,
             'action': 'store_true',
             'help': "If you give this flag, Anvi'o will not open new browser to show Contigs database statistics and write all stats \
                      to TAB separated file and you should also give --output-file with this flag otherwise Anvi'o will complain."}
                ),
    'dump-dir': (
            ['--dump-dir'],
            {'required': False,
             'help': "Modelling and annotating structures requires a lot of moving parts, each which have \
                      their own outputs. The output of this program is a structure database containing the \
                      pertinent results of this computation, however a lot of stuff doesn't make the cut. \
                      By providing a directory for this parameter you will get, in addition to the structure \
                      database, a directory containing the raw output for everything."}
                ),
    'workflow': (
            ['-w', '--workflow'],
            {'required': False,
             'help': "You must specify a workflow name. To see a list of available workflows\
                      run --list-workflows."}
                ),
    'list-workflows': (
            ['--list-workflows'],
            {'required': False,
             'action': 'store_true',
             'help': "Print a list of available snakemake workflows"}
                ),
    'save-workflow-graph': (
            ['--save-workflow-graph'],
            {'required': False,
             'action': 'store_true',
             'help': "Save a graph representation of the workflow. If you are using this flag and if your\
                      system is unable to generate such graph outputs, you will hear anvi'o complaining\
                      (still, totally worth trying)."}
                ),
    'get-default-config': (
            ['--get-default-config'],
            {'metavar': 'OUTPUT_FILENAME',
             'type': str,
             'help': "Store a json formatted config file with all the default settings of the\
                      workflow. This is a good draft you could use in order to write your own\
                      config file. This config file contains all parameters that could be configured\
                      for this workflow. NOTICE: the config file is provided with default values\
                      only for parameters that are set by us in the workflow. The values for the rest\
                      of the parameters are determined by the relevant program."}
                ),
    'list-dependencies': (
            ['--list-dependencies'],
            {'required': False,
             'action': 'store_true',
             'help': "Print a list of the dependencies of this workflow. You must provide a workflow name\
                      and a config file. snakemake will figure out which rules need to be run according\
                      to your config file, and according to the files available on your disk. According\
                      to the rules that need to be run, we will let you know which programs are going to\
                      be used, so that you can make sure you have all of them installed and loaded."}
                ),
    'config-file': (
            ['-c', '--config-file'],
            {'required': False,
             'help': "TBD"}
                ),
    'additional-params': (
            ['-A', '--additional-params'],
            {'required': False,
             'nargs':'...', 'type':str,
             'help': "Additional snakemake parameters to add when running snakemake. NOTICE: --additional-params \
                      HAS TO BE THE LAST ARGUMENT THAT IS PASSED TO anvi-run-workflow, ANYTHING THAT \
                      FOLLOWS WILL BE CONSIDERED AS PART OF THE ADDITIONAL PARAMETERS THAT ARE PASSED TO SNAKEMAKE. \
                      Any parameter that is accepted by snakemake should be fair game here, but it is your \
                      responsibility to make sure that whatever you added makes sense. To see what parameters are \
                      available please refer to the snakemake documentation. For example, you could use this to set \
                      up cluster submission using --additional-params --cluster \"YOUR-CLUSTER-SUBMISSION-CMD\""}
                ),
    'self-key': (
            ['--self-key'],
            {'default': None,
             'type': str,
             'help': "The key you wish to set or change."}
                ),
    'self-value': (
            ['--self-value'],
            {'default': None,
             'type': str,
             'help': "The value you wish to set for the self key."}
                ),
    'no-variability': (
            ['--no-variability'],
            {'required': False,
             'action': 'store_true',
             'help': "If provided, no measures of sequence heterogeneity (from short read data) will be overlayed\
                      on structures."}
                ),
    'compute-gene-coverage-stats': (
            ['--compute-gene-coverage-stats'],
            {'required': False,
             'action': 'store_true',
             'help': "If provided, gene coverage statistics will be appended for each entry in variability report.\
                      This is very useful information, but will not be included by default because it is an expensive\
                      operation, and may take some additional time."}
                ),
    'repository': (
            ['--repository'],
            {'default': 'merenlab/anvio',
             'type': str,
             'help': "Source repository to download releases, currently only Github is supported. Enter in 'merenlab/anvio' format."}
                ),
    'inseq-stats': (
            ['--inseq-stats'],
            {'required': False,
             'action': 'store_true',
             'default': False,
             'help': "Provide if working with INSeq/Tn-Seq genomic data. With this, all gene level \
                      coverage stats will be calculated using INSeq/Tn-Seq statistical methods."}
                ),
}

# two functions that works with the dictionary above.
def A(param_id, exclude_param=None):
    if exclude_param:
        return [p for p in D[param_id][0] if p != exclude_param]
    else:
        return D[param_id][0]

def K(param_id, params_dict={}):
    kwargs = copy.deepcopy(D[param_id][1])
    for key in params_dict:
        kwargs[key] = params_dict[key]

    return kwargs

# The rest of this file is composed of code that responds to '-v' or '--version' calls from clients,
# and provides access to the database version numbers for all anvi'o modules.

import anvio.tables as t
from anvio.terminal import Run


run = Run()


def set_version():
    return anvio_version, \
           anvio_codename, \
           t.contigs_db_version, \
           t.pan_db_version, \
           t.profile_db_version, \
           t.genes_db_version, \
           t.auxiliary_data_version, \
           t.genomes_storage_vesion, \
           t.structure_db_version


def get_version_tuples():
    return [("Anvi'o version", "%s (v%s)" % (__codename__, __version__)),
            ("Profile DB version", __profile__version__),
            ("Contigs DB version", __contigs__version__),
            ("Genes DB version", __genes__version__),
            ("Auxiliary data storage version", __auxiliary_data_version__),
            ("Pan DB version", __pan__version__),
            ("Genome data storage version", __genomes_storage_version__),
            ("Structure DB version", __structure__version__)]


def print_version():
    run.info("Anvi'o version", "%s (v%s)" % (__codename__, __version__), mc='green')
    run.info("Profile DB version", __profile__version__)
    run.info("Contigs DB version", __contigs__version__)
    run.info("Pan DB version", __pan__version__)
    run.info("Genome data storage version", __genomes_storage_version__)
    run.info("Auxiliary data storage version", __auxiliary_data_version__)
    run.info("Structure DB version", __structure__version__)


__version__, \
__codename__, \
__contigs__version__, \
__pan__version__, \
__profile__version__, \
__genes__version__, \
__auxiliary_data_version__, \
__genomes_storage_version__ , \
__structure__version__ = set_version()


if '-v' in sys.argv or '--version' in sys.argv:
    print_version()
    sys.exit()
