# -*- coding: utf-8
"""Everything about anvi'o workflows and artifacts"""

# when defining a new anvi'o workflow, it is essential to document most critical
# artifacts produced by a given workflow. it is also important to mention the
# third party programs used (such as Bowtie2 or HMMER). there is no need to
# mention anvi'o programs here as it will be discovered from program tags. to
# associate a program with a workflow, please add the workflow name(s) to the
# meta tag `__anvio_workflows__`.
ANVIO_WORKFLOWS = {
    "contigs": {
        "authors": ['ShaiberAlon', 'ivagljiva', 'meren', 'mschecht'],
        "artifacts_produced": ['contigs-db'],
        "artifacts_accepted": ['fasta-txt'],
        "anvio_workflows_inherited": [],
        "third_party_programs_used": [
            ('Gene calling', ['prodigal']),
            ('HMM search', ['HMMER']),
            ('Gene taxonomy', ['krakenuniq', 'centrifuge']),
            ('Sequence search against various databases', ['DIAMOND'])
            ],
        "one_sentence_summary": "From FASTA files to annotated anvi'o contigs databases",
        "one_paragraph_summary": ("Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor "
            "incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco "
            "laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate "
            "velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, "
            "sunt in culpa qui officia deserunt mollit anim id est laborum")
    },

    "metagenomics": {
        "authors": ['ShaiberAlon'],
        "artifacts_produced": ['contigs-db', 'profile-db'],
        "artifacts_accepted": ['samples-txt', 'fasta-txt'],
        "anvio_workflows_inherited": ['contigs'],
        "third_party_programs_used": [
            ('Quality control of short reads', ['illumina-utils']),
            ('Assembly', ['IDBA-UD', 'metaSPAdes', 'MEGAHIT']),
            ('BAM file manipulations', ['samtools']),
            ('Gene calling', ['prodigal']),
            ('HMM search', ['HMMER']),
            ('Gene taxonomy', ['krakenuniq', 'centrifuge']),
            ('Read recruitment', ['Bowtie2'])
            ],
        "one_sentence_summary": "From FASTA and/or FASTQ files to anvi'o contigs and profile databases",
        "one_paragraph_summary": ("Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor "
            "incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco "
            "laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate "
            "velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, "
            "sunt in culpa qui officia deserunt mollit anim id est laborum")
    },
    "ecophylo": {
        "authors": ['mschecht'],
        "artifacts_accepted": ['samples-txt', 'hmm-list', 'external-genomes', 'metagenomes'],
        "artifacts_produced": ['contigs-db', 'profile-db'],
        "anvio_workflows_inherited": [],
        "third_party_programs_used": [
            ('Read recruitment', ['Bowtie2']),
            ('Cluster open reading frames', ['MMseqs2']),
            ('Align protein sequences', ['muscle']),
            ('Trim multiple sequence alignment', ['trimal']),
            ('Calculate phylogenetic tree', ['IQ-TREE']),
            ('Calculate phylogenetic tree', ['FastTree']),
            ('Search for homologous sequences', ['HMMER'])
            ],
        "one_sentence_summary": "Co-characterize the biogeography and phylogeny of any protein",
        "one_paragraph_summary": ("The ecophylo workflow explores the **eco**logical and **phylo**genetic relationships between individual genes and environments. "
            "Briefly, the workflow extracts a target gene from any set of FASTA files (e.g., isolate genomes, [MAGs](https://anvio.org/vocabulary/#metagenome-assembled-genome-mag), "
            "[SAGs](https://anvio.org/vocabulary/#single-amplified-genome-sag), or simply assembled metagenomes) "
            "using a user-defined [HMM](https://anvio.org/vocabulary/#hidden-markov-models-hmms), and offers an integrated access "
            "to the phylogenetics of matching genes, and their distribution across environments.")
    },
    "sra-download": {
        "authors": ['mschecht'],
        "artifacts_accepted": [],
        "artifacts_produced": ['paired-end-fastq'],
        "anvio_workflows_inherited": [],
        "third_party_programs_used": [
            ('Downloads SRA accessions', ['prefetch']),
            ('Extracts FASTQ files from SRA accessions', ['fasterq-dump']),
            ('Compresses FASTQ files in parallel', ['pigz']),
            ],
        "one_sentence_summary": "Download, extract, and gzip paired-end FASTQ files automatically from the NCBI short-read archive (SRA)",
        "one_paragraph_summary": ("The sra-download workflow automatizes the process of downloading paired-end FASTQ files "
            "for a given list of SRA-accessions using [NCBI sra-tools wiki](https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump) "
            "then gzips them using [pigz](https://zlib.net/pigz/).")
    },
}

# the purpose of this variable is to have a list of third-party programs used from
# within anvi'o workflows.
THIRD_PARTY_PROGRAMS = {
    'IDBA-UD': {'link': 'https://github.com/loneknightpy/idba'},
    'illumina-utils': {'link': 'https://github.com/merenlab/illumina-utils'},
    'metaSPAdes': {'link': "https://cab.spbu.ru/software/meta-spades/"},
    'MEGAHIT': {'link': 'https://github.com/voutcn/megahit'},
    'samtools': {'link': 'http://www.htslib.org/'},
    'prodigal': {'link': 'https://github.com/hyattpd/Prodigal'},
    'HMMER': {'link': 'http://hmmer.org/'},
    'Bowtie2': {'link': 'https://github.com/BenLangmead/bowtie2'},
    'krakenuniq': {'link': 'https://github.com/fbreitwieser/krakenuniq'},
    'centrifuge': {'link': 'https://github.com/DaehwanKimLab/centrifuge'},
    'DIAMOND': {'link': 'https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/diamond/'},
    'MMseqs2': {'link': 'https://github.com/soedinglab/MMseqs2'},
    'muscle': {'link': 'http://www.drive5.com/muscle/'},
    'FastTree': {'link': 'http://www.microbesonline.org/fasttree/'},
    'IQ-TREE': {'link': 'https://github.com/Cibiv/IQ-TREE'},
    'trimal': {'link': 'https://github.com/inab/trimal'},
    'pigz': {'link': 'https://zlib.net/pigz/'},
    'prefetch': {'link': 'https://github.com/ncbi/sra-tools'},
    'fasterq-dump': {'link': 'https://github.com/ncbi/sra-tools'}
    }

# the purpose of dictionaries in this module is to describes all anvi'o items and concepts
# that are referred from 'requires' and 'provides' statements in anvi'o programs
ANVIO_ARTIFACTS ={
    "pan-db": {
        "name": "PAN",
        "type": "DB",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "contigs-db": {
        "name": "CONTIGS",
        "type": "DB",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "trnaseq-db": {
        "name": "TRNASEQ",
        "type": "DB",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "trnaseq-contigs-db": {
        "name": "TRNASEQ CONTIGS",
        "type": "DB",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "trnaseq-profile-db": {
        "name": "TRNASEQ PROFILE",
        "type": "DB",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "modules-db": {
        "name": "MODULES",
        "type": "DB",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "fasta": {
        "name": "REGULAR FASTA",
        "type": "FASTA",
        "provided_by_anvio": False,
        "provided_by_user": True
    },
    "contigs-fasta": {
        "name": "CONTIGS",
        "type": "FASTA",
        "provided_by_anvio": True,
        "provided_by_user": True
    },
    "dna-sequence": {
        "name": "DNA SEQUENCE",
        "type": "SEQUENCE",
        "provided_by_anvio": True,
        "provided_by_user": True
    },
    "trnaseq-fasta": {
        "name": "TRNASEQ",
        "type": "FASTA",
        "provided_by_anvio": False,
        "provided_by_user": True
    },
    "configuration-ini": {
        "name": "CONFIGURATION FILE",
        "type": "TXT",
        "provided_by_anvio": False,
        "provided_by_user": True
    },
    "external-gene-calls": {
        "name": "EXTERNAL GENE CALLS",
        "type": "TXT",
        "provided_by_anvio": False,
        "provided_by_user": True
    },
    "external-structures": {
        "name": "EXTERNAL STRUCTURES",
        "type": "TXT",
        "provided_by_anvio": False,
        "provided_by_user": True
    },
    "concatenated-gene-alignment-fasta": {
        "name": "CONCATENATED GENE ALIGNMENT",
        "type": "FASTA",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "short-reads-fasta": {
        "name": "SHORT READS",
        "type": "FASTA",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "paired-end-fastq": {
        "name": "SHORT READS",
        "type": "FASTQ",
        "provided_by_anvio": True,
        "provided_by_user": True
    },
    "genes-fasta": {
        "name": "GENES",
        "type": "FASTA",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "bam-file": {
        "name": "BAM FILE",
        "type": "BAM",
        "provided_by_anvio": False,
        "provided_by_user": True
    },
    "bam-stats-txt": {
        "name": "BAM STATS TXT",
        "type": "TXT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "bams-and-profiles-txt": {
        "name": "BAMS AND PROFILES TXT",
        "type": "TXT",
        "provided_by_anvio": False,
        "provided_by_user": True
    },
    "markdown-txt": {
        "name": "MARKDOWN TXT",
        "type": "TXT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "protein-structure-txt": {
        "name": "PDB FILE",
        "type": "TXT",
        "provided_by_anvio": False,
        "provided_by_user":True
    },
    "samples-txt": {
        "name": "SAMPLES TXT",
        "type": "TXT",
        "provided_by_anvio": False,
        "provided_by_user":True
    },
    "primers-txt": {
        "name": "PRIMERS TXT",
        "type": "TXT",
        "provided_by_anvio": False,
        "provided_by_user":True
    },
    "fasta-txt": {
        "name": "FASTA TXT",
        "type": "TXT",
        "provided_by_anvio": False,
        "provided_by_user":True
    },
    "raw-bam-file": {
        "name": "RAW BAM FILE",
        "type": "BAM",
        "provided_by_anvio": False,
        "provided_by_user":True
    },
    "locus-fasta": {
        "name": "LOCUS",
        "type": "FASTA",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "structure-db": {
        "name": "PROTEIN STRUCTURES",
        "type": "DB",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "pdb-db": {
        "name": "PDB DB",
        "type": "DB",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "kegg-data": {
        "name": "KEGG MODULES DB",
        "type": "DB",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "user-modules-data": {
        "name": "USER MODULES DB",
        "type": "DB",
        "provided_by_anvio": True,
        "provided_by_user": True
    },
    "single-profile-db": {
        "name": "SINGLE PROFILE",
        "type": "DB",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "profile-db": {
        "name": "PROFILE",
        "type": "DB",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "genes-db": {
        "name": "GENES",
        "type": "DB",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "genomes-storage-db": {
        "name": "GENOMES STORAGE",
        "type": "DB",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "contigs-stats": {
        "name": "CONTIGS STATS",
        "type": "STATS",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "svg": {
        "name": "SVG",
        "type": "SVG",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "bin": {
        "name": "BIN",
        "type": "BIN",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "collection": {
        "name": "COLLECTION",
        "type": "COLLECTION",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "collection-txt": {
        "name": "COLLECTION",
        "type": "TXT",
        "provided_by_anvio": True,
        "provided_by_user":True
    },
    "hmm-source": {
        "name": "HMM SOURCE",
        "type": "HMM",
        "provided_by_anvio": False,
        "provided_by_user":True
    },
    "hmm-hits": {
        "name": "HMM PROFILE",
        "type": "CONCEPT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "completion": {
        "name": "COMPLETION",
        "type": "CONCEPT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "cogs-data": {
        "name": "COGs DATA",
        "type": "DATA",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "pfams-data": {
        "name": "PFAMs DATA",
        "type": "DATA",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "cazyme-data": {
        "name": "CAZymes HMM DATA",
        "type": "DATA",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "misc-data-items-txt": {
        "name": "ITEMS DATA",
        "type": "TXT",
        "provided_by_anvio": True,
        "provided_by_user":True
    },
    "misc-data-items": {
        "name": "ITEMS DATA",
        "type": "CONCEPT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "misc-data-layers-txt": {
        "name": "LAYERS DATA",
        "type": "TXT",
        "provided_by_anvio": True,
        "provided_by_user": True
    },
    "misc-data-layers": {
        "name": "LAYERS DATA",
        "type": "CONCEPT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "misc-data-nucleotides-txt": {
        "name": "NUCLEOTIDES DATA",
        "type": "TXT",
        "provided_by_anvio": True,
        "provided_by_user": True
    },
    "misc-data-nucleotides": {
        "name": "NUCLEOTIDES DATA",
        "type": "CONCEPT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "misc-data-amino-acids-txt": {
        "name": "AMINO_ACIDS DATA",
        "type": "TXT",
        "provided_by_anvio": True,
        "provided_by_user": True
    },
    "misc-data-amino-acids": {
        "name": "AMINO_ACIDS DATA",
        "type": "CONCEPT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "genome-similarity": {
        "name": "GENOME SIMILARITY",
        "type": "CONCEPT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "misc-data-layer-orders-txt": {
        "name": "LAYER ORDERS DATA",
        "type": "TXT",
        "provided_by_anvio": True,
        "provided_by_user": True
    },
    "misc-data-layer-orders": {
        "name": "LAYER ORDERS DATA",
        "type": "CONCEPT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "misc-data-items-order-txt": {
        "name": "LAYER ORDERS DATA",
        "type": "TXT",
        "provided_by_anvio": True,
        "provided_by_user": True
    },
    "misc-data-items-order": {
        "name": "ITEM ORDERS DATA",
        "type": "CONCEPT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "dendrogram": {
        "name": "DENDROGRAM",
        "type": "NEWICK",
        "provided_by_anvio": True,
        "provided_by_user": True
    },
    "metapangenome": {
        "name": "METAPANGENOME",
        "type": "CONCEPT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "oligotypes": {
        "name": "OLIGOTYPES",
        "type": "CONCEPT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "linkmers-txt": {
        "name": "LINKMERS",
        "type": "TXT",
        "provided_by_anvio": False,
        "provided_by_user": True
    },
    "palindromes-txt": {
        "name": "PALINDROMES",
        "type": "TXT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "inversions-txt": {
        "name": "INVERSIONS TXT",
        "type": "TXT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "phylogeny": {
        "name": "PHYLOGENY",
        "type": "NEWICK",
        "provided_by_anvio": True,
        "provided_by_user": True
    },
    "gene-calls-txt": {
        "name": "GENE CALLS",
        "type": "TXT",
        "provided_by_anvio": True,
        "provided_by_user": True
    },
    "binding-frequencies-txt": {
        "name": "AMINO ACID BINDING FREQUENCIES TEXT",
        "type": "TXT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "interacdome-data": {
        "name": "INTERACDOME DATA",
        "type": "DATA",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "functions": {
        "name": "GENE FUNCTIONS",
        "type": "CONCEPT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "functions-txt": {
        "name": "GENE FUNCTIONS",
        "type": "TXT",
        "provided_by_anvio": True,
        "provided_by_user": True
    },
    "functional-enrichment-txt": {
        "name": "ENRICHMENT SCORES",
        "type": "TXT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "kegg-functions": {
        "name": "KOFAM FUNCTIONS",
        "type": "CONCEPT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "functions-across-genomes-txt": {
        "name": "FUNCTIONS ACROSS GENOMES",
        "type": "TXT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "hmm-hits-across-genomes-txt": {
        "name": "HMM HITS ACROSS GENOMES",
        "type": "TXT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "interactive": {
        "name": "INTERACTIVE DISPLAY",
        "type": "DISPLAY",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "genome-view": {
        "name": "GENOME VIEW",
        "type": "DISPLAY",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "view-data": {
        "name": "VIEW DATA",
        "type": "TXT",
        "provided_by_anvio": True,
        "provided_by_user": True
    },
    "layer-taxonomy": {
        "name": "LAYER TAXONOMY",
        "type": "CONCEPT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "layer-taxonomy-txt": {
        "name": "LAYER TAXONOMY",
        "type": "TXT",
        "provided_by_anvio": True,
        "provided_by_user": True
    },
    "gene-taxonomy": {
        "name": "GENE TAXONOMY",
        "type": "CONCEPT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "gene-taxonomy-txt": {
        "name": "GENE TAXONOMY",
        "type": "TXT",
        "provided_by_anvio": True,
        "provided_by_user": True
    },
    "genome-taxonomy": {
        "name": "GENOME TAXONOMY",
        "type": "CONCEPT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "genome-taxonomy-txt": {
        "name": "GENOME TAXONOMY",
        "type": "TXT",
        "provided_by_anvio": False,
        "provided_by_user": False
    },
    "scgs-taxonomy-db": {
        "name": "SCG TAXONOMY DB",
        "type": "CONCEPT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "scgs-taxonomy": {
        "name": "SCG TAXONOMY",
        "type": "CONCEPT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "trna-taxonomy-db": {
        "name": "TRNA TAXONOMY DB",
        "type": "CONCEPT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "trna-taxonomy": {
        "name": "TRNA TAXONOMY",
        "type": "CONCEPT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "external-genomes": {
        "name": "EXTERNAL GENOMES",
        "type": "TXT",
        "provided_by_anvio": False,
        "provided_by_user": True
    },
    "internal-genomes": {
        "name": "INTERNAL GENOMES",
        "type": "TXT",
        "provided_by_anvio": False,
        "provided_by_user": True
    },
    "metagenomes": {
        "name": "METAGENOMES",
        "type": "TXT",
        "provided_by_anvio": False,
        "provided_by_user": True
    },
    "hmm-list": {
        "name": "HMM-LIST",
        "type": "TXT",
        "provided_by_anvio": False,
        "provided_by_user": True
    },
    "coverages-txt": {
        "name": "COVERAGES",
        "type": "TXT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "detection-txt": {
        "name": "DETECTIONS",
        "type": "TXT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "variability-profile": {
        "name": "VARIABILITY PROFILE",
        "type": "CONCEPT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "variability-profile-txt": {
        "name": "VARIABILITY PROFILE",
        "type": "TXT",
        "provided_by_anvio": True,
        "provided_by_user": True
    },
    "variability-profile-xml": {
        "name": "VARIABILITY NETWORK",
        "type": "XML",
        "provided_by_anvio": True,
        "provided_by_user": True
    },
    "codon-frequencies-txt": {
        "name": "CODON FREQUENCIES",
        "type": "TXT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "aa-frequencies-txt": {
        "name": "AA FREQUENCIES",
        "type": "TXT",
        "provided_by_anvio": True ,
        "provided_by_user": False
    },
    "fixation-index-matrix": {
        "name": "FIXATION INDEX MATRIX",
        "type": "TXT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "trnaseq-seed-txt": {
        "name": "TRNASEQ SEED SUMMARY",
        "type": "TXT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "seeds-specific-txt": {
        "name": "TRNASEQ SEED SUMMARY",
        "type": "TXT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "seeds-non-specific-txt": {
        "name": "TRNASEQ SEED SUMMARY",
        "type": "TXT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "modifications-txt": {
        "name": "TRNASEQ MODIFICATION SUMMARY",
        "type": "TXT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "trnaseq-plot": {
        "name": "TRNASEQ PLOT",
        "type": "DISPLAY",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "summary": {
        "name": "STATIC SUMMARY",
        "type": "SUMMARY",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "quick-summary": {
        "name": "QUICK SUMMARY",
        "type": "TXT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "split-bins": {
        "name": "SPLIT BINS",
        "type": "CONCEPT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "state": {
        "name": "INTERACTIVE STATE",
        "type": "CONCEPT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "ngrams": {
        "name": "NGRAM",
        "type": "CONCEPT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "state-json": {
        "name": "INTERACTIVE STATE",
        "type": "JSON",
        "provided_by_anvio": True,
        "provided_by_user": True
    },
    "kegg-metabolism": {
        "name": "KEGG METABOLISM ESTIMATES",
        "type": "TXT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "user-metabolism": {
        "name": "USER METABOLISM ESTIMATES",
        "type": "TXT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "augustus-gene-calls": {
        "name": "AUGUSTUS GENE CALLS",
        "type": "TXT",
        "provided_by_anvio": False,
        "provided_by_user": True
    },
    "genes-stats": {
        "name": "GENE STATS",
        "type": "STATS",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "vcf": {
        "name": "VCF",
        "type": "TXT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "blast-table": {
        "name": "BLAST TABLE",
        "type": "TXT",
        "provided_by_anvio": False,
        "provided_by_user": True
    },
    "splits-txt": {
        "name": "SPLITS",
        "type": "TXT",
        "provided_by_anvio": False,
        "provided_by_user": True
    },
    "genbank-file": {
        "name": "GENBANK FILE",
        "type": "TXT",
        "provided_by_anvio": False,
        "provided_by_user": True
    },
    "groups-txt": {
        "name": "GROUPS",
        "type": "TXT",
        "provided_by_anvio": False,
        "provided_by_user": True
    },
    "splits-taxonomy-txt": {
        "name": "SPLITS TAXONOMY",
        "type": "TXT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "pn-ps-data": {
        "name": "PN/PS OUTPUT",
        "type": "CONCEPT",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "clustering-configuration": {
        "name": "CLUSTERING CONFIG",
        "type": "TXT",
        "provided_by_anvio": False,
        "provided_by_user": True
    },
    "workflow-config": {
        "name": "WORKFLOW CONFIG",
        "type": "JSON",
        "provided_by_anvio": False,
        "provided_by_user": True
    },
    "contigs-workflow": {
        "name": "CONTIGS WORKFLOW",
        "type": "WORKFLOW",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "metagenomics-workflow": {
        "name": "METAGENOMICS WORKFLOW",
        "type": "WORKFLOW",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "pangenomics-workflow": {
        "name": "PANGENOMICS WORKFLOW",
        "type": "WORKFLOW",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "phylogenomics-workflow": {
        "name": "PHYLOGENOMICS WORKFLOW",
        "type": "WORKFLOW",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "trnaseq-workflow": {
        "name": "TRNASEQ WORKFLOW",
        "type": "WORKFLOW",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "ecophylo-workflow": {
        "name": "ECOPHYLO WORKFLOW",
        "type": "WORKFLOW",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "sra-download-workflow": {
        "name": "SRA-DOWNLOAD WORKFLOW",
        "type": "WORKFLOW",
        "provided_by_anvio": True,
        "provided_by_user": False
    },
    "contig-inspection" : {
        "name" : "CONTIG INSPECTION",
        "type" : "DISPLAY",
        "provided_by_anvio" : True,
        "provided_by_user" : False
    },
    "gene-cluster-inspection" : {
        "name" : "GENE CLUSTER INSPECTION",
        "type" : "DISPLAY",
        "provided_by_anvio" : True,
        "provided_by_user" : False
    },
    "enzymes-txt": {
        "name": "ENZYMES TXT",
        "type": "TXT",
        "provided_by_anvio": False,
        "provided_by_user":True
    },
    "enzymes-list-for-module": {
        "name": "ENZYMES LIST",
        "type": "TXT",
        "provided_by_anvio": False,
        "provided_by_user":True
    },
    "metabolic-independence-score": {
        "name": "METABOLIC INDEPENDENCE SCORE",
        "type": "CONCEPT",
        "provided_by_anvio": True,
        "provided_by_user": False
    }
}
