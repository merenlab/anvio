# -*- coding: utf-8
"""Everything about anvi'o artifacts"""

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
        "provided_by_user":True
    },
    "bam-stats-txt": {
        "name": "BAM STATS TXT",
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
        "name": "KEGG DB",
        "type": "DB",
        "provided_by_anvio": True,
        "provided_by_user": False
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
    "hmm-hits-matrix-txt": {
        "name": "HMM HITS MATRIX",
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
    }
}
