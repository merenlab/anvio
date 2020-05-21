# -*- coding: utf-8
"""Everything about anvi'o artifacts"""

# the purpose of dictionaries in this module is to describes all anvi'o items and concepts
# that are referred from 'requires' and 'provides' statements in anvi'o programs
ANVIO_ARTIFACTS ={
    "pan-db": {
        "name": "PAN",
        "type": "DB",
        "internal": True
    },
    "contigs-db": {
        "name": "CONTIGS",
        "type": "DB",
        "internal": True
    },
    "fasta": {
        "name": "REGULAR FASTA",
        "type": "FASTA",
        "internal": False
    },
    "contigs-fasta": {
        "name": "CONTIGS",
        "type": "FASTA",
        "internal": False
    },
    "external-gene-calls": {
        "name": "EXTERNAL GENE CALLS",
        "type": "TXT",
        "internal": False
    },
    "concatenated-gene-alignment-fasta": {
        "name": "CONCATENATED GENE ALIGNMENT",
        "type": "FASTA",
        "internal": False
    },
    "short-reads-fasta": {
        "name": "SHORT READS",
        "type": "FASTA",
        "internal": False
    },
    "genes-fasta": {
        "name": "GENES",
        "type": "FASTA",
        "internal": False
    },
    "bam-file": {
        "name": "BAM FILE",
        "type": "BAM",
        "internal": False
    },
    "protein-structure": {
        "name": "PDB FILE",
        "type": "TXT",
        "internal": False
    },
    "raw-bam-file": {
        "name": "RAW BAM FILE",
        "type": "BAM",
        "internal": False
    },
    "locus-fasta": {
        "name": "LOCUS",
        "type": "FASTA",
        "internal": False
    },
    "structure-db": {
        "name": "STRUCTURE",
        "type": "DB",
        "internal": True
    },
    "pdb-db": {
        "name": "PDB DB",
        "type": "CONCEPT",
        "internal": True
    },
    "kegg-db": {
        "name": "KEGG DB",
        "type": "DB",
        "internal": True
    },
    "single-profile-db": {
        "name": "SINGLE PROFILE",
        "type": "DB",
        "internal": True
    },
    "profile-db": {
        "name": "PROFILE",
        "type": "DB",
        "internal": True
    },
    "genes-db": {
        "name": "GENES",
        "type": "DB",
        "internal": True
    },
    "genomes-storage-db": {
        "name": "GENOMES STORAGE",
        "type": "DB",
        "internal": True
    },
    "contigs-stats": {
        "name": "CONTIGS STATS",
        "type": "STATS",
        "internal": False
    },
    "svg": {
        "name": "SVG",
        "type": "SVG",
        "internal": False
    },
    "bin": {
        "name": "BIN",
        "type": "BIN",
        "internal": True
    },
    "collection": {
        "name": "COLLECTION",
        "type": "COLLECTION",
        "internal": True
    },
    "collection-txt": {
        "name": "COLLECTION",
        "type": "TXT",
        "internal": False
    },
    "hmm-source": {
        "name": "HMM SOURCE",
        "type": "HMM",
        "internal": False
    },
    "hmm-profile": {
        "name": "HMM PROFILE",
        "type": "CONCEPT",
        "internal": True
    },
    "cogs-data": {
        "name": "COGs DATA",
        "type": "DATA",
        "internal": True
    },
    "pfams-data": {
        "name": "PFAMs DATA",
        "type": "DATA",
        "internal": True
    },
    "trna-genes": {
        "name": "TRNA GENES",
        "type": "CONCEPT",
        "internal": True
    },
    "misc-data-items-txt": {
        "name": "ITEMS DATA",
        "type": "TXT",
        "internal": False
    },
    "misc-data-items": {
        "name": "ITEMS DATA",
        "type": "CONCEPT",
        "internal": False
    },
    "misc-data-layers-txt": {
        "name": "LAYERS DATA",
        "type": "TXT",
        "internal": False
    },
    "misc-data-layers": {
        "name": "LAYERS DATA",
        "type": "CONCEPT",
        "internal": False
    },
    "misc-data-nucleotides-txt": {
        "name": "NUCLEOTIDES DATA",
        "type": "TXT",
        "internal": False
    },
    "misc-data-nucleotides": {
        "name": "NUCLEOTIDES DATA",
        "type": "CONCEPT",
        "internal": False
    },
    "misc-data-amino-acids-txt": {
        "name": "AMINO_ACIDS DATA",
        "type": "TXT",
        "internal": False
    },
    "misc-data-amino-acids": {
        "name": "AMINO_ACIDS DATA",
        "type": "CONCEPT",
        "internal": False
    },
    "misc-data-layers-category": {
        "name": "LAYERS DATA CATEGORY",
        "type": "CONCEPT",
        "internal": True
    },
    "genome-similarity": {
        "name": "GENOME SIMILARITY",
        "type": "CONCEPT",
        "internal": False
    },
    "misc-data-layer-orders": {
        "name": "LAYER ORDERS DATA",
        "type": "CONCEPT",
        "internal": False
    },
    "misc-data-item-orders": {
        "name": "ITEM ORDERS DATA",
        "type": "CONCEPT",
        "internal": False
    },
    "misc-data-layer-orders-txt": {
        "name": "LAYER ORDERS DATA",
        "type": "TXT",
        "internal": False
    },
    "misc-data-item-orders-txt": {
        "name": "LAYER ORDERS DATA",
        "type": "TXT",
        "internal": False
    },
    "dendrogram": {
        "name": "DENDROGRAM",
        "type": "NEWICK",
        "internal": False
    },
    "metapangenome": {
        "name": "METAPANGENOME",
        "type": "CONCEPT",
        "internal": True
    },
    "oligotypes": {
        "name": "OLIGOTYPES",
        "type": "CONCEPT",
        "internal": False
    },
    "linkmers-txt": {
        "name": "LINKMERS",
        "type": "TXT",
        "internal": False
    },
    "phylogeny": {
        "name": "PHYLOGENY",
        "type": "NEWICK",
        "internal": False
    },
    "gene-calls-txt": {
        "name": "GENE CALLS",
        "type": "TXT",
        "internal": False
    },
    "functions": {
        "name": "GENE FUNCTIONS",
        "type": "CONCEPT",
        "internal": True
    },
    "functions-txt": {
        "name": "GENE FUNCTIONS",
        "type": "TXT",
        "internal": False
    },
    "functional-enrichment-txt": {
        "name": "ENRICHED FUNCTIONS",
        "type": "TXT",
        "internal": False
    },
    "kegg-functions": {
        "name": "KOFAM FUNCTIONS",
        "type": "CONCEPT",
        "internal": True
    },
    "interactive": {
        "name": "INTERACTIVE DISPLAY",
        "type": "DISPLAY",
        "internal": True
    },
    "view-data": {
        "name": "VIEW DATA",
        "type": "TXT",
        "internal": False
    },
    "layer-taxonomy": {
        "name": "LAYER TAXONOMY",
        "type": "CONCEPT",
        "internal": True
    },
    "layer-taxonomy-txt": {
        "name": "LAYER TAXONOMY",
        "type": "TXT",
        "internal": False
    },
    "gene-taxonomy": {
        "name": "GENE TAXONOMY",
        "type": "CONCEPT",
        "internal": True
    },
    "gene-taxonomy-txt": {
        "name": "GENE TAXONOMY",
        "type": "TXT",
        "internal": False
    },
    "genome-taxonomy": {
        "name": "GENOME TAXONOMY",
        "type": "CONCEPT",
        "internal": True
    },
    "genome-taxonomy-txt": {
        "name": "GENOME TAXONOMY",
        "type": "TXT",
        "internal": False
    },
    "scgs-taxonomy-db": {
        "name": "SCG TAXONOMY DB",
        "type": "CONCEPT",
        "internal": True
    },
    "scgs-taxonomy": {
        "name": "SCG TAXONOMY",
        "type": "CONCEPT",
        "internal": True
    },
    "external-genomes": {
        "name": "EXTERNAL GENOMES",
        "type": "TXT",
        "internal": False
    },
    "internal-genomes": {
        "name": "INTERNAL GENOMES",
        "type": "TXT",
        "internal": False
    },
    "coverages-txt": {
        "name": "COVERAGES",
        "type": "TXT",
        "internal": False
    },
    "genome-distance-txt": {
        "name": "DISTANCE ESTIMATES",
        "type": "TXT",
        "internal": False
    },
    "genome-distance": {
        "name": "DISTANCE ESTIMATES",
        "type": "CONCEPT",
        "internal": True
    },
    "variability-profile": {
        "name": "VARIABILITY PROFILE",
        "type": "CONCEPT",
        "internal": False
    },
    "codon-frequencies-txt": {
        "name": "CODON FREQUENCIES",
        "type": "TXT",
        "internal": False
    },
    "aa-frequencies-txt": {
        "name": "AA FREQUENCIES",
        "type": "TXT",
        "internal": False
    },
    "fixation-index-matrix": {
        "name": "FIXATION INDEX MATRIX",
        "type": "TXT",
        "internal": False
    },
    "summary": {
        "name": "STATIC SUMMARY",
        "type": "SUMMARY",
        "internal": False
    },
    "split-bins": {
        "name": "SPLIT BINS",
        "type": "CONCEPT",
        "internal": False
    },
    "state": {
        "name": "INTERACTIVE STATE",
        "type": "CONCEPT",
        "internal": True
    },
    "ngrams": {
        "name": "NGRAM",
        "type": "CONCEPT",
        "internal": True
    },
    "state-json": {
        "name": "INTERACTIVE STATE",
        "type": "JSON",
        "internal": False
    },
    "kegg-metabolism": {
        "name": "KEGG METABOLISM ESTIMATES",
        "type": "TXT",
        "internal": False
    }
}
