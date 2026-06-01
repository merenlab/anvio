window.ANVIO_PRELOADED_WORKFLOWS = [
  {
    id: "contigs",
    title: "contigs workflow",
    source: "merenlab/anvio anvio/workflows/contigs/Snakefile",
    description: "From FASTA files to annotated anvi'o contigs databases.",
    steps: [
      { id: "reformat-prefix", rule: "anvi_script_reformat_fasta_prefix_only", command: "anvi-script-reformat-fasta", description: "Reformat FASTA headers with the sample/group prefix.", inputs: ["contigs-fasta"], outputs: ["contigs-fasta"], params: ["--prefix", "--keep-ids", "--exclude-ids", "--simplify-names", "--seq-type", "--report-file"] },
      { id: "reformat-length", rule: "anvi_script_reformat_fasta", command: "anvi-script-reformat-fasta", description: "Apply contig length filters after prefix formatting.", inputs: ["contigs-fasta"], outputs: ["contigs-fasta"], params: ["--min-len", "--max-len", "--seq-type"] },
      { id: "external-gene-calls", rule: "reformat_external_gene_calls_table", command: "reformat_external_gene_calls_table", description: "Reformat external gene calls so contig names match the reformatted FASTA.", inputs: ["external-gene-calls", "contigs-fasta"], outputs: ["external-gene-calls"], params: [] },
      { id: "contigs-db", rule: "anvi_gen_contigs_database", command: "anvi-gen-contigs-database", description: "Generate the anvi'o contigs database.", inputs: ["contigs-fasta", "external-gene-calls"], outputs: ["contigs-db"], params: ["--description", "--skip-gene-calling", "--external-gene-calls", "--ignore-internal-stop-codons", "--skip-predict-frame", "--skip-mindful-splitting", "--contigs-fasta", "--project-name", "--split-length", "--kmer-size"] },
      { id: "external-functions", rule: "import_external_functions", command: "anvi-import-functions", description: "Import external functional annotations into the contigs database.", inputs: ["contigs-db", "functions-txt"], outputs: ["functions"], params: [] },
      { id: "gene-calls", rule: "export_gene_calls_for_centrifuge", command: "anvi-get-sequences-for-gene-calls", description: "Export gene calls for downstream taxonomy with centrifuge.", inputs: ["contigs-db"], outputs: ["genes-fasta"], params: ["--get-aa-sequences"] },
      { id: "centrifuge", rule: "centrifuge", command: "centrifuge", description: "Run centrifuge on exported gene calls.", inputs: ["genes-fasta"], outputs: ["gene-taxonomy-txt"], params: ["-x", "--threads", "--report-file"] },
      { id: "gene-taxonomy", rule: "anvi_import_taxonomy_for_genes", command: "anvi-import-taxonomy-for-genes", description: "Import centrifuge taxonomy calls for genes.", inputs: ["contigs-db", "gene-taxonomy-txt"], outputs: ["gene-taxonomy"], params: ["--parser"] },
      { id: "hmms", rule: "anvi_run_hmms", command: "anvi-run-hmms", description: "Run HMM searches on the contigs database.", inputs: ["contigs-db"], outputs: ["hmm-hits"], params: ["--installed-hmm-profile", "--hmm-profile-dir", "--also-scan-trnas", "--add-to-functions-table", "--just-do-it"] },
      { id: "pfams", rule: "anvi_run_pfams", command: "anvi-run-pfams", description: "Annotate genes with Pfam hits.", inputs: ["contigs-db"], outputs: ["pfams-data"], params: ["--pfam-data-dir"] },
      { id: "kofams", rule: "anvi_run_kegg_kofams", command: "anvi-run-kegg-kofams", description: "Annotate genes with KEGG KOfam hits.", inputs: ["contigs-db"], outputs: ["kegg-functions"], params: ["--kegg-data-dir", "--hmmer-program", "--keep-all-hits", "--log-bitscores", "--just-do-it"] },
      { id: "cogs", rule: "anvi_run_ncbi_cogs", command: "anvi-run-ncbi-cogs", description: "Annotate genes with NCBI COGs.", inputs: ["contigs-db"], outputs: ["cogs-data"], params: ["--cog-data-dir", "--temporary-dir-path", "--search-with"] },
      { id: "scg-taxonomy", rule: "anvi_run_scg_taxonomy", command: "anvi-run-scg-taxonomy", description: "Populate the contigs database with SCG taxonomy.", inputs: ["contigs-db", "hmm-hits"], outputs: ["scgs-taxonomy"], params: ["--scgs-taxonomy-data-dir"] },
      { id: "trna-scan", rule: "anvi_run_trna_scan", command: "anvi-scan-trnas", description: "Scan tRNA genes in the contigs database.", inputs: ["contigs-db"], outputs: ["trna-gene-hits"], params: ["--trna-cutoff-score", "--trna-model"] },
      { id: "aa-sequences", rule: "anvi_get_sequences_for_gene_calls", command: "anvi-get-sequences-for-gene-calls", description: "Export amino-acid sequences for EggNOG mapper.", inputs: ["contigs-db"], outputs: ["genes-fasta"], params: ["--get-aa-sequences"] },
      { id: "emapper", rule: "emapper", command: "emapper.py", description: "Run EggNOG mapper on exported amino-acid sequences.", inputs: ["genes-fasta"], outputs: ["functions-txt"], params: ["--database", "--usemem", "--override"] },
      { id: "eggnog-import", rule: "anvi_script_run_eggnog_mapper", command: "anvi-script-run-eggnog-mapper", description: "Import EggNOG mapper annotations into the contigs database.", inputs: ["contigs-db", "functions-txt"], outputs: ["functions"], params: ["--use-version"] }
    ]
  },
  {
    id: "metagenomics",
    title: "metagenomics workflow",
    source: "merenlab/anvio anvio/workflows/metagenomics/Snakefile",
    description: "From FASTA and/or FASTQ files to anvi'o contigs and profile databases.",
    includes: ["contigs"],
    steps: [
      { id: "iu-configs", rule: "iu_gen_configs", command: "iu-gen-configs", description: "Generate illumina-utils QC configuration files from samples.txt.", inputs: ["samples-txt"], outputs: ["configuration-ini"], params: ["--r1-prefix", "--r2-prefix"] },
      { id: "qc", rule: "iu_filter_quality_minoche", command: "iu-filter-quality-minoche", description: "Quality-filter paired-end short reads.", inputs: ["configuration-ini", "paired-end-fastq"], outputs: ["paired-end-fastq"], params: ["--ignore-deflines", "--visualize-quality-curves", "--limit-num-pairs", "--print-qual-scores", "--store-read-fate"] },
      { id: "remove-refs", rule: "remove_short_reads_based_on_references", command: "iu-remove-ids-from-fastq", description: "Optionally remove reads that map to reference sequences.", inputs: ["paired-end-fastq", "raw-bam-file"], outputs: ["paired-end-fastq"], params: [] },
      { id: "qc-report", rule: "gen_qc_report", command: "gen_qc_report", description: "Generate a QC report from illumina-utils statistics.", inputs: ["paired-end-fastq"], outputs: ["summary"], params: [] },
      { id: "assembly", rule: "megahit/metaspades/idba_ud", command: "megahit", description: "Assemble quality-filtered reads into contigs using the configured assembler.", inputs: ["paired-end-fastq"], outputs: ["contigs-fasta"], params: ["--threads"] },
      { id: "bowtie-build", rule: "bowtie_build", command: "bowtie2-build", description: "Build a Bowtie2 index from contigs.", inputs: ["contigs-fasta"], outputs: ["short-reads-fasta"], params: [] },
      { id: "bowtie", rule: "bowtie", command: "bowtie2", description: "Map quality-filtered reads to contigs.", inputs: ["paired-end-fastq", "short-reads-fasta"], outputs: ["raw-bam-file"], params: ["--threads"] },
      { id: "init-bam", rule: "anvi_init_bam", command: "anvi-init-bam", description: "Sort and index BAM files for anvi'o profiling.", inputs: ["raw-bam-file"], outputs: ["bam-file"], params: [] },
      { id: "profile", rule: "anvi_profile", command: "anvi-profile", description: "Profile each BAM file against the contigs database.", inputs: ["bam-file", "contigs-db"], outputs: ["single-profile-db"], params: ["--sample-name", "--min-contig-length", "--profile-SCVs"] },
      { id: "merge", rule: "anvi_merge", command: "anvi-merge", description: "Merge individual profiles into a merged profile database.", inputs: ["single-profile-db", "contigs-db"], outputs: ["profile-db"], params: ["--sample-name", "--description", "--skip-hierarchical-clustering", "--enforce-hierarchical-clustering"] },
      { id: "cluster", rule: "anvi_cluster_contigs", command: "anvi-cluster-contigs", description: "Optionally run automatic binning drivers on the merged profile.", inputs: ["profile-db", "contigs-db"], outputs: ["collection", "bin"], params: ["--driver", "--collection-name", "--just-do-it"] },
      { id: "summarize", rule: "anvi_summarize", command: "anvi-summarize", description: "Summarize bins/collections from the merged profile.", inputs: ["profile-db", "contigs-db", "collection"], outputs: ["summary"], params: ["--collection-name"] },
      { id: "split", rule: "anvi_split", command: "anvi-split", description: "Optionally split profile databases by collection.", inputs: ["profile-db", "contigs-db", "collection"], outputs: ["split-bins"], params: [] },
      { id: "krakenuniq", rule: "krakenuniq", command: "krakenuniq", description: "Optionally profile short-read taxonomy with KrakenUniq.", inputs: ["paired-end-fastq"], outputs: ["layer-taxonomy-txt"], params: ["--db", "--threads", "--paired"] },
      { id: "import-layer-taxonomy", rule: "import_krakenuniq_taxonomy", command: "anvi-import-taxonomy-for-layers", description: "Import KrakenUniq taxonomy into profile layers.", inputs: ["profile-db", "layer-taxonomy-txt"], outputs: ["layer-taxonomy"], params: ["--parser", "--min-abundance"] }
    ]
  },
  {
    id: "pangenomics",
    title: "pangenomics workflow",
    source: "merenlab/anvio anvio/workflows/pangenomics/Snakefile",
    description: "Generate a genomes storage, pangenome, optional genome similarity data, and pangenomic phylogeny.",
    includes: ["contigs", "phylogenomics"],
    steps: [
      { id: "genomes-storage", rule: "anvi_gen_genomes_storage", command: "anvi-gen-genomes-storage", description: "Generate an anvi'o genomes storage database.", inputs: ["external-genomes", "internal-genomes", "contigs-db"], outputs: ["genomes-storage-db"], params: ["--external-genomes", "--internal-genomes", "--gene-caller"] },
      { id: "pan", rule: "anvi_pan_genome", command: "anvi-pan-genome", description: "Run the pangenome analysis.", inputs: ["genomes-storage-db"], outputs: ["pan-db"], params: ["--project-name", "--genome-names", "--skip-alignments", "--align-with", "--mcl-inflation", "--minbit"] },
      { id: "gc-seqs", rule: "anvi_get_sequences_for_gene_clusters", command: "anvi-get-sequences-for-gene-clusters", description: "Export gene-cluster sequences from the pangenome.", inputs: ["pan-db", "genomes-storage-db"], outputs: ["genes-fasta"], params: ["--gene-cluster-id", "--collection-name", "--concatenate-gene-clusters", "--align-with"] },
      { id: "similarity", rule: "anvi_compute_genome_similarity", command: "anvi-compute-genome-similarity", description: "Compute genome similarity metrics and optionally store them in the pan database.", inputs: ["external-genomes", "internal-genomes", "pan-db"], outputs: ["genome-similarity"], params: ["--external-genomes", "--internal-genomes"] },
      { id: "trimal", rule: "trimal", command: "trimal", description: "Trim pangenomic alignments.", inputs: ["genes-fasta"], outputs: ["concatenated-gene-alignment-fasta"], params: ["-gt"] },
      { id: "iqtree", rule: "iqtree", command: "iqtree", description: "Compute a pangenomic phylogenetic tree.", inputs: ["concatenated-gene-alignment-fasta"], outputs: ["phylogeny"], params: ["-m", "-bb"] },
      { id: "import-tree", rule: "import_phylogenetic_tree_to_pangenome", command: "anvi-import-misc-data", description: "Import the pangenomic tree as layer ordering data.", inputs: ["pan-db", "phylogeny"], outputs: ["misc-data-layer-orders"], params: ["--target-data-table"] }
    ]
  },
  {
    id: "phylogenomics",
    title: "phylogenomics workflow",
    source: "merenlab/anvio anvio/workflows/phylogenomics/Snakefile",
    description: "Export HMM-hit sequences, trim alignments, and infer a phylogenomic tree.",
    steps: [
      { id: "hmm-seqs", rule: "anvi_get_sequences_for_hmm_hits", command: "anvi-get-sequences-for-hmm-hits", description: "Export sequences for selected HMM hits.", inputs: ["external-genomes", "internal-genomes", "hmm-list"], outputs: ["genes-fasta"], params: ["--internal-genomes", "--external-genomes", "--return-best-hit", "--separator", "--align-with", "--concatenate-genes", "--get-aa-sequences", "--gene-names", "--hmm-sources"] },
      { id: "trimal", rule: "trimal", command: "trimal", description: "Remove gaps from the concatenated protein alignment.", inputs: ["genes-fasta"], outputs: ["concatenated-gene-alignment-fasta"], params: ["-gt"] },
      { id: "iqtree", rule: "iqtree", command: "iqtree", description: "Infer the phylogenomic tree with IQ-TREE.", inputs: ["concatenated-gene-alignment-fasta"], outputs: ["phylogeny"], params: ["-m", "-bb"] }
    ]
  },
  {
    id: "sra-download",
    title: "sra-download workflow",
    source: "merenlab/anvio anvio/workflows/sra_download/Snakefile",
    description: "Download, verify, extract, compress, and tabulate SRA reads.",
    steps: [
      { id: "prefetch", rule: "prefetch", command: "prefetch", description: "Fetch SRA accessions with the SRA toolkit.", inputs: ["sra-accessions"], outputs: ["sra-file"], params: ["--output-directory", "--max-size"] },
      { id: "md5", rule: "check_md5sum", command: "curl", description: "Fetch SRA metadata and verify MD5 checksums.", inputs: ["sra-file"], outputs: ["completion"], params: [] },
      { id: "dump", rule: "fasterq_dump", command: "fasterq-dump", description: "Extract FASTQ files from prefetched SRA files.", inputs: ["sra-file"], outputs: ["paired-end-fastq", "single-end-fastq"], params: ["--split-3", "--threads"] },
      { id: "pigz", rule: "pigz", command: "pigz", description: "Compress FASTQ files with parallel gzip.", inputs: ["paired-end-fastq", "single-end-fastq"], outputs: ["paired-end-fastq", "single-end-fastq"], params: ["--processes"] },
      { id: "samples", rule: "generate_samples_txt", command: "generate_samples_txt", description: "Generate samples.txt and samples_single_reads.txt descriptors.", inputs: ["paired-end-fastq", "single-end-fastq"], outputs: ["samples-txt"], params: [] }
    ]
  },
  {
    id: "trnaseq",
    title: "trnaseq workflow",
    source: "merenlab/anvio anvio/workflows/trnaseq/Snakefile",
    description: "Process tRNA-seq reads, build per-sample tRNA databases, merge them, and optionally assign tRNA taxonomy.",
    steps: [
      { id: "iu-input", rule: "make_iu_input", command: "make_iu_input", description: "Create the illumina-utils input table for tRNA-seq samples.", inputs: ["samples-txt"], outputs: ["samples-txt"], params: [] },
      { id: "iu-configs", rule: "iu_gen_configs", command: "iu-gen-configs", description: "Generate illumina-utils configuration files for paired-end tRNA-seq reads.", inputs: ["samples-txt"], outputs: ["configuration-ini"], params: ["--r1-prefix", "--r2-prefix"] },
      { id: "merge-pairs", rule: "iu_merge_pairs", command: "iu-merge-pairs", description: "Merge paired-end tRNA-seq reads before anvi'o tRNA profiling.", inputs: ["configuration-ini", "paired-end-fastq"], outputs: ["trnaseq-fasta"], params: ["--marker-gene-stringent", "--max-num-mismatches", "--report-r1-prefix", "--report-r2-prefix", "--num-threads"] },
      { id: "qc-report", rule: "gen_qc_report", command: "gen_qc_report", description: "Generate a quality-control report for merged tRNA-seq reads.", inputs: ["trnaseq-fasta"], outputs: ["summary"], params: [] },
      { id: "reformat-fasta", rule: "anvi_reformat_fasta", command: "anvi-script-reformat-fasta", description: "Simplify and reformat merged tRNA-seq FASTA records.", inputs: ["trnaseq-fasta"], outputs: ["trnaseq-fasta"], params: ["--simplify-names", "--report-file"] },
      { id: "profile-sample", rule: "anvi_trnaseq", command: "anvi-trnaseq", description: "Profile each tRNA-seq sample and create an anvi'o tRNA-seq database.", inputs: ["trnaseq-fasta"], outputs: ["trnaseq-db"], params: ["--trnaseq-fasta", "--sample-name", "--output-dir", "--treatment", "--overwrite-output-destinations", "--description", "--write-checkpoints", "--load-checkpoint", "--feature-param-file", "--threeprime-termini", "--min-length-long-fiveprime", "--min-trna-fragment-size", "--agglomeration-max-mismatch-freq", "--skip-INDEL-profiling", "--max-indel-freq", "--left-indel-buffer", "--right-indel-buffer", "--skip-fasta-check", "--profiling-chunk-size", "--alignment-target-chunk-size"] },
      { id: "merge-trnaseq", rule: "anvi_merge_trnaseq", command: "anvi-merge-trnaseq", description: "Merge tRNA-seq sample databases into contigs and profile databases.", inputs: ["trnaseq-db"], outputs: ["trnaseq-contigs-db", "trnaseq-profile-db"], params: ["--project-name", "--max-reported-trna-seeds", "--overwrite-output-destinations", "--description", "--feature-threshold", "--preferred-treatment", "--nonspecific-output", "--min-variation", "--min-third-fourth-nt", "--min-indel-fraction", "--distance", "--linkage"] },
      { id: "taxonomy", rule: "anvi_run_trna_taxonomy", command: "anvi-run-trna-taxonomy", description: "Assign taxonomy to tRNA sequences in the merged contigs database.", inputs: ["trnaseq-contigs-db"], outputs: ["trna-taxonomy"], params: ["--min-percent-identity", "--max-num-target-sequences", "--write-buffer-size", "--all-hits-output-file"] },
      { id: "tabulate", rule: "anvi_tabulate_trnaseq", command: "anvi-tabulate-trnaseq", description: "Tabulate tRNA-seq outputs for downstream inspection.", inputs: ["trnaseq-contigs-db", "trnaseq-profile-db"], outputs: ["summary"], params: ["--output-dir", "--overwrite-output-destinations"] }
    ]
  },
  {
    id: "ecophylo",
    title: "ecophylo workflow",
    source: "merenlab/anvio anvio/workflows/ecophylo/Snakefile and rules/*.smk",
    description: "Co-characterize protein biogeography and phylogeny across genomes and metagenomes.",
    includes: ["metagenomics"],
    steps: [
      { id: "hmms", rule: "anvi_run_hmms_hmmsearch", command: "anvi-run-hmms", description: "Run HMM searches for selected targets.", inputs: ["contigs-db", "hmm-source"], outputs: ["hmm-hits"], params: ["--hmmer-program", "--hmmer-output-dir", "--installed-hmm-profile", "--hmm-profile-dir", "--domain-hits-table"] },
      { id: "filter-hmms", rule: "filter_hmm_hits_by_model_coverage", command: "anvi-script-filter-hmm-hits-table", description: "Filter HMM hits by model coverage.", inputs: ["contigs-db", "hmm-hits"], outputs: ["hmm-hits"], params: ["--domain-hits-table", "--hmm-source", "--min-model-coverage", "--filter-out-partial-gene-calls"] },
      { id: "process-hmms", rule: "process_hmm_hits", command: "anvi-get-sequences-for-hmm-hits", description: "Export HMM-hit amino-acid and nucleotide sequences.", inputs: ["contigs-db", "hmm-hits"], outputs: ["genes-fasta", "external-gene-calls"], params: ["--hmm-sources", "--gene-names", "--get-aa-sequences", "--just-do-it"] },
      { id: "combine", rule: "combine_sequence_data", command: "combine_sequence_data", description: "Combine AA/NT FASTA and external gene calls from all samples.", inputs: ["genes-fasta", "external-gene-calls"], outputs: ["genes-fasta", "external-gene-calls"], params: [] },
      { id: "cluster", rule: "cluster_X_percent_sim_mmseqs", command: "mmseqs", description: "Cluster extracted protein sequences.", inputs: ["genes-fasta"], outputs: ["genes-fasta"], params: [] },
      { id: "align", rule: "align_sequences", command: "muscle", description: "Align representative protein sequences.", inputs: ["genes-fasta"], outputs: ["concatenated-gene-alignment-fasta"], params: [] },
      { id: "trim", rule: "trim_alignment", command: "trimal", description: "Trim the alignment.", inputs: ["concatenated-gene-alignment-fasta"], outputs: ["concatenated-gene-alignment-fasta"], params: ["-gt", "-gappyout"] },
      { id: "tree", rule: "iqtree/fasttree", command: "iqtree", description: "Build the EcoPhylo tree.", inputs: ["concatenated-gene-alignment-fasta"], outputs: ["phylogeny"], params: ["-m"] },
      { id: "rename-tree", rule: "rename_tree_tips", command: "rename_tree_tips", description: "Rename tree tips for interactive display compatibility.", inputs: ["phylogeny", "genes-fasta"], outputs: ["phylogeny"], params: [] },
      { id: "fasta-txt", rule: "make_fasta_txt", command: "make_fasta_txt", description: "Make a fasta.txt for the metagenomics subworkflow.", inputs: ["genes-fasta", "external-gene-calls"], outputs: ["fasta-txt"], params: [] },
      { id: "metagenomics-config", rule: "make_metagenomics_config_file", command: "anvi-run-workflow", description: "Generate a metagenomics workflow config for profile mode.", inputs: ["fasta-txt", "samples-txt"], outputs: ["workflow-config"], params: ["--get-default-config"] },
      { id: "run-metagenomics", rule: "run_metagenomics_workflow", command: "anvi-run-workflow", description: "Run the metagenomics workflow to profile EcoPhylo targets.", inputs: ["workflow-config"], outputs: ["profile-db", "contigs-db"], params: ["--workflow", "--config-file"] },
      { id: "summarize", rule: "anvi_summarize", command: "anvi-summarize", description: "Summarize target coverages.", inputs: ["profile-db", "contigs-db"], outputs: ["summary"], params: ["--init-gene-coverages", "--just-do-it"] },
      { id: "state", rule: "make_anvio_state_file", command: "make_anvio_state_file", description: "Generate an anvi'o state file for the EcoPhylo display.", inputs: ["summary", "misc-data-items"], outputs: ["state-json"], params: [] },
      { id: "import-display", rule: "anvi_import_everything_metagenome/tree", command: "anvi-import-state", description: "Import state, tree, and misc data into an interactive display profile.", inputs: ["profile-db", "state-json", "phylogeny", "misc-data-items"], outputs: ["interactive"], params: ["--name"] }
    ]
  }
];
