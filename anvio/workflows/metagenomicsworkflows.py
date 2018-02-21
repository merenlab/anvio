# This dictionary defines the parameters that could be changed for each rule
ALL_RULES = ['all', 'gen_configs', 'qc', 'gen_qc_report', 'gzip_fastqs',\
             'fq2fa', 'merge_fastas_for_co_assembly', 'megahit', 'reformat_fasta',\
             'gen_contigs_db', 'export_gene_calls', 'run_centrifuge',\
             'import_taxonomy', 'anvi_run_hmms', 'anvi_run_ncbi_cogs',\
             'bowtie_build', 'bowtie', 'samtools_view', 'anvi_init_bam',\
             'anvi_profile', 'annotate_contigs_database', 'anvi_merge']

accessible_params_dict = {}

# defining the accesible params per rule
accessible_params_dict['all'] = []
accessible_params_dict['gen_configs'] = []
accessible_params_dict['qc'] = []
accessible_params_dict['gen_qc_report'] = []
accessible_params_dict['gzip_fastqs'] = []
accessible_params_dict['fq2fa'] = []
accessible_params_dict['merge_fastas_for_co_assembly'] = []
accessible_params_dict['megahit'] = []
accessible_params_dict['reformat_fasta'] = []
accessible_params_dict['gen_contigs_db'] = []
accessible_params_dict['export_gene_calls'] = []
accessible_params_dict['run_centrifuge'] = []
accessible_params_dict['import_taxonomy'] = []
accessible_params_dict['anvi_run_hmms'] = []
accessible_params_dict['anvi_run_ncbi_cogs'] = []
accessible_params_dict['bowtie_build'] = []
accessible_params_dict['bowtie'] = []
accessible_params_dict['samtools_view'] = []
accessible_params_dict['anvi_init_bam'] = []
accessible_params_dict['anvi_profile'] = []
accessible_params_dict['annotate_contigs_database'] = []
accessible_params_dict['anvi_merge'] = []

