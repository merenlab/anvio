# profile-mode with read recruitment

rule make_metagenomics_config_file:
    """Make a METAGENOMICS WORKFLOW config.json customized for ECOPHYLO_WORKFLOW - PROFILE MODE"""

    version: 1.0
    log: os.path.join(dirs_dict['LOGS_DIR'], "make_metagenomics_config_file.log")
    input:
        rules.make_fasta_txt.output.fasta_txt
    output:
        config = os.path.join(dirs_dict['HOME'], "METAGENOMICS_WORKFLOW", "metagenomics_config.json")
    threads: M.T('make_metagenomics_config_file')
    run:

        shell('anvi-run-workflow -w metagenomics --get-default-config {output.config}')

        config = open(output.config)
        config_dict = json.load(config)
        config_dict['fasta_txt'] = 'fasta.txt'
        sample_txt_path = M.samples_txt_file
        config_dict['samples_txt'] = sample_txt_path
        config_dict['references_mode'] = True
        config_dict['anvi_run_hmms']['run'] = False
        config_dict["anvi_script_reformat_fasta"]['run'] = False
        config_dict['anvi_run_kegg_kofams']['run'] = False
        config_dict['anvi_run_ncbi_cogs']['run'] = False
        config_dict['anvi_run_scg_taxonomy']['run'] = False
        config_dict['anvi_run_trna_scan']['run'] = False
        config_dict['anvi_run_scg_taxonomy']['run'] = False
        config_dict['iu_filter_quality_minoche']['run'] = False
        config_dict['anvi_profile']['--min-contig-length'] = 0
        config_dict['anvi_profile']['--profile-SCVs'] = True
        config_dict['bowtie']['threads'] = 5
        config_dict['bowtie_build']['threads'] = 5
        config_dict['anvi_gen_contigs_database']['threads'] = 5
        config_dict['anvi_init_bam']['threads'] = 2
        config_dict['anvi_profile']['--profile-SCVs'] = True 

        if M.clusterize_metagenomics_workflow == True:
            config_dict['bowtie']['threads'] = 10
            config_dict['anvi_profile']['threads'] = 10
            config_dict['anvi_merge']['threads'] = 10

        if M.bowtie2_additional_params:
            config_dict['bowtie']['additional_params'] = " ".join(["--no-unal", M.bowtie2_additional_params]) 

        if M.anvi_profile_min_percent_identity:
            config_dict['anvi_profile']['--min-percent-identity'] = M.anvi_profile_min_percent_identity

        with open(output.config, "w") as outfile:
            json.dump(config_dict, outfile, indent=4)


rule run_metagenomics_workflow:
    """Run metagenomics workflow to profile hmm_hits"""

    version: 1.0
    log: "00_LOGS/run_metagenomics_workflow.log"
    input:
        config = rules.make_metagenomics_config_file.output.config,
    output:
        done = touch(os.path.join(dirs_dict['HOME'], "METAGENOMICS_WORKFLOW", "metagenomics_workflow.done"))
    threads: M.T('run_metagenomics_workflow')
    params:
    run:
        metagenomics_workflow_path = os.path.join(dirs_dict['HOME'], "METAGENOMICS_WORKFLOW") 

        # Convert r1 and r2 to absolute paths
        samples_txt_new_path = os.path.join(metagenomics_workflow_path, M.samples_txt_file)
        samples_txt_new = pd.read_csv(M.samples_txt_file, sep='\t', index_col=False)
        samples_txt_new['r1'] = samples_txt_new['r1'].apply(lambda x: os.path.abspath(str(x)))
        samples_txt_new['r2'] = samples_txt_new['r2'].apply(lambda x: os.path.abspath(str(x)))
        samples_txt_new.to_csv(samples_txt_new_path, sep="\t", index=False, header=True)

        log_path = os.path.join(dirs_dict['HOME'], "METAGENOMICS_WORKFLOW", "00_LOGS")
        log_file = os.path.join(log_path, "run_metagenomics_workflow.log")
        log = os.path.join("00_LOGS", "run_metagenomics_workflow.log")
        shell(f"mkdir -p {log_path} && touch {log_file}")

        if M.clusterize_metagenomics_workflow == True:
            # If we are using slurm and clusterize: https://github.com/ekiefl/clusterize
            shell('cd {metagenomics_workflow_path} && anvi-run-workflow -w metagenomics -c metagenomics_config.json --additional-params --cluster \"clusterize -j={{rule}} -o={{log}} -n={{threads}} {M.clusterize_metagenomics_submission_params} -x\" {M.metagenomics_workflow_snakemake_additional_params} --latency-wait 100 --keep-going --rerun-incomplete &> {log} && cd -')
        elif M.metagenomics_workflow_HPC_string:
            # User-defined --cluster string: https://snakemake.readthedocs.io/en/stable/executing/cluster.html
            shell('cd {metagenomics_workflow_path} && anvi-run-workflow -w metagenomics -c metagenomics_config.json --additional-params --cluster \"{M.metagenomics_workflow_HPC_string}\" {M.metagenomics_workflow_snakemake_additional_params} --rerun-incomplete --latency-wait 100 --keep-going &> {log} && cd -')
        else:
            # Running snakemake on local
            shell('cd {metagenomics_workflow_path} && anvi-run-workflow -w metagenomics -c metagenomics_config.json --additional-params {M.metagenomics_workflow_snakemake_additional_params} --rerun-incomplete --latency-wait 100 --keep-going &> {log} && cd -')


rule add_default_collection:
    """Make default collection for profile-db that contains all splits"""

    version: 1.0
    log: os.path.join(dirs_dict['LOGS_DIR'], "add_default_collection_{hmm}.log")
    input: metagenomics_workflow_done = rules.run_metagenomics_workflow.output.done
    params:
        contigsDB = ancient(os.path.join(dirs_dict['HOME'], "METAGENOMICS_WORKFLOW", "03_CONTIGS", "{hmm}.db")),
        profileDB = os.path.join(dirs_dict['HOME'], "METAGENOMICS_WORKFLOW", "06_MERGED", "{hmm}", "PROFILE.db")
    output: done = touch(os.path.join(dirs_dict['HOME'], "METAGENOMICS_WORKFLOW", "{hmm}_add_default_collection.done"))
    threads: M.T('add_default_collection')
    run:
        shell('anvi-script-add-default-collection -c {params.contigsDB} -p {params.profileDB}')


rule anvi_summarize:
    """Get coverage values for hmm_hits"""

    version: 1.0
    log: os.path.join(dirs_dict['LOGS_DIR'], "anvi_summarize_{hmm}.log")
    input: 
        done = rules.add_default_collection.output.done
    params:
        contigsDB = ancient(os.path.join(dirs_dict['HOME'], "METAGENOMICS_WORKFLOW", "03_CONTIGS", "{hmm}-contigs.db")),
        profileDB = os.path.join(dirs_dict['HOME'], "METAGENOMICS_WORKFLOW", "06_MERGED", "{hmm}", "PROFILE.db"),
        output_dir = os.path.join(dirs_dict['HOME'], "METAGENOMICS_WORKFLOW", "07_SUMMARY", "{hmm}")
    output: touch(os.path.join(dirs_dict['HOME'], "METAGENOMICS_WORKFLOW", "07_SUMMARY", "{hmm}_summarize.done"))
    threads: M.T('anvi_summarize')
    run: 
        shell('anvi-summarize -c {params.contigsDB} -p {params.profileDB} -o {params.output_dir} -C DEFAULT --init-gene-coverages --just-do-it >> {log} 2>&1')
        

rule make_anvio_state_file:
    """Make a state file customized for EcoPhylo workflow interactive interface"""

    version: 1.0
    log: os.path.join(dirs_dict['LOGS_DIR'], "make_anvio_state_file_{hmm}.log")
    input:
        num_tree_tips = rules.subset_DNA_reps_with_QCd_AA_reps_for_mapping.output.NT_for_mapping,
        done_scg = rules.anvi_scg_taxonomy.output.done
    params:
        tax_data_final = os.path.join(dirs_dict['MISC_DATA'], "{hmm}_scg_taxonomy_data.tsv"),
        misc_data_final = os.path.join(dirs_dict['MISC_DATA'], "{hmm}_misc.tsv"),
    output:
        state_file = os.path.join(dirs_dict['HOME'], "METAGENOMICS_WORKFLOW", "{hmm}_ECOPHYLO_WORKFLOW_state.json")
    threads: M.T('make_anvio_state_file')
    run:

        hmm_source = M.hmm_dict[wildcards.hmm]['source']

        # Read in misc data headers for layer_order
        if hmm_source in M.internal_hmm_sources:
            with open(params.tax_data_final) as f:
                lines = f.read()
                first = lines.split('\n', 1)[0]
            scg_taxonomy_layers_list = first.split("\t")

        with open(params.misc_data_final) as f:
            lines = f.read()
            first = lines.split('\n', 1)[0]
        misc_layers_list = first.split("\t")

        state_dict = {}

        # basics
        state_dict['version'] = '3'
        state_dict['tree-type'] = 'phylogram'
        state_dict['current-view'] = 'single'

        # height and width
        # FIXME: It's unclear to me how the interactive interface determines
        # height and width of a tree when the input value is 0. There has to 
        # be some kind of calculation to determine the tree shape in the backend
        # of the interface because even after I export a "default" state file
        # the height and width are still "0". However, if you change the height and width
        # values within the interface to "" the tree will disappear. I need to sort this 
        # out eventually to have a clean way of changing the tree shape to 
        # match the dimensions of the number of SCGs vs metagenomes. 
        # num_tree_tips = pd.read_csv(input.num_tree_tips, \
        #                             sep="\t", \
        #                             index_col=None)


        # layer-orders
        first_layers = ["__parent__", "length", "gc_content"]
        metagenomes = []

        for metagenome in M.sample_names_for_mapping_list:
            metagenomes.append(metagenome)

        if hmm_source in M.internal_hmm_sources:
            layer_order = first_layers + metagenomes + misc_layers_list + scg_taxonomy_layers_list 
        else:
            layer_order = first_layers + metagenomes + misc_layers_list 

        state_dict['layer-order'] = layer_order

        # layers
        layers_dict = {}

        metagenome_layers_dict = {}

        metagenome_attributes = {
            "color": "#000000",
            "height": "180",
            "margin": "15",
            "type": "bar",
            "color-start": "#FFFFFF"
            }

        for metagenome in metagenomes:
            metagenome_layers_dict[str(metagenome)] = metagenome_attributes

        layer_attributes_parent = {
            "color": "#000000",
            "height": "0",
            "margin": "15",
            "type": "color",
            "color-start": "#FFFFFF"
            }

        length = {
            "color": "#000000",
            "height": "0",
            "margin": "15",
            "type": "color",
            "color-start": "#FFFFFF"
            }

        gc_content = {
            "color": "#000000",
            "height": "0",
            "margin": "15",
            "type": "color",
            "color-start": "#FFFFFF"
            }

        identifier = {
            "color": "#000000",
            "height": "0",
            "margin": "15",
            "type": "color",
            "color-start": "#FFFFFF"
            }

        percent_identity = {
            "color": "#000000",
            "height": "180",
            "margin": "15",
            "type": "line",
            "color-start": "#FFFFFF"
            }
            
        layers_dict.update(metagenome_layers_dict)
        layers_dict['__parent__'] = layer_attributes_parent
        layers_dict['length'] = length
        layers_dict['gc_content'] = gc_content
        layers_dict['identifier'] = identifier
        layers_dict['percent_identity'] = percent_identity

        state_dict['layers'] = layers_dict

        # views
        views_dict = {}

        single_dict = {}

        mean_coverage_dict = {}

        percent_identity = {
            "normalization": "none",
            "min": {
                "value": "90",
                "disabled": "false"
                },
            "max": {
                "value": "100",
                "disabled": "false"
                }
        }

        cluster_size = {
            "normalization": "none"
        }

        single_dict['percent_identity'] = percent_identity 
        mean_coverage_dict['percent_identity'] = percent_identity 
        mean_coverage_dict['cluster_size'] = cluster_size 
        views_dict['single'] = single_dict
        views_dict['mean_coverage'] = mean_coverage_dict
        state_dict['views'] = views_dict

        # samples-layer-order
        samples_layers_dict = {
            "default": {
                "num_INDELs_reported": {
                    "height": 0,
                },
                "total_reads_kept": {
                    "height": 0,
                },
                "num_SCVs_reported": {
                    "height": 0,
                },
                "num_SNVs_reported": {
                    "height": 0,
                },
                "total_reads_mapped": {
                    "height": 0,
                },
            }
        }
        
        state_dict['samples-layers'] = samples_layers_dict

        with open(output.state_file, "w") as outfile:
                json.dump(state_dict, outfile, indent=4)

rule anvi_import_everything_metagenome:
    """Import state file, phylogenetic tree, AND misc data to interactive interface
    If samples.txt is NOT provided then we will make an Ad Hoc profileDB for the tree to import misc data
    """

    version: 1.0
    log: os.path.join(dirs_dict['LOGS_DIR'], "anvi_import_state_{hmm}.log")
    input:
        tree = rules.rename_tree_tips.output.tree,
        misc_data = rules.make_misc_data.output.misc_data_final,
        state = rules.make_anvio_state_file.output.state_file,
        done = rules.run_metagenomics_workflow.output.done
    params:
        tax_data_final = rules.anvi_scg_taxonomy.params.tax_data_final,
        profileDB = os.path.join(dirs_dict['HOME'], "METAGENOMICS_WORKFLOW", "06_MERGED", "{hmm}", "PROFILE.db"),
        tree_profileDB = os.path.join(dirs_dict['TREES'], "{hmm}", "{hmm}-PROFILE.db")
    output: 
        touch(os.path.join(dirs_dict['HOME'], "METAGENOMICS_WORKFLOW", "{hmm}_state_imported_profile.done"))
    threads: M.T('anvi_import_state')
    run:
        state = os.path.join(dirs_dict['HOME'], "METAGENOMICS_WORKFLOW", f"{wildcards.hmm}_ECOPHYLO_WORKFLOW_state.json")

        shell("echo -e 'Step 1: anvi-import-state:\n' >> {log}")
        shell("anvi-import-state -p {params.profileDB} -s {state} -n default >> {log} 2>&1")
        shell("echo -e '' >> {log}")

        shell("echo -e 'Step 2: anvi-import-items-order:\n' >> {log}")
        shell("anvi-import-items-order -p {params.profileDB} -i {input.tree} --name {wildcards.hmm}_tree >> {log} 2>&1")
        shell("echo -e '' >> {log}")

        shell("echo -e 'Step 3: anvi-import-misc-data:\n' >> {log}")
        shell("anvi-import-misc-data -p {params.profileDB} --target-data-table items {input.misc_data} --just-do-it >> {log} 2>&1")
        shell("echo -e '' >> {log}")

        hmm_source = M.hmm_dict[wildcards.hmm]['source']
        
        if hmm_source in M.internal_hmm_sources:
            shell("echo -e 'Step 4: anvi-import-misc-data:\n' >> {log}")
            shell("anvi-import-misc-data -p {params.profileDB} --target-data-table items {params.tax_data_final} --just-do-it >> {log} 2>&1")
            shell("echo -e '' >> {log}")
