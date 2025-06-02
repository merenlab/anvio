
rule make_anvio_state_file_tree:
    """Make a state file customized for EcoPhylo workflow interactive interface - TREE MODE"""

    version: 1.0
    log: os.path.join(dirs_dict['LOGS_DIR'], "make_anvio_state_file_{hmm_source}_{hmm_name}.log")
    input:
        M.get_target_files_make_anvio_state_file("{hmm_source}", "{hmm_name}")
    params:
        tax_data_final = os.path.join(dirs_dict['MISC_DATA'], "{hmm_source}", "{hmm_name}", "{hmm_name}_scg_taxonomy_data.tsv"),
        misc_data_final = os.path.join(dirs_dict['MISC_DATA'], "{hmm_source}", "{hmm_name}", "{hmm_name}_misc.tsv")
    output:
        state_file = os.path.join(dirs_dict['MISC_DATA'], "{hmm_source}", "{hmm_name}", "{hmm_name}_TREE_state.json")
    threads: M.T('make_anvio_state_file')
    run:

        # Read in misc data headers for layer_order
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

        layer_order = first_layers + misc_layers_list

        #if anvi-estimate-scg-taxonomy was run:
        if os.path.isfile(params.tax_data_final):
            with open(params.tax_data_final) as f:
                lines = f.read()
                first = lines.split('\n', 1)[0]
                scg_taxonomy_layers_list = first.split("\t")
            layer_order.extend(scg_taxonomy_layers_list)

        state_dict['layer-order'] = layer_order

        # layers
        layers_dict = {}

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

        names = {
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

        layers_dict['__parent__'] = layer_attributes_parent
        layers_dict['length'] = length
        layers_dict['gc_content'] = gc_content
        layers_dict['identifier'] = identifier
        layers_dict['percent_identity'] = percent_identity
        layers_dict['names'] = names

        state_dict['layers'] = layers_dict

        # views
        views_dict = {}

        single_dict = {}

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

        single_dict['percent_identity'] = percent_identity
        views_dict['single'] = single_dict
        state_dict['views'] = views_dict

        with open(output.state_file, "w") as outfile:
                json.dump(state_dict, outfile, indent=4)

rule anvi_import_everything_tree:
    """Import state file, phylogenetic tree, AND misc data to interactive interface
    If samples.txt is NOT provided then we will make an Ad Hoc profileDB for the tree to import misc data
    """

    version: 1.0
    log: os.path.join(dirs_dict['LOGS_DIR'], "anvi_import_state_{hmm_source}_{hmm_name}.log")
    input:
        tree = rules.rename_tree_tips.output.tree,
        state = rules.make_anvio_state_file_tree.output.state_file
    params:
        misc_data = rules.make_misc_data.output.misc_data_final,
        tax_data_final = rules.anvi_estimate_scg_taxonomy.params.tax_data_final,
        tree_profileDB = os.path.join(dirs_dict['TREES'], "{hmm_source}", "{hmm_name}", "{hmm_name}-PROFILE.db")
    output:
        touch(os.path.join(dirs_dict['TREES'], "{hmm_source}", "{hmm_name}", "state_imported_tree.done"))

    threads: M.T('anvi_import_state')
    run:
        # Make place holder profileDB for tree
        import anvio
        import anvio.utils as utils
        import anvio.filesnpaths as filesnpaths
        from anvio.dbops import ProfileSuperclass, ProfileDatabase

        filesnpaths.is_file_exists(input.tree)
        newick_tree_text = ''.join([l.strip() for l in open(os.path.abspath(input.tree)).readlines()])

        p_meta = {}
        views = {}

        p_meta['output_dir'] = None
        p_meta['views'] = {}
        p_meta['db_type'] = 'profile'
        p_meta['merged'] = True
        p_meta['blank'] = True
        p_meta['default_view'] = 'single'

        p_meta['item_orders'] = {}
        p_meta['available_item_orders'] = []
        p_meta['default_item_order'] = []

        item_order_name = '%s:unknown:unknown' % filesnpaths.get_name_from_file_path(input.tree)
        p_meta['default_item_order'] = item_order_name
        p_meta['available_item_orders'].append(item_order_name)
        p_meta['item_orders'][item_order_name] = {'type': 'newick', 'data': newick_tree_text}

        p_meta['default_view'] = "single"
        p_meta['samples'] = "names"
        p_meta['sample_id'] = 'AD HOC DISPLAY'

        profile_db = ProfileDatabase(params.tree_profileDB)
        profile_db.create({'db_type': 'profile',
                            'db_variant': 'ad-hoc-display',
                            'blank': True,
                            'merged': True,
                            'contigs_db_hash': None,
                            'items_ordered': False,
                            'samples': ', '.join(p_meta['samples']),
                            'sample_id': p_meta['sample_id']})

        shell("anvi-import-items-order -p {params.tree_profileDB} -i {input.tree} --name {wildcards.hmm_name}_tree >> {log} 2>&1")

        shell("anvi-import-misc-data -p {params.tree_profileDB} --target-data-table items {params.misc_data} --just-do-it >> {log} 2>&1")

        shell("anvi-import-state -p {params.tree_profileDB} -s {input.state} -n default >> {log} 2>&1")

        if os.path.isfile(params.tax_data_final):
            shell("anvi-import-misc-data -p {params.tree_profileDB} --target-data-table items {params.tax_data_final} --just-do-it >> {log} 2>&1")
