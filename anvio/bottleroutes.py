# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Common routes for bottle web server.

    Functions defined here are used from wihtin programs such as
    anvi-interactive, or anvi-refine.
"""

import os
import json
import datetime

from bottle import static_file

import anvio
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.summarizer as summarizer

from anvio.errors import RefineError


__author__ = "Ozcan Esen"
__copyright__ = "Copyright 2016, The anvio Project"
__credits__ = ["A. Murat Eren"]
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


run = terminal.Run()
progress = terminal.Progress()


def set_default_headers(response):
    response.set_header('Content-Type', 'application/json')
    response.set_header('Pragma', 'no-cache')
    response.set_header('Cache-Control', 'no-cache, no-store, max-age=0, must-revalidate')
    response.set_header('Expires', 'Thu, 01 Dec 1994 16:00:00 GMT')


def state_autoload(d, response):
    # see --state parameter.
    set_default_headers(response)

    return json.dumps(d.state_autoload)


def state_all(d, response):
    set_default_headers(response)

    return json.dumps(d.states_table.states)


def save_state(args, d, request, response):
    if args.read_only:
        return json.dumps({'status_code': '0'})

    name = request.forms.get('name')
    content = request.forms.get('content')
    last_modified = datetime.datetime.now().strftime("%d.%m.%Y %H:%M:%S")

    d.states_table.store_state(name, content, last_modified)

    return json.dumps({'status_code': '1'})


def get_state(d, request, response):
    set_default_headers(response)

    name = request.forms.get('name')

    if name in d.states_table.states:
        state = d.states_table.states[name]
        return json.dumps(state['content'])

    return json.dumps("")


def get_protein_clusters_summary(d, request):
    protein_cluster_ids = json.loads(request.forms.get('split_names'))
    bin_name = json.loads(request.forms.get('bin_name'))

    summary = d.get_summary_for_PCs_list(protein_cluster_ids)

    run.info_single('PC info has been requested for %d items in %s' % (len(protein_cluster_ids), bin_name))

    return json.dumps({'functions': summary['functions'],
                       'num_PCs': summary['num_PCs'],
                       'genomes_contributing': summary['genomes_contributing'],
                       'num_gene_calls': summary['num_gene_calls']})


def completeness(d, request):
    completeness_sources = {}
    completeness_averages = {}
    if not d.completeness:
        return json.dumps(completeness_sources)

    split_names = json.loads(request.forms.get('split_names'))
    bin_name = json.loads(request.forms.get('bin_name'))

    run.info_single('Completeness info has been requested for %d splits in %s' % (len(split_names), bin_name))

    p_completion, p_redundancy, scg_domain, domain_confidence, results_dict = d.completeness.get_info_for_splits(set(split_names))

    # convert results_dict (where domains are the highest order items) into a dict that is compatible with the
    # previous format of the dict (where hmm scg source names are the higher order items).
    for domain in results_dict:
        for source in results_dict[domain]:
            completeness_sources[source] = results_dict[domain][source]

    completeness_averages['percent_completion'] = p_completion
    completeness_averages['percent_redundancy'] = p_redundancy
    completeness_averages['domain'] = scg_domain
    completeness_averages['domain_confidence'] = domain_confidence

    return json.dumps({'stats': completeness_sources, 'averages': completeness_averages, 'refs': d.completeness.http_refs})


def get_index_total_previous_and_next_items(d, item_name):
    previous_item_name = None
    next_item_name = None
    index = None
    total = None

    index_of_item = d.displayed_item_names_ordered.index(item_name)
    if index_of_item:
        previous_item_name = d.displayed_item_names_ordered[index_of_item - 1]
    if (index_of_item + 1) < len(d.displayed_item_names_ordered):
        next_item_name = d.displayed_item_names_ordered[index_of_item + 1]

    index = index_of_item + 1
    total = len(d.displayed_item_names_ordered)

    return index, total, previous_item_name, next_item_name 


def get_AA_sequences_for_PC(d, pc_name):
    data = {}

    if pc_name not in d.protein_clusters:
        return data

    if not d.genomes_storage_is_available:
        return data

    # add the list of gene caller ids associated with this protein cluster into `data`:
    for genome_name in d.protein_clusters[pc_name]:
        for gene_callers_id in d.protein_clusters[pc_name][genome_name]:
            data['%s_%s' % (genome_name, str(gene_callers_id))] = d.genomes_storage.get_gene_sequence(genome_name, gene_callers_id)

    return json.dumps(data)


def inspect_pc(d, pc_name):
    data = {'pc_name': pc_name,
            'genomes': [],
            'index': None,
            'gene_caller_ids': [],
            'gene_caller_ids_in_genomes': {},
            'gene_aa_sequences_for_gene_caller_ids': {},
            'previous_pc_name': None,
            'next_pc_name': None,
            'index': None,
            'total': None
            }

    if pc_name not in d.protein_clusters:
        return data

    if not d.genomes_storage_is_available:
        return data

    # add the list of gene caller ids associated with this protein cluster into `data`:
    for genome_name in d.protein_clusters[pc_name]:
        for gene_callers_id in d.protein_clusters[pc_name][genome_name]:
            data['gene_caller_ids'].append((gene_callers_id, genome_name), )

    # the dict that explains the distribution of genes in genomes:
    data['gene_caller_ids_in_genomes'] = d.protein_clusters[pc_name]

    # add the list of genomes into data:
    data['genomes'] = sorted(data['gene_caller_ids_in_genomes'].keys())

    # get some contextual stuff
    data['index'], data['total'], data['previous_pc_name'], data['next_pc_name'] = get_index_total_previous_and_next_items(d, pc_name)

    for gene_callers_id, genome_name in data['gene_caller_ids']:
        sequence = d.genomes_storage.get_gene_sequence(genome_name, gene_callers_id)

        # restore alignment from the alignment summary if possible
        if d.protein_clusters_gene_alignments_available:
            alignment_summary = d.protein_clusters_gene_alignments[genome_name][gene_callers_id]
            sequence = utils.restore_alignment(sequence, alignment_summary) if alignment_summary else sequence

        data['gene_aa_sequences_for_gene_caller_ids'][gene_callers_id] = sequence

    return json.dumps(data)


def charts(d, split_name, show_outlier_SNVs=False):
    data = {'layers': [],
             'index': None,
             'total': None,
             'coverage': [],
             'variability': [],
             'competing_nucleotides': [],
             'previous_contig_name': None,
             'next_contig_name': None,
             'genes': [],
             'outlier_SNVs_shown': show_outlier_SNVs}

    if split_name not in d.split_names:
        return data

    if not d.auxiliary_profile_data_available:
        return data

    data['index'], data['total'], data['previous_contig_name'], data['next_contig_name'] = get_index_total_previous_and_next_items(d, split_name)

    layers = sorted(d.p_meta['samples'])

    coverage_values_dict = d.split_coverage_values.get(split_name)
    data['coverage'] = [coverage_values_dict[layer].tolist() for layer in layers]

    ## get the variability information dict for split:
    progress.new('Variability')
    progress.update('Collecting info for "%s"' % split_name)
    split_variability_info_dict = d.get_variability_information_for_split(split_name, return_outliers=show_outlier_SNVs)

    for layer in layers:
        progress.update('Formatting variability data: "%s"' % layer)
        data['layers'].append(layer)
        data['competing_nucleotides'].append(split_variability_info_dict[layer]['competing_nucleotides'])
        data['variability'].append(split_variability_info_dict[layer]['variability'])

    levels_occupied = {1: []}
    for entry_id in d.split_name_to_gene_caller_ids_dict[split_name]:
        gene_callers_id = d.genes_in_splits[entry_id]['gene_callers_id']
        p = d.genes_in_splits[entry_id]
        # p looks like this at this point:
        #
        # {'percentage_in_split': 100,
        #  'start_in_split'     : 16049,
        #  'stop_in_split'      : 16633}
        #  'prot'               : u'prot2_03215',
        #  'split'              : u'D23-1contig18_split_00036'}
        #
        # we will add a bit more attributes:
        p['source'] = d.genes_in_contigs_dict[gene_callers_id]['source']
        p['direction'] = d.genes_in_contigs_dict[gene_callers_id]['direction']
        p['start_in_contig'] = d.genes_in_contigs_dict[gene_callers_id]['start']
        p['stop_in_contig'] = d.genes_in_contigs_dict[gene_callers_id]['stop']
        p['complete_gene_call'] = 'No' if d.genes_in_contigs_dict[gene_callers_id]['partial'] else 'Yes'
        p['length'] = p['stop_in_contig'] - p['start_in_contig']
        p['functions'] = d.gene_function_calls_dict[gene_callers_id] if gene_callers_id in d.gene_function_calls_dict else None

        for level in levels_occupied:
            level_ok = True
            for gene_tuple in levels_occupied[level]:
                if (p['start_in_split'] >= gene_tuple[0] - 100 and p['start_in_split'] <= gene_tuple[1] + 100) or\
                            (p['stop_in_split'] >= gene_tuple[0] - 100 and p['stop_in_split'] <= gene_tuple[1] + 100):
                    level_ok = False
                    break
            if level_ok:
                levels_occupied[level].append((p['start_in_split'], p['stop_in_split']), )
                p['level'] = level
                break
        if not level_ok:
            levels_occupied[level + 1] = [(p['start_in_split'], p['stop_in_split']), ]
            p['level'] = level + 1

        data['genes'].append(p)

    progress.end()

    return json.dumps(data)


def store_collections_dict(args, d, request, response):
    if args.read_only:
        return json.dumps("Sorry! This is a read-only instance.")

    source = request.forms.get('source')
    data = json.loads(request.forms.get('data'))
    colors = json.loads(request.forms.get('colors'))

    if not len(source):
        run.info_single('Lousy attempt from the user to store their collection under an empty source identifier name :/')
        return json.dumps("Error: Collection name cannot be empty.")

    num_splits = sum(len(l) for l in data.values())
    if not num_splits:
        run.info_single('The user to store 0 splits as a collection :/')
        return json.dumps("Error: There are no selections to store (you haven't selected anything).")

    if source in d.collections.collections_dict:
        e = d.collections.collections_dict[source]
        if e['read_only']:
            run.info_single('Lousy attempt from the user to store their collection under "%s" :/' % source)
            return json.dumps("Well, '%s' is a read-only collection, so you need to come up with a different name... Sorry!" % source)

    run.info_single('A request to store %d bins that describe %d splits under the collection id "%s"\
                     has been made.' % (len(data), num_splits, source), cut_after=None)

    bins_info_dict = {}
    for bin_name in data:
        bins_info_dict[bin_name] = {'html_color': colors[bin_name], 'source': "anvi-interactive"}

    # the db here is either a profile db, or a pan db, but it can't be both:
    db_path = d.pan_db_path or d.profile_db_path
    collections = dbops.TablesForCollections(db_path)
    collections.append(source, data, bins_info_dict)

    # a new collection is stored in the database, but the interactive object
    # does not know about that and needs updatin'
    d.collections.populate_collections_dict(db_path)

    msg = "New collection '%s' with %d bin%s been stored." % (source, len(data), 's have' if len(data) > 1 else ' has')
    run.info_single(msg)
    return json.dumps(msg)


def store_refined_bins(args, r, request, response):
    data = json.loads(request.forms.get('data'))
    colors = json.loads(request.forms.get('colors'))

    bins_info_dict = {}
    for bin_name in data:
        bins_info_dict[bin_name] = {'html_color': colors[bin_name], 'source': "anvi-refine"}

    try:
        r.store_refined_bins(data, bins_info_dict)
    except RefineError as e:
        return json.dumps({'status': -1, 'message': e.clear_text()})

    message = 'Done! Collection %s is updated in the database. You can close your browser window (or continue updating).' % (r.collection_name)
    return json.dumps({'status': 0, 'message': message})


def gen_summary(args, d, request, response, collection_name):
    set_default_headers(response)

    if args.read_only:
        return json.dumps({'error': "Sorry! This is a read-only instance."})

    if d.manual_mode:
        return json.dumps({'error': "Creating summaries is only possible with proper anvi'o runs at the moment :/"})

    run.info_single('A summary of collection "%s" has been requested.' % collection_name)

    # get a dummy args instance, and fill it down below
    summarizer_args = summarizer.ArgsTemplateForSummarizerClass()

    # common params. we will set pan/profile specific params a bit later:
    summarizer_args.collection_name = collection_name
    summarizer_args.taxonomic_level = d.taxonomic_level

    if d.mode == 'pan':
        summarizer_args.pan_db = d.pan_db_path
        summarizer_args.genomes_storage = d.genomes_storage_path
        summarizer_args.output_dir = os.path.join(os.path.dirname(summarizer_args.pan_db), 'SUMMARY_%s' % collection_name)
    elif d.mode == 'full':
        summarizer_args.profile_db = d.profile_db_path
        summarizer_args.contigs_db = d.contigs_db_path
        summarizer_args.output_dir = os.path.join(os.path.dirname(summarizer_args.profile_db), 'SUMMARY_%s' % collection_name)
    else:
        return json.dumps({'error': 'We do not know anything about this mode: "%s"' % d.mode})

    # call the summary:
    try:
        summary = summarizer.PanSummarizer(summarizer_args, r=run, p=progress) if d.mode == 'pan' else summarizer.ProfileSummarizer(summarizer_args, r=run, p=progress)
        summary.process()
    except Exception as e:
        return json.dumps({'error': 'Something failed in the "%s" summary mode. This is what we know: %s' % (d.mode, e)})

    run.info_single('HTML output for summary is ready: %s' % summary.index_html)

    path = "summary/%s/index.html" % (collection_name)
    return json.dumps({'path': path})


def send_summary_static(args, d, request, response, collection_name, filename):
    set_default_headers(response)

    if d.mode == 'pan':
        return static_file(filename, root=os.path.join(os.path.dirname(d.pan_db_path), 'SUMMARY_%s' % collection_name))
    elif d.mode == 'full':
        return static_file(filename, root=os.path.join(os.path.dirname(d.profile_db_path), 'SUMMARY_%s' % collection_name))
    else:
        return json.dumps({'error': 'Something failed. This is what we know: %s' % e})


def get_collection_dict(args, d, request, response, collection_name):
    run.info_single('Data for collection "%s" has been requested.' % collection_name)
    set_default_headers(response)

    collection_dict = d.collections.get_collection_dict(collection_name)
    bins_info_dict = d.collections.get_bins_info_dict(collection_name)

    colors_dict = {}
    for bin_name in bins_info_dict:
        colors_dict[bin_name] = bins_info_dict[bin_name]['html_color']

    return json.dumps({'data': collection_dict,
                       'colors': colors_dict})


def get_collections(args, d, request, response):
    csd = d.collections.collections_dict
    run.info_single('Collection sources has been requested (info dict with %d item(s) has been returned).' % len(csd), cut_after=None)
    set_default_headers(response)
    return json.dumps(csd)


def get_tree(args, d, request, response, tree_id):
    set_default_headers(response)

    if tree_id in d.p_meta['clusterings']:
        run.info_single("Clustering of '%s' has been requested" % (tree_id))
        return json.dumps(d.p_meta['clusterings'][tree_id]['newick'])

    return json.dumps("")


def get_sequence_for_gene_call(args, d, request, response, gene_callers_id):
    set_default_headers(response)

    try:
        gene_callers_id = int(gene_callers_id)
    except:
        return json.dumps({'error': "Gene caller id does not seem to be 'integerable'. Not good :/"})

    try:
        gene_calls_tuple = d.get_sequences_for_gene_callers_ids([gene_callers_id])
    except Exception as e:
        return json.dumps({'error': "Something went wrong when I tried to access to that gene: '%s' :/" % e})

    entry = gene_calls_tuple[1][gene_callers_id]
    sequence = entry['sequence']
    header = '%d|' % (gene_callers_id) + '|'.join(['%s:%s' % (k, str(entry[k])) for k in ['contig', 'start', 'stop', 'direction', 'rev_compd', 'length']])

    return json.dumps({'sequence': sequence, 'header': header})


def get_sequence_for_split(args, d, request, response, split_name):
    set_default_headers(response)

    try:
        sequence = d.split_sequences[split_name]
        header = split_name
    except Exception as e:
        return json.dumps({'error': "Something went wrong when I tried to access that split sequence: '%s' :/" % e})

    return json.dumps({'sequence': sequence, 'header': header})


def get_hmm_hit_from_bin(args, d, request, response, bin_name, gene_name):
    set_default_headers(response)

    if d.mode != 'collection':
        return json.dumps({'error': "HMM hits from bins can only be requested in 'collection' mode. You are doing something wrong..."})

    if not d.collection:
        return json.dumps({'error': "You are in 'collection' mode, but your collection is empty. You are killing me."})

    hmm_sequences_dict = d.hmm_access.get_sequences_dict_for_hmm_hits_in_splits({bin_name: set(d.collection[bin_name])})
    gene_sequences = utils.get_filtered_dict(hmm_sequences_dict, 'gene_name', set([gene_name]))

    if not gene_sequences:
        return json.dumps({'error': "Sorry. It seems %s does not have a hit for %s." % (bin_name, gene_name)})

    unique_id_for_longest_hit = sorted([(gene_sequences[gene_id]['length'], gene_id) for gene_id in gene_sequences], reverse=True)[0][1]

    header, sequence = d.hmm_access.get_FASTA_header_and_sequence_for_gene_unique_id(gene_sequences, unique_id_for_longest_hit)

    return json.dumps({'sequence': sequence, 'header': header})


def get_view_data(args, d, request, response, view_id):
    set_default_headers(response)

    return json.dumps(d.views[view_id])
