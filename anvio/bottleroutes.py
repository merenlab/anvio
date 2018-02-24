# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Common routes for bottle web server.

    Functions defined here are used from wihtin programs such as
    anvi-interactive, or anvi-refine.
"""

import os
import re
import io
import sys
import copy
import time
import json
import random
import argparse
import requests
import datetime
from multiprocessing import Process

from bottle import Bottle
from bottle import BaseRequest
from bottle import redirect, static_file

import anvio
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.drivers as drivers
import anvio.terminal as terminal
import anvio.summarizer as summarizer
import anvio.filesnpaths as filesnpaths
import anvio.auxiliarydataops as auxiliarydataops

from anvio.serverAPI import AnviServerAPI
from anvio.errors import RefineError, ConfigError
from anvio.tables.miscdata import TableForLayerOrders
from anvio.tables.collections import TablesForCollections


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = ["A. Murat Eren"]
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Ozcan Esen"
__email__ = "ozcanesen@gmail.com"


run = terminal.Run()
progress = terminal.Progress()

# increase maximum size of form data to 100 MB
BaseRequest.MEMFILE_MAX = 1024 * 1024 * 100


class BottleApplication(Bottle):
    def __init__(self, interactive, args, mock_request=None, mock_response=None):
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        super(BottleApplication, self).__init__()
        self.interactive = interactive
        self.args = args
        self.read_only = A('read_only')
        self.browser_path = A('browser_path')
        self.export_svg = A('export_svg')
        self.server_only = A('server_only')

        self.unique_session_id = random.randint(0,9999999999)
        self.static_dir = os.path.join(os.path.dirname(utils.__file__), 'data/interactive')

        self.register_hooks()
        self.register_routes()

        # this part is required to inject request and responses from anvi server
        if mock_request or mock_response:
            global request
            global response
            request = mock_request
            response = mock_response
        else:
            from bottle import response, request


    def register_hooks(self):
        self.add_hook('before_request', self.set_default_headers)


    def set_default_headers(self):
        response.set_header('Content-Type', 'application/json')
        response.set_header('Pragma', 'no-cache')
        response.set_header('Cache-Control', 'no-cache, no-store, max-age=0, must-revalidate')
        response.set_header('Expires', 'Thu, 01 Dec 1994 16:00:00 GMT')


    def register_routes(self):
        self.route('/',                                        callback=self.redirect_to_app)
        self.route('/app/:filename#.*#',                       callback=self.send_static)
        self.route('/data/news',                               callback=self.get_news)
        self.route('/data/<name>',                             callback=self.send_data)
        self.route('/data/view/<view_id>',                     callback=self.get_view_data)
        self.route('/tree/<items_order_id>',                   callback=self.get_items_order)
        self.route('/state/autoload',                          callback=self.state_autoload)
        self.route('/state/all',                               callback=self.state_all)
        self.route('/state/get/<state_name>',                  callback=self.get_state)
        self.route('/state/save/<state_name>',                 callback=self.save_state, method='POST')
        self.route('/data/charts/<order_name>/<item_name>',    callback=self.charts, method='POST')
        self.route('/data/completeness',                       callback=self.completeness, method='POST')
        self.route('/data/collections',                        callback=self.get_collections)
        self.route('/data/collection/<collection_name>',       callback=self.get_collection_dict)
        self.route('/store_collection',                        callback=self.store_collections_dict, method='POST')
        self.route('/store_description',                       callback=self.store_description, method='POST')
        self.route('/upload_project',                          callback=self.upload_project, method='POST')
        self.route('/data/contig/<split_name>',                callback=self.get_sequence_for_split)
        self.route('/summarize/<collection_name>',             callback=self.gen_summary)
        self.route('/summary/<collection_name>/:filename#.*#', callback=self.send_summary_static)
        self.route('/data/gene/<gene_callers_id>',             callback=self.get_sequence_for_gene_call)
        self.route('/data/hmm/<bin_name>/<gene_name>',         callback=self.get_hmm_hit_from_bin)
        self.route('/data/geneclusterssummary',             callback=self.get_gene_clusters_summary, method='POST')
        self.route('/data/get_AA_sequences_for_gene_cluster/<gene_cluster_name>',  callback=self.get_AA_sequences_for_gene_cluster)
        self.route('/data/pan_gene_popup/<gene_callers_id>/<genome_name>',         callback=self.get_gene_popup_for_pan)
        self.route('/data/geneclusters/<order_name>/<gene_cluster_name>',          callback=self.inspect_gene_cluster)
        self.route('/data/charts_for_single_gene/<order_name>/<item_name>',        callback=self.charts_for_single_gene, method='POST')
        self.route('/data/store_refined_bins',                 callback=self.store_refined_bins, method='POST')
        self.route('/data/phylogeny/aligners',                 callback=self.get_available_aligners)
        self.route('/data/phylogeny/programs',                 callback=self.get_available_phylogeny_programs)
        self.route('/data/phylogeny/generate_tree',            callback=self.generate_tree, method='POST')
        self.route('/data/search_functions',                   callback=self.search_functions, method='POST')
        self.route('/data/get_contigs_stats',                  callback=self.get_contigs_stats)
        self.route('/data/filter_gene_clusters',               callback=self.filter_gene_clusters, method='POST')


    def run_application(self, ip, port):
        try:
            server_process = Process(target=self.run, kwargs={'host': ip, 'port': port, 'quiet': True, 'server': 'cherrypy'})
            server_process.start()

            if self.export_svg:
                try:
                    utils.run_selenium_and_export_svg("http://%s:%d/app/index.html" % (ip, port),
                                                      self.args.export_svg,
                                                      browser_path=self.browser_path,
                                                      run=run)
                except Exception as e:
                    print(e)
                finally:
                    server_process.terminate()
                    sys.exit(0)

            if not self.server_only:
                # Sometimes browser opens before web server actually starts so we see
                # message like "Website can not be reached" and user needs to refresh
                # I have added sleep below to delay web browser little bit.
                time.sleep(1.5)

                utils.open_url_in_browser(url="http://%s:%d" % (ip, port),
                                          browser_path=self.browser_path,
                                          run=run)

            run.info_single('The server is now listening the port number "%d". When you are finished, press CTRL+C to terminate the server.' % port, 'green', nl_before = 1, nl_after=1)
            server_process.join()
        except KeyboardInterrupt:
            run.warning('The server is being terminated.', header='Please wait...')
            server_process.terminate()
            sys.exit(0)


    def redirect_to_app(self):
        homepage = 'index.html' 
        if self.interactive.mode == 'contigs':
            homepage = 'contigs.html'

        redirect('/app/%s?rand=%s' % (homepage, self.random_hash(8)))


    def send_static(self, filename):
        ret = static_file(filename, root=self.static_dir)
        ret.set_header('Pragma', 'no-cache')
        ret.set_header('Cache-Control', 'no-cache, no-store, max-age=0, must-revalidate')
        ret.set_header('Expires', 'Thu, 01 Dec 1994 16:00:00 GMT')

        # cache killer, it adds random query string to .js, .css source urls.
        if filename.endswith('.html'):
            pattern = re.compile(b".*(<script|<link).*(href|src)\=[\'\"]((?!http\:\/\/).+?)\".*", re.MULTILINE)

            buff = io.BytesIO()
            index = 0
            for result in re.finditer(pattern, ret.body.read()):
                pos = result.end(3)
                suffix = b'?rand=' + self.random_hash(32).encode() 

                # read chunk from original file and write to buffer, 
                # then store pos to index, next iteration we are going 
                # to read from that position
                ret.body.seek(index)
                buff.write(ret.body.read(pos - index))
                buff.write(suffix)
                index = pos

            # write rest of the file
            ret.body.seek(index)
            buff.write(ret.body.read())
            ret.body = buff
            ret.body.seek(0)
            ret.headers['Content-Length'] = buff.getbuffer().nbytes

        return ret


    def get_news(self):
        ret = []
        try:
            news_markdown = requests.get('https://raw.githubusercontent.com/merenlab/anvio/master/NEWS.md')
            news_items = news_markdown.text.split("***")
            
            """ FORMAT
            # Title with spaces (01.01.1970) #
            Lorem ipsum, dolor sit amet
            ***
            # Title with spaces (01.01.1970) #
            Lorem ipsum, dolor sit amet
            ***
            # Title with spaces (01.01.1970) #
            Lorem ipsum, dolor sit amet
            """
            for news_item in news_items:
                if len(news_item) < 5:
                    # too short to parse, just skip it
                    continue

                ret.append({
                        'date': news_item.split("(")[1].split(")")[0].strip(),
                        'title': news_item.split("#")[1].split("(")[0].strip(),
                        'content': news_item.split("#\n")[1].strip()
                    })
        except:
            ret.append({
                    'date': '',
                    'title': 'Something has failed',
                    'content': 'Anvi\'o failed to retrieve any news for you, maybe you do not have internet connection or something :('
                })

        return json.dumps(ret)


    def random_hash(self, size=8):
        r = random.getrandbits(size * 4)
        return '{1:0{0}x}'.format(size, r)


    def send_data(self, name):
        if name == "init":
            bin_prefix = "Bin_"
            if self.interactive.mode == 'refine':
                bin_prefix = list(self.interactive.bins)[0] + "_" if len(self.interactive.bins) == 1 else "Refined_",

            return json.dumps( { "title":                              self.interactive.title,
                                 "description":                        (self.interactive.p_meta['description']),
                                 "item_orders":                        (self.interactive.p_meta['default_item_order'], self.interactive.p_meta['item_orders']),
                                 "views":                              (self.interactive.default_view, dict(list(zip(list(self.interactive.views.keys()), list(self.interactive.views.keys()))))),
                                 "contigLengths":                      dict([tuple((c, self.interactive.splits_basic_info[c]['length']),) for c in self.interactive.splits_basic_info]),
                                 "defaultView":                        self.interactive.views[self.interactive.default_view],
                                 "mode":                               self.interactive.mode,
                                 "server_mode":                        False,
                                 "readOnly":                           self.read_only,
                                 "binPrefix":                          bin_prefix,
                                 "sessionId":                          self.unique_session_id,
                                 "samplesOrder":                       self.interactive.layers_order_data_dict,
                                 "sampleInformation":                  self.interactive.layers_additional_data_dict,
                                 "sampleInformationDefaultLayerOrder": self.interactive.layers_additional_data_keys,
                                 "stateAutoload":                      self.interactive.state_autoload,
                                 "collectionAutoload":                 self.interactive.collection_autoload,
                                 "noPing":                             False,
                                 "inspectionAvailable":                self.interactive.auxiliary_profile_data_available,
                                 "sequencesAvailable":                 True if (self.interactive.split_sequences or self.interactive.mode == 'gene') else False,
                                 "functions_initialized":              self.interactive.gene_function_calls_initiated })

        elif name == "session_id":
            return json.dumps(self.unique_session_id)


    def get_view_data(self, view_id):
        return json.dumps(self.interactive.views[view_id])


    def get_items_order(self, items_order_id):
        if items_order_id in self.interactive.p_meta['item_orders']:
            items_order = self.interactive.p_meta['item_orders'][items_order_id]

            if items_order['type'] == 'newick':
                run.info_single("The newick order '%s' has been requested" % (items_order_id))
            elif items_order['type'] == 'basic':
                run.info_single("The basic order '%s' has been requested" % (items_order_id))
            else:
                return json.dumps({'error': "The interface requested something anvi'o doesn't know about. Item orders\
                                             can only be in the form of 'newick' or 'basic'. But the interface requested\
                                             a '%s'. We are all confused here :/" % items_order_id})

            return json.dumps(items_order['data'])

        return json.dumps("")


    def state_autoload(self):
        return json.dumps(self.interactive.state_autoload)


    def state_all(self):
        return json.dumps(self.interactive.states_table.states)


    def save_state(self, state_name):
        if self.read_only:
            return json.dumps({'status_code': '0'})

        content = request.forms.get('content')
        last_modified = datetime.datetime.now().strftime("%d.%m.%Y %H:%M:%S")

        self.interactive.states_table.store_state(state_name, content, last_modified)

        return json.dumps({'status_code': '1'})


    def get_state(self, state_name):
        if state_name in self.interactive.states_table.states:
            state = self.interactive.states_table.states[state_name]
            return json.dumps(state['content'])

        return json.dumps("")


    def charts(self, order_name, item_name):
        title = None
        if self.interactive.mode == 'gene':
            split_name = self.interactive.gene_callers_id_to_split_name_dict[int(item_name)]
            title = "Gene '%d' in split '%s'" % (int(item_name), split_name)
        else:
            split_name = item_name
            title = split_name

        state = json.loads(request.forms.get('state'))

        data = {'layers': [],
                 'title': title,
                 'index': None,
                 'total': None,
                 'coverage': [],
                 'variability': [],
                 'competing_nucleotides': [],
                 'previous_contig_name': None,
                 'next_contig_name': None,
                 'genes': [],
                 'outlier_SNVs_shown': not self.args.hide_outlier_SNVs}

        if split_name not in self.interactive.split_names:
            return data

        if not self.interactive.auxiliary_profile_data_available:
            return data

        data['index'], data['total'], data['previous_contig_name'], data['next_contig_name'] = self.get_index_total_previous_and_next_items(order_name, item_name)

        layers = [layer for layer in sorted(self.interactive.p_meta['samples']) if float(state['layers'][layer]['height']) > 0]

        auxiliary_coverages_db = auxiliarydataops.AuxiliaryDataForSplitCoverages(self.interactive.auxiliary_data_path,
                                                                                 self.interactive.p_meta['contigs_db_hash'])
        coverages = auxiliary_coverages_db.get(split_name)
        auxiliary_coverages_db.close()

        data['coverage'] = [coverages[layer].tolist() for layer in layers]

        ## get the variability information dict for split:
        progress.new('Variability', discard_previous_if_exists=True)
        progress.update('Collecting info for "%s"' % split_name)
        split_variability_info_dict = self.interactive.get_variability_information_for_split(split_name, skip_outlier_SNVs=self.args.hide_outlier_SNVs)

        for layer in layers:
            progress.update('Formatting variability data: "%s"' % layer)
            data['layers'].append(layer)
            data['competing_nucleotides'].append(split_variability_info_dict[layer]['competing_nucleotides'])
            data['variability'].append(split_variability_info_dict[layer]['variability'])

        levels_occupied = {1: []}
        for entry_id in  self.interactive.split_name_to_genes_in_splits_entry_ids[split_name]:
            gene_callers_id =  self.interactive.genes_in_splits[entry_id]['gene_callers_id']
            p =  self.interactive.genes_in_splits[entry_id]
            # p looks like this at this point:
            #
            # {'percentage_in_split': 100,
            #  'start_in_split'     : 16049,
            #  'stop_in_split'      : 16633}
            #  'prot'               : u'prot2_03215',
            #  'split'              : u'D23-1contig18_split_00036'}
            #
            # we will add a bit more attributes:
            p['source'] =  self.interactive.genes_in_contigs_dict[gene_callers_id]['source']
            p['direction'] =  self.interactive.genes_in_contigs_dict[gene_callers_id]['direction']
            p['start_in_contig'] =  self.interactive.genes_in_contigs_dict[gene_callers_id]['start']
            p['stop_in_contig'] =  self.interactive.genes_in_contigs_dict[gene_callers_id]['stop']
            p['complete_gene_call'] = 'No' if  self.interactive.genes_in_contigs_dict[gene_callers_id]['partial'] else 'Yes'
            p['length'] = p['stop_in_contig'] - p['start_in_contig']
            p['functions'] =  self.interactive.gene_function_calls_dict[gene_callers_id] if gene_callers_id in  self.interactive.gene_function_calls_dict else None

            for level in levels_occupied:
                level_ok = True
                for gene_tuple in levels_occupied[level]:
                    if (p['start_in_split'] >= gene_tuple[0] - 100 and p['start_in_split'] <= gene_tuple[1] + 100) or\
                                (p['stop_in_split'] >= gene_tuple[0] - 100 and p['stop_in_split'] <= gene_tuple[1] + 100) or \
                                (p['start_in_split'] <= gene_tuple[0] - 100 and p['stop_in_split'] >= gene_tuple[1] + 100):
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


    def charts_for_single_gene(self, order_name, item_name):
        gene_callers_id = int(item_name)
        split_name = self.interactive.gene_callers_id_to_split_name_dict[gene_callers_id]
        gene_info = [e for e in self.interactive.genes_in_splits.values() if e['gene_callers_id'] == gene_callers_id][0]

        focus_region_start, focus_region_end = max(0, gene_info['start_in_split'] - 100), min(self.interactive.split_lengths_info[split_name], gene_info['stop_in_split'] + 100)

        state = json.loads(request.forms.get('state'))
        data = {'layers': [],
                 'title': "Gene '%d' in split '%s'" % (gene_callers_id, split_name),
                 'index': None,
                 'total': None,
                 'coverage': [],
                 'variability': [],
                 'competing_nucleotides': [],
                 'previous_contig_name': None,
                 'next_contig_name': None,
                 'genes': [],
                 'outlier_SNVs_shown': not self.args.hide_outlier_SNVs}

        data['index'], data['total'], data['previous_contig_name'], data['next_contig_name'] = self.get_index_total_previous_and_next_items(order_name, str(gene_callers_id))

        layers = [layer for layer in sorted(self.interactive.p_meta['samples']) if float(state['layers'][layer]['height']) > 0]

        auxiliary_coverages_db = auxiliarydataops.AuxiliaryDataForSplitCoverages(self.interactive.auxiliary_data_path,
                                                                                 self.interactive.p_meta['contigs_db_hash'])
        coverages = auxiliary_coverages_db.get(split_name)
        auxiliary_coverages_db.close()

        data['coverage'] = []
        for layer in layers:
            coverage_list = coverages[layer].tolist()
            # gene -+ 100 gap if possible
            data['coverage'].append(coverage_list[focus_region_start:focus_region_end])

        ## get the variability information dict for split:
        progress.new('Variability')
        progress.update('Collecting info for "%s"' % split_name)
        split_variability_info_dict = self.interactive.get_variability_information_for_split(split_name, skip_outlier_SNVs=self.args.hide_outlier_SNVs)

        for layer in layers:
            progress.update('Formatting variability data: "%s"' % layer)
            data['layers'].append(layer)
            data['competing_nucleotides'].append(split_variability_info_dict[layer]['competing_nucleotides'])
            data['variability'].append(split_variability_info_dict[layer]['variability'])

        levels_occupied = {1: []}
        for entry_id in self.interactive.split_name_to_genes_in_splits_entry_ids[split_name]:
            gene_callers_id = self.interactive.genes_in_splits[entry_id]['gene_callers_id']
            p =  self.interactive.genes_in_splits[entry_id]
            # p looks like this at this point:
            #
            # {'percentage_in_split': 100,
            #  'start_in_split'     : 16049,
            #  'stop_in_split'      : 16633}
            #  'prot'               : u'prot2_03215',
            #  'split'              : u'D23-1contig18_split_00036'}
            #

            if p['start_in_split'] <= focus_region_start and p['stop_in_split'] <= focus_region_start:
                continue
            if p['start_in_split'] >= focus_region_end and p['stop_in_split'] >= focus_region_end:
                continue

            # because Python. when we don't do this, the organization of genes in the interface split pages
            # gets all screwed up in gene view due the permanence of the changes in the dictionary.
            p = copy.deepcopy(p)

            # add offset
            p['start_in_split'] -= focus_region_start
            p['stop_in_split'] -= focus_region_start

            # we will add a bit more attributes:
            p['source'] =  self.interactive.genes_in_contigs_dict[gene_callers_id]['source']
            p['direction'] =  self.interactive.genes_in_contigs_dict[gene_callers_id]['direction']
            p['start_in_contig'] =  self.interactive.genes_in_contigs_dict[gene_callers_id]['start']
            p['stop_in_contig'] =  self.interactive.genes_in_contigs_dict[gene_callers_id]['stop']
            p['complete_gene_call'] = 'No' if  self.interactive.genes_in_contigs_dict[gene_callers_id]['partial'] else 'Yes'
            p['length'] = p['stop_in_contig'] - p['start_in_contig']
            p['functions'] =  self.interactive.gene_function_calls_dict[gene_callers_id] if gene_callers_id in  self.interactive.gene_function_calls_dict else None

            for level in levels_occupied:
                level_ok = True
                for gene_tuple in levels_occupied[level]:
                    if (p['start_in_split'] >= gene_tuple[0] - 100 and p['start_in_split'] <= gene_tuple[1] + 100) or\
                                (p['stop_in_split'] >= gene_tuple[0] - 100 and p['stop_in_split'] <= gene_tuple[1] + 100) or \
                                (p['start_in_split'] <= gene_tuple[0] - 100 and p['stop_in_split'] >= gene_tuple[1] + 100):
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


    def get_index_total_previous_and_next_items(self, order_name, item_name):
        previous_item_name = None
        next_item_name = None
        index = None
        total = None

        # FIXME: improve performance here
        items_order_entry = self.interactive.p_meta['item_orders'][order_name]
        items_order = None
        if items_order_entry['type'] == 'newick':
            items_order = utils.get_names_order_from_newick_tree(items_order_entry['data'])
        else:
            items_order = items_order_entry['data']

        index_of_item = items_order.index(item_name)
        if index_of_item:
            previous_item_name = items_order[index_of_item - 1]
        if (index_of_item + 1) < len(items_order):
            next_item_name = items_order[index_of_item + 1]

        index = index_of_item + 1
        total = len(items_order)

        return index, total, previous_item_name, next_item_name


    def completeness(self):
        completeness_sources = {}
        completeness_averages = {}
        if not self.interactive.completeness:
            return json.dumps(completeness_sources)

        split_names = json.loads(request.forms.get('split_names'))
        bin_name = json.loads(request.forms.get('bin_name'))

        run.info_single('Completeness info has been requested for %d splits in %s' % (len(split_names), bin_name))

        p_completion, p_redundancy, scg_domain, domain_confidence, results_dict = self.interactive.completeness.get_info_for_splits(set(split_names))

        # convert results_dict (where domains are the highest order items) into a dict that is compatible with the
        # previous format of the dict (where hmm scg source names are the higher order items).
        for domain in results_dict:
            for source in results_dict[domain]:
                completeness_sources[source] = results_dict[domain][source]

        completeness_averages['percent_completion'] = p_completion
        completeness_averages['percent_redundancy'] = p_redundancy
        completeness_averages['domain'] = scg_domain
        completeness_averages['domain_confidence'] = domain_confidence

        return json.dumps({'stats': completeness_sources, 'averages': completeness_averages, 'refs': self.interactive.completeness.http_refs})


    def get_collections(self):
        csd = self.interactive.collections.collections_dict
        run.info_single('Collection sources has been requested (info dict with %d item(s) has been returned).' % len(csd), cut_after=None)
        return json.dumps(csd)


    def get_collection_dict(self, collection_name):
        run.info_single('Data for collection "%s" has been requested.' % collection_name)
        #set_default_headers(response)

        collection_dict = self.interactive.collections.get_collection_dict(collection_name)
        bins_info_dict = self.interactive.collections.get_bins_info_dict(collection_name)

        colors_dict = {}
        for bin_name in bins_info_dict:
            colors_dict[bin_name] = bins_info_dict[bin_name]['html_color']

        return json.dumps({'data': collection_dict, 'colors': colors_dict})


    def store_collections_dict(self):
        if self.read_only:
            return json.dumps("Sorry! This is a read-only instance.")

        source = request.forms.get('source')
        data = json.loads(request.forms.get('data'))
        colors = json.loads(request.forms.get('colors'))

        if not len(source):
            run.info_single('Lousy attempt from the user to store their collection under an empty source identifier name :/')
            return json.dumps("Error: Collection name cannot be empty.")

        num_splits = sum(len(l) for l in list(data.values()))
        if not num_splits:
            run.info_single('The user to store 0 splits as a collection :/')
            return json.dumps("Error: There are no selections to store (you haven't selected anything).")

        if source in self.interactive.collections.collections_dict:
            e = self.interactive.collections.collections_dict[source]
            if e['read_only']:
                run.info_single('Lousy attempt from the user to store their collection under "%s" :/' % source)
                return json.dumps("Well, '%s' is a read-only collection, so you need to come up with a different name... Sorry!" % source)

        run.info_single('A request to store %d bins that describe %d splits under the collection id "%s"\
                         has been made.' % (len(data), num_splits, source), cut_after=None)

        bins_info_dict = {}
        for bin_name in data:
            bins_info_dict[bin_name] = {'html_color': colors[bin_name], 'source': "anvi-interactive"}

        # the db here is either a profile db, or a pan db, but it can't be both:
        db_path = self.interactive.pan_db_path or self.interactive.profile_db_path
        collections = TablesForCollections(db_path)
        try:
            collections.append(source, data, bins_info_dict)
        except ConfigError as e:
            return json.dumps(e.clear_text())

        # a new collection is stored in the database, but the interactive object
        # does not know about that and needs updatin'
        self.interactive.collections.populate_collections_dict(db_path)

        msg = "New collection '%s' with %d bin%s been stored." % (source, len(data), 's have' if len(data) > 1 else ' has')
        run.info_single(msg)
        return json.dumps(msg)


    def store_description(self):
        if self.read_only:
            return

        description = request.forms.get('description')

        db_path = self.interactive.pan_db_path or self.interactive.profile_db_path
        dbops.update_description_in_db(db_path, description)
        self.interactive.p_meta['description'] = description
        return json.dumps("")


    def get_sequence_for_split(self, split_name):
        try:
            sequence = self.interactive.split_sequences[split_name]
            header = split_name
        except Exception as e:
            return json.dumps({'error': "Something went wrong when I tried to access that split sequence: '%s' :/" % e})

        return json.dumps({'sequence': sequence, 'header': header})


    def gen_summary(self, collection_name):
        #set_default_headers(response)

        if self.read_only:
            return json.dumps({'error': "Sorry! This is a read-only instance."})

        if self.interactive.mode == 'manual':
            return json.dumps({'error': "Creating summaries is only possible with proper anvi'o runs at the moment :/"})

        run.info_single('A summary of collection "%s" has been requested.' % collection_name)

        # get a dummy args instance, and fill it down below
        summarizer_args = summarizer.ArgsTemplateForSummarizerClass()

        # common params. we will set pan/profile specific params a bit later:
        summarizer_args.collection_name = collection_name
        summarizer_args.taxonomic_level = self.interactive.taxonomic_level

        if self.interactive.mode == 'pan':
            summarizer_args.pan_db = self.interactive.pan_db_path
            summarizer_args.genomes_storage = self.interactive.genomes_storage_path
            summarizer_args.output_dir = os.path.join(os.path.dirname(summarizer_args.pan_db), 'SUMMARY_%s' % collection_name)
        elif self.interactive.mode == 'full':
            summarizer_args.profile_db = self.interactive.profile_db_path
            summarizer_args.contigs_db = self.interactive.contigs_db_path
            summarizer_args.output_dir = os.path.join(os.path.dirname(summarizer_args.profile_db), 'SUMMARY_%s' % collection_name)
        else:
            return json.dumps({'error': 'We do not know anything about this mode: "%s"' % self.interactive.mode})

        # call the summary:
        try:
            summary = summarizer.PanSummarizer(summarizer_args, r=run, p=progress) if self.interactive.mode == 'pan' else summarizer.ProfileSummarizer(summarizer_args, r=run, p=progress)
            summary.process()
        except Exception as e:
            return json.dumps({'error': 'Something failed in the "%s" summary mode. This is what we know: %s' % (self.interactive.mode, e)})

        run.info_single('HTML output for summary is ready: %s' % summary.index_html)

        path = "summary/%s/index.html" % (collection_name)
        return json.dumps({'path': path})


    def send_summary_static(self, collection_name, filename):
        if self.interactive.mode == 'pan':
            ret = static_file(filename, root=os.path.join(os.path.dirname(self.interactive.pan_db_path), 'SUMMARY_%s' % collection_name))
            ret.set_header('Pragma', 'no-cache')
            ret.set_header('Cache-Control', 'no-cache, no-store, max-age=0, must-revalidate')
            ret.set_header('Expires', 'Thu, 01 Dec 1994 16:00:00 GMT')
            return ret
        elif self.interactive.mode == 'full':
            ret = static_file(filename, root=os.path.join(os.path.dirname(self.interactive.profile_db_path), 'SUMMARY_%s' % collection_name))
            ret.set_header('Pragma', 'no-cache')
            ret.set_header('Cache-Control', 'no-cache, no-store, max-age=0, must-revalidate')
            ret.set_header('Expires', 'Thu, 01 Dec 1994 16:00:00 GMT')
            return ret
        else:
            return json.dumps({'error': 'The server has no idea how to handle the mode "%s" :/' % self.interactive.mode})


    def get_sequence_for_gene_call(self, gene_callers_id):
        try:
            gene_callers_id = int(gene_callers_id)
        except:
            return json.dumps({'error': "Gene caller id does not seem to be 'integerable'. Not good :/"})

        try:
            gene_calls_tuple = self.interactive.get_sequences_for_gene_callers_ids([gene_callers_id])
        except Exception as e:
            return json.dumps({'error': "Something went wrong when I tried to access to that gene: '%s' :/" % e})

        entry = gene_calls_tuple[1][gene_callers_id]
        sequence = entry['sequence']
        header = '%d|' % (gene_callers_id) + '|'.join(['%s:%s' % (k, str(entry[k])) for k in ['contig', 'start', 'stop', 'direction', 'rev_compd', 'length']])

        return json.dumps({'sequence': sequence, 'header': header})


    def get_gene_popup_for_pan(self, gene_callers_id, genome_name):
        if not self.interactive.genomes_storage_is_available:
            return json.dumps({'error': 'Genome storage does not seem to be available :/ So that button will not work..'})

        gene_callers_id = int(gene_callers_id)

        if genome_name not in self.interactive.genomes_storage.gene_info:
            return json.dumps({'error': "Your request contains a genome name anvi'o genomes storage does not know about. What are you doing?"})

        if gene_callers_id not in self.interactive.genomes_storage.gene_info[genome_name]:
            return json.dumps({'error': "Your gene caller id does not work for anvi'o :("})

        return json.dumps({'status': 0, 'gene_info': self.interactive.genomes_storage.gene_info[genome_name][gene_callers_id]})


    def get_hmm_hit_from_bin(self, bin_name, gene_name):
        if self.interactive.mode != 'collection':
            return json.dumps({'error': "HMM hits from bins can only be requested in 'collection' mode. You are doing something wrong..."})

        if not self.interactive.collection:
            return json.dumps({'error': "You are in 'collection' mode, but your collection is empty. You are killing me."})

        hmm_sequences_dict = self.interactive.hmm_access.get_sequences_dict_for_hmm_hits_in_splits({bin_name: set(self.interactive.collection[bin_name])})
        gene_sequences = utils.get_filtered_dict(hmm_sequences_dict, 'gene_name', set([gene_name]))

        if not gene_sequences:
            return json.dumps({'error': "Sorry. It seems %s does not have a hit for %s." % (bin_name, gene_name)})

        unique_id_for_longest_hit = sorted([(gene_sequences[gene_id]['length'], gene_id) for gene_id in gene_sequences], reverse=True)[0][1]

        header, sequence = self.interactive.hmm_access.get_FASTA_header_and_sequence_for_gene_unique_id(gene_sequences, unique_id_for_longest_hit)

        return json.dumps({'sequence': sequence, 'header': header})


    def get_gene_clusters_summary(self):
        gene_cluster_ids = json.loads(request.forms.get('split_names'))
        bin_name = json.loads(request.forms.get('bin_name'))

        summary = self.interactive.get_summary_for_gene_clusters_list(gene_cluster_ids)

        run.info_single('Gene cluster info has been requested for %d items in %s' % (len(gene_cluster_ids), bin_name))

        return json.dumps({'functions': summary['functions'],
                           'num_gene_clusters': summary['num_gene_clusters'],
                           'genomes_contributing': summary['genomes_contributing'],
                           'num_gene_calls': summary['num_gene_calls']})


    def get_AA_sequences_for_gene_cluster(self, gene_cluster_name):
        data = {}

        if gene_cluster_name not in self.interactive.gene_clusters:
            return data

        if not self.interactive.genomes_storage_is_available:
            return data

        # add the list of gene caller ids associated with this pootein cluster into `data`:
        for genome_name in self.interactive.gene_clusters[gene_cluster_name]:
            for gene_callers_id in self.interactive.gene_clusters[gene_cluster_name][genome_name]:
                data['%s_%s' % (genome_name, str(gene_callers_id))] = self.interactive.genomes_storage.get_gene_sequence(genome_name, gene_callers_id)

        return json.dumps(data)


    def inspect_gene_cluster(self, order_name, gene_cluster_name):
        data = {'gene_cluster_name': gene_cluster_name,
                'genomes': [],
                'index': None,
                'gene_caller_ids': [],
                'gene_caller_ids_in_genomes': {},
                'aa_sequences_in_gene_cluster': {},
                'previous_gene_cluster_name': None,
                'next_gene_cluster_name': None,
                'index': None,
                'total': None
                }

        if gene_cluster_name not in self.interactive.gene_clusters:
            return data

        if not self.interactive.genomes_storage_is_available:
            return data

        AA_sequences = self.interactive.get_sequences_for_gene_clusters(gene_cluster_names=set([gene_cluster_name]))

        # add the list of gene caller ids associated with this gene cluster into `data`:
        for genome_name in self.interactive.gene_clusters[gene_cluster_name]:
            data['aa_sequences_in_gene_cluster'][genome_name] = {}
            for gene_callers_id in self.interactive.gene_clusters[gene_cluster_name][genome_name]:
                data['gene_caller_ids'].append((gene_callers_id, genome_name), )
                data['aa_sequences_in_gene_cluster'][genome_name][gene_callers_id] = AA_sequences[gene_cluster_name][genome_name][gene_callers_id]

        # the dict that explains the distribution of genes in genomes:
        data['gene_caller_ids_in_genomes'] = self.interactive.gene_clusters[gene_cluster_name]

        # add the list of genomes into data:
        data['genomes'] = sorted(data['gene_caller_ids_in_genomes'].keys())

        # get some contextual stuff
        data['index'], data['total'], data['previous_gene_cluster_name'], data['next_gene_cluster_name'] = self.get_index_total_previous_and_next_items(order_name, gene_cluster_name)

        return json.dumps(data)


    def search_functions(self):
        try:
            full_report = self.interactive.search_for_functions(request.forms.get('terms'))
            return json.dumps({'status': 0, 'results': full_report})
        except Exception as e:
            message = str(e.clear_text()) if hasattr(e, 'clear_text') else str(e)
            return json.dumps({'status': 1, 'message': message})


    def store_refined_bins(self):
        data = json.loads(request.forms.get('data'))
        colors = json.loads(request.forms.get('colors'))

        bins_info_dict = {}
        for bin_name in data:
            bins_info_dict[bin_name] = {'html_color': colors[bin_name], 'source': "anvi-refine"}

        try:
            self.interactive.store_refined_bins(data, bins_info_dict)
        except RefineError as e:
            return json.dumps({'status': -1, 'message': e.clear_text()})

        message = 'Done! Collection %s is updated in the database. You can close your browser window (or continue updating).' % (self.interactive.collection_name)
        return json.dumps({'status': 0, 'message': message})


    def get_available_phylogeny_programs(self):
        return json.dumps(list(drivers.driver_modules['phylogeny'].keys()))

    def get_available_aligners(self):
        return json.dumps(list(drivers.Aligners().aligners.keys()))

    def generate_tree(self):
        gene_clusters = set(request.forms.getall('gene_clusters[]'))
        name = request.forms.get('name')
        program = request.forms.get('program')
        aligner = request.forms.get('aligner')
        store_tree = request.forms.get('store_tree')

        temp_fasta_file = filesnpaths.get_temp_file_path()
        temp_tree_file = filesnpaths.get_temp_file_path()
        tree_text = None

        try:
            self.interactive.write_sequences_in_gene_clusters_for_phylogenomics(gene_cluster_names=gene_clusters, output_file_path=temp_fasta_file, align_with=aligner)
            drivers.driver_modules['phylogeny'][program]().run_command(temp_fasta_file, temp_tree_file)
            tree_text = open(temp_tree_file,'rb').read().decode()

            if store_tree:
                TableForLayerOrders(self.interactive.args).add({name: {'data_type': 'newick', 'data_value': tree_text}})

                # TO DO: instead of injecting new newick tree, we can use TableForLayerOrders.get()
                self.interactive.layers_order_data_dict[name] = {'newick': tree_text, 'basic': None}
        except Exception as e:
            message = str(e.clear_text()) if 'clear_text' in dir(e) else str(e)
            return json.dumps({'status': 1, 'message': message})

        return json.dumps({'status': 0, 'tree': tree_text})


    def upload_project(self):
        try:
            args = argparse.Namespace()
            args.user = request.forms.get('username')
            args.password = request.forms.get('password')
            args.api_url = anvio.D['api-url'][1]['default']
            args.project_name = request.forms.get('project_name')
            args.delete_if_exists = True if request.forms.get('delete_if_exists') == "true" else False

            view_name = request.forms.get('view')
            if view_name in self.interactive.views:
                view_path = filesnpaths.get_temp_file_path()
                utils.store_array_as_TAB_delimited_file(self.interactive.views[view_name][1:], view_path, self.interactive.views[view_name][0])
                args.view_data = view_path

            item_order_name = request.forms.get('ordering')
            if item_order_name in self.interactive.p_meta['item_orders']:
                ordering_path = filesnpaths.get_temp_file_path()
                items_order = self.interactive.p_meta['item_orders'][item_order_name]

                f = open(ordering_path, 'w')
                if items_order['type'] == 'newick':
                    f.write(items_order['data'])
                    args.tree = ordering_path
                elif items_order['type'] == 'basic':
                    f.write("\n".join(items_order['data']))
                    args.items_order = ordering_path
                f.close()

            state_name = request.forms.get('state')
            if state_name in self.interactive.states_table.states:
                state_path = filesnpaths.get_temp_file_path()
                f = open(state_path, 'w')
                f.write(self.interactive.states_table.states[state_name]['content'])
                f.close()

                args.state = state_path

            if request.forms.get('include_description') == "true":
                description_path = filesnpaths.get_temp_file_path()
                f = open(description_path, 'w')
                f.write(self.interactive.p_meta['description'])
                f.close()

                args.description = description_path

            if request.forms.get('include_samples') == "true":
                # FIXME: this will break
                if len(self.interactive.layers_order_data_dict):
                    layers_order_data_path = filesnpaths.get_temp_file_path()
                    utils.store_dict_as_TAB_delimited_file(self.interactive.layers_order_data_dict, layers_order_data_path, headers=['attributes', 'basic', 'newick'])
                    args.layers_order_data_path = layers_order_data_path

                if len(self.interactive.layers_additional_data_dict):
                    layers_additional_data_path = filesnpaths.get_temp_file_path()
                    utils.store_dict_as_TAB_delimited_file(self.interactive.layers_additional_data_dict, layers_additional_data_path)
                    args.layers_additional_data_file = layers_additional_data_path

            collection_name = request.forms.get('collection')
            if collection_name in self.interactive.collections.collections_dict:
                collection_path_prefix = filesnpaths.get_temp_file_path()
                self.interactive.collections.export_collection(collection_name, output_file_prefix=collection_path_prefix)

                args.bins = collection_path_prefix + '.txt'
                args.bins_info = collection_path_prefix + '-info.txt'

            server = AnviServerAPI(args)
            server.login()
            server.push()
            return json.dumps({'status': 0})
        except Exception as e:
            message = str(e.clear_text()) if hasattr(e, 'clear_text') else str(e)
            return json.dumps({'status': 1, 'message': message})


    def get_contigs_stats(self):
        return json.dumps({'stats': self.interactive.contigs_stats,
                           'tables': self.interactive.tables,
                           'human_readable_keys': self.interactive.human_readable_keys})


    def filter_gene_clusters(self):
        try:
            parameters = {}
            for key in request.forms:
                parameters[key] = int(request.forms.get(key))

            gene_clusters_dict, _ = self.interactive.filter_gene_clusters_from_gene_clusters_dict(copy.deepcopy(self.interactive.gene_clusters), **parameters)
            return json.dumps({'status': 0, 'gene_clusters_list': list(gene_clusters_dict.keys())})
        except Exception as e:
            message = str(e.clear_text()) if hasattr(e, 'clear_text') else str(e)
            return json.dumps({'status': 1, 'message': message})
