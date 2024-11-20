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
import math
import copy
import time
import json
import base64
import random
import getpass
import argparse
import datetime
import importlib

from hashlib import md5
from ete3 import Tree
from bottle import Bottle
from bottle import BaseRequest
from bottle import redirect, static_file

# multiprocess is a fork of multiprocessing that uses the dill serializer instead of pickle
# using the multiprocessing module directly results in a pickling error in Python 3.10 which
# goes like this:
#
#   >>> AttributeError: Can't pickle local object 'SOMEFUNCTION.<locals>.<lambda>' multiprocessing
#
import multiprocess as multiprocessing

import anvio
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.drivers as drivers
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.summarizer as summarizer
import anvio.filesnpaths as filesnpaths
import anvio.taxonomyops.scg as scgtaxonomyops
import anvio.auxiliarydataops as auxiliarydataops

from anvio.serverAPI import AnviServerAPI
from anvio.errors import RefineError, ConfigError
from anvio.tables.miscdata import TableForLayerOrders
from anvio.tables.collections import TablesForCollections


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
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
    def __init__(self, interactive, mock_request=None, mock_response=None):
        super(BottleApplication, self).__init__()
        self.interactive = interactive

        # WSGI for bottle to use
        self._wsgi_for_bottle = "paste"

        self.additional_gc_data = None

        A = lambda x: self.args.__dict__[x] if x in self.args.__dict__ else None

        if self.interactive:
            self.args = self.interactive.args
            self.read_only = A('read_only')
            self.browser_path = A('browser_path')
            self.export_svg = A('export_svg')
            self.server_only = A('server_only')
            self.user_server_shutdown = A('user_server_shutdown')
            self.password_protected = A('password_protected')
            self.password = ''
            self.authentication_secret = ''
            if self.password_protected:
                print('')
                self.password = getpass.getpass('Enter password to secure interactive interface: ').encode('utf-8')
                salt = 'using_md5_in_2018_'.encode('utf-8')

                self.authentication_secret = md5(salt + self.password).hexdigest()

        self.session_id = random.randint(0,9999999999)
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

        # if there is a contigs database, and scg taxonomy was run on it get an instance
        # of the SCG Taxonomy class early on:
        if A('contigs_db') and dbops.ContigsDatabase(A('contigs_db')).meta['scg_taxonomy_was_run']:
            self.scg_taxonomy = scgtaxonomyops.SCGTaxonomyEstimatorSingle(argparse.Namespace(contigs_db=self.interactive.contigs_db_path))
        else:
            self.scg_taxonomy = None


    def set_password(self, password):
        self.password_protected = True
        self.password = password.encode('utf-8')
        salt = 'using_md5_in_2018_'.encode('utf-8')

        self.authentication_secret = md5(salt + self.password).hexdigest()


    def register_hooks(self):
        self.add_hook('before_request', self.before_request)


    def before_request(self):
        # /app/ contains static files and not password protected.
        if self.password_protected and not request.path.startswith('/app/'):
            if not self.authentication_secret == request.get_cookie('authentication_secret'):
                redirect('/app/login.html')

        response.set_header('Content-Type', 'application/json')
        response.set_header('Pragma', 'no-cache')
        response.set_header('Cache-Control', 'no-cache, no-store, max-age=0, must-revalidate')
        response.set_header('Expires', 'Thu, 01 Dec 1994 16:00:00 GMT')


    def register_routes(self):
        self.route('/',                                        callback=self.redirect_to_app)
        self.route('/app/:filename#.*#',                       callback=self.send_static)
        self.route('/app/shutdown',                            callback=self.server_shutdown)
        self.route('/data/news',                               callback=self.get_news)
        self.route('/data/<name>',                             callback=self.send_data)
        self.route('/data/view/<view_id>',                     callback=self.get_view_data)
        self.route('/tree/<items_order_id>',                   callback=self.get_items_order)
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
        self.route('/summarize/<collection_name>',             callback=self.gen_summary, method='POST')
        self.route('/summary/<collection_name>/:filename#.*#', callback=self.send_summary_static)
        self.route('/data/gene/<gene_callers_id>',             callback=self.get_sequence_for_gene_call)
        self.route('/data/hmm/<bin_name>/<gene_name>',         callback=self.get_hmm_hit_from_bin)
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
        self.route('/data/get_initial_data',                   callback=self.get_initial_data)
        self.route('/data/get_column_info',                    callback=self.get_column_info, method='POST')
        self.route('/data/get_structure/<gene_callers_id:int>',callback=self.get_structure)
        self.route('/data/get_variability',                    callback=self.get_variability, method='POST')
        self.route('/data/store_variability',                  callback=self.store_variability, method='POST')
        self.route('/data/store_structure_as_pdb',             callback=self.store_structure_as_pdb, method='POST')
        self.route('/data/get_gene_function_info/<gene_callers_id:int>',             callback=self.get_gene_function_info)
        self.route('/data/get_model_info/<gene_callers_id:int>',             callback=self.get_model_info)
        self.route('/data/filter_gene_clusters',               callback=self.filter_gene_clusters, method='POST')
        self.route('/data/reroot_tree',                        callback=self.reroot_tree, method='POST')
        self.route('/data/save_tree',                          callback=self.save_tree, method='POST')
        self.route('/data/check_homogeneity_info',             callback=self.check_homogeneity_info, method='POST')
        self.route('/data/get_additional_gc_data/<gc_id>/<gc_key>',    callback=self.get_additional_gc_data, method='POST')
        self.route('/data/search_items',                       callback=self.search_items_by_name, method='POST')
        self.route('/data/get_taxonomy',                       callback=self.get_taxonomy, method='POST')
        self.route('/data/get_functions_for_gene_clusters',    callback=self.get_functions_for_gene_clusters, method='POST')
        self.route('/data/get_gene_info/<gene_callers_id>',    callback=self.get_gene_info)
        self.route('/data/get_metabolism',                     callback=self.get_metabolism)
        self.route('/data/get_scale_bar',                      callback=self.get_scale_bar, method='POST')
        self.route('/data/get_psgc_data/<psgc_name>',          callback=self.get_psgc_data)


    def run_application(self, ip, port):
        # check for the wsgi module bottle will use.
        if not importlib.util.find_spec(self._wsgi_for_bottle):
            raise ConfigError("Anvi'o uses `%(wsgi)s` as a web server gateway interface, and you don't seem to have it. Which "
                              "means bad news. But the good news is that you can actually install it very easily. If you are "
                              "in a conda environment, try 'conda install %(wsgi)s'. If you are in a Python environment "
                              "try 'pip install %(wsgi)s'. If you are not sure, start with conda, if it doesn't work, try pip." \
                                    % {'wsgi': self._wsgi_for_bottle})

        try:
            # allow output to terminal when debugging
            if anvio.DEBUG:
                server_process = multiprocessing.Process(target=self.run, kwargs={'host': ip, 'port': port, 'quiet': False, 'server': self._wsgi_for_bottle})
                server_process.start()
            else:
                with terminal.SuppressAllOutput():
                    server_process = multiprocessing.Process(target=self.run, kwargs={'host': ip, 'port': port, 'quiet': True, 'server': self._wsgi_for_bottle})
                    server_process.start()

            url = "http://%s:%d" % (ip, port)

            if self.export_svg:
                try:
                    utils.run_selenium_and_export_svg("/".join([url, "app/index.html"]),
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

                utils.open_url_in_browser(url=url, browser_path=self.browser_path, run=run)

                run.info_single("The server is up and running 🎉", mc='green', nl_before = 1)

                run.warning("If you are using OSX and if the server terminates prematurely before you can see anything in your browser, "
                            "try running the same command by putting 'sudo ' at the beginning of it (you will be prompted to enter your "
                            "password if sudo requires super user credentials on your system). If your browser does not show up, try "
                            "manually entering the URL shown below into the address bar of your favorite browser. *cough* CHROME *cough*.")

            run.info('Server address', url, mc="green", nl_before=1, nl_after=1)

            run.info_single("When you are ready, press CTRL+C once to terminate the server and go back to the command line.", nl_after=1)

            server_process.join()
        except KeyboardInterrupt:
            run.info_single("The server is being terminated...", mc="red", nl_before=1)
            server_process.terminate()
            sys.exit(0)


    def redirect_to_app(self):
        homepage = 'index.html'
        if self.interactive.mode == 'contigs':
            homepage = 'contigs.html'
        elif self.interactive.mode == 'structure':
            homepage = 'structure.html'
        elif self.interactive.mode == 'metabolism':
            homepage = 'metabolism.html'
        elif self.interactive.mode == 'inspect':
            redirect('/app/charts.html?id=%s&show_snvs=true&rand=%s' % (self.interactive.inspect_split_name, self.random_hash(8)))

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


    def server_shutdown(self, **kwd):
        if self.user_server_shutdown:
            run.info_single('User Requested shutdown via web.', nl_after=1)
            # Could do sys.exit(0) instead, but raising KeyboardInterrupt will force consistent shutdown process
            raise KeyboardInterrupt
        return json.dumps({'error': "The server cannot be shutdown by a web user.", 'status_code': 0})


    def get_news(self):
        if self.interactive.anvio_news:
            return json.dumps(self.interactive.anvio_news)
        else:
            return json.dumps([{'date': '',
                                'title': 'No news for you :(',
                                'content': "Anvi'o couldn't bring any news for you. You can bring yourself to the news by clicking [here](%s)." % constants.anvio_news_url}])


    def random_hash(self, size=8):
        r = random.getrandbits(size * 4)
        return '{1:0{0}x}'.format(size, r)


    def send_data(self, name):
        if name == "init":
            bin_prefix = "Bin_"
            if self.interactive.mode == 'refine':
                bin_prefix = list(self.interactive.bin_names_of_interest)[0] + "_" if len(self.interactive.bin_names_of_interest) == 1 else "Refined_",

            default_view = self.interactive.default_view
            default_order = self.interactive.p_meta['default_item_order']
            autodraw = False
            state_dict = None

            if self.interactive.state_autoload:
                state_dict = json.loads(self.interactive.states_table.states[self.interactive.state_autoload]['content'])

                if 'current-view' in state_dict and state_dict['current-view'] in self.interactive.views:
                    default_view = state_dict['current-view']

                if 'order-by' in state_dict and state_dict['order-by'] in self.interactive.p_meta['item_orders']:
                    default_order = state_dict['order-by']

                autodraw = True

            collection_dict = None
            if self.interactive.mode != 'collection' and self.interactive.mode != 'refine' and self.interactive.collection_autoload:
                collection_dict = json.loads(self.get_collection_dict(self.interactive.collection_autoload))
                autodraw = True

            item_lengths = {}
            if self.interactive.mode == 'full' or self.interactive.mode == 'refine':
                item_lengths = dict([tuple((c, self.interactive.splits_basic_info[c]['length']),) for c in self.interactive.splits_basic_info])
            elif self.interactive.mode == 'pan':
                for gene_cluster in self.interactive.gene_clusters:
                    item_lengths[gene_cluster] = 0
                    for genome in self.interactive.gene_clusters[gene_cluster]:
                        item_lengths[gene_cluster] += len(self.interactive.gene_clusters[gene_cluster][genome])

            functions_sources = []
            if self.interactive.mode == 'full' or self.interactive.mode == 'gene' or self.interactive.mode == 'refine':
                functions_sources = list(self.interactive.gene_function_call_sources)
            elif self.interactive.mode == 'pan':
                functions_sources = list(self.interactive.gene_clusters_function_sources)

            inspection_available = self.interactive.auxiliary_profile_data_available

            return json.dumps( { "version":                            anvio.anvio_version,
                                 "title":                              self.interactive.title,
                                 "description":                        self.interactive.p_meta['description'],
                                 "item_orders":                        (default_order, self.interactive.p_meta['item_orders'][default_order], list(self.interactive.p_meta['item_orders'].keys())),
                                 "views":                              (default_view, self.interactive.views[default_view], list(self.interactive.views.keys())),
                                 "item_lengths":                       item_lengths,
                                 "mode":                               self.interactive.mode,
                                 "server_mode":                        False,
                                 "read_only":                          self.read_only,
                                 "bin_prefix":                         bin_prefix,
                                 "session_id":                         self.session_id,
                                 "layers_order":                       self.interactive.layers_order_data_dict,
                                 "layers_information":                 self.interactive.layers_additional_data_dict,
                                 "layers_information_default_order":   self.interactive.layers_additional_data_keys,
                                 "check_background_process":           True,
                                 "autodraw":                           autodraw,
                                 "inspection_available":               inspection_available,
                                 "sequences_available":                True if (self.interactive.split_sequences or self.interactive.mode == 'gene') else False,
                                 "functions_initialized":              self.interactive.gene_function_calls_initiated,
                                 "functions_sources":                  functions_sources,
                                 "state":                              (self.interactive.state_autoload, state_dict),
                                 "collection":                         collection_dict,
                                 "samples":                            self.interactive.p_meta['samples'] if self.interactive.mode in ['full', 'refine'] else [],
                                 "load_full_state":                    self.interactive.load_full_state })

        elif name == "session_id":
            return json.dumps(self.session_id)


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

            return json.dumps(items_order)

        return json.dumps("")


    def save_tree(self):
        try:
            order_full_name = request.forms.get('name')
            order_data = request.forms.get('data')
            tree_type = request.forms.get('tree_type')
            additional = request.forms.get('additional')

            if tree_type == 'samples':
                order_name = order_full_name
                distance = 'NA'
                linkage = 'NA'

                if order_name in self.interactive.layers_order_data_dict:
                    raise ConfigError("Tree name '%s' already exists, overwriting currently not supported." % order_name)

                self.interactive.layers_order_data_dict[order_name] = {'newick': order_data, 'basic': ''}
                TableForLayerOrders(self.interactive.args).add({order_name: {'data_type': 'newick', 'data_value': order_data}})
            else:
                self.interactive.p_meta['item_orders'][order_full_name] = {'type': 'newick', 'data': order_data, 'additional': additional}

                order_name, distance, linkage = order_full_name.split(':')
                anvio_db_path = self.interactive.pan_db_path or self.interactive.profile_db_path

                dbops.add_items_order_to_db(anvio_db_path, order_name, order_data, order_data_type_newick=True, distance=distance, linkage=linkage, additional_data=additional, dont_overwrite=True)

            return json.dumps({'status': 0, 'message': 'New order "%s (D: %s; L: %s)" successfully saved to the database.' % (order_name, distance, linkage)})

        except ConfigError as e:
            message = str(e.clear_text()) if hasattr(e, 'clear_text') else str(e)
            return json.dumps({'status': 1, 'message': message})


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
            state_dict = json.loads(state['content'])

            if self.interactive.mode == 'structure':
                return json.dumps({'content': state['content']})
            else:

                default_view = self.interactive.default_view
                default_order = self.interactive.p_meta['default_item_order']

                if 'current-view' in state_dict and state_dict['current-view'] in self.interactive.views: # extra checks to accomodate minified states without current-view/order-by data
                    default_view = state_dict['current-view']

                if 'order-by' in state_dict and state_dict['order-by'] in self.interactive.p_meta['item_orders']:
                    default_order = state_dict['order-by']

                return json.dumps((state_dict, self.interactive.p_meta['item_orders'][default_order], self.interactive.views[default_view]))

        return json.dumps("")


    def charts(self, order_name, item_name):
        title = None
        state = {}

        if self.interactive.mode == 'gene':
            split_name = self.interactive.gene_callers_id_to_split_name_dict[int(item_name)]
            title = "Gene '%d' in split '%s'" % (int(item_name), split_name)
        else:
            split_name = item_name
            title = split_name

        if self.interactive.mode == 'inspect':
            order_name = 'alphabetical'

        data = {'layers': [],
                 'title': title,
                 'index': None,
                 'total': None,
                 'coverage': [],
                 'variability': [],
                 'indels': [],
                 'competing_nucleotides': [],
                 'previous_contig_name': None,
                 'next_contig_name': None,
                 'genes': [],
                 'outlier_SNVs_shown': not self.args.hide_outlier_SNVs,
                 'state': {}}

        if split_name not in self.interactive.split_names:
            return data

        if not self.interactive.auxiliary_profile_data_available:
            return data

        data['index'], data['total'], data['previous_contig_name'], data['next_contig_name'] = self.get_index_total_previous_and_next_items(order_name, item_name)

        if self.interactive.mode == 'inspect':
            if self.interactive.state_autoload:
                state = json.loads(self.interactive.states_table.states[self.interactive.state_autoload]['content'])
                layers = [layer for layer in sorted(self.interactive.p_meta['samples']) if (layer not in state['layers'] or float(state['layers'][layer]['height']) > 0)]
            else:
                layers = [layer for layer in sorted(self.interactive.p_meta['samples'])]

                # anvi-inspect is called so there is no state stored in localstorage written by main anvio plot
                # and there is no default state in the database, we are going to generate a mock state.
                # only the keys we need is enough.
                state['layer-order'] = layers
                state['layers'] = {}
                for layer in layers:
                    state['layers'][layer] = {'height': 1, 'color': '#00000'}

        else:
            state = json.loads(request.forms.get('state'))
            layers = [layer for layer in sorted(self.interactive.p_meta['samples']) if (layer not in state['layers'] or float(state['layers'][layer]['height']) > 0)]

        data['state'] = state

        db_variant = str(self.interactive.contigs_db_variant)
        try:
            auxiliary_coverages_db = auxiliarydataops.AuxiliaryDataForSplitCoverages(self.interactive.auxiliary_data_path,
                                                                                     self.interactive.p_meta['contigs_db_hash'],
                                                                                     db_variant=db_variant)
            coverages = auxiliary_coverages_db.get(split_name)
            auxiliary_coverages_db.close()

            data['coverage'] = [coverages[layer].tolist() for layer in layers]
        except:
            data['coverage'] = [[0] * self.interactive.splits_basic_info[split_name]['length']]

        data['sequence'] = self.interactive.split_sequences[split_name]['sequence']

        ## get the variability information dict for split:
        split_variability_info_dict = self.interactive.get_variability_information_for_split(split_name, skip_outlier_SNVs=self.args.hide_outlier_SNVs)

        ## get the indels information dict for split:
        split_indels_info_dict = self.interactive.get_indels_information_for_split(split_name)

        # building layer data
        for layer in layers:
            data['layers'].append(layer)
            data['competing_nucleotides'].append(split_variability_info_dict[layer]['competing_nucleotides'])
            data['variability'].append(split_variability_info_dict[layer]['variability'])
            data['indels'].append(split_indels_info_dict[layer]['indels'])

        levels_occupied = {1: []}
        gene_entries_in_split = self.interactive.split_name_to_genes_in_splits_entry_ids[split_name]

        # we get all the amino acid sequences for genes in this split here to avoid
        # multiple database calls for each gene later. this dictionary will be used
        # below as we go through each gene call.
        gene_aa_sequences_dict = self.interactive.get_gene_amino_acid_sequence([self.interactive.genes_in_splits[e]['gene_callers_id'] for e in gene_entries_in_split])

        for entry_id in gene_entries_in_split:
            gene_callers_id =  self.interactive.genes_in_splits[entry_id]['gene_callers_id']

            # this is a CRAZY case where a gene caller id is found in a split, but
            # it is not occurring in the genes table. ABSOLUTELY CRAZY, BUT FLORIAN
            # MANAGED TO DO IT, SO THERE WE GO.
            if gene_callers_id not in self.interactive.genes_in_contigs_dict:
                progress.reset()
                run.info_single(f"Gene caller id {gene_callers_id} is missing from the contigs db. But the "
                                f"split {split_name} thinks it has it :/ Anvi'o ignores this. Anvi'o is too old "
                                f"for stuff like this.")
                continue

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
            p['call_type'] =  self.interactive.genes_in_contigs_dict[gene_callers_id]['call_type']
            p['complete_gene_call'] = 'No' if  self.interactive.genes_in_contigs_dict[gene_callers_id]['partial'] else 'Yes'
            p['length'] = p['stop_in_contig'] - p['start_in_contig']
            p['functions'] =  self.interactive.gene_function_calls_dict[gene_callers_id] if gene_callers_id in  self.interactive.gene_function_calls_dict else None

            # get amino acid sequence for the gene call:
            p['aa_sequence'] = gene_aa_sequences_dict[gene_callers_id]

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

        return json.dumps(data)


    def get_gene_info(self, gene_callers_id):
        # TO DO: there are three functions returns gene info dict for different purposes
        # this needs to be organized.

        gene_callers_id = int(gene_callers_id)
        split_name = self.interactive.gene_callers_id_to_split_name_dict[gene_callers_id]

        # we need eentry id
        entry_id = None
        for candidate_entry_id in self.interactive.split_name_to_genes_in_splits_entry_ids[split_name]:
            if int(gene_callers_id) == int(self.interactive.genes_in_splits[candidate_entry_id]['gene_callers_id']):
                entry_id = candidate_entry_id

        if not entry_id:
            raise ConfigError("Can not find this gene_callers_id in any splits.")

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
        p['call_type'] =  self.interactive.genes_in_contigs_dict[gene_callers_id]['call_type']
        p['complete_gene_call'] = 'No' if  self.interactive.genes_in_contigs_dict[gene_callers_id]['partial'] else 'Yes'
        p['length'] = p['stop_in_contig'] - p['start_in_contig']
        p['functions'] =  self.interactive.gene_function_calls_dict[gene_callers_id] if gene_callers_id in  self.interactive.gene_function_calls_dict else None

        # get amino acid sequence for the gene call:
        p['aa_sequence'] = self.interactive.get_gene_amino_acid_sequence([gene_callers_id])

        return json.dumps(p)


    def search_items_by_name(self):
        items_per_page = 30

        query = request.forms.get('search-query')
        page = int(request.forms.get('page') or 0)

        if query and len(query) > 0:
            query = query.lower()
            results = []
            for name in self.interactive.displayed_item_names_ordered:
                if query in name.lower():
                    results.append(name)
        else:
            results = self.interactive.displayed_item_names_ordered

        page_start = max(0, page * items_per_page)
        page_end = min(len(results), (page + 1) * items_per_page)

        total_page = math.ceil(len(results) / items_per_page)

        results = results[page_start:page_end]

        return json.dumps({
            'search-query': query,
            'results': results,
            'page': page,
            'total_page': total_page
            })


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
                 'indels': [],
                 'competing_nucleotides': [],
                 'previous_contig_name': None,
                 'next_contig_name': None,
                 'genes': [],
                 'outlier_SNVs_shown': not self.args.hide_outlier_SNVs,
                 'state': state}

        data['index'], data['total'], data['previous_contig_name'], data['next_contig_name'] = self.get_index_total_previous_and_next_items(order_name, str(gene_callers_id))

        layers = [layer for layer in sorted(self.interactive.p_meta['samples']) if (layer not in state['layers'] or float(state['layers'][layer]['height']) > 0)]

        db_variant = str(self.interactive.contigs_db_variant)
        auxiliary_coverages_db = auxiliarydataops.AuxiliaryDataForSplitCoverages(self.interactive.auxiliary_data_path,
                                                                                 self.interactive.p_meta['contigs_db_hash'],
                                                                                 db_variant=db_variant)
        coverages = auxiliary_coverages_db.get(split_name)
        auxiliary_coverages_db.close()

        data['coverage'] = []
        for layer in layers:
            coverage_list = coverages[layer].tolist()
            data['coverage'].append(coverage_list[focus_region_start:focus_region_end])

        data['sequence'] = self.interactive.split_sequences[split_name]['sequence'][focus_region_start:focus_region_end]

        ## get the variability information dict for split:
        progress.new('Variability')
        progress.update('Collecting info for "%s"' % split_name)
        split_variability_info_dict = self.interactive.get_variability_information_for_split(split_name, skip_outlier_SNVs=self.args.hide_outlier_SNVs)

        for layer in layers:
            progress.update('Formatting variability data: "%s"' % layer)
            data['layers'].append(layer)

            # filter and substract offset variability, and competing nucleotide information.
            variability_dict_original = copy.deepcopy(split_variability_info_dict[layer]['variability'])
            variability_dict = {}
            for nucleotide_pos_in_codon in variability_dict_original:
                variability_dict[nucleotide_pos_in_codon] = {}
                for pos in variability_dict_original[nucleotide_pos_in_codon]:
                    if pos < focus_region_start or pos > focus_region_end:
                        continue

                    variability_dict[nucleotide_pos_in_codon][pos - focus_region_start] = variability_dict_original[nucleotide_pos_in_codon][pos]

            competing_nucleotides_dict_original = copy.deepcopy(split_variability_info_dict[layer]['competing_nucleotides'])
            competing_nucleotides_dict = {}
            for pos in competing_nucleotides_dict_original:
                if pos < focus_region_start or pos > focus_region_end:
                    continue

                entry = competing_nucleotides_dict_original[pos]
                entry['pos_in_split'] = pos
                entry['pos'] = pos - focus_region_start

                competing_nucleotides_dict[entry['pos']] = entry

            data['competing_nucleotides'].append(competing_nucleotides_dict)
            data['variability'].append(variability_dict)

        progress.end()

        ## get the indels information dict for split:
        progress.new('Indels')
        progress.update('Collecting info for "%s"' % split_name)
        split_indels_info_dict = self.interactive.get_indels_information_for_split(split_name)

        for layer in layers:
            progress.update('Formatting indels data: "%s"' % layer)

            indels_dict_original = copy.deepcopy(split_indels_info_dict[layer]['indels'])
            indels_dict = {}
            for indel_entry_id in indels_dict_original:
                pos = indels_dict_original[indel_entry_id]['pos']
                if pos < focus_region_start or pos > focus_region_end:
                    continue
                else:
                    indels_dict[indel_entry_id] = indels_dict_original[indel_entry_id]
                    indels_dict[indel_entry_id]['pos'] = indels_dict[indel_entry_id]['pos'] - focus_region_start

            data['indels'].append(indels_dict)

        levels_occupied = {1: []}
        gene_entries_in_split = self.interactive.split_name_to_genes_in_splits_entry_ids[split_name]

        # we get all the amino acid sequences for genes in this split here to avoid
        # multiple database calls for each gene later. this dictionary will be used
        # below as we go through each gene call.
        gene_aa_sequences_dict = self.interactive.get_gene_amino_acid_sequence([self.interactive.genes_in_splits[e]['gene_callers_id'] for e in gene_entries_in_split])

        for entry_id in gene_entries_in_split:
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

            # add amino acid sequence for the gene call:
            p['aa_sequence'] = gene_aa_sequences_dict[gene_callers_id]

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
            names_with_only_digits_ok = self.interactive.mode == 'gene'
            items_order = utils.get_names_order_from_newick_tree(items_order_entry['data'], names_with_only_digits_ok=names_with_only_digits_ok)
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
        completeness_data = {}
        if not self.interactive.completeness:
            return json.dumps(completeness_sources)

        split_names = json.loads(request.forms.get('split_names'))
        bin_name = json.loads(request.forms.get('bin_name'))

        run.info_single('Completeness info has been requested for %d splits in %s' % (len(split_names), bin_name))

        p_completion, p_redundancy, scg_domain, domain_probabilities, info_text, results_dict = self.interactive.completeness.get_info_for_splits(set(split_names))

        # convert results_dict (where domains are the highest order items) into a dict that is compatible with the
        # previous format of the dict (where hmm scg source names are the higher order items).
        for domain in results_dict:
            for source in results_dict[domain]:
                completeness_sources[source] = results_dict[domain][source]

        completeness_data['percent_completion'] = p_completion
        completeness_data['percent_redundancy'] = p_redundancy
        completeness_data['domain'] = scg_domain
        completeness_data['info_text'] = info_text
        completeness_data['domain_probabilities'] = domain_probabilities

        # FIXME: We need to look into what we are sending and sort out what needs to be shown:
        return json.dumps({'stats': completeness_sources, 'averages': completeness_data, 'refs': self.interactive.completeness.http_refs})


    def get_collections(self):
        csd = self.interactive.collections.collections_dict
        run.info_single('Collection sources has been requested (info dict with %d item(s) has been returned).' % len(csd), cut_after=None)
        return json.dumps(csd)


    def get_collection_dict(self, collection_name):
        run.info_single('Data for collection "%s" has been requested.' % collection_name)

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

        run.info_single('A request to store %d bins that describe %d splits under the collection id "%s" '
                        'has been made.' % (len(data), num_splits, source), cut_after=None)

        bins_info_dict = {}
        for bin_name in data:
            bins_info_dict[bin_name] = {'html_color': colors[bin_name], 'source': "anvi-interactive"}

        # the db here is either a profile db, or a pan db, but it can't be both:
        if self.interactive.mode == 'gene':
            db_path = self.interactive.genes_db_path
        else:
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
            sequence = self.interactive.split_sequences[split_name]['sequence']
            header = split_name
        except Exception as e:
            return json.dumps({'error': "Something went wrong when I tried to access that split sequence: '%s' :/" % e})

        return json.dumps({'sequence': sequence, 'header': header})


    def gen_summary(self, collection_name):
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
        init_gene_coverages = json.loads(request.forms.get('init_gene_coverages'))

        if self.interactive.mode == 'pan':
            summarizer_args.pan_db = self.interactive.pan_db_path
            summarizer_args.genomes_storage = self.interactive.genomes_storage_path
            summarizer_args.output_dir = os.path.join(os.path.dirname(summarizer_args.pan_db), 'SUMMARY_%s' % collection_name)
        elif self.interactive.mode == 'full':
            summarizer_args.profile_db = self.interactive.profile_db_path
            summarizer_args.contigs_db = self.interactive.contigs_db_path
            if init_gene_coverages:
                summarizer_args.init_gene_coverages = True

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
            gene_calls_tuple = self.interactive.get_sequences_for_gene_callers_ids([gene_callers_id], include_aa_sequences=True)
        except Exception as e:
            return json.dumps({'error': "Something went wrong when I tried to access to that gene: '%s' :/" % e})

        entry = gene_calls_tuple[1][gene_callers_id]

        sequence = entry['sequence']
        aa_sequence = entry['aa_sequence']

        header = '%d|' % (gene_callers_id) + '|'.join(['%s:%s' % (k, str(entry[k])) for k in ['contig', 'start', 'stop', 'direction', 'rev_compd', 'length']])

        return json.dumps({'sequence': sequence, 'aa_sequence': aa_sequence, 'header': header})


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

        if self.interactive.hmm_access is None:
            return json.dumps({'error': "HMMs for single-copy core genes were not run for this contigs database. "})

        hmm_sequences_dict = self.interactive.hmm_access.get_sequences_dict_for_hmm_hits_in_splits({bin_name: set(self.interactive.collection[bin_name])})
        gene_sequences = utils.get_filtered_dict(hmm_sequences_dict, 'gene_name', set([gene_name]))

        if not gene_sequences:
            return json.dumps({'error': "Sorry. It seems %s does not have a hit for %s." % (bin_name, gene_name)})

        unique_id_for_longest_hit = sorted([(gene_sequences[gene_id]['length'], gene_id) for gene_id in gene_sequences], reverse=True)[0][1]

        header, sequence = self.interactive.hmm_access.get_FASTA_header_and_sequence_for_gene_unique_id(gene_sequences, unique_id_for_longest_hit)

        return json.dumps({'sequence': sequence, 'header': header})


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
            requested_sources = request.forms.getall('sources[]')

            if not len(requested_sources):
                requested_sources = None

            items, full_report = self.interactive.search_for_functions(request.forms.get('terms'), requested_sources)

            items_unique = set([])
            for search_term in items:
                items_unique = items_unique.union(set(items[search_term]))

            return json.dumps({'status': 0, 'results': full_report, 'item_count': len(items_unique)})
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
        gene_cluster_names = set(request.forms.getall('gene_clusters[]'))
        gene_clusters = self.interactive.filter_gene_clusters_dict(argparse.Namespace(gene_clusters_names_of_interest=gene_cluster_names))
        name = request.forms.get('name')
        program = request.forms.get('program')
        aligner = request.forms.get('aligner')
        store_tree = request.forms.get('store_tree')

        temp_fasta_file = filesnpaths.get_temp_file_path()
        temp_tree_file = filesnpaths.get_temp_file_path()
        tree_text = None

        try:
            self.interactive.write_sequences_in_gene_clusters_for_phylogenomics(gene_clusters_dict=gene_clusters, output_file_path=temp_fasta_file, align_with=aligner)
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

            # FIX ME: this broke
            # if request.forms.get('include_samples') == "true":
            #     # FIXME: this will break
            #     if len(self.interactive.layers_order_data_dict):
            #         layers_order_data_path = filesnpaths.get_temp_file_path()
            #         print(layers_order_data_path)
            #         utils.store_dict_as_TAB_delimited_file(self.interactive.layers_order_data_dict, layers_order_data_path, headers=['attributes', 'basic', 'newick'])
            #         args.layers_order_data_path = layers_order_data_path

            #     if len(self.interactive.layers_additional_data_dict):
            #         layers_additional_data_path = filesnpaths.get_temp_file_path()
            #         print(layers_additional_data_path)
            #         utils.store_dict_as_TAB_delimited_file(self.interactive.layers_additional_data_dict, layers_additional_data_path)
            #         args.layers_additional_data_file = layers_additional_data_path

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


    def get_initial_data(self):
        return json.dumps(self.interactive.get_initial_data())


    def get_column_info(self):
        gene_callers_id = int(request.forms.get('gene_callers_id'))
        engine = request.forms.get('engine')

        return json.dumps(self.interactive.get_column_info(gene_callers_id, engine))


    def get_metabolism(self):
        return json.dumps(self.interactive.get_metabolism_data())


    def get_structure(self, gene_callers_id):
        return json.dumps(self.interactive.get_structure(gene_callers_id))


    def get_variability(self):
        options = json.loads(request.forms.get('options'))
        return self.interactive.get_variability(options)


    def store_variability(self):
        options = json.loads(request.forms.get('options'))
        return self.interactive.store_variability(options)


    def store_structure_as_pdb(self):
        options = json.loads(request.forms.get('options'))
        return self.interactive.store_structure_as_pdb(options)


    def get_gene_function_info(self, gene_callers_id):
        return json.dumps(self.interactive.get_gene_function_info(gene_callers_id))


    def get_model_info(self, gene_callers_id):
        return json.dumps(self.interactive.get_model_info(gene_callers_id))


    def filter_gene_clusters(self):
        try:
            parameters = {}
            for key in request.forms:
                parameters[key] = float(request.forms.get(key))

            gene_clusters_dict, _ = self.interactive.filter_gene_clusters_from_gene_clusters_dict(copy.deepcopy(self.interactive.gene_clusters), **parameters)
            return json.dumps({'status': 0, 'gene_clusters_list': list(gene_clusters_dict.keys())})
        except Exception as e:
            message = str(e.clear_text()) if hasattr(e, 'clear_text') else str(e)
            return json.dumps({'status': 1, 'message': message})


    def check_homogeneity_info(self):
        try:
            return json.dumps({'status': 0,
                               'functional_homogeneity_info_is_available': self.interactive.functional_homogeneity_info_is_available,
                               'geometric_homogeneity_info_is_available': self.interactive.geometric_homogeneity_info_is_available,
                               'combined_homogeneity_info_is_available': self.interactive.combined_homogeneity_info_is_available})
        except:
            return json.dumps({'status': 1})


    def get_additional_gc_data(self, gc_id, gc_key):
        try:
            return json.dumps({
                'status': 0,
                'gene_cluster_data': self.interactive.items_additional_data_dict[gc_id][gc_key]
            })
        except:
            return json.dumps({'status': 1})


    def get_psgc_data(self, psgc_name):
        """Gets PSGC data for the inspection page"""

        try:
            psgc_data = {
                psgc_name: []
            }

            associated_gcs = []
            for gc_id, psgc_id in self.interactive.gc_psgc_associations.items():
                if psgc_id == psgc_name:
                    associated_gcs.append(gc_id)
            
            for gene_id, gene_entries in self.interactive.gc_tracker.items():
                for gene_data in gene_entries:
                    gc_id = gene_data['gene_cluster_id']

                    # Skip if this gene's cluster is not associated with our PSGC
                    # Its important to do this before we add the gene info to psgc_data
                    if gc_id not in associated_gcs:
                        continue

                    gene_info = {
                        'gene_callers_id': gene_id,
                        'gene_cluster_id': gc_id,
                        'genome_name': gene_data['genome_name'],
                        'alignment_summary': gene_data['alignment_summary']
                    }
                    psgc_data[psgc_name].append(gene_info)

            return json.dumps({'status': 0, 'data': psgc_data})

        except Exception as e:
            return json.dumps({'status': 1, 'message': f"Error getting PSGC data: {str(e)}"})


    def reroot_tree(self):
        # Get the Newick tree string from the form data
        newick = request.forms.get('newick')
        internal_node = request.forms.get('internal_node')
        tree = Tree(newick, format=1)

        branch_support_value = {1:''}
        unique_index = 1

        if not internal_node:
        # for every node that is not leaf or root, associate with a unique index
            for node in tree.traverse():
                if not node.is_leaf() and not node.is_root():
                    unique_index += 1
                    branch_support_value[unique_index] = node.name
                    node.support = unique_index

        # Find the leftmost and rightmost nodes based on the provided names
        left_most = tree.search_nodes(name=request.forms.get('left_most'))[0]
        right_most = tree.search_nodes(name=request.forms.get('right_most'))[0]

        # Find the new root by identifying the common ancestor of the leftmost and rightmost nodes
        new_root = tree.get_common_ancestor(left_most, right_most)

        # Set the new root as the outgroup
        tree.set_outgroup(new_root)

        if not internal_node:
            # Assign support values as node names for non-leaf and non-root nodes
            for node in tree.traverse():
                if not node.is_leaf() and not node.is_root():
                    node.name = branch_support_value[node.support]

        # Encode node names using base32 encoding
        for node in tree.traverse('preorder'):
            node.name = 'base32' + base64.b32encode(node.name.encode('utf-8')).decode('utf-8')

        # Serialize the tree to Newick format
        new_newick = tree.write(format=1)

        # Decode base32 encoded node names back
        new_newick = re.sub(r"base32(\w*)", lambda m: base64.b32decode(m.group(1).replace('_', '=')).decode('utf-8'), new_newick)

        # Return the modified Newick format tree string as JSON
        return json.dumps({'newick': new_newick})


    def get_taxonomy(self):
        collection = json.loads(request.forms.get('collection'))

        if not self.scg_taxonomy:
            message = "You first need to run `anvi-run-scg-taxonomy` on your contigs database for this to work :("
            run.warning(message)
            return json.dumps({'status': 1, 'message': message})

        output = {}
        try:
            for bin_name in collection:
                output[bin_name] = self.scg_taxonomy.estimate_for_list_of_splits(collection[bin_name], bin_name=bin_name)

            run.info_single('Taxonomy estimation has been requested for bin(s) "%s".' % (", ".join(collection.keys())))
        except Exception as e:
            message = str(e.clear_text()) if hasattr(e, 'clear_text') else str(e)
            return json.dumps({'status': 1, 'message': message})

        return json.dumps(output)


    def get_functions_for_gene_clusters(self):
        if not len(self.interactive.gene_clusters_function_sources):
            message = "Gene cluster functions seem to have not been initialized, so that button has nothing to show you :/ Please carry on."
            run.warning(message)
            return json.dumps({'status': 1, 'message': message})

        gene_cluster_names = json.loads(request.forms.get('gene_clusters'))

        d = {}
        for gene_cluster_name in gene_cluster_names:
            if gene_cluster_name not in self.interactive.gene_clusters_functions_summary_dict:
                message = (f"At least one of the gene clusters in your list (e.g., {gene_cluster_name}) is missing in "
                           f"the functions summary dict :/")
                return json.dumps({'status': 1, 'message': message})
                
            d[gene_cluster_name] = self.interactive.gene_clusters_functions_summary_dict[gene_cluster_name]

        return json.dumps({'functions': d, 'sources': list(self.interactive.gene_clusters_function_sources)})


    def get_scale_bar(self):
        try:
            newick = request.json.get('newick')
            tree = Tree(newick, format=1)

            total_branch_length = tree.get_farthest_leaf()[1]

        except Exception as e:
            message = str(e.clear_text()) if hasattr(e, 'clear_text') else str(e)
            return json.dumps({'status': 1, 'message': message})
            
        return json.dumps({'scale_bar_value': total_branch_length})