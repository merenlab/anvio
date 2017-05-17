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
import json
import random
import requests
import datetime
import webbrowser
from multiprocessing import Process

from bottle import Bottle
from bottle import BaseRequest
from bottle import redirect, static_file

import anvio
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.summarizer as summarizer

from anvio.errors import RefineError, ConfigError


__author__ = "Ozcan Esen"
__copyright__ = "Copyright 2016, The anvio Project"
__credits__ = ["A. Murat Eren"]
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


run = terminal.Run()
progress = terminal.Progress()

# increase maximum size of form data to 100 MB
BaseRequest.MEMFILE_MAX = 1024 * 1024 * 100

class BottleApplication(Bottle):
    def __init__(self, interactive, args, mock_request=None, mock_response=None):
        super(BottleApplication, self).__init__()
        self.interactive = interactive
        self.args = args
        self.read_only = args.read_only

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
        self.route('/tree/<items_ordering_id>',                callback=self.get_items_ordering)
        self.route('/state/autoload',                          callback=self.state_autoload)
        self.route('/state/all',                               callback=self.state_all)
        self.route('/state/get/<state_name>',                  callback=self.get_state)
        self.route('/state/save/<state_name>',                 callback=self.save_state, method='POST')
        self.route('/data/charts/<split_name>',                callback=self.charts)
        self.route('/data/completeness',                       callback=self.completeness, method='POST')
        self.route('/data/collections',                        callback=self.get_collections)
        self.route('/data/collection/<collection_name>',       callback=self.get_collection_dict)
        self.route('/store_collection',                        callback=self.store_collections_dict, method='POST')
        self.route('/store_description',                       callback=self.store_description, method='POST')
        self.route('/data/contig/<split_name>',                callback=self.get_sequence_for_split)
        self.route('/summarize/<collection_name>',             callback=self.gen_summary)
        self.route('/summary/<collection_name>/:filename#.*#', callback=self.send_summary_static)
        self.route('/data/gene/<gene_callers_id>',             callback=self.get_sequence_for_gene_call)
        self.route('/data/hmm/<bin_name>/<gene_name>',         callback=self.get_hmm_hit_from_bin)
        self.route('/data/proteinclusterssummary',             callback=self.get_protein_clusters_summary, method='POST')
        self.route('/data/get_AA_sequences_for_PC/<pc_name>',  callback=self.get_AA_sequences_for_PC)
        self.route('/data/proteinclusters/<pc_name>',          callback=self.inspect_pc)
        self.route('/data/store_refined_bins',                 callback=self.store_refined_bins, method='POST')

    def run_application(self, ip, port):
        try:
            server_process = Process(target=self.run, kwargs={'host': ip, 'port': port, 'quiet': True, 'server': 'cherrypy'})
            server_process.start()

            if not self.args.server_only:
                webbrowser.open_new("http://%s:%d" % (ip, port))

            run.info_single('The server is now listening the port number "%d". When you are finished, press CTRL+C to terminate the server.' % port, 'green', nl_before = 1, nl_after=1)
            server_process.join()
        except KeyboardInterrupt:
            run.warning('The server is being terminated.', header='Please wait...')
            server_process.terminate()
            sys.exit(0)

    def redirect_to_app(self):
        redirect('/app/index.html?rand=' + self.random_hash(8))

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

            return json.dumps( { "title":                               self.interactive.title,
                                 "description":                        (self.interactive.p_meta['description']),
                                 "clusterings":                        (self.interactive.p_meta['default_clustering'], self.interactive.p_meta['clusterings']),
                                 "views":                              (self.interactive.default_view, dict(list(zip(list(self.interactive.views.keys()), list(self.interactive.views.keys()))))),
                                 "contigLengths":                      dict([tuple((c, self.interactive.splits_basic_info[c]['length']),) for c in self.interactive.splits_basic_info]),
                                 "defaultView":                        self.interactive.views[self.interactive.default_view],
                                 "mode":                               self.interactive.mode,
                                 "readOnly":                           self.read_only,
                                 "binPrefix":                          bin_prefix,
                                 "sessionId":                          self.unique_session_id,
                                 "samplesOrder":                       self.interactive.samples_order_dict,
                                 "sampleInformation":                  self.interactive.samples_information_dict,
                                 "sampleInformationDefaultLayerOrder": self.interactive.samples_information_default_layer_order,
                                 "stateAutoload":                      self.interactive.state_autoload,
                                 "collectionAutoload":                 self.interactive.collection_autoload,
                                 "noPing":                             False,
                                 "inspectionAvailable":                self.interactive.auxiliary_profile_data_available,
                                 "sequencesAvailable":                 True if self.interactive.split_sequences else False})

        elif name == "session_id":
            return json.dumps(self.unique_session_id)

    def get_view_data(self, view_id):
        return json.dumps(self.interactive.views[view_id])

    def get_items_ordering(self, items_ordering_id):
        if items_ordering_id in self.interactive.p_meta['clusterings']:
            items_ordering = self.interactive.p_meta['clusterings'][items_ordering_id]

            if 'newick' in items_ordering:
                run.info_single("The newick order '%s' has been requested" % (items_ordering_id))
                return json.dumps(items_ordering['newick'])
            elif 'basic' in items_ordering:
                run.info_single("The list order '%s' has been requested" % (items_ordering_id))
                return json.dumps(items_ordering['basic'])
            else:
                return json.dumps({'error': "The interface requested something anvi'o doesn't know about. Item orderings\
                                             can only be in the form of 'newick' or 'basic'. But the interface requested\
                                             a '%s'. We are all confused here :/" % items_ordering_id})

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


    def charts(self, split_name):
        data = {'layers': [],
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

        data['index'], data['total'], data['previous_contig_name'], data['next_contig_name'] = self.get_index_total_previous_and_next_items(split_name)

        layers = sorted(self.interactive.p_meta['samples'])

        coverage_values_dict = self.interactive.split_coverage_values.get(split_name)
        data['coverage'] = [coverage_values_dict[layer].tolist() for layer in layers]

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

    def get_index_total_previous_and_next_items(self, item_name):
        previous_item_name = None
        next_item_name = None
        index = None
        total = None

        index_of_item = self.interactive.displayed_item_names_ordered.index(item_name)
        if index_of_item:
            previous_item_name = self.interactive.displayed_item_names_ordered[index_of_item - 1]
        if (index_of_item + 1) < len(self.interactive.displayed_item_names_ordered):
            next_item_name = self.interactive.displayed_item_names_ordered[index_of_item + 1]

        index = index_of_item + 1
        total = len(self.interactive.displayed_item_names_ordered)

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
        collections = dbops.TablesForCollections(db_path)
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
            return json.dumps({'error': 'Something failed. This is what we know: %s' % e})


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


    def get_protein_clusters_summary(self):
        protein_cluster_ids = json.loads(request.forms.get('split_names'))
        bin_name = json.loads(request.forms.get('bin_name'))

        summary = self.interactive.get_summary_for_PCs_list(protein_cluster_ids)

        run.info_single('PC info has been requested for %d items in %s' % (len(protein_cluster_ids), bin_name))

        return json.dumps({'functions': summary['functions'],
                           'num_PCs': summary['num_PCs'],
                           'genomes_contributing': summary['genomes_contributing'],
                           'num_gene_calls': summary['num_gene_calls']})


    def get_AA_sequences_for_PC(self, pc_name):
        data = {}

        if pc_name not in self.interactive.protein_clusters:
            return data

        if not self.interactive.genomes_storage_is_available:
            return data

        # add the list of gene caller ids associated with this protein cluster into `data`:
        for genome_name in self.interactive.protein_clusters[pc_name]:
            for gene_callers_id in self.interactive.protein_clusters[pc_name][genome_name]:
                data['%s_%s' % (genome_name, str(gene_callers_id))] = self.interactive.genomes_storage.get_gene_sequence(genome_name, gene_callers_id)

        return json.dumps(data)


    def inspect_pc(self, pc_name):
        data = {'pc_name': pc_name,
                'genomes': [],
                'index': None,
                'gene_caller_ids': [],
                'gene_caller_ids_in_genomes': {},
                'aa_sequences_in_pc': {},
                'previous_pc_name': None,
                'next_pc_name': None,
                'index': None,
                'total': None
                }

        if pc_name not in self.interactive.protein_clusters:
            return data

        if not self.interactive.genomes_storage_is_available:
            return data

        AA_sequences = self.interactive.get_AA_sequences_for_PCs(pc_names=set([pc_name]))

        # add the list of gene caller ids associated with this protein cluster into `data`:
        for genome_name in self.interactive.protein_clusters[pc_name]:
            data['aa_sequences_in_pc'][genome_name] = {}
            for gene_callers_id in self.interactive.protein_clusters[pc_name][genome_name]:
                data['gene_caller_ids'].append((gene_callers_id, genome_name), )
                data['aa_sequences_in_pc'][genome_name][gene_callers_id] = AA_sequences[pc_name][genome_name][gene_callers_id]

        # the dict that explains the distribution of genes in genomes:
        data['gene_caller_ids_in_genomes'] = self.interactive.protein_clusters[pc_name]

        # add the list of genomes into data:
        data['genomes'] = sorted(data['gene_caller_ids_in_genomes'].keys())

        # get some contextual stuff
        data['index'], data['total'], data['previous_pc_name'], data['next_pc_name'] = self.get_index_total_previous_and_next_items(pc_name)

        return json.dumps(data)


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
