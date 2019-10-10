# -*- coding: utf-8
# pylint: disable=line-too-long
"""A library to build vignettes and networks of anvi'o programs"""

import os
import sys
import glob
import json
import importlib

from collections import Counter

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.summaryhtml import SummaryHTMLOutput
from anvio.errors import ConfigError

__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


G = lambda d: [p for p in glob.glob(os.path.join(d, 'anvi-*')) if utils.is_program_exists(p, dont_raise=True)]
M = lambda m: [x for x in G(os.path.dirname(utils.is_program_exists(m)))]
S = lambda s: [x for x in G(os.path.dirname(utils.is_program_exists(s)))]
J = lambda x: '\n'.join(x) if x else ''

# this dictionary describes all anvi'o items that are referred from 'requires' and
# 'provudes' statements written in anvi'o programs
ANVIO_ITEMS = {'pan-db': {'name': 'PAN', 'type': 'DB', 'internal': True},
               'contigs-db': {'name': 'CONTIGS', 'type': 'DB', 'internal': True},
               'contigs-fasta': {'name': 'CONTIGS', 'type': 'FASTA', 'internal': False},
               'concatenated-gene-alignment-fasta': {'name': 'CONCATENATED GENE ALIGNMENT', 'type': 'FASTA', 'internal': False},
               'short-reads-fasta': {'name': 'SHORT READS', 'type': 'FASTA', 'internal': False},
               'genes-fasta': {'name': 'GENES', 'type': 'FASTA', 'internal': False},
               'bam-file': {'name': 'BAM FILE', 'type': 'BAM', 'internal': False},
               'protein-structure': {'name': 'PDB FILE', 'type': 'TXT', 'internal': False},
               'raw-bam-file': {'name': 'RAW BAM FILE', 'type': 'BAM', 'internal': False},
               'locus-fasta': {'name': 'LOCUS', 'type': 'FASTA', 'internal': False},
               'structure-db': {'name': 'STRUCTURE', 'type': 'DB', 'internal': True},
               'single-profile-db': {'name': 'SINGLE PROFILE', 'type': 'DB', 'internal': True},
               'profile-db': {'name': 'PROFILE', 'type': 'DB', 'internal': True},
               'genes-db': {'name': 'GENES', 'type': 'DB', 'internal': True},
               'genomes-storage-db': {'name': 'GENOMES STORAGE', 'type': 'DB', 'internal': True},
               'contigs-stats': {'name': 'CONTIGS STATS', 'type': 'STATS', 'internal': False},
               'svg': {'name': 'SVG', 'type': 'SVG', 'internal': False},
               'bin': {'name': 'BIN', 'type': 'BIN', 'internal': True},
               'collection': {'name': 'COLLECTION', 'type': 'COLLECTION', 'internal': True},
               'collection-txt': {'name': 'COLLECTION', 'type': 'TXT', 'internal': False},
               'hmm-source': {'name': 'HMM SOURCE', 'type': 'HMM', 'internal': False},
               'hmm-profile': {'name': 'HMM PROFILE', 'type': 'CONCEPT', 'internal': True},
               'cogs-data': {'name': 'COGs DATA', 'type': 'DATA', 'internal': True},
               'pfams-data': {'name': 'PFAMs DATA', 'type': 'DATA', 'internal': True},
               'misc-data-items-txt': {'name': 'ITEMS DATA', 'type': 'TXT', 'internal': False},
               'misc-data-items': {'name': 'ITEMS DATA', 'type': 'CONCEPT', 'internal': False},
               'misc-data-layers-txt': {'name': 'LAYERS DATA', 'type': 'TXT', 'internal': False},
               'misc-data-layers': {'name': 'LAYERS DATA', 'type': 'CONCEPT', 'internal': False},
               'misc-data-layers-category': {'name': 'LAYERS DATA CATEGORY', 'type': 'CONCEPT', 'internal': True},
               'genome-similarity': {'name': 'GENOME SIMILARITY', 'type': 'CONCEPT', 'internal': False},
               'misc-data-layer-orders': {'name': 'LAYER ORDERS DATA', 'type': 'CONCEPT', 'internal': False},
               'misc-data-item-orders': {'name': 'ITEM ORDERS DATA', 'type': 'CONCEPT', 'internal': False},
               'misc-data-layer-orders-txt': {'name': 'LAYER ORDERS DATA', 'type': 'TXT', 'internal': False},
               'misc-data-item-orders-txt': {'name': 'LAYER ORDERS DATA', 'type': 'TXT', 'internal': False},
               'dendrogram': {'name': 'DENDROGRAM', 'type': 'NEWICK', 'internal': False},
               'metapangenome': {'name': 'METAPANGENOME', 'type': 'CONCEPT', 'internal': True},
               'oligotypes': {'name': 'OLIGOTYPES', 'type': 'CONCEPT', 'internal': False},
               'linkmers-txt': {'name': 'LINKMERS', 'type': 'TXT', 'internal': False},
               'phylogeny': {'name': 'PHYLOGENY', 'type': 'NEWICK', 'internal': False},
               'gene-calls-txt': {'name': 'GENE CALLS', 'type': 'TXT', 'internal': False},
               'functions': {'name': 'GENE FUNCTIONS', 'type': 'CONCEPT', 'internal': True},
               'functions-txt': {'name': 'GENE FUNCTIONS', 'type': 'TXT', 'internal': False},
               'functional-enrichment-txt': {'name': 'ENRICHED FUNCTIONS', 'type': 'TXT', 'internal': False},
               'interactive': {'name': 'INTERACTIVE DISPLAY', 'type': 'DISPLAY', 'internal': True},
               'view-data': {'name': 'VIEW DATA', 'type': 'TXT', 'internal': False},
               'layer-taxonomy': {'name': 'LAYER TAXONOMY', 'type': 'CONCEPT', 'internal': True},
               'layer-taxonomy-txt': {'name': 'LAYER TAXONOMY', 'type': 'TXT', 'internal': False},
               'gene-taxonomy': {'name': 'GENE TAXONOMY', 'type': 'CONCEPT', 'internal': True},
               'gene-taxonomy-txt': {'name': 'GENE TAXONOMY', 'type': 'TXT', 'internal': False},
               'genome-taxonomy': {'name': 'GENOME TAXONOMY', 'type': 'CONCEPT', 'internal': True},
               'genome-taxonomy-txt': {'name': 'GENOME TAXONOMY', 'type': 'TXT', 'internal': False},
               'scgs-taxonomy-db': {'name': 'SCG TAXONOMY DB', 'type': 'CONCEPT', 'internal': True},
               'scgs-taxonomy': {'name': 'SCG TAXONOMY', 'type': 'CONCEPT', 'internal': True},
               'external-genomes': {'name': 'EXTERNAL GENOMES', 'type': 'TXT', 'internal': False},
               'internal-genomes': {'name': 'INTERNAL GENOMES', 'type': 'TXT', 'internal': False},
               'coverages-txt': {'name': 'COVERAGES', 'type': 'TXT', 'internal': False},
               'genome-distance-txt': {'name': 'DISTANCE ESTIMATES', 'type': 'TXT', 'internal': False},
               'genome-distance': {'name': 'DISTANCE ESTIMATES', 'type': 'CONCEPT', 'internal': True},
               'variability-profile': {'name': 'VARIABILITY PROFILE', 'type': 'CONCEPT', 'internal': False},
               'codon-frequencies-txt': {'name': 'CODON FREQUENCIES', 'type': 'TXT', 'internal': False},
               'aa-frequencies-txt': {'name': 'AA FREQUENCIES', 'type': 'TXT', 'internal': False},
               'fixation-index-matrix': {'name': 'FIXATION INDEX MATRIX', 'type': 'TXT', 'internal': False},
               'summary': {'name': 'STATIC SUMMARY', 'type': 'SUMMARY', 'internal': False},
               'split-bins': {'name': 'SPLIT BINS', 'type': 'CONCEPT', 'internal': False},
               'state': {'name': 'INTERACTIVE STATE', 'type': 'CONCEPT', 'internal': True},
               'state-json': {'name': 'INTERACTIVE STATE', 'type': 'JSON', 'internal': False}}

ANVIO_CONCEPTS = {'functions': {'goes_in': ['contigs_db', 'genomes-storage-db'],
                               'used_by': ['anvi-search-functions']}
                 }

run = terminal.Run()
progress = terminal.Progress()


def get_until_blank(output):
    section = []
    while 1:
        if output[0] == '':
            break
        else:
            section.append(output.pop(0))

    return section


def get_meta_information_from_file(file_path, meta_tag):
    all_lines = [l.strip() for l in open(file_path, 'rU').readlines()]

    meta_tag_content = ''

    while 1:
        if not len(all_lines):
            return []

        if not all_lines[0].startswith(meta_tag):
            all_lines.pop(0)
        else:
            break

    meta_tag_content = all_lines.pop(0)

    while 1:
        line = all_lines.pop(0)
        if line == '' or line.startswith('__'):
            break
        else:
            meta_tag_content += line

    if meta_tag_content:
        return eval(meta_tag_content.split('=')[1].strip())
    else:
        return []


def get_param_set(output):
    if output[0] in ['optional arguments:', 'positional arguments:']:
        section = output.pop(0)
        desc = ''
        _params = J([p for p in get_until_blank(output) if not p.startswith('  -h, --help')])
    else:
        section = output.pop(0)
        if output[0].startswith('  -'):
            # no description, goes into params immediately (someone did a crappy job)
            desc = ''
        else:
            desc = get_until_blank(output)
            output.pop(0)

        _params = J(get_until_blank(output))

    return section, desc, _params


def skip_until_usage(output):
    while 1:
        if not len(output):
            return

        if output[0].startswith('usage:'):
            return

        output.pop(0)


def parse_help_output(output):
    skip_until_usage(output)

    if not len(output):
        raise ConfigError("This is not the help menu output we are looking for.")

    if not output[0].startswith('usage:'):
        raise ConfigError("This output does not seem to have the proper usage statement.")

    usage = J([l[7:] for l in get_until_blank(output)])

    if output.pop(0) != '':
        raise ConfigError("This output is missing the description start marker.")

    description = J(get_until_blank(output))

    params = {}
    while 1:
        if output.pop(0) != '':
            raise ConfigError("The params section does not seem to be where this script expects to find it.")

        if not len(output):
            break

        section, desc, _params = get_param_set(output)
        if _params == '':
            pass
        else:
            params[section] = {'description': J(desc),
                               'params': _params}

    return usage, description, params, output


class AnvioPrograms:
    def __init__(self, args, r=terminal.Run(), p=terminal.Progress()):
        self.args = args
        self.run = r
        self.progress = p

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.program_names_to_focus = A("program_names_to_focus")

        try:
            self.main_program_filepaths = M('anvi-interactive')
            self.script_filepaths = S('anvi-script-gen-programs-vignette')

            self.all_program_filepaths = sorted(list(set(self.main_program_filepaths + self.script_filepaths)))
        except:
            raise ConfigError("Something is wrong. Either your installation or anvi'o setup on this computer is missing some of\
                               the fundamental programs, or your configuration is broken :/")

        if not len(self.main_program_filepaths) or not len(self.script_filepaths):
            raise ConfigError("Somethings fishy is happening. This script is unable to find things that want to be found :(")

        self.run.info("Main anvi'o programs found", len(self.main_program_filepaths))
        self.run.info("Anvi'o ad hoc scripts found", len(self.script_filepaths))

        if self.program_names_to_focus:
            self.program_names_to_focus = [p.strip() for p in self.program_names_to_focus.split(',')]
            run.info("Program names to focus", len(self.program_names_to_focus))

            self.all_program_filepaths = [p for p in self.all_program_filepaths if os.path.basename(p) in self.program_names_to_focus]

            if not len(self.all_program_filepaths):
                raise ConfigError("No anvi'o programs left to analyze after changing the focus to your list of program names.\
                                   Probably there is a typo or something :/")


    def create_program_classes(self, okay_if_no_meta=False, quiet=False):
        programs_dict = {}
        num_all_programs = len(self.all_program_filepaths)

        meta_count = 0
        self.programs = []
        self.progress.new('Characterizing program', progress_total_items=num_all_programs)

        for program_filepath in self.all_program_filepaths:
            self.progress.update(os.path.basename(program_filepath))
            self.progress.increment()

            program = Program(program_filepath)

            if program.meta_info['provides']['value'] or program.meta_info['requires']['value']:
                meta_count += 1

            if not (program.meta_info['provides']['value'] or program.meta_info['requires']['value']) and not okay_if_no_meta:
                pass
            else:
                self.programs.append(program)

        self.progress.end()

        if not meta_count and not okay_if_no_meta:
            raise ConfigError("None of the %d anvi'o programs found contained any provides or\
                               requires statements :/" % len(self.all_program_filepaths))

        if not quiet:
            self.run.info_single("Of %d programs found, %d did contain provides and/or requires \
                                  statements." % (len(self.all_program_filepaths), meta_count),
                                  nl_after=1, nl_before=1)
        if anvio.DEBUG:
            absentees = ', '.join(list(set([os.path.basename(p) for p in self.all_program_filepaths]) - set([p.name for p in self.programs])))
            self.run.info_single("Here is a list of programs that do not contain any information\
                                  about themselves: %s" % (absentees), nl_after=1, nl_before=1, mc="red")



class Program:
    def __init__(self, program_path):
        self.program_path = program_path
        self.name = os.path.basename(program_path)

        self.meta_info = {
            'requires': {
                'object_name': '__requires__',
                'null_object': []
            },
            'provides': {
                'object_name': '__provides__',
                'null_object': []
            },
            'tags': {
                'object_name': '__tags__',
                'null_object': []
            },
            'resources': {
                'object_name': '__resources__',
                'null_object': []
            },
            'description': {
                'object_name': '__description__',
                'null_object': ''
            },
            'status': {
                'object_name': '__status__',
                'null_object': ''
            },
        }

        self.module = self.load_as_module(self.program_path)
        self.get_meta_info()


    def get_meta_info(self):
        for info_type in self.meta_info.keys():
            try:
                info = getattr(self.module, self.meta_info[info_type]['object_name'])
            except AttributeError:
                info = self.meta_info[info_type]['null_object']

            if info_type == 'requires' or info_type == 'provides':
                # these info_types have their items cast as Item types
                info = [Item(item_name) for item_name in info]

            if type(info) == str:
                info = info.replace('\n', ' ')

            self.meta_info[info_type]['value'] = info


    def load_as_module(self, path):
        """
        Importing the program as a module has the advantage of grabbing the meta info as python
        objects directly instead of parsing the file lines as a string object. if is not a Python
        file, self.module is None.

        Taken from stackoverflow user Ciro Santilli:
        https://stackoverflow.com/questions/2601047/import-a-python-module-without-the-py-extension/56090741#56090741
        """
        try:
            module_name = os.path.basename(path).replace('-', '_')
            spec = importlib.util.spec_from_loader(
                module_name,
                importlib.machinery.SourceFileLoader(module_name, path)
            )
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            sys.modules[module_name] = module
            return module
        except:
            return None


class Item:
    def __init__(self, item_id, internal=True, optional=True, single=True):
        if item_id not in ANVIO_ITEMS:
            raise ConfigError("Ehem. Anvi'o does not know about item '%s'. There are two was this could happen:\
                               one, you made a type (easy to fix), two, you just added a new program into anvi'o\
                               but have not yet updated `anvio/programs.py`." % item_id)

        item = ANVIO_ITEMS[item_id]
        self.id = item_id
        self.name = item['name']
        self.type = item['type']
        self.internal = item['internal']

        # attributes set by the context master
        self.single = single
        self.optional = optional


class ProgramsNetwork(AnvioPrograms):
    def __init__(self, args, r=terminal.Run(), p=terminal.Progress()):
        self.args = args
        self.run = r
        self.progress = p

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.output_file_path = A("output_file") or 'NETWORK.json'

        filesnpaths.is_output_file_writable(self.output_file_path)

        AnvioPrograms.__init__(self, args, r=self.run, p=self.progress)


    def generate(self):
        self.create_program_classes()
        self.report_network()


    def report_network(self):
        item_names_seen = set([])
        items_seen = Counter({})
        all_items = []
        for program in self.programs:
            for item in program.meta_info['provides']['value'] + program.meta_info['requires']['value']:
                items_seen[item.id] += 1
                if not item.id in item_names_seen:
                    all_items.append(item)
                    item_names_seen.add(item.id)

        programs_seen = Counter({})
        for item in all_items:
            for program in self.programs:
                for program_item in program.meta_info['provides']['value'] + program.meta_info['requires']['value']:
                    if item.name == program_item.name:
                        programs_seen[program.name] += 1

        network_dict = {"graph": [], "nodes": [], "links": [], "directed": False, "multigraph": False}

        node_indices = {}

        index = 0
        types_seen = set(["PROGRAM"])
        for item in all_items:
            types_seen.add(item.type)
            network_dict["nodes"].append({"size": items_seen[item.id],
                                          "score": 0.5 if item.internal else 1,
                                          "color": '#00AA00' if item.internal else "#AA0000",
                                          "id": item.id,
                                          "name": item.name,
                                          "internal": True if item.internal else False,
                                          "type": item.type})
            node_indices[item.id] = index
            index += 1

        for program in self.programs:
            network_dict["nodes"].append({"size": programs_seen[program.name],
                                          "score": 0.1,
                                          "color": "#AAAA00",
                                          "id": program.name,
                                          "name": program.name,
                                          "type": "PROGRAM"})
            node_indices[program.name] = index
            index += 1

        for item in all_items:
            for program in self.programs:
                for item_provided in program.meta_info['provides']['value']:
                    if item_provided.id == item.id:
                        network_dict["links"].append({"source": node_indices[program.name], "target": node_indices[item.id]})
                for item_needed in program.meta_info['requires']['value']:
                    if item_needed.id == item.id:
                        network_dict["links"].append({"target": node_indices[program.name], "source": node_indices[item.id]})

        open(self.output_file_path, 'w').write(json.dumps(network_dict, indent=2))

        self.run.info('JSON description of network', self.output_file_path)
        self.run.info('Item types seen', ', '.join(sorted(list(types_seen))))


class ProgramsVignette(AnvioPrograms):
    def __init__(self, args, r=terminal.Run(), p=terminal.Progress()):
        self.args = args
        self.run = r
        self.progress = p

        self.programs_to_skip = ['anvi-script-gen-programs-vignette']

        AnvioPrograms.__init__(self, args, r=self.run, p=self.progress)

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.output_file_path = A("output_file")


    def generate(self):
        self.create_program_classes(okay_if_no_meta = True, quiet = True)

        d = {}
        log_file = filesnpaths.get_temp_file_path()
        for i, program in enumerate(self.programs):
            if program.name in self.programs_to_skip:
                run.warning("Someone doesn't want %s to be in the output :/ Fine. Skipping." % (program.name))

            progress.new('Bleep bloop')
            progress.update('%s (%d of %d)' % (program.name, i+1, len(self.programs)))

            output = utils.run_command_STDIN('%s --help' % (program.program_path), log_file, '').split('\n')

            if anvio.DEBUG:
                    usage, description, params, output = parse_help_output(output)
            else:
                try:
                    usage, description, params, output = parse_help_output(output)
                except Exception as e:
                    progress.end()
                    run.warning("The program '%s' does not seem to have the expected help menu output. Skipping to the next.\
                                 For the curious, this was the error message: '%s'" % (program.name, str(e).strip()))
                    continue

            d[program.name] = {'usage': usage,
                               'description': description,
                               'params': params,
                               'tags': program.meta_info['tags']['value'],
                               'resources': program.meta_info['resources']['value']}

            progress.end()

        os.remove(log_file)

        # generate output
        program_names = sorted([p for p in d if not p.startswith('anvi-script-')])
        script_names = sorted([p for p in d if p.startswith('anvi-script-')])
        vignette = {'vignette': d,
                    'program_names': program_names,
                    'script_names': script_names,
                    'all_names': program_names + script_names,
                    'meta': {'summary_type': 'vignette',
                             'version': '\n'.join(['|%s|%s|' % (t[0], t[1]) for t in anvio.get_version_tuples()]),
                             'date': utils.get_date()}}

        if anvio.DEBUG:
            run.warning(None, 'THE OUTPUT DICT')
            import json
            print(json.dumps(d, indent=2))

        open(self.output_file_path, 'w').write(SummaryHTMLOutput(vignette, r=run, p=progress).render())

        run.info('Output file', os.path.abspath(self.output_file_path))
