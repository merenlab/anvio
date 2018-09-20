# -*- coding: utf-8
# pylint: disable=line-too-long
"""A library to build vignettes and networks of anvi'o programs"""

import os
import glob
import json

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
               'raw-bam-file': {'name': 'RAW BAM FILE', 'type': 'BAM', 'internal': False},
               'locus-fasta': {'name': 'LOCUS', 'type': 'FASTA', 'internal': False},
               'structure-db': {'name': 'STRUCTURE', 'type': 'DB', 'internal': True},
               'single-profile-db': {'name': 'SINGLE PROFILE', 'type': 'DB', 'internal': True},
               'profile-db': {'name': 'PROFILE', 'type': 'DB', 'internal': True},
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
               'misc-data-layer-orders': {'name': 'LAYER ORDERS DATA', 'type': 'CONCEPT', 'internal': False},
               'misc-data-layer-orders-txt': {'name': 'LAYER ORDERS DATA', 'type': 'TXT', 'internal': False},
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
               'external-genomes': {'name': 'EXTERNAL GENOMES', 'type': 'TXT', 'internal': False},
               'internal-genomes': {'name': 'INTERNAL GENOMES', 'type': 'TXT', 'internal': False},
               'ani-txt': {'name': 'ANI ESTIMATES', 'type': 'TXT', 'internal': False},
               'ani': {'name': 'ANI ESTIMATES', 'type': 'CONCEPT', 'internal': True},
               'variability-profile': {'name': 'VARIABILITY PROFILE', 'type': 'CONCEPT', 'internal': False},
               'codon-frequencies-txt': {'name': 'CODON FREQUENCIES', 'type': 'TXT', 'internal': False},
               'aa-frequencies-txt': {'name': 'AA FREQUENCIES', 'type': 'TXT', 'internal': False},
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
            self.main_programs = M('anvi-interactive')
            self.scripts = S('anvi-script-gen-programs-vignette')

            self.all_programs = sorted(list(set(self.main_programs + self.scripts)))
        except:
            raise ConfigError("Something is wrong. Either your installation or anvi'o setup on this computer is missing some of\
                               the fundamental programs, or your configuration is broken :/")

        if not len(self.main_programs) or not len(self.scripts):
            raise ConfigError("Somethings fishy is happening. This script is unable to find things that want to be found :(")

        self.run.info("Main anvi'o programs found", len(self.main_programs))
        self.run.info("Anvi'o ad hoc scripts found", len(self.scripts))

        if self.program_names_to_focus:
            self.program_names_to_focus = [p.strip() for p in self.program_names_to_focus.split(',')]
            run.info("Program names to focus", len(self.program_names_to_focus))

            self.all_programs = [p for p in self.all_programs if os.path.basename(p) in self.program_names_to_focus]

            if not len(self.all_programs):
                raise ConfigError("No anvi'o programs left to analyze after changing the focus to your list of program names.\
                                   Probably there is a typo or something :/")


class Program:
    def __init__(self, name, program_data):
        self.name = name

        self.provides = []
        self.requires = []

        factory = {'requires': self.requires,
                   'provides': self.provides}

        for x in factory:
            [factory[x].append(Item(item_name)) for item_name in program_data[x]]


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
        programs_dict = {}
        num_all_programs = len(self.all_programs)

        self.progress.new('Bleep bloop')
        for i in range(0, num_all_programs):
            program_path = self.all_programs[i]
            program_name = os.path.basename(program_path)

            self.progress.update('%s (%d of %d)' % (program_name, i+1, num_all_programs))

            requires = get_meta_information_from_file(program_path, '__requires__')
            provides = get_meta_information_from_file(program_path, '__provides__')

            if requires or provides:
                programs_dict[program_name] = {'requires': requires,
                                               'provides': provides}

        progress.end()

        if len(programs_dict):
            self.run.info_single("Of %d programs found, %d did contain provides and requires \
                                  statements." % (len(self.all_programs), len(programs_dict)),
                                  nl_after=1, nl_before=1)
            if anvio.DEBUG:
                absentees = ', '.join(list(set([os.path.basename(p) for p in self.all_programs]) - set(programs_dict.keys())))
                self.run.info_single("Here is a list of programs that do not contain any information\
                                      about themselves: %s" % (absentees), nl_after=1, nl_before=1, mc="red")
        else:
            raise ConfigError("None of the %d anvi'o programs found contained any provides or\
                               requires statements :/" % len(self.all_programs))

        self.programs = [Program(p, programs_dict[p]) for p in programs_dict]

        self.report_network()


    def report_network(self):
        item_names_seen = set([])
        items_seen = Counter({})
        all_items = []
        for program in self.programs:
            for item in program.provides + program.requires:
                items_seen[item.id] += 1
                if not item.id in item_names_seen:
                    all_items.append(item)
                    item_names_seen.add(item.id)

        programs_seen = Counter({})
        for item in all_items:
            for program in self.programs:
                for program_item in program.provides + program.requires:
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
                for item_provided in program.provides:
                    if item_provided.id == item.id:
                        network_dict["links"].append({"source": node_indices[program.name], "target": node_indices[item.id]})
                for item_needed in program.requires:
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
        d = {}

        log_file = filesnpaths.get_temp_file_path()
        num_all_programs = len(self.all_programs)
        for i in range(0, num_all_programs):
            program_path = self.all_programs[i]
            program_name = os.path.basename(program_path)

            if program_name in self.programs_to_skip:
                run.warning("Someone doesn't want %s to be in the output :/ Fine. Skipping." % (program_name))

            progress.new('Bleep bloop')
            progress.update('%s (%d of %d)' % (program_name, i+1, num_all_programs))

            output = utils.run_command_STDIN('%s --help' % (program_path), log_file, '').split('\n')

            if anvio.DEBUG:
                    usage, description, params, output = parse_help_output(output)
            else:
                try:
                    usage, description, params, output = parse_help_output(output)
                except Exception as e:
                    progress.end()
                    run.warning("The program '%s' does not seem to have the expected help menu output. Skipping to the next.\
                                 For the curious, this was the error message: '%s'" % (program_name, str(e).strip()))
                    continue

            d[program_name] = {'usage': usage,
                               'description': description,
                               'params': params,
                               'tags': get_meta_information_from_file(program_path, '__tags__'),
                               'resources': get_meta_information_from_file(program_path, '__resources__')}

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
