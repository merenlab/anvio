# -*- coding: utf-8
# pylint: disable=line-too-long
"""A library to help anvi'o desribe itself"""

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

from anvio.errors import ConfigError
from anvio.docs import ANVIO_ARTIFACTS
from anvio.summaryhtml import SummaryHTMLOutput


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
            self.all_program_names = [os.path.basename(p) for p in self.all_program_filepaths]
        except:
            raise ConfigError("Something is wrong. Either your installation or anvi'o setup on this computer is missing some of "
                              "the fundamental programs, or your configuration is broken :/")

        if not len(self.main_program_filepaths) or not len(self.script_filepaths):
            raise ConfigError("Somethings fishy is happening. This script is unable to find things that want to be found :(")

        self.run.info("Main anvi'o programs found", len(self.main_program_filepaths))
        self.run.info("Anvi'o ad hoc scripts found", len(self.script_filepaths))

        if self.program_names_to_focus:
            self.program_names_to_focus = [p.strip() for p in self.program_names_to_focus.split(',')]
            run.info("Program names to focus", len(self.program_names_to_focus))

            self.all_program_filepaths = [p for p in self.all_program_filepaths if os.path.basename(p) in self.program_names_to_focus]

            if not len(self.all_program_filepaths):
                raise ConfigError("No anvi'o programs left to analyze after changing the focus to your list of program names. "
                                  "Probably there is a typo or something :/")


    def init_programs(self, okay_if_no_meta=False, quiet=False):
        """Initializes the `self.programs` dictionary."""

        num_all_programs = len(self.all_program_filepaths)

        self.programs = {}
        self.progress.new('Characterizing program', progress_total_items=num_all_programs)

        programs_with_usage_info = set([])
        programs_without_usage_info = set([])
        programs_with_provides_requires_info = set([])
        programs_without_provides_requires_info = set([])

        for program_filepath in self.all_program_filepaths:
            self.progress.update(os.path.basename(program_filepath), increment=True)

            program = Program(program_filepath, r=self.run, p=self.progress)

            program_usage_information_path = os.path.join(self.docs_path, 'programs/%s.md' % (program.name))

            if program.meta_info['provides']['value'] or program.meta_info['requires']['value']:
                programs_with_provides_requires_info.add(program.name)
            else:
                programs_without_provides_requires_info.add(program.name)

            # learn about the usage statement of the program
            if os.path.exists(program_usage_information_path):
                program.usage = self.read_anvio_markdown(program_usage_information_path)
                programs_with_usage_info.add(program.name)
            else:
                programs_without_usage_info.add(program.name)

            if not (program.meta_info['provides']['value'] or program.meta_info['requires']['value']) and not okay_if_no_meta:
                try:
                    programs_with_usage_info.remove(program.name)
                    programs_without_usage_info.remove(program.name)
                except:
                    pass
            else:
                self.programs[program.name] = program

        self.progress.end()

        # report missing provides/requires information
        self.run.info_single("Of %d programs found, %d did not contain PROVIDES AND/OR REQUIRES "
                             "statements :/ This may be normal for some programs, but here is the "
                             "complete list of those that are missing __provides__ and __requires__ "
                             "tags in their code in case you see something you can complete: '%s'." % \
                                        (len(self.all_program_filepaths),
                                         len(programs_without_provides_requires_info),
                                         ', '.join(programs_without_provides_requires_info)),
                             nl_after=1, nl_before=1)

        # report missing provides/requires information
        self.run.info_single("Of %d programs found, %d did not have any USAGE INFORMATION. You can "
                             "help by adding usage information for programs by creating markdown "
                             "formatted files under the directory '%s'. Please see examples in anvi'o "
                             "codebase: https://github.com/merenlab/anvio/tree/master/anvio/docs. "
                             "Here is a complete list of programs that are missing usage statements: %s " % \
                                        (len(self.all_program_filepaths),
                                         len(programs_without_provides_requires_info),
                                         self.docs_path, 
                                         ', '.join(programs_without_provides_requires_info)),
                             nl_after=1, nl_before=1)


class Program:
    def __init__(self, program_path, r=terminal.Run(), p=terminal.Progress()):
        self.run = r
        self.progress = p

        self.program_path = program_path
        self.name = os.path.basename(program_path)
        self.usage = None

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
                # these info_types have their items cast as Artifact types
                info = [Artifact(artifact_name) for artifact_name in info]

            if type(info) == str:
                if 'read_as_is' in self.meta_info[info_type] and self.meta_info[info_type]['read_as_is']:
                    info = info
                else:
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


    def __str__(self):
        self.run.warning(None, header='%s' % self.name, lc='green')
        for info_type in self.meta_info:
            self.run.info(info_type, self.meta_info[info_type]['value'])

        return ''


    def __repr__(self):
        return "PROG::%s" % self.name


class Artifact:
    """A class to describe an anvi'o artifact"""

    def __init__(self, artifact_id, internal=True, optional=True, single=True):
        if artifact_id not in ANVIO_ARTIFACTS:
            raise ConfigError("Ehem. Anvi'o does not know about artifact '%s'. There are two was this could happen: "
                              "one, you made a type (easy to fix), two, you just added a new program into anvi'o "
                              "but have not yet updated `anvio/programs.py`." % artifact_id)

        artifact = ANVIO_ARTIFACTS[artifact_id]
        self.id = artifact_id
        self.name = artifact['name']
        self.type = artifact['type']
        self.internal = artifact['internal']

        # attributes set by the context master
        self.single = single
        self.optional = optional


    def __repr__(self):
        return "ARTIFACT::%s" % self.id


class AnvioArtifacts:
    """Information on anvi'o artifacts"""

    def __init__(self, args, r=terminal.Run(), p=terminal.Progress()):
        self.args = args
        self.run = r
        self.progress = p

        self.artifacts_info = {}

        if not hasattr(self, 'programs'):
            raise ConfigError("AnvioArtifacts class is upset. You need to treat this class as a base class, and initialize "
                              "it from within another class that has already initialized AnvioPrograms class. If this is "
                              "confusing, take a look at the AnvioDocs class.")

        if not len(self.programs):
            raise ConfigError("AnvioArtifacts is initialized with a blank `self.programs` variable. HOW CUTE.")


    def init_artifacts(self):
        """Generate `required_by` and `provided_by` statements.

        Returns
        =======
        artifacts_info, dict:
            Running this function will fill in the dictionary `self.artifacts_info`
        """

        artifacts_with_descriptions = set([])
        artifacts_without_descriptions = set([])
    
        for artifact in ANVIO_ARTIFACTS:
            self.artifacts_info[artifact] = {'required_by': [], 'provided_by': [], 'description': None}

            # learn about the description of the artifact
            artifact_description_path = os.path.join(self.docs_path, 'artifacts/%s.md' % (artifact))
            if os.path.exists(artifact_description_path):
                self.artifacts_info[artifact]['description'] = self.read_anvio_markdown(artifact_description_path)
                artifacts_with_descriptions.add(artifact)
            else:
                artifacts_without_descriptions.add(artifact)

            # learn about what provides or requires this artifact
            for program in self.programs.values():
                if artifact in [a.id for a in program.meta_info['requires']['value']]:
                    self.artifacts_info[artifact]['required_by'].append(program.name)

                if artifact in [a.id for a in program.meta_info['provides']['value']]:
                    self.artifacts_info[artifact]['provided_by'].append(program.name)

        if len(artifacts_without_descriptions):
            self.run.info_single("Of %d artifacts found, %d did not contain any DESCRIPTION. If you would like to "
                                 "see examples and add new descriptions, please see the directory '%s'. Here is the "
                                 "full list of artifacts that are not yet explained: %s." \
                                        % (len(ANVIO_ARTIFACTS),
                                           len(artifacts_without_descriptions),
                                           self.docs_path,
                                           ', '.join(artifacts_without_descriptions)), nl_after=1, nl_before=1)


class AnvioDocs(AnvioPrograms, AnvioArtifacts):
    """Generate a docs output.

       The purpose of this class is to generate a static HTML output with
       interlinked files that serve as the primary documentation for anvi'o
       programs, input files they expect, and output files the generate.

       The default client of this class is `anvi-script-gen-help-docs`.
    """

    def __init__(self, args, r=terminal.Run(), p=terminal.Progress()):
        self.args = args
        self.run = r
        self.progress = p

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.output_directory_path = A("output_dir") or 'ANVIO-HELP'

        self.docs_path = os.path.join(os.path.dirname(anvio.__file__), 'docs')
        if not os.path.exists(self.docs_path):
            raise ConfigError("The anvi'o docs path is not where it should be :/ Something funny is going on.")

        filesnpaths.gen_output_directory(self.output_directory_path, delete_if_exists=True, dont_warn=True)

        self.artifacts_output_dir = filesnpaths.gen_output_directory(os.path.join(self.output_directory_path, 'artifacts'))
        self.programs_output_dir = filesnpaths.gen_output_directory(os.path.join(self.output_directory_path, 'programs'))

        self.base_url = "/software/anvio/help"
        self.anvio_markdown_variables_conversion_dict = {}

        AnvioPrograms.__init__(self, args, r=self.run, p=self.progress)
        self.init_programs()

        AnvioArtifacts.__init__(self, args, r=self.run, p=self.progress)
        self.init_artifacts()

        if not len(self.programs):
            raise ConfigError("AnvioDocs is asked ot process the usage statements of some programs, but the "
                              "`self.programs` dictionary seems to be empty :/")

        self.images_source_directory = os.path.join(os.path.dirname(anvio.__file__), 'docs/images/png')

        self.sanity_check()


    def sanity_check(self):
        """Quick sanity checks to ensure things are working"""

        if not os.path.exists(self.images_source_directory):
            raise ConfigError("AnvioDocs speaking: the images source directory does not seem to be "
                              "where it should have been :/")

        # make sure each artifact type has an icon
        A_PNG = lambda x: os.path.exists(os.path.join(self.images_source_directory, 'icons', ANVIO_ARTIFACTS[x]['type'] + '.png'))
        missing_images_for_artifact_types = [ANVIO_ARTIFACTS[artifact]['type'] for artifact in self.artifacts_info if not A_PNG(artifact)]
        if len(missing_images_for_artifact_types):
            raise ConfigError("Some artifacts do not have matching images. If you just added a new artifact type, you "
                              "also need to add a corresponding PNG icon for the type under the directory '%s'. See "
                              "examples in that directory, and if they are not enough, get in touch with a developer. "
                              "Regardless. These are the artifact types missing images: %s." \
                                                                % (os.path.join(self.images_source_directory, 'icons'),
                                                                   ', '.join(missing_images_for_artifact_types)))


    def generate(self):
        self.copy_images()

        self.generate_pages_for_artifacts()

        self.generate_pages_for_programs()

        self.generate_index_page()


    def copy_images(self):
        """Copies images from the codebase to the output directory"""

        utils.shutil.copytree(self.images_source_directory, os.path.join(self.output_directory_path, 'images'))


    def init_anvio_markdown_variables_conversion_dict(self):
        for program_name in self.all_program_names:
            self.anvio_markdown_variables_conversion_dict[program_name] = """<span class="artifact-n">[%s](%s/programs/%s)</span>""" % (program_name, self.base_url, program_name)

        for artifact_name in ANVIO_ARTIFACTS:
            self.anvio_markdown_variables_conversion_dict[artifact_name] = """<span class="artifact-n">[%s](%s/artifacts/%s)</span>""" % (artifact_name, self.base_url, artifact_name)


    def read_anvio_markdown(self, file_path):
        """Reads markdown descriptions filling in anvi'o variables"""

        filesnpaths.is_file_plain_text(file_path)

        if not len(self.anvio_markdown_variables_conversion_dict):
            self.init_anvio_markdown_variables_conversion_dict()

        markdown_content = open(file_path).read()

        try:
            markdown_content = markdown_content % self.anvio_markdown_variables_conversion_dict
        except KeyError as e:
            self.progress.end()
            raise ConfigError("One of the variables, %s, in '%s' is not yet described anywhere :/ If it is not a typo but "
                              "a new artifact, you can add it to the file `anvio/docs/__init__.py`. After which everything "
                              "should work. But please also remember to update provides / requires statements of programs "
                              "for everything to be linked together." % (e, file_path))
        except Exception as e:
            self.progress.end()
            raise ConfigError("Something went wrong while working with '%s' :/ This is what we know: '%s'." % (file_path, e))

        markdown_content = markdown_content.replace("""{{ codestart }}""", """<div class="codeblock" markdown="1">""")
        markdown_content = markdown_content.replace("""{{ codestop }}""", """</div>""")
        markdown_content = markdown_content.replace("""-""", """&#45;""")


        return markdown_content


    def get_program_requires_provides_dict(self, prefix="../../"):
        d = {}

        for program_name in self.programs:
            d[program_name] = {}

            program = self.programs[program_name]
            d[program_name]['requires'] = [(r.id, '%sartifacts/%s' % (prefix, r.id)) for r in program.meta_info['requires']['value']]
            d[program_name]['provides'] = [(r.id, '%sartifacts/%s' % (prefix, r.id)) for r in program.meta_info['provides']['value']]

        return d


    def generate_pages_for_artifacts(self):
        """Generates static pages for artifacts in the output directory"""

        for artifact in ANVIO_ARTIFACTS:
            d = {'artifact': ANVIO_ARTIFACTS[artifact],
                 'meta': {'summary_type': 'artifact',
                          'version': '\n'.join(['|%s|%s|' % (t[0], t[1]) for t in anvio.get_version_tuples()]),
                          'date': utils.get_date()}
                }

            d['artifact']['name'] = artifact
            d['artifact']['required_by'] = [(r, '../../programs/%s' % r) for r in self.artifacts_info[artifact]['required_by']]
            d['artifact']['provided_by'] = [(r, '../../programs/%s' % r) for r in self.artifacts_info[artifact]['provided_by']]
            d['artifact']['description'] = self.artifacts_info[artifact]['description']
            d['artifact']['icon'] = '../../images/icons/%s.png' % ANVIO_ARTIFACTS[artifact]['type']

            if anvio.DEBUG:
                run.warning(None, 'THE OUTPUT DICT')
                import json
                print(json.dumps(d, indent=2))

            artifact_output_dir = filesnpaths.gen_output_directory(os.path.join(self.artifacts_output_dir, artifact))
            output_file_path = os.path.join(artifact_output_dir, 'index.md')
            open(output_file_path, 'w').write(SummaryHTMLOutput(d, r=run, p=progress).render())


    def generate_pages_for_programs(self):
        """Generates static pages for programs in the output directory"""

        program_provides_requires_dict = self.get_program_requires_provides_dict()

        for program_name in self.programs:
            program = self.programs[program_name]
            d = {'program': {},
                 'meta': {'summary_type': 'program',
                          'version': '\n'.join(['|%s|%s|' % (t[0], t[1]) for t in anvio.get_version_tuples()]),
                          'date': utils.get_date()}
                }

            d['program']['name'] = program_name
            d['program']['usage'] = program.usage
            d['program']['description'] = program.meta_info['description']['value']
            d['program']['resources'] = program.meta_info['resources']['value']
            d['program']['requires'] = program_provides_requires_dict[program_name]['requires']
            d['program']['provides'] = program_provides_requires_dict[program_name]['provides']
            d['program']['icon'] = '../../images/icons/%s.png' % 'PROGRAM'

            if anvio.DEBUG:
                run.warning(None, 'THE OUTPUT DICT')
                import json
                print(json.dumps(d, indent=2))

            program_output_dir = filesnpaths.gen_output_directory(os.path.join(self.programs_output_dir, program_name))
            output_file_path = os.path.join(program_output_dir, 'index.md')
            open(output_file_path, 'w').write(SummaryHTMLOutput(d, r=run, p=progress).render())


    def generate_index_page(self):
        d = {'programs': [(p, 'programs/%s' % p, self.programs[p].meta_info['description']['value']) for p in self.programs],
             'artifacts': [(a, 'artifacts/%s' % a) for a in self.artifacts_info],
             'meta': {'summary_type': 'programs_and_artifacts_index',
                      'version': '%s (%s)' % (anvio.anvio_version, anvio.anvio_codename),
                      'date': utils.get_date()}
            }

        d['program_provides_requires'] = self.get_program_requires_provides_dict(prefix='')

        output_file_path = os.path.join(self.output_directory_path, 'index.md')
        open(output_file_path, 'w').write(SummaryHTMLOutput(d, r=run, p=progress).render())


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
        self.init_programs()
        self.report_network()


    def report_network(self):
        artifact_names_seen = set([])
        artifacts_seen = Counter({})
        all_artifacts = []
        for program in self.programs.values():
            for artifact in program.meta_info['provides']['value'] + program.meta_info['requires']['value']:
                artifacts_seen[artifact.id] += 1
                if not artifact.id in artifact_names_seen:
                    all_artifacts.append(artifact)
                    artifact_names_seen.add(artifact.id)

        programs_seen = Counter({})
        for artifact in all_artifacts:
            for program in self.programs.values():
                for program_artifact in program.meta_info['provides']['value'] + program.meta_info['requires']['value']:
                    if artifact.name == program_artifact.name:
                        programs_seen[program.name] += 1

        network_dict = {"graph": [], "nodes": [], "links": [], "directed": False, "multigraph": False}

        node_indices = {}

        index = 0
        types_seen = set(["PROGRAM"])
        for artifact in all_artifacts:
            types_seen.add(artifact.type)
            network_dict["nodes"].append({"size": artifacts_seen[artifact.id],
                                          "score": 0.5 if artifact.internal else 1,
                                          "color": '#00AA00' if artifact.internal else "#AA0000",
                                          "id": artifact.id,
                                          "name": artifact.name,
                                          "internal": True if artifact.internal else False,
                                          "type": artifact.type})
            node_indices[artifact.id] = index
            index += 1

        for program in self.programs.values():
            network_dict["nodes"].append({"size": programs_seen[program.name],
                                          "score": 0.1,
                                          "color": "#AAAA00",
                                          "id": program.name,
                                          "name": program.name,
                                          "type": "PROGRAM"})
            node_indices[program.name] = index
            index += 1

        for artifact in all_artifacts:
            for program in self.programs.values():
                for artifact_provided in program.meta_info['provides']['value']:
                    if artifact_provided.id == artifact.id:
                        network_dict["links"].append({"source": node_indices[program.name], "target": node_indices[artifact.id]})
                for artifact_needed in program.meta_info['requires']['value']:
                    if artifact_needed.id == artifact.id:
                        network_dict["links"].append({"target": node_indices[program.name], "source": node_indices[artifact.id]})

        open(self.output_file_path, 'w').write(json.dumps(network_dict, indent=2))

        self.run.info('JSON description of network', self.output_file_path)
        self.run.info('Artifacts seen', ', '.join(sorted(list(types_seen))))


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
        self.init_programs(okay_if_no_meta = True, quiet = True)

        d = {}
        log_file = filesnpaths.get_temp_file_path()
        for i, program_name in enumerate(self.programs):
            program = self.programs[program_name]

            if program_name in self.programs_to_skip:
                run.warning("Someone doesn't want %s to be in the output :/ Fine. Skipping." % (program.name))

            progress.new('Bleep bloop')
            progress.update('%s (%d of %d)' % (program_name, i+1, len(self.programs)))

            output = utils.run_command_STDIN('%s --help --quiet' % (program.program_path), log_file, '').split('\n')

            if anvio.DEBUG:
                    usage, description, params, output = parse_help_output(output)
            else:
                try:
                    usage, description, params, output = parse_help_output(output)
                except Exception as e:
                    progress.end()
                    run.warning("The program '%s' does not seem to have the expected help menu output. Skipping to the next. "
                                "For the curious, this was the error message: '%s'" % (program.name, str(e).strip()))
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
