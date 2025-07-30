# -*- coding: utf-8
# pylint: disable=line-too-long
"""A library to help anvi'o describe itself"""

import os
import sys
import json
import copy
import argparse
import importlib

from collections import Counter

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.authors import AnvioAuthors
from anvio.docs import ANVIO_ARTIFACTS, ANVIO_WORKFLOWS, THIRD_PARTY_PROGRAMS
from anvio.summaryhtml import SummaryHTMLOutput


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


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
    all_lines = [l.strip() for l in open(file_path, 'r').readlines()]

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
        output.pop(0)
        section = output.pop(0)

        if section.startswith('━'):
            return None, None, None

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

    usage = J([l for l in get_until_blank(output)])

    if output.pop(0) != '':
        raise ConfigError("This output is missing the description start marker.")

    params = {}
    while 1:
        if not len(output):
            break

        section, desc, _params = get_param_set(output)

        if section == None:
            break

        if _params == '':
            pass
        else:
            params[section] = {'description': J(desc),
                               'params': _params}

    return usage,  params, output


class AnvioPrograms(AnvioAuthors):
    def __init__(self, args=None, r=terminal.Run(), p=terminal.Progress()):
        self.args = args
        self.run = r
        self.progress = p

        if not self.args:
             args = type('Args', (), {})()

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.program_names_to_focus = A("program_names_to_focus")

        # initiate `self.authors`
        AnvioAuthors.__init__(self, r=terminal.Run(verbose=False), p=self.progress)

        self.program_names_and_paths = self.get_anvio_program_names_and_their_paths()

        self.run.info("Anvi'o programs found and located", len(self.program_names_and_paths))

        if self.program_names_to_focus:
            self.program_names_to_focus = [p.strip() for p in self.program_names_to_focus.split(',')]

            self.run.info(" - Program names to focus", len(self.program_names_to_focus))

            # figure out what do we not want
            program_names_to_exclude = [p for p in self.program_names_and_paths if p not in self.program_names_to_focus]

            # remove them from the main dictionary
            [self.program_names_and_paths.pop(p) for p in program_names_to_exclude]

            if not len(self.program_names_and_paths):
                raise ConfigError("No anvi'o programs left to analyze after changing the focus to your list of program names. "
                                  "Probably there is a typo or something :/")

            self.run.info(" - Final number of programs kept", len(self.program_names_to_focus), mc="red")


    def sanity_check(self):
        """Check whether known programs and available programs make sense"""
        available_programs_according_to_python_environment = utils.get_available_program_names_in_active_environment(prefix='anvi-')
        available_programs_according_to_anvio = set(list(self.program_names_and_paths.keys()))

        programs_only_environment_knows_about = available_programs_according_to_python_environment - available_programs_according_to_anvio
        programs_only_anvio_knows_about = available_programs_according_to_anvio - available_programs_according_to_python_environment

        if programs_only_environment_knows_about or programs_only_anvio_knows_about:
                    self.run.warning("Please read the following lines carefully, since you may need to act on "
                                     "this information. There is a mismatch between the anvi'o programs the "
                                     "active anvi'o codebase knows about (through the entry points described "
                                     "in the `pyproject.toml`), and the anvi'o programs your Python environment "
                                     "knows about (through the list of programs accessible via $PATH).",
                                     header="FRIENDLY WARNING: ANVIO ENVIRONMENT IS CONFUSE", overwrite_verbose=True, lc='yellow')

                    if programs_only_environment_knows_about:
                        self.run.info_single("There are some anvi'o programs that are accessible in your Python environment, "
                                             "but your active codebase does not know about them. Here is a list of such "
                                             "programs:", overwrite_verbose=True, nl_after=1, level=0)

                        for program_name in programs_only_environment_knows_about:
                            self.run.info_single(program_name, overwrite_verbose=True, mc='red')

                        self.run.info_single(f"This can happen if you at some point had switched to an anvi'o branch where "
                                             f"these programs are described in the `pyproject.toml`, and ran `pip install -e .` "
                                             f"to install them to your environment, and then you switched to another branch with "
                                             f"a version of `pyproject.toml` that does not include these program names. This means, "
                                             f"if you were to run, let's say, '{list(programs_only_environment_knows_about)[0]}' "
                                             f"in your terminal right now, you would not get a 'command not found' error from your "
                                             f"shell, but a 'ModuleNotFoundError' error from Python.",
                                             overwrite_verbose=True, nl_after=1, level=0, nl_before=1)

                    if programs_only_anvio_knows_about:
                        self.run.info_single("There are some anvi'o programs that are known to your active anvi'o codebase, "
                                             "but they are not accessible to you in your Python environment. Here is a "
                                             "list of such programs:", overwrite_verbose=True, nl_after=1, level=0)

                        for program_name in programs_only_anvio_knows_about:
                            self.run.info_single(program_name, overwrite_verbose=True, mc='red')

                        self.run.info_single(f"This happens when you switch to a branch where there are new anvi'o programs described "
                                             f"in the `pyproject.toml` file, but they are not yet installed in the Python environment. "
                                             f"Which means, if you were to run, let's say, '{list(programs_only_anvio_knows_about)[0]}' "
                                             f"in your terminal right now, you would get a 'command not found' error from your shell "
                                             f"(rather than a 'ModuleNotFoundError' error from Python).",
                                             overwrite_verbose=True, nl_before=1, nl_after=1, level=0)

                    self.run.info_single("The universal solution here is to run the following command right now in your anvi'o "
                                         "source code directory:", overwrite_verbose=True, nl_after=1, level=0)
                    self.run.info_single("    pip install -e . --force-reinstall --upgrade",
                                         overwrite_verbose=True, nl_after=1, level=0, pretty_indentation=False)
                    self.run.info_single("This will synchronize your anvi'o codebase with its installed version in your active "
                                         "Python environment. This is indeed very annoying, since you will likely have to do it "
                                         "again when you go back to another branch, but this is how it goes. It is also a "
                                         "viable alternative to ignore this message, if you think this mismatch is not a concern "
                                         "for you at this stage.", overwrite_verbose=True, level=0)


    def get_anvio_program_names_and_their_paths(self):
        """Parses the package pyproject.toml file and returns a dictionary that links program names to
           program Python files under `anvio/cli`
        """

        program_names_and_paths = {}

        try:
            import tomli as tomllib
        except ImportError:
            raise ConfigError("The AnvioPrograms class needs `tomli` to be available in this Python environment. You can "
                              "simply install it by running `pip install tomli` (and hope for the best).")

        anvio_dir = os.path.dirname(anvio.__file__)
        pyproject_path = os.path.join(os.path.dirname(anvio_dir), 'pyproject.toml')

        if not os.path.exists(pyproject_path):
            raise ConfigError("The pyproject.toml for the anvi'o package does not seem to be where it is expected :/ "
                              "Without that file, this class cannot associate anvi'o program names to the actual "
                              "Python files that implement them :(")

        # read the contensts of the pyproject.toml
        with open(pyproject_path, 'rb') as f:
            pyproject_data = tomllib.load(f)

        ########################################
        # figure out entry points in the file
        ########################################
        entry_points = pyproject_data.get('project', {}).get('scripts', {})

        if not entry_points:
            raise ConfigError("The pyproject.toml is there, but it does not seem to contain any entry points. This "
                              "function needs an adult to figure this out :(")

        # turn entry point entries into program name / absolute path pairs
        for entry_point in entry_points:
            program_name = entry_point
            program_path = entry_points[entry_point]

            program_relative_path = program_path.split(':')[0].replace('.', '/') + '.py'
            program_path = os.path.abspath(os.path.join(anvio_dir, '..', program_relative_path))

            if not os.path.exists(program_path):
                raise ConfigError(f"Parsing the entry points in the pyproject.toml file did not lead to an actual "
                                  f"program path for `{program_name}` as there was nothing at `{program_path}` :/ "
                                  f"This should have never happened, but must be solved before this program can "
                                  f"continue doing its job.")

            program_names_and_paths[program_name] = program_path

        ########################################
        # figure out non-python scripts
        ########################################
        non_python_scripts = pyproject_data.get('tool', {}).get('setuptools', {})['script-files']

        for non_python_script in non_python_scripts:
            program_name = os.path.basename(non_python_script)
            program_path = os.path.abspath(os.path.join(anvio_dir, '..', non_python_script))

            program_names_and_paths[program_name] = program_path

        return program_names_and_paths


    def init_programs(self, okay_if_no_meta=False, always_include_those_with_docs=True, quiet=False):
        """Initializes the `self.programs` dictionary."""

        num_all_programs = len(self.program_names_and_paths)

        self.programs = {}
        self.progress.new('Characterizing program', progress_total_items=num_all_programs)

        programs_with_usage_info = set([])
        programs_without_usage_info = set([])
        programs_with_provides_requires_info = set([])
        programs_without_provides_requires_info = set([])

        for program_name in self.program_names_and_paths:
            program_filepath = self.program_names_and_paths[program_name]
            self.progress.update(os.path.basename(program_name), increment=True)

            program = Program(program_name, program_filepath, r=self.run, p=self.progress)

            program_usage_information_path = os.path.join(anvio.DOCS_PATH, 'programs/%s.md' % (program.name))

            if program.meta_info['provides']['value'] or program.meta_info['requires']['value']:
                programs_with_provides_requires_info.add(program.name)
            else:
                programs_without_provides_requires_info.add(program.name)

            # learn about the usage statement of the program if you have access to reading
            # markdown reader function:
            if hasattr(self, 'read_anvio_markdown'):
                if os.path.exists(program_usage_information_path):
                    program.usage = self.read_anvio_markdown(program_usage_information_path)
                    programs_with_usage_info.add(program.name)
                else:
                    programs_without_usage_info.add(program.name)

            keep_program = True
            if not (program.meta_info['provides']['value'] or program.meta_info['requires']['value']):
                # if we are here, it means the program is missing both provides AND requires statements.
                # If the user hasn't set `okay_if_no_meta=True`, we're going to get rid of them, and
                # will NOT include them in `self.programs`
                if not okay_if_no_meta:
                    keep_program = False

                # BUT, there are programs that have no provides/requires statements, such as anvi-self-test,
                # but have a usage statement under docs already, we may want to keep them in the list
                # regardless. so here we test that:
                if program.name in programs_with_usage_info and always_include_those_with_docs:
                    keep_program = True

            if keep_program:
                # include the program in our final list
                self.programs[program.name] = program
            else:
                # forget all about it
                try:
                    programs_with_usage_info.remove(program.name)
                    programs_without_usage_info.remove(program.name)
                except:
                    pass

        self.progress.end()

        # here we will go through each program, and see if there are any with no author information
        programs_with_no_authors = [p for p in self.programs if not len(self.programs[p].meta_info['authors']['value'])]
        if len(programs_with_no_authors) and anvio.DEBUG:
            self.run.warning(f"The following programs have no `__authors__` tag: {', '.join(programs_with_no_authors)}.")

        # here we will go through each program, and see if there are any with authors that
        # are not described in the AUTHORS.yaml file:
        programs_with_unknown_authors = set([])
        unknown_authors = set([])
        known_authors_in_programs = set([])
        for program in self.programs:
            for author in self.programs[program].meta_info['authors']['value']:
                if author not in self.authors:
                    programs_with_unknown_authors.add(program)
                    unknown_authors.add(author)
                else:
                    known_authors_in_programs.add(author)

        if len(programs_with_unknown_authors):
            raise ConfigError(f"The following programs have authors who are not defined in the authors YAML file "
                              f"({', '.join(unknown_authors)}) -- every author listed in anvi'o programs must have "
                              f"an entry in the authors YAML file: {', '.join(programs_with_unknown_authors)}.")

        # report missing provides/requires information
        self.run.info_single("Of %d programs found, %d did not contain PROVIDES AND/OR REQUIRES "
                             "statements :/ This may be normal for some programs, but here is the "
                             "complete list of those that are missing __provides__ and __requires__ "
                             "tags in their code in case you see something you can complete: '%s'." % \
                                        (len(self.program_names_and_paths),
                                         len(programs_without_provides_requires_info),
                                         ', '.join(programs_without_provides_requires_info)),
                             nl_after=1, nl_before=1)

        # report missing provides/requires information
        self.run.info_single("Of %d programs found, %d did not have any PROVIDES/REQUIRES statements. You can "
                             "help by adding usage information for programs by creating markdown "
                             "formatted files under the directory '%s'. Please see examples in anvi'o "
                             "codebase: https://github.com/merenlab/anvio/tree/master/anvio/docs. "
                             "Here is a complete list of programs that are missing usage statements: %s " % \
                                        (len(self.program_names_and_paths),
                                         len(programs_without_provides_requires_info),
                                         anvio.DOCS_PATH,
                                         ', '.join(programs_without_provides_requires_info)),
                             nl_after=1, nl_before=1)


class Program:
    def __init__(self, program_name, program_path, r=terminal.Run(), p=terminal.Progress()):
        self.run = r
        self.progress = p

        self.name = program_name
        self.program_path = program_path
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
            'authors': {
                'object_name': '__authors__',
                'null_object': []
            },
            'anvio_workflows': {
                'object_name': '__anvio_workflows__',
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

            # Lower case the github usernames
            if info_type == "authors":
                info = [a.lower() for a in info]

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

    def __init__(self, artifact_id, provided_by_anvio=True, optional=True, single=True):
        if artifact_id not in ANVIO_ARTIFACTS:
            progress.reset()
            raise ConfigError("Ehem. Anvi'o does not know about artifact '%s'. There are two was this could happen: "
                              "one, you've made a typo (easy to fix), two, you've just updated __provides__ or __requires__ "
                              "statement in an anvi'o program with an artifact that does not exist and have not yet updated "
                              "`anvio/docs/__init__.py` (which is also easy to fix). Please consider also adding a description of "
                              "this artifact under anvio/docs/artifacts while you are at it :)" % artifact_id)

        artifact = ANVIO_ARTIFACTS[artifact_id]
        self.id = artifact_id
        self.name = artifact['name']
        self.type = artifact['type']
        self.provided_by_anvio = artifact['provided_by_anvio']
        self.provided_by_user = artifact['provided_by_user']

        # attributes set by the context master
        self.single = single
        self.optional = optional


    def __repr__(self):
        return "ARTIFACT::%s" % self.id


class AnvioWorkflows:
    """Information on anvi'o workflows"""

    def __init__(self, args, r=terminal.Run(), p=terminal.Progress()):
        self.args = args
        self.run = r
        self.progress = p

        self.workflows = {}

        if not hasattr(self, 'programs') or not hasattr(self, 'artifacts_info') or not hasattr(self, 'authors'):
            raise ConfigError("AnvioWorkflows class is upset. You need to treat this class as a base class, and initialize "
                              "it from within another class that has already initialized AnvioPrograms, AnvioArtifacts, AND"
                              "AnvioAuthors classes. If this is confusing, take a look at the AnvioDocs class.")

        if not len(self.programs):
            raise ConfigError("AnvioWorkflows is being initialized with a blank `self.programs` variable :/")

        if not len(self.artifacts_info):
            raise ConfigError("AnvioWorkflows is being initialized with a blank `self.artifacts_info` variable :/")


    def init_workflows(self):
        """Learn all about anvi'o workflows and initiate the class.

        Returns
        =======
        workflows, dict:
            Running this function will fill in the dictionary `self.workflows`
        """

        self.workflows= {}

        expected_keys = ['authors', 'artifacts_produced', 'artifacts_accepted', 'anvio_workflows_inherited', 'third_party_programs_used', 'one_sentence_summary', 'one_paragraph_summary']

        workflows_without_descriptions = set([])

        for workflow in ANVIO_WORKFLOWS:
            self.workflows[workflow] = copy.deepcopy(ANVIO_WORKFLOWS[workflow])

            for key in expected_keys:
                if key not in self.workflows[workflow]:
                    raise ConfigError(f"Every workflow must contain the keys \"{', '.join(expected_keys)}\". But "
                                      f"it is not the case for the workflow '{workflow}' :(")

            self.workflows[workflow]['name'] = workflow
            self.workflows[workflow]['anvio_programs_used'] = []

            # a workflow description includes the list of third party programs that are
            # optionally used from whithin the workflow. here we will sanity check that
            # they all have descriptions in `THIRD_PARTY_PROGRAMS`
            # dictionary.
            for purpose, program_names in self.workflows[workflow]['third_party_programs_used']:
                for program_name in program_names:
                    if program_name not in THIRD_PARTY_PROGRAMS:
                        raise ConfigError(f"The workflow {workflow} lists the program '{program_name}' in its "
                                          f"description for third-party programs that are used from within, "
                                          f"however, there is no entry for this program in the variable "
                                          f"'THIRD_PARTY_PROGRAMS' in the file "
                                          f"'anvio/docs/__init__.py'. Please add a necessary description for "
                                          f"this program into that dict, and try this again.")

            # learn about the description of the workflow
            workflow_description_path = os.path.join(anvio.DOCS_PATH, 'workflows/%s.md' % (workflow))
            if os.path.exists(workflow_description_path):
                self.workflows[workflow]['description'] = self.read_anvio_markdown(workflow_description_path)
            else:
                workflows_without_descriptions.add(workflow)

            # add anvi'o programs used by the workflow
            for program in self.programs.values():
                    if workflow in [a for a in program.meta_info['anvio_workflows']['value']]:
                        self.workflows[workflow]['anvio_programs_used'].append(program.name)

            # make sure 'workflow-config' artifact is not included in 'artifacts_accepted'
            if 'workflow-config' in self.workflows[workflow]['artifacts_accepted']:
                raise ConfigError(f"The 'artifacts_accepted' description for the workflow '{workflow}' includes "
                                  f"`workflow-config`, but this artifact is automatically added to each workflow "
                                  f"on-the-fly, thus, it shouldn't be listed in the workflow description :/ Sorry!")

            # every workflow should accept the artifact `workflow-config` by default, so add it here:
            self.workflows[workflow]['artifacts_accepted'] = ['workflow-config'] + self.workflows[workflow]['artifacts_accepted']

            # sanity check artifacts accepted:
            for artifact_name in self.workflows[workflow]['artifacts_accepted']:
                if artifact_name not in self.artifacts_info:
                    raise ConfigError(f"The artifact '{artifact_name}' that is listed as one of the artifacts the workflow "
                                      f"{workflow} accepts does not seem to be an artifact anvi'o knows about :/ If this is "
                                      f"a new artifact for workflow, please first describe it in the dictionary `ANVIO_ARTIFACTS` "
                                      f"in anvio/docs/__init__.py")

        # sanity check of author names
        author_names_appear_in_workflows = set([])
        [author_names_appear_in_workflows.update(w['authors']) for w in self.workflows.values()]
        author_names_missing_in_authors_file = [a for a in author_names_appear_in_workflows if a not in self.authors]
        if len(author_names_missing_in_authors_file):
            self.run.warning(None, header="SOME SNAFU TOOK PLACE [poop emoji]")
            self.run.info("Author names anvi'o knows about", ', '.join(self.authors), mc='green')
            self.run.info("Author names anvi'o does not know about", ', '.join(author_names_missing_in_authors_file), mc='red')
            raise ConfigError("Some author names in anvi'o workflows defined under `anvio/docs/__init__.py` do not "
                              "appear in the DEVELOPERS.yaml file. If there is no typo here, please update the "
                              "contents of the DEVELOPERS.yaml file with the GitHub username of the developer you "
                              "wish to associate with a workflow. The problematic authors are shown above.")

        # make sure every workflow has at least one author
        workflows_missing_authors = set([])
        [workflows_missing_authors.add(w) for w in self.workflows if not len(self.workflows[w]['authors'])]
        if len(workflows_missing_authors):
            raise ConfigError(f"One or more workflows defined under `anvio/docs/__init__.py` do not have "
                              f"any authors. Every workflow must have at least one :/ Here is the list of those that "
                              f"are missing any authors: {', '.join(workflows_missing_authors)}")

        # note the workflows that are missing descriptions.
        if len(workflows_without_descriptions):
            self.run.info_single("Of %d workflows found, %d did not contain any DESCRIPTION. If you would like to "
                                 "see examples and add new descriptions, please see the directory '%s'. Here is the "
                                 "full list of workflows that are not yet explained: %s." \
                                        % (len(ANVIO_WORKFLOWS),
                                           len(workflows_without_descriptions),
                                           anvio.DOCS_PATH,
                                           ', '.join(workflows_without_descriptions)), nl_after=1, nl_before=1)

        # makes ure every workflow mentioned in programs in fact are explained
        # as a workflow
        workflows_mentioned_in_programs = set([])
        [workflows_mentioned_in_programs.update(p.meta_info['anvio_workflows']['value']) for p in self.programs.values()]
        unknown_workflows_mentioned_in_programs = [w for w in workflows_mentioned_in_programs if w not in self.workflows]
        if len(unknown_workflows_mentioned_in_programs):
            raise ConfigError(f"Some anvi'o programs include `__anvio_workflows__` tags with workflow names anvi'o "
                              f"dees not recognize :/ Here is the missing workflow names so you can either fix some "
                              f"typos, or add entries for these workflows in `anvio/docs/__init.py__`: "
                              f"{', '.join(unknown_workflows_mentioned_in_programs)}")


class AnvioArtifacts:
    """Information on anvi'o artifacts"""

    def __init__(self, args, r=terminal.Run(), p=terminal.Progress()):
        self.args = args
        self.run = r
        self.progress = p

        self.artifacts_info = {}
        self.artifact_types = {}

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
            self.artifacts_info[artifact] = {'required_by': [], 'provided_by': [], 'description': None, 'type': ANVIO_ARTIFACTS[artifact]['type']}

            # learn about the description of the artifact
            artifact_description_path = os.path.join(anvio.DOCS_PATH, 'artifacts/%s.md' % (artifact))
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

            # register artifact type
            artifact_type = self.artifacts_info[artifact]['type']

            if artifact_type not in self.artifact_types:
                self.artifact_types[artifact_type] = []

            self.artifact_types[artifact_type].append(artifact)

        if len(artifacts_without_descriptions):
            self.run.info_single("Of %d artifacts found, %d did not contain any DESCRIPTION. If you would like to "
                                 "see examples and add new descriptions, please see the directory '%s'. Here is the "
                                 "full list of artifacts that are not yet explained: %s." \
                                        % (len(ANVIO_ARTIFACTS),
                                           len(artifacts_without_descriptions),
                                           anvio.DOCS_PATH,
                                           ', '.join(artifacts_without_descriptions)), nl_after=1, nl_before=1)


class AnvioDocs(AnvioPrograms, AnvioArtifacts, AnvioWorkflows):
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

        if not os.path.exists(anvio.DOCS_PATH):
            raise ConfigError("The anvi'o docs path is not where it should be :/ Something funny is going on.")

        filesnpaths.gen_output_directory(self.output_directory_path, delete_if_exists=True, dont_warn=True)

        self.artifacts_output_dir = filesnpaths.gen_output_directory(os.path.join(self.output_directory_path, 'artifacts'))
        self.programs_output_dir = filesnpaths.gen_output_directory(os.path.join(self.output_directory_path, 'programs'))
        self.workflows_output_dir = filesnpaths.gen_output_directory(os.path.join(self.output_directory_path, 'workflows'))

        self.version_short_identifier = 'm' if anvio.anvio_version_for_help_docs == 'main' else anvio.anvio_version_for_help_docs
        self.base_url = os.path.join("/help", anvio.anvio_version_for_help_docs)
        self.anvio_markdown_variables_conversion_dict = {}

        AnvioPrograms.__init__(self, args, r=self.run, p=self.progress)
        self.init_programs()

        AnvioArtifacts.__init__(self, args, r=self.run, p=self.progress)
        self.init_artifacts()

        AnvioWorkflows.__init__(self, args, r=self.run, p=self.progress)
        self.init_workflows()

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

        self.generate_pages_for_workflows()

        self.generate_index_page()


    def copy_images(self):
        """Copies images from the codebase to the output directory"""

        utils.shutil.copytree(self.images_source_directory, os.path.join(self.output_directory_path, 'images'))

        os.makedirs(os.path.join(self.output_directory_path, 'images/authors'))

        for author in self.authors:
            utils.shutil.copy(self.authors[author]['avatar'], os.path.join(self.output_directory_path, 'images/authors', os.path.basename(self.authors[author]['avatar'])))


    def init_anvio_markdown_variables_conversion_dict(self):
        for program_name in self.program_names_and_paths:
            self.anvio_markdown_variables_conversion_dict[program_name] = """<span class="artifact-p">[%s](%s/programs/%s)</span>""" % (program_name, self.base_url, program_name)

        for artifact_name in ANVIO_ARTIFACTS:
            self.anvio_markdown_variables_conversion_dict[artifact_name] = """<span class="artifact-n">[%s](%s/artifacts/%s)</span>""" % (artifact_name, self.base_url, artifact_name)


    def read_anvio_markdown(self, file_path):
        """Reads markdown descriptions filling in anvi'o variables.

        Basically a lot of l_l83Я 1337 Я0XX0ЯZ stuff's going on down there, so you better run while you can.
        """

        filesnpaths.is_file_plain_text(file_path)

        if not len(self.anvio_markdown_variables_conversion_dict):
            self.init_anvio_markdown_variables_conversion_dict()

        markdown_content = open(file_path).read()

        # this is quite a big deal thing to do here:
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
            additional_info = ("If you're stumped by that message, here are some common errors and their solutions: "
                               "(1) 'unsupported format character' could mean that one of your tags specified with "
                               "'%(tag)s' did not have the appended 's'. (2) 'not enough arguments for format string' "
                               "could mean that your document has a '%' sign used in natural language, i.e. '85% similar'. "
                               "This must be replaced with '85%% similar'.")
            raise ConfigError("Something went wrong while working with '%s' :/ This is what we know: '%s'. %s" % (file_path, e, additional_info))

        # now we have replaced anvi'o variables with markdown links, it is time to replace
        # hyphens in anvi'o codeblocks with HTML hyphens so markdown does not freakout when it is
        # time to visualize these and replace -- characters with en dash.
        markdwon_lines = markdown_content.split('\n')
        line_nums_for_codestart_tags = [i for i in range(0, len(markdwon_lines)) if markdwon_lines[i].strip() == "{{ codestart }}"]
        line_nums_for_codestop_tags = [i for i in range(0, len(markdwon_lines)) if markdwon_lines[i].strip() == "{{ codestop }}"]

        if len(line_nums_for_codestart_tags) != len(line_nums_for_codestop_tags):
            raise ConfigError("In %s, the number of {{ codestart }} tags do not match to the number of {{ codestop }} tags :/" % file_path)


        for line_start, line_end in list(zip(line_nums_for_codestart_tags, line_nums_for_codestop_tags)):
            for line_num in range(line_start + 1, line_end):
                markdwon_lines[line_num] = markdwon_lines[line_num].replace("-", "&#45;").replace("*", "&#42;").replace("==", "&#61;&#61;")

        # all lines are processed: merge them back into a single text:
        markdown_content = '\n'.join(markdwon_lines)

        # now we have a proper markdown, it is time to remove anvi'o {{ codestart }} and {{ codestop }} blocks.
        markdown_content = markdown_content.replace("""{{ codestart }}""", """<div class="codeblock" markdown="1">""")
        markdown_content = markdown_content.replace("""{{ codestop }}""", """</div>""")

        # return it like a pro.
        return markdown_content


    def get_program_requires_provides_dict(self, prefix="../../"):
        d = {}

        for program_name in self.programs:
            d[program_name] = {}

            program = self.programs[program_name]
            d[program_name]['requires'] = [(r.id, '%sartifacts/%s' % (prefix, r.id)) for r in program.meta_info['requires']['value']]
            d[program_name]['provides'] = [(r.id, '%sartifacts/%s' % (prefix, r.id)) for r in program.meta_info['provides']['value']]
            d[program_name]['anvio_workflows'] = [(w, '%sworkflows/%s' % (prefix, w)) for w in program.meta_info['anvio_workflows']['value']]

        return d


    def get_workflow_produced_artifacts_list(self, workflow_name, prefix="../../"):
        return [(r, '%sartifacts/%s' % (prefix, r)) for r in self.workflows[workflow_name]['artifacts_produced']]


    def get_workflow_accepted_artifacts_list(self, workflow_name, prefix="../../"):
        return [(r, '%sartifacts/%s' % (prefix, r)) for r in self.workflows[workflow_name]['artifacts_accepted']]


    def generate_pages_for_artifacts(self):
        """Generates static pages for artifacts in the output directory"""

        self.progress.new("Rendering artifact pages", progress_total_items=len(ANVIO_ARTIFACTS))
        self.progress.update('...')

        for artifact in ANVIO_ARTIFACTS:
            self.progress.update(f"'{artifact}' ...", increment=True)

            d = {'artifact': ANVIO_ARTIFACTS[artifact],
                 'meta': {'summary_type': 'artifact',
                          'version': '\n'.join(['|%s|%s|' % (t[0], t[1]) for t in anvio.get_version_tuples()]),
                          'date': utils.get_date(),
                          'version_short_identifier': self.version_short_identifier}
                }

            d['artifact']['name'] = artifact
            d['artifact']['required_by'] = [(r, '../../programs/%s' % r) for r in self.artifacts_info[artifact]['required_by']]
            d['artifact']['provided_by'] = [(r, '../../programs/%s' % r) for r in self.artifacts_info[artifact]['provided_by']]
            d['artifact']['description'] = self.artifacts_info[artifact]['description']
            d['artifact']['icon'] = '../../images/icons/%s.png' % ANVIO_ARTIFACTS[artifact]['type']

            if anvio.DEBUG:
                self.progress.reset()
                run.warning(None, 'THE OUTPUT DICT')
                import json
                print(json.dumps(d, indent=2))

            self.progress.update(f"'{artifact}' ... rendering ...", increment=False)
            artifact_output_dir = filesnpaths.gen_output_directory(os.path.join(self.artifacts_output_dir, artifact))
            output_file_path = os.path.join(artifact_output_dir, 'index.md')
            open(output_file_path, 'w').write(SummaryHTMLOutput(d, r=run, p=progress).render())

        self.progress.end()


    def get_HTML_formatted_authors_data(self, authors):
        """for a given program, returns HTML-formatted authors data"""

        d = ""

        for author in authors:
            d += '''<div class="anvio-person"><div class="anvio-person-info">'''
            d += f'''<div class="anvio-person-photo"><img class="anvio-person-photo-img" src="../../images/authors/{os.path.basename(self.authors[author]['avatar'])}" /></div>'''
            d += '''<div class="anvio-person-info-box">'''
            d += f'''<a href="/people/{self.authors[author]['github']}" target="_blank"><span class="anvio-person-name">{self.authors[author]['name']}</span></a>'''
            d += '''<div class="anvio-person-social-box">'''

            if 'web' in self.authors[author]:
                d += f'''<a href="{self.authors[author]['web']}" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a>'''

            d += f'''<a href="mailto:{self.authors[author]['email']}" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a>'''

            if 'twitter' in self.authors[author]:
                d += f'''<a href="http://twitter.com/{self.authors[author]['twitter']}" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a>'''

            d += f'''<a href="http://github.com/{self.authors[author]['github']}" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a>'''

            d += '''</div></div></div></div>\n\n'''

        return d


    def get_HTML_formatted_authors_data_mini(self, authors):
        """for a given list of authors, returns a tiny version of the HTML-formatted authors data"""

        d = ""

        for author in authors:
            d += '''<div class="anvio-person-mini"><div class="anvio-person-photo-mini">'''
            d += f'''<a href="/people/{self.authors[author]['github']}" target="_blank"><img class="anvio-person-photo-img-mini" title="{self.authors[author]['name']}" src="images/authors/{os.path.basename(self.authors[author]['avatar'])}" /></a>'''
            d += '''</div></div>\n'''

        return d


    def get_HTML_formatted_third_party_programs(self, workflow_name):
        """Get a template-friendly list of third-party programs used from within a workflow"""

        d = []

        for purpose, program_names in self.workflows[workflow_name]['third_party_programs_used']:
            for program_name in program_names:
                d.append(f'''<a href="{THIRD_PARTY_PROGRAMS[program_name]['link']}" target="_blank">{program_name}</a> ({purpose})''')

        return d


    def generate_pages_for_workflows(self):
        """Generate static pages for anvi'o workflows in the output directory"""

        self.progress.new("Rendering workflow pages", progress_total_items=len(self.workflows))
        self.progress.update('...')

        for workflow_name in self.workflows:
            self.progress.update(f"'{workflow_name}' ...", increment=True)

            d = {'workflow': self.workflows[workflow_name],
                 'meta': {'summary_type': 'workflow',
                          'version': '\n'.join(['|%s|%s|' % (t[0], t[1]) for t in anvio.get_version_tuples()]),
                          'date': utils.get_date(),
                          'version_short_identifier': self.version_short_identifier}
                 }

            d['workflow']['artifacts_produced'] = self.get_workflow_produced_artifacts_list(workflow_name)
            d['workflow']['artifacts_accepted'] = self.get_workflow_accepted_artifacts_list(workflow_name)
            d['workflow']['third_party_programs_used'] = self.get_HTML_formatted_third_party_programs(workflow_name)
            d['workflow']['authors'] = self.get_HTML_formatted_authors_data(d['workflow']['authors'])

            # also add information regarding the artifacts
            d['artifacts'] = self.artifacts_info

            if anvio.DEBUG:
                self.progress.reset()
                run.warning(None, 'THE WORKFLOW OUTPUT DICT')
                import json
                print(json.dumps(d, indent=2))

            self.progress.update(f"'{workflow_name}' ... rendering ...", increment=False)
            workflow_output_dir = filesnpaths.gen_output_directory(os.path.join(self.workflows_output_dir, workflow_name))
            output_file_path = os.path.join(workflow_output_dir, 'index.md')
            open(output_file_path, 'w').write(SummaryHTMLOutput(d, r=run, p=progress).render())

        self.progress.end()



    def generate_pages_for_programs(self):
        """Generates static pages for programs in the output directory"""

        self.progress.new("Rendering program pages", progress_total_items=len(self.programs))
        self.progress.update('...')

        program_provides_requires_dict = self.get_program_requires_provides_dict()

        for program_name in self.programs:
            self.progress.update(f"'{program_name}' ...", increment=True)

            program = self.programs[program_name]
            d = {'program': {},
                 'meta': {'summary_type': 'program',
                          'version': '\n'.join(['|%s|%s|' % (t[0], t[1]) for t in anvio.get_version_tuples()]),
                          'date': utils.get_date(),
                          'version_short_identifier': self.version_short_identifier}
                }

            d['program']['name'] = program_name
            d['program']['usage'] = program.usage
            d['program']['description'] = program.meta_info['description']['value']
            d['program']['resources'] = program.meta_info['resources']['value']
            d['program']['requires'] = program_provides_requires_dict[program_name]['requires']
            d['program']['provides'] = program_provides_requires_dict[program_name]['provides']
            d['program']['icon'] = '../../images/icons/%s.png' % 'PROGRAM'
            d['program']['authors'] = self.get_HTML_formatted_authors_data(program.meta_info['authors']['value'])
            d['artifacts'] = self.artifacts_info
            d['workflows'] = self.workflows

            if anvio.DEBUG:
                self.progress.reset()
                run.warning(None, 'THE OUTPUT DICT')
                import json
                print(json.dumps(d, indent=2))

            self.progress.update(f"'{program_name}' ... rendering ...", increment=False)
            program_output_dir = filesnpaths.gen_output_directory(os.path.join(self.programs_output_dir, program_name))
            output_file_path = os.path.join(program_output_dir, 'index.md')
            open(output_file_path, 'w').write(SummaryHTMLOutput(d, r=run, p=progress).render())

            # create the program network, too
            self.progress.update(f"'{program_name}' ... rendering ... network json ...", increment=False)
            program_output_dir = filesnpaths.gen_output_directory(os.path.join(self.programs_output_dir, program_name))
            program_network = ProgramsNetwork(argparse.Namespace(output_file=os.path.join(program_output_dir, "network.json"), program_names_to_focus=program_name), r=terminal.Run(verbose=False))
            program_network.generate()

        self.progress.end()


    def generate_index_page(self):
        """Generates the index page for help where all programs and artifacts are listed"""

        self.progress.new("Index page")
        self.progress.update('...')

        # let's add the 'path' for each artifact to simplify
        # access from the template:
        for artifact in self.artifacts_info:
            self.artifacts_info[artifact]['path'] = f"artifacts/{artifact}"

        # quick update of the author information in workflows so they contain nice HTML
        # code instad of a list of author names
        for workflow in self.workflows:
            self.workflows[workflow]['authors'] = self.get_HTML_formatted_authors_data_mini(ANVIO_WORKFLOWS[workflow]['authors'])

        # please note that artifacts get a fancy dictionary with everything, while programs get a crappy tuples list.
        # if we need to improve the functionality of the help index page, we may need to update programs
        # to a fancy dictionary, too.
        d = {'programs': [(p, 'programs/%s' % p, self.programs[p].meta_info['description']['value'], self.get_HTML_formatted_authors_data_mini(self.programs[p].meta_info['authors']['value'])) for p in self.programs],
             'workflows': self.workflows,
             'artifacts': self.artifacts_info,
             'artifact_types': self.artifact_types,
             'meta': {'summary_type': 'programs_and_artifacts_index',
                      'version': '%s (%s)' % (anvio.anvio_version, anvio.anvio_codename),
                      'date': utils.get_date()}
            }

        d['program_provides_requires'] = self.get_program_requires_provides_dict(prefix='')

        self.progress.update('Rendering...')
        output_file_path = os.path.join(self.output_directory_path, 'index.md')

        self.progress.update('Writing...')
        open(output_file_path, 'w').write(SummaryHTMLOutput(d, r=run, p=progress).render())

        self.progress.end()


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


    def report_network(self, program_name=None):
        """Reports an association network for anvi'o programs and artifacs

        By default this function will report a network for all programs, unless the user
        set args.program_names_to_focus to a list of programs, in which case the function
        will focus only on that program and its artifats, reporting a sub-network.
        """

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
                                          "score": 0.5 if artifact.provided_by_anvio else 1,
                                          "color": '#00AA00' if artifact.provided_by_anvio else "#AA0000",
                                          "id": artifact.id,
                                          "name": artifact.name,
                                          "provided_by_anvio": True if artifact.provided_by_anvio else False,
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
                    usage, params, output = parse_help_output(output)
            else:
                try:
                    usage, params, output = parse_help_output(output)
                except Exception as e:
                    progress.end()
                    run.warning("The program '%s' does not seem to have the expected help menu output. Skipping to the next. "
                                "For the curious, this was the error message: '%s'" % (program.name, str(e).strip()))
                    continue

            d[program.name] = {'usage': usage,
                               'description': program.meta_info['description']['value'],
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


