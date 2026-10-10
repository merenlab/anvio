"""Render outputs from a documentation dataset.

HelpPagesRenderer builds the help website and LLMsTextRenderer builds llms-full.txt.
Rendering uses only the dataset, its image assets, and installed website templates.
It does not discover CLI programs or consult documentation source registries.
"""

from __future__ import annotations

import os
import re
import json
import shutil
from collections import Counter
from dataclasses import asdict
from pathlib import Path
from typing import Any

import anvio
import anvio.filesnpaths as filesnpaths
import anvio.terminal as terminal
import anvio.utils as utils

from anvio.docs.data import (
    DATASET_FILENAME, ArtifactDocumentation, DocumentationData, DocumentRecord,
    ProgramDocumentation, WorkflowDocumentation,
)
from anvio.errors import ConfigError, FilesNPathsError
from anvio.summaryhtml import SummaryHTMLOutput


run = terminal.Run()
progress = terminal.Progress()


class HelpPagesRenderer:
    """Generate a docs output.

    The purpose of this class is to generate a static HTML output with
    interlinked files that serve as the primary documentation for anvi'o
    programs, input files they expect, and output files the generate.

    The default client of this class is `anvi-script-gen-help-docs`.
    """


    def __init__(self, dataset: DocumentationData, output_directory: str, r=run, p=progress):
        self.dataset = dataset
        self.output_directory = Path(output_directory).resolve()
        self.run = r
        self.progress = p
        self.programs = dataset.programs
        self.workflows = dataset.workflows
        self.authors = dataset.authors
        self.program_names_and_paths = dataset.program_source_paths
        self.version_short_identifier = dataset.meta.version_short_identifier
        self.base_url = dataset.meta.base_url
        self.anvio_markdown_variables_conversion_dict = {}

        # Resolve every body before replacing any existing output.
        self.bodies = {record.id: self.render_body(record)
                       for records in (dataset.programs, dataset.artifacts, dataset.workflows)
                       for record in records.values()}

        # The existing templates consume dictionaries rather than documentation records.
        self.artifacts_info = {}
        self.artifact_types = {}
        for artifact, record in dataset.artifacts.items():
            self.artifact_types.setdefault(record.type, []).append(artifact)
            self.artifacts_info[artifact] = {'type': record.type, 'description': self.bodies[record.id]}
            for relation, inverse in (('requires', 'required_by'), ('provides', 'provided_by'),
                                      ('can_use', 'can_used_by'), ('can_provide', 'can_provided_by')):
                self.artifacts_info[artifact][inverse] = [name for name, program in self.programs.items()
                                                        if artifact in getattr(program, relation)]


    def programs_network(self, program_names: list[str] | None = None, artifact_names_as_ids: bool = False) -> dict[str, Any]:
        """Returns an association network for anvi'o programs and artifacs

        By default this function will report a network for all programs, unless the user
        passed program_names as a list of programs, in which case the function
        will focus only on that program and its artifats, reporting a sub-network.
        """

        programs = {name: program for name, program in self.programs.items()
                    if (program_names is None or name in program_names) and (program.requires or program.provides)}
        artifact_names = {name: name if artifact_names_as_ids else record.name
                          for name, record in self.dataset.artifacts.items()}

        artifact_names_seen = set([])
        artifacts_seen = Counter({})
        all_artifacts = []
        for program in programs.values():
            all_program_artifacts = (program.provides +
                                     program.requires +
                                     program.can_use +
                                     program.can_provide)
            for artifact in all_program_artifacts:
                artifacts_seen[artifact] += 1
                if not artifact in artifact_names_seen:
                    all_artifacts.append(artifact)
                    artifact_names_seen.add(artifact)

        programs_seen = Counter({})
        for artifact in all_artifacts:
            for program in programs.values():
                all_program_artifacts = (program.provides +
                                         program.requires +
                                         program.can_use +
                                         program.can_provide)
                for program_artifact in all_program_artifacts:
                    if artifact_names[artifact] == artifact_names[program_artifact]:
                        programs_seen[program.name] += 1

        network_dict = {"graph": [], "nodes": [], "links": [], "directed": False, "multigraph": False}

        node_indices = {}

        index = 0
        for artifact in all_artifacts:
            network_dict["nodes"].append({"size": artifacts_seen[artifact],
                                          "score": 0.5 if self.dataset.artifacts[artifact].provided_by_anvio else 1,
                                          "color": '#00AA00' if self.dataset.artifacts[artifact].provided_by_anvio else "#AA0000",
                                          "id": artifact,
                                          "name": artifact_names[artifact],
                                          "provided_by_anvio": True if self.dataset.artifacts[artifact].provided_by_anvio else False,
                                          "type": self.dataset.artifacts[artifact].type})
            node_indices[artifact] = index
            index += 1

        for program in programs.values():
            network_dict["nodes"].append({"size": programs_seen[program.name],
                                          "score": 0.1,
                                          "color": "#AAAA00",
                                          "id": program.name,
                                          "name": program.name,
                                          "type": "PROGRAM"})
            node_indices[program.name] = index
            index += 1

        for artifact in all_artifacts:
            for program in programs.values():
                for artifact_provided in program.provides:
                    if artifact_provided == artifact:
                        network_dict["links"].append({"source": node_indices[program.name], "target": node_indices[artifact], "type": "provides"})
                for artifact_needed in program.requires:
                    if artifact_needed == artifact:
                        network_dict["links"].append({"target": node_indices[program.name], "source": node_indices[artifact], "type": "requires"})
                for artifact_can_use in program.can_use:
                    if artifact_can_use == artifact:
                        network_dict["links"].append({"target": node_indices[program.name], "source": node_indices[artifact], "type": "can_use"})
                for artifact_can_provide in program.can_provide:
                    if artifact_can_provide == artifact:
                        network_dict["links"].append({"source": node_indices[program.name], "target": node_indices[artifact], "type": "can_provide"})

        return network_dict


    def render_body(self, record: DocumentRecord) -> str | None:
        if record.source_markdown is None:
            return None
        return self.render_anvio_markdown(record.source_markdown, record.documentation_path)


    def init_anvio_markdown_variables_conversion_dict(self):
        for program_name in self.program_names_and_paths:
            self.anvio_markdown_variables_conversion_dict[program_name] = """<span class="artifact-p">[%s](%s/programs/%s)</span>""" % (program_name, self.base_url, program_name)

        for artifact_name in self.dataset.artifacts:
            self.anvio_markdown_variables_conversion_dict[artifact_name] = """<span class="artifact-n">[%s](%s/artifacts/%s)</span>""" % (artifact_name, self.base_url, artifact_name)


    def render_anvio_markdown(self, content: str, source_path: str) -> str:
        """Renders markdown descriptions filling in anvi'o variables.

        Basically a lot of l_l83Я 1337 Я0XX0ЯZ stuff's going on down there, so you better run while you can.
        """

        markdown_content = content
        file_path = source_path

        if not len(self.anvio_markdown_variables_conversion_dict):
            self.init_anvio_markdown_variables_conversion_dict()

        # this is quite a big deal thing to do here:
        try:
            markdown_content = markdown_content % self.anvio_markdown_variables_conversion_dict
        except KeyError as e:
            self.progress.end()
            raise ConfigError("One of the variables, %s, in '%s' is not yet described anywhere :/ If it is not a typo but "
                              "a new artifact, you can add it to the file `anvio/programsdata.py`. After which everything "
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
            self.progress.end()
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


    def generate(self) -> None:
        """Validate the source bundle before preparing and writing the help output."""

        asset_sources = self.dataset.asset_sources
        source_paths = list(asset_sources.values())
        if self.dataset.dataset_path is not None:
            source_paths.append(self.dataset.dataset_path)
        for source in source_paths:
            if source.resolve() == self.output_directory or self.output_directory in source.resolve().parents:
                raise ConfigError(f"The output directory '{self.output_directory}' contains a source file needed to "
                                  f"render the documentation: '{source}'. Generating help pages replaces the output "
                                  f"directory, so please choose a separate output directory to keep the source intact.")

        try:
            for asset in self.dataset.assets:
                source = asset_sources.get(asset.path)
                if source is None or not source.is_file():
                    raise ConfigError(f"Anvi'o cannot find the documentation image '{asset.path}'. Please keep the "
                                      f"exported images directory beside documentation.json when moving a bundle, "
                                      f"or regenerate the documentation from its sources.")
                if utils.get_file_md5(source) != asset.md5:
                    raise ConfigError(f"The documentation image '{source}' no longer matches the checksum recorded "
                                      f"in the dataset. Please restore the image from the same export, or regenerate "
                                      f"the documentation so the dataset and its images agree.")

            filesnpaths.check_output_directory(self.output_directory, ok_if_exists=True)
            filesnpaths.gen_output_directory(self.output_directory, progress=self.progress, run=self.run,
                                            delete_if_exists=True, dont_warn=True)
            self.artifacts_output_dir = filesnpaths.gen_output_directory(self.output_directory / 'artifacts', progress=self.progress, run=self.run)
            self.programs_output_dir = filesnpaths.gen_output_directory(self.output_directory / 'programs', progress=self.progress, run=self.run)
            self.workflows_output_dir = filesnpaths.gen_output_directory(self.output_directory / 'workflows', progress=self.progress, run=self.run)

            self.copy_images()

            self.generate_pages_for_artifacts()

            self.generate_pages_for_programs()

            self.generate_pages_for_workflows()

            self.generate_index_page()

            self.dataset.write(self.output_directory / DATASET_FILENAME)
        except OSError as error:
            raise FilesNPathsError(f"Anvi'o could not finish generating documentation in '{self.output_directory}'. "
                                   f"Please check that the source images are readable and the output location has "
                                   f"enough disk space and write permissions. Here is the operating system error: "
                                   f"{error}") from error
        finally:
            self.progress.end()


    def copy_images(self):
        """Copies images from the documentation dataset to the output directory"""

        for asset in self.dataset.assets:
            destination = self.output_directory / asset.path
            filesnpaths.gen_output_directory(destination.parent, progress=self.progress, run=self.run)
            shutil.copyfile(self.dataset.asset_sources[asset.path], destination)


    def get_workflow_produced_artifacts_list(self, workflow_name, prefix="../../"):
        return [(r, '%sartifacts/%s' % (prefix, r)) for r in self.workflows[workflow_name].artifacts_produced]


    def get_workflow_accepted_artifacts_list(self, workflow_name, prefix="../../"):
        return [(r, '%sartifacts/%s' % (prefix, r)) for r in self.workflows[workflow_name].artifacts_accepted]


    def get_program_requires_provides_dict(self, prefix="../../"):
        d = {}

        for program_name in self.programs:
            d[program_name] = {}

            program = self.programs[program_name]
            d[program_name]['requires'] = [(r, '%sartifacts/%s' % (prefix, r)) for r in program.requires]
            d[program_name]['provides'] = [(r, '%sartifacts/%s' % (prefix, r)) for r in program.provides]
            d[program_name]['can_use'] = [(r, '%sartifacts/%s' % (prefix, r)) for r in program.can_use]
            d[program_name]['can_provide'] = [(r, '%sartifacts/%s' % (prefix, r)) for r in program.can_provide]
            d[program_name]['anvio_workflows'] = [(w, '%sworkflows/%s' % (prefix, w)) for w in program.anvio_workflows]

        return d


    def generate_pages_for_artifacts(self):
        """Generates static pages for artifacts in the output directory"""

        self.progress.new("Rendering artifact pages", progress_total_items=len(self.dataset.artifacts))
        self.progress.update('...')

        for artifact in self.dataset.artifacts:
            self.progress.update(f"'{artifact}' ...", increment=True)

            d = {'artifact': asdict(self.dataset.artifacts[artifact]),
                 'meta': {'summary_type': 'artifact',
                          'version': '\n'.join(['|%s|%s|' % (t[0], t[1]) for t in self.dataset.meta.versions]),
                          'date': self.dataset.meta.date,
                          'version_short_identifier': self.version_short_identifier}
                }

            d['artifact']['name'] = artifact
            d['artifact']['required_by'] = [(r, '../../programs/%s' % r) for r in self.artifacts_info[artifact]['required_by']]
            d['artifact']['provided_by'] = [(r, '../../programs/%s' % r) for r in self.artifacts_info[artifact]['provided_by']]
            d['artifact']['can_used_by'] = [(r, '../../programs/%s' % r) for r in self.artifacts_info[artifact]['can_used_by']]
            d['artifact']['can_provided_by'] = [(r, '../../programs/%s' % r) for r in self.artifacts_info[artifact]['can_provided_by']]
            d['artifact']['description'] = self.artifacts_info[artifact]['description']
            d['artifact']['icon'] = '../../images/icons/%s.png' % self.dataset.artifacts[artifact].type

            if anvio.DEBUG:
                self.progress.reset()
                self.run.warning(None, 'THE OUTPUT DICT')
                import json
                print(json.dumps(d, indent=2))

            self.progress.update(f"'{artifact}' ... rendering ...", increment=False)
            artifact_output_dir = filesnpaths.gen_output_directory(os.path.join(self.artifacts_output_dir, artifact), progress=self.progress, run=self.run)
            output_file_path = os.path.join(artifact_output_dir, 'index.md')
            open(output_file_path, 'w', encoding='utf-8').write(SummaryHTMLOutput(d, r=self.run, p=self.progress).render())

        self.progress.end()


    def get_HTML_formatted_authors_data(self, authors):
        """for a given program, returns HTML-formatted authors data"""

        d = ""

        for author in authors:
            d += '''<div class="anvio-person"><div class="anvio-person-info">'''
            d += f'''<div class="anvio-person-photo"><img class="anvio-person-photo-img" src="../../images/authors/{os.path.basename(self.authors[author].avatar)}" /></div>'''
            d += '''<div class="anvio-person-info-box">'''
            d += f'''<a href="/people/{self.authors[author].github}" target="_blank"><span class="anvio-person-name">{self.authors[author].name}</span></a>'''
            d += '''<div class="anvio-person-social-box">'''

            if self.authors[author].web is not None:
                d += f'''<a href="{self.authors[author].web}" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a>'''

            d += f'''<a href="mailto:{self.authors[author].email}" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a>'''

            if self.authors[author].twitter is not None:
                d += f'''<a href="http://twitter.com/{self.authors[author].twitter}" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a>'''

            d += f'''<a href="http://github.com/{self.authors[author].github}" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a>'''

            d += '''</div></div></div></div>\n\n'''

        return d


    def get_HTML_formatted_authors_data_mini(self, authors):
        """for a given list of authors, returns a tiny version of the HTML-formatted authors data"""

        d = ""

        for author in authors:
            d += '''<div class="anvio-person-mini"><div class="anvio-person-photo-mini">'''
            d += f'''<a href="/people/{self.authors[author].github}" target="_blank"><img class="anvio-person-photo-img-mini" title="{self.authors[author].name}" src="images/authors/{os.path.basename(self.authors[author].avatar)}" /></a>'''
            d += '''</div></div>\n'''

        return d


    def get_HTML_formatted_third_party_programs(self, workflow_name):
        """Get a template-friendly list of third-party programs used from within a workflow"""

        d = []

        for purpose, program_names in self.workflows[workflow_name].third_party_programs_used:
            for program_name in program_names:
                d.append(f'''<a href="{self.dataset.third_party_programs[program_name].link}" target="_blank">{program_name}</a> ({purpose})''')

        return d


    def generate_pages_for_workflows(self):
        """Generate static pages for anvi'o workflows in the output directory"""

        self.progress.new("Rendering workflow pages", progress_total_items=len(self.workflows))
        self.progress.update('...')

        for workflow_name in self.workflows:
            self.progress.update(f"'{workflow_name}' ...", increment=True)

            d = {'workflow': asdict(self.workflows[workflow_name]),
                 'meta': {'summary_type': 'workflow',
                          'version': '\n'.join(['|%s|%s|' % (t[0], t[1]) for t in self.dataset.meta.versions]),
                          'date': self.dataset.meta.date,
                          'version_short_identifier': self.version_short_identifier}
                 }

            d['workflow']['description'] = self.bodies[self.workflows[workflow_name].id]
            d['workflow']['artifacts_produced'] = self.get_workflow_produced_artifacts_list(workflow_name)
            d['workflow']['artifacts_accepted'] = self.get_workflow_accepted_artifacts_list(workflow_name)
            d['workflow']['third_party_programs_used'] = self.get_HTML_formatted_third_party_programs(workflow_name)
            d['workflow']['authors'] = self.get_HTML_formatted_authors_data(d['workflow']['authors'])

            # also add information regarding the artifacts
            d['artifacts'] = self.artifacts_info

            if anvio.DEBUG:
                self.progress.reset()
                self.run.warning(None, 'THE WORKFLOW OUTPUT DICT')
                import json
                print(json.dumps(d, indent=2))

            self.progress.update(f"'{workflow_name}' ... rendering ...", increment=False)
            workflow_output_dir = filesnpaths.gen_output_directory(os.path.join(self.workflows_output_dir, workflow_name), progress=self.progress, run=self.run)
            output_file_path = os.path.join(workflow_output_dir, 'index.md')
            open(output_file_path, 'w', encoding='utf-8').write(SummaryHTMLOutput(d, r=self.run, p=self.progress).render())

        self.progress.end()


    def generate_pages_for_programs(self):
        """Generates static pages for programs in the output directory"""

        self.progress.new("Rendering program pages", progress_total_items=len(self.programs))
        self.progress.update('...')

        program_provides_requires_dict = self.get_program_requires_provides_dict()

        resources_example_program = 'anvi-interactive'
        resources_example_path = None
        if resources_example_program in self.program_names_and_paths:
            resources_example_path = self.program_names_and_paths[resources_example_program]
            resources_example_path = resources_example_path.replace(os.sep, '/')

        for program_name in self.programs:
            self.progress.update(f"'{program_name}' ...", increment=True)

            program = self.programs[program_name]
            program_source_path = program.source_path
            d = {'program': {},
                 'meta': {'summary_type': 'program',
                          'version': '\n'.join(['|%s|%s|' % (t[0], t[1]) for t in self.dataset.meta.versions]),
                          'date': self.dataset.meta.date,
                          'anvio_version': self.dataset.meta.version,
                          'version_short_identifier': self.version_short_identifier}
                }

            d['program']['name'] = program_name
            d['program']['usage'] = self.bodies[program.id]
            d['program']['description'] = program.description
            d['program']['resources'] = program.resources
            d['program']['source_path'] = program_source_path
            d['program']['resources_example_source_path'] = resources_example_path or program_source_path
            d['program']['requires'] = program_provides_requires_dict[program_name]['requires']
            d['program']['provides'] = program_provides_requires_dict[program_name]['provides']
            d['program']['can_use'] = program_provides_requires_dict[program_name]['can_use']
            d['program']['can_provide'] = program_provides_requires_dict[program_name]['can_provide']
            d['program']['icon'] = '../../images/icons/%s.png' % 'PROGRAM'
            d['program']['authors'] = self.get_HTML_formatted_authors_data(program.authors)
            d['program']['command_line_usage'] = program.command_line_usage
            d['program']['parameter_groups'] = [asdict(group) for group in program.parameter_groups]
            d['artifacts'] = self.artifacts_info
            d['workflows'] = {name: asdict(workflow) for name, workflow in self.workflows.items()}

            if anvio.DEBUG:
                self.progress.reset()
                self.run.warning(None, 'THE OUTPUT DICT')
                print(json.dumps(d, indent=2))

            self.progress.update(f"'{program_name}' ... rendering ...", increment=False)
            program_output_dir = filesnpaths.gen_output_directory(os.path.join(self.programs_output_dir, program_name), progress=self.progress, run=self.run)
            output_file_path = os.path.join(program_output_dir, 'index.md')
            open(output_file_path, 'w', encoding='utf-8').write(SummaryHTMLOutput(d, r=self.run, p=self.progress).render())

            # create the program network, too
            self.progress.update(f"'{program_name}' ... rendering ... network json ...", increment=False)
            program_output_dir = filesnpaths.gen_output_directory(os.path.join(self.programs_output_dir, program_name), progress=self.progress, run=self.run)
            program_network = self.programs_network([program_name], artifact_names_as_ids=True)
            with open(os.path.join(program_output_dir, "network.json"), "w", encoding="utf-8") as output:
                json.dump(program_network, output, indent=2)

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
        workflows = {name: asdict(workflow) for name, workflow in self.workflows.items()}
        for workflow in workflows:
            workflows[workflow]['authors'] = self.get_HTML_formatted_authors_data_mini(self.workflows[workflow].authors)

        # please note that artifacts get a fancy dictionary with everything, while programs get a crappy tuples list.
        # if we need to improve the functionality of the help index page, we may need to update programs
        # to a fancy dictionary, too.
        d = {'programs': [(p, 'programs/%s' % p, self.programs[p].description, self.get_HTML_formatted_authors_data_mini(self.programs[p].authors)) for p in self.programs],
             'workflows': workflows,
             'artifacts': self.artifacts_info,
             'artifact_types': self.artifact_types,
             'meta': {'summary_type': 'programs_and_artifacts_index',
                      'version': self.dataset.meta.version,
                      'date': self.dataset.meta.date}
            }

        d['program_provides_requires'] = self.get_program_requires_provides_dict(prefix='')

        self.progress.update('Rendering...')
        output_file_path = os.path.join(self.output_directory, 'index.md')

        self.progress.update('Writing...')
        open(output_file_path, 'w', encoding='utf-8').write(SummaryHTMLOutput(d, r=self.run, p=self.progress).render())

        self.progress.end()


class LLMsTextRenderer:
    """Flatten a documentation dataset into llms-full.txt.

    The output follows the llms-full.txt convention: every program, artifact, and workflow
    becomes one self-contained section with a predictable heading, its metadata, and its
    authored Markdown. Website-specific markup is converted to plain Markdown so the file
    can be read directly or split into chunks for retrieval.
    """

    def __init__(self, dataset: DocumentationData, output_file: str = "llms-full.txt", r=run, p=progress):
        """Prepare lookups shared by all sections.

        Args:
            dataset: The documentation dataset to flatten.
            output_file: Destination of llms-full.txt.
            r: Reports the output location.
            p: Reports rendering progress.
        """
        self.dataset = dataset
        self.output_file = Path(output_file)
        self.run = r
        self.progress = p
        self.site_url = f"https://anvio.org{dataset.meta.base_url}"

        self.body_top_heading_level = 3
        self.max_heading_level = 6
        self._artifact_relations = {"requires": "Required by", "provides": "Provided by",
                                    "can_use": "Can be used by", "can_provide": "Can be provided by"}


    def render(self) -> str:
        """Render the full document.

        Returns:
            The Markdown text of llms-full.txt.
        """
        renderers = [(self.dataset.workflows, self.render_workflow),
                     (self.dataset.programs, self.render_program),
                     (self.dataset.artifacts, self.render_artifact)]
        total = sum(len(records) for records, _ in renderers)

        self.progress.new("Rendering llms-full.txt", progress_total_items=total)
        sections = [self.render_header()]
        for records, render_record in renderers:
            for name, record in sorted(records.items()):
                self.progress.update(f"{record.kind} '{name}' ...", increment=True)
                sections.append(render_record(record))
        self.progress.end()

        return "\n\n".join(sections) + "\n"


    def write(self) -> None:
        """Render the document and write it to the output file."""
        filesnpaths.is_output_file_writable(self.output_file)
        self.output_file.write_text(self.render(), encoding="utf-8")
        self.run.info("llms-full.txt", self.output_file)


    def render_header(self) -> str:
        meta = self.dataset.meta
        versions = "\n".join(f"- {label}: {value}" for label, value in meta.versions)

        return (f"# anvi'o documentation\n\n"
                f"> anvi'o is an open-source platform for microbial 'omics. This file contains the help pages "
                f"for every anvi'o workflow, program, and artifact (the files and concepts programs use and "
                f"produce). Each entry is a level-two section titled `<name> (<kind>)`.\n\n"
                f"Generated {meta.date} from anvi'o {meta.version}. Online version: {self.site_url}/\n\n"
                f"{versions}")


    def render_program(self, record: ProgramDocumentation) -> str:
        fields = [("URL", record.html_url),
                  ("Description", record.description),
                  ("Requires", self.format_names(record.requires)),
                  ("Can use", self.format_names(record.can_use)),
                  ("Provides", self.format_names(record.provides)),
                  ("Can provide", self.format_names(record.can_provide)),
                  ("Workflows", self.format_names(record.anvio_workflows)),
                  ("Tags", ", ".join(record.tags)),
                  ("Authors", self.format_authors(record.authors))]

        resources = "\n".join(f"- [{title}]({url})" for title, url in record.resources)
        if resources:
            resources = f"Resources:\n\n{resources}"

        return self.render_section(record.name, record.kind, fields, [self.render_body(record), resources])


    def render_artifact(self, record: ArtifactDocumentation) -> str:
        fields = [("URL", record.html_url),
                  ("Display name", record.name),
                  ("Type", record.type),
                  ("Provided by anvi'o", "yes" if record.provided_by_anvio else "no"),
                  ("Provided by the user", "yes" if record.provided_by_user else "no")]

        artifact_id = record.id.removeprefix("artifact:")
        for relation, label in self._artifact_relations.items():
            programs = [name for name, program in sorted(self.dataset.programs.items())
                        if artifact_id in getattr(program, relation)]
            fields.append((label, self.format_names(programs)))

        return self.render_section(artifact_id, record.kind, fields, [self.render_body(record)])


    def render_workflow(self, record: WorkflowDocumentation) -> str:
        third_party = "; ".join(f"{purpose}: {', '.join(programs)}" for purpose, programs in record.third_party_programs_used)

        fields = [("URL", record.html_url),
                  ("Summary", record.one_sentence_summary),
                  ("Accepts", self.format_names(record.artifacts_accepted)),
                  ("Produces", self.format_names(record.artifacts_produced)),
                  ("anvi'o programs used", self.format_names(record.anvio_programs_used)),
                  ("Third-party programs used", third_party),
                  ("Inherits workflows", self.format_names(record.anvio_workflows_inherited)),
                  ("Authors", self.format_authors(record.authors))]

        return self.render_section(record.name, record.kind, fields, [record.one_paragraph_summary, self.render_body(record)])


    def render_section(self, name: str, kind: str, fields: list[tuple[str, str]], blocks: list[str]) -> str:
        metadata = "\n".join(f"- {label}: {value}" for label, value in fields if value)
        parts = [f"## {name} ({kind})", metadata] + [block for block in blocks if block]
        return "\n\n".join(parts)


    def render_body(self, record: ProgramDocumentation | ArtifactDocumentation | WorkflowDocumentation) -> str:
        """Nest the authored Markdown under its section heading.

        The Markdown is kept as authored except for headings, which are shifted so the
        shallowest heading in the body sits directly below the section heading. This keeps
        the main sections of every entry at the same level regardless of how the source
        Markdown was authored. Lines inside code blocks are never treated as headings.

        Args:
            record: The documentation record whose authored Markdown is nested.

        Returns:
            The nested Markdown, or an empty string if the record has no authored Markdown.
        """
        if record.source_markdown is None:
            return ""

        lines = record.source_markdown.strip().split("\n")
        heading_indices = []
        in_code = False
        for index, line in enumerate(lines):
            stripped = line.strip()
            if stripped in ("{{ codestart }}", "{{ codestop }}"):
                in_code = stripped == "{{ codestart }}"
            elif stripped.startswith("```"):
                in_code = not in_code
            elif not in_code and re.match(r"#{1,6} ", line):
                heading_indices.append(index)

        if heading_indices:
            top_level = min(self.heading_level(lines[index]) for index in heading_indices)
            offset = self.body_top_heading_level - top_level
            for index in heading_indices:
                lines[index] = self.shift_heading(lines[index], offset)

        return "\n".join(lines)


    def heading_level(self, line: str) -> int:
        return len(line) - len(line.lstrip("#"))


    def shift_heading(self, line: str, offset: int) -> str:
        level = self.heading_level(line)
        return "#" * min(level + offset, self.max_heading_level) + line[level:]


    def format_names(self, names: list[str]) -> str:
        return ", ".join(f"`{name}`" for name in names)


    def format_authors(self, authors: list[str]) -> str:
        return ", ".join(f"{self.dataset.authors[author].name} (@{author})" for author in authors)
