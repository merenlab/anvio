"""Render the existing anvi'o help pages."""

import os
import argparse

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.programs import AnvioPrograms, AnvioArtifacts, AnvioWorkflows, ProgramsNetwork, run, progress
from anvio.programsdata import ANVIO_ARTIFACTS, ANVIO_WORKFLOWS, THIRD_PARTY_PROGRAMS
from anvio.summaryhtml import SummaryHTMLOutput


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
        self.repo_root = os.path.abspath(os.path.join(os.path.dirname(anvio.__file__), '..'))

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
                              "Regardless. These are the artifact types missing images: %s."
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
            d[program_name]['can_use'] = [(r.id, '%sartifacts/%s' % (prefix, r.id)) for r in program.meta_info['can_use']['value']]
            d[program_name]['can_provide'] = [(r.id, '%sartifacts/%s' % (prefix, r.id)) for r in program.meta_info['can_provide']['value']]
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
            d['artifact']['can_used_by'] = [(r, '../../programs/%s' % r) for r in self.artifacts_info[artifact]['can_used_by']]
            d['artifact']['can_provided_by'] = [(r, '../../programs/%s' % r) for r in self.artifacts_info[artifact]['can_provided_by']]
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

        resources_example_program = 'anvi-interactive'
        resources_example_path = None
        if resources_example_program in self.program_names_and_paths:
            resources_example_path = os.path.relpath(self.program_names_and_paths[resources_example_program], self.repo_root)
            resources_example_path = resources_example_path.replace(os.sep, '/')

        for program_name in self.programs:
            self.progress.update(f"'{program_name}' ...", increment=True)

            program = self.programs[program_name]
            program_source_path = os.path.relpath(program.program_path, self.repo_root).replace(os.sep, '/')
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
            d['program']['source_path'] = program_source_path
            d['program']['resources_example_source_path'] = resources_example_path or program_source_path
            d['program']['requires'] = program_provides_requires_dict[program_name]['requires']
            d['program']['provides'] = program_provides_requires_dict[program_name]['provides']
            d['program']['can_use'] = program_provides_requires_dict[program_name]['can_use']
            d['program']['can_provide'] = program_provides_requires_dict[program_name]['can_provide']
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
