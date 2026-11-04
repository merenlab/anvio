"""Collect documentation from CLI metadata, source registries, Markdown, and images."""

import os
import sys
import argparse
from pathlib import Path

import anvio
import anvio.filesnpaths as filesnpaths
import anvio.terminal as terminal
import anvio.utils as utils

from anvio.programsdata import ANVIO_ARTIFACTS, THIRD_PARTY_PROGRAMS
from anvio.docs.data import (
    ArtifactDocumentation,
    AuthorDocumentation,
    BuildMetadata,
    DocumentationAsset,
    DocumentationData,
    ParameterGroup,
    ProgramDocumentation,
    ProgramParameter,
    ThirdPartyProgram,
    WorkflowDocumentation,
)
from anvio.errors import ConfigError
from anvio.argparse import get_parser_for_program
from anvio.programs import AnvioPrograms, AnvioArtifacts, AnvioWorkflows


class AnvioDocs(AnvioPrograms, AnvioArtifacts, AnvioWorkflows):
    """Collect source documentation and metadata into a serializable dataset."""


    def __init__(self, args, r=terminal.Run(), p=terminal.Progress()):
        self.args = args
        self.run = r
        self.progress = p

        self.repo_root = os.path.abspath(os.path.join(os.path.dirname(anvio.__file__), '..'))

        if not os.path.exists(anvio.DOCS_PATH):
            raise ConfigError("The anvi'o docs path is not where it should be :/ Something funny is going on.")

        self.version_short_identifier = 'm' if anvio.anvio_version_for_help_docs == 'main' else anvio.anvio_version_for_help_docs
        self.base_url = os.path.join("/help", anvio.anvio_version_for_help_docs)

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


    def read_anvio_markdown(self, file_path):
        """Collect authored Markdown; website formatting belongs to the renderer."""
        filesnpaths.is_file_plain_text(file_path)
        with open(file_path, encoding='utf-8') as source:
            return source.read()


    def get_JSON_value(self, value):
        """Turns parameter defaults and choices into values JSON can hold and read back as they were"""

        if value is None or isinstance(value, (bool, int, float, str)):
            return value
        elif isinstance(value, (set, frozenset)):
            return sorted([self.get_JSON_value(v) for v in value], key=str)
        elif isinstance(value, (list, tuple, range, type({}.keys()))):
            return [self.get_JSON_value(v) for v in value]
        else:
            return str(value)


    def get_command_line_parameters(self, program_name):
        """Returns the usage line and the parameter groups of a program, as its argument parser describes them.

        The parser comes from `get_parser_for_program`, which learns it without running the program.
        Programs without a parser we can reach (such as the non-Python scripts, or retired programs
        that refuse to run) get no usage line and an empty list of parameter groups.
        """

        entry_point = self.program_entry_points.get(program_name)
        parser = get_parser_for_program(program_name, entry_point) if entry_point else None

        if not parser:
            return None, []

        # the width is set to the maximum so the usage line does not depend on the terminal of the
        # person who happens to generate the documentation
        formatter = argparse.HelpFormatter(parser.prog, width=sys.maxsize)
        formatter.add_usage(parser.usage, parser._actions, parser._mutually_exclusive_groups, prefix='')
        command_line_usage = formatter.format_help().strip()

        # programs with modes get one more group for each mode, titled with the flag and the mode that
        # enable its parameters, since the groups of mode parsers carry titles formatted for the terminal
        groups = [(group, group.title, None) for group in parser._action_groups]
        for mode, mode_parser in parser.mode_parsers.items():
            groups += [(group, f"{parser.mode_flag} {mode}", mode) for group in mode_parser._action_groups]

        parameter_groups = []
        for group, title, mode in groups:
            parameters = []
            for action in group._group_actions:
                # the help flag is the same for every program, and the parameters whose help is
                # suppressed are hidden from users on purpose
                if isinstance(action, (argparse._HelpAction, argparse._VersionAction)) or action.help == argparse.SUPPRESS:
                    continue

                # this fills in things like `%(default)s` the same way `--help` would
                try:
                    help_text = formatter._expand_help(action) if action.help else ''
                except (KeyError, TypeError, ValueError):
                    help_text = action.help

                metavar = ' '.join(action.metavar) if isinstance(action.metavar, tuple) else action.metavar

                parameters.append(ProgramParameter(flags=list(action.option_strings),
                                                   dest=action.dest,
                                                   help=' '.join(help_text.split()),
                                                   metavar=metavar,
                                                   default=None if action.default is argparse.SUPPRESS else self.get_JSON_value(action.default),
                                                   choices=None if action.choices is None else self.get_JSON_value(action.choices),
                                                   required=action.required,
                                                   takes_value=action.nargs != 0,
                                                   nargs=action.nargs))

            if parameters:
                parameter_groups.append(ParameterGroup(title=title,
                                                       description=' '.join(group.description.split()) if group.description else None,
                                                       parameters=parameters,
                                                       mode=mode))

        return command_line_usage, parameter_groups


    def collect(self) -> DocumentationData:
        """Capture help records without presentation changes or inferred provenance."""
        programs: dict[str, ProgramDocumentation] = {}
        programs_without_parameters = []
        self.progress.new('Learning program parameters', progress_total_items=len(self.programs))
        for name, program in self.programs.items():
            self.progress.update(name, increment=True)
            command_line_usage, parameter_groups = self.get_command_line_parameters(name)
            if command_line_usage is None:
                programs_without_parameters.append(name)

            metadata = program.meta_info
            programs[name] = ProgramDocumentation(
                id=f'program:{name}', kind='program', name=name,
                source_markdown=program.usage,
                documentation_path=f'anvio/docs/programs/{name}.md',
                html_url=f'https://anvio.org{self.base_url}/programs/{name}/',
                source_path=Path(program.program_path).relative_to(self.repo_root).as_posix(),
                description=metadata['description']['value'],
                authors=list(metadata['authors']['value']),
                resources=[(title, url) for title, url in metadata['resources']['value']],
                tags=list(metadata['tags']['value']),
                anvio_workflows=list(metadata['anvio_workflows']['value']),
                requires=[artifact.id for artifact in metadata['requires']['value']],
                provides=[artifact.id for artifact in metadata['provides']['value']],
                can_use=[artifact.id for artifact in metadata['can_use']['value']],
                can_provide=[artifact.id for artifact in metadata['can_provide']['value']],
                command_line_usage=command_line_usage,
                parameter_groups=parameter_groups,
            )
        self.progress.end()

        if programs_without_parameters:
            self.run.info_single(f"Anvi'o could not learn the command line parameters of {len(programs_without_parameters)} "
                                 f"of {len(self.programs)} programs, so their help pages will not list any. These are "
                                 f"usually programs that are not written in Python, or retired programs that only tell "
                                 f"their users what to use instead: {', '.join(sorted(programs_without_parameters))}.",
                                 nl_after=1, nl_before=1)

        artifacts = {
            name: ArtifactDocumentation(
                id=f'artifact:{name}', kind='artifact', name=artifact['name'],
                source_markdown=self.artifacts_info[name]['description'],
                documentation_path=f'anvio/docs/artifacts/{name}.md',
                html_url=f'https://anvio.org{self.base_url}/artifacts/{name}/',
                type=artifact['type'], provided_by_anvio=artifact['provided_by_anvio'],
                provided_by_user=artifact['provided_by_user'],
            ) for name, artifact in ANVIO_ARTIFACTS.items()
        }
        workflows = {
            name: WorkflowDocumentation(
                id=f'workflow:{name}', kind='workflow', name=name,
                source_markdown=workflow.get('description'),
                documentation_path=f'anvio/docs/workflows/{name}.md',
                html_url=f'https://anvio.org{self.base_url}/workflows/{name}/',
                authors=list(workflow['authors']),
                artifacts_produced=list(workflow['artifacts_produced']),
                artifacts_accepted=list(workflow['artifacts_accepted']),
                anvio_workflows_inherited=list(workflow['anvio_workflows_inherited']),
                third_party_programs_used=[(purpose, list(names)) for purpose, names in workflow['third_party_programs_used']],
                one_sentence_summary=workflow['one_sentence_summary'],
                one_paragraph_summary=workflow['one_paragraph_summary'],
                anvio_programs_used=list(workflow['anvio_programs_used']),
            ) for name, workflow in self.workflows.items()
        }
        images = Path(self.images_source_directory)
        asset_sources = {f'images/{image.relative_to(images).as_posix()}': image
                         for image in sorted(images.rglob('*')) if image.is_file()}
        authors: dict[str, AuthorDocumentation] = {}
        for name, entry in self.authors.items():
            source = Path(entry['avatar'])
            destination = f'images/authors/{source.name}'
            authors[name] = AuthorDocumentation(
                github=entry['github'], name=entry['name'], email=entry['email'],
                avatar=destination, web=entry.get('web'), twitter=entry.get('twitter'),
            )
            asset_sources[destination] = source
        metadata = BuildMetadata(
            docs_version=anvio.anvio_version_for_help_docs,
            version=f'{anvio.anvio_version} ({anvio.anvio_codename})',
            versions=anvio.get_version_tuples(),
            date=utils.get_date(),
            base_url=self.base_url,
            version_short_identifier=self.version_short_identifier,
        )
        dataset = DocumentationData(
            meta=metadata,
            program_source_paths={name: Path(path).relative_to(self.repo_root).as_posix()
                                  for name, path in self.program_names_and_paths.items()},
            programs=programs, artifacts=artifacts, workflows=workflows, authors=authors,
            third_party_programs={name: ThirdPartyProgram(link=record['link']) for name, record in THIRD_PARTY_PROGRAMS.items()},
            assets=[DocumentationAsset(path=path, md5=utils.get_file_md5(source)) for path, source in asset_sources.items()],
        )
        dataset.asset_sources = asset_sources
        return dataset
