"""Collect documentation from CLI metadata, source registries, Markdown, and images."""

import argparse
import os
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
    ProgramDocumentation,
    ThirdPartyProgram,
    WorkflowDocumentation,
)
from anvio.errors import ConfigError
from anvio.programs import AnvioPrograms, AnvioArtifacts, AnvioWorkflows


class AnvioDocs(AnvioPrograms, AnvioArtifacts, AnvioWorkflows):
    """Collect source documentation and metadata into a serializable dataset."""

    def __init__(self, args: argparse.Namespace, r: terminal.Run = terminal.Run(),
                 p: terminal.Progress = terminal.Progress()) -> None:
        self.args = args
        self.run = r
        self.progress = p

        self.repo_root = os.path.abspath(os.path.join(os.path.dirname(anvio.__file__), '..'))

        if not os.path.exists(anvio.DOCS_PATH):
            raise ConfigError("The anvi'o docs path is not where it should be :/ Something funny is going on.")

        self.version_short_identifier = 'm' if anvio.anvio_version_for_help_docs == 'main' else anvio.anvio_version_for_help_docs
        self.base_url = f'/help/{anvio.anvio_version_for_help_docs}'

        AnvioPrograms.__init__(self, args, r=self.run, p=self.progress)
        self.init_programs()

        AnvioArtifacts.__init__(self, args, r=self.run, p=self.progress)
        self.init_artifacts()

        AnvioWorkflows.__init__(self, args, r=self.run, p=self.progress)
        self.init_workflows()

        if not len(self.programs):
            raise ConfigError("AnvioDocs is asked to process the usage statements of some programs, but the "
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


    def read_anvio_markdown(self, file_path: str) -> str:
        """Collect authored Markdown; website formatting belongs to the renderer."""
        filesnpaths.is_file_plain_text(file_path)
        with open(file_path, encoding='utf-8') as source:
            return source.read()


    def collect(self) -> DocumentationData:
        """Capture help records without presentation changes or inferred provenance."""
        programs: dict[str, ProgramDocumentation] = {}
        for name, program in self.programs.items():
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
            )
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
