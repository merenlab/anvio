"""Render the help website from a documentation dataset.

Rendering uses only the dataset, its image assets, and installed website templates.
It does not discover CLI programs or consult documentation source registries.
"""

from __future__ import annotations

import json
import shutil
from collections import Counter
from pathlib import Path
from typing import Any

import anvio
import anvio.filesnpaths as filesnpaths
import anvio.terminal as terminal
import anvio.utils as utils

from anvio.docs.data import (
    DATASET_FILENAME, RELATIONSHIPS, DocumentationData, DocumentRecord,
    WorkflowDocumentation,
)
from anvio.errors import ConfigError, FilesNPathsError
from anvio.summaryhtml import SummaryHTMLOutput


ProgramLinks = dict[str, list[tuple[str, str]]]
run = terminal.Run()
progress = terminal.Progress()


class HelpPagesRenderer:
    """Render the Jekyll help pages from typed documentation records."""

    def __init__(self, dataset: DocumentationData, output_directory: str, r=run, p=progress):
        self.dataset = dataset
        self.output_directory = Path(output_directory).resolve()
        self.run = r
        self.progress = p
        self.markdown_links = {
            name: f'<span class="artifact-p">[{name}]({dataset.meta.base_url}/programs/{name})</span>'
            for name in dataset.program_source_paths
        }
        self.markdown_links.update(
            {
                name: f'<span class="artifact-n">[{name}]({dataset.meta.base_url}/artifacts/{name})</span>'
                for name in dataset.artifacts
            }
        )
        # Resolve every body before replacing any existing output.
        self.bodies = {
            record.id: self.render_body(record)
            for records in (dataset.programs, dataset.artifacts, dataset.workflows)
            for record in records.values()
        }
        # These dictionaries are adapters for the existing template lookup filter.
        self.artifact_context = {
            name: {"type": artifact.type, "path": f"artifacts/{name}"}
            for name, artifact in dataset.artifacts.items()
        }

    def programs_network(
        self,
        program_names: list[str] | None = None,
        artifact_names_as_ids: bool = False,
    ) -> dict[str, Any]:
        """Format dataset relationships for the website visualization."""
        programs = {
            name: record
            for name, record in self.dataset.programs.items()
            if (program_names is None or name in program_names)
            and (record.requires or record.provides)
        }

        artifact_order: dict[str, int] = {}
        for program in programs.values():
            for relation in ("provides", "requires", "can_use", "can_provide"):
                for artifact in getattr(program, relation):
                    artifact_order[artifact] = artifact_order.get(artifact, 0) + 1

        network: dict[str, Any] = {
            "graph": [],
            "nodes": [],
            "links": [],
            "directed": False,
            "multigraph": False,
        }

        indices: dict[str, int] = {}
        for artifact, count in artifact_order.items():
            record = self.dataset.artifacts[artifact]
            indices[artifact] = len(network["nodes"])
            network["nodes"].append(
                {
                    "size": count,
                    "score": 0.5 if record.provided_by_anvio else 1,
                    "color": "#00AA00" if record.provided_by_anvio else "#AA0000",
                    "id": artifact,
                    "name": artifact if artifact_names_as_ids else record.name,
                    "provided_by_anvio": record.provided_by_anvio,
                    "type": record.type,
                }
            )

        display_names = {
            artifact: artifact
            if artifact_names_as_ids
            else self.dataset.artifacts[artifact].name
            for artifact in artifact_order
        }

        name_counts = Counter(display_names.values())
        for name, program in programs.items():
            indices[name] = len(network["nodes"])
            # The legacy visualization counts matching display names, including aliases.
            network["nodes"].append(
                {
                    "size": sum(
                        name_counts[display_names[artifact]]
                        for relation in RELATIONSHIPS
                        for artifact in getattr(program, relation)
                    ),
                    "score": 0.1,
                    "color": "#AAAA00",
                    "id": name,
                    "name": name,
                    "type": "PROGRAM",
                }
            )

        for artifact in artifact_order:
            for name, program in programs.items():
                for relation in ("provides", "requires", "can_use", "can_provide"):
                    for referenced_artifact in getattr(program, relation):
                        if referenced_artifact != artifact:
                            continue
                        if relation in ("provides", "can_provide"):
                            edge = {
                                "source": indices[name],
                                "target": indices[artifact],
                                "type": relation,
                            }
                        else:
                            edge = {
                                "target": indices[name],
                                "source": indices[artifact],
                                "type": relation,
                            }
                        network["links"].append(edge)

        return network

    def render_body(self, record: DocumentRecord) -> str | None:
        if record.source_markdown is None:
            return None
        return self.render_anvio_markdown(
            record.source_markdown, record.documentation_path
        )

    def render_anvio_markdown(self, content: str, source_path: str) -> str:
        """Resolve anvi'o references and preserve the website's code-block markup."""
        try:
            # This is the authored documentation's substitution syntax, not Python message formatting.
            content = content % self.markdown_links
        except (KeyError, TypeError, ValueError) as error:
            raise ConfigError(
                f"Could not resolve documentation references in '{source_path}': {error}. "
                "Use %(name)s for known programs/artifacts and %% for literal percent signs."
            ) from error
        lines = content.split("\n")
        starts = [
            number
            for number, line in enumerate(lines)
            if line.strip() == "{{ codestart }}"
        ]
        stops = [
            number
            for number, line in enumerate(lines)
            if line.strip() == "{{ codestop }}"
        ]
        if len(starts) != len(stops):
            raise ConfigError(f"Unmatched code-block markers in '{source_path}'.")
        for start, stop in zip(starts, stops):
            for number in range(start + 1, stop):
                lines[number] = (
                    lines[number]
                    .replace("-", "&#45;")
                    .replace("*", "&#42;")
                    .replace("==", "&#61;&#61;")
                )
        return (
            "\n".join(lines)
            .replace("{{ codestart }}", '<div class="codeblock" markdown="1">')
            .replace("{{ codestop }}", "</div>")
        )

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
            for directory in ("artifacts", "programs", "workflows"):
                filesnpaths.gen_output_directory(self.output_directory / directory, progress=self.progress, run=self.run)
            for asset in self.dataset.assets:
                destination = self.output_directory / asset.path
                filesnpaths.gen_output_directory(destination.parent, progress=self.progress, run=self.run)
                shutil.copyfile(asset_sources[asset.path], destination)

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

    def render_page(
        self, relative_path: str, summary_type: str, context: dict[str, Any]
    ) -> None:
        metadata = self.dataset.meta
        context = {
            **context,
            "meta": {
                "summary_type": summary_type,
                "version": metadata.version
                if summary_type == "programs_and_artifacts_index"
                else "\n".join(f"|{name}|{version}|" for name, version in metadata.versions),
                "date": metadata.date,
                "version_short_identifier": metadata.version_short_identifier,
            },
        }
        if anvio.DEBUG:
            self.progress.reset()
            self.run.warning(None, "THE OUTPUT DICT")
            print(json.dumps(context, indent=2))
        path = self.output_directory / relative_path
        filesnpaths.gen_output_directory(path.parent, progress=self.progress, run=self.run)
        path.write_text(
            SummaryHTMLOutput(context, r=self.run, p=self.progress).render(),
            encoding="utf-8",
        )

    def get_program_requires_provides_dict(
        self, prefix: str = "../../"
    ) -> dict[str, ProgramLinks]:
        return {
            name: {
                "requires": [
                    (artifact, f"{prefix}artifacts/{artifact}")
                    for artifact in program.requires
                ],
                "provides": [
                    (artifact, f"{prefix}artifacts/{artifact}")
                    for artifact in program.provides
                ],
                "can_use": [
                    (artifact, f"{prefix}artifacts/{artifact}")
                    for artifact in program.can_use
                ],
                "can_provide": [
                    (artifact, f"{prefix}artifacts/{artifact}")
                    for artifact in program.can_provide
                ],
                "anvio_workflows": [
                    (workflow, f"{prefix}workflows/{workflow}")
                    for workflow in program.anvio_workflows
                ],
            }
            for name, program in self.dataset.programs.items()
        }

    def generate_pages_for_artifacts(self) -> None:
        self.progress.new(
            "Rendering artifact pages", progress_total_items=len(self.dataset.artifacts)
        )
        for name, artifact in self.dataset.artifacts.items():
            self.progress.update(name, increment=True)
            context: dict[str, Any] = {
                "name": name,
                "type": artifact.type,
                "provided_by_anvio": artifact.provided_by_anvio,
                "provided_by_user": artifact.provided_by_user,
                "description": self.bodies[artifact.id],
                "icon": f"../../images/icons/{artifact.type}.png",
            }
            for relation, inverse in (
                ("requires", "required_by"),
                ("provides", "provided_by"),
                ("can_use", "can_used_by"),
                ("can_provide", "can_provided_by"),
            ):
                context[inverse] = [
                    (program_name, f"../../programs/{program_name}")
                    for program_name, program in self.dataset.programs.items()
                    if name in getattr(program, relation)
                ]
            self.render_page(
                f"artifacts/{name}/index.md", "artifact", {"artifact": context}
            )
        self.progress.end()

    def get_HTML_formatted_authors_data(self, authors: list[str]) -> str:
        result = ""
        for name in authors:
            author = self.dataset.authors[name]
            result += '<div class="anvio-person"><div class="anvio-person-info">'
            result += f'<div class="anvio-person-photo"><img class="anvio-person-photo-img" src="../../images/authors/{Path(author.avatar).name}" /></div>'
            result += '<div class="anvio-person-info-box">'
            result += f'<a href="/people/{author.github}" target="_blank"><span class="anvio-person-name">{author.name}</span></a>'
            result += '<div class="anvio-person-social-box">'
            if author.web is not None:
                result += f'<a href="{author.web}" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a>'
            result += f'<a href="mailto:{author.email}" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a>'
            if author.twitter is not None:
                result += f'<a href="http://twitter.com/{author.twitter}" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a>'
            result += f'<a href="http://github.com/{author.github}" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a>'
            result += "</div></div></div></div>\n\n"
        return result

    def get_HTML_formatted_authors_data_mini(self, authors: list[str]) -> str:
        result = ""
        for name in authors:
            author = self.dataset.authors[name]
            result += (
                '<div class="anvio-person-mini"><div class="anvio-person-photo-mini">'
            )
            result += f'<a href="/people/{author.github}" target="_blank"><img class="anvio-person-photo-img-mini" title="{author.name}" src="images/authors/{Path(author.avatar).name}" /></a>'
            result += "</div></div>\n"
        return result

    def get_HTML_formatted_third_party_programs(
        self, workflow: WorkflowDocumentation
    ) -> list[str]:
        return [
            f'<a href="{self.dataset.third_party_programs[name].link}" target="_blank">{name}</a> ({purpose})'
            for purpose, names in workflow.third_party_programs_used
            for name in names
        ]

    def generate_pages_for_workflows(self) -> None:
        self.progress.new(
            "Rendering workflow pages", progress_total_items=len(self.dataset.workflows)
        )
        for name, workflow in self.dataset.workflows.items():
            self.progress.update(name, increment=True)
            context = {
                "name": name,
                "one_sentence_summary": workflow.one_sentence_summary,
                "one_paragraph_summary": workflow.one_paragraph_summary,
                "description": self.bodies[workflow.id],
                "authors": self.get_HTML_formatted_authors_data(workflow.authors),
                "artifacts_produced": [
                    (artifact, f"../../artifacts/{artifact}")
                    for artifact in workflow.artifacts_produced
                ],
                "artifacts_accepted": [
                    (artifact, f"../../artifacts/{artifact}")
                    for artifact in workflow.artifacts_accepted
                ],
                "third_party_programs_used": self.get_HTML_formatted_third_party_programs(
                    workflow
                ),
            }
            self.render_page(
                f"workflows/{name}/index.md",
                "workflow",
                {"workflow": context, "artifacts": self.artifact_context},
            )
        self.progress.end()

    def generate_pages_for_programs(self) -> None:
        self.progress.new(
            "Rendering program pages", progress_total_items=len(self.dataset.programs)
        )
        links = self.get_program_requires_provides_dict()
        example_path = self.dataset.program_source_paths.get("anvi-interactive")
        for name, program in self.dataset.programs.items():
            self.progress.update(name, increment=True)
            context = {
                "name": name,
                "usage": self.bodies[program.id],
                "description": program.description,
                "resources": program.resources,
                "source_path": program.source_path,
                "resources_example_source_path": example_path or program.source_path,
                "authors": self.get_HTML_formatted_authors_data(program.authors),
                **links[name],
            }
            self.render_page(
                f"programs/{name}/index.md",
                "program",
                {"program": context, "artifacts": self.artifact_context},
            )
            network_path = self.output_directory / "programs" / name / "network.json"
            with network_path.open("w", encoding="utf-8") as network_file:
                json.dump(
                    self.programs_network([name], artifact_names_as_ids=True),
                    network_file,
                    indent=2,
                )
        self.progress.end()

    def generate_index_page(self) -> None:
        artifact_types: dict[str, list[str]] = {}
        for name, artifact in self.dataset.artifacts.items():
            artifact_types.setdefault(artifact.type, []).append(name)
        context = {
            "programs": [
                (
                    name,
                    f"programs/{name}",
                    program.description,
                    self.get_HTML_formatted_authors_data_mini(program.authors),
                )
                for name, program in self.dataset.programs.items()
            ],
            "workflows": {
                name: {
                    "one_sentence_summary": workflow.one_sentence_summary,
                    "authors": self.get_HTML_formatted_authors_data_mini(
                        workflow.authors
                    ),
                }
                for name, workflow in self.dataset.workflows.items()
            },
            "artifacts": self.artifact_context,
            "artifact_types": artifact_types,
            "program_provides_requires": self.get_program_requires_provides_dict(
                prefix=""
            ),
        }
        self.render_page("index.md", "programs_and_artifacts_index", context)
