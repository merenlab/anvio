"""Serializable documentation records collected from anvi'o source metadata.

The dataset preserves authored Markdown, declared relationships, canonical URLs,
and an image inventory. documentation.json and its adjacent images form a portable
bundle that renderers can consume without consulting the source registries.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from functools import cached_property
from pathlib import Path

import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError, FilesNPathsError


SCHEMA_VERSION = 1
DATASET_FILENAME = "documentation.json"
RELATIONSHIPS = ("requires", "provides", "can_use", "can_provide")


@dataclass(frozen=True, kw_only=True)
class BuildMetadata:
    docs_version: str
    version: str
    versions: list[tuple[str, str]]
    date: str
    base_url: str
    version_short_identifier: str


@dataclass(frozen=True, kw_only=True)
class DocumentRecord:
    id: str
    kind: str
    source_markdown: str | None
    documentation_path: str
    html_url: str


@dataclass(frozen=True, kw_only=True)
class ProgramDocumentation(DocumentRecord):
    name: str
    source_path: str
    description: str
    authors: list[str]
    resources: list[tuple[str, str]]
    tags: list[str]
    anvio_workflows: list[str]
    requires: list[str]
    provides: list[str]
    can_use: list[str]
    can_provide: list[str]


@dataclass(frozen=True, kw_only=True)
class ArtifactDocumentation(DocumentRecord):
    name: str
    type: str
    provided_by_anvio: bool
    provided_by_user: bool


@dataclass(frozen=True, kw_only=True)
class WorkflowDocumentation(DocumentRecord):
    name: str
    authors: list[str]
    artifacts_produced: list[str]
    artifacts_accepted: list[str]
    anvio_workflows_inherited: list[str]
    third_party_programs_used: list[tuple[str, list[str]]]
    one_sentence_summary: str
    one_paragraph_summary: str
    anvio_programs_used: list[str]


@dataclass(frozen=True, kw_only=True)
class AuthorDocumentation:
    github: str
    name: str
    email: str
    avatar: str
    web: str | None = None
    twitter: str | None = None


@dataclass(frozen=True, kw_only=True)
class ThirdPartyProgram:
    link: str


@dataclass(frozen=True, kw_only=True)
class DocumentationAsset:
    path: str
    md5: str


@dataclass(kw_only=True)
class DocumentationData:
    """Typed snapshot shared by serialization, page rendering, and network generation."""

    meta: BuildMetadata
    program_source_paths: dict[str, str]
    programs: dict[str, ProgramDocumentation]
    artifacts: dict[str, ArtifactDocumentation]
    workflows: dict[str, WorkflowDocumentation]
    authors: dict[str, AuthorDocumentation]
    third_party_programs: dict[str, ThirdPartyProgram]
    assets: list[DocumentationAsset]
    schema_version: int = SCHEMA_VERSION
    dataset_path: Path | None = field(default=None, repr=False, compare=False)


    @cached_property
    def asset_sources(self) -> dict[str, Path]:
        if self.dataset_path is None:
            raise ConfigError("DocumentationData does not know where to find its image assets. Please load an exported "
                              "dataset with DocumentationData.load(), collect one with AnvioDocs.collect(), or provide "
                              "an asset_sources dictionary that maps each asset path to its source file.")

        return {asset.path: self.dataset_path.parent / asset.path for asset in self.assets}


    @classmethod
    def load(cls, path: Path) -> DocumentationData:
        """Load an exported documentation dataset and locate its adjacent image assets."""

        path = Path(path)
        filesnpaths.is_file_exists(path)

        try:
            with path.open(encoding="utf-8") as source:
                data = json.load(source)
        except OSError as error:
            raise FilesNPathsError(f"Anvi'o could not read the documentation dataset at '{path}'. Please check that this "
                                   f"is a readable file. Here is the error from the operating system: {error}") from error
        except (json.JSONDecodeError, UnicodeError) as error:
            raise ConfigError(f"The documentation dataset at '{path}' could not be read as UTF-8 JSON. Please use the "
                              f"documentation.json produced by anvi-script-gen-help-pages, or regenerate the export. "
                              f"Here is what the reader reported: {error}") from error

        try:
            if data["schema_version"] != SCHEMA_VERSION:
                raise ConfigError(f"The documentation dataset at '{path}' uses schema version '{data['schema_version']}', "
                                  f"but this version of anvi'o reads schema version {SCHEMA_VERSION}. Please regenerate "
                                  f"the dataset with this version of anvi'o, or use a compatible anvi'o installation.")

            # Restore tuple fields from JSON arrays.
            data["meta"]["versions"] = [(name, version) for name, version in data["meta"]["versions"]]
            for record in data["programs"].values():
                record["resources"] = [(label, url) for label, url in record["resources"]]
            for record in data["workflows"].values():
                record["third_party_programs_used"] = [(operation, programs) for operation, programs in record["third_party_programs_used"]]

            return cls(schema_version=data["schema_version"],
                       meta=BuildMetadata(**data["meta"]),
                       program_source_paths=data["program_source_paths"],
                       programs={name: ProgramDocumentation(**record) for name, record in data["programs"].items()},
                       artifacts={name: ArtifactDocumentation(**record) for name, record in data["artifacts"].items()},
                       workflows={name: WorkflowDocumentation(**record) for name, record in data["workflows"].items()},
                       authors={name: AuthorDocumentation(**record) for name, record in data["authors"].items()},
                       third_party_programs={name: ThirdPartyProgram(**record) for name, record in data["third_party_programs"].items()},
                       assets=[DocumentationAsset(**record) for record in data["assets"]],
                       dataset_path=path)
        except (KeyError, TypeError, ValueError, AttributeError) as error:
            raise ConfigError(f"The JSON file at '{path}' does not contain the documentation records anvi'o expected. "
                              f"Please use a complete documentation.json export from anvi-script-gen-help-pages. "
                              f"The record reader reported: {error}") from error


    def write(self, path: Path) -> None:
        """Write the dataset as JSON; image assets are copied separately by the renderer."""

        path = Path(path)
        filesnpaths.is_output_file_writable(path)
        data = asdict(self)
        # Local input locations are runtime state, not part of the portable dataset.
        del data["dataset_path"]
        try:
            with path.open("w", encoding="utf-8") as output:
                json.dump(data, output, ensure_ascii=False, indent=2)
                output.write("\n")
        except OSError as error:
            raise FilesNPathsError(f"Anvi'o could not write the documentation dataset to '{path}'. Please check the "
                                   f"available disk space and your write permissions, or choose another output path. "
                                   f"Here is the error from the operating system: {error}") from error
