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
        assert self.dataset_path is not None, "No assets because dataset_path not provided"
        return {
            asset.path: self.dataset_path.parent / asset.path for asset in self.assets
        }

    @classmethod
    def load(cls, path: Path) -> DocumentationData:
        with path.open(encoding="utf-8") as source:
            data = json.load(source)

        # Restore tuple fields from JSON arrays.
        data["meta"]["versions"] = [(name, version) for name, version in data["meta"]["versions"]]
        for record in data["programs"].values():
            record["resources"] = [(label, url) for label, url in record["resources"]]
        for record in data["workflows"].values():
            record["third_party_programs_used"] = [
                (operation, programs)
                for operation, programs in record["third_party_programs_used"]
            ]

        return cls(
            schema_version=data["schema_version"],
            meta=BuildMetadata(**data["meta"]),
            program_source_paths=data["program_source_paths"],
            programs={
                name: ProgramDocumentation(**record)
                for name, record in data["programs"].items()
            },
            artifacts={
                name: ArtifactDocumentation(**record)
                for name, record in data["artifacts"].items()
            },
            workflows={
                name: WorkflowDocumentation(**record)
                for name, record in data["workflows"].items()
            },
            authors={
                name: AuthorDocumentation(**record)
                for name, record in data["authors"].items()
            },
            third_party_programs={
                name: ThirdPartyProgram(**record)
                for name, record in data["third_party_programs"].items()
            },
            assets=[DocumentationAsset(**record) for record in data["assets"]],
            dataset_path=path,
        )

    def write(self, path: Path) -> None:
        data = asdict(self)
        # Local input locations are runtime state, not part of the portable dataset.
        del data["dataset_path"]
        with path.open("w", encoding="utf-8") as output:
            json.dump(data, output, ensure_ascii=False, indent=2)
            output.write("\n")
