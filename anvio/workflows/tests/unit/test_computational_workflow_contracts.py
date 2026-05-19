"""
    Unit tests for computationally meaningful workflow contracts.

    These tests focus on static Snakemake hygiene: parseable rule structure,
    log/resource declarations, reproducibility-oriented environment hooks, and
    avoiding deprecated rule directives.
"""

import re
import unittest
from pathlib import Path


WORKFLOWS_ROOT = Path(__file__).resolve().parents[2]

SNAKEMAKE_FILES = [
    WORKFLOWS_ROOT / 'contigs' / 'Snakefile',
    WORKFLOWS_ROOT / 'ecophylo' / 'Snakefile',
    WORKFLOWS_ROOT / 'ecophylo' / 'rules' / 'profile_mode.smk',
    WORKFLOWS_ROOT / 'ecophylo' / 'rules' / 'tree_mode.smk',
    WORKFLOWS_ROOT / 'metagenomics' / 'Snakefile',
    WORKFLOWS_ROOT / 'pangenomics' / 'Snakefile',
    WORKFLOWS_ROOT / 'phylogenomics' / 'Snakefile',
    WORKFLOWS_ROOT / 'sra_download' / 'Snakefile',
    WORKFLOWS_ROOT / 'trnaseq' / 'Snakefile',
]

WORKFLOW_MODULES = [
    WORKFLOWS_ROOT / 'contigs' / '__init__.py',
    WORKFLOWS_ROOT / 'ecophylo' / '__init__.py',
    WORKFLOWS_ROOT / 'metagenomics' / '__init__.py',
    WORKFLOWS_ROOT / 'pangenomics' / '__init__.py',
    WORKFLOWS_ROOT / 'phylogenomics' / '__init__.py',
    WORKFLOWS_ROOT / 'sra_download' / '__init__.py',
    WORKFLOWS_ROOT / 'trnaseq' / '__init__.py',
]

TARGET_RULE_NAMES = {
    'contigs_workflow_target_rule',
    'ecophylo_workflow_target_rule',
    'metagenomics_workflow_target_rule',
    'pangenomics_workflow_target_rule',
    'phylogenomics_workflow_target_rule',
    'sra_download_workflow_target_rule',
    'trnaseq_workflow_target_rule',
}


def is_target_rule(rule_name):
    return rule_name in TARGET_RULE_NAMES or rule_name.lower().endswith('target_rule')


def read_text(path):
    return path.read_text()


def rule_blocks(path):
    text = read_text(path)
    matches = list(re.finditer(r'^(?P<indent>\s*)rule\s+(?P<name>[A-Za-z_][A-Za-z0-9_]*)\s*:', text, re.MULTILINE))

    for index, match in enumerate(matches):
        end = matches[index + 1].start() if index + 1 < len(matches) else len(text)
        yield match.group('name'), match.group('indent'), text[match.start():end]


def directive_names(rule_block):
    return set(re.findall(r'^\s{4,}([A-Za-z_][A-Za-z0-9_]*)\s*:', rule_block, re.MULTILINE))


class ComputationalWorkflowContractTestCase(unittest.TestCase):
    def test_rule_version_directives_are_static_when_present(self):
        for path in SNAKEMAKE_FILES:
            with self.subTest(path=path.relative_to(WORKFLOWS_ROOT)):
                versions = re.findall(r'^\s*version\s*:\s*(.+)$', read_text(path), re.MULTILINE)

                for version in versions:
                    self.assertRegex(version, r'^(1\.0|anvio\.__[A-Za-z0-9_]+__version__)(\s*#.*)?$')


    def test_rules_have_declared_outputs_or_are_explicit_target_rules(self):
        for path in SNAKEMAKE_FILES:
            for rule_name, _, block in rule_blocks(path):
                with self.subTest(path=path.relative_to(WORKFLOWS_ROOT), rule=rule_name):
                    directives = directive_names(block)
                    if is_target_rule(rule_name):
                        self.assertIn('input', directives)
                    else:
                        self.assertIn('output', directives)


    def test_non_target_rules_have_logs_for_reproducible_diagnostics(self):
        allowed_without_logs = {
            'all_reformatting_done',
            'cat_anvi_profile_blitz',
            'copy_contigs_db',
            'import_percent_of_reads_mapped_for_single_profiles',
        }

        for path in SNAKEMAKE_FILES:
            for rule_name, _, block in rule_blocks(path):
                with self.subTest(path=path.relative_to(WORKFLOWS_ROOT), rule=rule_name):
                    if is_target_rule(rule_name) or rule_name in allowed_without_logs:
                        continue

                    self.assertIn('log', directive_names(block))


    def test_declared_node_resources_are_tied_to_workflow_thread_controls(self):
        for path in SNAKEMAKE_FILES:
            for rule_name, _, block in rule_blocks(path):
                with self.subTest(path=path.relative_to(WORKFLOWS_ROOT), rule=rule_name):
                    if not re.search(r'^\s+resources\s*:', block, re.MULTILINE):
                        continue

                    self.assertIn('nodes', block)
                    self.assertRegex(block, r'nodes\s*=\s*(M\.T\(|w\.T\(|1)')


    def test_shell_rules_redirect_stdout_and_stderr_to_logs(self):
        for path in SNAKEMAKE_FILES:
            for rule_name, _, block in rule_blocks(path):
                with self.subTest(path=path.relative_to(WORKFLOWS_ROOT), rule=rule_name):
                    if 'shell:' not in directive_names(block):
                        continue

                    self.assertIn('{log}', block)
                    self.assertRegex(block, r'(2>&1|2>>\s*\{log\}|2>>\{log\})')


    def test_metagenomics_external_tool_rules_support_conda_or_existing_envs(self):
        metagenomics_module = read_text(WORKFLOWS_ROOT / 'metagenomics' / '__init__.py')
        metagenomics_snakefile = read_text(WORKFLOWS_ROOT / 'metagenomics' / 'Snakefile')

        for tool in ['flye', 'minimap2', 'bowtie', 'megahit', 'metaspades', 'idba_ud']:
            with self.subTest(tool=tool):
                self.assertIn('"conda_yaml"', metagenomics_module)
                self.assertIn('"conda_env"', metagenomics_module)
                self.assertIn(f"self.ensure_tool_in_path_or_conda('{tool}'", metagenomics_module)

        self.assertIn('def get_conda_yaml_path(workflow, tool):', read_text(WORKFLOWS_ROOT / '__init__.py'))
        self.assertIn('def get_conda_env_prefix(workflow, tool):', read_text(WORKFLOWS_ROOT / '__init__.py'))
        self.assertIn('w.get_conda_yaml_path(M,', metagenomics_snakefile)
        self.assertIn('w.get_conda_env_prefix(M,', metagenomics_snakefile)
        self.assertIn('conda:', metagenomics_snakefile)


    def test_workflow_modules_define_default_max_threads_or_rule_threads(self):
        for path in WORKFLOW_MODULES:
            with self.subTest(path=path.relative_to(WORKFLOWS_ROOT)):
                text = read_text(path)
                self.assertTrue("'threads'" in text or '"threads"' in text or "'max_threads'" in text or '"max_threads"' in text,
                                'workflow module should expose thread controls in default config')


if __name__ == '__main__':
    unittest.main()
