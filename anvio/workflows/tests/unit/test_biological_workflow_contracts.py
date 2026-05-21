"""
    Unit tests for biologically meaningful workflow contracts.

    These tests intentionally avoid importing the workflow modules. Importing
    anvi'o requires the exact runtime environment and many external programs,
    while these contracts only need to make sure each workflow still advertises
    the biological steps that define its purpose.
"""

import re
import unittest
from pathlib import Path


WORKFLOWS_ROOT = Path(__file__).resolve().parents[2]


METAGENOMICS_RULE_FILES = sorted((WORKFLOWS_ROOT / 'metagenomics' / 'rules').glob('*.smk'))


WORKFLOWS = {
    'contigs': {
        'module': WORKFLOWS_ROOT / 'contigs' / '__init__.py',
        'snakefile': WORKFLOWS_ROOT / 'contigs' / 'Snakefile',
        'rules': {
            'anvi_gen_contigs_database',
            'anvi_run_hmms',
            'anvi_run_scg_taxonomy',
            'anvi_run_ncbi_cogs',
            'anvi_run_kegg_kofams',
        },
        'tokens': {
            "'--also-scan-trnas': True",
            '"--simplify-names": True',
            '"--project-name": "{group}"',
        },
    },
    'metagenomics': {
        'module': WORKFLOWS_ROOT / 'metagenomics' / '__init__.py',
        'snakefile': WORKFLOWS_ROOT / 'metagenomics' / 'Snakefile',
        'extra_snakefiles': [WORKFLOWS_ROOT / 'read_recruitment' / 'rules' / 'main.smk'],
        'rules': {
            'iu_filter_quality_minoche',
            'megahit',
            'idba_ud',
            'metaspades',
            'flye',
            'bowtie',
            'minimap2',
            'anvi_profile',
            'anvi_merge',
            'krakenuniq',
            'anvi_cluster_contigs',
        },
        'tokens': {
            'min_contig_length_for_assembly = 1000',
            '"--min-contig-len": min_contig_length_for_assembly',
            '"--min_contig": min_contig_length_for_assembly',
            "'--gzip-compressed': True",
        },
    },
    'pangenomics': {
        'module': WORKFLOWS_ROOT / 'pangenomics' / '__init__.py',
        'snakefile': WORKFLOWS_ROOT / 'pangenomics' / 'Snakefile',
        'rules': {
            'anvi_gen_genomes_storage',
            'anvi_pan_genome',
            'anvi_get_sequences_for_gene_clusters',
            'import_phylogenetic_tree_to_pangenome',
        },
        'tokens': {
            "self.valid_sequence_sources_for_phylogeny = ['gene_clusters', 'hmm']",
            "'--concatenate-gene-clusters'",
            "'--align-with'",
        },
    },
    'phylogenomics': {
        'module': WORKFLOWS_ROOT / 'phylogenomics' / '__init__.py',
        'snakefile': WORKFLOWS_ROOT / 'phylogenomics' / 'Snakefile',
        'rules': {
            'anvi_get_sequences_for_hmm_hits',
            'trimal',
            'iqtree',
        },
        'tokens': {
            "'--return-best-hit': True",
            "'--concatenate-genes': True",
            "'--get-aa-sequences': True",
            "'--hmm-sources': 'Bacteria_71'",
            "'-bb': 1000",
        },
    },
    'sra_download': {
        'module': WORKFLOWS_ROOT / 'sra_download' / '__init__.py',
        'snakefile': WORKFLOWS_ROOT / 'sra_download' / 'Snakefile',
        'rules': {
            'prefetch',
            'fasterq_dump',
            'pigz',
            'generate_samples_txt',
        },
        'tokens': {
            "'SRR', 'ERR', 'DRR'",
            "'--split-files'",
            "'Remove_unzipped_SRA_files': True",
        },
    },
    'trnaseq': {
        'module': WORKFLOWS_ROOT / 'trnaseq' / '__init__.py',
        'snakefile': WORKFLOWS_ROOT / 'trnaseq' / 'Snakefile',
        'rules': {
            'iu_merge_pairs',
            'anvi_reformat_fasta',
            'anvi_trnaseq',
            'anvi_merge_trnaseq',
            'anvi_run_trna_taxonomy',
            'anvi_tabulate_trnaseq',
        },
        'tokens': {
            "KNOWN_TREATMENTS = ['untreated', 'demethylase']",
            "'--skip-fasta-check': True",
            "'--min-percent-identity': 90",
            "'--max-num-target-sequences': 100",
        },
    },
    'ecophylo': {
        'module': WORKFLOWS_ROOT / 'ecophylo' / '__init__.py',
        'snakefile': WORKFLOWS_ROOT / 'ecophylo' / 'Snakefile',
        'rules': {
            'anvi_run_hmms_hmmsearch',
            'filter_hmm_hits_by_model_coverage',
            'cluster_X_percent_sim_mmseqs',
            'align_sequences',
            'trim_alignment',
            'fasttree',
            'iqtree',
            'run_metagenomics_workflow',
            'anvi_estimate_scg_taxonomy',
        },
        'tokens': {
            "'--min-model-coverage': 0.8",
            "'--filter-out-partial-gene-calls': True",
            "'--min-seq-id': 0.97",
            "'clustering_threshold_for_OTUs': [0.99, 0.98]",
            "'scg_taxonomy_database_version': \"GTDB: v214.1; Anvi'o: v1\"",
        },
    },
}


def read_text(path):
    return path.read_text()


def snakefile_rule_names(paths):
    rule_names = set([])

    for path in paths:
        text = read_text(path)
        rule_names.update(re.findall(r'^\s*rule\s+([A-Za-z_][A-Za-z0-9_]*)\s*:', text, re.MULTILINE))

    return rule_names


def workflow_snakefiles(contract):
    paths = [contract['snakefile']]
    rules_dir = contract['snakefile'].parent / 'rules'
    if rules_dir.exists():
        paths.extend(sorted(rules_dir.glob('*.smk')))
    paths.extend(contract.get('extra_snakefiles', []))
    return paths


class BiologicalWorkflowContractTestCase(unittest.TestCase):
    def test_each_workflow_keeps_domain_defining_rules(self):
        for workflow_name, contract in WORKFLOWS.items():
            with self.subTest(workflow=workflow_name):
                rule_names = snakefile_rule_names(workflow_snakefiles(contract))
                self.assertTrue(contract['rules'].issubset(rule_names),
                                f"{workflow_name} is missing expected biological rules: "
                                f"{sorted(contract['rules'] - rule_names)}")


    def test_each_workflow_keeps_biologically_relevant_defaults_and_guards(self):
        for workflow_name, contract in WORKFLOWS.items():
            with self.subTest(workflow=workflow_name):
                module_text = read_text(contract['module']).replace('"', "'")

                for token in contract['tokens']:
                    self.assertIn(token.replace('"', "'"), module_text,
                                  f"{workflow_name} is missing biological contract token: {token}")


    def test_short_and_long_read_metagenomics_paths_remain_distinct(self):
        snakefile_text = read_text(WORKFLOWS['metagenomics']['snakefile'])
        read_recruitment_text = read_text(WORKFLOWS_ROOT / 'read_recruitment' / 'rules' / 'main.smk')
        combined_text = snakefile_text + read_recruitment_text + ''.join(
            read_text(path) for path in METAGENOMICS_RULE_FILES
        )

        self.assertIn('SR_RS_RE = w.regex_from_ids(SR_READSETS)', snakefile_text)
        self.assertIn('LR_RS_RE = w.regex_from_ids(LR_READSETS)', snakefile_text)
        self.assertRegex(combined_text, r'group\s*=\s*SR_GRP_RE')
        self.assertRegex(combined_text, r'group\s*=\s*LR_GRP_RE')
        self.assertIn('bowtie2', combined_text)
        self.assertIn('minimap2', combined_text)


    def test_ecophylo_preserves_profile_and_tree_modes(self):
        profile_mode_text = read_text(WORKFLOWS_ROOT / 'ecophylo' / 'rules' / 'profile_mode.smk')
        tree_mode_text = read_text(WORKFLOWS_ROOT / 'ecophylo' / 'rules' / 'tree_mode.smk')

        self.assertIn('run_metagenomics_workflow', profile_mode_text)
        self.assertIn('add_default_collection', profile_mode_text)
        self.assertIn('anvi_import_everything_metagenome', profile_mode_text)
        self.assertIn('anvi_import_everything_tree', tree_mode_text)


if __name__ == '__main__':
    unittest.main()
