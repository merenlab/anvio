"""Exercise the actual workflow read_csv expressions without importing Snakemake.

Extracting expressions lets these parsing regressions run on both scientific
stacks without bypassing anvi'o's interpreter guard or executing external tools.
"""
import ast
import io
from pathlib import Path
import re
import runpy
import tempfile
import tokenize
from types import SimpleNamespace
import unittest

import pandas as pd
from pandas.api.types import is_numeric_dtype

WORKFLOWS = Path(__file__).resolve().parents[3] / 'workflows'


def reader_expressions(relative_path):
    source = (WORKFLOWS / relative_path).read_text()
    for match in re.finditer(r'pd\.read_csv\(', source):
        tail = source[match.start():]
        lines = tail.splitlines(keepends=True)
        depth = 0
        for token in tokenize.generate_tokens(io.StringIO(tail).readline):
            if token.type != tokenize.OP:
                continue
            if token.string == '(':
                depth += 1
            elif token.string == ')':
                depth -= 1
                if depth == 0:
                    row, column = token.end
                    expression = ast.parse(''.join(lines[:row-1]) + lines[row-1][:column], mode='eval')
                    yield expression
                    break


def reader_expression(relative_path, argument):
    for expression in reader_expressions(relative_path):
        if ast.unparse(expression.body.args[0]) == argument:
            return expression
    raise AssertionError(f'Reader not found: {relative_path}: {argument}')


class IdentifierParsingTests(unittest.TestCase):
    def parse(self, path, argument, text):
        expression = reader_expression(path, argument)
        with tempfile.TemporaryDirectory() as directory:
            fixture = Path(directory) / 'input.tsv'
            fixture.write_text(text)
            # Keep every actual production keyword; replace only the input path.
            expression.body.args[0] = ast.Constant(str(fixture))
            ast.fix_missing_locations(expression)
            return eval(compile(expression, path, 'eval'), {'pd': pd, 'str': str})

    def test_changed_reader_expressions_compile_with_unique_keywords(self):
        paths = ['contigs/rules/contigs_database.smk', 'contigs/rules/functional_annotation.smk',
                 'ecophylo/__init__.py', 'ecophylo/rules/clustering_alignment.smk',
                 'ecophylo/rules/hmm_processing.smk', 'ecophylo/rules/misc_taxonomy.smk',
                 'ecophylo/rules/sequence_data.smk', 'ecophylo/scripts/subset_external_gene_calls_file.py',
                 'sra_download/__init__.py']
        for path in paths:
            for expression in reader_expressions(path):
                argument = ast.unparse(expression.body.args[0])
                with self.subTest(path=path, argument=argument):
                    keywords = [keyword.arg for keyword in expression.body.keywords]
                    self.assertEqual(len(keywords), len(set(keywords)))
                    compile(expression, path, 'eval')

    def test_header_and_cluster_identifiers(self):
        cases = [
            ('ecophylo/rules/hmm_processing.smk', 'reformat_file'),
            ('ecophylo/rules/clustering_alignment.smk', 'input.mmseqs_cluster_rep_index'),
            ('ecophylo/rules/clustering_alignment.smk', 'input.reformat_report'),
            ('ecophylo/rules/misc_taxonomy.smk', 'mmseqs_cluster_rep_index'),
            ('ecophylo/rules/misc_taxonomy.smk', 'params.reformat_file'),
        ]
        for path, argument in cases:
            with self.subTest(path=path, argument=argument):
                result = self.parse(path, argument, 'NA\t001\nNULL\t002\n001\t003\n')
                self.assertEqual(result.iloc[:, 0].tolist(), ['NA', 'NULL', '001'])
                self.assertEqual(result.iloc[:, 1].tolist(), ['001', '002', '003'])

    def test_single_header_identifiers(self):
        cases = [
            ('ecophylo/rules/misc_taxonomy.smk', 'final_sequences_headers'),
            ('ecophylo/rules/misc_taxonomy.smk', 'input.final_list_of_sequences_for_mapping_headers'),
            ('ecophylo/scripts/subset_external_gene_calls_file.py', 'snakemake.input.headers'),
        ]
        for path, argument in cases:
            with self.subTest(path=path, argument=argument):
                result = self.parse(path, argument, 'NA\nNULL\n001\n')
                self.assertEqual(result.iloc[:, 0].tolist(), ['NA', 'NULL', '001'])

    def test_contigs_reformat_report(self):
        result = self.parse('contigs/rules/contigs_database.smk', 'input.reformat_report[0]',
                            '001\tNA\n002\tNULL\n003\t001\n')
        self.assertEqual(result.index.tolist(), ['NA', 'NULL', '001'])
        self.assertEqual(result.iloc[:, 0].tolist(), ['001', '002', '003'])

    def test_sra_accessions_remain_strings(self):
        result = self.parse('sra_download/__init__.py', 'self.SRA_accession_list',
                            'NA\nNULL\n001\n')
        self.assertEqual(result['accessions'].tolist(), ['NA', 'NULL', '001'])

    def test_hmm_list_identifiers_and_blank_group_fallback(self):
        tree = ast.parse((WORKFLOWS / 'ecophylo/__init__.py').read_text())
        statements = [node for node in ast.walk(tree)
                      if (isinstance(node, ast.Assign) and ast.unparse(node.targets[0]) == "hmm_df['id']")
                      or (isinstance(node, ast.If) and ast.unparse(node.test) == "'group' not in hmm_df")]
        self.assertEqual(len(statements), 2)
        statements.sort(key=lambda node: node.lineno)
        program = compile(ast.Module(body=statements, type_ignores=[]), 'actual HMM group fallback', 'exec')
        for group_column in [True, False]:
            with self.subTest(group_column=group_column):
                text = 'name\tsource\tpath' + ('\tgroup' if group_column else '') + '\n'
                rows = [('NA', 'NULL', '001', ''), ('NULL', '001', '002', 'NULL'), ('001', 'NA', '003', '001')]
                text += ''.join('\t'.join(row if group_column else row[:3]) + '\n' for row in rows)
                result = self.parse('ecophylo/__init__.py', 'self.hmm_list_path', text)
                self.assertEqual(result['name'].tolist(), ['NA', 'NULL', '001'])
                self.assertEqual(result['source'].tolist(), ['NULL', '001', 'NA'])
                self.assertEqual(result['path'].tolist(), ['001', '002', '003'])
                exec(program, {'hmm_df': result})
                self.assertEqual(result['id'].tolist(), ['NULL_NA', '001_NULL', 'NA_001'])
                expected = ['NULL_NA', 'NULL', '001'] if group_column else result['id'].tolist()
                self.assertEqual(result['group'].tolist(), expected)

    def test_external_table_identifiers(self):
        for argument in ['self.metagenomes', 'self.external_genomes']:
            with self.subTest(argument=argument):
                result = self.parse('ecophylo/__init__.py', argument,
                                    'name\tpath\nNA\t001\nNULL\t002\n001\t003\n')
                self.assertEqual(result['name'].tolist(), ['NA', 'NULL', '001'])
                self.assertEqual(result['path'].tolist(), ['001', '002', '003'])

    def test_function_tables_keep_numeric_gene_ids(self):
        for argument in ['gene_functional_annotation_file', 'external_gene_calls_file']:
            with self.subTest(argument=argument):
                result = self.parse('contigs/rules/functional_annotation.smk', argument,
                                    'gene_callers_id\taccession\tannotation\n1\tNA\tNULL\n2\tNULL\tNA\n')
                self.assertEqual(result['accession'].tolist(), ['NA', 'NULL'])
                self.assertEqual(result['annotation'].tolist(), ['NULL', 'NA'])
                self.assertEqual(result.index.tolist(), [1, 2])
                self.assertTrue(is_numeric_dtype(result.index.dtype))

    def test_hmm_metadata_keeps_numeric_measurements(self):
        cases = [('ecophylo/rules/hmm_processing.smk', 'input.hmm_hits'),
                 ('ecophylo/rules/hmm_processing.smk', 'params.hmm_hits'),
                 ('ecophylo/rules/sequence_data.smk', 'hmm_hit')]
        for path, argument in cases:
            with self.subTest(path=path, argument=argument):
                result = self.parse(path, argument, 'source\tgene_name\tscore\nNA\t001\t2.5\nNULL\tNULL\t3.5\n001\tNA\t4.5\n')
                self.assertEqual(result['source'].tolist(), ['NA', 'NULL', '001'])
                self.assertEqual(result['gene_name'].tolist(), ['001', 'NULL', 'NA'])
                self.assertTrue(is_numeric_dtype(result['score']))
                self.assertEqual(result['score'].tolist(), [2.5, 3.5, 4.5])

    def test_external_gene_calls_keep_numeric_coordinates(self):
        cases = [('ecophylo/rules/hmm_processing.smk', 'external_gene_calls'),
                 ('ecophylo/scripts/subset_external_gene_calls_file.py', 'snakemake.params.external_gene_calls_all')]
        for path, argument in cases:
            with self.subTest(path=path, argument=argument):
                result = self.parse(path, argument, 'contig source version start stop\nNA NULL 001 1 5\nNULL NA 002 2 6\n001 001 NULL 3 7\n')
                self.assertEqual(result['contig'].tolist(), ['NA', 'NULL', '001'])
                self.assertEqual(result['source'].tolist(), ['NULL', 'NA', '001'])
                self.assertEqual(result['version'].tolist(), ['001', '002', 'NULL'])
                self.assertTrue(is_numeric_dtype(result['start']))
                self.assertTrue(is_numeric_dtype(result['stop']))

    def test_coverage_keys_keep_numeric_coverage_and_missing_values(self):
        cases = [('ecophylo/rules/clustering_alignment.smk', 'input.coverages', ['contig', 'gene_callers_id']),
                 ('ecophylo/rules/misc_taxonomy.smk', 'params.taxonomy_long', ['metagenome_name', 'gene_name', 'gene_callers_id'])]
        for path, argument, columns in cases:
            with self.subTest(path=path, argument=argument):
                text = '\t'.join(columns + ['coverage']) + '\n'
                text += '\t'.join(['NA'] * len(columns) + ['2.5']) + '\n'
                text += '\t'.join(['NULL'] * len(columns) + ['']) + '\n'
                text += '\t'.join(['001'] * len(columns) + ['4.5']) + '\n'
                result = self.parse(path, argument, text)
                for column in columns:
                    self.assertEqual(result[column].tolist(), ['NA', 'NULL', '001'])
                self.assertTrue(is_numeric_dtype(result['coverage']))
                self.assertEqual(result['coverage'].iloc[0], 2.5)
                self.assertTrue(pd.isna(result['coverage'].iloc[1]))

    def test_subset_script_retains_matching_identifier_rows(self):
        with tempfile.TemporaryDirectory() as directory:
            calls, headers, output = [Path(directory) / name for name in ['calls.tsv', 'headers.tsv', 'output.tsv']]
            calls.write_text('gene_callers_id contig source version start stop\n1 NA NULL 001 1 5\n2 NULL NA 002 2 6\n3 001 NULL 003 3 7\n4 other NA 004 4 8\n')
            headers.write_text('NA\nNULL\n001\n')
            snakemake = SimpleNamespace(params=SimpleNamespace(external_gene_calls_all=str(calls)),
                                        input=SimpleNamespace(headers=str(headers)),
                                        output=SimpleNamespace(external_gene_calls_subset=str(output)))
            runpy.run_path(str(WORKFLOWS / 'ecophylo/scripts/subset_external_gene_calls_file.py'), init_globals={'snakemake': snakemake})
            result = pd.read_csv(output, sep='\t', keep_default_na=False, dtype={'contig': str, 'version': str})
            self.assertEqual(result['contig'].tolist(), ['NA', 'NULL', '001'])
            self.assertEqual(result['version'].tolist(), ['001', '002', '003'])
            self.assertEqual(result['gene_callers_id'].tolist(), [0, 1, 2])
            self.assertEqual(result['start'].tolist(), [1, 2, 3])


if __name__ == '__main__':
    unittest.main()
