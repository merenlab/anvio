"""
    Unit tests for workflow dependency extraction helpers.
"""

import types
import unittest
import importlib.util
from pathlib import Path


WORKFLOWS_ROOT = Path(__file__).resolve().parents[2]
DEPENDENCIES_PATH = WORKFLOWS_ROOT / 'scripts' / 'dependencies.py'


def load_dependencies_module():
    spec = importlib.util.spec_from_file_location('test_workflow_dependencies', DEPENDENCIES_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    return module


class WorkflowDependencyParserTestCase(unittest.TestCase):
    def setUp(self):
        self.dependencies = load_dependencies_module()


    def test_parser_skips_comments_shell_builtins_and_templated_prefixes(self):
        shellcmd = r"""
            # Run an assembler through an optional conda/env prefix.
            mkdir -p {params.outdir}

            {params.env_prefix} flye {params.meta} {params.read_type} {input.reads} \
                 -o {params.outdir} -t {threads} \
                 >> {log} 2>&1

            if [ -f "{params.outdir}/assembly.fasta" ]; then
                cp "{params.outdir}/assembly.fasta" {output.raw_fasta} && rm "{params.outdir}/assembly.fasta"
            else
                echo "Could not find assembly fasta" >> {log}
                exit 1
            fi
        """

        self.assertEqual(self.dependencies.get_programs_from_shell_command(shellcmd), ['flye'])


    def test_parser_extracts_multiple_real_commands_from_multiline_shell(self):
        shellcmd = r"""
            {params.env_prefix} bowtie2 --threads {threads} \
                -x {params.index_prefix} \
                -1 {params.r1} -2 {params.r2} \
                -S {output.sam} >> {log} 2>&1
            samtools view -bS {input} -o {output} >> {log} 2>&1
            anvi-init-bam {input} -o {output.bam} -T {threads} >> {log} 2>&1
        """

        self.assertEqual(self.dependencies.get_programs_from_shell_command(shellcmd),
                         ['bowtie2', 'samtools', 'anvi-init-bam'])


    def test_parser_ignores_flat_file_utilities_that_are_not_workflow_dependencies(self):
        shellcmd = r"""
            echo "message" >> {log}
            gzip {input.fastq} >> {log} 2>&1
            mv {params.summary_temp_dir} {output} >> {log} 2>&1
            sed 's/^[0-9]/g&/' {input} > {output}
            krakenuniq --db {params.db} {input} > {output}
        """

        self.assertEqual(self.dependencies.get_programs_from_shell_command(shellcmd), ['krakenuniq'])


    def test_parser_keeps_multiline_quoted_echo_messages_together(self):
        shellcmd = r"""
            echo -e 'The group {wildcards.group} has only one sample. Hence, there is nothing to merge, but you can find\n\
                     the profile database here: {input.profiles}. Also, just so you know, profile was done using\n\
                     --cluster-contigs so you can visualize this profile database using anvi-interactive.' > {output} 2>>{log}

            anvi-merge {input.profiles} -o {output.merged}
        """

        self.assertEqual(self.dependencies.get_programs_from_shell_command(shellcmd), ['anvi-merge'])


    def test_workflow_extraction_deduplicates_programs_while_preserving_order(self):
        workflow = types.SimpleNamespace(rules=[
            types.SimpleNamespace(shellcmd='anvi-script-reformat-fasta {input} -o {output}'),
            types.SimpleNamespace(shellcmd='anvi-script-reformat-fasta {input} -o {output}'),
            types.SimpleNamespace(shellcmd='{params.env_prefix} bowtie2 --threads {threads}'),
        ])

        self.assertEqual(self.dependencies.get_programs_from_snakemake_workflow(workflow),
                         ['anvi-script-reformat-fasta', 'bowtie2'])


if __name__ == '__main__':
    unittest.main()
