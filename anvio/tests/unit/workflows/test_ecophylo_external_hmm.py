import os
import shutil
import tempfile
import unittest

from anvio.errors import ConfigError
from anvio.workflows.ecophylo import EcoPhyloWorkflow


class TestEcoPhyloExternalHMM(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.external_hmm_fixture = os.path.join(
            os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
            'sandbox',
            'workflows',
            'ecophylo',
            'Ribosomal_L16_external_hmm',
        )

    def make_workflow(self, hmm_name, group=None):
        temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(temporary_directory.cleanup)

        source_name = 'Ribosomal_L16_external_hmm'
        hmm_path = os.path.join(temporary_directory.name, source_name)
        shutil.copytree(self.external_hmm_fixture, hmm_path)

        hmm_list_path = os.path.join(temporary_directory.name, 'hmm_list.txt')
        with open(hmm_list_path, 'w') as hmm_list:
            columns = ['name', 'source', 'path']
            values = [hmm_name, source_name, hmm_path]
            if group is not None:
                columns.append('group')
                values.append(group)
            hmm_list.write('\t'.join(columns) + '\n')
            hmm_list.write('\t'.join(values) + '\n')

        workflow = EcoPhyloWorkflow.__new__(EcoPhyloWorkflow)
        workflow.hmm_list_path = hmm_list_path
        workflow.run_scg_taxonomy = False
        workflow.init_hmm_list_txt()
        return workflow

    def test_external_hmm_uses_default_and_explicit_groups(self):
        for group in (None, 'ribosomal_group'):
            with self.subTest(group=group):
                workflow = self.make_workflow('Ribosomal_L16', group)
                hmm_entry = workflow.hmm_dict['Ribosomal_L16_external_hmm_Ribosomal_L16']
                self.assertEqual(
                    hmm_entry['group'],
                    group if group is not None else 'Ribosomal_L16_external_hmm_Ribosomal_L16',
                )

    def test_external_hmm_wrong_gene_name_is_config_error(self):
        with self.assertRaisesRegex(ConfigError, r'change the gene name Wrong_gene to\s+this: Ribosomal_L16'):
            self.make_workflow('Wrong_gene')


if __name__ == '__main__':
    unittest.main()
