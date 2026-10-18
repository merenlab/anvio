""" Classes to define and work with anvi'o SRA_downloads workflow. """

import os
import anvio
import hashlib
import pandas as pd

import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.workflows import WorkflowSuperClass


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = ['mschecht']
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Matthew S. Schechter"
__email__ = "mschechter@uchicago.edu"


run = terminal.Run()

class SRADownloadWorkflow(WorkflowSuperClass):
    def __init__(self, args=None, run=terminal.Run(), progress=terminal.Progress()):
        self.init_workflow_super_class(args, workflow_name='sra_download')

        # check that NCBI SRA Toolkit and other programs are installed
        NCBI_sra_tool_programs = ['prefetch', 'fasterq-dump']
        other_programs = ['pigz']

        for program in NCBI_sra_tool_programs:
            if not utils.is_program_exists(program, dont_raise=True):
                raise ConfigError(f"The program {program} is not installed in your anvi'o conda environment. "
                                  f"'prefetch' and 'fasterq-dump' are from the NCBI SRA toolkit and must be installed for the "
                                  f"sra_download workflow to work. Please check out the installation instructions here: "
                                  f"https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit")
        for program in other_programs:
            if not utils.is_program_exists(program, dont_raise=True):
                raise ConfigError(f"The program {program} is not installed in your anvi'o conda environment. Please "
                                  f"double check you installed all of the programs listed in the anvi'o installation tutorial: https://anvio.org/install/")

        # Snakemake rules
        self.rules.extend(['prefetch',
                           'fasterq_dump',
                           'check_md5sum',
                           'pigz'])


        # Directory structure for Snakemake workflow
        self.dirs_dict.update({"SRA_prefetch": "01_NCBI_SRA"})
        self.dirs_dict.update({"FASTQ_DIR": "02_FASTQ"})


    def init(self):
        """Load the SRA_accession_list and creates a list of target files for the Snakemake workflow."""

        self.migrate_legacy_fastas_dir_name()

        super().init()

        self.warn_about_legacy_fastq_directory()

        # Load SRA_accession_list
        self.SRA_accession_list = self.get_param_value_from_config(['SRA_accession_list'])

        if not self.SRA_accession_list:
            raise ConfigError('Please provide a list of SRA accessions')

        if self.SRA_accession_list:
            filesnpaths.is_file_exists(self.SRA_accession_list)
            try:
                self.SRA_accession_list_df = pd.read_csv(self.SRA_accession_list, sep='\t', index_col=False, header=None, names=['accessions'])
                self.accessions_list = self.SRA_accession_list_df['accessions'].tolist()
            except IndexError as e:
                raise ConfigError(f"Looks like your SRA accession list file, {self.SRA_accession_list}, is not properly formatted. "
                                  f"This is what we know: {e}")

        for accession in self.accessions_list:
            if not accession.startswith(('SRR', 'ERR', 'DRR')):
                if accession.startswith('SAMEA'):
                    raise ConfigError(f"anvi'o found an NCBI BioSample in your {self.SRA_accession_list}: {accession}. "
                                      f"The anvi'o sra-download workflow only processes sequencing accessions that start with the prefix: ERR, SRR, or DRR. "
                                      f"Search for the BioSample accession '{accession}' on the [NCBI SRA website](https://www.ncbi.nlm.nih.gov/sra) "
                                      f"and find the sequencing accessions.")
                else:
                    raise ConfigError(f"Looks like one of your \"SRA accessions\", {accession}, is not an SRA accession :( "
                                      f"anvi'o asks that you kindly double check your SRA_accession_list.txt ({self.SRA_accession_list}) to confirm you "
                                      f"are using the correct accessions. Hint: SRA accessions start with the prefix: ERR, SRR, or DRR")

        self.target_files = self.get_target_files()


    def migrate_legacy_fastas_dir_name(self):
        """Accept the old `FASTAS` key in a config file's `output_dirs` section.

        The directory holding the FASTQ files this workflow downloads used to be declared with
        the key `FASTAS`. Config files written before the rename still use it, and the sanity
        check for `output_dirs` would reject them as an unknown directory, so we quietly
        translate the old key into the current one and let the user know."""

        output_dirs = (self.config or {}).get('output_dirs') or {}

        if 'FASTAS' not in output_dirs:
            return

        legacy_value = output_dirs.pop('FASTAS')
        output_dirs.setdefault('FASTQ_DIR', legacy_value)

        self.run.warning(f"The `output_dirs` section of your config file asks for a directory called "
                         f"`FASTAS`. That key is now called `FASTQ_DIR`, since what this workflow puts "
                         f"there are FASTQ files rather than FASTA files. Anvi'o will use your value "
                         f"('{legacy_value}') for `FASTQ_DIR` and carry on, but please do rename the key "
                         f"in your config file when you have a moment.")


    def warn_about_legacy_fastq_directory(self):
        """Point out output from an earlier run that lives under the old directory name.

        Snakemake decides what work is left to do by looking for output files, so a run that
        downloaded its reads into `02_FASTA` before the rename would look entirely unfinished
        now and every accession would be downloaded a second time. Downloads are the expensive
        part of this workflow, so it is worth saying something before that happens."""

        fastq_dir = self.dirs_dict['FASTQ_DIR']
        legacy_dir = '02_FASTA'

        if fastq_dir == legacy_dir or not os.path.exists(legacy_dir) or os.path.exists(fastq_dir):
            return

        self.run.warning(f"Anvi'o found a directory called '{legacy_dir}' here, which is where this workflow "
                         f"used to keep the FASTQ files it downloads. It now uses '{fastq_dir}' instead. If "
                         f"'{legacy_dir}' holds reads from an earlier run, anvi'o will not see them and will "
                         f"download everything again from scratch. To avoid that, either rename '{legacy_dir}' "
                         f"to '{fastq_dir}', or add \"output_dirs\": {{\"FASTQ_DIR\": \"{legacy_dir}\"}} to your "
                         f"config file to keep using the old location.")


    def get_target_files(self):
        """Get list of target files for snakemake target rule"""

        target_files = [os.path.join(self.dirs_dict['FASTQ_DIR'], "generate_samples_txt.done")]

        return target_files


    def calculate_md5(self, file_path):
        """Calculate the md5sum of a file"""
        hash_md5 = hashlib.md5(usedforsecurity=False)
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(65536), b""):
                hash_md5.update(chunk)
        return hash_md5.hexdigest()
