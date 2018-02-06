import os
import anvio

from snakemake.shell import shell

__author__ = "Alon Shaiber"
__copyright__ = "Copyright 2017, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Alon Shaiber"
__email__ = "alon.shaiber@gmail.com"


shell(
    "anvi-script-reformat-fasta {snakemake.input.contigs} -o {snakemake.output.contigs} -r"
    "{snakemake.output.report} --simplify-names --prefix {snakemake.params.prefix} &>> {snakemake.log}")
