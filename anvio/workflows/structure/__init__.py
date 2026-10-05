"""
    StructureModule — a mixin providing protein structure prediction steps for anvi'o workflows.
"""

import os
import shlex

from snakemake.io import ancient

import anvio
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.workflows import WorkflowSuperClass


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__


class StructureModule(WorkflowSuperClass):
    """Mixin providing structure database generation with `anvi-gen-structure-database`.

    A structure database is predicted either in a single step (rule `anvi_gen_structure_database`,
    for MODELLER or ColabFold), or, for ColabFold with a local database, in two steps that split the
    CPU-heavy MSA generation (rule `colabfold_msa`, `--only-msa`) from the GPU-heavy prediction (rule
    `colabfold_predict`, `--only-predict`), connected by a `--dump-dir` checkpoint. Every prediction
    parameter lives under `anvi_gen_structure_database` in the config, so both modes read them from the
    same place; `colabfold_msa` and `colabfold_predict` only carry their own resources.

    Rules that need a GPU declare it as the `gpu` resource, which users can map onto their own scheduler.

    Subclasses must set before calling `init_structure()`:
      - self.dirs_dict: includes "STRUCTURE_DIR"
      - self.fasta_information: dict[str, dict] (may be empty)
      - self.get_contigs_db_path(): returns the contigs database path with a `{group}` wildcard
    """

    valid_structure_engines = ['modeller', 'colabfold']

    # parameters of `anvi_gen_structure_database` that are passed verbatim to the program (see
    # get_structure_engine_params). `--colabfold-additional-parameters` is not in this list because its
    # value starts with dashes and must be attached with an equals sign
    structure_engine_params = ['--num-models', '--skip-DSSP',
                               '--deviation', '--modeller-database', '--scoring-method', '--very-fast',
                               '--percent-cutoff', '--alignment-fraction-cutoff', '--max-number-templates',
                               '--modeller-executable', '--offline-mode', '--pdb-db',
                               '--colabfold-conda-env', '--num-recycle', '--amber', '--write-buffer-size']

    def __init__(self):
        # the accepted parameters, types, and default values of these rules are declared in
        # anvio/workflows/structure/params.json, which WorkflowSuperClass picks up for any workflow
        # that inherits this module
        self.rules.extend(['anvi_gen_structure_database', 'colabfold_msa', 'colabfold_predict'])

        self.run_structure = False
        self.structure_engine = 'modeller'
        self.split_msa_and_predict = False


    def init_structure(self):
        """Read and sanity check the structure prediction configuration."""
        self.run_structure = self.get_param_value_from_config(['anvi_gen_structure_database', 'run']) == True
        self.structure_engine = self.get_param_value_from_config(['anvi_gen_structure_database', '--engine']) or 'modeller'
        self.split_msa_and_predict = self.get_param_value_from_config(['anvi_gen_structure_database', 'split_msa_and_predict']) == True

        if not self.run_structure:
            return

        if self.structure_engine not in self.valid_structure_engines:
            raise ConfigError(f"The `--engine` for `anvi_gen_structure_database` must be one of "
                              f"{', '.join(self.valid_structure_engines)}, but your config says '{self.structure_engine}'.")

        colabfold_db = self.get_param_value_from_config(['anvi_gen_structure_database', '--colabfold-db'])
        colabfold_msa_server = self.get_param_value_from_config(['anvi_gen_structure_database', '--colabfold-msa-server']) == True

        if self.structure_engine == 'colabfold':
            if colabfold_db and colabfold_msa_server:
                raise ConfigError("Your config for `anvi_gen_structure_database` asks ColabFold to use both a local "
                                  "database (`--colabfold-db`) and the public MSA server (`--colabfold-msa-server`). "
                                  "Please pick only one.")
            if not colabfold_db and not colabfold_msa_server:
                raise ConfigError("ColabFold needs to know how to generate the multiple sequence alignments. Please "
                                  "set either `--colabfold-db` (a local ColabFold database, recommended for many "
                                  "sequences) or `--colabfold-msa-server` (the public MMseqs2 server) for "
                                  "`anvi_gen_structure_database` in your config.")
        elif colabfold_db or colabfold_msa_server:
            raise ConfigError(f"Your config for `anvi_gen_structure_database` includes ColabFold parameters, but the "
                              f"`--engine` is '{self.structure_engine}'. Please set `--engine` to 'colabfold', or "
                              f"drop the ColabFold parameters.")

        if self.split_msa_and_predict:
            if self.structure_engine != 'colabfold':
                raise ConfigError("`split_msa_and_predict` splits ColabFold's MSA and prediction steps into two jobs, "
                                  "so it only works with `--engine colabfold`.")
            if not colabfold_db:
                raise ConfigError("`split_msa_and_predict` only works with a local ColabFold database "
                                  "(`--colabfold-db`): the public MSA server generates the MSA and predicts the "
                                  "structure in a single step that cannot be split.")

        for group, info in self.fasta_information.items():
            genes_of_interest = info.get('structure_genes_of_interest')
            if genes_of_interest and not filesnpaths.is_file_exists(genes_of_interest, dont_raise=True):
                raise ConfigError(f"The `structure_genes_of_interest` file for '{group}' in your fasta_txt does not "
                                  f"exist: '{genes_of_interest}'.")

        self.warn_about_structure_gene_selection()


    def warn_about_structure_gene_selection(self):
        """Let the user know which groups will get a structure for every single gene."""
        groups_with_all_genes = [g for g in self.get_structure_group_names()
                                 if not self.fasta_information.get(g, {}).get('structure_genes_of_interest')]

        if not groups_with_all_genes:
            return

        self.run.warning(f"Anvi'o will predict a structure for EVERY gene in the contigs databases of the following "
                         f"(it found no `structure_genes_of_interest` for them in your fasta_txt): "
                         f"{', '.join(groups_with_all_genes)}. A single genome has thousands of genes, and a "
                         f"metagenomic assembly can have hundreds of thousands, so this may take a very long time "
                         f"(days, if not weeks, with ColabFold). If that is not what you want, you can let the "
                         f"workflow build your contigs databases first, pick your genes of interest, add them to "
                         f"your fasta_txt under the column `structure_genes_of_interest`, and then run the "
                         f"workflow again: it will only run the missing structure steps.",
                         header="STRUCTURE PREDICTION FOR ALL GENES", lc='yellow')


    def get_structure_group_names(self):
        """The groups for which a structure database is generated."""
        return self.group_names


    def get_structure_db_path(self, group='{group}'):
        return os.path.join(self.dirs_dict["STRUCTURE_DIR"], f"{group}-STRUCTURE.db")


    def get_colabfold_msa_dir(self, group='{group}'):
        return os.path.join(self.dirs_dict["STRUCTURE_DIR"], f"{group}-COLABFOLD-MSA")


    def get_structure_target_files(self):
        if not self.run_structure:
            return []

        return [self.get_structure_db_path(group) for group in self.get_structure_group_names()]


    def get_input_for_structure_rules(self, wildcards):
        # the contigs database is ancient(): annotation rules keep writing to it after the structure
        # database exists, and that must not trigger another (potentially days-long) prediction
        d = {}
        d['contigs_db'] = ancient(self.get_contigs_db_path().format(group=wildcards.group))

        genes_of_interest = self.fasta_information.get(wildcards.group, {}).get('structure_genes_of_interest')
        if genes_of_interest:
            d['genes_of_interest'] = genes_of_interest

        return d


    def get_structure_genes_of_interest_param(self, wildcards, input):
        if 'genes_of_interest' in input.keys():
            return f"--genes-of-interest {input.genes_of_interest}"
        return ''


    def get_structure_engine_params(self, msa_source=True):
        """The prediction parameters shared by every structure rule, as a command line string.

        `msa_source` is False for `--only-predict`, which reads its MSAs from the checkpoint and takes
        no MSA source."""
        rule = 'anvi_gen_structure_database'

        params = [f"--engine {self.structure_engine}"]
        params.extend([self.get_rule_param(rule, p) for p in self.structure_engine_params])

        if msa_source:
            params.append(self.get_rule_param(rule, '--colabfold-db'))
            params.append(self.get_rule_param(rule, '--colabfold-msa-server'))

        additional_parameters = self.get_param_value_from_config([rule, '--colabfold-additional-parameters'])
        if additional_parameters:
            params.append('--colabfold-additional-parameters=' + shlex.quote(additional_parameters))

        return ' '.join(p for p in params if p)


    def get_structure_gpu_resource(self, rule):
        """The number of GPUs a structure rule declares as its `gpu` resource.

        Rules that predict structures with ColabFold need a GPU; the MSA step and MODELLER do not. Users
        can override it with the `gpu` parameter of a rule."""
        gpu = self.get_param_value_from_config([rule, 'gpu'])
        if gpu not in [None, '']:
            return int(gpu)

        if rule == 'colabfold_msa':
            return 0
        elif rule == 'colabfold_predict':
            return 1
        else:
            return 1 if self.structure_engine == 'colabfold' else 0
