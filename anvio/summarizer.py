# coding: utf-8
# pylint: disable=line-too-long
"""Summarizes information for a collection.

It also gives access to bin data that may be useful. For instance, did you know
could totally do this to access to gene coverage and detection dicts for a given
bin in a given collection:

import anvio.summarizer as summarizer

    >>> class Args: None
    >>> args = Args()
    >>> args.profile_db = 'SAMPLES-MERGED/PROFILE.db'
    >>> args.contigs_db = 'CONTIGS.db'
    >>> args.collection_name = "CONCOCT"

    >>> summary = summarizer.ProfileSummarizer(args)
    >>> summary.init()
    >>> _bin = summarizer.Bin(summary, 'Bin_1')

    # now you have these:
    >>> _bin.gene_coverages
    >>> _bin.gene_detection

"""

import os
import sys
import gzip
import numpy
import mistune
import argparse
import textwrap
import pandas as pd

from collections import Counter

import anvio
import anvio.utils as utils
import anvio.hmmops as hmmops
import anvio.sequence as seqlib
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.filesnpaths as filesnpaths
import anvio.ccollections as ccollections
import anvio.completeness as completeness
import anvio.taxonomyops.scg as scgtaxonomyops

from anvio.errors import ConfigError
from anvio.dbops import DatabasesMetaclass, ContigsSuperclass, PanSuperclass
from anvio.hmmops import SequencesForHMMHits
from anvio.summaryhtml import SummaryHTMLOutput, humanize_n, pretty
from anvio.tables.miscdata import TableForLayerAdditionalData, MiscDataTableFactory


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


pp = terminal.pretty_print
run = terminal.Run()
progress = terminal.Progress()
progress_quiet = terminal.Progress(verbose=False)
P = lambda x, y: float(x) * 100 / float(y)


class ArgsTemplateForSummarizerClass:
    """Just a dummy args template for ad hoc use of summary classes.

    You can use it like this:

        >>> args = summarizer.ArgsTemplateForSummarizerClass()
        >>> args.profile_db = profile_db_path
        >>> args.contigs_db = contigs_db_path
        >>> args.collection_name = collection_name
        >>> args.output_dir = output_dir

        >>> summarizer.ProfileSummarizer(args)
    """

    def __init__(self):
        self.profile_db = None
        self.pan_db = None
        self.contigs_db = None
        self.collection_name = None
        self.list_collections = None
        self.list_bins = None
        self.debug = None
        self.quick_summary = False
        self.init_gene_coverages = False
        self.calculate_Q2Q3_carefully = False
        self.skip_check_collection_name = False
        self.skip_init_functions = False
        self.cog_data_dir = None
        self.output_dir = filesnpaths.get_temp_directory_path()
        self.report_aa_seqs_for_gene_calls = False
        self.reformat_contig_names = False


class SummarizerSuperClass(object):
    def __init__(self, args, r=run, p=progress):
        self.summary = {}
        self.collection_name = None

        self.collections = ccollections.Collections()

        if self.summary_type == "pan":
            self.collections.populate_collections_dict(self.pan_db_path)
        else:
            self.collections.populate_collections_dict(
                self.pan_db_path if self.summary_type == "pan" else self.profile_db_path
            )
            (
                self.collections.populate_collections_dict(self.contigs_db_path)
                if self.contigs_db_path
                else None
            )

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None

        if A("list_collections"):
            self.collections.list_collections()
            sys.exit()

        self.collection_name = A("collection_name")

        if A("list_bins"):
            if not self.collection_name:
                raise ConfigError(
                    "It may come across as a surprise, but you really can't list bins in a collection "
                    "without providing a collection name. Bioinformatics. Never understands what you "
                    "need and all :/"
                )
            self.collections.list_bins_in_collection(
                collection_name=self.collection_name
            )
            sys.exit()

        self.skip_check_collection_name = A("skip_check_collection_name")
        self.skip_init_functions = A("skip_init_functions")
        self.init_gene_coverages = A("init_gene_coverages")
        self.calculate_Q2Q3_carefully = A("calculate_Q2Q3_carefully")
        self.output_directory = A("output_dir")
        self.quick = A("quick_summary")
        self.debug = A("debug")
        self.cog_data_dir = A("cog_data_dir")
        self.report_aa_seqs_for_gene_calls = A("report_aa_seqs_for_gene_calls")
        self.report_DNA_sequences = A("report_DNA_sequences")
        self.delete_output_directory_if_exists = (
            False
            if A("delete_output_directory_if_exists") == None
            else A("delete_output_directory_if_exists")
        )
        self.just_do_it = A("just_do_it")
        self.reformat_contig_names = A("reformat_contig_names")

        if not self.lazy_init:
            self.sanity_check()

        if self.output_directory:
            self.output_directory = filesnpaths.check_output_directory(
                self.output_directory,
                ok_if_exists=self.delete_output_directory_if_exists or self.just_do_it,
            )
            filesnpaths.gen_output_directory(
                self.output_directory,
                delete_if_exists=self.delete_output_directory_if_exists
                or self.just_do_it,
            )
        else:
            self.output_directory = "SUMMARY"

    def report_misc_data_files(self, target_table="layers"):
        if target_table == "layer_orders":
            raise ConfigError(
                "Report misc data files do not know how to work with layer orders yet :/"
            )

        run_obj = terminal.Run(verbose=False)

        db_path = (
            self.pan_db_path if self.summary_type == "pan" else self.profile_db_path
        )
        additional_data = MiscDataTableFactory(
            argparse.Namespace(
                pan_or_profile_db=db_path, target_data_table=target_table
            ),
            r=run_obj,
        )

        data_groups, data_dict = additional_data.get_all()
        for data_group in data_groups:
            output_file_obj = self.get_output_file_handle(
                sub_directory="misc_data_%s" % target_table,
                prefix="%s.txt" % data_group,
            )
            output_file_path = output_file_obj.name
            output_file_obj.close()
            MiscDataTableFactory(
                argparse.Namespace(
                    pan_or_profile_db=db_path,
                    target_data_table=target_table,
                    target_data_group=data_group,
                ),
                r=run_obj,
            ).export(output_file_path=output_file_path)

        if "misc_data" not in self.summary:
            self.summary["misc_data"] = {}

        data_groups_reported = list(data_groups.keys())
        self.summary["misc_data"][target_table] = data_groups
        self.run.info(
            "Misc data reported for %s" % target_table,
            ", ".join(data_groups_reported) if data_groups else "None",
            nl_before=1,
        )

    def sanity_check(self):
        if not self.skip_check_collection_name:
            if not self.collection_name:
                raise ConfigError("You must specify a collection id :/")

            if self.collection_name not in self.collections.collections_dict:
                raise ConfigError(
                    "%s is not a valid collection ID. See a list of available ones with '--list-collections' flag"
                    % self.collection_name
                )

    def get_output_file_handle(
        self,
        sub_directory=None,
        prefix="output.txt",
        overwrite=False,
        within=None,
        compress_output=False,
        add_project_name=False,
    ):
        if sub_directory:
            output_directory = os.path.join(self.output_directory, sub_directory)
        else:
            output_directory = self.output_directory

        if not os.path.exists(output_directory):
            filesnpaths.gen_output_directory(output_directory)

        key = prefix.split(".")[0].replace("-", "_")

        if add_project_name:
            prefix = "%s_%s" % (self.p_meta["project_name"].replace(" ", "_"), prefix)

        if within:
            file_path = os.path.join(output_directory, "%s_%s" % (within, prefix))
        else:
            file_path = os.path.join(output_directory, "%s" % (prefix))

        if compress_output:
            file_path += ".gz"

        if os.path.exists(file_path) and not overwrite:
            raise ConfigError(
                'get_output_file_handle: well, this file already exists: "%s"'
                % file_path
            )

        if within:
            if within not in self.summary["files"]:
                self.summary["files"][within] = {}
            self.summary["files"][within][key] = file_path[
                len(self.output_directory) :
            ].strip("/")
        else:
            self.summary["files"][key] = file_path[len(self.output_directory) :].strip(
                "/"
            )

        if compress_output:
            return gzip.open(file_path, "wb")
        else:
            return open(file_path, "w")


class PanSummarizer(PanSuperclass, SummarizerSuperClass):
    """Creates a dictionary of summary for anvi'o pan profiles"""

    def __init__(self, args=None, lazy_init=False, r=run, p=progress):
        self.summary_type = "pan"
        self.debug = False
        self.quick = False
        self.pan_db_path = None
        self.lazy_init = lazy_init
        self.output_directory = None
        self.genomes_storage_path = None

        self.run = r
        self.progress = p

        PanSuperclass.__init__(self, args, self.run, self.progress)
        if not self.genomes_storage_is_available:
            raise ConfigError("No genomes storage no summary. Yes. Very simple stuff.")

        SummarizerSuperClass.__init__(self, args, self.run, self.progress)

        # init gene clusters and functions from Pan super.
        self.init_gene_clusters()

        # init items additional data.
        self.init_items_additional_data()

        if not self.skip_init_functions:
            self.init_gene_clusters_functions()

        self.collection_dict, self.bins_info_dict, self.bin_ids = None, None, None
        if self.collection_name:
            self.collection_dict, self.bins_info_dict = self.init_collection_profile(
                self.collection_name
            )
            self.bin_ids = sorted(self.collection_dict.keys())

        # see if COG functions or categories are available
        self.cog_functions_are_called = (
            "COG_FUNCTION" in self.gene_clusters_function_sources
        )
        self.cog_categories_are_called = (
            "COG_CATEGORY" in self.gene_clusters_function_sources
        )

    def get_occurrence_of_functions_in_pangenome(
        self, gene_clusters_functions_summary_dict
    ):
        """
        For each function we will create a fake merged gc, with an occurrence vector
        which is the "or" product of all gcs that match this function.

        We also keep track of all the GCs annotated with the function.
        """
        occurrence_of_functions_in_pangenome_dict = {}

        self.progress.new("Computing presence absence for functions in genomes")
        self.progress.update("Creating a dictionary")

        for gene_cluster_id in gene_clusters_functions_summary_dict:
            gene_cluster_function = gene_clusters_functions_summary_dict[
                gene_cluster_id
            ]["gene_cluster_function"]
            gene_cluster_function_accession = gene_clusters_functions_summary_dict[
                gene_cluster_id
            ]["gene_cluster_function_accession"]
            if gene_cluster_function:
                if (
                    gene_cluster_function
                    not in occurrence_of_functions_in_pangenome_dict
                ):
                    occurrence_of_functions_in_pangenome_dict[gene_cluster_function] = (
                        {}
                    )
                    occurrence_of_functions_in_pangenome_dict[gene_cluster_function][
                        "gene_clusters_ids"
                    ] = []
                    occurrence_of_functions_in_pangenome_dict[gene_cluster_function][
                        "occurrence"
                    ] = None
                    occurrence_of_functions_in_pangenome_dict[gene_cluster_function][
                        "accession"
                    ] = gene_cluster_function_accession
                occurrence_of_functions_in_pangenome_dict[gene_cluster_function][
                    "gene_clusters_ids"
                ].append(gene_cluster_id)

        from anvio.dbops import PanDatabase

        pan_db = PanDatabase(self.pan_db_path)

        view_data, _ = pan_db.db.get_view_data("gene_cluster_frequencies")
        gene_cluster_frequencies_dataframe = pd.DataFrame.from_dict(
            view_data, orient="index"
        )

        self.progress.update(
            "Merging presence/absence of gene clusters with the same function"
        )

        D = {}
        for gene_cluster_function in occurrence_of_functions_in_pangenome_dict:
            v = None
            for gene_cluster_id in occurrence_of_functions_in_pangenome_dict[
                gene_cluster_function
            ]["gene_clusters_ids"]:
                if v is None:
                    v = gene_cluster_frequencies_dataframe.loc[gene_cluster_id,].astype(
                        int
                    )
                else:
                    v = v.add(gene_cluster_frequencies_dataframe.loc[gene_cluster_id,])
            D[gene_cluster_function] = {}
            for genome in v.index:
                D[gene_cluster_function][genome] = v[genome]

        self.progress.end()

        return pd.DataFrame.from_dict(D), occurrence_of_functions_in_pangenome_dict

    def functional_enrichment_stats(self):
        """
        Compute functional enrichment.

        To learn more refer to the documentation:
            anvi-compute-functional-enrichment-in-pan -h
            anvi-compute-functional-enrichment-across-genomes -h
            anvi-compute-metabolic-enrichment -h
        """
        A = lambda x: self.args.__dict__[x] if x in self.args.__dict__ else None
        output_file_path = A("output_file")
        tmp_functional_occurrence_file = filesnpaths.get_temp_file_path()

        enrichment_file_path = output_file_path
        if not enrichment_file_path:
            # if no output was requested it means a programmer is calling this function
            # in that case, we will use a tmp file for the enrichment output
            enrichment_file_path = filesnpaths.get_temp_file_path()

        # set the tmp output_file for the functional occurrence.
        # this is a little hacky, but allows for the functional occurrence to be used without
        # the functional enrichment (if we will want that at some point)
        self.args.output_file = tmp_functional_occurrence_file

        # Call functional occurrence. Output is saved to the tmp file
        self.functional_occurrence_stats()

        # run enrichment script. this uses the output saved from the previous step.
        enrichment_stats = utils.run_functional_enrichment_stats(
            tmp_functional_occurrence_file,
            output_file_path,
            run=self.run,
            progress=self.progress,
        )

        return enrichment_stats

    def functional_occurrence_stats(self):
        """
        Compute the functional occurrence stats for a pangenome

        returns the functional_occurrence_summary_dict, which has the form:
            value = functional_occurrence_summary_dict[category_name][function_name][value_name]

        If self.args.output_file_path exists then an output file is created.

        To learn more about how this works refer to the documentation:
            anvi-get-functional-occurrence-summary-per-pan-group -h
        """

        A = lambda x: self.args.__dict__[x] if x in self.args.__dict__ else None
        output_file_path = A("output_file")
        category_variable = A("category_variable")
        functional_annotation_source = A("annotation_source")
        list_functional_annotation_sources = A("list_annotation_sources")
        functional_occurrence_table_output = A("functional_occurrence_table_output")
        include_ungrouped = A("include_ungrouped")

        if output_file_path:
            filesnpaths.is_output_file_writable(output_file_path)

        if functional_occurrence_table_output:
            filesnpaths.is_output_file_writable(functional_occurrence_table_output)

        if not self.functions_initialized:
            raise ConfigError(
                "For some reason funtions are not initialized for this pan class instance. We "
                "can't summarize functional occurrence stats without that :/"
            )

        if not len(self.gene_clusters_functions_dict):
            raise ConfigError(
                "The gene clusters functions dict seems to be empty. We assume this error makes "
                "zero sense to you, and it probably will not help you to know that it also makes "
                "zero sense to anvi'o too :/ Maybe you forgot to provide a genomes storage?"
            )

        if list_functional_annotation_sources:
            self.run.info(
                "Available functional annotation sources",
                ", ".join(self.gene_clusters_function_sources),
            )
            sys.exit()

        if not category_variable:
            raise ConfigError(
                "For this to work, you must provide a category variable .. and it better be in "
                "the misc additional layer data table, too. If you don't have any idea what is "
                "available, try `anvi-show-misc-data`."
            )

        if not functional_annotation_source:
            raise ConfigError(
                "You haven't provided a functional annotation source to make sense of functional "
                "occurrence stats as defined by the categorical variable %s. These are the functions "
                "that are available, so pick one: %s."
                % (category_variable, ", ".join(self.gene_clusters_function_sources))
            )

        if functional_annotation_source not in self.gene_clusters_function_sources:
            raise ConfigError(
                "Your favorite functional annotation source '%s' does not seem to be among one of the sources "
                "that are available to you. Here are the ones you should choose from: %s."
                % (
                    functional_annotation_source,
                    ", ".join(self.gene_clusters_function_sources),
                )
            )

        keys, categories_dict = TableForLayerAdditionalData(
            argparse.Namespace(pan_db=self.pan_db_path)
        ).get(additional_data_keys_requested=[category_variable])

        values_that_are_not_none = [
            s for s in categories_dict.values() if s[category_variable] is not None
        ]
        if not values_that_are_not_none:
            raise ConfigError(
                "The variable '%s' contains only values of type None,\
                               this is probably a mistake, surely you didn't mean to provide an empty category.\
                               Do you think this is a mistake on our part? Let us know."
                % category_variable
            )

        type_category_variable = type(values_that_are_not_none[0][category_variable])
        if type_category_variable != str:
            raise ConfigError(
                "The variable '%s' does not seem to resemble anything that could be a category. "
                "Anvi'o expects these variables to be of type string, yet yours is type %s :/ "
                "Do you think this is a mistake on our part? Let us know."
                % (category_variable, type_category_variable)
            )

        gene_clusters_functions_summary_dict = (
            self.get_gene_clusters_functions_summary_dict(functional_annotation_source)
        )

        self.run.info("Category", category_variable)
        self.run.info("Functional annotation source", functional_annotation_source)
        self.run.info("Include ungrouped", include_ungrouped)

        (
            occurrence_frequency_of_functions_in_pangenome_dataframe,
            occurrence_of_functions_in_pangenome_dict,
        ) = self.get_occurrence_of_functions_in_pangenome(
            gene_clusters_functions_summary_dict
        )

        if functional_occurrence_table_output:
            occurrence_frequency_of_functions_in_pangenome_dataframe.astype(
                int
            ).transpose().to_csv(functional_occurrence_table_output, sep="\t")
            self.run.info(
                "Occurrence frequency of functions:", functional_occurrence_table_output
            )

        # Get the presence/absence info for functions, which we will use for the comparisson between groups
        occurrence_of_functions_in_pangenome_dataframe = (
            occurrence_frequency_of_functions_in_pangenome_dataframe.astype(
                bool
            ).astype(int)
        )

        self.progress.new("Functional occurrence analysis")
        self.progress.update("Creating a dictionary")

        # get a list of unique function names
        functions_names = set(occurrence_of_functions_in_pangenome_dataframe.columns)

        # add a category column to the dataframe
        occurrence_of_functions_in_pangenome_dataframe["category"] = (
            occurrence_of_functions_in_pangenome_dataframe.index.map(
                lambda x: str(categories_dict[x][category_variable])
            )
        )

        # the sum of occurrences of each function in each category
        functions_in_categories = (
            occurrence_of_functions_in_pangenome_dataframe.groupby("category").sum()
        )

        # unique names of categories
        categories = set(
            [
                str(categories_dict[g][category_variable])
                for g in categories_dict.keys()
                if (
                    categories_dict[g][category_variable] is not None
                    or include_ungrouped
                )
            ]
        )

        categories_to_genomes_dict = {}
        for c in categories:
            categories_to_genomes_dict[c] = set(
                [
                    genome
                    for genome in categories_dict.keys()
                    if str(categories_dict[genome][category_variable]) == c
                ]
            )

        functional_occurrence_summary_dict = {}
        function_occurrence_table = {}

        # populate the number of genomes per category once
        categories_few_genomes = []
        for c in categories:
            function_occurrence_table[c] = {}
            function_occurrence_table[c]["N"] = len(categories_to_genomes_dict[c])
            if function_occurrence_table[c]["N"] < 8:
                categories_few_genomes.append(c)

        # warn user if they have a low number of genomes per group
        if categories_few_genomes:
            self.progress.reset()
            categories_string = ", ".join(categories_few_genomes)
            self.run.warning(
                "Some of your groups have very few genomes in them, so if you are running functional enrichment, the statistical test may not be very reliable. "
                "The minimal number of genomes in a group for the test to be reliable depends on a number of factors, "
                "but we recommend proceeding with great caution because the following groups have fewer than 8 genomes: "
                f"{categories_string}."
            )

        self.progress.update(
            "Generating the input table for functional enrichment analysis"
        )
        for f in functions_names:
            for c in categories:
                function_occurrence_table[c]["p"] = (
                    functions_in_categories.loc[c, f]
                    / function_occurrence_table[c]["N"]
                )
            function_occurrence_table_df = pd.DataFrame.from_dict(
                function_occurrence_table, orient="index"
            )

            functional_occurrence_summary_dict[f] = {}
            functional_occurrence_summary_dict[f]["gene_clusters_ids"] = (
                occurrence_of_functions_in_pangenome_dict[f]["gene_clusters_ids"]
            )
            functional_occurrence_summary_dict[f]["accession"] = (
                occurrence_of_functions_in_pangenome_dict[f]["accession"]
            )
            for c in categories:
                functional_occurrence_summary_dict[f]["p_" + c] = (
                    function_occurrence_table[c]["p"]
                )
            for c in categories:
                functional_occurrence_summary_dict[f]["N_" + c] = (
                    function_occurrence_table[c]["N"]
                )
            enriched_groups_vector = utils.get_enriched_groups(
                function_occurrence_table_df["p"].values,
                function_occurrence_table_df["N"].values,
            )
            c_dict = dict(
                zip(
                    function_occurrence_table_df["p"].index,
                    range(len(function_occurrence_table_df["p"].index)),
                )
            )
            associated_groups = [
                c for c in categories if enriched_groups_vector[c_dict[c]]
            ]
            functional_occurrence_summary_dict[f][
                "associated_groups"
            ] = associated_groups

        if output_file_path:
            self.progress.update("Generating the output file")
            functional_occurrence_summary_data_frame = (
                self.get_functional_occurrence_summary_dict_as_dataframe(
                    functional_occurrence_summary_dict, functional_annotation_source
                )
            )

            # Sort the columns the way we want them
            columns = [
                functional_annotation_source,
                "accession",
                "gene_clusters_ids",
                "associated_groups",
            ]
            columns.extend([s + c for s in ["p_", "N_"] for c in categories])
            functional_occurrence_summary_data_frame.to_csv(
                output_file_path,
                sep="\t",
                index=False,
                float_format="%.4f",
                columns=columns,
            )

        self.progress.end()

        if output_file_path:
            self.run.info("Functional occurrence summary", output_file_path)

        return functional_occurrence_summary_dict

    def get_functional_occurrence_summary_dict_as_dataframe(
        self, functional_occurrence_summary_dict, functional_annotation_source
    ):
        # convert dictionary to pandas
        # we want to deal with sequences (see below) and so this is easier than using pandas from_dict
        # we first convert it to a dict of dicts and then convert to pandas
        # because this is faster than alternatives
        i = 0
        D = {}
        for f in functional_occurrence_summary_dict:
            D[i] = {}
            D[i][functional_annotation_source] = f
            for key, value in functional_occurrence_summary_dict[f].items():
                if type(value) == str:
                    D[i][key] = value
                else:
                    try:
                        # if there is a sequence of values
                        # merge them with commas for nicer printing
                        D[i][key] = ", ".join(iter(value))
                    except:
                        # if it is not a sequnce, it is a single value
                        D[i][key] = value
            i += 1
        functional_occurrence_summary_data_frame = pd.DataFrame.from_dict(
            D, orient="index"
        )

        return functional_occurrence_summary_data_frame

    def process(self):
        # let bin names known to all
        bin_ids = list(self.collection_profile.keys())

        genome_names = ", ".join(list(self.gene_clusters.values())[0].keys())

        # set up the initial summary dictionary
        self.summary["meta"] = {
            "quick": self.quick,
            "cog_functions_are_called": self.cog_functions_are_called,
            "cog_categories_are_called": self.cog_categories_are_called,
            "output_directory": self.output_directory,
            "summary_type": self.summary_type,
            "collection": bin_ids,
            "num_bins": len(bin_ids),
            "collection_name": self.collection_name,
            "total_num_genes_in_collection": 0,
            "anvio_version": __version__,
            "pan": self.p_meta,
            "genomes": {
                "num_genomes": self.genomes_storage.num_genomes,
                "functions_available": (
                    True if len(self.gene_clusters_function_sources) else False
                ),
                "function_sources": self.gene_clusters_function_sources,
            },
            "percent_of_genes_collection": 0.0,
            "genome_names": genome_names,
        }

        # I am not sure whether this is the best place to do this,
        self.summary["basics_pretty"] = {
            "pan": [
                ("Created on", self.p_meta["creation_date"]),
                ("Version", anvio.__pan__version__),
                (
                    "Number of genes",
                    pretty(int(self.p_meta["num_genes_in_gene_clusters"])),
                ),
                (
                    "Number of gene clusters",
                    pretty(int(self.p_meta["num_gene_clusters"])),
                ),
                (
                    "Partial genes excluded",
                    "Yes" if self.p_meta["exclude_partial_gene_calls"] else "No",
                ),
                ("Minbit parameter", self.p_meta["minbit"]),
                (
                    "Gene cluster min occurrence parameter",
                    pretty(int(self.p_meta["gene_cluster_min_occurrence"])),
                ),
                ("MCL inflation parameter", self.p_meta["mcl_inflation"]),
                (
                    "NCBI blastp or DIAMOND?",
                    "NCBI blastp" if self.p_meta["use_ncbi_blast"] else ("DIAMOND"),
                ),
                (
                    "Additional parameters for sequence search",
                    self.p_meta["additional_params_for_seq_search"],
                ),
                ("Number of genomes used", pretty(int(self.p_meta["num_genomes"]))),
                (
                    "Items aditional data keys",
                    (
                        "--"
                        if not self.items_additional_data_keys
                        else ", ".join(self.items_additional_data_keys)
                    ),
                ),
            ],
            "genomes": [
                ("Created on", "Storage DB knows nothing :("),
                ("Version", anvio.__genomes_storage_version__),
                (
                    "Number of genomes described",
                    pretty(self.genomes_storage.num_genomes),
                ),
                (
                    "Functional annotation",
                    (
                        "Available"
                        if len(self.gene_clusters_function_sources)
                        else "Not available :/"
                    ),
                ),
                (
                    "Functional annotation sources",
                    (
                        "--"
                        if not len(self.gene_clusters_function_sources)
                        else ", ".join(self.gene_clusters_function_sources)
                    ),
                ),
            ],
        }

        self.summary["files"] = {}
        self.summary["collection_profile"] = (
            self.collection_profile
        )  # reminder; collection_profile comes from the superclass!

        self.generate_gene_clusters_file(self.collection_dict)

        self.report_misc_data_files(target_table="layers")
        self.report_misc_data_files(target_table="items")

        if self.debug:
            import json

            print(json.dumps(self.summary, sort_keys=True, indent=4))

        self.index_html = SummaryHTMLOutput(
            self.summary, r=self.run, p=self.progress
        ).generate(quick=self.quick)

    def generate_gene_clusters_file(self, collection_dict, compress_output=True):
        """Generates the gene summary file"""

        # generate a dict of gene cluster ~ bin id relationships
        gene_cluster_name_to_bin_name = dict(
            list(
                zip(
                    self.gene_clusters_in_pan_db_but_not_binned,
                    [None] * len(self.gene_clusters_in_pan_db_but_not_binned),
                )
            )
        )
        for bin_id in collection_dict:
            for gene_cluster_name in collection_dict[bin_id]:
                gene_cluster_name_to_bin_name[gene_cluster_name] = bin_id

        ###############################################
        # generate an output file for gene clusters.
        ###############################################
        output_file_obj = self.get_output_file_handle(
            prefix="gene_clusters_summary.txt",
            compress_output=compress_output,
            add_project_name=True,
        )

        # standard headers
        header = [
            "unique_id",
            "gene_cluster_id",
            "bin_name",
            "genome_name",
            "gene_callers_id",
        ]

        # extend the header with items additional data keys
        for items_additional_data_key in self.items_additional_data_keys:
            header.append(items_additional_data_key)

        # extend the header with functions if there are any
        for function_source in self.gene_clusters_function_sources:
            if self.quick:
                header.append(function_source + "_ACC")
            else:
                header.append(function_source + "_ACC")
                header.append(function_source)

        # if this is not a quick summary, have sequences in the output
        sequences = None
        if not self.quick:
            header.append(
                "dna_sequence" if self.report_DNA_sequences else "aa_sequence"
            )
            sequences = self.get_sequences_for_gene_clusters(
                gene_cluster_names=self.gene_cluster_names,
                report_DNA_sequences=self.report_DNA_sequences,
            )

        # write the header
        output_file_obj.write(("\t".join(header) + "\n").encode("utf-8"))

        self.progress.new("Gene clusters summary file")
        self.progress.update("...")

        # uber loop for the file content
        unique_id = 1
        for gene_cluster_name in self.gene_clusters:
            for genome_name in self.gene_clusters[gene_cluster_name]:
                for gene_caller_id in self.gene_clusters[gene_cluster_name][
                    genome_name
                ]:
                    entry = [
                        unique_id,
                        gene_cluster_name,
                        gene_cluster_name_to_bin_name[gene_cluster_name],
                        genome_name,
                        gene_caller_id,
                    ]

                    # populate the entry with item aditional data
                    for items_additional_data_key in self.items_additional_data_keys:
                        if gene_cluster_name in self.items_additional_data_dict:
                            entry.append(
                                self.items_additional_data_dict[gene_cluster_name][
                                    items_additional_data_key
                                ]
                            )
                        else:
                            entry.append("")

                    # populate the entry with functions.
                    for function_source in self.gene_clusters_function_sources:
                        annotations_dict = self.gene_clusters_functions_dict[
                            gene_cluster_name
                        ][genome_name][gene_caller_id]
                        if function_source in annotations_dict:
                            annotation_blob = self.gene_clusters_functions_dict[
                                gene_cluster_name
                            ][genome_name][gene_caller_id][function_source]

                            # FIXME: this is an artifact from Py2 to Py3 swtich. DBs generated in Py2 and used from Py3 will have type
                            # bytes for annotation_blob. so we will convert them to str.. If the db is generated with Py3, there is no
                            # such problem, so we can remove this extra step around July 2017.
                            if isinstance(annotation_blob, bytes):
                                annotation_blob = annotation_blob.decode("utf-8")

                            accessions, annotations = [
                                l.split("!!!") for l in annotation_blob.split("|||")
                            ]
                            entry.append("|".join(accessions))
                            entry.append("|".join(annotations))
                        else:
                            entry.append("")
                            entry.append("")

                    if not self.quick:
                        entry.append(
                            sequences[gene_cluster_name][genome_name][gene_caller_id]
                        )

                    output_file_obj.write(
                        (
                            "\t".join(
                                [
                                    str(e) if e not in [None, "UNKNOWN"] else ""
                                    for e in entry
                                ]
                            )
                            + "\n"
                        ).encode("utf-8")
                    )
                    unique_id += 1

        # we're done here.
        output_file_obj.close()

        self.progress.end()


class ProfileSummarizer(DatabasesMetaclass, SummarizerSuperClass):
    """Creates an über dictionary of 'summary' for anvi'o profiles."""

    def __init__(self, args=None, lazy_init=False, r=run, p=progress):
        self.args = args

        self.run = r
        self.progress = p
        self.lazy_init = lazy_init

        self.summary = {}
        self.summary_type = "profile"
        self.debug = False
        self.quick = False
        self.profile_db_path = None
        self.contigs_db_path = None
        self.output_directory = None
        self.split_names_per_bin = None
        self.completeness_data_available = False
        self.gene_level_coverage_stats_available = False
        self.non_single_copy_gene_hmm_data_available = False
        self.reformat_contig_names = False

        DatabasesMetaclass.__init__(self, self.args, self.run, self.progress)
        SummarizerSuperClass.__init__(self, self.args, self.run, self.progress)

        # databases initiated, let's make sure we have gene covereges data avaialable.
        if self.gene_level_coverage_stats_dict:
            self.gene_level_coverage_stats_available = True

        self.collection_dict = {}
        self.bins_info_dict = {}
        self.initialized = False

    def init(self):
        # init profile data for colletion.
        self.collection_dict, self.bins_info_dict = self.init_collection_profile(
            self.collection_name, calculate_Q2Q3_carefully=self.calculate_Q2Q3_carefully
        )

        # let bin names known to all
        self.bin_ids = list(self.collection_profile.keys())

        # load completeness information if available
        self.completeness = completeness.Completeness(self.contigs_db_path)
        if len(self.completeness.sources):
            self.completeness_data_available = True

        # load HMM sources for non-single-copy genes if available
        if self.non_singlecopy_gene_hmm_sources and not self.quick:
            self.init_non_singlecopy_gene_hmm_sources()
            self.non_single_copy_gene_hmm_data_available = True

        # load gene functions from contigs db superclass
        self.init_functions()

        # get an instance of SCG taxonomy:
        if self.a_meta["scg_taxonomy_was_run"]:
            self.scg_taxonomy = scgtaxonomyops.SCGTaxonomyEstimatorSingle(
                argparse.Namespace(contigs_db=self.contigs_db_path)
            )
        else:
            self.scg_taxonomy = None

        # set up the initial summary dictionary
        self.summary["meta"] = {
            "quick": self.quick,
            "summary_type": self.summary_type,
            "output_directory": self.output_directory,
            "collection": self.bin_ids,
            "num_bins": len(self.bin_ids),
            "collection_name": self.collection_name,
            "total_nts_in_collection": 0,
            "num_contigs_in_collection": 0,
            "anvio_version": __version__,
            "profile": self.p_meta,
            "contigs": self.a_meta,
            "gene_level_coverage_stats_available": self.gene_level_coverage_stats_available,
            "completeness_data_available": self.completeness_data_available,
            "non_single_copy_gene_hmm_data_available": self.non_single_copy_gene_hmm_data_available,
            "percent_contigs_nts_described_by_collection": 0.0,
            "percent_profile_nts_described_by_collection": 0.0,
            "percent_contigs_nts_described_by_profile": P(
                self.p_meta["total_length"], self.a_meta["total_length"]
            ),
            "percent_contigs_contigs_described_by_profile": P(
                self.p_meta["num_contigs"], self.a_meta["num_contigs"]
            ),
            "percent_contigs_splits_described_by_profile": P(
                self.p_meta["num_splits"], self.a_meta["num_splits"]
            ),
        }

        # FIXME: mistune had an API change, and is not behaving correctly
        # this is the part where we try to initialize
        markdown = None
        try:
            # this is for mistune `0.8.4` and earlier
            renderer = mistune.Renderer(escape=False)
            markdown = mistune.Markdown(renderer=renderer)
        except AttributeError:
            # this is for mistune `2.0.0rc1` and later
            try:
                markdown = mistune.create_markdown(escape=False)
            except Exception as e:
                self.run.warning(
                    f"Well :( Anvi'o failed to initialize `mistune`. This is the error "
                    f"we got from the downstream library: '{e}'. Probably this needs a developer "
                    f"to take a look. Meanwhile, we are turning off the rendering of markdown "
                    f"descriptions. Things will look ugly in some places in the interface and "
                    f"summary, but at least things will work from a functional standpoint :/"
                )

        # I am not sure whether this is the best place to do this,
        T = lambda x: "True" if x else "False"

        J = lambda x: x in self.p_meta and self.p_meta[x] != None

        self.summary["basics_pretty"] = {
            "profile": [
                ("Created on", self.p_meta["creation_date"]),
                ("Version", self.p_meta["version"]),
                ("Number of contigs", pretty(int(self.p_meta["num_contigs"]))),
                ("Number of splits", pretty(int(self.p_meta["num_splits"]))),
                ("Contig length cutoff min", pretty(self.p_meta["min_contig_length"])),
                ("Contig length cutoff max", pretty(self.p_meta["max_contig_length"])),
                ("Samples in profile", ", ".join(self.p_meta["samples"])),
                ("Total nucleotides", humanize_n(int(self.p_meta["total_length"]))),
                ("SNVs profiled", T(self.p_meta["SNVs_profiled"])),
                ("SCVs profiled", T(self.p_meta["SCVs_profiled"])),
                ("IN/DELs profiled", T(self.p_meta["INDELs_profiled"])),
                ("Report variability full", T(self.p_meta["report_variability_full"])),
                (
                    "Min coverage for variability",
                    (
                        humanize_n(int(self.p_meta["min_coverage_for_variability"]))
                        if J("min_coverage_for_variability")
                        else "NA"
                    ),
                ),
            ],
            "contigs": [
                ("Project name", self.a_meta["project_name"]),
                ("Created on", self.a_meta["creation_date"]),
                ("Version", self.a_meta["version"]),
                ("Total nucleotides", humanize_n(int(self.a_meta["total_length"]))),
                ("Number of contigs", pretty(int(self.a_meta["num_contigs"]))),
                ("Number of splits", pretty(int(self.a_meta["num_splits"]))),
                ("Genes are called", T(self.a_meta["genes_are_called"])),
                ("External gene calls", T(self.a_meta["external_gene_calls"])),
                (
                    "External amino acid sequences",
                    T(self.a_meta["external_gene_amino_acid_seqs"]),
                ),
                ("K-mer size", self.a_meta["kmer_size"]),
                ("Split length", pretty(int(self.a_meta["split_length"]))),
                (
                    "Splits consider gene calls",
                    T(self.a_meta["splits_consider_gene_calls"]),
                ),
                ("SCG taxonomy was run", T(self.a_meta["scg_taxonomy_was_run"])),
                (
                    "Gene function sources",
                    (
                        ", ".join(self.gene_function_call_sources)
                        if self.gene_function_call_sources
                        else "None :("
                    ),
                ),
                ("Summary reformatted contig names", self.reformat_contig_names),
            ],
            "description": (
                markdown(self.p_meta["description"])
                if markdown
                else self.p_meta["description"]
            ),
        }

        self.summary["max_shown_header_items"] = 10
        self.summary["slice_header_items_tmpl"] = (
            "0:%d" % self.summary["max_shown_header_items"]
        )
        self.summary["num_not_shown_samples"] = (
            len(self.p_meta["samples"]) - self.summary["max_shown_header_items"]
        )
        self.summary["num_not_shown_hmm_items"] = dict(
            [
                (
                    hmm_search_source,
                    len(self.hmm_sources_info[hmm_search_source]["genes"])
                    - self.summary["max_shown_header_items"],
                )
                for hmm_search_type, hmm_search_source in self.hmm_searches_header
            ]
        )

        self.summary["files"] = {}
        self.summary["misc_data"] = {}
        self.summary["collection"] = {}
        self.summary["collection_profile"] = (
            self.collection_profile
        )  # reminder; collection_profile comes from ProfileSuperclass!
        self.summary["collection_profile_items"] = (
            []
            if not len(list(self.collection_profile.values()))
            else list(self.collection_profile.values())[0].keys()
        )

        # add hmm items for each seach type:
        if self.non_single_copy_gene_hmm_data_available:
            self.summary["meta"]["hmm_items"] = dict(
                [
                    (
                        hmm_search_source,
                        self.hmm_sources_info[hmm_search_source]["genes"],
                    )
                    for hmm_search_type, hmm_search_source in self.hmm_searches_header
                ]
            )

        # yay
        self.initialized = True

    def process(self):
        if not self.initialized:
            self.init()

        if not self.output_directory:
            raise ConfigError(
                "It seems the summarizer class have been inherited without an `output_directory` argument :/ Show stopper "
                "mistake stopped the show. Bye!"
            )

        # summarize bins:
        self.progress.new("Summarizing ...", progress_total_items=len(self.bin_ids))
        for i in range(0, len(self.bin_ids)):
            bin_id = self.bin_ids[i]

            self.progress.increment()
            self.progress.update_pid(
                "Summarizing %d of %d: '%s'" % (i + 1, len(self.bin_ids), bin_id)
            )
            self.progress.update("...")

            bin = Bin(self, bin_id, self.run, self.progress)
            bin.output_directory = os.path.join(
                self.output_directory, "bin_by_bin", bin_id
            )
            bin.bin_profile = self.collection_profile[bin_id]

            self.summary["collection"][bin_id] = bin.create()
            self.summary["collection"][bin_id]["color"] = (
                self.bins_info_dict[bin_id]["html_color"] or "#212121"
            )
            self.summary["collection"][bin_id]["source"] = (
                self.bins_info_dict[bin_id]["source"] or "unknown_source"
            )
            self.summary["meta"]["total_nts_in_collection"] += self.summary[
                "collection"
            ][bin_id]["total_length"]
            self.summary["meta"]["num_contigs_in_collection"] += self.summary[
                "collection"
            ][bin_id]["num_contigs"]

        self.progress.end()

        # bins are computed, add some relevant meta info:
        self.summary["meta"]["percent_contigs_nts_described_by_collection"] = "%.2f" % (
            self.summary["meta"]["total_nts_in_collection"]
            * 100.0
            / int(self.a_meta["total_length"])
        )
        self.summary["meta"]["percent_profile_nts_described_by_collection"] = "%.2f" % (
            self.summary["meta"]["total_nts_in_collection"]
            * 100.0
            / int(self.p_meta["total_length"])
        )
        self.summary["meta"]["bins"] = self.get_bins_ordered_by_completeness_and_size()

        if not self.quick:
            # generate a TAB-delimited text output file for bin summaries
            summary_of_bins = {}
            properties = ["total_length", "num_contigs", "N50", "GC_content"]

            if self.completeness_data_available:
                properties += ["percent_completion", "percent_redundancy"]

            if self.scg_taxonomy:
                properties += constants.levels_of_taxonomy

            for bin_name in self.summary["collection"]:
                summary_of_bins[bin_name] = dict(
                    [
                        (prop, self.summary["collection"][bin_name][prop])
                        for prop in properties
                    ]
                )

            output_file_obj = self.get_output_file_handle(prefix="bins_summary.txt")
            utils.store_dict_as_TAB_delimited_file(
                summary_of_bins,
                None,
                headers=["bins"] + properties,
                file_obj=output_file_obj,
            )

            self.report_misc_data_files(target_table="layers")
            self.report_misc_data_files(target_table="items")

            # save merged matrices for bins x samples
            for table_name in self.summary["collection_profile_items"]:
                d = {}
                for bin_id in self.collection_profile:
                    d[bin_id] = self.collection_profile[bin_id][table_name]

                output_file_obj = self.get_output_file_handle(
                    sub_directory="bins_across_samples", prefix="%s.txt" % table_name
                )
                utils.store_dict_as_TAB_delimited_file(
                    d,
                    None,
                    headers=["bins"] + sorted(self.p_meta["samples"]),
                    file_obj=output_file_obj,
                )

            # merge and store matrices for hmm hits
            if self.non_single_copy_gene_hmm_data_available:
                for hmm_search_source in self.summary["meta"]["hmm_items"]:
                    # this is to keep numbers per hmm item:
                    d = {}

                    for bin_id in self.summary["meta"]["bins"]:
                        d[bin_id] = self.summary["collection"][bin_id]["hmms"][
                            hmm_search_source
                        ]

                    output_file_obj = self.get_output_file_handle(
                        sub_directory="bins_across_samples",
                        prefix="%s.txt" % hmm_search_source,
                        within="hmms",
                    )
                    utils.store_dict_as_TAB_delimited_file(
                        d,
                        None,
                        headers=["bins"]
                        + sorted(self.summary["meta"]["hmm_items"][hmm_search_source]),
                        file_obj=output_file_obj,
                    )

                # this is to keep number of hmm hits per bin:
                n = dict([(bin_id, {}) for bin_id in self.summary["meta"]["bins"]])
                for hmm_search_source in self.summary["meta"]["hmm_items"]:
                    for bin_id in self.summary["meta"]["bins"]:
                        n[bin_id][hmm_search_source] = sum(
                            self.summary["collection"][bin_id]["hmms"][
                                hmm_search_source
                            ].values()
                        )

                output_file_obj = self.get_output_file_handle(
                    sub_directory="bins_across_samples", prefix="hmm_hit_totals.txt"
                )
                utils.store_dict_as_TAB_delimited_file(
                    n,
                    None,
                    headers=["bins"] + sorted(self.summary["meta"]["hmm_items"]),
                    file_obj=output_file_obj,
                )

            # store percent abundance of each bin
            if self.p_meta["blank"]:
                self.summary["bin_percent_recruitment"] = None
                self.summary["bin_percent_abundance_items"] = None
            else:
                self.summary["bin_percent_recruitment"] = (
                    self.bin_percent_recruitment_per_sample
                )
                self.summary["bin_percent_abundance_items"] = sorted(
                    list(self.bin_percent_recruitment_per_sample.values())[0].keys()
                )
                output_file_obj = self.get_output_file_handle(
                    sub_directory="bins_across_samples",
                    prefix="bins_percent_recruitment.txt",
                )
                utils.store_dict_as_TAB_delimited_file(
                    self.bin_percent_recruitment_per_sample,
                    None,
                    headers=["samples"]
                    + sorted(self.collection_profile.keys())
                    + ["__splits_not_binned__"],
                    file_obj=output_file_obj,
                )

        if self.reformat_contig_names:
            self.run.warning(
                "You have asked anvi'o to reformat contig names for bins in the summary output. Which means, the original names "
                "found in the contigs database and BAM files are no longer there in FASTA files (and hopefully all other relevant "
                "files) for your bins. Instead, they are replaced to include the bin name, and they look very neat. Just to make "
                "sure you have an idea how the name conversion looked like, anvi'o kept a copy of the conversion map for each bin "
                "you can find under directoris stored under `bin_by_bin/` directory. Please be extra careful for your downstream "
                "analyses to make sure this change will not break things."
            )

        if self.debug:
            import json

            print(json.dumps(self.summary, sort_keys=True, indent=4))

        self.index_html = SummaryHTMLOutput(
            self.summary, r=self.run, p=self.progress
        ).generate(quick=self.quick)

    def get_bins_ordered_by_completeness_and_size(self):
        if self.completeness_data_available:
            return [
                t[2]
                for t in sorted(
                    [
                        (
                            self.summary["collection"][bin]["percent_completion"],
                            self.summary["collection"][bin]["total_length"],
                            bin,
                        )
                        for bin in self.summary["collection"]
                    ],
                    reverse=True,
                )
            ]
        else:
            return sorted(self.summary["collection"].keys())


class ContigSummarizer(SummarizerSuperClass):
    def __init__(self, contigs_db_path, run=run, progress=progress):
        self.contigs_db_path = contigs_db_path
        self.run = run
        self.progress = progress

    def get_contigs_db_info_dict(
        self,
        run=run,
        progress=progress,
        include_AA_counts=False,
        split_names=None,
        gene_caller_to_use=None,
    ):
        """Returns an info dict for a given contigs db.

        Please note that this function will only return gene calls made by `gene_caller_to_use`,
        but it will report other gene callers found in the contigs database, and how many genes
        were not reported. The client side should check for those to report to the user.
        """

        if not gene_caller_to_use:
            gene_caller_to_use = constants.default_gene_caller

        args = argparse.Namespace(contigs_db=self.contigs_db_path)

        if split_names:
            args.split_names_of_interest = set(split_names)

        run = terminal.Run()
        progress = terminal.Progress()
        run.verbose = False
        progress.verbose = False

        c = ContigsSuperclass(args, r=run, p=progress)

        info_dict = {
            "path": self.contigs_db_path,
            "gene_caller_ids": set([]),
            "gene_caller": gene_caller_to_use,
        }

        for key in c.a_meta:
            info_dict[key] = c.a_meta[key]

        gene_calls_from_other_gene_callers = Counter()
        impossible_gene_calls_missing_from_contigs_db = set([])

        # Two different strategies here depending on whether we work with a given set if split ids or
        # everything in the contigs database.
        def process_gene_call(g):
            if g not in c.genes_in_contigs_dict:
                impossible_gene_calls_missing_from_contigs_db.add(g)
                return

            gene_caller = c.genes_in_contigs_dict[g]["source"]
            if gene_caller == gene_caller_to_use:
                info_dict["gene_caller_ids"].add(g)
            else:
                gene_calls_from_other_gene_callers[gene_caller] += 1

        if split_names:
            split_names = set(split_names)
            c.init_split_sequences()
            seq = "".join(
                [
                    c.split_sequences[split_name]["sequence"]
                    for split_name in split_names
                ]
            )
            for e in list(c.genes_in_splits.values()):
                if e["split"] in split_names:
                    process_gene_call(e["gene_callers_id"])

            info_dict["num_splits"] = len(split_names)
            info_dict["num_contigs"] = len(
                set(
                    [
                        c.splits_basic_info[split_name]["parent"]
                        for split_name in split_names
                    ]
                )
            )
        else:
            c.init_contig_sequences()
            seq = "".join([e["sequence"] for e in list(c.contig_sequences.values())])

            for g in c.genes_in_contigs_dict:
                process_gene_call(g)

        if len(gene_calls_from_other_gene_callers):
            run.info_single(
                "PLEASE READ CAREFULLY. Contigs db info summary will not include %d gene calls that were "
                'not identified by "%s", the default gene caller. Other gene calls found in this contigs '
                "database include, %s. If you are more interested in gene calls in any of those, you should "
                "indicate that through the `--gene-caller` parameter in your program."
                % (
                    sum(gene_calls_from_other_gene_callers.values()),
                    gene_caller_to_use,
                    ", ".join(
                        [
                            "%d gene calls by %s" % (tpl[1], tpl[0])
                            for tpl in gene_calls_from_other_gene_callers.items()
                        ]
                    ),
                )
            )

        if len(impossible_gene_calls_missing_from_contigs_db):
            self.progress.reset()
            self.run.warning(
                f"SOMETHING IMPOSSIBLE HAS HAPPENED: splits in your contigs database contains gene calls that "
                f"are not found in the gene calls table. This really should not happen under any circumstance "
                f"yet there were {len(impossible_gene_calls_missing_from_contigs_db)} of them in your database."
                f"Anvi'o ignored those gene calls and is going to continue processing your data, but what you "
                f"should be doing is to remove this contigs database, and start everything from scratch. There is "
                f"something quite wrong with it. If you are curious, here are the gene calls in question: "
                f"{', '.join([str(g) for g in impossible_gene_calls_missing_from_contigs_db])}",
                overwrite_verbose=True,
            )

        info_dict["gene_calls_from_other_gene_callers"] = (
            gene_calls_from_other_gene_callers
        )
        info_dict["gc_content"] = seqlib.Composition(seq).GC_content
        info_dict["total_length"] = len(seq)

        info_dict["partial_gene_calls"] = set([])
        for gene_caller_id in info_dict["gene_caller_ids"]:
            if c.genes_in_contigs_dict[gene_caller_id]["partial"]:
                info_dict["partial_gene_calls"].add(gene_caller_id)

        info_dict["num_genes"] = len(info_dict["gene_caller_ids"])
        if info_dict["num_genes"]:
            info_dict["avg_gene_length"] = numpy.mean(
                [
                    c.gene_lengths[gene_caller_id]
                    for gene_caller_id in info_dict["gene_caller_ids"]
                ]
            )
            info_dict["num_genes_per_kb"] = (
                info_dict["num_genes"] * 1000.0 / info_dict["total_length"]
            )
        else:
            info_dict["avg_gene_length"], info_dict["num_genes_per_kb"] = 0.0, 0

        # get completeness / contamination estimates
        (
            p_completion,
            p_redundancy,
            domain,
            domain_probabilities,
            info_text,
            results_dict,
        ) = completeness.Completeness(self.contigs_db_path).get_info_for_splits(
            split_names if split_names else set(c.splits_basic_info.keys())
        )
        domain_confidence = domain_probabilities[domain] if domain else 0.0

        info_dict["hmm_sources_info"] = c.hmm_sources_info
        info_dict["percent_completion"] = p_completion
        info_dict["percent_redundancy"] = p_redundancy
        info_dict["scg_domain"] = domain
        info_dict["scg_domain_confidence"] = domain_confidence

        info_dict["hmms_for_scgs_were_run"] = True if len(results_dict) else False

        # lets get all amino acids used in all complete gene calls:
        if include_AA_counts:
            if split_names:
                AA_counts_dict = c.get_AA_counts_dict(split_names=split_names)
            else:
                AA_counts_dict = c.get_AA_counts_dict()

            info_dict["AA_counts"] = AA_counts_dict["AA_counts"]
            info_dict["total_AAs"] = AA_counts_dict["total_AAs"]

        missing_keys = [
            key for key in constants.essential_genome_info if key not in info_dict
        ]
        if len(missing_keys):
            raise ConfigError(
                "We have a big problem. I am reporting from get_contigs_db_info_dict. This function must "
                "produce a dictionary that meets the requirements defined in the constants module of anvi'o "
                "for 'essential genome info'. But when I look at the resulting dictionary this function is "
                "about to return, I can see it is missing some stuff :/ This is not a user error, but it needs "
                "the attention of an anvi'o developer. Here are the keys that should have been in the results "
                "but missing: '%s'" % (", ".join(missing_keys))
            )

        return info_dict

    def get_summary_dict_for_assembly(self, gene_caller_to_use=None):
        """Returns a simple summary dict for a given contigs database"""
        if not gene_caller_to_use:
            gene_caller_to_use = constants.default_gene_caller

        self.progress.new("Generating contigs db summary")

        self.progress.update(
            "Initiating contigs super for %s ..." % self.contigs_db_path
        )
        run, progress = terminal.Run(), terminal.Progress()
        run.verbose, progress.verbose = False, False
        c = ContigsSuperclass(
            argparse.Namespace(contigs_db=self.contigs_db_path), r=run, p=progress
        )

        self.progress.update("Recovering info about %s ..." % self.contigs_db_path)
        num_genes = len(
            [
                True
                for v in c.genes_in_contigs_dict.values()
                if v["source"] == gene_caller_to_use
            ]
        )
        project_name = c.a_meta["project_name"]
        contig_lengths = sorted(
            [e["length"] for e in c.contigs_basic_info.values()], reverse=True
        )
        total_length = sum(contig_lengths)
        num_contigs = len(contig_lengths)

        self.progress.update("Figuring out HMM hits in %s ..." % self.contigs_db_path)
        hmm = hmmops.SequencesForHMMHits(self.contigs_db_path, run=self.run)

        self.progress.update("Summarizing %s ..." % self.contigs_db_path)
        summary = {}
        summary["project_name"] = project_name
        summary["total_length"] = total_length
        summary["num_genes"] = num_genes
        summary["gene_caller"] = gene_caller_to_use
        summary["num_contigs"] = num_contigs
        summary["n_values"] = self.calculate_N_values(
            contig_lengths, total_length, N=100
        )
        summary["contig_lengths"] = contig_lengths
        summary["gene_hit_counts_per_hmm_source"] = (
            hmm.get_gene_hit_counts_per_hmm_source()
        )
        summary["hmm_sources_for_SCGs"] = sorted(
            [
                s
                for s in hmm.hmm_hits_info
                if hmm.hmm_hits_info[s]["search_type"] == "singlecopy"
            ]
        )
        summary["num_genomes_per_SCG_source_dict"] = (
            hmm.get_num_genomes_from_SCG_sources_dict()
        )

        self.progress.end()

        return summary

    def calculate_N_values(self, contig_lengths, total_length, N=100):
        results = []

        temp_length = 0
        contigs_index = 0
        n_index = 1

        while n_index <= N:
            if temp_length >= int(((total_length / N) * n_index)):
                results.append(
                    {
                        "num_contigs": contigs_index,
                        "length": contig_lengths[contigs_index - 1],
                    }
                )
                n_index += 1
            else:
                temp_length += contig_lengths[contigs_index]
                contigs_index += 1

        return results


class PanBin:
    def __init__(self, summary, bin_id, r=run, p=progress):
        self.summary = summary
        self.progress = p
        self.run = r

        if bin_id not in self.summary.bin_ids:
            raise ConfigError(
                "Bin '%s' does not seem to be in this summary :/ These are the ones in it: %s."
                % (bin_id, ", ".join(self.summary.bin_ids))
            )

        self.bin_id = bin_id
        self.split_names = summary.collection_dict[self.bin_id]
        self.views = {}
        self.genomes = []

        self.num_gene_clusters = None
        self.num_genes_in_gene_clusters = None

        self.num_splits = len(self.split_names)

        # here we are subsetting views based on what gene clusters are available
        # in this particular bin.
        for view in self.summary.views:
            self.views[view] = {}
            self.views[view]["table_name"] = self.summary.views[view]["table_name"]
            self.views[view]["header"] = self.summary.views[view]["header"]
            self.views[view]["dict"] = {}

            for split_name in self.split_names:
                self.views[view]["dict"][split_name] = self.summary.views[view]["dict"][
                    split_name
                ]

            if not self.genomes:
                self.genomes = self.views[view]["header"]

        if "gene_cluster_frequencies" in self.views:
            self.num_gene_clusters = len(self.views["gene_cluster_frequencies"]["dict"])
            self.num_genes_in_gene_clusters = sum(
                [
                    sum(x.values())
                    for x in self.views["gene_cluster_frequencies"]["dict"].values()
                ]
            )


class Bin:
    def __init__(self, summary, bin_id, r=run, p=progress):
        self.summary = summary

        if not self.summary.initialized:
            raise ConfigError(
                "The summary object you sent to the `Bin` class to make sense of the bin '%s' does "
                "not seem to have been initialized. Anvi'o could have taken care of it for you, but "
                "it will not (not only because anvi'o is implemented by mean people, but also it kinda "
                "likes to be explicit about this kind of stuff). Please initialize your summary object "
                "first." % (bin_id)
            )

        if bin_id not in self.summary.bin_ids:
            raise ConfigError(
                "Bin '%s' does not seem to be in this summary :/ These are the ones in it: %s."
                % (bin_id, ", ".join(self.summary.bin_ids))
            )

        self.bin_id = bin_id
        self.split_names = summary.collection_dict[self.bin_id]
        self.progress = p
        self.run = r

        self.across_samples = {}
        self.bin_profile = {}
        self.bin_info_dict = {"files": {}}

        # this dictionary will keep new contig names if the user asked contig names in the summary output
        # to be reformatted nicely
        self.contig_name_conversion_dict = {}

        # this will quickly populate self.contig_names, self.total_length, and self.num_contigs variables
        self.process_contigs(quick=True)

        self.output_directory = None
        self.contig_lengths = []

        # make sure all split_names in the collection is actually in the contigs database.
        # in collections stored in the contigs database, split_names that are not in the
        # oritinal contigs used to generate contigs database *may* end up in the
        # collections table. we gotta make sure we deal with them properly:
        missing_ids = [
            split_id
            for split_id in self.split_names
            if split_id not in self.summary.split_sequences
        ]
        if len(missing_ids):
            for missing_id in missing_ids:
                self.split_names.remove(missing_id)

            self.run.warning(
                '%d split id(s) in bin "%s" reported by collection "%s" is not found in the '
                "contigs database and removed from the bin summary. If this does not make "
                "any sense, you may need make sure everything is in order. The thing is, "
                "sometimes external clustering results that are added to the contigs via "
                "`anvi-populate-collections-table` may include split names that are not used "
                "while the contigs database was generated."
                % (len(missing_ids), bin_id, self.summary.collection_name)
            )

        self.gene_caller_ids = self.get_gene_caller_ids()
        self.num_splits = len(self.split_names)

        # make these dicts avilable:
        self.gene_level_coverage_stats_dict = {}
        self.split_coverage_values_per_nt_dict = {}

        A = lambda x: self.summary.gene_level_coverage_stats_dict[gene_callers_id][
            sample_name
        ][x]

        # populate gene coverage and detection dictionaries by subsetting them from the parent summary object
        if self.summary.gene_level_coverage_stats_dict:
            for gene_callers_id in self.gene_caller_ids:
                self.gene_level_coverage_stats_dict[gene_callers_id] = {}

                for sample_name in self.summary.p_meta["samples"]:
                    self.gene_level_coverage_stats_dict[gene_callers_id][
                        sample_name
                    ] = {
                        "mean_coverage": A("mean_coverage"),
                        "detection": A("detection"),
                        "non_outlier_mean_coverage": A("non_outlier_mean_coverage"),
                        "non_outlier_coverage_std": A("non_outlier_coverage_std"),
                    }

                    if (
                        "gene_coverage_values_per_nt"
                        in self.summary.gene_level_coverage_stats_dict[gene_callers_id][
                            sample_name
                        ]
                    ):
                        self.gene_level_coverage_stats_dict[gene_callers_id][
                            sample_name
                        ]["gene_coverage_values_per_nt"] = A(
                            "gene_coverage_values_per_nt"
                        )
                        self.gene_level_coverage_stats_dict[gene_callers_id][
                            sample_name
                        ]["non_outlier_positions"] = A("non_outlier_positions")

        # populate coverage values per nucleutide for the bin.
        if self.summary.split_coverage_values_per_nt_dict:
            for split_name in self.split_names:
                self.split_coverage_values_per_nt_dict[split_name] = (
                    self.summary.split_coverage_values_per_nt_dict[split_name]
                )

    def create(self):
        self.create_bin_dir()

        self.process_contigs()

        self.store_sequences_for_hmm_hits()

        if self.summary.completeness_data_available:
            self.access_completeness_scores()

        if self.summary.non_single_copy_gene_hmm_data_available:
            self.summarize_hmm_hits()

        self.compute_basic_stats()

        self.recover_scg_taxonomy()

        self.store_genes_basic_info()

        self.store_gene_level_coverage_stats()

        self.store_profile_data()

        self.report_contig_name_conversion_dict()

        return self.bin_info_dict

    def create_bin_dir(self):
        self.progress.update("Creating the output directory ...")

        if not self.output_directory:
            self.progress.end()
            raise ConfigError(
                'You called Bin.create() before setting an output directory. Anvio says "nope, thanks".'
            )

        filesnpaths.gen_output_directory(self.output_directory)

    def access_completeness_scores(self):
        self.progress.update("Accessing completeness scores ...")

        (
            p_completion,
            p_redundancy,
            domain,
            domain_probabilities,
            info_text,
            results_dict,
        ) = self.summary.completeness.get_info_for_splits(set(self.split_names))
        domain_confidence = domain_probabilities[domain] if domain else 0.0

        self.bin_info_dict["completeness"] = results_dict

        self.bin_info_dict["percent_redundancy"] = p_redundancy
        self.bin_info_dict["percent_completion"] = p_completion
        self.bin_info_dict["scg_domain"] = domain
        self.bin_info_dict["scg_domain_confidence"] = domain_confidence

        for k in ["percent_redundancy", "percent_completion"]:
            self.store_data_in_file("%s.txt" % k, "%.4f" % self.bin_info_dict[k])

        self.store_data_in_file(
            "scg_domain.txt", "%s" % self.bin_info_dict["scg_domain"]
        )
        self.store_data_in_file(
            "scg_domain_confidence.txt",
            "%.2f" % self.bin_info_dict["scg_domain_confidence"],
        )

    def store_profile_data(self):
        if self.summary.quick:
            return

        self.progress.update("Storing profile data ...")

        for table_name in self.bin_profile:
            output_file_obj = self.get_output_file_handle("%s.txt" % table_name)
            utils.store_dict_as_TAB_delimited_file(
                {table_name: self.bin_profile[table_name]},
                None,
                headers=["bin"] + self.summary.p_meta["samples"],
                file_obj=output_file_obj,
            )

    def report_contig_name_conversion_dict(self):
        if not self.summary.reformat_contig_names:
            return

        self.progress.update("Storing contig name conversion dict ...")

        output_file_obj = self.get_output_file_handle("contig-name-conversion-map.txt")
        utils.store_dict_as_TAB_delimited_file(
            self.contig_name_conversion_dict,
            None,
            headers=["original_contig_name", "reformatted_contig_name"],
            file_obj=output_file_obj,
        )

    def summarize_hmm_hits(self):
        """Make sense of everything there is to make sense of regarding hmm hits.

        Unfortunately this is *VERY* complicated. Here we try to make sense of any
        HMM collection with respect to nubmer of hits that happens to be in splits
        associated with this bin, and split - hit associations. This function fills
        all the information into self.bin_mm_profile_dict, and the process function
        up above later makes sense of all to generate files and matrices, as well as
        dictionaries to diplay part of this information in the interface.
        """

        if self.summary.quick:
            return

        info_dict = {}

        # lets limit our interest space into splits that are in our bin and have hmm hits from the get go:
        split_names_with_hmm_hits = [
            split_id
            for split_id in self.split_names
            if split_id in self.summary.hmm_searches_dict
        ]

        for hmm_search_type, hmm_search_source in self.summary.hmm_searches_header:
            hmm_items = self.summary.hmm_sources_info[hmm_search_source]["genes"]
            info_dict[hmm_search_source] = dict(
                [(hmm_item, 0) for hmm_item in hmm_items]
            )

            hits_in_splits = []
            # keep track of unique identifiers of hmm hits to not count a single hit that spans across multiple splits:
            unique_identifiers_seen = set([])

            for split_id in split_names_with_hmm_hits:
                for hmm_item, unique_identifier in self.summary.hmm_searches_dict[
                    split_id
                ][hmm_search_source]:
                    hits_in_splits.append(
                        (split_id, hmm_item, unique_identifier),
                    )

                    if unique_identifier in unique_identifiers_seen:
                        continue

                    unique_identifiers_seen.add(unique_identifier)
                    info_dict[hmm_search_source][hmm_item] += 1

            output_file_obj = self.get_output_file_handle(
                "%s-hmm-hits.txt" % hmm_search_source
            )
            output_file_obj.write("contigs\thmm_profile\tunique_identifier\n")
            for item in hits_in_splits:
                output_file_obj.write("%s\n" % "\t".join(item))
            output_file_obj.close()

        self.bin_info_dict["hmms"] = info_dict

    def get_gene_caller_ids(self):
        """Returns a set of gene caller ids in the bin.

        Because splits can be cut from arbitrary locations, we may have partial hits of genes in a bin.
        we don't want genes to appear in a bin more than once due to this, or end up appearing in
        two different bins just because one bin has a fraction of a gene. here we will build the
        genes_dict, which will contain every gene hit in all splits that are found in this genome
        bin.
        """

        genes_dict = {}

        for split_name in self.split_names:
            if split_name not in self.summary.split_name_to_genes_in_splits_entry_ids:
                continue

            for gene_entry_id in self.summary.split_name_to_genes_in_splits_entry_ids[
                split_name
            ]:
                gene_call_in_split = self.summary.genes_in_splits[gene_entry_id]
                gene_callers_id = gene_call_in_split["gene_callers_id"]

                if gene_callers_id in genes_dict:
                    genes_dict[gene_callers_id].append(gene_call_in_split)
                else:
                    genes_dict[gene_callers_id] = [gene_call_in_split]

        # here we have every gene hit in this bin stored in genes_dict. what we will do is to find gene
        # call ids for genes more than 90% of which apper to be in this bin (so nothing wil be reported for
        # a gene where only like 20% of it ended up in this bin).
        gene_callers_ids_for_complete_genes = set([])
        for gene_caller_id in genes_dict:
            if sum([x["percentage_in_split"] for x in genes_dict[gene_caller_id]]) > 90:
                gene_callers_ids_for_complete_genes.add(gene_caller_id)

        del genes_dict

        return gene_callers_ids_for_complete_genes

    def store_gene_level_coverage_stats(self):
        if self.summary.quick:
            return

        if not self.summary.gene_level_coverage_stats_dict:
            return

        self.progress.update("Storing gene coverages ...")

        headers = ["gene_callers_id"] + self.summary.p_meta["samples"]

        for key, file_name in [
            ("mean_coverage", "gene_coverages.txt"),
            ("detection", "gene_detection.txt"),
            ("non_outlier_mean_coverage", "gene_non_outlier_coverages.txt"),
            ("non_outlier_coverage_std", "gene_non_outlier_coverage_stds.txt"),
        ]:
            # we will create a new dictionary here by subestting values of `key` from self.gene_level_coverage_stats_dict,
            # so we can store that information into `file_name`. magical stuff .. by us .. level 3000 wizards who can summon
            # inefficiency at most random places. SHUT UP.

            d = utils.get_values_of_gene_level_coverage_stats_as_dict(
                self.gene_level_coverage_stats_dict, key
            )

            utils.store_dict_as_TAB_delimited_file(
                d,
                None,
                headers=headers,
                file_obj=self.get_output_file_handle(file_name),
            )

    def store_genes_basic_info(self):
        if self.summary.quick:
            return

        self.progress.update("Sorting out gene calls ...")

        d = {}

        headers = ["contig", "start", "stop", "direction"]
        header_items_for_gene_sequences = ["dna_sequence"]
        if self.summary.report_aa_seqs_for_gene_calls:
            header_items_for_gene_sequences.append("aa_sequence")

        for gene_callers_id in self.gene_caller_ids:
            d[gene_callers_id] = {}
            # add sample independent information into `d`;
            for header in headers:
                if gene_callers_id not in self.summary.genes_in_contigs_dict:
                    progress.reset()
                    raise ConfigError(
                        "Bad news :( A very very rare error has occurred. A gene caller id found in your splits table "
                        "is not occurring in the genes in contigs table of your contigs database. This is due to a "
                        "rare migration error we discovered thanks to the following issue Jon Sanders filed: "
                        "https://github.com/merenlab/anvio/issues/1596. If you are seeing this error please "
                        "get in touch with us so we can help you recover from this."
                    )

                d[gene_callers_id][header] = self.summary.genes_in_contigs_dict[
                    gene_callers_id
                ][header]

            # add functions if there are any:
            if len(self.summary.gene_function_call_sources):
                for source in self.summary.gene_function_call_sources:
                    if gene_callers_id not in self.summary.gene_function_calls_dict:
                        # this gene did not get any functional annotation
                        d[gene_callers_id][source] = ""
                        d[gene_callers_id][source + " (ACCESSION)"] = ""
                        continue

                    if self.summary.gene_function_calls_dict[gene_callers_id][source]:
                        d[gene_callers_id][source + " (ACCESSION)"] = (
                            self.summary.gene_function_calls_dict[gene_callers_id][
                                source
                            ][0]
                        )
                        d[gene_callers_id][source] = (
                            self.summary.gene_function_calls_dict[gene_callers_id][
                                source
                            ][1]
                        )
                    else:
                        d[gene_callers_id][source + " (ACCESSION)"] = ""
                        d[gene_callers_id][source] = ""

            # finally add the dna and amino acid sequence for gene calls:
            contig = self.summary.genes_in_contigs_dict[gene_callers_id]["contig"]
            start = self.summary.genes_in_contigs_dict[gene_callers_id]["start"]
            stop = self.summary.genes_in_contigs_dict[gene_callers_id]["stop"]

            dna_sequence = self.summary.contig_sequences[contig]["sequence"][start:stop]
            if self.summary.genes_in_contigs_dict[gene_callers_id]["direction"] == "r":
                dna_sequence = utils.rev_comp(dna_sequence)

            d[gene_callers_id]["dna_sequence"] = dna_sequence

            # if the user asked for it, report amino acid sequences as well
            if self.summary.report_aa_seqs_for_gene_calls:
                try:
                    d[gene_callers_id]["aa_sequence"] = (
                        utils.get_translated_sequence_for_gene_call(
                            dna_sequence, gene_callers_id
                        )
                    )
                except:
                    d[gene_callers_id]["aa_sequence"] = ""

        output_file_obj = self.get_output_file_handle("gene_calls.txt")

        if self.summary.gene_function_call_sources:
            sources = [
                [source, source + " (ACCESSION)"]
                for source in self.summary.gene_function_call_sources
            ]
            headers = (
                ["gene_callers_id"]
                + headers
                + [item for sublist in sources for item in sublist]
                + header_items_for_gene_sequences
            )
        else:
            headers = ["gene_callers_id"] + headers + header_items_for_gene_sequences

        if self.summary.reformat_contig_names:
            for gene_callers_id in d:
                reformatted_contig_name = self.contig_name_conversion_dict[
                    d[gene_callers_id]["contig"]
                ]["reformatted_contig_name"]
                d[gene_callers_id]["contig"] = reformatted_contig_name

        self.progress.update("Storing genes basic info ...")
        utils.store_dict_as_TAB_delimited_file(
            d, None, headers=headers, file_obj=output_file_obj
        )

        self.bin_info_dict["genes"] = {"num_genes_found": len(self.gene_caller_ids)}

    def store_sequences_for_hmm_hits(self):
        if self.summary.quick:
            return

        s = SequencesForHMMHits(
            self.summary.contigs_db_path,
            split_names_of_interest=self.split_names,
            progress=progress_quiet,
            bin_name=self.bin_id,
        )
        hmm_sequences_dict = s.get_sequences_dict_for_hmm_hits_in_splits(
            {self.bin_id: self.split_names}
        )

        if self.summary.reformat_contig_names:
            for entry_id in hmm_sequences_dict:
                reformatted_contig_name = self.contig_name_conversion_dict[
                    hmm_sequences_dict[entry_id]["contig"]
                ]["reformatted_contig_name"]
                hmm_sequences_dict[entry_id]["contig"] = reformatted_contig_name

        single_copy_gene_hmm_sources = [
            hmm_search_source
            for hmm_search_type, hmm_search_source in self.summary.hmm_searches_header
        ]
        non_single_copy_gene_hmm_sources = self.summary.completeness.sources

        for hmm_search_source in (
            single_copy_gene_hmm_sources + non_single_copy_gene_hmm_sources
        ):
            filtered_hmm_sequences_dict = utils.get_filtered_dict(
                hmm_sequences_dict, "source", set([hmm_search_source])
            )

            output_file_obj = self.get_output_file_handle(
                "%s-hmm-sequences.txt" % hmm_search_source, key=hmm_search_source
            )

            for gene_unique_id in filtered_hmm_sequences_dict:
                header, sequence = s.get_FASTA_header_and_sequence_for_gene_unique_id(
                    hmm_sequences_dict, gene_unique_id
                )
                output_file_obj.write(">%s\n%s\n" % (header, sequence))

    def process_contigs(self, quick=False):
        """Storing contig sequences.

        This is not an easy problem. We split contigs into smaller sequences at the beginning. Only
        a part of a given contig may be used during the binning process. On the other hand we can't
        simply store sequences of splits, whenever possible, we must store the entire sequence of
        the contig (only if all splits are selected from a contig in to the same bin). So, this
        function first identifies all splits coming from the same parent, then identifies sequential
        blocks of splits (see `SequentialBlocks` class), then checks whether all splits of a given
        contig is included in the bin. If that is the case, it puts the contig as a single entry,
        witht he identical FASTA id to the original contigs in the assembly file. Otherwise it appends
        `_partial_X_Y` to the FASTA id, X and Y being the start and stop positions.
        """

        self.contig_names = set(
            [
                self.summary.splits_basic_info[split_name]["parent"]
                for split_name in self.summary.splits_basic_info
                if split_name in self.split_names
            ]
        )
        self.total_length = sum(
            [
                self.summary.splits_basic_info[split_name]["length"]
                for split_name in self.summary.splits_basic_info
                if split_name in self.split_names
            ]
        )
        self.num_contigs = len(self.contig_names)

        self.bin_info_dict["total_length"] = self.total_length
        self.bin_info_dict["contig_names"] = self.contig_names
        self.bin_info_dict["num_contigs"] = len(self.contig_names)

        if self.summary.quick or quick:
            return

        self.progress.update("Creating the FASTA file ...")

        # store original split names:
        self.store_data_in_file("original_split_names.txt", "\n".join(self.split_names))

        fasta_file = self.get_output_file_handle("contigs.fa")
        fasta_file.write(self.get_bin_sequence())
        fasta_file.close()

        self.store_data_in_file(
            "num_contigs.txt", "%d" % self.bin_info_dict["num_contigs"]
        )
        self.store_data_in_file(
            "total_length.txt", "%d" % self.bin_info_dict["total_length"]
        )

    def get_bin_sequence(self):
        output = ""

        # this dict will keep all the contig ids found in this bin with split names ordered:
        contigs_represented = utils.get_contigs_splits_dict(
            self.split_names, self.summary.splits_basic_info
        )

        # now it is time to go through each contig found in contigs_represented to
        # figure out what fraction of the contig is in fact in this bin
        contig_name_counter = 1
        for contig_name in contigs_represented:
            splits_order = list(contigs_represented[contig_name].keys())

            # this is critical: sequential_blocks is a list of one ore more lists, where each item of this list
            # describes a range of splits that follow each other to represent a coherent
            # chunk of the parent sequence (if all splits from a contig is selected into this bin,
            # then there would be one list item that spans across the entire contig):
            sequential_blocks = ccollections.GetSequentialBlocksOfSplits(
                splits_order
            ).process()

            for sequential_block in sequential_blocks:
                first_split = contigs_represented[contig_name][sequential_block[0]]
                last_split = contigs_represented[contig_name][sequential_block[-1]]

                contig_sequence_start_in_splits = self.summary.splits_basic_info[
                    first_split
                ]["start"]
                contig_sequence_end_in_splits = self.summary.splits_basic_info[
                    last_split
                ]["end"]

                # so this much of the contig is represented by its splits:
                total_contig_length_in_splits = (
                    contig_sequence_end_in_splits - contig_sequence_start_in_splits
                )

                # and this is is actual length:
                contig_sequence_length = self.summary.contigs_basic_info[contig_name][
                    "length"
                ]

                if contig_sequence_length == total_contig_length_in_splits:
                    # the entireity of the contig is represented!
                    appendix = ""
                else:
                    appendix = "_partial_%d_%d" % (
                        contig_sequence_start_in_splits,
                        contig_sequence_end_in_splits,
                    )

                sequence = ""
                for split_order in sequential_block:
                    sequence += self.summary.split_sequences[
                        contigs_represented[contig_name][split_order]
                    ]["sequence"]

                if self.summary.reformat_contig_names:
                    reformatted_contig_name = "%s_contig_%06d" % (
                        self.bin_id,
                        contig_name_counter,
                    )
                    self.contig_name_conversion_dict[contig_name] = {
                        "reformatted_contig_name": reformatted_contig_name
                    }
                    contig_name = reformatted_contig_name

                fasta_id = contig_name + appendix
                self.contig_lengths.append(len(sequence))

                output += ">%s\n" % fasta_id
                output += "%s\n" % textwrap.fill(sequence, 80, break_on_hyphens=False)

            contig_name_counter += 1

        return output

    def recover_scg_taxonomy(self):
        self.bin_info_dict["scg_taxonomy"] = None
        self.bin_info_dict["scg_taxonomy_simple"] = None

        if self.summary.quick or not self.summary.a_meta["scg_taxonomy_was_run"]:
            return

        self.progress.update("Filling in taxonomy info ...")

        scg_taxonomy_dict = self.summary.scg_taxonomy.estimate_for_list_of_splits(
            self.split_names, bin_name=self.bin_id
        )
        self.bin_info_dict["scg_taxonomy_dict"] = scg_taxonomy_dict

        self.bin_info_dict["scg_taxonomy_simple"] = "N/A"
        for level in constants.levels_of_taxonomy[::-1]:
            if scg_taxonomy_dict["consensus_taxonomy"][level]:
                self.bin_info_dict["scg_taxonomy_simple"] = scg_taxonomy_dict[
                    "consensus_taxonomy"
                ][level]
                break

        for level in constants.levels_of_taxonomy:
            self.bin_info_dict[level] = scg_taxonomy_dict["consensus_taxonomy"][level]

        scg_taxonomy_output_headers = [
            "gene_callers_id",
            "gene_name",
            "percent_identity",
            "supporting_consensus",
        ] + constants.levels_of_taxonomy
        scg_taxonomy_output = self.get_output_file_handle(
            "scg_taxonomy_details.txt", just_the_path=True
        )
        utils.store_dict_as_TAB_delimited_file(
            scg_taxonomy_dict["scgs"],
            scg_taxonomy_output,
            headers=scg_taxonomy_output_headers,
        )

    def compute_basic_stats(self):
        self.progress.update("Computing basic stats ...")

        self.bin_info_dict["N50"] = utils.get_N50(self.contig_lengths)
        self.bin_info_dict["GC_content"] = (
            numpy.mean(
                [
                    self.summary.splits_basic_info[split_id]["gc_content"]
                    for split_id in self.split_names
                ]
            )
            * 100
        )

        self.store_data_in_file(
            "N50.txt",
            "%d" % self.bin_info_dict["N50"] if self.bin_info_dict["N50"] else "NA",
        )
        self.store_data_in_file(
            "GC_content.txt", "%.4f" % self.bin_info_dict["GC_content"]
        )

    def get_output_file_handle(
        self, prefix="output.txt", overwrite=False, key=None, just_the_path=False
    ):
        file_path = os.path.join(self.output_directory, "%s-%s" % (self.bin_id, prefix))

        if os.path.exists(file_path) and not overwrite:
            raise ConfigError(
                'get_output_file_handle: well, this file already exists: "%s"'
                % file_path
            )

        if not key:
            key = prefix.split(".")[0].replace("-", "_")

        self.bin_info_dict["files"][key] = file_path[
            len(self.summary.output_directory) :
        ].strip("/")

        return file_path if just_the_path else open(file_path, "w")

    def store_data_in_file(self, output_file_name_posfix, content):
        output_file_obj = self.get_output_file_handle(output_file_name_posfix)
        output_file_obj.write("%s\n" % content)
        output_file_obj.close()
