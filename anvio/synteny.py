#!/usr/bin/env python # -*- coding: utf-8
# pylint: disable=line-too-long

"""Classes to work with ngrams of contig functions.

These are classes to deconstruct loci into ngrams. They are used
to analyze conserved genes and synteny structures across loci.
"""

import sys
import pandas as pd

from collections import Counter

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import anvio.genomestorage as genomestorage

from anvio.dbops import PanDatabase
from anvio.errors import ConfigError

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Matthew Schechter"
__email__ = "mschechter@uchicago.edu"


class NGram(object):
    """class for counting NGrams

    anvi-analyze-synteny is designed to work with a group of similar loci, where each locus is a contig (which
    can lie in any number of contigs dbs)

    Parameters
    ==========
    args : argparse.Namespace
        For examples, arguments accepted by anvi-analyze-syntenty

    skip_sanity_check : bool, false
        If True, sanity_check will not be called.

    Notes
    =====
    - Currently the design assumes that each locus is a contig. In the future we have plans to expand this to
      compare genomes, not just loci, to each other. If that behavior is desired in the current design, each
      genome should be a single contig.
    """

    def __init__(
        self,
        args,
        run=terminal.Run(),
        progress=terminal.Progress(),
        skip_sanity_check=False,
    ):
        """Parses arguments and run sanity_check"""

        self.args = args
        self.run = run
        self.progress = progress

        # Parse arguments
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.annotation_source = A("annotation_source")
        self.window_range = A("ngram_window_range") or "2:3"
        self.is_in_unknowns_mode = A("analyze_unknown_functions")
        self.output_file = A("output_file")
        self.skip_init_functions = A("skip_init_functions")
        self.genome_names_to_focus = A("genome_names")
        self.ngram_source = A("ngram_source")
        self.first_functional_hit_only = A("first_functional_hit_only")
        self.annotation_source_dict = {}

        self.pan_db_path = A("pan_db")

        if self.annotation_source and self.pan_db_path:
            self.annotation_sources = [self.annotation_source, "gene_clusters"]

        if self.pan_db_path:
            filesnpaths.is_file_exists(self.pan_db_path)
            self.pan_db = PanDatabase(self.pan_db_path)

            self.p_meta = self.pan_db.meta

            self.p_meta["creation_date"] = (
                utils.get_time_to_date(self.p_meta["creation_date"])
                if "creation_date" in self.p_meta
                else "unknown"
            )
            self.p_meta["genome_names"] = sorted(
                [
                    s.strip()
                    for s in self.p_meta["external_genome_names"].split(",")
                    + self.p_meta["internal_genome_names"].split(",")
                    if s
                ]
            )
            self.p_meta["num_genomes"] = len(self.p_meta["genome_names"])
            self.genome_names = self.p_meta["genome_names"]
            self.gene_clusters_gene_alignments_available = self.p_meta[
                "gene_alignments_computed"
            ]
        else:
            self.pan_db = None

        self.genomes_storage_path = A("genomes_storage")

        # confirm genome-storage and pangenome hashes match of pangenome is provided
        if self.pan_db:
            self.genomes_storage = genomestorage.GenomeStorage(
                self.genomes_storage_path,
                self.p_meta["genomes_storage_hash"],
                genome_names_to_focus=self.p_meta["genome_names"],
                skip_init_functions=self.skip_init_functions,
                run=self.run,
                progress=self.progress,
            )
        else:
            self.genomes_storage = genomestorage.GenomeStorage(
                self.genomes_storage_path,
                skip_init_functions=self.skip_init_functions,
                run=self.run,
                progress=self.progress,
            )

        # list-annotation-resources
        self.list_annotation_sources = A("list_annotation_sources")
        self.gene_function_source_set = self.genomes_storage.db.get_table_as_dataframe(
            "gene_function_calls"
        ).source.unique()
        if self.list_annotation_sources:
            self.run.info(
                "Available functional annotation sources",
                ", ".join(self.gene_function_source_set),
            )
            sys.exit()

        # This houses the ngrams' data
        self.ngram_attributes_list = []

        # Focus on specfic set of genomes
        if self.genome_names_to_focus:
            if filesnpaths.is_file_exists(self.genome_names_to_focus, dont_raise=True):
                self.genome_names_to_focus = utils.get_column_data_from_TAB_delim_file(
                    self.genome_names_to_focus,
                    column_indices=[0],
                    expected_number_of_fields=1,
                )[0]
            else:
                self.genome_names_to_focus = [
                    g.strip() for g in self.genome_names_to_focus.split(",")
                ]

            self.run.warning(
                "A subset of genome names is found, and anvi'o will focus only on to those."
            )

            self.genomes_storage = genomestorage.GenomeStorage(
                self.genomes_storage_path,
                storage_hash=None,
                genome_names_to_focus=self.genome_names_to_focus,
            )
            self.genomes = self.genomes_storage.get_genomes_dict()

            self.external_genome_names = [
                g for g in self.genomes if self.genomes[g]["external_genome"]
            ]
            self.internal_genome_names = [
                g for g in self.genomes if not self.genomes[g]["external_genome"]
            ]

            self.hash_to_genome_name = {}
            for genome_name in self.genomes:
                self.hash_to_genome_name[self.genomes[genome_name]["genome_hash"]] = (
                    genome_name
                )

            # number of genomes in genome-storage
            self.num_contigs_in_external_genomes_with_genes = len(self.genomes)

        # number of genomes in genome-storage
        if self.genome_names_to_focus:
            self.num_contigs_in_external_genomes_with_genes = len(
                self.genome_names_to_focus
            )
        else:
            self.num_contigs_in_external_genomes_with_genes = len(
                self.genomes_storage.get_all_genome_names()
            )

        if not skip_sanity_check:
            self.sanity_check()

        # unless we are in debug mode, let's keep things quiet.
        if anvio.DEBUG:
            self.run_object = terminal.Run()
        else:
            self.run_object = terminal.Run(verbose=False)

    def sanity_check(self):
        """Sanity_check will confirm input for NGram class"""

        # checking if the annotation source is common across all contigs databases
        if (
            self.annotation_source
            and self.annotation_source not in self.gene_function_source_set
        ):
            raise ConfigError(
                "The annotation source you requested does not appear to be in all of "
                "the contigs databases from the external-genomes file. "
                "Please confirm your annotation-source and that all contigs databases have it :)"
            )

        if (self.annotation_source and self.pan_db) and not self.ngram_source:
            raise ConfigError(
                "anvi-analyze-synteny needs to know which annotation source to slice Ngrams with. "
                "Please use the --ngram-source flag to declare one :)"
            )

        if not self.args.output_file:
            raise ConfigError("You should provide an output file name.")

        # checking window-range input
        if self.window_range.count(":") != 1:
            raise ConfigError(
                "anvi'o would love to slice and dice your loci, but the "
                "Format of window_range must be x:y (e.g. Window sizes 2 to 4 would be denoted as: 2:4)"
            )

        try:
            self.window_range = [int(n) for n in self.window_range.split(":")]
        except ValueError:
            raise ConfigError(
                "anvi'o would love to slice and dice your loci, but the "
                "window-ranges need to be integers :)"
            )

        if self.window_range[0] > self.window_range[1]:
            raise ConfigError(
                "anvi'o would love to slice and dice your loci, but the "
                "window-range needs to be from small to big :)"
            )

        # Window-range must contain 2 integers for window
        if (
            len(self.window_range) > 2
            or not isinstance(self.window_range[0], int)
            or not isinstance(self.window_range[1], int)
        ):
            raise ConfigError(
                "anvi'o would love to slice and dice your loci, but... the "
                "window_range must only contain 2 integers and be formated as x:y (e.g. Window sizes 2 to 4 would be denoted as: 2:4)"
            )

        # Loop through each contigs db, test that each contig contains at least as many genes as max window size and confirm every contig has annotations
        for contigs_db_name in self.genomes_storage.get_genomes_dict():

            gene_caller_ids = self.genomes_storage.get_gene_caller_ids(contigs_db_name)
            num_genes = len(gene_caller_ids)

            if self.window_range[1] > num_genes:
                raise ConfigError(
                    "The largest window size you requested (%d) is larger than the number of genes found on this genome: %s"
                    % (self.window_range[1], contigs_db_name)
                )

    def populate_ngram_attributes(self):
        """Iterates through all contigs and use self.count_synteny to count all ngrams in that contig.

        This populates the self.ngram_attributes_list, where each element looks like:
        (ngram, count, contigs_db_name, contig_name, N)
        """

        # Get gene cluster info from panDB
        if self.pan_db:
            gene_cluster_frequencies_dataframe = self.pan_db.db.get_table_as_dataframe(
                "gene_clusters"
            )

        self.run.warning(
            "Anvi'o is now looking for Ngrams in your contigs!", lc="green"
        )

        self.run.info_single(
            "What do we say to loci that appear to have no coherent synteny patterns...? Not today! ⚔️",
            nl_before=1,
            nl_after=1,
        )

        for contigs_db_name in self.genomes_storage.get_genomes_dict():

            # Get list of genes-callers-ids
            gene_caller_ids_list = list(
                self.genomes_storage.get_gene_caller_ids(contigs_db_name)
            )

            # Use gene_caller_ids_list to get functions table
            gene_function_call_df = self.genomes_storage.db.get_table_as_dataframe(
                "gene_function_calls"
            )

            # Create dict for annotate Ngrams
            if self.annotation_source:
                self.gene_caller_id_to_accession_dict = (
                    self.get_genes_and_functions_dict(
                        contigs_db_name, gene_function_call_df
                    )
                )
                self.annotation_source_dict[self.annotation_source] = (
                    self.gene_caller_id_to_accession_dict
                )

            if self.pan_db:
                self.gene_caller_id_to_gene_cluster_dict = self.get_gene_cluster_dict(
                    contigs_db_name, gene_cluster_frequencies_dataframe
                )
                self.annotation_source_dict["gene_clusters"] = (
                    self.gene_caller_id_to_gene_cluster_dict
                )

            # Iterate over range of window sizes and run synteny algorithm to count occurrences of ngrams in a contig
            for N in range(*self.window_range):
                ngram_counts_dict, annotations_dict = self.count_synteny(
                    N, gene_caller_ids_list
                )

                for ngram, count in ngram_counts_dict.items():

                    for annotation_source, annotation in annotations_dict[
                        ngram
                    ].items():
                        if annotation_source == self.ngram_source:
                            pass
                        else:
                            self.ngram_attributes_list.append(
                                [
                                    ngram,
                                    count,
                                    annotations_dict[ngram][annotation_source],
                                    contigs_db_name,
                                    N,
                                ]
                            )

    def count_synteny(self, N, gene_caller_ids_list):
        """This method counts synteny patterns of size N on a contig

        This method counts synteny patterns of size N on a contig by taking a window of
        gene-callers-ids, annotating it, ordering it consistently, and then counting its occurence
        on the contig.

        Parameters
        ==========
        n : int
            A window size to extract a ngram

        gene_caller_ids_list: list
            list of all gene-callers-ids  on a contig

        Returns
        =======
        ngram_counts_dict : dict
            A dict of ngram counts on a contig A tuple of annotations {ngram:count}

        Notes
        =====
        This function assumes that the input list of gene-callers-ids are in the order of which they
        occur on the contig
        """

        ngram_counts_dict = Counter({})
        annotations_dict = {}
        gene_callers_id_windows = self.get_windows(N, gene_caller_ids_list)

        for window in gene_callers_id_windows:

            annotated_window_dict = self.annotate_window(window)
            if self.ngram_source:
                ngram = self.order_window(annotated_window_dict[self.ngram_source])
            elif self.annotation_source:
                ngram = self.order_window(annotated_window_dict[self.annotation_source])
            else:
                ngram = self.order_window(annotated_window_dict["gene_clusters"])

            # flip annotation if the ngram was flipped
            if ngram[1] == True:
                annotated_window_dict_ordered = {}
                for annotation_source, annotation in annotated_window_dict.items():
                    annotated_window_flipped = annotation[::-1]
                    annotated_window_dict_ordered[annotation_source] = (
                        annotated_window_flipped
                    )
                    annotations_dict[ngram[0]] = (
                        annotated_window_dict_ordered  # record flipped version of annotation
                    )
            else:
                annotations_dict[ngram[0]] = annotated_window_dict

            ngram_counts_dict[ngram[0]] += 1

        return ngram_counts_dict, annotations_dict

    def annotate_window(self, window):
        """This method annotates a gene-callers-id window

        This method will annotate a gene-callers-id window based using annotation sources provided
        by the user (e.g. COGs, pan_db). If the user provided the `first_functional_hit_only` flag, the COG
        annotation will be split by "!!!" and the first (best hit) item will be used.

        Parameters
        ==========
        window : tuple
            A tuple of gene gene-callers-ids that represents an unannotated Ngram

        Returns
        =======
        ngram_gene_clusters : tuple
            A tuple of annotations
        """

        # Annotate window based on user input
        gene_annotation_dict = {}
        for annotation_source, annotations_dict in self.annotation_source_dict.items():
            if self.first_functional_hit_only:
                ngram_annotation = []
                for g in window:
                    annotation = annotations_dict[g]
                    if "!!!" in annotation:
                        annotation_first = annotation.split("!!!")[0]
                        ngram_annotation.append(annotation_first)
                    else:
                        ngram_annotation.append(annotation)

                ngram_annotation = tuple(ngram_annotation)

            else:
                ngram_annotation = tuple([annotations_dict[g] for g in window])

            gene_annotation_dict[annotation_source] = ngram_annotation

        return gene_annotation_dict

    def order_window(self, annotated_window):
        """This method orients an annotated_window in a consistent order

        This method is to make sure that all Ngrams are ordered in the same direction. For example,
        we want to count A-B-C and C-B-A as the same Ngram.

        Parameters
        ==========
        annotated_window : tuple
            A tuple of annotations accessions

        Returns
        =======
        ngram : tuple
            A tuple Ngram that is correctly oriented
        """

        original_order = annotated_window
        flipped_order = annotated_window[::-1]

        if original_order[0] < flipped_order[0]:
            ngram = original_order
            flipped = False
        else:
            ngram = flipped_order
            flipped = True

        ngram = [ngram, flipped]

        return ngram

    def convert_to_df(self):
        """Takes self.ngram_attributes_list and returns a pandas dataframe"""

        ngram_count_df_list = []
        for ngram_attribute in self.ngram_attributes_list:
            ngram = "::".join(map(str, list(ngram_attribute[0])))
            annotation = "::".join(map(str, list(ngram_attribute[2])))
            if self.pan_db and self.annotation_source:
                df = pd.DataFrame(
                    columns=[
                        "ngram",
                        "count",
                        "annotation",
                        "contig_db_name",
                        "N",
                        "number_of_loci",
                    ]
                )
                df = df.append(
                    {
                        "ngram": ngram,
                        "count": ngram_attribute[1],
                        "annotation": annotation,
                        "contig_db_name": ngram_attribute[3],
                        "N": ngram_attribute[4],
                        "number_of_loci": self.num_contigs_in_external_genomes_with_genes,
                    },
                    ignore_index=True,
                )
            elif self.pan_db and not self.annotation_source:
                ngram = "::".join(map(str, list(ngram_attribute[0])))
                df = pd.DataFrame(
                    columns=["ngram", "count", "contig_db_name", "N", "number_of_loci"]
                )
                df = df.append(
                    {
                        "ngram": ngram,
                        "count": ngram_attribute[1],
                        "contig_db_name": ngram_attribute[3],
                        "N": ngram_attribute[4],
                        "number_of_loci": self.num_contigs_in_external_genomes_with_genes,
                    },
                    ignore_index=True,
                )
            else:
                ngram = "::".join(map(str, list(ngram_attribute[0])))
                df = pd.DataFrame(
                    columns=["ngram", "count", "contig_db_name", "N", "number_of_loci"]
                )
                df = df.append(
                    {
                        "ngram": ngram,
                        "count": ngram_attribute[1],
                        "contig_db_name": ngram_attribute[3],
                        "N": ngram_attribute[4],
                        "number_of_loci": self.num_contigs_in_external_genomes_with_genes,
                    },
                    ignore_index=True,
                )

            ngram_count_df_list.append(df)

        ngram_count_df_final = pd.concat(ngram_count_df_list)

        if not self.is_in_unknowns_mode:
            ngram_count_df_final = ngram_count_df_final[
                ~ngram_count_df_final["ngram"].str.contains(
                    "unknown-function" or "no-gene-cluster-annotation"
                )
            ]

        return ngram_count_df_final

    def report_ngrams_to_user(self):
        """Counts ngrams per contig and reports as tab-delimited file"""

        self.populate_ngram_attributes()
        df = self.convert_to_df()
        df.to_csv(self.output_file, sep="\t", index=False)
        self.run.info("Ngram table", self.output_file)

    def get_genes_and_functions_dict(self, contigs_db_name, gene_function_call_df):
        """This method will extract a list of gene attributes from each contig within a contigsDB.

        Returns
        =======
        output : list of lists
            first element is gene_caller_id, second is function accession, third is the contig name
        """

        # get contigsDB
        gene_function_call_df_filtered = gene_function_call_df[
            (gene_function_call_df["genome_name"] == contigs_db_name)
            & (gene_function_call_df["source"] == self.annotation_source)
        ]
        gene_callers_id_to_accession_dict = (
            gene_function_call_df_filtered[["gene_callers_id", "accession"]]
            .set_index("gene_callers_id")["accession"]
            .to_dict()
        )

        gene_caller_ids_list = self.genomes_storage.get_gene_caller_ids(contigs_db_name)

        # Make list of lists containing gene attributes. If there is not annotation add one in!
        genes_and_functions_list = (
            []
        )  # List of lists [gene-caller-id, accessions, contig-name]
        counter = 0
        for gene_callers_id in gene_caller_ids_list:
            list_of_gene_attributes = []

            if gene_callers_id in gene_callers_id_to_accession_dict:
                accession = gene_callers_id_to_accession_dict[gene_callers_id]
                accession = accession.replace(" ", "")
                list_of_gene_attributes.extend((gene_callers_id, accession))
                genes_and_functions_list.append(list_of_gene_attributes)
            else:
                # adding in "unknown annotation" if there is none
                accession = "unknown-function"
                list_of_gene_attributes.extend((counter, accession))
                genes_and_functions_list.append(list_of_gene_attributes)
            counter += 1

        gene_caller_id_to_accession_dict = {}
        for entry in genes_and_functions_list:
            gene_caller_id_to_accession_dict[entry[0]] = entry[1]

        return gene_caller_id_to_accession_dict

    def get_gene_cluster_dict(
        self, contigs_db_name, gene_cluster_frequencies_dataframe
    ):
        gene_cluster_frequencies_dataframe_filtered = (
            gene_cluster_frequencies_dataframe[
                gene_cluster_frequencies_dataframe["genome_name"] == contigs_db_name
            ]
        )
        gene_callers_id_to_gene_cluster_id_dict = (
            gene_cluster_frequencies_dataframe_filtered[
                ["gene_caller_id", "gene_cluster_id"]
            ]
            .set_index("gene_caller_id")["gene_cluster_id"]
            .to_dict()
        )

        gene_caller_ids_list = self.genomes_storage.get_gene_caller_ids(contigs_db_name)

        # Make list of lists containing gene cluster attributes. If there is not annotation add one in!
        genes_cluster_list = (
            []
        )  # List of lists [gene-caller-id, gene-cluster-id, contig-name]
        counter = 0
        for gene_callers_id in gene_caller_ids_list:
            list_of_gene_attributes = []

            if gene_callers_id in gene_callers_id_to_gene_cluster_id_dict:
                gene_cluster_id = gene_callers_id_to_gene_cluster_id_dict[
                    gene_callers_id
                ]
                gene_cluster_id = gene_cluster_id.replace(" ", "")
                list_of_gene_attributes.extend((gene_callers_id, gene_cluster_id))
                genes_cluster_list.append(list_of_gene_attributes)
            else:
                # adding in "unknown annotation" if there is none
                gene_cluster_id = "no-gene-cluster-annotation"
                list_of_gene_attributes.extend((counter, gene_cluster_id))
                genes_cluster_list.append(list_of_gene_attributes)
            counter += 1

        gene_caller_id_to_gene_cluster_dict = {}
        for entry in genes_cluster_list:
            gene_caller_id_to_gene_cluster_dict[entry[0]] = entry[1]

        return gene_caller_id_to_gene_cluster_dict

    def get_windows(self, N, gene_caller_ids_list):
        """This method will count NGrams in contigs

        This method will use a sliding window of size N to extract gene-callers-id windows.
        The final output will be a list of windows of size N.

        Parameters
        ==========
        gene_caller_ids_list : list
            A list of gene gene-callers-ids as they appear in the contig

        n : int
            A window size to extract a ngram

        Returns
        =======
        gene_callers_id_windows : list
            A list of n sized windows of gene-callers-ids from a contig
        """

        gene_callers_id_windows = []
        for i in range(0, len(gene_caller_ids_list) - N + 1):

            # extract window
            window = tuple(gene_caller_ids_list[i : i + N])

            gene_callers_id_windows.append(window)

        return gene_callers_id_windows
