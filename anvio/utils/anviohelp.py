import shutil

import numpy as np
import pandas as pd

from collections import Counter

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.constants as constants
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.utils.algorithms import get_list_of_outliers
from anvio.dbinfo import is_contigs_db, is_pan_or_profile_db
from anvio.utils.database import get_required_version_for_db

# for full output
pd.options.display.max_columns=100
pd.options.display.max_rows=100

class CoverageStats:
    """A class to return coverage stats for an array of nucleotide level coverages.

    FIXME: This class should replace `coverage_c` function in bamops to avoid redundancy.
    """

    def __init__(self, coverage, skip_outliers=False):
        self.min = np.amin(coverage)
        self.max = np.amax(coverage)
        self.median = np.median(coverage)
        self.mean = np.mean(coverage)
        self.std = np.std(coverage)
        self.detection = np.sum(coverage > 0) / len(coverage)

        if coverage.size < 4:
            self.mean_Q2Q3 = self.mean
        else:
            sorted_c = np.sort(coverage)
            Q = int(coverage.size * 0.25)
            Q2Q3 = sorted_c[Q:-Q]
            self.mean_Q2Q3 = np.mean(Q2Q3)

        if skip_outliers:
            self.is_outlier = None
        else:
            self.is_outlier = get_list_of_outliers(coverage, median=self.median) # this is an array not a list


def ununique_BLAST_tabular_output(tabular_output_path, names_dict):
    """FIXME <A one line descriptor here>

    Notes
    =====
    - Assumes outfmt has `qseqid` and `sseqid` as 1st and 2nd columns, respectively
    """

    new_search_output_path = tabular_output_path + '.ununiqued'
    new_tabular_output = open(new_search_output_path, 'w')

    for line in open(tabular_output_path):
        fields = line.strip().split('\t')
        for query_id in names_dict[fields[0]]:
            for subject_id in names_dict[fields[1]]:
                new_tabular_output.write('%s\t%s\t%s\n' % (query_id, subject_id, '\t'.join(fields[2:])))


    new_tabular_output.close()

    shutil.move(tabular_output_path, tabular_output_path + '.unique')
    shutil.move(new_search_output_path, tabular_output_path)

    return tabular_output_path



def convert_SSM_to_single_accession(matrix_data):
    """
    The substitution scores from the SSM dictionaries created in anvio.data.SSMs are accessed via a dictionary of
    dictionaries, e.g.  data["Ala"]["Trp"]. This returns a new dictionary accessed via the concatenated sequence element
    pair, e.g. data["AlaTrp"], data["AT"], etc.  where they are ordered alphabetically.
    """
    items = matrix_data.keys()
    new_data = {}

    for row in items:
        for column in items:

            if row > column:
                continue
            new_data[''.join([row, column])] = matrix_data[row][column]
    return new_data



def get_default_gene_caller(contigs_db_path):
    """Returns the default gene caller, but in a smart way.

    Well. Smart is not a very accurate way to say this. Here is some history. For the longest time,
    anvi'o used `prodigal` as its default gene caller. So we had this one entry in anvio/constants.py
    that said `default_gene_caller = 'prodigal'`, and life was simple. Then it became clear that
    we needed to switch to pyrodigal-gv, which was an effort to maintain the stale repository of
    prodigal, but also included some models for giant viruses. But we couldn't simply change
    `default_gene_caller` to `pyrodigal-gv`, since there were many contigs-db files out there with
    prodigal gene calls, which are equally fine. But why is this a problem?

    It is a problem when the user wants to build a pangenome, or export amino acid sequences for
    a given contigs-db file. In those cases, anvi'o needs to know which default gene caller it
    should use, and only if it can't find that should it ask the user to choose a gene caller
    explicitly. Removing `prodigal` from the constants file was going to make users' life very
    difficult when they were using a newer version of anvi'o (with `pyrodigal-gv` as the default
    gene caller in the constants.py) but using contigs-db files generated with older versions.

    The purpose of this function is to solve that problem by retrieving the 'default' gene caller
    by considering all gene callers listed in constants.py, and taking a look at the most
    frequent source of gene calls in the contigs-db file. If `prodigal` is the most frequent
    gene caller, and if it is still in the list of default gene callers in constants.py, then
    this function would return `prodigal`. The same for `pyrodigal-gv`. Only after checking
    for all recognized default gene callers would the relevant context ask the user to provide
    a gene caller for downstream operations.

    Parameters
    ==========
    contigs_db_path : str
        Path to an anvi'o contigs-db file

    Returns
    =======
    value : str
        None if most frequent gene caller source in contigs_db is not in constants.py, else
        the gene caller source as default
    """

    is_contigs_db(contigs_db_path)

    contigs_db = db.DB(contigs_db_path, anvio.__contigs__version__)

    gene_call_sources_in_contigs_db = contigs_db.get_single_column_from_table(t.genes_in_contigs_table_name, 'source')

    try:
        most_frequent_gene_caller = Counter(gene_call_sources_in_contigs_db).most_common(1)[0][0]
    except IndexError:
        most_frequent_gene_caller = None

    contigs_db.disconnect()

    if most_frequent_gene_caller in constants.default_gene_callers:
        return most_frequent_gene_caller
    else:
        return None



def get_split_and_contig_names_of_interest(contigs_db_path, gene_caller_ids):
    """Takes a set of gene caller ids, returns all split and contig names in a
       contigs database that are affiliated with them.
    """

    if not isinstance(gene_caller_ids, set):
        raise ConfigError("`gene_caller_ids` must be of type `set`.")

    is_contigs_db(contigs_db_path)

    contigs_db = db.DB(contigs_db_path, anvio.__contigs__version__)

    where_clause_genes = "gene_callers_id in (%s)" % ', '.join(['%d' % g for g in gene_caller_ids])
    genes_in_contigs = contigs_db.get_some_rows_from_table_as_dict(t.genes_in_contigs_table_name, where_clause=where_clause_genes)
    contig_names_of_interest = set([e['contig'] for e in genes_in_contigs.values()])

    where_clause_contigs = "parent in (%s)" % ', '.join(['"%s"' % c for c in contig_names_of_interest])
    splits_info = contigs_db.get_some_rows_from_table_as_dict(t.splits_info_table_name, where_clause=where_clause_contigs)
    split_names_of_interest = set(splits_info.keys())

    contigs_db.disconnect()

    return (split_names_of_interest, contig_names_of_interest)



def get_contigs_splits_dict(split_ids, splits_basic_info):
    """
    For a given list of split ids, create a dictionary of contig names
    that represents all parents as keys, and ordered splits as items.

    split_ids is a set of split IDs, splits_basic_info comes from the contigs database:

     >>> contigs_db = dbops.ContigsDatabase(contigs_db_path)
     >>> splits_basic_info = contigs_db.db.get_table_as_dict(t.splits_info_table_name)
     >>> znnotation_db.disconnect()
     >>> x = get_contigs_splits_dict(set([contig_A_split_00001, contig_A_split_00002, contig_A_split_00004,
                                         contig_C_split_00003, contig_C_split_00004, contig_C_split_00005]),
                                    splits_basic_info)
     >>> print x
         {
             'contig_A': {
                             0: 'contig_A_split_00001',
                             1: 'contig_A_split_00002',
                             4: 'contig_A_split_00004'
                         },
             'contig_C': {
                             3: 'contig_C_split_00003',
                             4: 'contig_C_split_00004',
                             5: 'contig_C_split_00005'
                         }
         }
    """

    contigs_splits_dict = {}

    for split_id in split_ids:
        s = splits_basic_info[split_id]
        if s['parent'] in contigs_splits_dict:
            contigs_splits_dict[s['parent']][s['order_in_parent']] = split_id
        else:
            contigs_splits_dict[s['parent']] = {s['order_in_parent']: split_id}

    return contigs_splits_dict



def get_variabile_item_frequencies(e, engine='NT'):
    """
    e is a row from variable_nucleotide_positions table defined in tables.
    this function extends dictionary with consensus and departure from consensus.
    """

    items = constants.nucleotides if engine=='NT' else constants.amino_acids
    frequency_dict = Counter(dict([(item, e[item]) for item in items]))
    return frequency_dict.most_common()



def get_consensus_and_departure_data(variable_item_frequencies):
    """Make sense of `variable_item_frequencies`.

       The format of `variable_item_frequencies` follows this:

           >>> [('A', 45), ('T', 5), ('G', 0), ('N', 0), ('C', 0)]

       For a given entry of the variable_XX_frequencies table, the `variable_item_frequencies`
       tuple can be obtained via `get_variabile_item_frequencies`.
    """

    frequency_of_consensus = variable_item_frequencies[0][1]
    total_frequency_of_all_but_the_consensus = sum([tpl[1] for tpl in variable_item_frequencies[1:]])
    coverage = total_frequency_of_all_but_the_consensus + frequency_of_consensus

    n2n1ratio = variable_item_frequencies[1][1] / frequency_of_consensus if frequency_of_consensus else -1
    consensus = variable_item_frequencies[0][0]
    departure_from_consensus = total_frequency_of_all_but_the_consensus / coverage if coverage else -1

    return (n2n1ratio, consensus, departure_from_consensus)



def get_bin_name_from_item_name(anvio_db_path, item_name, collection_name=None):
    is_pan_or_profile_db(anvio_db_path, genes_db_is_also_accepted=True)
    database = db.DB(anvio_db_path, None, ignore_version=True)

    if t.collections_splits_table_name not in database.get_table_names():
        raise ConfigError("The database %s does not contain a collections table :/")

    if collection_name:
        where_clause = 'split = "%s" and collection_name = "%s"' % (item_name, collection_name)
    else:
        where_clause = 'split = "%s"' % (item_name)

    rows = database.get_some_rows_from_table(t.collections_splits_table_name, where_clause=where_clause)

    database.disconnect()

    return rows



def get_contig_name_to_splits_dict(contigs_db_path):
    """Returns a dict for contig name to split name conversion"""

    is_contigs_db(contigs_db_path)

    contigs_db = db.DB(contigs_db_path, get_required_version_for_db(contigs_db_path))

    contig_name_to_split_names_dict = {}
    for split_name, contig_name in contigs_db.get_some_columns_from_table(t.splits_info_table_name, "split, parent"):
        if contig_name not in contig_name_to_split_names_dict:
            contig_name_to_split_names_dict[contig_name] = set([])

        contig_name_to_split_names_dict[contig_name].add(split_name)

    return contig_name_to_split_names_dict



def get_variability_table_engine_type(table_path, dont_raise=False):
    """A non-extensive test to determine if a file was generated by anvi-gen-variability-profile,
       and if it was, what engine (NT, CDN, or AA) was used.
    """
    filesnpaths.is_file_tab_delimited(table_path)
    columns_names = set(pd.read_csv(table_path, sep="\t", nrows = 0).columns)

    if set(constants.nucleotides) < columns_names:
        return "NT"

    elif set(constants.codons) < columns_names:
        return "CDN"

    elif set(constants.amino_acids) < columns_names:
        return "AA"

    else:
        if dont_raise:
            return ""
        raise ConfigError("anvi'o does not recognize %s as being a variability table generated by "
                          "anvi-gen-variability-profile." % table_path)
