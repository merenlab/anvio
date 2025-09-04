import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.terminal import Run, Progress, SuppressAllOutput


def get_names_order_from_newick_tree(newick_tree, newick_format=1, reverse=False, names_with_only_digits_ok=False):
    # import ete3
    with SuppressAllOutput():
        from ete3 import Tree

    filesnpaths.is_proper_newick(newick_tree, names_with_only_digits_ok=names_with_only_digits_ok)

    tree = Tree(newick_tree, format=newick_format)

    names = [n.name for n in tree.get_leaves()]

    return list(reversed(names)) if reverse else names



def gen_NEXUS_format_partition_file_for_phylogenomics(partition_file_path, sequence_lengths, separator='', run=Run(), progress=Progress()):
    """ Generates a NEXUS-formatted partition file for phylogenomics. See
        https://github.com/merenlab/anvio/issues/1333 for details

    Parameters
    ==========
    partition_file_path: `str`
        File path to be generated.
    sequence_lengths: `list` of `tuples`
        A list that contins sequence names and lenghts as tuples. I.e.,
        [('seq_1', 100), ('seq_2', 42), ...]
    separator: `str`
        Characters used to separate sequences from each other in a multi-alignment
        file.
    run: `object`
        Anvi'o run object
    run: `progress`
        Anvi'o progress object

    Returns
    =======
    None
    """

    filesnpaths.is_output_file_writable(partition_file_path)

    if not isinstance(sequence_lengths, list):
        raise ConfigError("Sequence lengths must be passed as a list of tuples.")

    if not isinstance(sequence_lengths[0], tuple):
        raise ConfigError("Sequence lengths must be passed as a list of tuples.")

    with open(partition_file_path, 'w') as partition_file:
        partition_file.write("#nexus\nbegin sets;\n")
        index = 1
        for sequence_name, sequence_length in sequence_lengths:
            partition_file.write("    charset %s = %d-%d;\n" % (sequence_name, index, index + sequence_length - 1))
            index += (sequence_length + len(separator))
        partition_file.write("end;\n")

    progress.reset()
    run.info("Partition file", partition_file_path, mc='yellow')
    run.info_single("Your partition file is ready. Please do not forget to replace placeholders for model names ('[MODEL]') "
                    "in this file with appropriate model names prior to your phylogenomic analysis.", nl_before=1, nl_after=1)

