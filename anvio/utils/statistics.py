import os

import numpy as np

import anvio
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.terminal import Run, Progress

from anvio.utils.commandline import run_command
from anvio.utils.system import is_program_exists
from anvio.utils.files import get_TAB_delimited_file_as_dictionary


def run_functional_enrichment_stats(functional_occurrence_stats_input_file_path, enrichment_output_file_path=None, run=Run(), progress=Progress()):
    """This function runs the enrichment analysis implemented by Amy Willis.

    Since the enrichment analysis is an R script, we interface with that program by
    producing a compatible input file first, and then calling this function from various
    places in the anvi'o code.

    Parameters
    ==========
    functional_occurrence_stats_input_file_path, str file path
        This is the primary input file for the R script, `anvi-script-enrichment-stats`.
        For the most up-do-date file header, please see the header section of the R
        script.
    enrichment_output_file_path, str file path
        An optional output file path for the enrichment analysis.

    Returns
    =======
    enrichment_output: dict
        The enrichment analysis results
    """

    run.warning("This program will compute enrichment scores using an R script developed by Amy Willis. "
                "You can find more information about it in the following paper: Shaiber, Willis et al "
                "(https://doi.org/10.1186/s13059-020-02195-w). When you publish your findings, please "
                "do not forget to properly credit this work. :)", lc='green', header="CITATION")

    # sanity check for R packages
    package_dict = get_required_packages_for_enrichment_test()
    check_R_packages_are_installed(package_dict)

    # make sure the input file path is a TAB delmited file that exists.
    filesnpaths.is_file_tab_delimited(functional_occurrence_stats_input_file_path)

    if not enrichment_output_file_path:
        enrichment_output_file_path = filesnpaths.get_temp_file_path()
    elif filesnpaths.is_file_exists(enrichment_output_file_path, dont_raise=True):
        raise ConfigError(f"The file {enrichment_output_file_path} already exists and anvi'o doesn't like to overwrite it :/ "
                          f"Please either delete the existing file, or provide another file path before re-running this "
                          f"program again.")

    log_file_path = filesnpaths.get_temp_file_path()

    run.warning(None, header="AMY's ENRICHMENT ANALYSIS 🚀", lc="green")
    run.info("Functional occurrence stats input file path: ", functional_occurrence_stats_input_file_path)
    run.info("Functional enrichment output file path: ", enrichment_output_file_path)
    run.info("Temporary log file (use `--debug` to keep): ", log_file_path, nl_after=2)

    # run enrichment script
    progress.new('Functional enrichment analysis')
    progress.update("Running Amy's enrichment")
    run_command(['anvi-script-enrichment-stats',
                 '--input', f'{functional_occurrence_stats_input_file_path}',
                 '--output', f'{enrichment_output_file_path}'], log_file_path)
    progress.end()

    if not filesnpaths.is_file_exists(enrichment_output_file_path, dont_raise=True):
        raise ConfigError(f"Something went wrong during the functional enrichment analysis :( We don't "
                          f"know what happened, but this log file could contain some clues: {log_file_path}")

    if filesnpaths.is_file_empty(enrichment_output_file_path):
        raise ConfigError(f"Something went wrong during the functional enrichment analysis :( "
                          f"An output file was created, but it was empty... We hope that this "
                          f"log file offers some clues: {log_file_path}")

    # if everything went okay, we remove the log file
    if anvio.DEBUG:
        run.warning(f"Due to the `--debug` flag, anvi'o keeps the log file at '{log_file_path}'.", lc='green', header="JUST FYI")
    else:
        os.remove(log_file_path)

    enrichment_stats = get_TAB_delimited_file_as_dictionary(enrichment_output_file_path)

    # here we will naively try to cast every column that matches `p_*` to float, and every
    # column that matches `N_*` to int.
    column_names = list(enrichment_stats.values())[0].keys()
    column_names_to_cast = [(c, float) for c in ['unadjusted_p_value', 'adjusted_q_value', 'enrichment_score']] + \
                           [(c, float) for c in column_names if c.startswith('p_')] + \
                           [(c, int) for c in column_names if c.startswith('N_')]
    for entry in enrichment_stats:
        for column_name, to_cast in column_names_to_cast:
            try:
                enrichment_stats[entry][column_name] = to_cast(enrichment_stats[entry][column_name])
            except:
                raise ConfigError(f"Something sad happened :( Anvi'o expects the functional enrichment output to contain "
                                  f"values for the column name `{column_name}` that can be represented as `{to_cast}`. Yet, the "
                                  f"entry `{entry}` in your output file contained a value of `{enrichment_stats[entry][column_name]}`. "
                                  f"We have no idea how this happened, but it is not good :/ If you would like to mention this "
                                  f"to someone, please attach to your inquiry the following file: '{enrichment_output_file_path}'.")

    return enrichment_stats


# # FIXME
# def is_external_genomes_compatible_with_pan_database(pan_db_path, external_genomes_path):


def get_enriched_groups(props, reps):
    '''
        Accepts a vector of proportions and number of replicates per group and
        returns a boolean vector where each group that has proportion above
        the "expected" (i.e. the overall proportion) is True and the rest are False.
    '''
    # if the function doesn't occur at all then test_statistic is zero and p-value is 1
    if not np.count_nonzero(props):
        return np.zeros(len(props))
    overall_portion = np.sum(np.multiply(props, reps)) / np.sum(reps)

    return props > overall_portion



def get_required_packages_for_enrichment_test():
    ''' Return a dict with the packages as keys and installation instrucstions as values'''
    packages = ["tidyverse", "stringi", "magrittr", "qvalue", "optparse"]

    installation_instructions = ["conda install -c r r-tidyverse",
                                 "conda install -c r r-stringi",
                                 "conda install -c bioconda r-magrittr",
                                 "conda install -c bioconda bioconductor-qvalue",
                                 "conda install -c conda-forge r-optparse"]

    return dict(zip(packages,installation_instructions))



def check_R_packages_are_installed(required_package_dict):
    """Checks if R and the provided R packages are installed on the user's system.
    If not, raises an error with installation instructions for any missing packages.

    Credits to Ryan Moore (https://github.com/mooreryan) for this solution!
    (https://github.com/merenlab/anvio/commit/91f9cf1531febdbf96feb74c3a68747b91e868de#r35353982)

    Parameters
    ==========
    required_package_dict, dictionary
        keys should be R package names, values should be the corresponding installation instruction for the package
        See get_required_packages_for_enrichment_test() for an example
    """

    is_program_exists('Rscript')

    missing_packages = []
    log_file = filesnpaths.get_temp_file_path()
    for lib in required_package_dict:
        ret_val = run_command(["Rscript", "-e", "library('%s')" % lib], log_file)
        if ret_val != 0:
            missing_packages.append(lib)

    if missing_packages:
        if len(missing_packages) == 1 and 'qvalue' in missing_packages:
            raise ConfigError("It seems you're struggling with the R package `qvalue`. It can be a pain to install. In our experience "
                              "best way to install this package is to do it through Bioconductor directly. For that, please "
                              "copy-paste this command as a single line into your terminal and run it: "
                              "Rscript -e 'install.packages(\"BiocManager\", repos=\"https://cran.rstudio.com\"); BiocManager::install(\"qvalue\")'")
        else:
            raise ConfigError("The following R packages are required in order to run this, but seem to be missing or broken: '%(missing)s'. "
                              "If you have installed anvi'o through conda, BEFORE ANYTHING ELSE we would suggest you to run the command "
                              "Rscript -e \"update.packages(repos='https://cran.rstudio.com')\" in your terminal. This will try to update "
                              "all R libraries on your conda environment and will likely solve this problem. If it doesn't work, then you "
                              "will need to try a bit harder, so here are some pointers: if you are using conda, in an ideal world you"
                              "should be able to install these packages by running the following commands: %(conda)s. But if this option "
                              "doesn't seem to be working for you, then you can also try to install the problem libraries directly through R, "
                              "for instance by typing in your terminal, Rscript -e 'install.packages(\"%(example)s\", "
                              "repos=\"https://cran.rstudio.com\")' and see if it will address the installation issue. UNFORTUNATELY, in "
                              "some cases you may continue to see this error despite the fact that you have these packages installed :/ It "
                              "would most likely mean that some other issues interfere with their proper usage during run-time. If you have "
                              "these packages installed but you continue seeing this error, please run in your terminal Rscript -e "
                              "\"library(%(example)s)\" to see what is wrong with %(example)s on your system. Running this on your "
                              "terminal will test whether the package is properly loading or not and the resulting error messages will likely "
                              "be much more helpful solving the issue. If none of the solutions offered here worked for you, feel free to "
                              "come to anvi'o Discord and ask around -- others may already have a solution for it already. Apologies for the "
                              "frustration. R frustrates everyone." % {'missing': ', '.join(missing_packages),
                                                                       'conda': ', '.join(['"%s"' % required_package_dict[i] for i in missing_packages]),
                                                                       'example': missing_packages[0]})
    else:
        os.remove(log_file)



def get_values_of_gene_level_coverage_stats_as_dict(gene_level_coverage_stats_dict, key, genes_of_interest=None, samples_of_interest=None, as_pandas=False):
    """
        This function takes the gene_level_coverage_stats_dict and return one of the values
        as a matrix-like dict of dicts.
        THIS FUNCTION IS IN utils AND NOT IN summarizer, or dbops, because it used to be in summarizer
        and why should it be in summarizer?!? that makes no sense. And also mcg-classifier doesn't want
        to initialize summarizer, it wants to be able to just get the gene_level_coverage_stats_dict as
        input and then deal with it.

        There is also an option to as to get the data back as a pandas dataframe.
    """
    legal_keys = {'mean_coverage', 'detection', 'non_outlier_mean_coverage', 'non_outlier_coverage_std'}
    if key not in legal_keys and as_pandas:
        raise ConfigError("%s is not a valid key for creating a pandas dataframe of values of gene_level_coverage_stats_dict. "
                           "Here is a list of the valid keys: %s" % (key, list(legal_keys)))

    gene_callers_ids = set(gene_level_coverage_stats_dict.keys())
    samples = set(next(iter(gene_level_coverage_stats_dict.values())).keys())

    if genes_of_interest is not None:
        missing_genes = [g for g in genes_of_interest if g not in gene_callers_ids]
        if len(missing_genes):
            raise ConfigError("The following genes are not in the gene_level_coverage_stats_dict, and yet you are asking for them: %s" % missing_genes)
    else:
        genes_of_interest = gene_callers_ids

    if samples_of_interest is not None:
        missing_samples = [s for s in samples_of_interest if s not in samples]
        if len(missing_samples):
            raise ConfigError("The following samples are not in the gene_level_coverage_stats_dict, and yet you are asking for them: %s" % missing_samples)
    else:
        samples_of_interest = samples

    d = {}

    for gene_callers_id in genes_of_interest:
        d[gene_callers_id] = {}
        for sample_name in samples_of_interest:
            d[gene_callers_id][sample_name] = gene_level_coverage_stats_dict[gene_callers_id][sample_name][key]

    if as_pandas:
        # This option is used by the mcg-classifier.
        import pandas as pd
        return pd.DataFrame.from_dict(d, orient='index')
    else:
        return d

