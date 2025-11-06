import pandas as pd
import numpy as np

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError

from anvio.metabolism.context import KeggContext


class KeggModuleEnrichment(KeggContext):
    """This class is a driver for anvi-script-enrichment-stats for modules input.

    It takes in the modules mode output from anvi-estimate-metabolism, formats it for the enrichment script,
    and runs the script.

    ==========
    args: Namespace object
        All the arguments supplied by user to anvi-compute-functional-enrichment
    """

    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.modules_txt = A('modules_txt')
        self.groups_txt = A('groups_txt')
        self.sample_header_in_modules_txt = A('sample_header') or 'db_name'
        self.module_completion_threshold = A('module_completion_threshold') or 0.75
        self.output_file_path = A('output_file')
        self.include_missing = True if A('include_samples_missing_from_groups_txt') else False
        self.use_stepwise_completeness = A('use_stepwise_completeness')
        self.qlambda = A('qlambda')

        # init the base class
        KeggContext.__init__(self, self.args)

        # if necessary, assign 0 completion threshold, which evaluates to False above
        if A('module_completion_threshold') == 0:
            self.module_completion_threshold = 0.0
            if not self.quiet:
                self.run.warning("Your completion threshold is set to 0, which will make the enrichment results MEANINGLESS. Why? Because "
                                 "with a threshold this low, every module will be considered present in every single sample and therefore will "
                                 "be equally present in every single group. So you should stop what you are doing RIGHT NOW.")
            if not self.just_do_it:
                raise ConfigError("We are stopping you right there because your completion threshold is 0 and that will make the enrichment results "
                                  "meaningless (see the warnings above, if you haven't suppressed them with --quiet). But if you really really really "
                                  "want to do it, you can run again with --just-do-it and then we won't stop you. You have the right to be without meaning.")

        # sanity checkses my precious
        if not self.modules_txt:
            raise ConfigError("To compute module enrichment, you must provide a modules-txt file (aka modules mode output from "
                              "`anvi-estimate-metabolism`).")
        if not self.groups_txt:
            raise ConfigError("To compute module enrichment, you must provide a groups-txt file mapping each sample to a group.")

        filesnpaths.is_file_exists(self.modules_txt)
        filesnpaths.is_file_plain_text(self.modules_txt)
        filesnpaths.is_file_exists(self.groups_txt)
        filesnpaths.is_file_plain_text(self.groups_txt)

        # make sure we are not overwriting anything
        filesnpaths.is_output_file_writable(self.output_file_path, ok_if_exists=False)


        if not self.quiet:
            self.run.info("modules-txt input file", self.modules_txt)
            self.run.info("groups-txt input file", self.groups_txt)
            self.run.info("sample column in modules-txt", self.sample_header_in_modules_txt)
            self.run.info("module completion threshold", self.module_completion_threshold)


    def get_enrichment_input(self, output_file_path):
        """This function converts modules mode output into input for anvi-script-enrichment-stats

        The input format for anvi-script-enrichment-stats is described in a comment at the top of that script, and here is
        how we get the values for each column:
        The first column, 'KEGG_MODULE', and second column 'accession', are already in the modules mode output as 'module_name'
        and 'module', respectively.
        The 'N_*' columns are the total number of samples in each group.
        For each module, this function determines which samples the module is 'present' in according to the specified completion threshold.
        This determines the list of samples for the 'sample_ids' column as well as the 'p_*' proportions for each group of samples.
        Finally, the fourth column, 'associated_groups', is computed from the 'p_*' proportions and 'N_*' totals.

        PARAMETERS
        ==========
        output_file_path : str
            a file path where we will store the (temporary) input file for the enrichment script
        """

        filesnpaths.is_output_file_writable(output_file_path)

        # read the files into dataframes
        modules_df = pd.read_csv(self.modules_txt, sep='\t')

        completeness_header = 'pathwise_module_completeness'
        if self.use_stepwise_completeness:
            completeness_header = 'stepwise_module_completeness'
        self.progress.reset()
        self.run.info("Completeness score being used for determining sample presence", completeness_header)

        # make sure we have all the columns we need in modules mode output, since this output can be customized
        required_modules_txt_headers = ['module', 'module_name', completeness_header]
        missing_headers = []
        for h in required_modules_txt_headers:
            if h not in modules_df.columns:
                missing_headers.append(h)
        if missing_headers:
            missing_string = ", ".join(missing_headers)
            self.progress.reset()
            raise ConfigError("We cannot go on! *dramatic sweep*   We trust that you have provided us with "
                              "modules mode output, but unfortunately the modules-txt input does not contain "
                              f"the following required headers: {missing_string}   Please re-generate your "
                              "modules-txt to include these before trying again.")

        if 'unique_id' in modules_df.columns:
            modules_df = modules_df.drop(columns=['unique_id'])

        # samples column sanity check - this column will become the index
        if self.sample_header_in_modules_txt not in modules_df.columns:
            col_list = ", ".join(modules_df.columns)
            self.progress.reset()
            raise ConfigError(f"You have specified that your sample names are in the column with header '{self.sample_header_in_modules_txt}' "
                               "in the modules-txt file, but that column does not exist. :( Please figure out which column is right and submit "
                               "it using the --sample-header parameter. Just so you know, the columns in modules-txt that you can choose from "
                               f"are: {col_list}")

        samples_to_groups_dict, groups_to_samples_dict = utils.get_groups_txt_file_as_dict(self.groups_txt, include_missing_samples_is_true=self.include_missing)

        # make sure the samples all have a group
        samples_with_none_group = []
        for s,g in samples_to_groups_dict.items():
            if not g:
                samples_with_none_group.append(s)

        for s in samples_with_none_group:
            samples_to_groups_dict.pop(s)

        if samples_with_none_group:
            self.progress.reset()
            none_group_str = ", ".join(samples_with_none_group)
            self.run.warning("Some samples in your groups-txt did not have a group, and we will ignore those samples. If you "
                                 "want them to be included in the analysis, you need to fix the groups-txt to have a group for "
                                 "these samples. Anyway. Here are the samples we will be ignoring: "
                                 f"{none_group_str}")

        # sanity check for mismatch between modules-txt and groups-txt
        sample_names_in_modules_txt = set(modules_df[self.sample_header_in_modules_txt].unique())
        sample_names_in_groups_txt = set(samples_to_groups_dict.keys())
        samples_missing_in_groups_txt = sample_names_in_modules_txt.difference(sample_names_in_groups_txt)
        samples_missing_in_modules_txt = sample_names_in_groups_txt.difference(sample_names_in_modules_txt)
        if anvio.DEBUG:
            self.run.info("Samples in modules-txt", ", ".join(list(sample_names_in_modules_txt)))
            self.run.info("Samples in groups-txt", ", ".join(list(sample_names_in_groups_txt)))
            self.run.info("Missing samples from groups-txt", ", ".join(list(samples_missing_in_groups_txt)))
            self.run.info("Missing samples from modules-txt", ", ".join(list(samples_missing_in_modules_txt)))

        if samples_missing_in_groups_txt:
            missing_samples_str = ", ".join(samples_missing_in_groups_txt)
            if not self.include_missing:
                self.progress.reset()
                self.run.warning(f"Your groups-txt file does not contain some samples present in your modules-txt ({self.sample_header_in_modules_txt} "
                                "column). Since you have not elected to --include-samples-missing-from-groups-txt, we are not going to take these samples into consideration at all. "
                                "Here are the samples that we will be ignoring: "
                                f"{missing_samples_str}")
                # drop the samples that are not in groups-txt
                modules_df = modules_df[~modules_df[self.sample_header_in_modules_txt].isin(list(samples_missing_in_groups_txt))]
                if anvio.DEBUG:
                    self.run.info("Samples remaining in modules-txt dataframe after removing ungrouped", ", ".join(modules_df[self.sample_header_in_modules_txt].unique()))

            else:
                self.progress.reset()
                self.run.warning(f"Your groups-txt file does not contain some samples present in your modules-txt ({self.sample_header_in_modules_txt} "
                                "column). Since you have chosen to --include-samples-missing-from-groups-txt, for the purposes of this analysis we will now consider all of "
                                "these samples to belong to one group called 'UNGROUPED'."
                                "Here are the {len(samples_missing_in_groups_txt)} UNGROUPED samples that we will consider as one big happy family: "
                                f"{missing_samples_str}")
                # add those samples to the UNGROUPED group
                ungrouped_samples = list(samples_missing_in_groups_txt)
                for s in ungrouped_samples:
                    samples_to_groups_dict[s] = 'UNGROUPED'

        if samples_missing_in_modules_txt:
            missing_samples_str = ", ".join(samples_missing_in_modules_txt)
            if not self.just_do_it:
                self.progress.reset()
                raise ConfigError(f"Your modules-txt file ({self.sample_header_in_modules_txt} column) does not contain some samples that "
                                 "are present in your groups-txt. This is not necessarily a huge deal, it's just that those samples will "
                                 "not be included in the enrichment analysis because, well, you don't have any module information for them. "
                                 "If all of the missing samples belong to groups you don't care about at all, then feel free to ignore this "
                                 "message and re-run using --just-do-it. But if you do care about those groups, you'd better fix this because "
                                 "the enrichment results for those groups will be wrong. Here are the samples in question: "
                                  f"{missing_samples_str}")
            else:
                self.progress.reset()
                self.run.warning(f"Your modules-txt file ({self.sample_header_in_modules_txt} column) does not contain some samples that "
                                 "are present in your groups-txt. This is not necessarily a huge deal, it's just that those samples will "
                                 "not be included in the enrichment analysis because, well, you don't have any module information for them. "
                                 "Since you have used the --just-do-it parameter, we assume you don't care about this and are going to keep "
                                 "going anyway. We hope you know what you are doing :) Here are the samples in question: "
                                  f"{missing_samples_str}")
                # drop the samples that are not in modules-txt
                for s in list(samples_missing_in_modules_txt):
                    samples_to_groups_dict.pop(s)
                if anvio.DEBUG:
                    self.run.info("Samples remaining in groups-txt dataframe after removing ungrouped", ", ".join(samples_to_groups_dict.keys()))


        modules_df.set_index(self.sample_header_in_modules_txt, inplace=True)
        sample_groups_df = pd.DataFrame.from_dict(samples_to_groups_dict, orient="index", columns=['group'])

        # convert modules mode output to enrichment input
        N_values = sample_groups_df['group'].value_counts()
        group_list = N_values.keys()
        module_list = modules_df['module'].unique()

        output_dict = {}
        header_list = ['MODULE', 'accession', 'sample_ids', 'associated_groups']
        for c in group_list:
            header_list.append(f"p_{c}")
            header_list.append(f"N_{c}")

        for mod_num in module_list:
            query_string = f"module == '{mod_num}' and {completeness_header} >= {self.module_completion_threshold}"
            samples_with_mod_df = modules_df.query(query_string)
            if samples_with_mod_df.shape[0] == 0:
                continue
            # if we are working with module data from metagenomes, we may have multiple complete copies of the module in
            # the same sample. We drop these duplicates before proceeding.
            duplicates = samples_with_mod_df.index.duplicated()
            samples_with_mod_df = samples_with_mod_df[~duplicates]

            # we need to explicitly ignore samples without a group here, because they were taken out of sample_groups_df
            # and if only ungrouped samples end up having this module, we will get an index error
            samples_with_mod_list = list(samples_with_mod_df.index)
            for s in samples_with_none_group:
                if s in samples_with_mod_list:
                    samples_with_mod_list.remove(s)
            if len(samples_with_mod_list) == 0:
                continue

            mod_name = samples_with_mod_df['module_name'][0]
            output_dict[mod_name] = {}
            output_dict[mod_name]['MODULE'] = mod_name
            output_dict[mod_name]['accession'] = mod_num
            output_dict[mod_name]['sample_ids'] = ','.join(samples_with_mod_list)
            sample_group_subset = sample_groups_df.loc[samples_with_mod_list]
            p_values = sample_group_subset['group'].value_counts()

            # we need the categories p and N values to be in the same order for finding associated groups
            p_vector = np.array([])
            N_vector = np.array([])
            for c in group_list:
                if c not in p_values.index:
                    p_values[c] = 0
                p_vector = np.append(p_vector, p_values[c]/N_values[c])
                N_vector = np.append(N_vector, N_values[c])

            # compute associated groups for functional enrichment
            enriched_groups_vector = utils.get_enriched_groups(p_vector, N_vector)

            associated_groups = [c for i,c in enumerate(group_list) if enriched_groups_vector[i]]
            output_dict[mod_name]['associated_groups'] = ','.join(associated_groups)

            for c in group_list:
                output_dict[mod_name]["p_%s" % c] = p_values[c]/N_values[c]
                output_dict[mod_name]["N_%s" % c] = N_values[c]

        utils.store_dict_as_TAB_delimited_file(output_dict, output_file_path, key_header='accession', headers=header_list)


    def run_enrichment_stats(self):
        """This function is the driver for running the enrichment script on the modules data."""

        self.progress.new('Enrichment analysis')

        self.progress.update('Converting modules mode output into input for enrichment script')
        enrichment_input_path = filesnpaths.get_temp_file_path()

        if anvio.DEBUG:
            self.progress.reset()
            self.run.info("Temporary input file for enrichment script", enrichment_input_path)

        self.get_enrichment_input(enrichment_input_path)

        self.progress.end()

        # run the enrichment analysis
        enrichment_stats = utils.run_functional_enrichment_stats(enrichment_input_path,
                                                                 self.output_file_path,
                                                                 qlambda=self.qlambda,
                                                                 run=self.run,
                                                                 progress=self.progress)

        return enrichment_stats
