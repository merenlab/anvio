#!/usr/bin/env python
# -*- coding: utf-8
"""Core estimation algorithms for the metabolism module. NO BD ACCESS FROM HERE."""

import re
import statistics

from scipy import stats

import anvio
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.errors import ConfigError
from anvio.metabolism.constants import STRAY_KO_ANVIO_SUFFIX


class KeggEstimationAlgorithms:
    """Core estimation algorithms for KEGG metabolism."""

    def __init__(self, run=terminal.Run(), progress=terminal.Progress()):
        self.run = run
        self.progress = progress


    def split_module_path_into_individual_essential_components(self, path):
        """Given a list of atomic steps in a module, this function returns a list of each essential individual enzyme.

        When an atomic step is an enzyme complex (ie K01657+K01658), we need to split the complex into its individual components
        When an atomic step contains non-ssential components (ie K00765-K02502), we need to remove the nonessential components from the list
        When there are both nonessential and essential components, we need to remove the non-essential ones first and then split the essential ones

        PARAMETERS
        ==========
        path : list
            Each element in list is an atomic step, which can include enzyme complexes and non-essential components

        RETURNS
        ==========
        new_path : list
            Each element in list is a single enzyme
        """

        new_path = []
        for x in path:
            if '+' and '-' in x:
                # first remove the nonessentials, then add in the essential components
                idx_of_nonessential = x.index('-')
                new_x = x[:idx_of_nonessential]
                individual_components = new_x.split('+')
                new_path.extend(individual_components)
            elif '+' in x:
                # split essential components and add individually to list
                individual_components = x.split('+')
                new_path.extend(individual_components)
            elif '-' in x and x != '--':
                # remove nonessentials
                idx_of_nonessential = x.index('-')
                new_x = x[:idx_of_nonessential]
                new_path.append(new_x)
            else:
                new_path.append(x)

        return new_path


    def get_enzymes_from_module_definition_in_order(self, mod_definition, all_modules_in_db):
        """Given a module DEFINITION string, this function parses out the enzyme accessions in order of appearance.

        PARAMETERS
        ==========
        mod_definition : a string or list of strings containing the module DEFINITION lines

        RETURNS
        ==========

        """

        if isinstance(mod_definition, list):
            mod_definition = " ".join(mod_definition)

        acc_list = module_definition_to_enzyme_accessions(mod_definition)

        # remove anything that is not an enzyme and sanity check for weird characters
        mods_to_remove = set()
        for a in acc_list:
            if a in all_modules_in_db:
                mods_to_remove.add(a)
            if re.match('[^a-zA-Z0-9_\.]', a):
                raise ConfigError(f"The get_enzymes_from_module_definition_in_order() function found an enzyme accession that looks a bit funny. "
                                  f"Possibly this is a failure of our parsing strategy, or maybe the enzyme accession just has unexpected characters "
                                  f"in it. We don't know what module it is, but the weird enzyme is {a}. If you think that accession looks perfectly "
                                  f"fine, you should reach out to the developers and have them fix this function to accomodate the accession. Or, you "
                                  f"could just rename the enzyme?")
        if mods_to_remove:
            for m in mods_to_remove:
                acc_list.remove(m)

        return acc_list


    def remove_nonessential_enzymes_from_module_step(self, step_string, exclude_dashed_reactions=True):
        """This functions removes nonessential enzymes from a step definition string and returns the resulting definition.

        It is intended to be called on top-level steps of a module definition (not on the full module definition).

        A nonessential enzyme is any accession (ie, '-K11024') or group of accessions (ie, -(K00242,K18859,K18860))
        that is marked with a leading '-' character. This function finds the '-' characters in the given string and
        removes any subsequent accessions. The resulting definition with only essential components is returned.

        If a step does not contain nonessential enzymes, the original definition is returned.

        Likewise, '--' is a special case containing the '-' character which actually indicates a step that has no enzyme profile.
        It seems to always be present as its own step (ie, '--' is the entire step definition string), so we return the original
        definition in this case. However, in case there is an internal '--' within a more complicated definition, this function
        ignores the part of the string that includes it and processes the remainder of the string before re-joining the two parts.
        It is not able to do this for steps with more than one internal '--', which would require multiple splits and joins, so
        this case results in an error. Note that when self.exclude_dashed_reactions is True, we instead remove '--' entirely.

        PARAMETERS
        ==========
        step_string: str
            A string containing the definition of one step from a module

        RETURNS
        =======
        step_string: str
            The same string, with nonessential enzyme accessions (if any) removed.
        """

        if step_string == '--' and exclude_dashed_reactions:
            return ""
        if step_string != '--' and '-' in step_string:
            saw_double_dash = False             # a Boolean to indicate if we found '--' within the step definition
            str_prior_to_double_dash = None     # if we find '--', this variable stores the string that occurs prior to and including this '--'
            while '-' in step_string:
                idx = step_string.index('-')
                if step_string[idx+1] == '-': # internal '--' case
                    if saw_double_dash:
                        raise ConfigError("Unfortunately, a step containing >1 internal instances of '--' was passed to the function "
                                          "remove_nonessential_enzymes_from_module_step(). This function is not currently able to handle this "
                                          "situation. Please contact a developer and ask them to turn this into a smarter function. :) ")
                    saw_double_dash = True
                    if exclude_dashed_reactions: # remove the internal '--'
                        str_prior_to_double_dash = step_string[:idx]
                        step_string = step_string[idx+3:] # also remove the space after it
                    else:
                        str_prior_to_double_dash = step_string[:idx+2]
                        step_string = step_string[idx+2:] # continue processing the remainder of the string
                    continue
                elif step_string[idx+1] == '(': # group to eliminate
                    parens_index = idx+1
                    while step_string[parens_index] != ')':
                        parens_index += 1
                    step_string = step_string[:idx] + step_string[parens_index+1:]
                else: # single non-essential enzyme
                    punctuation_index = idx+1
                    while punctuation_index < len(step_string) and step_string[punctuation_index] not in [')','(','+',' ',',']:
                        punctuation_index += 1
                    step_string = step_string[:idx] + step_string[punctuation_index:]

            # if we found an internal '--', we put the two pieces of the definition back together
            if saw_double_dash:
                step_string = str_prior_to_double_dash + step_string

        return step_string.strip()


    def mark_kos_present_for_list_of_splits(self, kofam_hits_in_splits, all_modules_in_db, all_kos_in_db,
                                          ko_dict, stray_ko_dict=None, include_stray_kos=False,
                                          exclude_kos_no_threshold=True, ignore_unknown_kos=False,
                                          split_list=None, bin_name=None):
        """This function generates two bin-level dictionaries of dictionaries to store metabolism data.

        The first dictionary of dictionaries is the module completeness dictionary, which associates modules with the KOs
        that are present in the bin for each module.

        The structure of the dictionary is like this example:
        {mnum: {"gene_caller_ids" : set([132, 133, 431, 6777]),
                "kofam_hits" : {'K00033' : [431, 6777],
                                'K01057' : [133],
                                'K00036' : [132] },
                "genes_to_contigs": {132: 0,
                                     133: 0,
                                     431: 2,
                                    6777: 1 },
                "contigs_to_genes": { 0: set([132, 133]),
                                      1: set(6777),
                                      2: set(431) },}}
        This dictionary will be expanded later by other functions.

        The second dictionary of dictionaries is the KOfam hit dictionary, which stores all of the KOfam hits in the bin
        regardless of whether they are in a KEGG module or not.

        The structure of the dictionary is like this example:
        {ko: {"gene_caller_ids" : set([431, 6777]),
              "modules" : ["M00001", "M00555"],                 **Can be None if KO does not belong to any KEGG modules
              "genes_to_contigs": { 431: 2,
                                   6777: 1 },
              "contigs_to_genes": { 1: set(6777),
                                    2: set(431) }}}


        PARAMETERS
        ==========
        kofam_hits_in_splits : list
            (ko_num, gene_call_id, split, contig) tuples, one per KOfam hit in the splits we are considering
        split_list : list
            splits we are considering, this is only for debugging output
        bin_name : str
            name of the bin containing these splits, this is only for debugging output

        RETURNS
        =======
        bin_level_module_dict : dictionary of dictionaries
            initialized module completeness dictionary for the list of splits (genome, metagenome, or bin) provided
        bin_level_ko_dict : dictionary of dictionaries
            dictionary of ko hits within the list of splits provided
        """

        bin_level_module_dict = {}
        bin_level_ko_dict = {}

        if anvio.DEBUG:
            self.run.info("Marking KOs present for bin", bin_name)
            if split_list:
                num_splits = len(split_list)
            else:
                num_splits = "None"
            self.run.info("Number of splits", num_splits)

        # initialize all modules with empty lists and dicts for kos, gene calls
        modules = all_modules_in_db.keys()
        all_kos = all_kos_in_db.keys()
        for mnum in modules:
            bin_level_module_dict[mnum] = {"gene_caller_ids" : set(),
                                           "kofam_hits" : {},
                                           "genes_to_contigs" : {},
                                           "contigs_to_genes" : {},
                                           "unique_to_this_module": set(),
                                           "warnings" : set()
                                          }

        for knum in all_kos:
            """
            We can only add warnings about missing KOfam profiles because for other annotation sources, we don't
            have a way to know if profiles are missing. But for KOfams with missing profiles, this step is necessary
            so that we don't add the enzyme to the bin_level_ko_dict, because later this will cause problems since
            the enzyme is not in self.ko_dict

            Furthermore, this can only be done when we are using both KEGG data and user data (ie, not --only-user-modules)
            because we need access to the self.ko_dict
            """

            if knum.startswith("M"):
                continue

            if (ko_dict is not None
                and all_kos_in_db[knum]['annotation_source'] == 'KOfam'
                and knum not in ko_dict
                and (not stray_ko_dict or (knum not in stray_ko_dict and f"{knum}{STRAY_KO_ANVIO_SUFFIX}" not in stray_ko_dict))
                and exclude_kos_no_threshold):

                mods_it_is_in = all_kos_in_db[knum]['modules']
                if mods_it_is_in:
                    if anvio.DEBUG:
                        mods_str = ", ".join(mods_it_is_in)
                        self.run.warning(f"Oh dear. We do not appear to have a KOfam profile for {knum}. This means "
                                         f"that any modules this KO belongs to can never be fully complete (this includes "
                                         f"{mods_str}). ")
                    for m in mods_it_is_in:
                        bin_level_module_dict[m]["warnings"].add(f"No KOfam profile for {knum}")

                if anvio.DEBUG:
                    if exclude_kos_no_threshold:
                        self.run.warning(f"We cannot find an entry for KO {knum} in the `ko_list.txt` file downloaded "
                                         f"from KEGG. What this means is that you are somehow using KOfam annotations "
                                         f"that are different from the current version of KOfam on your computer (this can "
                                         f"happen with --enzymes-txt input). Because we are not considering these annotations, "
                                         f"you may get KeyErrors downstream. You can force the inclusion of these KOfams by "
                                         f"re-running this program with the --include-kos-not-in-kofam flag.")

                continue

            bin_level_ko_dict[knum] = {"gene_caller_ids" : set(),
                                     "modules" : None,
                                     "genes_to_contigs" : {},
                                     "contigs_to_genes" : {},
                                     "warnings" : set()
                                     }

        kos_not_in_modules = []

        for ko, gene_call_id, split, contig in kofam_hits_in_splits:
            # make sure we can count annotations to anvi'o versions of stray KO models by using KEGG's original accession
            is_anvio_version = False
            if ko.endswith(STRAY_KO_ANVIO_SUFFIX):
                is_anvio_version = True
                ko = ko.replace(STRAY_KO_ANVIO_SUFFIX, "")

            if ko not in all_kos_in_db:

                kos_not_in_modules.append(ko)

                # KOs that are not in modules will not be initialized above in the ko hit dictionary, so we add them here if we haven't already
                if ko not in bin_level_ko_dict:
                    bin_level_ko_dict[ko] = {"gene_caller_ids" : set(),
                                             "modules" : None,
                                             "genes_to_contigs" : {},
                                             "contigs_to_genes" : {},
                                             "warnings" : set()
                                             }
            else:

                # if we are missing the KO from the dictionary at this point, we should fail nicely instead of with a KeyError
                if ko not in bin_level_ko_dict:
                    if ignore_unknown_kos:
                        continue
                    raise ConfigError(f"We cannot find the KEGG enzyme {ko} in the dictionary of enzyme hits, even though this enzyme is "
                                      f"annotated in your data. There are 3 main ways this can happen: (1) you are using --enzymes-txt input "
                                      f"that includes KOs that are different from the set used for annotation with `anvi-run-kegg-kofams`, "
                                      f"(2) your contigs database was annotated with `anvi-run-kegg-kofams --include-stray-KOs` but you didn't "
                                      f"use the `--include-stray-KOs` flag for `anvi-estimate-metabolism`, or (3) you imported external annotations "
                                      f"with the source name `KOfam` that include KOs not in the KEGG data directory currently being used. "
                                      f"You have a few options to get around this error depending on which case applies to your situation. "
                                      f"If this is case (2) and you want to include these enzymes in the analysis, then re-run `anvi-estimate-metabolism` "
                                      f"with the flag `--include-stray-KOs`. If this is case (1) or (3) and you want to include these enzymes in the "
                                      f"analysis, then re-run `anvi-estimate-metabolism` with the flag `--include-kos-not-in-kofam`. And no matter "
                                      f"what the situation is, if you want to IGNORE these unknown annotations for the purposes of estimating "
                                      f"metabolism, you can re-run `anvi-estimate-metablism` with the flag `--ignore-unknown-KOs`. If this message "
                                     f"made you worry, you could also re-do your annotations or remove these unknown enzymes from your input "
                                     f"--enzymes-txt file to be on the safe side.")

                present_in_mods = all_kos_in_db[ko]['modules']
                bin_level_ko_dict[ko]["modules"] = present_in_mods

                # keep track of enzymes unique to this module
                is_unique = False
                if len(present_in_mods) == 1:
                    is_unique = True

                for m in present_in_mods:
                    bin_level_module_dict[m]["gene_caller_ids"].add(gene_call_id)
                    if ko in bin_level_module_dict[m]["kofam_hits"] and gene_call_id not in bin_level_module_dict[m]["kofam_hits"][ko]:
                        bin_level_module_dict[m]["kofam_hits"][ko].append(gene_call_id)
                    else:
                        bin_level_module_dict[m]["kofam_hits"][ko] = [gene_call_id]
                    bin_level_module_dict[m]["genes_to_contigs"][gene_call_id] = contig
                    if contig in bin_level_module_dict[m]["contigs_to_genes"]:
                        bin_level_module_dict[m]["contigs_to_genes"][contig].add(gene_call_id)
                    else:
                        bin_level_module_dict[m]["contigs_to_genes"][contig] = set([gene_call_id])

                    # make a special list for the enzymes that are unique
                    if is_unique:
                        bin_level_module_dict[m]["unique_to_this_module"].add(ko)
                    # warn the user if this enzyme is shared between multiple modules
                    else:
                        mod_str = "/".join(present_in_mods)
                        bin_level_module_dict[m]["warnings"].add(f"{ko} is present in multiple modules: {mod_str}")

                    # point out use of anvi'o-specific KO models
                    if is_anvio_version:
                        warning = f"used '{ko}{STRAY_KO_ANVIO_SUFFIX}' model to annotate {ko}"
                        bin_level_module_dict[m]["warnings"].add(warning)
                        bin_level_ko_dict[ko]["warnings"].add(warning)


            bin_level_ko_dict[ko]["gene_caller_ids"].add(gene_call_id)
            bin_level_ko_dict[ko]["genes_to_contigs"][gene_call_id] = contig
            if contig in bin_level_ko_dict[ko]["contigs_to_genes"]:
                bin_level_ko_dict[ko]["contigs_to_genes"][contig].add(gene_call_id)
            else:
                bin_level_ko_dict[ko]["contigs_to_genes"][contig] = set([gene_call_id])

        if anvio.DEBUG:
            self.run.info("Gene calls processed", "%d in bin" % len(kofam_hits_in_splits))
            if kos_not_in_modules:
                self.run.warning(f"Just so you know, the following enzymes did not belong to any modules in the MODULES.db: {', '.join(kos_not_in_modules)}")

        return bin_level_module_dict, bin_level_ko_dict


    def compute_stepwise_module_completeness_for_bin(self, mnum, meta_dict_for_bin, all_modules_in_db, module_completion_threshold, exclude_dashed_reactions=True):
        """This function calculates the stepwise completeness of the specified module within the given bin dictionary.

        It uses only the "top-level" steps of the module definition, which are the steps that you get when you first
        split the module definition on a space. Each "top-level" step is comprised of one or more enzymes that either
        work together or serve as alternatives to each other. In this calculation, we ignore the possible combinations
        of enzymes and simply decide whether or not a "top-level" step is complete or not. Then the module completeness
        is computed as the number of complete "top-level" steps divided by the total number of steps.

        PARAMETERS
        ==========
        mnum : string
            module number to work on
        meta_dict_for_bin : dictionary of dictionaries
            metabolism completeness dict for the current bin, to be modified in-place

        NEW KEYS ADDED TO METABOLISM COMPLETENESS DICT
        =======
        "stepwise_completeness"         the stepwise completeness of the module
        "stepwise_is_complete"          whether the module completeness falls over the completeness threshold
        "top_level_step_info"           a dictionary of each top level step
                                            keyed by integer from 0 to # of top level steps.
                                            inner dict contains the following keys:
                                            'step_definition' (string)
                                            'complete' (Boolean)
                                            'includes_modules' (Boolean)
                                            'included_module_list' (list of strings)

        RETURNS
        =======
        over_complete_threshold : boolean
            whether or not the module is considered "complete" overall based on the threshold fraction of completeness
        """

        top_level_steps = all_modules_in_db[mnum]['top_level_steps']
        num_steps = len(top_level_steps)
        num_complete = 0
        num_nonessential_steps = 0

        present_list_for_mnum = meta_dict_for_bin[mnum]["kofam_hits"].keys()

        meta_dict_for_bin[mnum]['top_level_step_info'] = {}

        for i, step in enumerate(top_level_steps):
            step_is_present_condition_statement = ""
            cur_index = 0  # current position in the step DEFINITION
            step_includes_modules = False
            included_module_list = []

            while cur_index < len(step):
                if step[cur_index] in ['(',')']:
                    step_is_present_condition_statement += step[cur_index]
                    cur_index += 1

                elif step[cur_index] == ",":
                    step_is_present_condition_statement += " or "
                    cur_index += 1

                elif step[cur_index] == "+" or step[cur_index] == ' ':
                    step_is_present_condition_statement += " and "
                    cur_index += 1

                elif step[cur_index] == "-":
                    # '--' no associated enzyme case, by default False (assumed incomplete)
                    if step[cur_index+1] == "-":
                        if exclude_dashed_reactions: # skip it instead
                            cur_index += 3 # skip over both '-' AND the following space
                        else:
                            step_is_present_condition_statement += "False"
                            cur_index += 2 # skip over both '-', the next character should be a space or end of DEFINITION line

                        if anvio.DEBUG:
                            self.run.warning(f"While estimating the stepwise completeness of KEGG module {mnum}, anvi'o saw "
                                             f"'--' in the module DEFINITION. This indicates a step in the pathway that has no "
                                             f"associated enzyme. By default, anvi'o marks steps like these incomplete, *unless* "
                                             f"you are using the flag --exclude-dashed-reactions. But if you aren't using that flag, "
                                             f"it is possible that this module might be falsely considered incomplete. So it may be in your "
                                             f"interest to take a closer look at module {mnum}.")
                        if cur_index < len(step) and step[cur_index] != " ":
                            raise ConfigError(f"Serious, serious parsing sadness is happening. We just processed a '--' in "
                                              f"a DEFINITION line for module {mnum} but did not see a space afterwards. Instead, "
                                              f"we found {step[cur_index+1]}. WHAT DO WE DO NOW?")

                    # a whole set of nonessential KOs - skip all of them
                    elif step[cur_index+1] == "(":
                        while step[cur_index] != ")":
                            cur_index += 1
                        cur_index += 1 # skip over the ')'

                    # anything else that follows a '-' should be an enzyme or enzyme component, and should be skipped
                    else:
                        # find the next space or '-' or the end of the step
                        while cur_index+1 < len(step) and (step[cur_index+1] not in [' ', ',', '+', '-', '(', ')']):
                            cur_index += 1
                        cur_index += 1
                        # if we found a non-accession character, the next iteration of the loop will take care of it
                        # if we reached the end of the step and the condition statement is empty, then the entire
                        #    step is nonessential so we need to avoid counting it (taken care of later)

                else: # enzyme or module accession
                    enzyme_start_index = cur_index
                    while cur_index+1 < len(step) and step[cur_index+1] not in [' ', ',', '+', '-', '(', ')']:
                        cur_index += 1
                    enzyme_end_index = cur_index

                    accession = step[enzyme_start_index : enzyme_end_index+1]
                    # module
                    if accession in all_modules_in_db:
                        step_includes_modules = True
                        included_module_list.append(accession)
                        # store the module accession in the condition string to be replaced later
                        step_is_present_condition_statement += accession
                    # enzyme
                    elif accession in present_list_for_mnum:
                        step_is_present_condition_statement += "True"
                    else:
                        step_is_present_condition_statement += "False"

                    cur_index += 1

            meta_dict_for_bin[mnum]['top_level_step_info'][i] = {}
            meta_dict_for_bin[mnum]['top_level_step_info'][i]['step_definition'] = step
            meta_dict_for_bin[mnum]['top_level_step_info'][i]['includes_modules'] = step_includes_modules
            meta_dict_for_bin[mnum]['top_level_step_info'][i]['included_module_list'] = included_module_list

            # entire step was nonessential, do not count it
            if step_is_present_condition_statement == "":
                num_nonessential_steps += 1
                meta_dict_for_bin[mnum]['top_level_step_info'][i]['complete'] = "nonessential"
            elif step_includes_modules:
                # we'll eval this condition statement in a later function once all other modules have stepwise completeness
                meta_dict_for_bin[mnum]['top_level_step_info'][i]['complete'] = step_is_present_condition_statement
            else:
                step_is_present = eval(step_is_present_condition_statement)
                meta_dict_for_bin[mnum]['top_level_step_info'][i]['complete'] = step_is_present
                if step_is_present:
                    num_complete += 1

        # compute stepwise completeness as number of complete (essential) steps / number of total (essential) steps
        mod_stepwise_completeness = num_complete / (num_steps - num_nonessential_steps)
        meta_dict_for_bin[mnum]["stepwise_completeness"] = mod_stepwise_completeness

        over_complete_threshold = True if meta_dict_for_bin[mnum]["stepwise_completeness"] >= module_completion_threshold else False
        meta_dict_for_bin[mnum]["stepwise_is_complete"] = over_complete_threshold

        return over_complete_threshold


    def compute_pathwise_module_completeness_for_bin(self, mnum, meta_dict_for_bin, module_paths_dict, module_completion_threshold,
                                                     exclude_dashed_reactions=True, all_modules_in_db=None):
        """This function calculates the pathwise completeness of the specified module within the given bin metabolism dictionary.

        To do this, it works with the unrolled module definition: a list of all possible paths, where each path is a list of atomic steps.
        Atomic steps include singular KOs, protein complexes, modules, non-essential steps, and steps without associated KOs.
        An atomic step (or parts of a protein complex) can be considered 'present' if the corresponding KO(s) has a hit in the bin.
        For each path, the function computes the path completeness as the number of present (essential) steps divided by the number of total steps in the path.
        The module completeness is simply the highest path completeness.

        There are some special cases to consider here.
        1) Non-essential steps. These are steps that are marked with a preceding "-" to indicate that they are not required for the module to
           be considered complete. They often occur in pathways with multiple forks. What we do with these is save and count them separately as
           non-essential steps, but we do not use them in our module completeness calculations. Another thing we do is continue parsing the rest
           of the module steps as normal, even though some of them may affect steps after the non-essential one. That may eventually change.
           See comments in the code below.
        2) Steps without associated KOs. These are steps marked as "--". They may require an enzyme, but if so that enzyme is not in the KOfam
           database, so we can't know whether they are complete or not from our KOfam hits. Therefore, we assume these steps are incomplete, and
           warn the user to go back and check the module manually.
        3) Steps defined by entire modules. These steps have module numbers instead of KOs, so they require an entire module to be complete in
           order to be complete. We can't figure this out until after we've evaluated all modules, so we simply parse these steps without marking
           them complete, and later will go back to adjust the completeness score once all modules have been marked complete or not.


        PARAMETERS
        ==========
        mnum : string
            module number to work on
        meta_dict_for_bin : dictionary of dictionaries
            metabolism completeness dict for the current bin, to be modified in-place

        NEW KEYS ADDED TO METABOLISM COMPLETENESS DICT
        =======
        "pathway_completeness"          a list of the completeness of each pathway
        "present_nonessential_kos"      a list of non-essential KOs in the module that were found to be present
        "most_complete_paths"           a list of the paths with maximum completeness
        "pathwise_percent_complete"     the completeness of the module, which is the maximum pathway completeness
        "pathwise_is_complete"          whether the module completeness falls over the completeness threshold

        RETURNS
        =======
        over_complete_threshold : boolean
            whether or not the module is considered "complete" overall based on the threshold fraction of completeness
        has_nonessential_step : boolean
            whether or not the module contains non-essential steps. Used for warning the user about these.
        has_no_ko_step : boolean
            whether or not the module contains steps without associated KOs. Used for warning the user about these.
        defined_by_modules : boolean
            whether or not the module contains steps defined by other modules. Used for going back to adjust completeness later.
        """

        present_list_for_mnum = meta_dict_for_bin[mnum]["kofam_hits"].keys()
        if not present_list_for_mnum:
            # no KOs in this module are present
            if anvio.DEBUG:
                self.run.warning(f"No KOs present for module {mnum}. Parsing for completeness is still being done to obtain module information.")

        # stuff to put in the module's dictionary
        module_nonessential_kos = [] # KOs that are present but unnecessary for module completeness

        # stuff that will be returned
        over_complete_threshold = False
        has_nonessential_step = False
        has_no_ko_step = False
        defined_by_modules = False

        meta_dict_for_bin[mnum]["pathway_completeness"] = []
        meta_dict_for_bin[mnum]["num_complete_copies_of_all_paths"] = []
        meta_dict_for_bin[mnum]["num_complete_copies_of_most_complete_paths"] = []

        for p in module_paths_dict[mnum]:
            num_complete_steps_in_path = 0
            num_nonessential_steps_in_path = 0 # so that we don't count nonessential steps when computing completeness
            atomic_step_copy_number = []

            for atomic_step in p:
                # there are 5 types of atomic steps to take care of
                if any(x in atomic_step for x in ['-','+']):
                    # 1) steps without associated enzymes, ie --
                    if atomic_step == "--":
                        # when '--' in a DEFINITION line happens, it signifies a reaction step that has no associated enzyme.
                        # by default, we assume that such steps are not complete
                        has_no_ko_step = True
                        if exclude_dashed_reactions:
                            warning_str = "'--' step was ignored in the calculation"
                            num_nonessential_steps_in_path += 1 # this is to ensure we fix the denominator later
                        else:
                            warning_str = "'--' steps are assumed incomplete"
                            atomic_step_copy_number.append(0)

                        meta_dict_for_bin[mnum]["warnings"].add(warning_str)
                    # 2) non-essential KOs, ie -Kxxxxx
                    elif atomic_step[0] == "-" and not any(x in atomic_step[1:] for x in ['-','+']):
                        """
                        OKAY, SO HERE WE HAVE SOME POOPINESS THAT MAY NEED TO BE FIXED EVENTUALLY.
                        Basically, some DEFINITION lines have KOs that seem to be marked non-essential;
                        ie, "-K11024" in "K11023 -K11024 K11025 K11026 K11027".
                        It was difficult to decide whether we should consider only K11024, or K11024 and all following KOs, to be non-essential.
                        For instance, the module M00778 is a complex case that gave us pause - see Fiesta issue 955.
                        But for now, we have decided to just track only the one KO as a 'non-essential step', and to not include such steps in
                        the module completeness estimate.
                        """
                        ko = atomic_step[1:]
                        if ko not in module_nonessential_kos:
                            module_nonessential_kos.append(ko)
                        num_nonessential_steps_in_path += 1
                        has_nonessential_step = True

                    # 3) protein complexes, ie Kxxxxx+Kyyyyy-Kzzzzz (2 types of complex components - essential and nonessential)
                    else:
                        # split on '+' or '-'
                        pattern = re.compile('\+|\-')
                        match_idxs = []
                        for match in re.finditer(pattern, atomic_step):
                            match_idxs.append(match.start())

                        essential_components = []
                        num_matches_processed = 0
                        for i, match_idx in enumerate(match_idxs):
                            # if this is the first match, we need to handle the initial component in the complex
                            if num_matches_processed == 0:
                                essential_components.append(atomic_step[0:match_idx])

                            # handle the component after the match character
                            if i < len(match_idxs)-1:
                                next_idx = match_idxs[i+1]
                            else:
                                next_idx = len(atomic_step)
                            component_ko = atomic_step[match_idx+1:next_idx]

                            # essential component after  +
                            if atomic_step[match_idx] == '+':
                                essential_components.append(component_ko)
                            # non-essential component after '-'
                            else:
                                has_nonessential_step = True
                                if component_ko not in module_nonessential_kos:
                                    module_nonessential_kos.append(component_ko)

                            num_matches_processed += 1

                        # after processing all components of the enzyme complex, we compute the complex completeness and copy number
                        num_present_components = 0
                        component_copy_number = []
                        for c in essential_components:
                            if c in present_list_for_mnum:
                                num_present_components += 1
                                num_copies = len(meta_dict_for_bin[mnum]["kofam_hits"][c])
                            else:
                                num_copies = 0
                            component_copy_number.append(num_copies)
                        component_completeness = num_present_components / len(essential_components)
                        num_complete_steps_in_path += component_completeness

                        if component_completeness >= module_completion_threshold:
                            atomic_step_copy_number.append(min(component_copy_number))
                        else:
                            atomic_step_copy_number.append(0)
                else:
                    # atomic step is a single enzyme or module
                    # 4) Module numbers, ie Mxxxxx
                    if all_modules_in_db and atomic_step in all_modules_in_db:
                        """
                        This happens when a module is defined by other modules. For example, photosynthesis module M00611 is defined as
                        (M00161,M00163) M00165 === (photosystem II or photosystem I) and calvin cycle

                        We need all the modules to have been evaluated before we can determine completeness of steps with module numbers.
                        So what we will do here is to use a flag variable to keep track of the modules that have this sort of definition
                        in a list so we can go back and evaluate completeness of steps with module numbers later.
                        """
                        defined_by_modules = True
                    # 5) regular old single enzymes, ie Kxxxxx (for KOs), COGyyyyy (for COGs), etc
                    else:
                        if atomic_step in present_list_for_mnum:
                            num_complete_steps_in_path += 1
                            num_copies = len(meta_dict_for_bin[mnum]["kofam_hits"][atomic_step])
                        else:
                            num_copies = 0

                        atomic_step_copy_number.append(num_copies)

            path_completeness = num_complete_steps_in_path / (len(p) - num_nonessential_steps_in_path)
            meta_dict_for_bin[mnum]["pathway_completeness"].append(path_completeness)

            # compute path copy number
            if defined_by_modules:
                path_copy_number = atomic_step_copy_number # save list with atomic step copy numbers to use when adjusting module copy number later
            else:
                path_copy_number = self.compute_num_complete_copies_of_path(atomic_step_copy_number, module_completion_threshold)
            meta_dict_for_bin[mnum]["num_complete_copies_of_all_paths"].append(path_copy_number)

        # once all paths have been evaluated, we find the path(s) of maximum completeness and set that as the overall module completeness
        # this is not very efficient as it takes two passes over the list but okay
        meta_dict_for_bin[mnum]["pathwise_percent_complete"] = max(meta_dict_for_bin[mnum]["pathway_completeness"])
        if meta_dict_for_bin[mnum]["pathwise_percent_complete"] > 0:
            meta_dict_for_bin[mnum]["most_complete_paths"] = [module_paths_dict[mnum][i] for i, pc in enumerate(meta_dict_for_bin[mnum]["pathway_completeness"]) if pc == meta_dict_for_bin[mnum]["pathwise_percent_complete"]]
            if not defined_by_modules:
                meta_dict_for_bin[mnum]["num_complete_copies_of_most_complete_paths"] = [meta_dict_for_bin[mnum]["num_complete_copies_of_all_paths"][i] for i, pc in enumerate(meta_dict_for_bin[mnum]["pathway_completeness"]) if pc == meta_dict_for_bin[mnum]["pathwise_percent_complete"]]
        else:
            meta_dict_for_bin[mnum]["most_complete_paths"] = []
            meta_dict_for_bin[mnum]["num_complete_copies_of_most_complete_paths"] = []

        # set module copy number as the maximum copy number of the path(s) of maximum completeness
        if meta_dict_for_bin[mnum]["num_complete_copies_of_most_complete_paths"]:
            meta_dict_for_bin[mnum]["pathwise_copy_number"] = max(meta_dict_for_bin[mnum]["num_complete_copies_of_most_complete_paths"])
        else:
            meta_dict_for_bin[mnum]["pathwise_copy_number"] = 'NA'

        # compute proportion of unique enzymes in the module (regardless of which path(s) enzyme is in or whether enzyme is essential)
        if meta_dict_for_bin[mnum]["unique_to_this_module"]:
            num_unique_enzymes_present = 0
            num_unique_enzymes_in_mod = len(meta_dict_for_bin[mnum]["unique_to_this_module"])
            for ko in present_list_for_mnum:
                if ko in meta_dict_for_bin[mnum]["unique_to_this_module"]:
                    num_unique_enzymes_present += 1

            meta_dict_for_bin[mnum]["proportion_unique_enzymes_present"] = num_unique_enzymes_present / num_unique_enzymes_in_mod
            meta_dict_for_bin[mnum]["unique_enzymes_context_string"] = f"{num_unique_enzymes_present} of {num_unique_enzymes_in_mod} unique enzymes in module"
        else:
            meta_dict_for_bin[mnum]["proportion_unique_enzymes_present"] = "NA"
            meta_dict_for_bin[mnum]["unique_enzymes_context_string"] = "NA"


        if anvio.DEBUG and len(meta_dict_for_bin[mnum]["most_complete_paths"]) > 1:
            self.run.warning("Found %d complete paths for module %s with completeness %s. " % (len(meta_dict_for_bin[mnum]["most_complete_paths"]), mnum, meta_dict_for_bin[mnum]["pathwise_percent_complete"]),
                             header='DEBUG OUTPUT', lc='yellow')
        over_complete_threshold = True if meta_dict_for_bin[mnum]["pathwise_percent_complete"] >= module_completion_threshold else False
        meta_dict_for_bin[mnum]["pathwise_is_complete"] = over_complete_threshold
        meta_dict_for_bin[mnum]["present_nonessential_kos"] = module_nonessential_kos

        return over_complete_threshold, has_nonessential_step, has_no_ko_step, defined_by_modules


    def adjust_stepwise_completeness_for_bin(self, mnum, meta_dict_for_bin, module_completion_threshold):
        """This function adjusts stepwise completeness of modules that are defined by other modules.

        This can only be done after all other modules have been evaluated for completeness.
        The function goes through the top-level steps established by compute_stepwise_module_completeness_for_bin()
        and re-assesses whether steps including other modules are complete. It updates the metabolism completess dictionary accordingly.

        PARAMETERS
        ==========
        mnum : string
            the module number to adjust
        meta_dict_for_bin : dictionary of dictionaries
            metabolism completeness dictionary for the current bin

        RETURNS
        =======
        now_complete : boolean
            whether or not the module is NOW considered "complete" overall based on the threshold fraction of completeness
        """

        num_steps = len(meta_dict_for_bin[mnum]['top_level_step_info'].keys())
        num_complete = 0
        num_nonessential_steps = 0

        for step_num, step_dict in meta_dict_for_bin[mnum]['top_level_step_info'].items():
            if step_dict['includes_modules']:
                # this condition statement has module accessions in it, we need to replace those with True/False and eval
                step_is_present_condition_statement = step_dict['complete']
                module_accessions_to_replace = step_dict['included_module_list']
                for m in module_accessions_to_replace:
                    mod_completeness = meta_dict_for_bin[m]['stepwise_completeness']
                    if mod_completeness >= module_completion_threshold:
                        step_is_present_condition_statement = step_is_present_condition_statement.replace(m, "True")
                    else:
                        step_is_present_condition_statement = step_is_present_condition_statement.replace(m, "False")

                # now evaluate to see if this step is complete
                step_is_present = eval(step_is_present_condition_statement)
                meta_dict_for_bin[mnum]['top_level_step_info'][step_num]['complete'] = step_is_present
                if step_is_present:
                    num_complete += 1

            else:
                if step_dict['complete'] == "nonessential":
                    num_nonessential_steps += 1
                elif step_dict['complete']:
                    num_complete += 1

        mod_stepwise_completeness = num_complete / (num_steps - num_nonessential_steps)
        meta_dict_for_bin[mnum]["stepwise_completeness"] = mod_stepwise_completeness

        now_complete = True if mod_stepwise_completeness >= module_completion_threshold else False
        meta_dict_for_bin[mnum]["stepwise_is_complete"] = now_complete

        return now_complete


    def adjust_pathwise_completeness_for_bin(self, mod, meta_dict_for_bin, module_paths_dict, module_completion_threshold):
        """This function adjusts pathwise completeness of modules that are defined by other modules.

        This can only be done after all other modules have been evaluated for completeness.
        The function uses similar logic as compute_pathwise_module_completeness_for_bin() to re-assess whether steps defined
        by other modules are complete, and updates the metabolism completess dictionary accordingly.

        PARAMETERS
        ==========
        mod : string
            the module number to adjust
        meta_dict_for_bin : dictionary of dictionaries
            metabolism completeness dictionary for the current bin

        RETURNS
        =======
        now_complete : boolean
            whether or not the module is NOW considered "complete" overall based on the threshold fraction of completeness
        """

        for i in range(len(module_paths_dict[mod])):
            p = module_paths_dict[mod][i]
            num_essential_steps_in_path = 0  # note that the len(p) will include nonessential steps; we should count only essential ones
            num_complete_module_steps = 0

            # take previously computed step copy numbers. This list includes all steps except those defined by modules.
            atomic_step_copy_numbers_in_path = meta_dict_for_bin[mod]["num_complete_copies_of_all_paths"][i]
            module_copy_num_should_be_NA = False # flag to indicate whether component modules have a copy number of NA

            for atomic_step in p:
                # module step; we need to count these based on previously computed module completeness
                if atomic_step in meta_dict_for_bin:  # This is a module
                    num_complete_module_steps += meta_dict_for_bin[atomic_step]["pathwise_percent_complete"]
                    num_essential_steps_in_path += 1

                    if meta_dict_for_bin[atomic_step]["pathwise_copy_number"] == "NA":
                        module_copy_num_should_be_NA = True
                    else:
                        atomic_step_copy_numbers_in_path.append(meta_dict_for_bin[atomic_step]["pathwise_copy_number"])
                # non-essential KO, don't count as a step in the path
                elif atomic_step[0] == '-' and not atomic_step == "--":
                    pass
                # single enzymes, protein complexes and '--' steps; were already counted as complete by previous function
                else:
                    num_essential_steps_in_path += 1

            # now we adjust the previous pathway completeness
            old_complete_steps_in_path = meta_dict_for_bin[mod]["pathway_completeness"][i] * num_essential_steps_in_path
            adjusted_num_complete_steps_in_path = old_complete_steps_in_path + num_complete_module_steps
            meta_dict_for_bin[mod]["pathway_completeness"][i] = adjusted_num_complete_steps_in_path / num_essential_steps_in_path

            # now we adjust the path copy number
            if module_copy_num_should_be_NA:
                path_copy_number = "NA"
            else:
                path_copy_number = self.compute_num_complete_copies_of_path(atomic_step_copy_numbers_in_path, module_completion_threshold)
            meta_dict_for_bin[mod]["num_complete_copies_of_all_paths"][i] = path_copy_number

        # after adjusting for all paths, adjust overall module completeness
        meta_dict_for_bin[mod]["pathwise_percent_complete"] = max(meta_dict_for_bin[mod]["pathway_completeness"])
        if meta_dict_for_bin[mod]["pathwise_percent_complete"] > 0:
            meta_dict_for_bin[mod]["most_complete_paths"] = [module_paths_dict[mod][i] for i, pc in enumerate(meta_dict_for_bin[mod]["pathway_completeness"]) if pc == meta_dict_for_bin[mod]["pathwise_percent_complete"]]
        else:
            meta_dict_for_bin[mod]["most_complete_paths"] = []

        now_complete = True if meta_dict_for_bin[mod]["pathwise_percent_complete"] >= module_completion_threshold else False
        meta_dict_for_bin[mod]["pathwise_is_complete"] = now_complete

        # and adjust overall module copy number
        if meta_dict_for_bin[mod]["num_complete_copies_of_most_complete_paths"]:
            meta_dict_for_bin[mod]["pathwise_copy_number"] = max(meta_dict_for_bin[mod]["num_complete_copies_of_most_complete_paths"])
        else:
            meta_dict_for_bin[mod]["pathwise_copy_number"] = 'NA'

        return now_complete


    def add_module_coverage(self, mod, meta_dict_for_bin, profile_db=None, enzymes_of_interest_df=None):
        """This function updates the metabolism dictionary with coverage values for the given module.

        For profile DB input, this function must be called after dbaccess.init_gene_coverage() so that the 
        relevant gene coverage values are initialized.

        NEW KEYS ADDED TO METABOLISM COMPLETENESS DICT
        =======
        "genes_to_coverage"             dictionary of mean coverage in each sample for each gene
                                        coverage = meta_dict_for_bin[module]["genes_to_coverage"][sample][gcid]
        "genes_to_detection"            dictionary of detection in each sample for each gene
                                        detection = meta_dict_for_bin[module]["genes_to_detection"][sample][gcid]
        "average_coverage_per_sample"   dictionary of average mean coverage of all genes in module, per sample
                                        avg_coverage = meta_dict_for_bin[module]["average_coverage_per_sample"][sample]
        "average_detection_per_sample"  dictionary of average detection of all genes in module, per sample
                                        avg_detection = meta_dict_for_bin[module]["average_detection_per_sample"][sample]
        """

        meta_dict_for_bin[mod]["genes_to_coverage"] = {}
        meta_dict_for_bin[mod]["genes_to_detection"] = {}
        meta_dict_for_bin[mod]["average_coverage_per_sample"] = {}
        meta_dict_for_bin[mod]["average_detection_per_sample"] = {}

        if enzymes_of_interest_df.empty and not profile_db:
            raise ConfigError("The add_module_coverage() function cannot work without a properly initialized "
                              "profile database.")

        num_genes = len(meta_dict_for_bin[mod]["gene_caller_ids"])
        if not self.coverage_sample_list:
            # here we establish the same 'sample list' that we use in estimate.add_gene_coverage_to_headers_list()
            if profile_db:
                self.coverage_sample_list = profile_db.p_meta['samples']
            elif enzymes_of_interest_df:
                self.coverage_sample_list = [self.contigs_db_project_name] 

        for s in self.coverage_sample_list:
            meta_dict_for_bin[mod]["genes_to_coverage"][s] = {}
            meta_dict_for_bin[mod]["genes_to_detection"][s] = {}
            coverage_sum = 0
            detection_sum = 0
            for g in meta_dict_for_bin[mod]["gene_caller_ids"]:
                if enzymes_of_interest_df is not None:
                    cov = enzymes_of_interest_df[enzymes_of_interest_df['gene_id'] == g]['coverage'].values[0]
                    det = enzymes_of_interest_df[enzymes_of_interest_df['gene_id'] == g]['detection'].values[0]
                else:
                    cov = profile_db.gene_level_coverage_stats_dict[g][s]['mean_coverage']
                    det = profile_db.gene_level_coverage_stats_dict[g][s]['detection']
                coverage_sum += cov
                detection_sum += det
                meta_dict_for_bin[mod]["genes_to_coverage"][s][g] = cov
                meta_dict_for_bin[mod]["genes_to_detection"][s][g] = det

            if num_genes == 0:
                meta_dict_for_bin[mod]["average_coverage_per_sample"][s] = 0
                meta_dict_for_bin[mod]["average_detection_per_sample"][s] = 0
            else:
                meta_dict_for_bin[mod]["average_coverage_per_sample"][s] = coverage_sum / num_genes
                meta_dict_for_bin[mod]["average_detection_per_sample"][s] = detection_sum / num_genes


    def estimate_for_list_of_splits(self, metabolism_dict_for_list_of_splits, bin_name=None, all_modules_in_db=None,
                                    module_paths_dict=None, module_completion_threshold=0.75, exclude_dashed_reactions=True,
                                    add_coverage=False, quiet=False, genome_mode=False):
        """This is the atomic metabolism estimator function, which builds up the metabolism completeness dictionary for an arbitrary list of splits.

        For example, the list of splits may represent a bin, a single isolate genome, or an entire metagenome.

        The function takes in a metabolism completeness dictionary already initialized with the relevant KOfam hits per module, and updates it
        with the individual steps and completion estimates for each module.

        PARAMETERS
        ==========
        metabolism_dict_for_list_of_splits : dictionary of dictionaries
            the metabolism completeness dictionary of dictionaries for this list of splits. It contains
            one dictionary of module steps and completion information for each module (keyed by module number).
            Calling functions should assign this dictionary to a metabolism superdict with the bin name as a key.
        bin_name : str
            the name of the bin/genome/metagenome that we are working with
        """

        pathwise_complete_mods = set([])
        stepwise_complete_mods = set([])
        mods_def_by_modules = [] # a list of modules that have module numbers in their definitions
        # modules to warn about
        mods_with_unassociated_ko = [] # a list of modules that have "--" steps without an associated KO
        mods_with_nonessential_steps = [] # a list of modules that have nonessential steps like "-K11024"

        # estimate completeness of each module
        for mod in metabolism_dict_for_list_of_splits.keys():
            # pathwise
            mod_is_complete, has_nonessential_step, has_no_ko_step, defined_by_modules \
            = self.compute_pathwise_module_completeness_for_bin(mod, metabolism_dict_for_list_of_splits,
                                                              module_paths_dict, module_completion_threshold,
                                                              exclude_dashed_reactions, all_modules_in_db)

            if mod_is_complete:
                pathwise_complete_mods.add(mod)
            if has_nonessential_step:
                mods_with_nonessential_steps.append(mod)
            if has_no_ko_step:
                mods_with_unassociated_ko.append(mod)
            if defined_by_modules:
                mods_def_by_modules.append(mod)

            # stepwise
            mod_is_complete = self.compute_stepwise_module_completeness_for_bin(mod, metabolism_dict_for_list_of_splits,
                                                                              all_modules_in_db, module_completion_threshold,
                                                                              exclude_dashed_reactions)
            self.compute_stepwise_module_copy_number_for_bin(mod, metabolism_dict_for_list_of_splits,
                                                           all_modules_in_db, module_completion_threshold,
                                                           exclude_dashed_reactions)

            if mod_is_complete:
                stepwise_complete_mods.add(mod)

            if add_coverage:
                self.add_module_coverage(mod, metabolism_dict_for_list_of_splits, profile_db=self.profile_db, 
                                        enzymes_of_interest_df=self.enzymes_of_interest_df)

        # go back and adjust completeness/copy number of modules that are defined by other modules
        if mods_def_by_modules:
            for mod in mods_def_by_modules:
                # pathwise
                mod_is_complete = self.adjust_pathwise_completeness_for_bin(mod, metabolism_dict_for_list_of_splits,
                                                                          module_paths_dict, module_completion_threshold)
                if mod_is_complete:
                    pathwise_complete_mods.add(mod)
                # stepwise
                mod_is_complete = self.adjust_stepwise_completeness_for_bin(mod, metabolism_dict_for_list_of_splits,
                                                                          module_completion_threshold)
                if mod_is_complete:
                    stepwise_complete_mods.add(mod)
                self.adjust_stepwise_copy_number_for_bin(mod, metabolism_dict_for_list_of_splits, module_completion_threshold)


        # estimate redundancy of each module
        for mod in metabolism_dict_for_list_of_splits.keys():
            self.compute_module_redundancy_for_bin(mod, metabolism_dict_for_list_of_splits)


        # notify user of the modules that gave some fishy results -- but only for genome mode because it's too wordy otherwise
        if not quiet and genome_mode:
            if mods_with_nonessential_steps:
                self.run.warning(f"Please note that anvi'o found one or more non-essential steps in the following modules: "
                                 f"{', '.join(mods_with_nonessential_steps)}. At this time, we are not counting these steps "
                                 f"in our percent completion estimates.")

            if mods_with_unassociated_ko:
                self.run.warning(f"Just so you know, while estimating the completeness of some modules, anvi'o saw "
                                 f"'--' in the module DEFINITION. This indicates a step in the pathway that has no "
                                 f"associated enzyme. So we really cannot know just based on gene annotations whether or not this "
                                 f"step is present. By default, anvi'o marks these steps incomplete. But they may not be, "
                                 f"and as a result their modules may be falsely considered incomplete. So it may be in your "
                                 f"interest to go back and take a look at these individual modules to see if you can find the "
                                 f"missing enzyme in some other way. Best of luck to you. Here is the list of modules to check out: "
                                 f"{', '.join(mods_with_unassociated_ko)}")

        if anvio.DEBUG or genome_mode:
            self.run.info("Bin name", bin_name)
            self.run.info("Module completion threshold", module_completion_threshold)
            self.run.info("Number of complete modules (pathwise)", len(pathwise_complete_mods))
            self.run.info("Number of complete modules (stepwise)", len(stepwise_complete_mods))
            if pathwise_complete_mods:
                self.run.info("Pathwise complete modules", ", ".join(sorted(list(pathwise_complete_mods))))
            if stepwise_complete_mods:
                self.run.info("Stepwise complete modules", ", ".join(sorted(list(stepwise_complete_mods))))

        return metabolism_dict_for_list_of_splits

######### REDUNDANCY FUNCTIONS (UNUSED IN NON-JSON OUTPUT) #########

    def compute_naive_redundancy_for_path(self, num_ko_hits_in_path):
        """This function computes a naive redundancy measure for a module path, given the number of hits per KO in the path.

        naive redundancy = # extra hits / len(path) where a hit is "extra" if it is not the first hit to the KO.

        PARAMETERS
        ==========
        num_ko_hits_in_path : list
            stores the number of copies of each enzyme in path
        """

        extra_hits = [x - 1 if x > 1 else 0 for x in num_ko_hits_in_path]
        return sum(extra_hits)/len(num_ko_hits_in_path)


    def compute_copywise_redundancy_for_path(self, num_ko_hits_in_path, aggregation_measure="average"):
        """This function computes redundancy based on the completeness of each extra copy of a path.

        The 'base' redundancy score is determined by the number of extra copies with 100% completeness.
        The completeness measurements of all other extra copies are aggregated (using the aggregation_measure) and
        added to this 'base' redundancy to get the overall path redundancy.

        PARAMETERS
        ==========
        num_ko_hits_in_path : list
            stores the number of copies of each enzyme in path
        """

        accepted_aggregation_measures = ["average", "median", "weighted_sum", "geometric_mean"]
        extra_hits = [x - 1 if x > 1 else 0 for x in num_ko_hits_in_path]
        base_redundancy = min(extra_hits) # number of extra copies of path that are 100% complete
        extra_copy_completeness = []
        # here we get the completeness of every extra copy of the path
        for i in range((base_redundancy+1), max(extra_hits) + 1):
            num_present_kos_in_copy = len([num_hits for num_hits in extra_hits if num_hits >= i])
            extra_copy_completeness.append(num_present_kos_in_copy/len(num_ko_hits_in_path))

        aggregated_completeness = None
        if not extra_copy_completeness: # this handles the case when ALL extra copies are 100% complete
            aggregated_completeness = 0
        else:
            if aggregation_measure == "average":
                aggregated_completeness = statistics.mean(extra_copy_completeness)
            elif aggregation_measure == "median":
                aggregated_completeness = statistics.median(extra_copy_completeness)
            elif aggregation_measure == "weighted_sum":
                aggregated_completeness = 0
                for c in range(len(extra_copy_completeness)):
                    aggregated_completeness += 1/(c+1) * extra_copy_completeness[c]
            elif aggregation_measure == "geometric_mean":
                aggregated_completeness = stats.gmean(extra_copy_completeness)
            else:
                raise ConfigError("The function compute_copywise_redundancy_for_path() doesn't know how to handle the aggregation measure '%s'. "
                                  "Accepted aggregation measures include: %s " % (aggregation_measure, ", ".join(accepted_aggregation_measures)))

        return (base_redundancy + aggregated_completeness), extra_copy_completeness


    def compute_entropy_weighted_redundancy_for_bin(self, num_ko_hits_in_path):
        """This function computes naive redundancy but weights it by the entropy of the hit distribution.

        PARAMETERS
        ==========
        num_ko_hits_in_path : list
            stores the number of copies of each enzyme in path
        """

        extra_hits = [x - 1 if x > 1 else 0 for x in num_ko_hits_in_path]
        total_extra_hits = sum(extra_hits)
        num_kos = len(num_ko_hits_in_path)
        naive_redundancy = total_extra_hits/num_kos
        if all(e == 0 for e in extra_hits):
            return 0.0
        entropy = stats.entropy(extra_hits)
        max_entropy_distribution = [total_extra_hits // num_kos] * num_kos
        for i in range(total_extra_hits % num_kos):
            max_entropy_distribution[i] += 1
        max_entropy = stats.entropy(max_entropy_distribution)
        # avoid divide by 0
        max_entropy += 1e-20

        return naive_redundancy * entropy/max_entropy


    def compute_num_complete_copies_of_path(self, copy_num_of_atomic_steps, module_completion_threshold):
        """This function computes the number of copies of a path that are >= x% complete,
        where x is the module completeness threshold.

        It does this based on the provided list in which each entry is the number of copies of
        each atomic step in the path.
        - first, these copy numbers are ordered (descending order)
        - then, we compute N, the number of steps needed to make the path at least X complete, where
          X is the module completeness threshold
        - finally, we loop from i=1 to the maximum number of hits. Each time, if the number x of steps
          with hit count >= i is x >= N, we add 1 to our count of path copy numbers.
        - the final count of path copy numbers is returned.

        PARAMETERS
        ==========
        copy_num_of_atomic_steps : list
            stores the number of copies of each step in path

        RETURNS
        ==========
        copy_number : int
            number of copies of path which are at least X complete, where X is module completeness threshold
        """

        import math
        path_length = len(copy_num_of_atomic_steps)
        num_enzymes_needed = math.ceil(module_completion_threshold * path_length)  # N
        copy_num_of_atomic_steps.sort(reverse=True)

        copy_number = 0
        for i in range(1, copy_num_of_atomic_steps[0]+1):
            x = len([h for h in copy_num_of_atomic_steps if h >= i])
            if x >= num_enzymes_needed:
                copy_number += 1

        return copy_number


    def compute_module_redundancy_for_bin(self, mnum, meta_dict_for_bin):
        """This function calculates the redundancy of the specified module within the given bin metabolism dictionary.

        Each module can have multiple paths, but (in most cases) we only compute redundancy on the paths with the highest completeness
        (stored under the "most_complete_paths" key). If there are no paths in this list (which only happens when there
        are 0 KOfam hits to the module), then we do not compute redundancy.

        PARAMETERS
        ==========
        mnum : string
            module number to work on
        meta_dict_for_bin : dictionary of dictionaries
            metabolism completeness dict for the current bin, to be modified in-place

        """

        meta_dict_for_bin[mnum]["naive_redundancy"] = []
        meta_dict_for_bin[mnum]["copywise_average"] = []
        meta_dict_for_bin[mnum]["copywise_completeness_distributions"] = []
        meta_dict_for_bin[mnum]["copywise_median"] = []
        meta_dict_for_bin[mnum]["copywise_weighted-sum"] = []
        meta_dict_for_bin[mnum]["copywise_geometric-mean"] = []
        meta_dict_for_bin[mnum]["entropy_weighted"] = []

        paths_of_highest_completeness = meta_dict_for_bin[mnum]["most_complete_paths"]
        if not paths_of_highest_completeness:
            return

        for p in paths_of_highest_completeness:
            p = self.split_module_path_into_individual_essential_components(p)
            num_hits_per_kofam = [len(meta_dict_for_bin[mnum]["kofam_hits"][k]) if k in meta_dict_for_bin[mnum]["kofam_hits"] else 0 for k in p]

            # for now, we will try a bunch of different redundancy calculations and put them all into the dictionary until we find the ones we like
            meta_dict_for_bin[mnum]["naive_redundancy"].append(self.compute_naive_redundancy_for_path(num_hits_per_kofam))
            cw_avg_redundancy, copy_completeness_distribution = self.compute_copywise_redundancy_for_path(num_hits_per_kofam, aggregation_measure="average")
            meta_dict_for_bin[mnum]["copywise_average"].append(cw_avg_redundancy)
            meta_dict_for_bin[mnum]["copywise_completeness_distributions"].append(copy_completeness_distribution)
            cw_med_redundancy, copy_completeness_distribution = self.compute_copywise_redundancy_for_path(num_hits_per_kofam, aggregation_measure="median")
            meta_dict_for_bin[mnum]["copywise_median"].append(cw_med_redundancy)
            cw_ws_redundancy, copy_completeness_distribution = self.compute_copywise_redundancy_for_path(num_hits_per_kofam, aggregation_measure="weighted_sum")
            meta_dict_for_bin[mnum]["copywise_weighted-sum"].append(cw_ws_redundancy)
            cw_gm_redundancy, copy_completeness_distribution = self.compute_copywise_redundancy_for_path(num_hits_per_kofam, aggregation_measure="geometric_mean")
            meta_dict_for_bin[mnum]["copywise_geometric-mean"].append(cw_gm_redundancy)
            meta_dict_for_bin[mnum]["entropy_weighted"].append(self.compute_entropy_weighted_redundancy_for_bin(num_hits_per_kofam))

        return


    def get_step_copy_number(self, step_string, enzyme_hit_counts, all_modules_in_db, exclude_dashed_reactions=True):
        """This function recursively calculates the copy number of a step in a module.

        It parses the definition string of a step and recurses as needed to compute copy number of
        substeps. Copy numbers of any substeps are mathematically combined to obtain a copy number for
        the step as a whole.

        The key base case in the recursion is an individual enzyme accession, for which the copy number
        is simply the number of times it is annotated in the sample (which we obtain from the enzyme_hit_counts dictionary).

        Combining copy numbers works as follows: If two enzymes (or substeps) are connected by an AND, then we need both, so
        we take the minimum copy number of both of them. If they are connected by an OR, then we can use either, so we can sum
        their copy numbers. In doing this, we follow correct order of operations, as established by any parentheses in the step definition.

        In short, this function accomplishes the same thing as modifying the step definition by replacing spaces and '+' signs with min()
        operations, replacing commas with + operations, and replacing enzyme accessions with their corresponding hit counts; then returning
        the value obtained by evaluating the resulting arithmetic expression.

        Some steps are defined by other modules. When module accessions are found, we initially treat them as having a copy number of 0, but
        we re-compute the copy number of the module later once we have the overall copy number of all other modules (and then we use the
        component module's copy number in the calculation instead).

        PARAMETERS
        ==========
        step_string : str
            A string containing the definition of one step (or substep) from a module
        enzyme_hit_counts : dict
            Keys are enzyme accessions, values are the number of times the enzyme was annotated in the current sample

        RETURNS
        ==========
        The copy number (int) of the given step/substep
        """

        # first, eliminate non-essential KOs from the step definition so they won't be considered
        step_string = self.remove_nonessential_enzymes_from_module_step(step_string, exclude_dashed_reactions)

        # sometimes a step will have commas outside parentheses. In this case, we need to split those first for proper order of operations
        step_list = utils.split_by_delim_not_within_parens(step_string, ",")
        if len(step_list) > 1:
            added_step_count = 0
            # call recursively on each split
            for s in step_list:
                added_step_count += self.get_step_copy_number(s, enzyme_hit_counts, all_modules_in_db, exclude_dashed_reactions)
            # combine results using addition and return
            return added_step_count

        # complex case - parentheses surround substeps, which need to be counted recursively and appropriately combined
        if '(' in step_string:
            open_parens_idx = step_string.index('(') # first (outermost) parenthesis
            close_parens_idx = None

            # find matching parenthesis
            i = open_parens_idx + 1
            parens_level = 1
            while not close_parens_idx:
                if step_string[i] == '(':
                    parens_level += 1
                if step_string[i] == ')':
                    parens_level -= 1

                    if parens_level == 0:
                        close_parens_idx = i
                i += 1

            # call recursively on string within outermost parentheses
            sub_step = step_string[open_parens_idx+1:close_parens_idx]
            sub_copy_num = self.get_step_copy_number(sub_step, enzyme_hit_counts, all_modules_in_db, exclude_dashed_reactions)

            # parse the rest of the string and combine with the copy number of the stuff within parentheses
            step_copy_num = None
            # handle anything prior to parentheses
            if open_parens_idx > 0:
                previous_str = step_string[:open_parens_idx]

                previous_steps = previous_str[:-1]
                prev_copy = self.get_step_copy_number(previous_steps, enzyme_hit_counts, all_modules_in_db, exclude_dashed_reactions)

                combo_element = previous_str[-1]
                if combo_element == ',': # OR
                    step_copy_num = (prev_copy + sub_copy_num)
                if combo_element == ' ' or combo_element == '+': # AND
                    step_copy_num = min(prev_copy,sub_copy_num)

            # handle anything following parentheses
            if close_parens_idx < len(step_string) - 1:
                post_steps = step_string[close_parens_idx+2:]
                post_copy = self.get_step_copy_number(post_steps, enzyme_hit_counts, all_modules_in_db, exclude_dashed_reactions)

                combo_element = step_string[close_parens_idx+1]
                if step_copy_num is None:
                    # no previous clause, so we only combine the parenthetical clause and what comes after
                    if combo_element == ',': # OR
                        step_copy_num = (sub_copy_num + post_copy)
                    if combo_element == ' ' or combo_element == '+': # AND
                        step_copy_num = min(sub_copy_num,post_copy)
                else:
                    # we have to combine the post clause with the already-combined previous clause
                    # and parenthetical clause
                    if combo_element == ',': # OR
                        step_copy_num += post_copy
                    if combo_element == ' ' or combo_element == '+': # AND
                        step_copy_num = min(step_copy_num,post_copy)

            # handle edge case where parentheses circles entire step
            if (open_parens_idx == 0) and (close_parens_idx == len(step_string) - 1):
                step_copy_num = sub_copy_num

            return step_copy_num

        # simple case - no substeps within parentheses
        else:
            if ',' in step_string: # OR - combine copy numbers using addition
                or_splits = step_string.split(',')
                added_step_count = 0
                for s in or_splits:
                    added_step_count += self.get_step_copy_number(s, enzyme_hit_counts, all_modules_in_db, exclude_dashed_reactions)
                return added_step_count

            elif ' ' in step_string or '+' in step_string: # AND - combine copy numbers using min()
                and_splits = step_string.replace('+', ' ').split(' ')
                min_step_count = None
                for s in and_splits:
                    s_count = self.get_step_copy_number(s, enzyme_hit_counts, all_modules_in_db, exclude_dashed_reactions)
                    if min_step_count is None:
                        min_step_count = s_count # make first step the minimum
                    min_step_count = min(min_step_count, s_count)
                return min_step_count

            # base cases
            elif '-' in step_string:
                if step_string == '--': # no KO profile => no copy number (unless user wants to ignore these, in which case
                        return 0        # they were already removed by remove_nonessential_enzymes_from_module_step() above)
                else: # contains non-essential KO, should never happen because we eliminated them above
                    raise ConfigError(f"Something is very wrong, because the get_step_copy_number() function found a nonessential "
                                      f"enzyme in the step definition {step_string}")
            elif step_string == '': # entire step was nonessential KO, skip computation
                return None
            else: # accession
                if all_modules_in_db and step_string in all_modules_in_db: # module
                    if step_string in enzyme_hit_counts: # we are currently adjusting, and know the module copy number
                        return enzyme_hit_counts[step_string]
                    else: # return 0 for now, will be adjusted later
                        return 0
                else: # enzyme
                    if step_string not in enzyme_hit_counts:
                        return 0
                    return enzyme_hit_counts[step_string]


    def are_enzymes_indirect_alternatives_within_step(self, enzyme_list, step):
        """An overly simplistic function to determine whether the relationship between the provided alternative
        enzymes in the given step is indirect.

        To do this, it simply walks through the step definition string to determine whether each pair of enzymes is separated by
        a character symbolizing a more complex relationship. That is, they are not separated only by commas and other enzymes (which
        indicates a direct relationship, as in the two enzymes are synonymous in the context of the metabolic pathway).

        For example, within the step (((K01657+K01658,K13503,K13501,K01656) K00766),K13497), the direct alternatives include
        K13503, K13501, and K01656. K01657 and K01658 are indirect alternatives to each other because they are two
        components of the same enzyme, while K01658 and K00766 are indirect because they catalyze two separate reactions in
        an alternative branch of the step.

        This algorithm is not perfect at identifying all indirect relationships - for instance, given K01658 and K13503 it will
        wrongly suggest they are direct alternatives. However, it is meant to be used only for identifying putative edge cases
        for the `get_dereplicated_enzyme_hits_for_step_in_module()` function, and it works well enough for that.

        PARAMETERS
        ==========
        enzyme_list : list of enzyme accessions
            the alternative enzymes to process
        step : string
            the definition string of the relevant step

        RETURNS
        =======
        contains_indirect : Boolean
            True if the list of provided enzymes contains those that are indirect alternatives within the given step.
        """

        enzyme_data = {e : {'index': step.index(e),
                                    'direct_alts': [],
                                    'indirect_alts': []} for e in enzyme_list}

        contains_indirect = False
        # get enzyme-specific list of alternatives
        for e in enzyme_list:
            for z in enzyme_list:
                if e != z:
                    e_index = enzyme_data[e]['index']
                    z_index = enzyme_data[z]['index']
                    indirect_alternatives = False

                    # indirect alts have a space, parentheses, or plus/minus sign between them
                    for c in step[min(e_index, z_index):max(e_index, z_index)]:
                        if c in [' ', '(', ')', '+', '-']:
                            indirect_alternatives = True

                    if indirect_alternatives:
                        enzyme_data[e]['indirect_alts'].append(z)
                        contains_indirect = True
                    else:
                        enzyme_data[e]['direct_alts'].append(z)

        return contains_indirect


    def get_dereplicated_enzyme_hits_for_step_in_module(self, meta_dict_for_mnum, step_to_focus_on, mnum,
                                                       add_copy_number=False):
        """This function returns a dictionary of enzyme accessions matched to the number of hits, with duplicate hits to the
        same gene removed, for the provided step in a metabolic pathway.

        Depreplicating the gene calls is necessary because the same gene can be annotated with multiple alternative enzymes for the
        same reaction, and we don't want these annotations to be double-counted in the stepwise copy number calculation.

        PARAMETERS
        ==========
        meta_dict_for_mnum : dictionary of dictionaries
            metabolism completeness dict for the current bin and metabolic module
        step_to_focus_on : string
            which step in the module to resolve alternative enzymes for, passed as a definition string for the step.
        mnum : string
            module ID (used only for warning output)

        RETURNS
        =======
        derep_enzyme_hits : dictionary
            matches enzyme accession to number of hits to unique genes
        """

        derep_enzyme_hits = {k : len(meta_dict_for_mnum["kofam_hits"][k]) for k in meta_dict_for_mnum["kofam_hits"] if k in step_to_focus_on}

        # map gene caller IDs to enzyme accessions
        gene_calls_to_enzymes = {gcid : [] for gcid in meta_dict_for_mnum['gene_caller_ids']}
        for enzyme, gene_list in meta_dict_for_mnum['kofam_hits'].items():
            for g in gene_list:
                if enzyme in step_to_focus_on:
                    gene_calls_to_enzymes[g].append(enzyme)

        for gcid, enzymes in gene_calls_to_enzymes.items():
            if len(enzymes) > 1:
                # simple solution (only works well for enzymes that are direct alternatives)
                # for each duplicated gene, we arbitrarily keep only the hit to the first enzyme
                # and for all other annotations, we reduce the count of hits by one
                for acc in enzymes[1:]:
                    derep_enzyme_hits[acc] -= 1

                if self.are_enzymes_indirect_alternatives_within_step(enzymes, step_to_focus_on) and add_copy_number:
                    enz_str = ", ".join(enzymes)
                    self.run.warning(f"The gene call {gcid} has multiple annotations to alternative enzymes "
                                     f"within the same step of a metabolic pathway ({enz_str}), and these enzymes "
                                     f"unfortunately have a complex relationship. The affected module is {mnum}, and "
                                     f"here is the step in question: {step_to_focus_on}. We arbitrarily kept only one of "
                                     f"the annotations to this gene in order to avoid inflating the step's copy number, "
                                     f"but due to the complex relationship between these alternatives, this could mean "
                                     f"that the copy number for this step is actually too low. Please heed this warning "
                                     f"and double check the stepwise copy number results for {mnum} and other pathways "
                                     f"containing gene call {gcid}.")

        return derep_enzyme_hits


    def compute_stepwise_module_copy_number_for_bin(self, mnum, meta_dict_for_bin, all_modules_in_db,
                                                   module_completion_threshold, exclude_dashed_reactions=True):
        """This function calculates the copy number of the specified module within the given bin metabolism dictionary.

        It goes through the top-level steps established by compute_stepwise_module_completeness_for_bin() and determines the
        copy number of each step. Then, the overall module copy number is calculated as the minimum copy number of all steps.

        PARAMETERS
        ==========
        mnum : string
            module number to work on
        meta_dict_for_bin : dictionary of dictionaries
            metabolism completeness dict for the current bin, to be modified in-place

        NEW KEYS ADDED TO METABOLISM COMPLETENESS DICT
        =======
        "stepwise_copy_number"         the stepwise copy number of the module

        [keys added in "top_level_step_info" dictionary]
            "copy_number"              the copy number of an individual step
        """

        all_step_copy_nums = []
        for key in meta_dict_for_bin[mnum]["top_level_step_info"]:
            if not meta_dict_for_bin[mnum]["top_level_step_info"][key]["includes_modules"]:
                step_string = meta_dict_for_bin[mnum]["top_level_step_info"][key]["step_definition"]
                enzyme_hits_dict = self.get_dereplicated_enzyme_hits_for_step_in_module(meta_dict_for_bin[mnum], step_string, mnum)

                step_copy_num = self.get_step_copy_number(step_string, enzyme_hits_dict, all_modules_in_db, exclude_dashed_reactions)
                meta_dict_for_bin[mnum]["top_level_step_info"][key]["copy_number"] = step_copy_num
                if step_copy_num is not None: # avoid taking minimum of None values (from non-essential steps)
                    all_step_copy_nums.append(step_copy_num)

        if all_step_copy_nums:
            module_stepwise_copy_num = min(all_step_copy_nums)
        else:
            module_stepwise_copy_num = None
        meta_dict_for_bin[mnum]["stepwise_copy_number"] = module_stepwise_copy_num


    def adjust_stepwise_copy_number_for_bin(self, mnum, meta_dict_for_bin, module_completion_threshold):
        """This function adjusts stepwise copy number of modules that are defined by other modules.

        This can only be done after all other modules have had their copy numbers calculated and added to the metabolism dictionary
        by the function compute_stepwise_module_copy_number_for_bin().

        The function goes through the top-level steps in the module and re-computes copy number for steps that include other modules.
        Then it re-calculates the overall module copy number as the minimum copy number of all steps. It updates the metabolism completess
        dictionary accordingly.

        PARAMETERS
        ==========
        mnum : string
            the module number to adjust
        meta_dict_for_bin : dictionary of dictionaries
            metabolism completeness dictionary for the current bin
        """

        enzyme_hits_dict = {k : len(meta_dict_for_bin[mnum]["kofam_hits"][k]) for k in meta_dict_for_bin[mnum]["kofam_hits"] }

        all_step_copy_nums = []
        for key in meta_dict_for_bin[mnum]["top_level_step_info"]:
            # re-calculate ONLY for steps with modules in definition
            if meta_dict_for_bin[mnum]["top_level_step_info"][key]["includes_modules"]:
                step_string = meta_dict_for_bin[mnum]["top_level_step_info"][key]["step_definition"]

                for included_module in meta_dict_for_bin[mnum]["top_level_step_info"][key]["included_module_list"]:
                    enzyme_hits_dict[included_module] = meta_dict_for_bin[included_module]["stepwise_copy_number"]

                step_copy_num = self.get_step_copy_number(step_string, enzyme_hits_dict, meta_dict_for_bin)
                meta_dict_for_bin[mnum]["top_level_step_info"][key]["copy_number"] = step_copy_num

            # take minimum over all steps, even those not defined by modules
            if meta_dict_for_bin[mnum]["top_level_step_info"][key]["copy_number"] is not None: # avoid taking minimum of None values (from non-essential steps)
                all_step_copy_nums.append(meta_dict_for_bin[mnum]["top_level_step_info"][key]["copy_number"])

        module_stepwise_copy_num = min(all_step_copy_nums)
        meta_dict_for_bin[mnum]["stepwise_copy_number"] = module_stepwise_copy_num


## STATIC FUNCTIONS
def module_definition_to_enzyme_accessions(mod_definition):
    """Parses a module definition string into a list of enzyme accessions."""

    # anything that is not (),-+ should be converted to spaces, then we can split on the spaces to get the accessions
    mod_definition = re.sub('[\(\)\+\-,]', ' ', mod_definition).strip()
    acc_list = re.split(r'\s+', mod_definition)

    return acc_list
