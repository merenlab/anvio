import os
import collections

import anvio
import anvio.utils as utils

from anvio.errors import ConfigError


class KeggContext(object):
    """Base class to define shared functions and file paths for all KEGG operations."""

    def __init__(self, args):
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        # default data directory will be called KEGG and will store the KEGG Module data as well
        self.default_kegg_dir = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/KEGG')
        self.kegg_data_dir = os.path.abspath(A('kegg_data_dir') or self.default_kegg_dir)
        self.user_input_dir = os.path.abspath(A('user_modules')) if A('user_modules') else None
        self.only_user_modules = A('only_user_modules')
        self.orphan_data_dir = os.path.join(self.kegg_data_dir, "orphan_data")
        self.kegg_module_data_dir = os.path.join(self.kegg_data_dir, "modules")
        self.kegg_hmm_data_dir = os.path.join(self.kegg_data_dir, "HMMs")
        self.pathway_data_dir = os.path.join(self.kegg_data_dir, "pathways")
        self.brite_data_dir = os.path.join(self.kegg_data_dir, "BRITE")
        self.binary_relation_data_dir = os.path.join(self.kegg_data_dir, "binary_relations")

        # The 'KEGG/map_images' directory has a structure of nested directories. 'map_images'
        # contains 'png' for image files and 'kgml' for XML mapping files. Within both 'png' and
        # 'kgml' are directories, '1x' and '2x', for lower and higher resolution maps. 'png/1x'
        # contains 5 directories of image files highlighting different things: 'map', 'ko', 'ec',
        # 'rn', and 'org'. 'png/2x' contains 1 directory, 'map', as higher resolution images are
        # only available for manually drawn maps. 'kgml/1x' and 'kgml/2x' each contain 4 directories
        # of XML files that allow modification of different lower and higher resolution maps: 'ko',
        # 'ec', 'rn', and 'org'.
        self.map_image_data_dir = os.path.join(self.kegg_data_dir, "map_images")
        self.png_dir = os.path.join(self.map_image_data_dir, "png")
        self.kgml_dir = os.path.join(self.map_image_data_dir, "kgml")
        self.png_1x_dir = os.path.join(self.png_dir, "1x")
        self.png_2x_dir = os.path.join(self.png_dir, "2x")
        self.png_1x_map_dir = os.path.join(self.png_1x_dir, "map")
        self.png_1x_ko_dir = os.path.join(self.png_1x_dir, "ko")
        self.png_1x_ec_dir = os.path.join(self.png_1x_dir, "ec")
        self.png_1x_rn_dir = os.path.join(self.png_1x_dir, "rn")
        self.png_1x_org_dir = os.path.join(self.png_1x_dir, "org")
        self.png_2x_map_dir = os.path.join(self.png_2x_dir, "map")
        self.kgml_1x_dir = os.path.join(self.kgml_dir, "1x")
        self.kgml_2x_dir = os.path.join(self.kgml_dir, "2x")
        self.kgml_1x_ko_dir = os.path.join(self.kgml_1x_dir, "ko")
        self.kgml_1x_ec_dir = os.path.join(self.kgml_1x_dir, "ec")
        self.kgml_1x_rn_dir = os.path.join(self.kgml_1x_dir, "rn")
        self.kgml_1x_org_dir = os.path.join(self.kgml_1x_dir, "org")
        self.kgml_2x_ko_dir = os.path.join(self.kgml_2x_dir, "ko")
        self.kgml_2x_ec_dir = os.path.join(self.kgml_2x_dir, "ec")
        self.kgml_2x_rn_dir = os.path.join(self.kgml_2x_dir, "rn")
        self.kgml_2x_org_dir = os.path.join(self.kgml_2x_dir, "org")

        self.quiet = A('quiet') or False
        self.just_do_it = A('just_do_it')

        # shared variables for all KEGG subclasses
        self.kofam_hmm_file_path = os.path.join(self.kegg_hmm_data_dir, "Kofam.hmm") # file containing concatenated KOfam hmms
        self.stray_ko_hmm_file_path = os.path.join(self.orphan_data_dir, "anvio_hmm_profiles_for_stray_KOs.hmm") # anvi'o-generated concatenated hmms for stray KOs
        self.stray_ko_hmms_from_kegg = os.path.join(self.orphan_data_dir, "hmm_profiles_with_kofams_with_no_threshold.hmm") # original concatentated hmms for stray KOs
        self.ko_list_file_path = os.path.join(self.kegg_data_dir, "ko_list.txt")
        self.stray_ko_thresholds_file = os.path.join(self.orphan_data_dir, "estimated_thresholds_for_stray_kos.txt")
        self.kegg_module_file = os.path.join(self.kegg_data_dir, "modules.keg")
        self.kegg_pathway_file = os.path.join(self.kegg_data_dir, "pathways.keg")
        self.kegg_brite_hierarchies_file = os.path.join(self.kegg_data_dir, "hierarchies.json")
        self.kegg_brite_pathways_file = os.path.join(self.kegg_data_dir, "br08901.json")
        self.kegg_modules_db_path = os.path.join(self.kegg_data_dir, "MODULES.db")
        self.kegg_binary_relation_files = {('KO', 'EC'): "ko2ec.xl", ('KO', 'RN'): "ko2rn.xl"}
        self.kegg_pathway_list_file = os.path.join(self.kegg_data_dir, "pathway_list.tsv")
        self.kegg_map_image_kgml_file = os.path.join(self.kegg_data_dir, "map_kgml.tsv")

        if self.user_input_dir:
            self.user_module_data_dir = os.path.join(self.user_input_dir, "modules")
            self.user_modules_db_path = os.path.join(self.user_input_dir, "USER_MODULES.db")

        # sanity check for incompatible arguments
        if A('kegg_data_dir') and A('only_user_modules'):
            raise ConfigError("The options --kegg-data-dir and --only-user-modules are incompatible. Please figure out which one you "
                              "want and try again :)")

        # sanity check to prevent automatic overwriting of non-default kegg data dir
        if self.__class__.__name__ in ['KeggSetup'] and not self.user_input_dir:
            if os.path.exists(self.kegg_data_dir) and self.kegg_data_dir != self.default_kegg_dir:
                raise ConfigError(f"You are attempting to set up KEGG in a non-default data directory ({self.kegg_data_dir}) which already exists. "
                                  f"To avoid automatically deleting a directory that may be important to you, anvi'o refuses to get rid of "
                                  f"directories that have been specified with --kegg-data-dir. If you really want to get rid of this "
                                  f"directory and replace it with the KEGG archive data, then please remove the directory yourself using "
                                  f"a command like `rm -r {self.kegg_data_dir}`. We are sorry to make you go through this extra trouble, but it really is "
                                  f"the safest way to handle things.")


    def setup_ko_dict(self, exclude_threshold=True, suppress_warnings=False):
        """The purpose of this function is to process the ko_list file into usable form by KEGG sub-classes.

        The ko_list file (which is downloaded along with the KOfam HMM profiles) contains important
        information for each KEGG Orthology number (KO, or knum), incuding pre-defined scoring thresholds
        for limiting HMM hits and annotation information.

        It looks something like this:

        knum    threshold    score_type    profile_type    F-measure    nseq    nseq_used    alen    mlen    eff_nseq    re/pos    definition
        K00001    329.57    domain    trim    0.231663    1473    1069    1798    371    17.12    0.590    alcohol dehydrogenase [EC:1.1.1.1]

        Since this information is useful for both the setup process (we need to know all the knums) and HMM process,
        all KEGG subclasses need to have access to this dictionary.

        This is a dictionary (indexed by knum) of dictionaries(indexed by column name).
        Here is an example of the dictionary structure:
        self.ko_dict["K00001"]["threshold"] = 329.57

        PARAMETERS
        ==========
        exclude_threshold : Boolean
            If this is true, we remove KOs without a bitscore threshold from the ko_dict
        suppress_warnings : Boolean
            If this is true, we don't print the warning message about stray KOs
        """

        self.ko_dict = utils.get_TAB_delimited_file_as_dictionary(self.ko_list_file_path)
        self.ko_skip_list, self.ko_no_threshold_list = self.get_ko_skip_list()

        # if we are currently setting up KOfams, we should generate a text file with the ko_list entries
        # of the KOs that have no scoring threshold
        if self.__class__.__name__ in ['KeggSetup', 'KOfamDownload']:
            stray_ko_dict = {ko:self.ko_dict[ko] for ko in self.ko_skip_list}
            stray_ko_dict.update({ko:self.ko_dict[ko] for ko in self.ko_no_threshold_list})

            if not os.path.exists(self.orphan_data_dir): # should not happen but we check just in case
                raise ConfigError(f"Hmm. Something is out of order. The orphan data directory {self.orphan_data_dir} does not exist "
                                  f"yet, but it needs to in order for the setup_ko_dict() function to work.")
            stray_ko_path = os.path.join(self.orphan_data_dir, "kofams_with_no_threshold.txt")
            stray_ko_headers = ["threshold","score_type","profile_type","F-measure","nseq","nseq_used","alen","mlen","eff_nseq","re/pos", "definition"]
            utils.store_dict_as_TAB_delimited_file(stray_ko_dict, stray_ko_path, key_header="knum", headers=stray_ko_headers)

        [self.ko_dict.pop(ko) for ko in self.ko_skip_list]
        if exclude_threshold:
            [self.ko_dict.pop(ko) for ko in self.ko_no_threshold_list]
        else:
            if not suppress_warnings:
                self.run.warning("FYI, we are including KOfams that do not have a bitscore threshold in the analysis.")


    def setup_stray_ko_dict(self, add_entries_to_regular_ko_dict=False):
        """This class sets up a dictionary of predicted bit score thresholds for stray KOs, if possible.

        Those predicted thresholds are generated during `anvi-setup-kegg-data --include-stray-KOs`
        (see KOfamDownload.process_all_stray_kos()), and are stored in a file that looks like this:

        knum	threshold	score_type	definition
        K11700	800.4	full	poly(A) RNA polymerase Cid12 [EC:2.7.7.19]
        K14747_anvio_version	1054.2	full	benzoylacetate-CoA ligase [EC:6.2.1.-]

        The dictionary structure is identical to that of self.ko_dict. Note that the `knum` column can contain
        normal KEGG Ortholog accessions (for KOs whose HMMs we haven't updated) and accessions that end with
        STRAY_KO_ANVIO_SUFFIX (for KOs that we created new models for).

        If thresholds have not been predicted, then this function throws an error.

        Parameters
        ==========
        add_entries_to_regular_ko_dict : Boolean
            If True, we don't create a separate self.stray_ko_dict but instead add the stray KOs to the
            regular self.ko_dict attribute. Useful if you don't need to keep the two sets separate.
        """

        if os.path.exists(self.stray_ko_thresholds_file):
            if add_entries_to_regular_ko_dict:
                stray_kos = utils.get_TAB_delimited_file_as_dictionary(self.stray_ko_thresholds_file)
                self.ko_dict.update(stray_kos)
                # initialize it to None so that things don't break if we try to access this downstream
                self.stray_ko_dict = None
            else:
                self.stray_ko_dict = utils.get_TAB_delimited_file_as_dictionary(self.stray_ko_thresholds_file)
        else:
            raise ConfigError(f"You've requested to include stray KO models in your analysis, but we cannot find the "
                              f"estimated bit score thresholds for these models, which can be generated during "
                              f"`anvi-setup-kegg-data` and stored at the following path: {self.stray_ko_thresholds_file}. This "
                              f"means that `anvi-setup-kegg-data` was run without the `--include-stray-KOs` flag for the KEGG "
                              f"data directory that you are using. You have two options: 1) give up on including these models and "
                              f"re-run your command without the `--include-stray-KOs` flag, or 2) change the KEGG data that you "
                              f"are using, which could be as simple as specifying a new `--kegg-data-dir` or as complex as re-running "
                              f"`anvi-setup-kegg-data` to obtain a dataset with predicted thresholds for stray KO models.")


    def get_ko_skip_list(self):
        """The purpose of this function is to determine which KO numbers have no associated data or just no score threshold in the ko_list file.

        That is, their ko_list entries look like this, with hypens in all but the first and last columns:

        K14936    -    -    -    -    -    -    -    -    -    -    small nucleolar RNA snR191
        K15035    -    -    -    -    -    -    -    -    -    -    transfer-messenger RNA
        K15841    -    -    -    -    -    -    -    -    -    -    small regulatory RNA GlmY
        K15851    -    -    -    -    -    -    -    -    -    -    quorum regulatory RNA Qrr
        K16736    -    -    -    -    -    -    -    -    -    -    bantam
        K16863    -    -    -    -    -    -    -    -    -    -    microRNA 21

        These are RNAs.

        Or, their ko_list entries look like this, with no score threshold (but the rest of the data is not completely blank):

        K23749 - - - - 1 1 2266 2266 0.39 0.592 spectinabilin polyketide synthase system NorC [EC:2.3.1.290]

        Returns:
        skip_list  list of strings, each string is a KO number that has no associated data (ie, RNAs)
        no_threshold_list   list of strings, each string is a KO number that has no scoring threshold
        """

        col_names_to_check = ["threshold","score_type","profile_type","F-measure","nseq","nseq_used","alen","mlen","eff_nseq","re/pos"]
        skip_list = []
        no_threshold_list = []
        for k in self.ko_dict.keys():
            should_skip = True
            no_threshold = False
            for c in col_names_to_check:
                if not self.ko_dict[k][c] == "-":
                    should_skip = False
                    break # here we stop checking this KO num because we already found a value in our columns of interest

                if c == "threshold":
                    no_threshold = True # if we got to this line of code, there is a '-' in the threshold column
            if should_skip: # should be True unless we found a value above
                skip_list.append(k)
            elif no_threshold:
                no_threshold_list.append(k)
        return skip_list, no_threshold_list


    @staticmethod
    def invert_brite_json_dict(brite_dict):
        """Invert a BRITE hierarchy dict loaded from a json file into a dict keyed by KEGG entries.

        There are only two keys expected in a BRITE json file, 'name' and 'children'. The value for
        'name' is a string and the value for 'children' is a list of dicts.

        Here is an example of what the beginning of the json dict looks like for 'br08902 BRITE
        Hierarchy Files', the 'hierarchy of all existing hierarchies':
           {
             "name": "br08902",
             "children": [
               {
                 "name": "Pathway and Brite",
                 "children": [
                   {
                     "name": "Pathway maps",
                     "children": [
                       {
                         "name": "br08901  KEGG pathway maps"
                       }
                     ]
                   }, ...
        Observe that innermost dicts only have a single entry keyed by 'name'.

        Here is the corresponding entry in the returned dict for the item in the example:
            'br08901  KEGG pathway maps':
                [['br08902', 'Pathway and Brite', 'Pathway maps']]
        The value is a list of lists because an item can occur multiple times in the same hierarchy.

        PARAMETERS
        ==========
        brite_dict : dict
            dict loaded from BRITE hierarchy json file

        RETURNS
        =======
        categorization_dict : dict
            dict of entry categorizations in BRITE hierarchy
        """

        children_stack = collections.deque()
        children_stack.append(([brite_dict['name']], brite_dict['children']))
        categorization_dict = {}
        while children_stack:
            hierarchy, children_list = children_stack.popleft()
            for child_dict in children_list:
                child_name = child_dict['name']
                if 'children' in child_dict:
                    children_stack.append((hierarchy + [child_name], child_dict['children']))
                else:
                    try:
                        categorization_dict[child_name].append(hierarchy)
                    except KeyError:
                        categorization_dict[child_name] = [hierarchy]

        return categorization_dict
