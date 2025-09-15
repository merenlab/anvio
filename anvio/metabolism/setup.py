import os
import re
import glob
import shutil

import anvio
import anvio.db as db
import anvio.utils as utils
import anvio.filesnpaths as filesnpaths
import anvio.terminal as terminal

from anvio.errors import ConfigError
from anvio.version import versions_for_db_types
from anvio.terminal import pluralize as P

from anvio.metabolism.context import KeggContext
from anvio.metabolism.modulesdb import ModulesDatabase


class KeggSetup(KeggContext):
    """Class for setting up KEGG Kofam HMM profiles and modules.

    It performs sanity checks and downloads, unpacks, and prepares the profiles for later use by `hmmscan`.
    It also downloads module files and creates the MODULES.db.

    Parameters
    ==========
    args: Namespace object
        All the arguments supplied by user to anvi-setup-kegg-data. If using this class through the API, please
        provide a Namespace object with the Boolean 'reset' parameter.
    skip_init: Boolean
        Developers can use this flag to skip the sanity checks and creation of directories when testing this class
    """

    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress(), skip_init=False):
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.args = args
        self.run = run
        self.progress = progress
        self.num_threads = 1 if not A('num_threads') else A('num_threads')
        self.kegg_archive_path = A('kegg_archive')
        self.kegg_snapshot = A('kegg_snapshot')
        self.download_from_kegg = True if A('download_from_kegg') else False
        self.only_download = True if A('only_download') else False
        self.only_processing = True if A('only_processing') else False
        self.skip_init = skip_init
        self.skip_brite_hierarchies = True if A('skip_brite_hierarchies') else False
        self.skip_binary_relations = True if A('skip_binary_relations') else False
        self.skip_map_images = True if A('skip_map_images') else False

        if self.kegg_archive_path and self.download_from_kegg:
            raise ConfigError("You provided two incompatible input options, --kegg-archive and --download-from-kegg. "
                              "Please pick either just one or none of these. ")
        if self.kegg_snapshot and self.download_from_kegg or self.kegg_snapshot and self.kegg_archive_path:
            raise ConfigError("You cannot request setup from an anvi'o KEGG snapshot at the same time as from KEGG directly or from one of your "
                              "KEGG archives. Please pick just one setup option and try again.")

        if (not self.download_from_kegg) and (self.only_download or self.only_processing):
            raise ConfigError("Erm. The --only-download and --only-processing options are only valid if you are also using the --download-from-kegg "
                              "option. Sorry.")
        if self.only_download and self.only_processing:
            raise ConfigError("The --only-download and --only-processing options are incompatible. Please choose only one. Or, if you want both "
                              "download AND database setup to happen, then use only the -D flag without providing either of these two options.")


        # initializing these to None here so that it doesn't break things downstream
        self.pathway_dict = None
        self.brite_dict = None

        # init the base class
        KeggContext.__init__(self, self.args)

        # get KEGG snapshot info for default setup
        self.target_snapshot = self.kegg_snapshot or 'v2025-08-07'
        self.target_snapshot_yaml = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/KEGG-SNAPSHOTS.yaml')
        self.snapshot_dict = utils.get_yaml_as_dict(self.target_snapshot_yaml)

        if self.target_snapshot not in self.snapshot_dict.keys():
            self.run.warning(None, header="AVAILABLE KEGG SNAPSHOTS", lc="yellow")
            available_snapshots = sorted(list(self.snapshot_dict.keys()))
            for snapshot_name in available_snapshots:
                self.run.info_single(f"{snapshot_name}\thash: {self.snapshot_dict[snapshot_name]['hash']}" + (' (latest)' if snapshot_name == available_snapshots[-1] else ''))

            raise ConfigError("Whoops. The KEGG snapshot you requested is not one that is known to anvi'o. Please try again, and "
                                "this time pick from the list shown above.")

        # default download path for KEGG snapshot
        self.default_kegg_data_url = self.snapshot_dict[self.target_snapshot]['url']
        self.default_kegg_archive_file = self.snapshot_dict[self.target_snapshot]['archive_name']

        # the KEGG API URL, in case its needed downstream
        self.kegg_rest_api_get = "http://rest.kegg.jp/get"

        if self.user_input_dir:
            self.run.warning(f"Just so you know, we will be setting up the metabolism data provided at the following "
                             f"location: '{self.user_input_dir}'. The success of this will be determined by how well you "
                             f"followed our formatting guidelines, so keep an eye out for errors below.")


        if not self.user_input_dir:

            # establish parent directory
            if self.download_from_kegg and not self.only_processing and not self.kegg_archive_path and not skip_init:
                filesnpaths.gen_output_directory(self.kegg_data_dir, delete_if_exists=args.reset)

        else: # user input setup
            filesnpaths.is_output_dir_writable(os.path.dirname(self.user_input_dir))

            self.check_user_input_dir_format()

            if not args.reset and not skip_init:
                self.is_user_database_exists()

            if args.reset:
                self.run.warning("Since you used the --reset flag, anvi'o will get rid of any existing user modules database. "
                                 "Now ye be warned.")
                paths_to_remove = [self.user_modules_db_path]
                for path in paths_to_remove:
                    if os.path.exists(path):
                        os.remove(path)
                        self.run.info("Successfully removed", path)


    def is_database_exists(self, files_to_check, fail_if_exists=True):
        """This function determines whether the user has already downloaded all required KEGG data.

        More specifically, it looks for the KEGG files that we use to learn what to download (as in
        the KEGG MODULE file) and for the existence of the data directories that are created by this
        program.

        PARAMETERS
        ==========
        files_to_check : list of file paths
            this list should contain the paths to all required KEGG data or directories. what those
            files are depends on the download mode.
        fail_if_exists : Boolean
            if this is True, this function will fail if the KEGG data already exists on the user's
            computer. If it is False, AND the user has already downloaded all required KEGG data,
            then this function will not fail. This is to enable the --only-processing option.
            Note that in this case we require all KEGG data to be pre-downloaded to avoid mixing
            older and newer KEGG data - so if this data is only partially downloaded, the function
            will raise an error even if this parameter is False.
        """

        if anvio.DEBUG:
            file_str = ", ".join(files_to_check)
            self.run.warning(f"We are looking for the following files to see if the KEGG data already "
                             f"exists on you computer: {file_str}")

        files_that_exist = []
        for f in files_to_check:
            if os.path.exists(f):
                if fail_if_exists:
                    raise ConfigError(f"It seems you already have data at {f}, please use the `--reset` flag "
                                      "or delete the KEGG data directory manually if you want to re-download KEGG data. "
                                      "See also the --only-processing option, which you can use if you already "
                                      "have all required KEGG data in that folder. (API users: skip this sanity "
                                      "check by initializing this class with `skip_init=True`)")
                else:
                    files_that_exist.append(f)

        if files_that_exist:
            exist_str = "\n".join(files_that_exist)
            # we require all data to be present. Otherwise we might produce chimeric KEGG data.
            if files_that_exist != files_to_check:
                raise ConfigError(f"We found some, but not all, required KEGG data on your computer in the KEGG "
                                  f"data directory. Since you don't have everything you need, we need you to re-download "
                                  f"everything from scratch. Please re-run this program using the --reset flag, and if "
                                  f"you were using the --only-processing option, remove that flag. :) HOWEVER, if you notice that "
                                  "KEGG BRITE data does not appear to be in the upcoming list, but you don't actually want "
                                  "to download BRITE data, then you can just add the --skip-brite-hierarchies to your previous "
                                  f"command and be on your way (ie, no --reset needed). Here is the KEGG data we found:\n{exist_str}")

            self.run.warning(f"We found already-downloaded KEGG data on your computer. Setup will continue using "
                             f"this data. However, if you think everything should be re-downloaded from scratch, please kill this program "
                             f"and restart it using the `--reset` flag. Here is the data you already have, in case you "
                             f"need to check it to make sure we are not using something that is too old:\n"
                             f"{exist_str}")

        if self.only_processing and not files_that_exist:
            raise ConfigError(f"We noticed that there is no KEGG data on your computer at {self.kegg_data_dir} even "
                              f"though you used the --only-processing option. If you don't actually have KEGG data already "
                              f"downloaded, you should get rid of the --only-processing flag and re-run this program. If you "
                              f"know that you DO have KEGG data, perhaps you gave us the wrong data directory?")


    def setup_from_archive(self):
        """This function sets up the KEGG data directory from an archive of a previously-setup KEGG data directory.

        To do so, it unpacks the archive and checks its structure and that all required components are there.
        """

        self.run.info("KEGG archive", self.kegg_archive_path)
        self.progress.new('Unzipping KEGG archive file...')
        if not self.kegg_archive_path.endswith("tar.gz"):
            self.progress.reset()
            raise ConfigError("The provided archive file %s does not appear to be an archive at all. Perhaps you passed "
                              "the wrong file to anvi'o?" % (self.kegg_archive_path))
        unpacked_archive_name = "KEGG_archive_unpacked"
        utils.tar_extract_file(self.kegg_archive_path, output_file_path=unpacked_archive_name, keep_original=True)

        self.progress.update('Checking KEGG archive structure and contents...')
        archive_is_ok = self.kegg_archive_is_ok(unpacked_archive_name)
        archive_contains_brite = self.check_archive_for_brite(unpacked_archive_name)
        archive_contains_binary_relations = self.check_archive_for_binary_relations(unpacked_archive_name)
        archive_contains_map_images = self.check_archive_for_map_files(unpacked_archive_name)
        self.progress.end()
        if archive_is_ok:
            if os.path.exists(self.kegg_data_dir):
                shutil.rmtree(self.kegg_data_dir)
            path_to_kegg_in_archive = os.path.join(unpacked_archive_name, "KEGG")
            shutil.move(path_to_kegg_in_archive, self.kegg_data_dir)
            shutil.rmtree(unpacked_archive_name)

            if not archive_contains_brite and not self.skip_brite_hierarchies:
                self.run.warning("The KEGG data archive does not contain the necessary files to set up BRITE hierarchy classification. "
                                 "This is not a problem, and KEGG set up proceeded without it. BRITE is guaranteed to be set up when "
                                 "downloading the latest version of KEGG with `anvi-setup-kegg-data`.")

            if not archive_contains_binary_relations and not self.skip_binary_relations:
                self.run.warning(
                    "The KEGG data archive does not contain the binary relation files needed for "
                    "`anvi-reaction-network`. This is not a problem, and KEGG setup proceeded "
                    "without it. Binary relation files are guaranteed to be set up when "
                    "downloading the latest version of KEGG with `anvi-setup-kegg-data`."
                )

            if not archive_contains_map_images and not self.skip_map_images:
                self.run.warning(
                    "The KEGG data archive does not contain the expected pathway map files used "
                    "for pathway visualization. This is not a problem, and KEGG setup proceeded "
                    "without it. Map files are guaranteed to be set up when downloading the latest "
                    "version of KEGG with `anvi-setup-kegg-data`."
                )

            # if necessary, warn user about migrating the modules db
            self.check_modules_db_version()

        else:
            debug_output = f"We kept the unpacked archive for you to take a look at it. It is at " \
                           f"{os.path.abspath(unpacked_archive_name)} and you may want " \
                           f"to delete it after you are done checking its contents."
            if not anvio.DEBUG:
                shutil.rmtree(unpacked_archive_name)
                debug_output = "The unpacked archive has been deleted, but you can re-run the script with the --debug " \
                               "flag to keep it if you want to see its contents."
            else:
                self.run.warning(f"The unpacked archive file {os.path.abspath(unpacked_archive_name)} was kept for "
                                 f"debugging purposes. You may want to clean it up after you are done looking through it.")

            raise ConfigError(f"SETUP FAILED. The provided archive file is missing some critical files, "
                              f"so anvi'o is unable to use it. {debug_output}")


    def check_modules_db_version(self):
        """This function checks if the MODULES.db is out of date and if so warns the user to migrate it"""

        # get current version of db
        db_conn = db.DB(self.kegg_modules_db_path, None, ignore_version=True)
        current_db_version = int(db_conn.get_meta_value('version'))
        db_conn.disconnect()

        # if modules.db is out of date, give warning
        target_version = int(versions_for_db_types['modules'])
        if current_db_version != target_version:
            self.run.warning(f"Just so you know, the KEGG archive that was just set up contains an outdated MODULES.db (version: "
                             f"{current_db_version}). You may want to run `anvi-migrate` on this database before you do anything else. "
                             f"Here is the path to the database: {self.kegg_modules_db_path}")


    def check_archive_for_brite(self, unpacked_archive_path):
        """Check the archive for the BRITE directory and 'hierarchy of hierarchies' json file.

        It is ok for archives not to have these present, but let the user know.
        """

        is_brite_included = True

        path_to_kegg_in_archive = os.path.join(unpacked_archive_path, "KEGG")
        brite_directories_and_files = [self.brite_data_dir,
                                       self.kegg_brite_hierarchies_file]
        for f in brite_directories_and_files:
            path_to_f_in_archive = os.path.join(path_to_kegg_in_archive, os.path.basename(f))
            if not os.path.exists(path_to_f_in_archive) and not self.skip_brite_hierarchies:
                is_brite_included = False
                if anvio.DEBUG:
                    self.run.warning(f"The KEGG archive does not contain the following optional BRITE file or directory: {path_to_f_in_archive}")

        return is_brite_included


    def check_archive_for_binary_relations(self, unpacked_archive_path):
        """
        Check the archive for the binary relations directory and files.

        It is ok for archives not to have these present, but let the user know.
        """
        path_to_kegg_in_archive = os.path.join(unpacked_archive_path, "KEGG")
        binary_relation_data_dir = os.path.join(
            path_to_kegg_in_archive, os.path.basename(self.binary_relation_data_dir)
        )
        if os.path.isdir(binary_relation_data_dir):
            is_binary_relation_dir_included = True
        else:
            is_binary_relation_dir_included = False
            if anvio.DEBUG and not self.skip_binary_relations:
                self.run.warning(
                    "The KEGG archive does not contain the following optional binary relations "
                    f"directory needed for `anvi-reaction-network`: {binary_relation_data_dir}"
                )

        if is_binary_relation_dir_included:
            missing_files = []
            for file in self.kegg_binary_relation_files.values():
                path = os.path.join(binary_relation_data_dir, file)
                if not os.path.isfile(path):
                    missing_files.append(file)
            if anvio.DEBUG and missing_files:
                self.run.warning(
                    "The following binary relation files expected in an up-to-date anvi'o KEGG "
                    f"installation are missing from the directory, '{binary_relation_data_dir}', "
                    f"in the archive: {', '.join(missing_files)}"
                )

        return is_binary_relation_dir_included


    def check_archive_for_map_files(self, unpacked_archive_path):
        """
        Check the archive for the pathway map directory containing image and KGML files, and for the
        BRITE json file classifying pathway maps.

        It is ok for archives not to have these present, but let the user know.
        """
        are_map_files_included = True
        path_to_kegg_in_archive = os.path.join(unpacked_archive_path, "KEGG")

        map_image_data_dir = os.path.join(
            path_to_kegg_in_archive, os.path.basename(self.map_image_data_dir)
        )
        if not os.path.isdir(map_image_data_dir):
            are_map_files_included = False
            if anvio.DEBUG and not self.skip_map_images:
                self.run.warning(
                    "The KEGG archive does not contain the following optional pathway map images "
                    f"directory, which is used in pathway visualization: {map_image_data_dir}"
                )

        brite_json_file = os.path.join(path_to_kegg_in_archive, os.path.basename(self.kegg_brite_pathways_file))
        if not os.path.isfile(brite_json_file):
            are_map_files_included = False
            if anvio.DEBUG and not self.skip_map_images:
                self.run.warning(
                    "The KEGG archive does not contain the following optional json file, a BRITE "
                    f"hierarchy classifying pathway maps: {self.kegg_brite_pathways_file}"
                )

        return are_map_files_included


    def setup_kegg_snapshot(self):
        """This is the default setup strategy in which we unpack a specific KEGG archive.

        We do this so that everyone who uses the same release of anvi'o will also have the same default KEGG
        data, which facilitates sharing and also means they do not have to continuously re-annotate their datasets
        when KEGG is updated.

        It is essentially a special case of setting up from an archive.
        """

        if anvio.DEBUG:
            self.run.info("Downloading from: ", self.default_kegg_data_url)
            self.run.info("Downloading to: ", self.default_kegg_archive_file)
        utils.download_file(self.default_kegg_data_url, self.default_kegg_archive_file, progress=self.progress, run=self.run)

        # a hack so we can use the archive setup function
        self.kegg_archive_path = self.default_kegg_archive_file
        self.setup_from_archive()

        # if all went well, let's get rid of the archive we used and the log file
        if not anvio.DEBUG:
            os.remove(self.default_kegg_archive_file)
        else:
            self.run.warning(f"Because you used the --debug flag, the KEGG archive file at {self.default_kegg_archive_file} "
                             "has been kept. You may want to remove it later.")


    def kegg_archive_is_ok(self, unpacked_archive_path):
        """This function checks the structure and contents of an unpacked KEGG archive and returns True if it is as expected.

        Please note that we check for existence of the files that are necessary to run KEGG scripts, but we don't check the file
        formats. This means that people could technically trick this function into returning True by putting a bunch of crappy files
        with the right names/paths into the archive file. But what would be the point of that?

        We also don't care about the contents of certain folders (ie modules) because they are not being directly used
        when running KEGG scripts. In the case of modules, all the information should already be in the MODULES.db so we don't
        waste our time checking that all the module files are there. We only check that the directory is there. If later changes
        to the implementation require the direct use of the files in these folders, then this function should be updated
        to check for those.

        Parameters
        ==========
        unpacked_archive_path : str
            Path to the unpacked archive directory
        """

        is_ok = True

        # check top-level files and folders
        path_to_kegg_in_archive = os.path.join(unpacked_archive_path, "KEGG")
        expected_directories_and_files = [self.orphan_data_dir,
                                          self.kegg_module_data_dir,
                                          self.kegg_hmm_data_dir,
                                          self.ko_list_file_path,
                                          self.kegg_module_file,
                                          self.kegg_modules_db_path]
        for f in expected_directories_and_files:
            path_to_f_in_archive = os.path.join(path_to_kegg_in_archive, os.path.basename(f))
            if not os.path.exists(path_to_f_in_archive):
                is_ok = False
                if anvio.DEBUG:
                    self.run.warning(f"The KEGG archive does not contain the following expected file or directory: "
                                     f"{path_to_f_in_archive}")

        # check hmm files
        path_to_hmms_in_archive = os.path.join(path_to_kegg_in_archive, os.path.basename(self.kegg_hmm_data_dir))
        kofam_hmm_basename = os.path.basename(self.kofam_hmm_file_path)
        expected_hmm_files = [kofam_hmm_basename]
        for h in expected_hmm_files:
            path_to_h_in_archive = os.path.join(path_to_hmms_in_archive, h)
            if not os.path.exists(path_to_h_in_archive):
                is_ok = False
                if anvio.DEBUG:
                    self.run.warning(f"The KEGG archive does not contain the following expected hmm file: "
                                     f"{path_to_h_in_archive}")
            expected_extensions = ['.h3f', '.h3i', '.h3m', '.h3p']
            for ext in expected_extensions:
                path_to_expected_hmmpress_file = path_to_h_in_archive + ext
                if not os.path.exists(path_to_expected_hmmpress_file):
                    is_ok = False
                    if anvio.DEBUG:
                        self.run.warning(f"The KEGG archive does not contain the following expected `hmmpress` output: "
                                         f"{path_to_expected_hmmpress_file}")

        return is_ok


    def setup_all_data_from_archive_or_snapshot(self):
        """This driver function controls whether we download one of our KEGG snapshots and set that up, or
        set up directly from an archive file already on the user's computer.
        """

        if os.path.exists(self.kegg_data_dir) and not self.args.reset:
            raise ConfigError(f"The directory {self.kegg_data_dir} already exists. Are you sure you want to "
                              f"overwrite it? If yes, feel free to restart this program with the --reset flag.")

        if self.kegg_archive_path:
            self.setup_from_archive()
        else:
            self.setup_kegg_snapshot()


    def check_user_input_dir_format(self):
        """This function checks whether the user input directory exists and contains the required subfolders

        The required subfolders are:
            modules : directory containing the user's metabolic pathway definitions (as text files)
        """

        for path in [self.user_input_dir, self.user_module_data_dir]:
            if not os.path.exists(path):
                raise ConfigError(f"There is a problem with the input directory you provided. The following path does not "
                                  f"exist: '{path}'. Please make sure that your input folder exists and that it follows the "
                                  f"formatting requirements. We're sorry for asking this of you, but it really helps us make "
                                  f"sure everything will go smoothly.")

            file_list = [f for f in glob.glob(os.path.join(path, '*'))]
            if not file_list:
                raise ConfigError(f"The folder '{path}' appears to be empty, so we have no data to work with. Please make "
                                  f"sure that you have provided the correct input directory and formatted it correctly so "
                                  f"that anvi'o can find your data.")


    def is_user_database_exists(self):
        """This function checks whether user data has already been set up in the provided input directory."""

        if os.path.exists(self.user_modules_db_path):
            raise ConfigError(f"It seems you already have a user modules database installed at '{self.user_modules_db_path}', "
                              f"please use the --reset flag or delete this file manually if you want to re-generate it.")


    def process_pathway_file(self):
        """This function reads the kegg pathway map file into a dictionary. It should be called during setup to get the KEGG pathway ids so the pathways can be downloaded.

        The structure of this file is like this:

        +C	Map number
        #<h2><a href="/kegg/kegg2.html"><img src="/Fig/bget/kegg3.gif" align="middle" border=0></a>&nbsp; KEGG Pathway Maps</h2>
        !
        A<b>Metabolism</b>
        B  Global and overview maps
        C    01100  Metabolic pathways
        C    01110  Biosynthesis of secondary metabolites
        C    01120  Microbial metabolism in diverse environments
        C    01200  Carbon metabolism
        C    01210  2-Oxocarboxylic acid metabolism

        Initial lines can be ignored and thereafter the line's information can be determined by the one-letter code at the start.
        A = Category of Pathway Map
        B = Sub-category of Pathway Map
        C = Pathway Map identifier number and name

        Note that not all Pathway Maps that we download will have ORTHOLOGY fields. We don't exclude these here, but processing later
        will have to be aware of the fact that not all pathways will have associated KOs.

        We do, however, exclude Pathway Maps that don't have existing `koXXXXX` identifiers (these yield 404 errors when attempting to
        download them). For instance, we exclude those that start with the code 010 (chemical structure maps) or with 07 (drug structure maps).
        """

        self.pathway_dict = {}

        filesnpaths.is_file_exists(self.kegg_pathway_file)
        filesnpaths.is_file_plain_text(self.kegg_pathway_file)

        f = open(self.kegg_pathway_file, 'r')
        self.progress.new("Parsing KEGG Pathway file")

        current_category = None
        current_subcategory = None


        for line in f.readlines():
            line = line.strip('\n')
            first_char = line[0]

            # garbage lines
            if first_char in ["+", "#", "!"]:
                continue
            else:
                # Category
                if first_char == "A":
                    fields = re.split('<[^>]*>', line) # we split by the html tag here
                    current_category = fields[1]
                # Sub-category
                elif first_char == "B":
                    fields = re.split('\s{2,}', line) # don't want to split the subcategory name, so we have to split at least 2 spaces
                    current_subcategory = fields[1]
                elif first_char == "C":
                    fields = re.split('\s{2,}', line)
                    konum = "ko" + fields[1]
                    if konum[:5] != "ko010" and konum[:4] != "ko07":
                        self.pathway_dict[konum] = {"name" : fields[2], "category" : current_category, "subcategory" : current_subcategory}
                # unknown code
                else:
                    raise ConfigError("While parsing the KEGG file %s, we found an unknown line code %s. This has "
                                      "made the file unparseable. It is likely that an update to KEGG has broken "
                                      "things such that anvi'o doesn't know what is going on anymore. Sad, we know. :( "
                                      "Please contact the developers to see if this is a fixable issue, and in the "
                                      "meantime use an older version of the KEGG data directory (if you have one). "
                                      "If we cannot fix it, we may be able to provide you with a legacy KEGG "
                                      "data archive that you can use to setup KEGG with the --kegg-archive flag." % (self.kegg_pathway_file, first_char))
        self.progress.end()


    def get_accessions_from_htext_file(self, htext_file):
        """This function can read generic KEGG htext files to get a list of accessions.

        Here is one example of the file structure, taken from a COMPOUND htext file:
        +D	Biochemical compound
        #<h2><a href="/kegg/brite.html"><img src="/Fig/bget/kegg3.gif" align="middle" border=0></a>&nbsp; Compounds with Biological Roles</h2>
        !
        A<b>Organic acids</b>
        B  Carboxylic acids [Fig]
        C    Monocarboxylic acids
        D      C00058  Formate; Methanoate
        D      C00033  Acetate; Ethanoate
        D      C00163  Propionate; Propanoate
        D      C00246  Butyrate; Butanoate

        The +(letter) at the start of the first line indicates how many levels the hierarchy contains (often C or D). For the purpose of
        downloading other files, we only need the accessions from the lowest level. All other information is skipped when reading the file.

        PARAMETERS
        ==========
        htext_file : str
            The filename of the hierachical text file downloaded from KEGG

        RETURNS
        =======
        accession_list : list of str
            Contains the KEGG identifiers contained in the lowest hierarchy of this file.
        """

        accession_list = []

        filesnpaths.is_file_exists(htext_file)
        filesnpaths.is_file_plain_text(htext_file)

        f = open(htext_file, 'r')
        self.progress.new(f"Parsing KEGG htext file: {htext_file}")

        target_level = None
        for line in f.readlines():
            line = line.strip('\n')
            first_char = line[0]

            if first_char == '+':  # first line of the file; second character gives us target level
                target_level = line[1]
            elif first_char == target_level: # need to extract the accession, which is second field (split on spaces)
                fields = re.split('\s+', line)
                accession_list.append(fields[1])
            else: # skip everything else
                continue

        self.progress.end()

        num_acc = len(accession_list)
        self.run.info("Number of accessions found in htext file", num_acc)
        return accession_list


    def download_generic_htext(self, h_accession, download_dir="./"):
        """Downloads the KEGG htext file for the provided accession.

        PARAMETERS
        ==========
        h_accession : str
            The accession for a KEGG hierarchy file
        download_dir : str
            Path to directory where file will be downloaded. Current working directory by default.
        """

        htext_url_prefix = "https://www.genome.jp/kegg-bin/download_htext?htext="
        htext_url_suffix = ".keg&format=htext&filedir="
        htext_url = htext_url_prefix+h_accession+htext_url_suffix

        htext_file = h_accession + ".keg"
        path_to_download_to = os.path.join(download_dir,htext_file)

        if filesnpaths.is_file_exists(path_to_download_to, dont_raise=True):
            if not self.args.reset:
                raise ConfigError(f"The file at {path_to_download_to} already exists. If you are "
                                  f"sure that you want to download it, you can avoid this error message "
                                  f"by using the 'reset' parameter. Make sure that won't erase your other "
                                  f"KEGG data, though.")

        try:
            utils.download_file(htext_url, path_to_download_to, progress=self.progress, run=self.run)
        except Exception as e:
            print(e)
            raise ConfigError(f"Anvi'o failed to download the KEGG htext file for {h_accession} from {htext_url}.")

        return path_to_download_to


    def download_generic_flat_file(self, accession, download_dir="./"):
        """Downloads the flat file for the given accession from the KEGG API.

        PARAMETERS
        ==========
        accession : str
            A KEGG identifier
        download_dir : str
            Path to the directory in which to download the file. Current working directory by default.
        """

        file_path = os.path.join(download_dir, accession)
        if filesnpaths.is_file_exists(file_path, dont_raise=True):
            if not self.args.reset:
                raise ConfigError(f"The file at {file_path} already exists. If you are "
                                  f"sure that you want to download it, you can avoid this error message "
                                  f"by using the 'reset' parameter. Make sure that won't erase your other "
                                  f"KEGG data, though.")

        utils.download_file(self.kegg_rest_api_get + '/' + accession,
            file_path, progress=self.progress, run=self.run)
        # verify entire file has been downloaded
        f = open(file_path, 'r')
        f.seek(0, os.SEEK_END)
        f.seek(f.tell() - 4, os.SEEK_SET)
        last_line = f.readline().strip('\n')
        if not last_line == '///':
            raise ConfigError(f"The KEGG flat file {file_path} was not downloaded properly. We were expecting the last line in the file "
                              f"to be '///', but instead it was {last_line}. Formatting of these files may have changed on the KEGG website. "
                              f"Please contact the developers to see if this is a fixable issue.")


    def download_kegg_files_from_hierarchy(self, h_accession, download_dir="./"):
        """Given the accession of a KEGG hierarchy, this function downloads all of its flat files.

        PARAMETERS
        ==========
        h_accession : str
            The accession for a KEGG hierarchy file
        download_dir : str
            Path to the directory in which to download the files. Current working directory by default.
            (a folder to store the hierarchy's flat files will be generated in this folder)
        """

        filesnpaths.is_output_dir_writable(download_dir)

        htext_filename = self.download_generic_htext(h_accession, download_dir)
        acc_list = self.get_accessions_from_htext_file(htext_filename)

        download_dir_name = os.path.join(download_dir,h_accession)
        filesnpaths.gen_output_directory(download_dir_name, delete_if_exists=self.args.reset)

        self.run.info("KEGG Module Database URL", self.kegg_rest_api_get)
        self.run.info("Number of KEGG files to download", len(acc_list))
        self.run.info("Directory to store files", download_dir_name)

        # download all modules
        for acc in acc_list:
            self.download_generic_flat_file(acc, download_dir_name)


    def download_pathways(self):
        """This function downloads the KEGG Pathways.

        To do so, it first processes a KEGG file containing pathway and map identifiers into a dictionary via the process_pathway_file()
        function. To verify that each file has been downloaded properly, we check that the last line is '///'.
        """

        # note that this is the same as the REST API for modules - perhaps at some point this should be printed elsewhere so we don't repeat ourselves.
        self.run.info("KEGG Pathway Database URL", self.kegg_rest_api_get)

        # download the kegg pathway file, which lists all modules
        try:
            utils.download_file(self.kegg_pathway_download_path, self.kegg_pathway_file, progress=self.progress, run=self.run)
        except Exception as e:
            print(e)
            raise ConfigError("Anvi'o failed to download the KEGG Pathway htext file from the KEGG website. Something "
                              "likely changed on the KEGG end. Please contact the developers to see if this is "
                              "a fixable issue. If it isn't, we may be able to provide you with a legacy KEGG "
                              "data archive that you can use to setup KEGG with the --kegg-archive flag.")

        # get pathway dict
        self.process_pathway_file()
        self.run.info("Number of KEGG Pathways", len(self.pathway_dict.keys()))

        # download all pathways
        for konum in self.pathway_dict.keys():
            file_path = os.path.join(self.pathway_data_dir, konum)
            utils.download_file(self.kegg_rest_api_get + '/' + konum,
                file_path, progress=self.progress, run=self.run)
            # verify entire file has been downloaded
            f = open(file_path, 'r')
            f.seek(0, os.SEEK_END)
            f.seek(f.tell() - 4, os.SEEK_SET)
            last_line = f.readline().strip('\n')
            if not last_line == '///':
                raise ConfigError("The KEGG pathway file %s was not downloaded properly. We were expecting the last line in the file "
                                  "to be '///', but instead it was %s. Formatting of these files may have changed on the KEGG website. "
                                  "Please contact the developers to see if this is a fixable issue. If it isn't, we may be able to "
                                  "provide you with a legacy KEGG data archive that you can use to setup KEGG with the --kegg-archive flag."
                                  % (file_path, last_line))


    def extract_data_field_from_kegg_file(self, file_path, target_field):
        """This function parses a KEGG file and returns the data value associated with the given target field.

        It can work on flat-text files obtained via the REST API (ie, self.kegg_rest_api_get).
        """

        data_to_return = []

        f = open(file_path, 'r')
        current_data_name = None

        for line in f.readlines():
            line = line.strip('\n')

            fields = re.split('\s{2,}', line)
            data_vals = None

            # when data name unknown, parse from first field
            if line[0] != ' ':
                current_data_name = fields[0]
            if line[0] == ' ' and not current_data_name:
                raise ConfigError(f"Uh oh. While trying to parse the KEGG file at {file_path}, we couldn't "
                "find the data field associated with the line '{line}'.")

            # note that if data name is known, first field still exists but is actually the empty string ''
            if len(fields) > 1:
                data_vals = fields[1]

            if (current_data_name == target_field) and (data_vals is not None):
                data_to_return.append(data_vals)

        f.close()

        return data_to_return


    def create_user_modules_dict(self):
        """This function establishes the self.module_dict parameter for user modules.

        It is essentially a replacement for the process_module_file() function.
        Since users will not have a modules file to process, we simply create the dictionary from the
        file names they provide for their module definitions. We don't add any dictionary values,
        but we won't need them (we hope).
        """

        user_module_list = [os.path.basename(k) for k in glob.glob(os.path.join(self.user_module_data_dir, '*'))]
        self.module_dict = {key: {} for key in user_module_list}

        # sanity check that they also have KEGG data since we need to compare module names
        if not os.path.exists(self.kegg_modules_db_path):
            raise ConfigError(f"Wait a second. We understand that you are setting up user-defined metabolism data, but "
                              f"unfortunately you need to FIRST have KEGG data set up on your computer. Why, you ask? "
                              f"Well, we need to make sure none of your module names overlap with those "
                              f"in the KEGG MODULES database. Long story short, we looked for KEGG data at "
                              f"{self.kegg_modules_db_path} but we couldn't find it. If this is the wrong place for us to be "
                              f"looking, please run this program again and use the --kegg-data-dir parameter to tell us where "
                              f"to find it.")

        # sanity check that user module names are distinct
        kegg_modules_db = ModulesDatabase(self.kegg_modules_db_path, args=self.args, quiet=True)
        kegg_mods = set(kegg_modules_db.get_all_modules_as_list())
        user_mods = set(user_module_list)
        bad_user_mods = kegg_mods.intersection(user_mods)
        if bad_user_mods:
            bad_mods_str = ", ".join(bad_user_mods)
            n = len(bad_user_mods)
            raise ConfigError(f"Hol'up a minute. You see, there {P('is a module', n, alt='are some modules')} "
                              f"in your user-defined modules data (at {self.user_module_data_dir}) which {P('has', n, alt='have')} "
                              f"the same name as an existing KEGG module. This is not allowed, for reasons. Please name {P('that module', n, alt='those modules')} "
                              f"differently. Append an underscore and your best friend's name to {P('it', n, alt='them')} or something. Just make sure it's "
                              f"unique. OK? ok. Here is the list of module names you should change: {bad_mods_str}")


    def setup_modules_db(self, db_path, module_data_directory, brite_data_directory=None, source='KEGG', skip_brite_hierarchies=False):
        """This function creates a Modules DB at the specified path."""

        if filesnpaths.is_file_exists(db_path, dont_raise=True):
            if self.overwrite_modules_db:
                os.remove(db_path)
            else:
                raise ConfigError(f"Woah there. There is already a modules database at {db_path}. If you really want to make a new modules database "
                                  f"in this folder, you should either delete the existing database yourself, or re-run this program with the "
                                  f"--overwrite-output-destinations flag. But the old database will go away forever in that case. Just making "
                                  f"sure you are aware of that, so that you have no regrets.")
        try:
            mod_db = ModulesDatabase(db_path, module_data_directory=module_data_directory, brite_data_directory=brite_data_directory, data_source=source,
                                     args=self.args, module_dictionary=self.module_dict, pathway_dictionary=self.pathway_dict, brite_dictionary=self.brite_dict,
                                     skip_brite_hierarchies=skip_brite_hierarchies, run=self.run, progress=self.progress)
            mod_db.create()
        except Exception as e:
            print(e)
            raise ConfigError("While attempting to build the MODULES.db, anvi'o encountered an error, which should be printed above. "
                              "If you look at that error and it seems like something you cannot handle, please contact the developers "
                              "for assistance. :) ")


    def setup_user_data(self):
        """This function sets up user metabolism data from the provided input directory.

        It processes the user's module files into the USER_MODULES.db.
        """

        self.create_user_modules_dict()
        self.setup_modules_db(db_path=self.user_modules_db_path, module_data_directory=self.user_module_data_dir, source='USER', skip_brite_hierarchies=True)
