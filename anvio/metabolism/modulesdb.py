import re
import os
import copy
import json
import time
import hashlib

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.errors import ConfigError
from anvio.terminal import pluralize as P

from anvio.metabolism.context import KeggContext


class ModulesDatabase(KeggContext):
    """To create or access a Modules DB.

    This DB should be created in the Kegg Data folder during KEGG setup, and will be populated with information from the
    Kegg Module files.

    If you want to load an existing database from the python terminal, all you need is the path to the database and an
    empty args object:
    ```
    >>> import argparse
    >>> from anvio.metabolism.modulesdb import ModulesDatabase
    >>> path_to_db = "YOUR/PATH/HERE/MODULES.db"
    >>> args = argparse.Namespace()
    >>> ModulesDatabase(path_to_db, args)
    ```
    """

    def __init__(self, db_path, args, module_data_directory=None, module_dictionary=None, pathway_dictionary=None, data_source='KEGG', run=terminal.Run(), progress=terminal.Progress(), quiet=False):
        self.db = None
        self.db_path = db_path
        self.module_data_directory = module_data_directory # only required for create()
        self.module_dict = module_dictionary
        self.pathway_dict = pathway_dictionary
        self.data_source = data_source
        self.run = run
        self.progress = progress
        self.quiet = quiet

        # keep track of functional annotation sources needed for the modules in this db
        self.annotation_sources = set()
        if self.data_source == 'KEGG':
            self.annotation_sources.add('KOfam')

        if anvio.DEBUG:
            self.run.info("Modules DB quiet param", self.quiet)

        # init the base class for access to shared paths and such
        KeggContext.__init__(self, args)

        # modules table info
        self.module_table_name = t.module_table_name
        self.module_table_structure = t.module_table_structure
        self.module_table_types = t.module_table_types

        # pathway maps table info
        self.pathway_table_name = t.pathway_table_name
        self.pathway_table_structure = t.pathway_table_structure
        self.pathway_table_types = t.pathway_table_types

        if os.path.exists(self.db_path):
            utils.is_kegg_modules_db(self.db_path)
            self.db = db.DB(self.db_path, anvio.__kegg_modules_version__, new_database=False)

            if not self.quiet:
                self.run.info("Modules database", f"An existing database, {self.db_path}, has been loaded.", quiet=self.quiet)
                self.run.info("Modules", f"{self.db.get_meta_value('num_modules')} found", quiet=self.quiet)
            
        else:
            # if self.module_dict is None, then we tried to initialize the DB outside of setup
            if not self.module_dict:
                raise ConfigError("ERROR - a new ModulesDatabase() cannot be initialized without providing a modules dictionary. This "
                                  "usually happens when you try to access a Modules DB before one has been setup. Running `anvi-setup-kegg-data` may fix this.")

######### DB GENERATION FUNCTIONS #########

    def touch(self):
        """Creates an empty Modules database on disk, and sets `self.db` to access to it.

        At some point self.db.disconnect() must be called to complete the creation of the new db.
        """

        # sanity check to avoid overriding previous Modules DB
        # this will probably never happen as long as this function is called through the setup script, but we check just in case
        if os.path.exists(self.db_path):
            raise ConfigError("A modules database at %s already exists. Please use the --reset flag when you restart the setup "
                              "if you really want to get rid of this one and make a new one." % (self.db_path))


        self.db = db.DB(self.db_path, anvio.__kegg_modules_version__, new_database=True)

        self.db.create_table(self.module_table_name, self.module_table_structure, self.module_table_types)


    def data_vals_sanity_check(self, data_vals, current_data_name, current_module_num):
        """This function checks if the data values were correctly parsed from a line in a KEGG module file.

        This is a sadly necessary step because some KEGG module file lines are problematic and don't follow the right format (ie, 2+ spaces
        between different fields). So here we check if the values that we parsed look like they are the right format, without any extra bits.
        Each data name (ORTHOLOGY, DEFINITION, etc) has a different format to check for.

        Note that we don't check the following data name types: NAME, CLASS, REFERENCE

        WARNING: The error checking and correction is by no means perfect and may well fail when KEGG is next updated. :(

        PARAMETERS
        ==========
        data_vals : str
            the data values field (split from the kegg module line)
        current_data_name : str
            which data name we are working on. It should never be None because we should have already figured this out by parsing the line.
        current_module_num : str
            which module we are working on. We need this to keep track of which modules throw parsing errors.

        RETURNS
        =======
        is_ok : bool
            whether the values look correctly formatted or not
        """

        is_ok = True
        is_corrected = False
        corrected_vals = None
        corrected_def = None

        if not current_data_name:
            raise ConfigError("data_vals_sanity_check() cannot be performed when the current data name is None. Something was not right "
                              "when parsing the KEGG module line.")
        elif current_data_name == "ENTRY":
            # example format: M00175
            if data_vals[0] != 'M' or len(data_vals) != 6:
                is_ok = False
                self.parsing_error_dict['bad_kegg_code_format'].append(current_module_num)
        elif current_data_name == "DEFINITION":
            # example format: (K01647,K05942) (K01681,K01682) (K00031,K00030) (K00164+K00658+K00382,K00174+K00175-K00177-K00176)
            # another example: (M00161,M00163) M00165
            knums = [x for x in re.split('\(|\)|,| |\+|-',data_vals) if x]
            for k in knums:
                if k[0] not in ['K','M'] or len(k) != 6:
                    is_ok = False
            if not is_ok: # this goes here to avoid counting multiple errors for the same line
                self.parsing_error_dict['bad_kegg_code_format'].append(current_module_num)
        elif current_data_name == "ORTHOLOGY":
            # example format: K00234,K00235,K00236,K00237
            # more complex example: (K00163,K00161+K00162)+K00627+K00382-K13997
            # another example:  (M00161         [ie, from (M00161  Photosystem II)]
            knums = [x for x in re.split('\(|\)|,|\+|-', data_vals) if x]
            for k in knums:
                if k[0] not in ['K','M'] or len(k) != 6:
                    is_ok = False
            # try to fix it by splitting on first space
            if not is_ok:
                self.parsing_error_dict['bad_line_splitting'].append(current_module_num)
                split_data_vals = data_vals.split(" ", maxsplit=1)
                corrected_vals = split_data_vals[0]
                corrected_def = split_data_vals[1]
                # double check that we don't have a knum in the new definition
                if re.match("K\d{5}",corrected_def):
                    corrected_vals = "".join([corrected_vals,corrected_def])
                    corrected_def = None
                is_corrected = True
        elif current_data_name == "PATHWAY":
            # example format: map00020
            if data_vals[0:3] != "map" or len(data_vals) != 8:
                is_ok = False
                self.parsing_error_dict['bad_line_splitting'].append(current_module_num)
                split_data_vals = data_vals.split(" ", maxsplit=1)
                corrected_vals = split_data_vals[0]
                corrected_def = split_data_vals[1]
                is_corrected = True
        elif current_data_name == "REACTION":
            # example format: R01899+R00268,R00267,R00709
            rnums = [x for x in re.split(',|\+', data_vals) if x]
            for r in rnums:
                if r[0] != 'R' or len(r) != 6:
                    is_ok = False
            if not is_ok:
                self.parsing_error_dict['bad_line_splitting'].append(current_module_num)
                split_data_vals = data_vals.split(" ", maxsplit=1)
                corrected_vals = split_data_vals[0]
                corrected_def = split_data_vals[1]
                is_corrected = True
        elif current_data_name == "COMPOUND":
            # example format: C00024
            if data_vals[0] not in ['C','G'] or len(data_vals) != 6:
                is_ok = False
                self.parsing_error_dict['bad_kegg_code_format'].append(current_module_num)
        elif current_data_name == "RMODULE":
            # example format: RM003
            if data_vals[0:2] != "RM" or len(data_vals) != 5:
                is_ok = False
                self.parsing_error_dict['bad_kegg_code_format'].append(current_module_num)


        if not is_ok and not is_corrected:
            self.num_uncorrected_errors += 1
            if self.just_do_it:
                self.progress.reset()
                self.run.warning("While parsing, anvi'o found an uncorrectable issue with a KEGG Module line in module %s, but "
                                 "since you used the --just-do-it flag, anvi'o will quietly ignore this issue and add the line "
                                 "to the MODULES.db anyway. Please be warned that this may break things downstream. In case you "
                                 "are interested, the line causing this issue has data name %s and data value %s."
                                 % (current_module_num, current_data_name, data_vals))
                is_ok = True # let's pretend that everything is alright so that the next function will take the original parsed values

            else:
                raise ConfigError("While parsing, anvi'o found an uncorrectable issue with a KEGG Module line in module %s. The "
                                  "current data name is %s, here is the incorrectly-formatted data value field: %s. If you think "
                                  "this is totally fine and want to ignore errors like this, please re-run the setup with the "
                                  "--just-do-it flag. But if you choose to do that of course we are obliged to inform you that things "
                                  "may eventually break as a result." % (current_module_num, current_data_name, data_vals))

        if is_corrected:
            self.num_corrected_errors += 1
            if anvio.DEBUG and not self.quiet:
                self.progress.reset()
                self.run.warning("While parsing a KEGG Module line, we found an issue with the formatting. We did our very best to parse "
                                 "the line correctly, but please check that it looks right to you by examining the following values.")
                self.run.info("Incorrectly parsed data value field", data_vals)
                self.run.info("Corrected data values", corrected_vals)
                self.run.info("Corrected data definition", corrected_def)

        return is_ok, corrected_vals, corrected_def


    def parse_kegg_modules_line(self, line, current_module, line_num=None, current_data_name=None, error_dictionary=None):
        """This function parses information from one line of a KEGG module file.

        These files have fields separated by 2 or more spaces. Fields can include data name (not always), data value (always), and data definition (not always).
        Lines for pathway module files can have between 1 and 4 fields, but in fact the only situation where there should be 4 lines is the ENTRY data,
        which for some inexplicable reason has multiple spaces between "Pathway" and "Module" in the data definition field. We can safely ignore this last "Module", I think.

        Some lines will have multiple entities in the data_value field (ie, multiple KOs or reaction numbers) and will be split into multiple db entries.

        PARAMETERS
        ==========
        line : str
            the line to parse
        current_module : str
            which module we are working on. We need this to keep track of which modules throw parsing errors
        line_num : int
            which line number we are working on. We need this to keep track of which entities come from the same line of the file.
        current_data_name : str
            which data name we are working on. If this is None, we need to parse this info from the first field in the line.

        RETURNS
        =======
        line_entries : list
            tuples, each containing information for one db entry, namely data name, data value, data definition, and line number.
            Not all parts of the db entry will be included (module num, for instance), so this information must be parsed and combined with
            the missing information before being added to the database.
        """

        if anvio.DEBUG:
            self.progress.reset()
            self.run.info("[DEBUG] Parsing line", line, mc='red', lc='yellow')
        fields = re.split('\s{2,}', line)
        data_vals = None
        data_def = None
        line_entries = []

        # when data name unknown, parse from first field
        if not current_data_name:
            # sanity check: if line starts with space then there is no data name field and we should have passed a current_data_name
            if line[0] == ' ':
                raise ConfigError("Oh, please. Some silly developer (you know who you are) has tried to call parse_kegg_modules_line() on "
                                  "a line without a data name field, and forgot to give it the current data name. Shame on you, go fix "
                                  "this. (For reference here is the line: %s)" % (line))

            current_data_name = fields[0]
        # note that if data name is known, first field still exists but is actually the empty string ''
        # so no matter the situation, data value is field 1 (index 0) and data definition (if any) is field 2 (index 1)
        # the only exception is that sometimes there is nothing in the data definition field (REFERENCE lines sometimes do this)
        if len(fields) > 1:
            data_vals = fields[1]

            if self.data_source == 'KEGG':
                # need to sanity check data value field because SOME modules don't follow the 2-space separation formatting
                vals_are_okay, corrected_vals, corrected_def = self.data_vals_sanity_check(data_vals, current_data_name, current_module)
            else:
                # TODO: USER data sanity check
                vals_are_okay = True # for now we always assume USER data is formatted properly

            if vals_are_okay and len(fields) > 2: # not all lines have a definition field
                data_def = fields[2]
            elif not vals_are_okay:
                data_vals = corrected_vals
                data_def = corrected_def
        else: # only the data name was in the line
            # these are the data types that we don't care if they have an empty line
            data_types_can_be_empty = ['REFERENCE', 'AUTHORS', 'TITLE', 'JOURNAL']
            if current_data_name in data_types_can_be_empty or self.just_do_it:
                if anvio.DEBUG:
                    self.run.warning(f"While parsing module {current_module} we found an empty {current_data_name} line. "
                                     "We think it is okay and it probably won't cause issues downstream.",
                                     header="DEBUG OUTPUT", lc='yellow')
            else:
                raise ConfigError(f"While parsing module {current_module} we found an empty {current_data_name} line. "
                                  "We are quitting here so you can check it, because this data type might be important. "
                                  "However, if you disagree, you can re-run the setup with --just-do-it and we will quietly "
                                  "incorporate this empty line into the MODULES.db (you may also need the --reset flag when you re-run). ")

        # some types of information may need to be split into multiple db entries
        data_types_to_split = ["ORTHOLOGY","REACTION"] # lines that fall under these categories need to have data_vals split on comma
        if current_data_name in data_types_to_split:
            # here we should NOT split on any commas within parentheses
            vals = [x for x in re.split('\(|\)|,|\+|-', data_vals) if x]
            for val in vals:
                line_entries.append((current_data_name, val, data_def, line_num))
        else:
            line_entries.append((current_data_name, data_vals, data_def, line_num))

        return line_entries


    def create(self):
        """Creates the Modules DB"""

        if not self.module_data_directory:
            raise ConfigError("Some dumb programmer forgot to provide a module_data_directory parameter value to the ModulesDatabase "
                              "class. The DB can't be created unless it knows where the modules are... Get yourself together.")

        self.touch()

        self.progress.new("Loading %s modules into Modules DB..." % len(self.module_dict.keys()))

        # sanity check that we setup the modules previously.
        # It shouldn't be a problem since this function should only be called during the setup process after modules download, but just in case.
        if not os.path.exists(self.module_data_directory) or len(self.module_dict.keys()) == 0:
            raise ConfigError("Appparently, the module data files were not correctly setup and now all sorts of things are broken. The "
                              "Modules DB cannot be created from broken things. BTW, this error is not supposed to happen to anyone "
                              "except maybe developers, so if you do not fall into that category you are likely in deep doo-doo. "
                              "Maybe re-running setup with --reset will work? (if not, you probably should email/Discord/telepathically "
                              "cry out for help to the developers). Anyway, if this helps make things any clearer, the number of modules "
                              "in the module dictionary is currently %s" % len(self.module_dict.keys()))

        # init the Modules table
        mod_table = ModulesTable(self.module_table_name)

        # keep track of errors encountered while parsing
        self.parsing_error_dict = {"bad_line_splitting" : [], "bad_kegg_code_format" : []}
        self.num_corrected_errors = 0
        self.num_uncorrected_errors = 0

        num_modules_parsed = 0
        line_number = 0
        for mnum in self.module_dict.keys():
            self.progress.update("Parsing Module %s" % mnum)
            mod_file_path = os.path.join(self.module_data_directory, mnum)
            f = open(mod_file_path, 'r')

            prev_data_name_field = None
            module_has_annotation_source = False
            orthology_to_annotation_source = {}  # for sanity check that each enzyme has an annotation source
            # for sanity check that each enzyme in definition has an orthology line
            mod_definition = []
            orth_list = set()
            for line in f.readlines():
                line = line.strip('\n')
                line_number += 1

                # check for last line ///. We don't want to send the last line to the parsing function because it will break.
                # we also check here that the line is not entirely blank (this happens sometimes in KEGG modules, inexplicably)
                if not line == '///' and re.search(r"\S+", line):
                    # parse the line into a tuple
                    entries_tuple_list = None
                    # here is the tricky bit about parsing these files. Not all lines start with the data_name field; those that don't start with a space.
                    # if this is the case, we need to tell the parsing function what the previous data_name field has been.
                    if line[0] == ' ':
                        entries_tuple_list = self.parse_kegg_modules_line(line, mnum, line_number, prev_data_name_field)
                    else:
                        entries_tuple_list = self.parse_kegg_modules_line(line, mnum, line_number)

                    prev_data_name_field = entries_tuple_list[0][0]

                    for name, val, definition, line in entries_tuple_list:
                        # there is one situation in which we want to ignore the entry, and that is Modules appearing in the ORTHOLOGY category, like so:
                        # (M00531  Assimilatory nitrate reduction, nitrate => ammonia)
                        if not (name == "ORTHOLOGY" and val[0] == '('):
                            # append_and_store will collect db entries and store every 10000 at a time
                            mod_table.append_and_store(self.db, mnum, name, val, definition, line)
                        else:
                            line -= 1

                        # save definition for later sanity check
                        if name == 'DEFINITION':
                            mod_definition.append(val)
                        if name == 'ORTHOLOGY':
                            orth_list.add(val)
                        # keep track of distinct annotation sources for user modules
                        if self.data_source != 'KEGG' and name == "ORTHOLOGY":
                            orthology_to_annotation_source[val] = None
                        if self.data_source != 'KEGG' and name == "ANNOTATION_SOURCE":
                            self.annotation_sources.add(definition)
                            module_has_annotation_source = True
                            if val not in orthology_to_annotation_source:
                                raise ConfigError(f"Woah. While parsing module {mnum} we found an ANNOTATION_SOURCE for "
                                                  f"an enzyme, {val}, that does not have an ORTHOLOGY line. Please check the "
                                                  f"module file and make sure that 1) ORTHOLOGY comes before ANNOTATION_SOURCE, "
                                                  f"and 2) each enzyme with an ORTHOLOGY also has a corresponding ANNOTATION_SOURCE.")
                            orthology_to_annotation_source[val] = definition

                        # sanity check for user-defined CLASS value
                        if self.data_source != 'KEGG' and name == "CLASS":
                            if len(val.split(";")) != 3:
                                raise ConfigError(f"The module {mnum} appears to have an invalid CLASS value. That value should be "
                                                  f"a string with a class, category, and sub-category separated by semi-colons (for a "
                                                  f"total of two semi-colons in the string). Instead, it is this: {val}")

                f.close()

            # every enzyme in the module definition needs an orthology line
            mod_definition = " ".join(mod_definition)
            # anything that is not (),-+ should be converted to spaces, then we can split on the spaces to get the accessions
            mod_definition = re.sub('[\(\)\+\-,]', ' ', mod_definition).strip()
            acc_list = re.split(r'\s+', mod_definition)
            accessions_in_def = set(acc_list)
            # remove any accession that is for a module (ie, not an enzyme)
            enzymes_in_def = set([acc for acc in accessions_in_def if acc not in self.module_dict.keys()])
            enzymes_without_orth = enzymes_in_def.difference(orth_list)
            if self.data_source != 'KEGG' and enzymes_without_orth:
                bad_list = ", ".join(enzymes_without_orth)
                n = len(enzymes_without_orth)
                raise ConfigError(f"So, there is a thing. And that thing is that there {P('is an enzyme', n, alt='are some enzymes')} "
                                  f"in the DEFINITION string for module {mnum} that {P('does', n, alt='do')} not have an ORTHOLOGY line, "
                                  f"and {P('it', n, alt='they')} really should have one. Here {P('it is', n, alt='they are')}: {bad_list}")

            # every user module needs at least one annotation source
            if self.data_source != 'KEGG' and not module_has_annotation_source:
                os.remove(self.db_path)
                raise ConfigError(f"While parsing user module {mnum}, we noticed that it does not have a single "
                                  f"'ANNOTATION_SOURCE' field. We are sorry to tell you that this is not okay, "
                                  f"because every user-defined module requires at least one of those fields to tell "
                                  f"anvi'o where to find its gene annotations. So you should go and take a look at "
                                  f"the module file at {mod_file_path} and add one 'ANNOTATION_SOURCE' line for each "
                                  f"gene in the module definition, before re-trying this setup program. Thank you!")
            # every enzyme needs an annotation source
            if self.data_source != 'KEGG' and not all(orthology_to_annotation_source.values()):
                nones = [k for k,v in orthology_to_annotation_source.items() if not v]
                nones_str = ", ".join(nones)
                n = len(nones)
                raise ConfigError(f"*Dalek noises* EXTERMINATE! EXTERMINATE! It seems that your module {mnum} contains "
                                  f"{P('an enzyme that is', n, alt='some enzymes that are')} missing an ANNOTATION_SOURCE line. "
                                  f"Please go back in time and fix this. Here {P('is the enzyme', n, alt='are the enzymes')} in "
                                  f"question here: {nones_str}")

            num_modules_parsed += 1
        # once we are done parsing all modules, we store whatever db entries remain in the db_entries list
        # this is necessary because append_and_store() above only stores every 10000 entries
        self.progress.update("Storing final batch of module entries into DB")
        mod_table.store(self.db)

        self.progress.end()

        # warn user about parsing errors
        if self.num_corrected_errors > 0 or self.num_uncorrected_errors > 0:
            if anvio.DEBUG:
                self.run.warning("Several parsing errors were encountered while building the Modules DB. "
                                 "Below you will see which modules threw each type of parsing error. Note that modules which "
                                 "threw multiple errors will occur in the list as many times as it threw each error.")
                self.run.info("Bad line splitting (usually due to rogue or missing spaces)", self.parsing_error_dict["bad_line_splitting"])
                self.run.info("Bad KEGG code format (not corrected; possibly problematic)", self.parsing_error_dict["bad_kegg_code_format"])
            else: # less verbose
                self.run.warning("First things first - don't panic. Several parsing errors were encountered while building the Modules DB. "
                                 "But that is probably okay, because if you got to this point it is likely that we already fixed all of them "
                                 "ourselves. So don't worry too much. Below you will see how many of each type of error was encountered. If "
                                 "you would like to see which modules threw these errors, please re-run the setup using the --debug flag (you "
                                 "will also probably need the --reset or --overwrite-output-destinations flag). When doing so, you will also "
                                 "see which lines caused issues; this can be a lot of output, so you can suppress the line-specific output with "
                                 "the `--quiet` flag if that makes things easier to read. So, in summary: You can probably ignore this warning. "
                                 "But if you want more info: run setup again with `--reset --debug --quiet` to see exactly which modules had "
                                 "issues, or run with `--reset --debug` to see exactly which lines in which modules had issues. Anvi'o developers "
                                 "thank you for your attention and patience 😇")
                self.run.info("Bad line splitting (usually due to rogue or missing spaces)", len(self.parsing_error_dict["bad_line_splitting"]))
                self.run.info("Bad KEGG code format (usually not correctable)", len(self.parsing_error_dict["bad_kegg_code_format"]))

        if not self.annotation_sources:
            os.remove(self.db_path)
            raise ConfigError("We're not sure how you made it this far without having any annotation sources defined in your module files, "
                              "because we should have noticed this while parsing them. But it happened, and here we are. You need to go add "
                              "'ANNOTATION_SOURCE' fields to those module files, and then re-do this setup.")
        annotation_source_list = ",".join(list(self.annotation_sources))


        # give some run info
        self.run.info('Modules database', 'A new database, %s, has been created.' % (self.db_path), quiet=self.quiet)
        self.run.info('Number of modules', num_modules_parsed, quiet=self.quiet)
        self.run.info('Number of module entries', mod_table.get_total_entries(), quiet=self.quiet)
        self.run.info('Number of module parsing errors (corrected)', self.num_corrected_errors, quiet=self.quiet)
        self.run.info('Number of module parsing errors (uncorrected)', self.num_uncorrected_errors, quiet=self.quiet)
        self.run.info('Annotation sources required for estimation', ", ".join(self.annotation_sources))

        # record some useful metadata
        self.db.set_meta_value('db_type', 'modules')
        self.db.set_meta_value('data_source', self.data_source)
        self.db.set_meta_value('annotation_sources', annotation_source_list)
        self.db.set_meta_value('num_modules', num_modules_parsed)
        self.db.set_meta_value('total_module_entries', mod_table.get_total_entries())
        self.db.set_meta_value('creation_date', time.time())
        self.db.set_meta_value('hash', self.get_db_content_hash())
        self.db.set_meta_value('version', anvio.__kegg_modules_version__)

        self.db.disconnect()


    def disconnect(self):
        self.db.disconnect()

######### SELF TABLE ACCESS FUNCTIONS #########

    def get_days_since_creation(self):
        """Returns the time (in days) since MODULES.db was created.

        The time units are seconds, and there are 60*60*24 = 86400 seconds per day,
        so we do the appropriate division to get the time in days.
        """

        return round((time.time() - float(self.db.get_meta_value('creation_date'))) / 86400)


    def get_db_content_hash(self):
        """Compute hash of all KOs and module numbers present in the db (used for tracking major changes to db content with future KEGG updates)"""
        mods = self.get_all_modules_as_list()
        mods.sort()
        orths = self.get_all_knums_as_list()
        orths.sort()
        mods_and_orths = mods + orths
        mods_and_orths = "".join(mods_and_orths)
        return str(hashlib.sha224(mods_and_orths.encode('utf-8')).hexdigest())[0:12]

######### MODULES TABLE ACCESS FUNCTIONS #########

    def get_modules_table_as_dict(self, data_names_of_interest=[]):
        """This function loads the modules table and returns it as a dictionary keyed by module.

        Every data_name for the module (NAME, DEFINITION, ORTHOLOGY, etc) will become a key in the
        inner dictionary, and the corresponding data_value will become its value. For data_names that
        have multiple values (ORTHOLOGY, COMPOUND, REACTION, etc), the value becomes yet another dictionary
        of data_value -> data_definition key-value pairs.

        The one exception is DEFINITION lines - there are occasionally multiple of these for one module but
        these must be returned as a list, not as a dictionary.

        PARAMETERS
        ==========
        data_names_of_interest : list of str
            the returned dictionary will contain only data names from this list. If the list is empty,
            all data names are returned.

        RETURNS
        =======
        module_dictionary : dict of dicts
            data for each module in the modules table. Outer dictionary is keyed by module number and
            inner dictionary is keyed by data name
        """

        if data_names_of_interest:
            data_names_list = [f"'{n}'" for n in data_names_of_interest]
            where_clause_string = f"data_name in ({','.join(data_names_list)})"
            # this WILL fail if you ask for a data name that doesn't exist, so know your data before you query
            dict_from_mod_table = self.db.get_some_rows_from_table_as_dict(self.module_table_name, where_clause_string, row_num_as_key=True)
        else:
            dict_from_mod_table = self.db.get_table_as_dict(self.module_table_name, row_num_as_key=True)
        # the returned dictionary is keyed by an arbitrary integer, and each value is a dict containing one row from the modules table
        # ex of one row in this dict: 0: {'module': 'M00001', 'data_name': 'ENTRY', 'data_value': 'M00001', 'data_definition': 'Pathway', 'line': 1}

        # now we convert this to a per-module dictionary
        module_dictionary = {}
        for entry in dict_from_mod_table:
            mod = dict_from_mod_table[entry]['module']
            data_name = dict_from_mod_table[entry]['data_name']
            data_value = dict_from_mod_table[entry]['data_value']
            data_definition = dict_from_mod_table[entry]['data_definition']

            if mod not in module_dictionary:
                module_dictionary[mod] = {}

            if data_name not in module_dictionary[mod]:
                if not data_definition:
                    module_dictionary[mod][data_name] = data_value
                else:
                    module_dictionary[mod][data_name] = {data_value : data_definition}
            else:
                # place multiple definition lines into list
                if data_name == "DEFINITION":
                    if isinstance(module_dictionary[mod][data_name], list):
                        module_dictionary[mod][data_name].append(data_value)
                    else:
                        existing_val = module_dictionary[mod][data_name]
                        module_dictionary[mod][data_name] = [existing_val, data_value]
                else:
                    # data_value -> data_definition dictionary
                    if isinstance(module_dictionary[mod][data_name], dict):
                        if data_value in module_dictionary[mod][data_name]:
                            module_dictionary[mod][data_name][data_value] += " / " + data_definition
                        else:
                            module_dictionary[mod][data_name][data_value] = data_definition
                    else:
                        existing_val = module_dictionary[mod][data_name]
                        existing_val_data_def = dict_from_mod_table[entry-1]['data_definition']
                        module_dictionary[mod][data_name] = {existing_val: existing_val_data_def, data_value: data_definition}

        return module_dictionary


    def get_data_value_entries_for_module_by_data_name(self, module_num, data_name):
        """This function returns data_value elements from the modules table for the specified module and data_name pair.

        All elements corresponding to the pair (ie, M00001 and ORTHOLOGY) will be returned.
        The function relies on the db.get_some_rows_from_table_as_dict() function to first fetch all rows corresponding \
        to a particular model, and then parses the resulting dictionary to find all the elements with the given data_name field.

        PARAMETERS
        ==========
        module_num : str
            the module to fetch data for
        data_name : str
            which data_name field we want

        RETURNS
        =======
        data_values_to_ret : list of str
            the data_values corresponding to the module/data_name pair
        """

        where_clause_string = "module = '%s'" % (module_num)
        dict_from_mod_table = self.db.get_some_rows_from_table_as_dict(self.module_table_name, where_clause_string, row_num_as_key=True)
        # the returned dictionary is keyed by an arbitrary integer, and each value is a dict containing one row from the modules table
        # ex of one row in this dict: 0: {'module': 'M00001', 'data_name': 'ENTRY', 'data_value': 'M00001', 'data_definition': 'Pathway', 'line': 1}
        data_values_to_ret = []
        for key in dict_from_mod_table.keys():
            if dict_from_mod_table[key]['data_name'] == data_name:
                data_values_to_ret.append(dict_from_mod_table[key]['data_value'])

        if not data_values_to_ret:
            self.run.warning("Just so you know, we tried to fetch data from the Modules database for the data_name field %s "
                             "and module %s, but didn't come up with anything, so an empty list is being returned. This may "
                             "cause errors down the line, and if so we're very sorry for that.")

        return data_values_to_ret


    def get_data_definition_entries_for_module_by_data_name(self, module_num, data_name):
        """This function returns data_definition elements from the modules table for the specified module and data_name pair.

        All elements corresponding to the pair (ie, M00001 and ORTHOLOGY) will be returned.
        The function relies on the db.get_some_rows_from_table_as_dict() function to first fetch all rows corresponding \
        to a particular model, and then parses the resulting dictionary to find all the elements with the given data_name field.

        PARAMETERS
        ==========
        module_num : str
            the module to fetch data for
        data_name : str
            which data_name field we want

        RETURNS
        =======
        data_defs_to_ret : list of str
            the data_definitions corresponding to the module/data_name pair
        """

        where_clause_string = "module = '%s'" % (module_num)
        dict_from_mod_table = self.db.get_some_rows_from_table_as_dict(self.module_table_name, where_clause_string, row_num_as_key=True)

        data_defs_to_ret = []
        for key in dict_from_mod_table.keys():
            if dict_from_mod_table[key]['data_name'] == data_name:
                data_defs_to_ret.append(dict_from_mod_table[key]['data_definition'])

        if not data_defs_to_ret and anvio.DEBUG:
            self.run.warning("Just so you know, we tried to fetch data definitions from the Modules database for the data_name field %s "
                             "and module %s, but didn't come up with anything, so an empty list is being returned. This may "
                             "cause errors down the line, and if so we're very sorry for that.")

        return data_defs_to_ret


    def get_all_modules_as_list(self):
        """This function returns a list of all modules in the DB."""
        return self.db.get_single_column_from_table(self.module_table_name, 'module', unique=True)


    def get_all_knums_as_list(self):
        """This function returns a list of all KO numbers in the DB."""
        where_clause_string = "data_name = 'ORTHOLOGY'"
        return self.db.get_single_column_from_table(self.module_table_name, 'data_value', unique=True, where_clause=where_clause_string)


    def get_ko_function_dict(self):
        """This function returns a 2-level dictionary keyed by KO number. The inner dict contains a 'definition' field
           that stores the KO's functional annotation.

           This method effectively builds a partial KO dict similar to the one produced by the setup_ko_dict() function,
           but containing only the KOs/gene annotations from the modules DB. When being used for user-defined metabolism,
           it should be expanded later with gene annotations from the contigs database(s) being worked on.
        """
        where_clause_string = "data_name = 'ORTHOLOGY'"
        kos_and_functions = self.db.get_some_columns_from_table(self.module_table_name, "data_value,data_definition", unique=True, where_clause=where_clause_string)
        ko_func_dict = {}
        for k,f in kos_and_functions:
            if k not in ko_func_dict:
                ko_func_dict[k] = {'definition': f }
        return ko_func_dict


    def get_modules_for_knum(self, knum):
        """This function returns a list of modules that the given KO belongs to."""
        where_clause_string = "data_value = '%s'" % (knum)
        return self.db.get_single_column_from_table(self.module_table_name, 'module', unique=True, where_clause=where_clause_string)


    def get_module_classes_for_knum_as_dict(self, knum):
        """This function returns the classes for the modules that a given KO belongs to in a dictionary of dictionaries keyed by module number."""
        mods = self.get_modules_for_knum(knum)
        all_mods_classes_dict = {}
        for mnum in mods:
            all_mods_classes_dict[mnum] = self.get_kegg_module_class_dict(mnum)
        return all_mods_classes_dict


    def get_module_classes_for_knum_as_list(self, knum):
        """This function returns the classes for the modules that a given KO belongs to as a list of strings."""
        mods = self.get_modules_for_knum(knum)
        all_mods_classes_list = []
        for mnum in mods:
            mod_class = self.get_data_value_entries_for_module_by_data_name(mnum, "CLASS")[0]
            all_mods_classes_list.append(mod_class)
        return all_mods_classes_list


    def get_module_name(self, mnum):
        """This function returns the name of the specified KEGG module."""

        # there should only be one NAME per module, so we return the first list element
        return self.get_data_value_entries_for_module_by_data_name(mnum, "NAME")[0]


    def get_module_names_for_knum(self, knum):
        """This function returns all names of each KEGG module that the given KO belongs to in a dictionary keyed by module number."""
        mods = self.get_modules_for_knum(knum)
        module_names = {}
        for mnum in mods:
            module_names[mnum] = self.get_module_name(mnum)
        return module_names


    def parse_kegg_class_value(self, class_data_val):
        """This function takes a data_value string for the CLASS field in the modules table and parses it into a dictionary.

        The data_value string of CLASS fields should look something like this: Pathway modules; Amino acid metabolism; Lysine metabolism
        so they can be parsed into 3 parts: class, category, and subcategory.
        """

        fields = class_data_val.split("; ")
        class_dict = {"class" : fields[0], "category" : fields[1], "subcategory" : fields[2] if len(fields) > 2 else None}
        return class_dict


    def get_kegg_module_class_dict(self, mnum, class_value=None):
        """This function returns a dictionary of values in the CLASS field for a specific module

        It really exists only for convenience to put together the data fetch and parsing functions.

        PARAMETERS
        ==========
        mnum : str
            the module number
        class_value : str
            The 'CLASS' string for the module. This parameter is optional, and if it is not provided,
            the 'CLASS' value will be queried from the modules DB.
        """

        if not class_value:
            # there should only be one CLASS line per module, so we extract the first list element
            class_value = self.get_data_value_entries_for_module_by_data_name(mnum, "CLASS")[0]
        return self.parse_kegg_class_value(class_value)


    def get_kegg_module_definition(self, mnum):
        """This function returns module DEFINITION fields as one string"""

        def_lines = self.get_data_value_entries_for_module_by_data_name(mnum, "DEFINITION")
        def_lines = [l.strip() for l in def_lines] # sometimes there are stray spaces at the end of the string that will mess us up later
        return " ".join(def_lines)


    def get_ko_definition_from_modules_table(self, ko_num):
        """This function returns the definition for the given KO from the modules data table.

        Note that the modules table will only contain information for KOs that belong to modules, so this
        function returns None for those KOs that are not in modules. If your use case depends on accessing
        definitions for all KOs, you are better off calling KeggContext.setup_ko_dict() and taking the
        definition from that dictionary.
        """

        where_clause_string = "data_name = 'ORTHOLOGY' AND data_value = '%s'" % (ko_num)
        dict_from_mod_table = self.db.get_some_rows_from_table_as_dict(self.module_table_name, where_clause_string, row_num_as_key=True, error_if_no_data=False)
        if not dict_from_mod_table:
            self.run.warning("get_ko_definition() speaking: No ORTHOLOGY entry found for KO %s - returning None."
                            % (ko_num))
            return None
        else:
            # there could be several rows for the same KO in different modules, but each definition should be
            # the same or similar, so we arbitrarily return the first one
            return dict_from_mod_table[0]['data_definition']


    def get_ko_reactions_from_modules_table(self, ko_num):
        """This function returns the KEGG reaction ID for the given KO from its data definition entry in the modules table.

        Reactions are indicated within brackets of the data definition entry, like these: [RN:R05339] or [RN:R01538 R03033].
        This function parses all reactions out of the entry and returns a list in which each reaction ID number is prefixed by
        the standard KEGG reaction indicator 'RN:', as in ["RN:R01538", "RN:R03033"].

        Note that the modules table will only contain information for KOs that belong to modules, so this function
        returns None for those KOs that are not in modules.
        """

        definition_line = self.get_ko_definition_from_modules_table(ko_num)
        if not definition_line: # this KO was not in the modules db
            return None

        def_fields = definition_line.split('[') # the last split should start with RN: and end with ]
        for f in def_fields:
            if f.startswith("RN:"):
                react_ids = f[3:-1].split(' ') # extract the ID numbers (without the initial RN: part or the final ])
                return ["RN:" + id for id in react_ids]



    def get_kos_in_module(self, mnum):
        """This function returns a list of KOs in the given module.

        It does this by parsing the ORTHOLOGY lines in the modules database. However,
        please note that these KOs are not always in the same order as the module
        definition, and may even contain duplicate entries for a KO. A good example
        of this is http://rest.kegg.jp/get/M00091 (K00551 is in two ORTHOLOGY lines)
        and http://rest.kegg.jp/get/M00176 (see KOs in the first top-level step). If
        this will be a problem, you should use the function get_kos_from_module_definition()
        instead.
        """

        return self.get_data_value_entries_for_module_by_data_name(mnum, "ORTHOLOGY")


    def get_kos_from_module_definition(self, mnum):
        """This function returns a list of KOs in the given module, in order of the DEFINITION.

        An alternative to get_kos_in_module().
        """

        mod_def = self.get_kegg_module_definition(mnum)
        ko_list = []
        k_indices = [x for x, v in enumerate(mod_def) if v == 'K']
        for idx in k_indices:
            ko_list.append(mod_def[idx:idx+6])

        return ko_list


    def get_kegg_module_compound_lists(self, mnum):
        """This function returns a list of substrates, a list of intermediates, and a list of products for the given module.

        We define 'substrate' to be any compound that is an input to but not an output from reactions in the module pathway.
        Likewise, a 'product' is any compound that is an output from but not an input to reactions in the module pathway.
        'Intermediate' is a compound that is both an input to and and output from reactions in the pathway.

        Note that this function refers to compounds by their KEGG identifier (format is 'C#####' where # is a digit).
        A separate function is used to convert these lists to human-readable compound names.

        RETURNS
        =======
        substrates : list
            Compounds that are only inputs to the module's metabolic pathway
        intermediates : list
            Compounds that are both outputs and inputs in the module's metabolic reactions
        products : list
            Compunds that are only outputs from the module's metabolic pathway
        """

        reactions_list = self.get_data_definition_entries_for_module_by_data_name(mnum, "REACTION")
        if not reactions_list:
            if anvio.DEBUG:
                self.run.warning(f"No REACTION entries found for module {mnum}, so no compounds will be returned by "
                                 "get_kegg_module_compound_lists()")

        inputs = set([])
        outputs = set([])

        for rxn_string in reactions_list:
            if '<->' in rxn_string:
                split_rxn = rxn_string.split('<->')
            else:
                split_rxn = rxn_string.split('->')
            if len(split_rxn) != 2:
                raise ConfigError(f"get_kegg_module_compound_lists('{mnum}') ran into an issue splitting the reaction {rxn_string}"
                                  "into 2 parts. Here is what the split looks like: {split_rxn}")
            rxn_inputs = [x.strip() for x in split_rxn[0].split('+')]
            rxn_outputs = [x.strip() for x in split_rxn[1].split('+')]
            # as of December 2023, some REACTION lines now include stoichiometry, like this:
            # C00390 + 2 C00125 -> C00399 + 2 C00126 + 2 C00080
            # to make sure we don't keep the stoichiometric numbers, we need to split each compound substring one more time
            # on a space and keep only the last element
            rxn_inputs = [x.split()[-1] for x in rxn_inputs]
            rxn_outputs = [x.split()[-1] for x in rxn_outputs]
            inputs = inputs.union(set(rxn_inputs))
            outputs = outputs.union(set(rxn_outputs))

        substrates = inputs.difference(outputs)
        products = outputs.difference(inputs)
        intermediates = inputs.intersection(outputs)

        return list(substrates), list(intermediates), list(products)


    def get_compound_dict_for_module(self, mnum, raise_error_if_no_data=False):
        """This function returns a dictionary mapping compound identifiers to their human-readable name for the given module

        If the module has no compounds, this function will either raise an error or return an empty dictionary depending on raise_error_if_no_data.
        If a compound doesn't have a human-readable name, then the compound identifier is used as the 'name'

        PARAMETERS
        ==========
        mnum : str
            module number to get compounds for
        raise_error_if_no_data : bool
            whether to quit all things if we don't get what we want
        """

        where_clause_string = "data_name = 'COMPOUND' AND module = '%s'" % (mnum)
        dict_from_mod_table = self.db.get_some_rows_from_table_as_dict(self.module_table_name, where_clause_string, row_num_as_key=True, error_if_no_data=raise_error_if_no_data)
        compound_dict = {}
        for key,row in dict_from_mod_table.items():
            compound = row['data_value']
            compound_name = row['data_definition']
            # if compound has no human-readable name in the database, we use the compound ID after all
            if not compound_name:
                compound_name = compound
            compound_dict[compound] = compound_name

        return compound_dict


    def get_human_readable_compound_lists_for_module(self, mnum):
        """This function returns a human-readable list of substrates, a list of intermediates, and a list of products for the given module.

        We define 'substrate' to be any compound that is an input to but not an output from reactions in the module pathway.
        Likewise, a 'product' is any compound that is an output from but not an input to reactions in the module pathway.
        'Intermediate' is a compound that is both an input to and and output from reactions in the pathway.

        RETURNS
        =======
        substrate_name_list : list of str
            List of substrate compounds
        intermediate_name_list : list of str
            List of intermediate compounds
        product_name_list : list of str
            List of product compounds
        """
        compound_to_name_dict = self.get_compound_dict_for_module(mnum)
        substrate_compounds, intermediate_compounds, product_compounds = self.get_kegg_module_compound_lists(mnum)

        substrate_name_list = [compound_to_name_dict[c] for c in substrate_compounds]
        intermediate_name_list = [compound_to_name_dict[c] for c in intermediate_compounds]
        product_name_list = [compound_to_name_dict[c] for c in product_compounds]

        return substrate_name_list, intermediate_name_list, product_name_list


######### MODULE DEFINITION UNROLLING FUNCTIONS #########

    def get_top_level_steps_in_module_definition(self, mnum):
        """This function access the DEFINITION line of a KEGG Module and returns the top-level steps as a list

        A 'top-level' step is one that you get by splitting on spaces (but not spaces in parentheses) just once -
        ie, the 'first layer' when unrolling the module.
        """

        def_string = self.get_kegg_module_definition(mnum)
        return utils.split_by_delim_not_within_parens(def_string, " ")


    def unroll_module_definition(self, mnum, def_lines = None):
        """This function accesses the DEFINITION line of a KEGG Module, unrolls it into all possible paths through the module, and
        returns the list of all paths.

        This is a driver for the recursive functions that do the actual unrolling of each definition line.

        PARAMETERS
        ==========
        mnum : str
            module number
        def_lines : list of str
            The DEFINITION lines for the module. This parameter is optional, and if it is not passed, the module
            definition will be looked up from the modules DB.
        """

        if not def_lines:
            def_lines = self.get_data_value_entries_for_module_by_data_name(mnum, "DEFINITION")
        combined_def_line = ""
        for d in def_lines:
            d = d.strip()
            combined_def_line += d + " "
        combined_def_line = combined_def_line.strip()
        
        try:
            def_line_paths = self.recursive_definition_unroller(combined_def_line)
        except RecursionError as re:
            raise ConfigError(f"Uh oh. While unrolling the definition of module {mnum}, we got an error. "
                              f"Unfortunately, we don't know much at this point, but here is the error message: '{re}'. "
                              f"At this point, our best guess is that something is wrong with the module defintion, so please "
                              f"take a look at the following DEFINITION lines and see if anything looks fishy (like stray "
                              f"spaces, missing spaces, or unbalanced parentheses): {combined_def_line}")

        return def_line_paths


    def recursive_definition_unroller(self, step):
        """This function recursively splits a module definition into its components.

        First, the definition is split into its component steps (separated by spaces).
        Each step is either an atomic step (a single KO, module number, '--', or nonessential KO starting with '-'),
        a protein complex, or a compound step.

        Atomic steps are used to extend each path that has been found so far. Protein complexes are split into
        their respective components, which may be split further by the split_paths() function to find all possible
        alternative complexes, before being used to extend each path. Compound steps are split and recursively processed
        by the split_paths() function before the resulting downstream paths are used to extend each path.

        PARAMETERS
        ==========
        step : str
            step definition to split into component steps as necessary

        RETURNS
        =======
        paths_list : list
            all paths that the input step has been unrolled into
        """

        split_steps = utils.split_by_delim_not_within_parens(step, " ")
        paths_list = [[]]  # list to save all paths, with initial empty path list to extend from
        for s in split_steps:
            # base case: step is a ko, mnum, non-essential step, or '--'
            if (len(s) == 6 and s[0] == "K") or (len(s) == 6 and s[0] == "M") or (s == "--") or (len(s) == 7 and s[0] == "-"):
                for p in paths_list:
                    p.extend([s])
            else:
                if s[0] == "(" and s[-1] == ")":
                    # here we try splitting to see if removing the outer parentheses will make the definition become unbalanced
                    # (the only way to figure this out is to try it because regex cannot handle nested parentheses)
                    comma_substeps = utils.split_by_delim_not_within_parens(s[1:-1], ",")
                    if not comma_substeps: # if it doesn't work, try without removing surrounding parentheses
                        comma_substeps = utils.split_by_delim_not_within_parens(s, ",")
                    space_substeps = utils.split_by_delim_not_within_parens(s[1:-1], " ")
                    if not space_substeps:
                        space_substeps = utils.split_by_delim_not_within_parens(s, " ")
                else:
                    comma_substeps = utils.split_by_delim_not_within_parens(s, ",")
                    space_substeps = utils.split_by_delim_not_within_parens(s, " ")

                # complex case: no commas OR spaces outside parentheses so this is a protein complex rather than a compound step
                if len(comma_substeps) == 1 and len(space_substeps) == 1:
                    complex_components, delimiters = utils.split_by_delim_not_within_parens(s, ["+","-"], return_delims=True)
                    complex_strs = [""]

                    # reconstruct the complex (and any alternate possible complexes) while keeping the +/- structure the same
                    for i in range(len(complex_components)):
                        c = complex_components[i]
                        if c[0] == '(':
                            alts = self.split_path(c)
                            new_complex_strs = []
                            for a in alts:
                                if len(a) > 1:
                                    raise ConfigError("Uh oh. recursive_definition_unroller() speaking. We found a protein complex with more "
                                                    "than one KO per alternative option here: %s" % s)
                                for cs in complex_strs:
                                    extended_complex = cs + a[0]
                                    new_complex_strs.append(extended_complex)
                            complex_strs = new_complex_strs
                        else:
                            for j in range(len(complex_strs)):
                                complex_strs[j] += c

                        if i < len(delimiters):
                            for j in range(len(complex_strs)):
                                complex_strs[j] += delimiters[i]

                    new_paths_list = []
                    for cs in complex_strs:
                        for p in paths_list:
                            p_copy = copy.copy(p)
                            p_copy.extend([cs])
                            new_paths_list.append(p_copy)
                    paths_list = new_paths_list

                # compound step case:
                else:
                    alts = self.split_path(s)
                    new_paths_list = []
                    for a in alts:
                        for p in paths_list:
                            p_copy = copy.copy(p)
                            p_copy.extend(a)
                            new_paths_list.append(p_copy)
                    paths_list = new_paths_list

        return paths_list


    def split_path(self, step):
        """This function handles compound steps that should be split into multiple alternative paths.

        It first splits the input step into substeps, and then since each substep could be its own mini-definition,
        it recursively calls the definition unrolling function to parse it. The list of all alternative paths
        that can be made from this step is returned.
        """

        if step[0] == "(" and step[-1] == ")":
            substeps = utils.split_by_delim_not_within_parens(step[1:-1], ",")
            if not substeps: # if it doesn't work, try without removing surrounding parentheses
                substeps = utils.split_by_delim_not_within_parens(step, ",")
        else:
            substeps = utils.split_by_delim_not_within_parens(step, ",")

        alt_path_list = []
        for s in substeps:
            alt_paths_from_substep = self.recursive_definition_unroller(s)
            for a in alt_paths_from_substep:
                alt_path_list.append(a)

        return alt_path_list


class ModulesTable:
    """This class defines operations for creating the KEGG Modules table in Modules.db"""

    def __init__(self, mod_table_name = None):
        """"""
        self.db_entries = []
        self.total_entries = 0

        if mod_table_name:
            self.module_table_name = mod_table_name
        else:
            raise ConfigError("Beep Beep. Warning. ModulesTable was initialized without knowing its own name.")


    def append_and_store(self, db, module_num, data_name, data_value, data_definition=None, line_num=None):
        """This function handles collects db entries (as tuples) into a list, and once we have 10,000 of them it stores that set into the Modules table.

        The db_entries list is cleared after each store so that future stores don't add duplicate entries to the table.
        """

        db_entry = tuple([module_num, data_name, data_value, data_definition, line_num])
        self.db_entries.append(db_entry)
        self.total_entries += 1

        # we can store chunks of 5000 at a time, so we don't want over 10,000 entries.
        if len(self.db_entries) >= 10000:
            self.store(db)
            self.db_entries = []


    def store(self, db):
        if len(self.db_entries):
            db._exec_many('''INSERT INTO %s VALUES (%s)''' % (self.module_table_name, (','.join(['?'] * len(self.db_entries[0])))), self.db_entries)


    def get_total_entries(self):
        return self.total_entries


