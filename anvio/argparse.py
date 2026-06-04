"""Overloading Python argparse for anvi'o purposes"""

import os
import sys
import argparse
import textwrap

from colored import fg, attr
from rich_argparse import RichHelpFormatter

import anvio
import anvio.docs as docs
import anvio.terminal as terminal

from anvio.programs import Program
from anvio.errors import ConfigError
from anvio.dbinfo import FindAnvioDBs


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


atty = sys.stdout.isatty()
P = terminal.pluralize


class ArgumentParser(argparse.ArgumentParser):
    num_threads_not_provided = object()

    def __init__(self, description="No description :/", epilog=None):
        super().__init__()

        self.description = description
        self.epilog = epilog or self.get_anvio_epilogue()

        self.anvio_allowed_ad_hoc_flags = ['--version', '--debug', '--force', '--fix-sad-tables',
                                           '--quiet', '--no-progress', '--as-markdown', '--display-db-calls',
                                           '--force-use-my-tree', '--force-overwrite', '--tmp-dir',
                                           '--I-know-this-is-not-a-good-idea', '--include-stray-KOs']


    def get_num_threads_action(self):
        for action in self._actions:
            if action.dest == 'num_threads':
                return action

        return None


    def get_num_threads_from_environment(self):
        try:
            num_threads = int(os.environ['ANVIO_THREADS'])
        except ValueError:
            raise ConfigError("The environmental variable `ANVIO_THREADS` must be a positive integer. "
                              f"But anvi'o found `{os.environ['ANVIO_THREADS']}` there :/")

        if num_threads <= 0:
            raise ConfigError(f"The environmental variable `ANVIO_THREADS` must be a positive integer. "
                              f"`{num_threads}` is not it :/")

        return num_threads


    def parse_known_args(self, args=None, namespace=None):
        num_threads_action = self.get_num_threads_action()

        if not num_threads_action:
            return super().parse_known_args(args, namespace)

        args_to_parse = sys.argv[1:] if args is None else args
        if '-h' in args_to_parse or '--help' in args_to_parse:
            return super().parse_known_args(args, namespace)

        default_num_threads = num_threads_action.default
        num_threads_action.default = self.num_threads_not_provided

        try:
            args, unknown = super().parse_known_args(args, namespace)
        finally:
            num_threads_action.default = default_num_threads

        if args.num_threads is self.num_threads_not_provided:
            if 'ANVIO_THREADS' in os.environ:
                args.num_threads = self.get_num_threads_from_environment()
            else:
                args.num_threads = default_num_threads

        return args, unknown


    def get_anvio_epilogue(self):
        """Function that formats the additional message that appears at the end of help."""

        # Here we intend to collect the requires and provides statements from anvi'o programs.
        # but the recovery of this 'epilogue' will yield an error for scripts that are not listed
        # in the $PATH, which can happen during the earlier stages of program development. so,
        # if we can't recover the epilogue, we go long hair don't care and return none.
        try:
            epilog = self.get_requires_and_provides_statements_for_program()
        except:
            return None

        version = anvio.anvio_version_for_help_docs

        general_help = f"https://anvio.org/help/{version}"
        program_help = f"{general_help}/programs/{self.prog}"

        if os.path.exists(os.path.join(os.path.dirname(docs.__file__), f"programs/{self.prog}.md")):
            if atty:
                epilog += f'''\n🍺 {attr('bold')}More on `{self.prog}`:{attr('reset')}\n\n   {fg('cyan') + program_help + attr('reset')}'''
            else:
                epilog += f'''\n🍺 More on `{self.prog}`:\n\n   {program_help}'''
        else:
            epilog = ""

        if atty:
            epilog += f'''\n\n🍻 {attr('bold')}All anvi'o programs and artifacts:{attr('reset')}\n\n   {fg('cyan') + general_help + attr('reset')}'''
        else:
            epilog += f'''\n\n🍻 All anvi'o programs and artifacts:\n\n   {general_help}'''

        if atty:
            return epilog + attr('reset')
        else:
            return epilog


    def get_requires_and_provides_statements_for_program(self):
        """Formats and returns requires and provides statements for the program"""

        requires_and_provides_statements = []

        def _get_values_for(statement):
            return [v.id for v in program.meta_info[statement]['value']]

        program = Program(self.prog)
        requires = _get_values_for('requires')
        provides = _get_values_for('provides')
        can_use = _get_values_for('can_use')
        can_provide = _get_values_for('can_provide')

        def get_block(statement, header=None, extra_items=None, extra_verb=''):
            blocks = []
            if len(statement):
                split = []
                for item in statement:
                    addition = [item, " / "] if statement[-1] != item else [item]
                    if len(' '.join(['   '] + split + addition)) > 90:
                        blocks.append(split)
                        split = addition
                    else:
                        split.extend(addition)

                blocks.append(split)
                blocks.append([''])

            # index of the last content block (before the trailing empty one)
            suffix_block_idx = len(blocks) - 2 if extra_items and len(blocks) >= 2 else None

            for i in range(0, len(blocks)):
                block = blocks[i]
                if atty:
                    block = [f"{fg('red') + b + attr('reset')}" if not b.strip() == '/' else b for b in block]
                else:
                    pass

                if len(block) > 1:
                    block.insert(0, '   ')

                # build suffix with artifact names colored, surrounding text plain
                if i == suffix_block_idx:
                    if atty:
                        colored = [fg('red') + item + attr('reset') for item in extra_items]
                    else:
                        colored = extra_items
                    block.append(f' (can also {extra_verb} {" / ".join(colored)})')

                blocks[i] = block

            if header:
                if atty:
                    blocks.insert(0, attr('bold') + header + attr('reset'))
                else:
                    blocks.insert(0, header)

            blocks.insert(1, "")

            blocks = [''.join(b) for b in blocks]

            return '\n'.join(blocks)


        requires_and_provides_statements.append(get_block(requires, """🧀 Requires:""", extra_items=can_use, extra_verb='use'))
        requires_and_provides_statements.append(get_block(provides, """🍕 Provides:""", extra_items=can_provide, extra_verb='generate'))

        return '\n' + '\n'.join(requires_and_provides_statements)


    def format_help(self):
        """Individual formatting of sections in the help text.

        When we use the same formatter for all, we either would lose the
        explicit spacing in the epilog, or lose the formatting in other
        sections. In this function we change the formatters, render
        different sections differently, and then concatenate everything
        into a single output.
        """

        RichHelpFormatter.styles["argparse.text"] = "italic"
        RichHelpFormatter.group_name_formatter = str.upper

        # we get our formatters here, fill them up down below, and finally render them at the end.
        if atty:
            usage_formatter = RichHelpFormatter(self.prog)
        else:
            usage_formatter = argparse.ArgumentDefaultsHelpFormatter(self.prog)

        description_formatter = argparse.RawDescriptionHelpFormatter(self.prog)
        epilog_formatter = argparse.RawDescriptionHelpFormatter(prog=self.prog)
        separator_formatter = argparse.RawDescriptionHelpFormatter(prog=self.prog)

        # usage
        usage_formatter.add_usage(self.usage, self._actions, self._mutually_exclusive_groups)

        # positionals, optionals and user-defined groups
        for action_group in self._action_groups:
            section_header = action_group.title

            usage_formatter.start_section(section_header)
            usage_formatter.add_text(action_group.description)
            usage_formatter.add_arguments(action_group._group_actions)
            usage_formatter.end_section()

        # separator
        separator_formatter.add_text('━' * 80 + '\n')

        # description
        if atty:
            description_text = [attr('bold') + '🔥 Program description:' + attr('reset'), '']
        else:
            description_text = ['🔥 Program description:', '']

        description_text.extend([textwrap.indent(l, '   ') for l in textwrap.wrap(" ".join(textwrap.dedent(self.description).split()), width=77)])
        description_formatter.add_text('\n'.join(description_text))

        # epilog
        epilog_formatter.add_text(self.epilog)

        # determine help from format above
        help_text = '\n'.join([usage_formatter.format_help().replace(":\n", "\n") if atty else usage_formatter.format_help(),
                               separator_formatter.format_help(),
                               description_formatter.format_help(),
                               epilog_formatter.format_help(),
                               separator_formatter.format_help()]) + '\n'

        return help_text


    def sanity_check(self, args):
        """Simple sanity checks for global arguments"""

        # make sure num threads are not 0 or negative.
        if args and 'num_threads' in args:
            try:
                args.num_threads = int(args.num_threads)
            except:
                raise ConfigError(f"Among your parameters anvi'o found one that sets the number of threads, "
                                  f"a value that should be a positive integer. But '{args.num_threads}' is not :/")

            if args.num_threads <= 0:
                raise ConfigError(f"The number of threads must be a positive integer. `{args.num_threads}` is not it :/")


    def get_args(self, parser, auto_fill_anvio_dbs=False):
        """A helper function to parse args anvi'o way.

        This function allows us to make sure some ad hoc parameters, such as `--debug`,
        can be used with any anvi'o program spontaneously even if they are not explicitly
        defined as an accepted argument, yet flags (or parameters) anvi'o does not expect
        to see can still be sorted out.
        """

        args, unknown = parser.parse_known_args()

        self.sanity_check(args)

        # A little historical luggage. We used to have the flag `--include-stray-KOs` in
        # metabolism framework associated programs such as `anvi-estimate-metabolism` or
        # `anvi-run-kegg-kofams`. Versions of these programs with this flag ended up in
        # various publications already, but in the meantime we realized that `--include-nt-KOs`
        # was a better name for this flag. So, here we are implementing a solution for backwards
        # compatibility.
        if '--include-stray-KOs' in unknown:
            args.include_nt_KOs = True

        if auto_fill_anvio_dbs:
            if anvio.DEBUG_AUTO_FILL_ANVIO_DBS:
                args = PopulateAnvioDBArgs(args, lazy_init=False).get_updated_args()
            else:
                args = PopulateAnvioDBArgs(args).get_updated_args()

        # if there are any args in the unknown that we do not expect to find
        # we we will make argparse complain about those.
        if len([f for f in unknown if f not in self.anvio_allowed_ad_hoc_flags]):
            for f in self.anvio_allowed_ad_hoc_flags:
                # handle non-boolean flags
                if f in ['--tmp-dir']:
                    parser.add_argument(f)
                else:
                    parser.add_argument(f, action='store_true')
            parser.parse_args()

        return args


class PopulateAnvioDBArgs(FindAnvioDBs):
    """A class to add missing anvi'o database arguments in an args instance

    This class does this by finding all anvi'o databases via `FindAnvioDBs`,
    and then reasoning with the information it generated.

    Parameters
    ==========
    args : Namespace instance
        Variable names matching to anvi'o databases defined below will be populated
        with paths in-place.
    search_path : str, default '.'
        The starting directory to search for all anvi'o databases underneath
    """

    def __init__(self, args, search_path='.', run=terminal.Run(), progress=terminal.Progress(), lazy_init=True):
        self.run = run
        self.progress = progress

        self.args = args
        self.search_path = search_path

        if lazy_init:
            self.anvio_dbs_found = False
            self.anvio_dbs = None
        else:
            FindAnvioDBs.__init__(self, run=self.run, progress=self.progress)

        self.__args_set = []
        self.__args_failed = []


    def set_arg(self, variable, value):
        self.__args_set.append((variable, value), )
        self.args.__dict__[variable] = value


    def fill_in_contigs_db(self, db_hash=None):
        if not 'contigs_db' in self.args:
            return

        if self.args.contigs_db:
            return

        FindAnvioDBs.__init__(self, run=self.run, progress=self.progress) if not self.anvio_dbs_found else None
        if 'contigs' not in self.anvio_dbs:
            return

        if db_hash:
            matching_contigs_dbs = [c for c in self.anvio_dbs['contigs'] if c.hash == db_hash]
        else:
            matching_contigs_dbs = self.anvio_dbs['contigs']

        if len(matching_contigs_dbs):
            self.set_arg('contigs_db', matching_contigs_dbs[0].path)
        else:
            self.__args_failed.append(('contigs_db', f'No matching contigs db (for hash "{db_hash}")'), )


    def fill_in_genomes_storage_db(self, db_hash=None):
        if not 'genomes_storage' in self.args:
            return

        if self.args.genomes_storage:
            return

        FindAnvioDBs.__init__(self, run=self.run, progress=self.progress) if not self.anvio_dbs_found else None
        if 'genomestorage' not in self.anvio_dbs:
            return

        if db_hash:
            matching_genomes_storage_dbs = [c for c in self.anvio_dbs['genomestorage'] if c.hash == db_hash]
        else:
            matching_genomes_storage_dbs = self.anvio_dbs['genomestorage']

        if len(matching_genomes_storage_dbs):
            self.set_arg('genomes_storage', matching_genomes_storage_dbs[0].path)
        else:
            self.__args_failed.append(('genomes_storage', f'No genomes storage around (for hash "{db_hash}")'), )


    def fill_in_profile_db(self):
        if 'profile_db' not in self.args:
            return

        FindAnvioDBs.__init__(self, run=self.run, progress=self.progress) if not self.anvio_dbs_found else None
        if not 'profile' in self.anvio_dbs:
            self.__args_failed.append(('profile_db', 'No profile databases around :/'))
            return

        profile_dbs = self.anvio_dbs['profile']

        if len(profile_dbs):
            # so we have some profile dbs. we can select the first one, or we can do something
            # slightly smarter and select the first one affiliated with a contigs database
            # if there are more than one
            if len(profile_dbs) > 1:
                try:
                    profile_db = [p for p in profile_dbs if p.hash][0]
                except:
                    # the exception here will come from the [0] in the previous line and will
                    # mean that although there were multiple profile databases, none was
                    # associated with a contigs db. FINE, we send back the first one of all
                    # these losers, then:
                    profile_db = profile_dbs[0]
            else:
                # there is only a single one. so it doesn't really matter
                profile_db = profile_dbs[0]
        else:
            # there is no profile db to be found around.
            self.__args_failed('profile_db', 'None around :/')
            return

        if not profile_db.hash:
            # the profile db is not associated with a contigs database
            # we shall try manual
            self.set_arg('profile_db', profile_db.path)
            self.set_arg('manual_mode', True)
        else:
            # it is associated with a contigs database. here we will set the
            # profile db, and next ask anvi'o to set the contigs db if it can
            # find one.
            self.set_arg('profile_db', profile_db.path)
            self.fill_in_contigs_db(db_hash=profile_db.hash)


    def fill_in_pan_db(self):
        if 'pan_db' not in self.args:
            return

        FindAnvioDBs.__init__(self, run=self.run, progress=self.progress) if not self.anvio_dbs_found else None
        if not 'pan' in self.anvio_dbs:
            return self.__args_failed.append(('pan_db', 'No pan databases around :/'))

        pan_dbs = self.anvio_dbs['pan']

        if len(pan_dbs):
            self.set_arg('pan_db', pan_dbs[0].path)
            self.fill_in_genomes_storage_db(db_hash=pan_dbs[0].hash)
        else:
            # there is no profile db to be found around.
            return self.__args_failed('pan_db', 'None around :/')


    def get_updated_args(self):
        if anvio.DEBUG_AUTO_FILL_ANVIO_DBS:
            self.run.warning(None, header="ANVI'O DBs FOUND", lc="yellow")
            if len(self.anvio_dbs):
                for db_type in self.anvio_dbs:
                    for anvio_db in self.anvio_dbs[db_type]:
                        print(anvio_db)
            else:
                self.run.info_single('lol no dbs around')

        if 'profile_db' in self.args and not self.args.profile_db:
            self.fill_in_profile_db()
        elif 'pan_db' in self.args and not self.args.pan_db:
            self.fill_in_pan_db()

        if len(self.__args_set):
            self.run.warning(None, header=f"ANVI'O FILLED IN THE {len(self.__args_set)} {P('ARG', len(self.__args_set), alt='ARGS')} BELOW AUTOMATICALLY", lc='yellow')
            for variable, value in self.__args_set:
                self.run.info(variable, value, nl_after= (1 if (variable, value) == self.__args_set[-1] else 0), lc="yellow")

        if len(self.__args_failed):
            self.run.warning(None, header="ARGS ANVI'O TRIED TO FILL IN BUT FAILED :(", lc='yellow')
            for variable, reason in self.__args_failed:
                self.run.info(variable, reason, nl_after= (1 if (variable, reason) == self.__args_failed[-1] else 0), lc="yellow")

        return self.args
