# -*- coding: utf-8
# pylint: disable=line-too-long
"""Overloading Python argparse for anvi'o purposes"""

import os
import sys
import argparse
import textwrap

from colored import fore, back, style
from colored import fg, bg, attr

import anvio
import anvio.docs as docs


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


atty = sys.stdout.isatty()


class ArgumentParser(argparse.ArgumentParser):
    def __init__(self, description="No description :/", epilog=None):
        super().__init__()

        self.description = description
        self.epilog = epilog or self.get_anvio_epilogue()


    def get_anvio_epilogue(self):
        version = "main" if anvio.anvio_version.endswith('master') else anvio.anvio_version

        general_help = f"https://merenlab.org/software/anvio/help/{version}"
        program_help = f"{general_help}/{self.prog}"
        separator = '━' * 80 + '\n'

        if os.path.exists(os.path.join(os.path.dirname(docs.__file__), f"programs/{self.prog}.md")):
            epilog = textwrap.dedent(f'''
                 🔥 Find more on {self.prog} here:

                    {program_help}
            ''')
        else:
            epilog = ""


        epilog += textwrap.dedent(f'''
             🌊 All anvi'o programs and artifacts:

                {general_help}
        ''')


        if atty:
            return separator + attr('bold') + epilog + attr('reset')
        else:
            return separator + epilog


    def format_help(self):
        usage_formatter = argparse.ArgumentDefaultsHelpFormatter(self.prog)

        # usage
        usage_formatter.add_usage(self.usage, self._actions, self._mutually_exclusive_groups)

        # description
        if atty:
            description_text =  attr('bold') + self.description + attr('reset')
        else:
            description_text =  self.description

        usage_formatter.add_text(description_text)

        # positionals, optionals and user-defined groups
        for action_group in self._action_groups:
            if atty:
                section_title = action_group.title + ' ' * (80 - len(action_group.title) if len(action_group.title) < 80 else 0)
                section_header = bg(250) + fg(236) + attr('bold') + section_title + attr('reset')
            else:
                section_header = action_group.title

            usage_formatter.start_section(section_header)
            usage_formatter.add_text(action_group.description)
            usage_formatter.add_arguments(action_group._group_actions)
            usage_formatter.end_section()

        # epilog
        epilog_formatter = argparse.RawDescriptionHelpFormatter(prog=self.prog)
        epilog_formatter.add_text(self.epilog)

        # determine help from format above
        if atty:
            help_text = usage_formatter.format_help().replace(":\n", "\n") + "\n" + epilog_formatter.format_help()
        else:
            help_text = usage_formatter.format_help() + "\n" + epilog_formatter.format_help()

        return help_text + "\n"


    def get_args(self, parser):
        """A helper function to parse args anvi'o way.

           This function allows us to make sure some ad hoc parameters such as `--debug`
           can be used with any anvi'o program spontaneously even if they are not explicitly
           defined as an accepted argument, yet flags (or parameters) anvi'o does not expect
           to see can still be sorted out.
        """

        allowed_ad_hoc_flags = ['--version', '--debug', '--force', '--fix-sad-tables', '--quiet', '--no-progress', '--as-markdown']

        args, unknown = parser.parse_known_args()

        # if there are any args in the unknown that we do not expect to find
        # we we will make argparse complain about those.
        if len([f for f in unknown if f not in allowed_ad_hoc_flags]):
            for f in allowed_ad_hoc_flags:
                parser.add_argument(f, action='store_true')
            parser.parse_args()

        return args



