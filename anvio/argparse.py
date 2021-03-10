# -*- coding: utf-8
# pylint: disable=line-too-long
"""Overloading Python argparse for anvi'o purposes"""

import os
import sys
import argparse
import textwrap

from colored import fg, bg, attr

import anvio
import anvio.docs as docs

from anvio.programs import Program
from anvio.utils import is_program_exists as get_program_path


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
        """Function that formats the additional message that appears at the end of help."""

        version = anvio.anvio_version_for_help_docs

        general_help = f"https://merenlab.org/software/anvio/help/{version}"
        program_help = f"{general_help}/programs/{self.prog}"

        # starting with the requires / provides statements
        epilog = self.get_requires_and_provides_statements_for_program()

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

        program = Program(get_program_path(self.prog))
        requires = [v.id for v in program.meta_info['requires']['value']]
        provides = [v.id for v in program.meta_info['provides']['value']]

        def get_block(statement, header):
            blocks = []
            if len(requires):
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

            for i in range(0, len(blocks)):
                block = blocks[i]
                if atty:
                    block = [f"{fg('red') + b + attr('reset')}" if not b.strip() == '/' else b for b in block]
                else:
                    pass

                block.insert(0, '   ')

                blocks[i] = block

            if atty:
                blocks.insert(0, attr('bold') + header + attr('reset'))
            else:
                blocks.insert(0, header)

            blocks.insert(1, "")

            blocks = [''.join(b) for b in blocks]

            return '\n'.join(blocks)


        requires_and_provides_statements.append(get_block(requires, """🧀 Can consume: """))
        requires_and_provides_statements.append(get_block(provides, """🍕 Can provide: """))

        return '\n' + '\n'.join(requires_and_provides_statements)


    def format_help(self):
        """Individual formatting of sections in the help text.

        When we use the same formatter for all, we either would lose the
        explicit spacing in the epilog, or lose the formatting in other
        sections. In this function we change the formatters, render
        different sections differently, and then concatenate everything
        into a single output.
        """

        # we get our formatters here, fill them up down bleow, and finally render them at the end:
        usage_formatter = argparse.ArgumentDefaultsHelpFormatter(self.prog)
        description_formatter = argparse.RawDescriptionHelpFormatter(self.prog)
        epilog_formatter = argparse.RawDescriptionHelpFormatter(prog=self.prog)
        separator_formatter = argparse.RawDescriptionHelpFormatter(prog=self.prog)

        # usage
        usage_formatter.add_usage(self.usage, self._actions, self._mutually_exclusive_groups)

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


    def get_args(self, parser):
        """A helper function to parse args anvi'o way.

        This function allows us to make sure some ad hoc parameters, such as `--debug`,
        can be used with any anvi'o program spontaneously even if they are not explicitly
        defined as an accepted argument, yet flags (or parameters) anvi'o does not expect
        to see can still be sorted out.
        """

        allowed_ad_hoc_flags = ['--version', '--debug', '--force', '--fix-sad-tables', '--quiet', '--no-progress', '--as-markdown', '--tmp-dir']

        args, unknown = parser.parse_known_args()

        # if there are any args in the unknown that we do not expect to find
        # we we will make argparse complain about those.
        if len([f for f in unknown if f not in allowed_ad_hoc_flags]):
            for f in allowed_ad_hoc_flags:
                # handle non-boolean flags
                if f in ['--tmp-dir']:
                    parser.add_argument(f)
                else:
                    parser.add_argument(f, action='store_true')
            parser.parse_args()

        return args
