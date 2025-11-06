# -*- coding: utf-8
# A library to search anvi'o programs based on keywords or input/output files
# see `anvi-help` for its default client

import sys
import textwrap
import pandas as pd

from itertools import groupby
from operator import itemgetter
from colored import Fore, Back, Style

import anvio
import anvio.programs as p
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.terminal import tabulate
from anvio.errors import ConfigError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['ekiefl']


run = terminal.Run()
progress = terminal.Progress()


TABLE_BG = Back.GREY_11
TABLE_FG = Fore.GREY_100
TABLE_STYLE = TABLE_BG + TABLE_FG
HIGHLIGHT_BG = Back.GREY_11
HIGHLIGHT_FG = Fore.GREEN
HIGHLIGHT_STYLE = HIGHLIGHT_BG + HIGHLIGHT_FG
COL_WIDTH = 20 # Does not apply to the description column, which is extended if possible to fill terminal width

T = lambda x, l=COL_WIDTH: textwrap.TextWrapper(width=l, break_on_hyphens=True, break_long_words=True).fill(text=x)
F = lambda x, l=COL_WIDTH: textwrap.TextWrapper(width=l, break_on_hyphens=False, break_long_words=False).fill(text=x)
D = lambda x, l=COL_WIDTH: textwrap.TextWrapper(width=l, break_on_hyphens=True, break_long_words=True).fill(text=x).replace('<>', '\n')

class ProgramSearch:
    def __init__(self, args):
        A = lambda x, t: t(args.__dict__[x]) if x in args.__dict__ else None
        null = lambda x: x
        self.search_term_user_input = A('search-term', null)
        self.requires = A('requires', bool)
        self.provides = A('provides', bool)
        self.name = A('name', bool)
        self.report = A('report', null)

        self.headers_to_report = []
        self.headers = ['Program', 'Description', 'Tags', 'Provides', 'Requires']
        self.table = pd.DataFrame([], columns = self.headers)

        self.handle_inputs()

        anvio_programs = p.AnvioPrograms(args)
        anvio_programs.init_programs(okay_if_no_meta=True, quiet=True)
        self.programs = anvio_programs.programs

        self.get_headers_to_report()


    def get_headers_to_report(self):
        if self.headers_to_report:
            return

        self.headers_to_report = ['Program', 'Description', 'Provides', 'Requires']


    def handle_inputs(self):
        if self.requires and self.provides:
            raise ConfigError("You can't provide --requires and --provides. Pick one.")
        if self.name and (self.requires or self.provides):
            raise ConfigError("You can't provide --name with --provides or --requires. Pick one.")

        self.search_terms = []
        if self.search_term_user_input:
            self.search_terms = [x.strip() for x in self.search_term_user_input.split(',')]

        for search_term in self.search_terms:
            if len(search_term) < 2:
                raise ConfigError("Okay, you can't search for something 1 character long. If you're "
                                  "trying to search for everything, use the search term 'ALL'")

        if self.report:
            self.headers_to_report = ['Program'] + self.report.split(',')
            for header in self.headers_to_report:
                if header not in self.headers:
                    raise ConfigError('%s isn\'t a valid option for --report. Here are your options (comma separate them): %s' \
                                           % (header, ', '.join(self.headers)))
        else:
            self.get_headers_to_report()


    def process(self):
        if sys.stdout.encoding.lower() != "utf-8":
            run.warning("Sorry :/ The encoding of your terminal is not UTF-8, and as a result you will "
                        "not see anything but gibberish from displaying the nice table anvi'o will "
                        "put together for you. Aborting mission.")
            return

        self.populate_table()
        self.search_table()
        self.display_table()


    def parse_description_from_program_call(self, program):
        log_file = filesnpaths.get_temp_file_path()
        output = utils.run_command_STDIN('%s --help' % (program.program_path), log_file, '').split('\n')

        try:
            _, description, _, _ = p.parse_help_output(output)
        except Exception:
            return ''

        return description.replace('\n', '')


    def populate_table(self):
        for program in self.programs.values():
            tags = program.meta_info['tags']['value']
            provides = [item.id for item in program.meta_info['provides']['value']]
            requires = [item.id for item in program.meta_info['requires']['value']]
            description = program.meta_info['description']['value']

            resources = program.meta_info['resources']['value']

            if not len(description):
                # Yes, it can be learned from parse_description_from_program_call,
                # No, it is not worth waiting ~1s for.
                pass

            # if there are resources listed for the program, add them at the end of the description
            if len(resources):
                resources_text = '<><>'.join(['%s: %s' % (r) for r in resources])

                description += "<><>RESOURCES:<><>%s<>" % (resources_text)

            row_info = [program.name, description, tags, provides, requires]
            row = dict(zip(self.headers, row_info))
            self.table = self.table.append(row, ignore_index = True)

        self.table.sort_values(by = 'Program').reset_index(inplace = True)


    def search_table(self):
        # The dataframe cells are converted to strings so that search terms can be searched within
        self.table_as_string = self.stringify_table(self.table, COL_WIDTH)

        if self.search_terms == ['ALL']:
            return

        indices_to_keep = []
        for idx, row in self.table_as_string.iterrows():
            match = False
            for header, string in row.items():
                if self.provides and header != 'Provides':
                    continue
                if self.requires and header != 'Requires':
                    continue
                if self.name and header != 'Program':
                    continue

                for search_term in self.search_terms:
                    if search_term in string:
                        indices_to_keep.append(idx)
                        match = True
                        break

                if match:
                    break

        self.table = self.table.loc[indices_to_keep, self.headers_to_report]
        self.table_as_string = self.table_as_string.loc[indices_to_keep, self.headers_to_report]


    def highlight_matches(self):
        # this approach handles overlapping search terms
        for i, row in self.table_as_string.iterrows():
            for j, cell in row.items():
                indices_to_color = []

                for search_term in self.search_terms:
                    start_index = cell.find(search_term)
                    if start_index > -1:
                        indices_to_color.extend(list(range(start_index, start_index + len(search_term))))

                indices_to_color = sorted(set(indices_to_color))
                blocks_to_color = []

                # https://stackoverflow.com/questions/2361945/detecting-consecutive-integers-in-a-list
                for k, g in groupby(enumerate(indices_to_color), lambda ix: ix[0] - ix[1]):
                    block = list(map(itemgetter(1), g))
                    blocks_to_color.append((block[0], block[-1]))

                new_cell = ''
                last_block_end = 0
                for block_start, block_end in blocks_to_color:
                    new_cell += cell[last_block_end:block_start]
                    new_cell += HIGHLIGHT_STYLE + cell[block_start:block_end+1] + TABLE_STYLE

                    last_block_end = block_end + 1
                new_cell += cell[last_block_end:]

                self.table_as_string.loc[i, j] = new_cell


    def display_table(self):
        self.highlight_matches()
        tabulated = tabulate(self.table_as_string, headers='keys', showindex=False, tablefmt="fancy_grid")

        terminal_width, _ = terminal.get_terminal_size()
        table_width = len(tabulated.__str__().split('\n')[0])
        diff = terminal_width - table_width

        if diff < 0:
            run.warning("The above table is wider than your window, so it looks weird. Try making "
                        "your window bigger, zooming out (press - while holding command), or "
                        "specifying fewer columns to report with --report.")
        else:
            # Use all of the terminal width
            if 'Description' in self.table.columns:
                self.table_as_string = self.stringify_table(self.table, COL_WIDTH + diff)
                self.highlight_matches()
                tabulated = tabulate(self.table_as_string, headers='keys', showindex=False, tablefmt="fancy_grid")

        print('\n' + TABLE_STYLE + tabulated + Style.RESET)

        run.info('Number of programs matching', self.table.shape[0], nl_before=1)


    def stringify_table(self, table, description_width=None):
        printable_table = self.table.copy(deep = True)

        for idx, row in printable_table.iterrows():
            printable_table.loc[idx, :] = self.make_row_printable(row, description_width)

        return printable_table


    def make_row_printable(self, row, description_width=None):
        for header, info in row.iteritems():
            if header == 'Program':
                formatted_info = F(info) # never linebreaks the program name

            elif header == 'Description':
                formatted_info = D(info, description_width) if description_width else D(info)

            elif type(info) == list:
                formatted_info = '\n\n'.join([T(str(item)) for item in info])

            elif type(info) == str:
                formatted_info = T(info)

            else:
                formatted_info = info

            row.loc[header] = formatted_info

        return row

