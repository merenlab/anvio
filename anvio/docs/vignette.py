"""Generate the existing programs vignette from CLI help output."""

import os

import anvio
import anvio.filesnpaths as filesnpaths
import anvio.terminal as terminal
import anvio.utils as utils

from anvio.programs import AnvioPrograms, parse_help_output
from anvio.summaryhtml import SummaryHTMLOutput


class ProgramsVignette:
    def __init__(self, args, r=terminal.Run(), p=terminal.Progress()):
        self.args = args
        self.run = r
        self.progress = p

        self.programs_to_skip = ['anvi-script-gen-programs-vignette']

        self.anvio_programs = AnvioPrograms(args, r=self.run, p=self.progress)

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.output_file_path = A("output_file")


    def generate(self):
        self.anvio_programs.init_programs(okay_if_no_meta = True, quiet = True)
        programs = self.anvio_programs.programs

        d = {}
        log_file = filesnpaths.get_temp_file_path()
        for i, program_name in enumerate(programs):
            program = programs[program_name]

            if program_name in self.programs_to_skip:
                self.run.warning("Someone doesn't want %s to be in the output :/ Fine. Skipping." % (program.name))

            self.progress.new('Bleep bloop')
            self.progress.update('%s (%d of %d)' % (program_name, i+1, len(programs)))

            output = utils.run_command_STDIN('%s --help --quiet' % (program.program_path), log_file, '').split('\n')

            if anvio.DEBUG:
                usage, params, output = parse_help_output(output)
            else:
                try:
                    usage, params, output = parse_help_output(output)
                except Exception as e:
                    self.progress.end()
                    self.run.warning("The program '%s' does not seem to have the expected help menu output. Skipping to the next. "
                                "For the curious, this was the error message: '%s'" % (program.name, str(e).strip()))
                    continue

            d[program.name] = {'usage': usage,
                               'description': program.meta_info['description']['value'],
                               'params': params,
                               'tags': program.meta_info['tags']['value'],
                               'resources': program.meta_info['resources']['value']}

            self.progress.end()

        os.remove(log_file)

        # generate output
        program_names = sorted([p for p in d if not p.startswith('anvi-script-')])
        script_names = sorted([p for p in d if p.startswith('anvi-script-')])
        vignette = {'vignette': d,
                    'program_names': program_names,
                    'script_names': script_names,
                    'all_names': program_names + script_names,
                    'meta': {'summary_type': 'vignette',
                             'version': '\n'.join(['|%s|%s|' % (t[0], t[1]) for t in anvio.get_version_tuples()]),
                             'date': utils.get_date()}}

        if anvio.DEBUG:
            self.run.warning(None, 'THE OUTPUT DICT')
            import json
            print(json.dumps(d, indent=2))

        open(self.output_file_path, 'w').write(SummaryHTMLOutput(vignette, r=self.run, p=self.progress).render())

        self.run.info('Output file', os.path.abspath(self.output_file_path))
