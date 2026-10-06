"""Generate the existing programs vignette from CLI help output."""

import json
import os
from argparse import Namespace

import anvio
import anvio.filesnpaths as filesnpaths
import anvio.terminal as terminal
import anvio.utils as utils

from anvio.programs import AnvioPrograms, parse_help_output
from anvio.summaryhtml import SummaryHTMLOutput


run = terminal.Run()
progress = terminal.Progress()


class ProgramsVignette(AnvioPrograms):
    def __init__(
        self,
        args: Namespace,
        r: terminal.Run = terminal.Run(),
        p: terminal.Progress = terminal.Progress(),
    ) -> None:
        self.args = args
        self.run = r
        self.progress = p

        self.programs_to_skip = ["anvi-script-gen-programs-vignette"]

        AnvioPrograms.__init__(self, args, r=self.run, p=self.progress)

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.output_file_path = A("output_file")

    def generate(self) -> None:
        self.init_programs(okay_if_no_meta=True, quiet=True)

        d = {}
        log_file = filesnpaths.get_temp_file_path()
        for i, program_name in enumerate(self.programs):
            program = self.programs[program_name]

            if program_name in self.programs_to_skip:
                run.warning(
                    f"Someone doesn't want {program.name} to be in the output :/ Fine. Skipping."
                )

            progress.new("Bleep bloop")
            progress.update(f"{program_name} ({i + 1} of {len(self.programs)})")

            output = utils.run_command_STDIN(
                f"{program.program_path} --help --quiet", log_file, ""
            ).split("\n")

            if anvio.DEBUG:
                usage, params, output = parse_help_output(output)
            else:
                try:
                    usage, params, output = parse_help_output(output)
                except Exception as e:
                    progress.end()
                    run.warning(
                        f"The program '{program.name}' does not seem to have the expected help menu output. Skipping to the next. "
                        f"For the curious, this was the error message: '{str(e).strip()}'"
                    )
                    continue

            d[program.name] = {
                "usage": usage,
                "description": program.meta_info["description"]["value"],
                "params": params,
                "tags": program.meta_info["tags"]["value"],
                "resources": program.meta_info["resources"]["value"],
            }

            progress.end()

        os.remove(log_file)

        # generate output
        program_names = sorted([p for p in d if not p.startswith("anvi-script-")])
        script_names = sorted([p for p in d if p.startswith("anvi-script-")])
        vignette = {
            "vignette": d,
            "program_names": program_names,
            "script_names": script_names,
            "all_names": program_names + script_names,
            "meta": {
                "summary_type": "vignette",
                "version": "\n".join(
                    f"|{name}|{version}|"
                    for name, version in anvio.get_version_tuples()
                ),
                "date": utils.get_date(),
            },
        }

        if anvio.DEBUG:
            run.warning(None, "THE OUTPUT DICT")
            print(json.dumps(d, indent=2))

        open(self.output_file_path, "w").write(
            SummaryHTMLOutput(vignette, r=run, p=progress).render()
        )

        run.info("Output file", os.path.abspath(self.output_file_path))
