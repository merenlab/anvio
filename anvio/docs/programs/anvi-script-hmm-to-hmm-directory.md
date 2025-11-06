This script converts a `.hmm` file as formatted by HMMER3 `hmmbuild` to an anvi'o compatible [hmm-source](https://anvio.org/help/main/artifacts/hmm-source/).

### Specifying the input HMM

To start you must point the command to a list of `.hmm` files you'd like to include in your file `hmm-source` using the `--hmm-list` parameter. You may specify one or more files

{{ codestart }}
anvi-script-hmm-to-hmm-directory --hmm-list FILE1.hmm FILE2.hmm
{{ codestop }}

### Specifying the input HMM source

Then you must specify a source for the HMM using `--hmm-source`. Here you have two options. If you specify exactly one source, that will be applied to all the `.hmm` files passed via `--hmm-list`.

{{ codestart }}
anvi-script-hmm-to-hmm-directory --hmm-list FILE1.hmm FILE2.hmm --hmm-source COMMON_SOURCE
{{ codestop }}

Otherwise you must specify one source per file passed to `--hmm-list`.
{{ codestart }}
anvi-script-hmm-to-hmm-directory --hmm-list FILE1.hmm FILE2.hmm --hmm-source SOURCE1 SOURCE2
{{ codestop }}

### Specifying the output directory
You must specify the output directory you'd like anvi'o to save your `hmm-source` into using `-o` or `--output-directory`. Anvi'o will create a folder at the path you specify if it does not exist, otherwise it will exit.

{{ codestart }}
anvi-script-hmm-to-hmm-directory -o path/to/output
{{ codestop }}


### Basic usage

Putting everything together, a basic command would look like this:

{{ codestart }}
anvi-script-hmm-to-hmm-directory --hmm-list FILE1.hmm --hmm-source SOURCE1 -o output_folder
{{ codestop }}