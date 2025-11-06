This script downloads all the HMM models made public by the MDM-Paris lab in their [public repository](https://github.com/mdmparis/defense-finder-models) and converts it to an anvi'o compatible [hmm-source](https://anvio.org/help/main/artifacts/hmm-source/).

### Basic usage

The command will create the `hmm-source` in a folder called `DefenseFinder_HMM` in the current directory, if no such folder exists.

{{ codestart }}
anvi-script-gen-defense-finder-models-to-hmm-directory
{{ codestop }}