This program computes the average nucleotide identity between reads in a single fasta file (using PyANI).

To compute the ANI (or other genome distance metrics) between two genomes in different fasta files, use %(anvi-compute-genome-similarity)s.

A default run of this program looks like this:

{{ codestart }}
anvi-script-compute-ani-for-fasta -f %(fasta)s \
                                  -o path/to/output \
                                  --method ANIb
{{ codestop }}

By default, anvi'o uses the pyANI-plus backend with ANIb (which aligns 1020 nt fragments of your sequences using BLASTN+). ANIm is also supported. pyANI-plus is installed as a required anvi'o dependency. Use `--ani-backend legacy` to select the legacy pyANI executable instead. ANIblastall and TETRA are retired in the Python 3.13 port; their results are not silently mapped to another method.

Use `--pyani-plus-program PATH` to select a specific pyANI-plus executable. This option is unavailable with `--ani-backend legacy`.

The pyANI-plus backend currently supports Linux only; it is unsupported on macOS and Windows. On Linux, anvi'o allocates the effective `--num-threads` count to pyANI-plus child processes with `taskset` (default 1, or `ANVIO_THREADS` when set). The pyANI-plus local workflow can use all CPUs visible within that allocation; requests above the CPUs currently available to anvi'o are capped and the actual CPU IDs are reported. This controls child-process CPU affinity, not every OS thread. The allocation requires Linux, `os.sched_getaffinity`, and `taskset`. The legacy backend keeps its existing thread behavior.

You also have the option to change the distance metric (from the default "euclidean") or the linkage method (from the default "ward") or provide a path to a log file for debug messages.
