The Distribution of Coverage (DisCov) score is designed to capture the presence-absence of a contig or genome ('bin') in a metagenomic sample, as an alternative to using arbitrary detection thresholds for this purpose. DisCov quantifies how well a sequence is covered by reads — not just how many bases have coverage, but whether that coverage is spread across the sequence and whether it is roughly even in depth. To achieve this, DisCov incorporates two independent metrics:

* **_S_** (spread): the proportion of non-overlapping windows across the sequence that include bases with at least 1X coverage. In other words, `S = # of windows with coverage / # of windows`.
* **_E_** (evenness): the proportion of covered bases whose coverage falls within an expected fold-range of the median nonzero coverage. In other words, `E = # covered bases that fall within some range of the median / # of covered bases`

These are combined into a single score using either a linear or geometric formula:

* **Linear**: DisCov = _αS_ + (1-_α_)_E_
* **Geometric**: DisCov = _S_^_α_ × _E_^(1-_α_)

where _α_ is a value between 0 and 1 that controls how much weight is given to spread (_S_) relative to evenness (_E_). By default _α_ = 0.5 (equal weighting) and the linear formula is used, so DisCov is simply the average of _S_ and _E_.

When DisCov values are reported -- for instance, in %(bam-stats-txt)s -- usually _S_, _E_, and the number of windows are individually reported in addition to the overall DisCov score.

## Technical details & parameters

Here you'll find a description on how DisCov is computed. You can modify most aspects of the DisCov calculation, and the parameters for each step are mentioned when relevant below.

### Programs that produce these stats

%(anvi-profile-blitz)s always computes DisCov scores as part of its output. %(anvi-summarize)s computes DisCov scores when the `--report-discov` flag is provided.

### Input data

_S_, _E_, and DisCov scores are computed based on per-nucleotide coverage arrays for individual contigs. When computing scores at the genome-level, the arrays for all contigs within a given genome (or more generically, a given bin) are concatenated. Note that the concatentation of per-contig coverage arrays means that the windows for the _S_ metric can span across contig boundaries within the genome.

The per-base coverage arrays are either obtained directly from a %(bam-file)s, or from the auxiliary database files (`AUXILIARY-DATA.db`) produced by %(anvi-profile)s.

### _S_ and window length

We divide each sequence (contig or genome) into non-overlapping windows to compute _S_. _S_ is therefore essentially a more coarse-grained version of the regular detection metric (proportion of bases with coverage), and the window size affects its sensitivity. We want the windows to be long enough to capture meaningful clusters of mapped reads (hence, longer than typical short-read lengths), but short enough to effectively assess the distribution of these mapped reads along the sequence.

![Schematic for Spread](../../images/discov-spread.png){:.center-img .width-50}

There are two strategies for setting the window length:

* **Fixed window length** (`--window-length`): a constant window size in bp applied to all sequences.
* **Percentage-based window length** (`--window-length-as-percentage`): the window size is set as a percentage of each sequence's length, optionally with a minimum length floor (`--min-window-length`).

If you don't explicitly set any of these parameters, the default sizing strategy depends on the context. For genome-level stats, the default is a fixed 1,000 bp window. For contig-level stats, a percentage-based window (1%% of contig length, with a minimum of 300 bp) is used to better accommodate contigs of varying sizes.

{:.warning}
If a contig (or genome) is _shorter_ than the window length, the sequence will not be divided; that is, the entire sequence will be a single window. In this case, the _S_ metric collapses to a binary value (does the sequence have any coverage, or not), and from there: (1) If the sequence has no coverage at all (_S_ = 0), DisCov will be 0 (as the evenness metric _E_ will also be 0). (2) If the sequence does have coverage (_S_ = 1), then the DisCov value will have a minimum value of _α_ and the increase from _α_ is determined by _E_: DisCov = _α_ + (1-_α_)_E_. DisCov may not be very meaningful for these short sequences. In most cases, anvi'o will warn you if sequences are shorter than the computed window length, and regardless you should be able to see the number of windows per sequence in the output files.

**Removal of very short 'remainder' windows**

Since the windows for computing _S_ are non-overlapping, the final window in the sequence is often shorter than the other windows (its length is the remainder when you divide the total sequence length by the window size). In general it is not a big problem to include one shorter window in the calculation, but it means that these final bases can sometimes have very disproportional representation in the spread metric. So when the remainder is _very small_ (<10%% of the window size), we exclude the final, shorter window entirely from the calculation of _S_.

What happens if we using a fixed window length and are dealing with very short contigs that by themselves are <10%% of that length? We don't throw anything away, and those contigs will have exactly 1 window.

### _E_ and the fold-range

Evenness is computed only based upon the bases with nonzero coverage in the input sequence; that is, we filter the coverage array to keep only the nonzero values. From these nonzero coverage values, we compute the median coverage depth and multiply this value by the fold-range lower and upper bound to establish a range of 'expected' coverage depth. Bases with coverage depth within this range are counted as "evenly covered".

![Schematic for Evenness](../../images/discov-evenness.png){:.center-img .width-50}

You can adjust the fold-range to accommodate more or less variation in coverage depth by setting `--foldrange-lower` and `--foldrange-upper`. Wider fold-range values are more permissive; narrower ones are stricter. The default lower bound is 0.5x the median nonzero coverage while the default upper bound is 2.0x the median nonzero coverage. Note that these values are not required to be inverses of each other; for example, you could set a lower bound of 0.2x and an upper bound of 3x.

{:.warning}
The reliance on nonzero coverage values means that _E_ can have extremely high values when few reads map to a sequence. Even if an input sequence is 1,000,000 bp long, if it has only one read of 300 bp mapping to it, _E_ will be computed based upon the subset of ~300 bp where that read maps. In this example, there is an even 1x coverage on every single base of the ~300bp, and therefore _E_ = 1. With the linear DisCov formula and the default _α_ of 0.5, the DisCov score in a situation like this would be 0.5 despite the relative lack of read mapping, due to the high evenness. If that seems like too high of a score for this situation, then you can adjust _α_ or use the geometric formula for DisCov instead to penalize cases where _S_ is very low.

### The weight parameter _α_

You can adjust the weights on _S_ and _E_ by setting the `--alpha` parameter. Increase _α_ to increase the importance of spread, and decrease _α_ to increase the importance of evenness. _α_ should always be between 0 and 1 (inclusive. In case you really want DisCov to be equivalent to one of _S_ or _E_ you can set _α_ = 1 or _α_ = 0, respectively).

### Choosing the DisCov formula

The parameter `--discov-formula` controls which equation is used to combine _S_ and _E_ into the final DisCov score. If it's `linear` (the default), we take a weighted average of _S_ and _E_, and if it's `geometric`, we compute a weighted geometric mean of the two.

Each formulation has pros and cons. The linear version is typically easier to interpret -- high numbers mean at least one of _S_ and _E_ is high and this often is a good evidence for presence in a metagenome sample. However, the exceptions to this are the edge cases mentioned in the warning boxes above -- _S_ will be 1.0 for very short input sequences with some coverage while _E_ will be 1.0 when a few isolated reads map to the sequence. Those cases yield mid-level DisCov scores even though the evidence for their presence is not very strong, and that means that linear DisCov values close to ~0.5 are rather ambiguous (they might result from similar yet low-ish values of _S_ and _E_, or they might result from these extreme edge cases where one takes a value near 1.0 and the other takes a value near 0.0). The alternative geometric formulation requires input sequences to have both high _S_ and high _E_ to yield a high DisCov score and can therefore more clearly distinguish ambiguous cases.

Still not sure which formulation is best for your data? You could always compute both and maybe even combine them :)

## Okay, but how do I use DisCov?

You can use DisCov to decide if a given sequence is present or absent in a metagenome sample by setting a threshold score above which it is considered present. Where you draw the line depends on how confident you want to be (or how much you want to avoid false positives). A rough guideline based on insights during testing: for [linear DisCov](https://github.com/merenlab/anvio/pull/2563#issuecomment-4422260917), don't use a threshold any lower than ~0.6, and for [geometric DisCov](https://github.com/merenlab/anvio/pull/2563#issuecomment-4544106000), don't use any threshold lower than ~0.5.

## More resources

- The [DisCov Pull Request](https://github.com/merenlab/anvio/pull/2563) explains some implementation decisions, shows results from the hyperparameter tuning used to set the default parameters, and also describes the outcome of various validation tests.
