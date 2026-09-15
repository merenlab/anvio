"""Downloading sequencing reads from the SRA as part of a workflow, without keeping them.

A workflow that maps public metagenomes against reference sequences does not need those
metagenomes afterwards, and on a machine with finite disk space it cannot afford to keep them:
a thousand metagenomes will not fit anywhere, no matter how patient you are. This module lets a
workflow download reads as it goes and hand them straight to snakemake's `temp()` machinery, so
each one is deleted the moment the last step that needed it is done.

Deleting reads promptly is only half the problem, though. Snakemake schedules whatever it can,
so left to its own devices it will happily download all thousand runs long before the first
mapping job finishes. Limiting how many downloads run *at once* does not help either, because a
download slot is released when the download ends, not when the reads are finally consumed. What
does work is making the download of a later sample depend on an earlier sample's reads having
been released: a sliding window built into the workflow graph itself. That is what the "release
unit" and "gate" machinery below is for.

A release unit is the set of samples that have to be on disk at the same time. Usually that is
one sample. When samples are co-assembled it is the whole co-assembly group, since an assembler
needs all of its input at once.
"""

import os
import math

import anvio.sra as sra
import anvio.terminal as terminal

from anvio.errors import ConfigError
from anvio.workflows import get_valid_lr_technologies


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__authors__ = ['FlorianTrigodet']


run = terminal.Run()

# The name of the release unit used when every sample has to be on disk simultaneously. It
# contains a character that anvi'o does not allow in sample or group names, so it can never be
# mistaken for one of them.
EVERY_SAMPLE_AT_ONCE = 'every-sample-at-once'


def suggested_disk_budget(gb):
    """Turn a predicted size into a `max_disk_gb` value that will actually hold it.

    A message that says "raise your budget to at least this much" has to round *up*: a suggestion
    that rounds down is still too small for what it is describing, and the user who follows it
    lands back on the very same error. Budgets are compared with a strict `>`, so a suggestion
    that falls exactly on the predicted size is accepted."""

    return f"{math.ceil(gb * 10) / 10:g}"


class SRAReadsModule:
    """Mixin that teaches a workflow to download its own reads from the SRA.

    This is a plain mixin in the spirit of `QCModule`: it contributes rules, directories and
    helpers, and it is deliberately not a workflow in its own right, so that including it never
    demands sra-tools from users who are not downloading anything."""

    def __init__(self):
        self.sra_mode = False
        self.sra_metadata = {}
        self.sra_accessions_by_sample = {}
        self.download_units = []
        self.download_unit_members = {}
        self.download_unit_peak_gb = {}
        self.download_unit_gates = {}
        self.max_disk_gb = None
        self.keep_reads = 'none'
        self.sra_metadata_cache_path = None
        self.derived_samples_txt = None
        self.derived_samples_txt_path = None

        self.rules.extend(['sra_prefetch', 'sra_fasterq_dump', 'sra_gather_reads'])

        self.dirs_dict.update({"SRA_DIR": "01_SRA"})


    def init_sra_reads(self, samples_txt):
        """Learn about the SRA runs in a samples-txt and predict where their reads will land.

        Returns a samples-txt with the `sra_accession` columns replaced by the paths of the FASTQ
        files this workflow is going to produce. Everything downstream — readsets, groups, QC,
        mapping — then works exactly as it does for reads that were already on disk, without
        knowing that anything was downloaded.

        If there are no accessions to speak of, the samples-txt is handed back untouched."""

        from anvio.artifacts.samples_txt import SamplesTxt

        self.sra_mode = samples_txt.has_sra_accessions()

        if not self.sra_mode:
            return samples_txt

        self.keep_reads = self.get_param_value_from_config(['download_reads', 'keep_reads']) or 'none'
        self.max_disk_gb = self.get_param_value_from_config(['download_reads', 'max_disk_gb'])

        self.sanity_check_sra_programs()

        if self.keep_reads not in ('none', 'raw', 'qc', 'both'):
            raise ConfigError(f"The `keep_reads` setting of the `download_reads` section of your config file must be "
                              f"one of 'none', 'raw', 'qc', or 'both'. You have '{self.keep_reads}' in there.")

        cache_path = self.get_param_value_from_config(['download_reads', 'metadata_cache']) or 'SRA-METADATA.txt'
        self.sra_metadata_cache_path = cache_path
        samples_txt_path = self.get_param_value_from_config(['samples_txt'])

        self.sra_metadata = sra.get_metadata_for_accessions(samples_txt.all_sra_accessions(),
                                                            cache_path=cache_path,
                                                            source_of_accessions=f"'{samples_txt_path}'",
                                                            run=self.run)

        sra.sanity_check_read_types_are_known(self.sra_metadata, cache_path)
        self.sanity_check_sra_read_types(samples_txt_path)
        self.sanity_check_samples_do_not_mix_sources(samples_txt, samples_txt_path)

        sra.warn_about_long_read_technologies(self.sra_metadata, run=self.run, cache_path=cache_path)

        if self.keep_reads != 'none':
            self.warn_about_kept_reads()

        if not self.max_disk_gb:
            self.run.warning(f"Anvi'o is going to download the reads of "
                             f"{terminal.pluralize('sample', len(samples_txt.samples()))} from the SRA, and you have "
                             f"not told it how much disk space it may use while doing so (that would be the "
                             f"`max_disk_gb` setting of the `download_reads` section of your config file). Without a "
                             f"budget anvi'o will download as fast as it can, which for a large set of samples can "
                             f"mean a great many of them sitting on your disk at the same time. Reads are still "
                             f"deleted as soon as nothing needs them anymore, so this is perfectly fine if you have "
                             f"the room.", header="NO DISK BUDGET SET", lc="yellow")

        derived = SamplesTxt(self.get_derived_samples_txt_data(samples_txt), expected_format="free",
                             skip_check_sanity=True, run=self.run)

        # Only the path is settled here. Writing it waits until the checks that can still refuse
        # this run have had their say (see `write_derived_samples_txt`).
        self.derived_samples_txt = derived
        self.derived_samples_txt_path = os.path.join(self.dirs_dict["SRA_DIR"], "samples-txt-with-downloaded-reads.txt")

        return derived


    def warn_about_kept_reads(self):
        """Say how much room the reads a user asked to keep are going to need.

        Kept reads are by definition not deleted as the workflow moves through the samples, so
        they pile up for the whole run. For a large set of samples that total is far larger than
        anything a disk budget ever allows on disk at one moment, and it is worth hearing before
        a thousand metagenomes are downloaded rather than after."""

        descriptions = {'raw': "the reads as they were downloaded",
                        'qc': "the quality-filtered reads",
                        'both': "both the reads as they were downloaded and the quality-filtered ones"}

        copies = 2 if self.keep_reads == 'both' else 1
        kept_gb = copies * sum(sra.predict_gzipped_fastq_bytes(e)
                               for e in self.sra_metadata.values()) / (1024 ** 3)

        budget_note = (" A disk budget cannot speak for them, either: `max_disk_gb` governs only the files "
                       "anvi'o is free to delete, and these are not among them." if self.max_disk_gb else "")

        self.run.warning(f"You have set `keep_reads` to '{self.keep_reads}', so anvi'o will hold on to "
                         f"{descriptions[self.keep_reads]} rather than deleting them once nothing needs them "
                         f"anymore. Unlike everything else this workflow downloads, these accumulate for the "
                         f"whole run: by the end of it anvi'o expects them to come to roughly "
                         f"{kept_gb:.1f} GB.{budget_note}",
                         header="THE READS YOU ARE KEEPING NEED ROOM OF THEIR OWN", lc="yellow")


    def sanity_check_sra_programs(self):
        """Make sure the tools that do the downloading are around, before anything else happens.

        These are only ever needed by users whose samples-txt names SRA accessions, which is why
        this is checked here rather than when the workflow class is instantiated."""

        import anvio.utils as utils

        missing = [p for p in ['prefetch', 'fasterq-dump', 'pigz'] if not utils.is_program_exists(p, dont_raise=True)]

        if missing:
            raise ConfigError(f"Your samples-txt file asks anvi'o to download reads from the SRA, but "
                              f"{terminal.pluralize('program', len(missing))} needed for that "
                              f"{'are' if len(missing) > 1 else 'is'} not installed in your anvi'o environment: "
                              f"{', '.join(missing)}. `prefetch` and `fasterq-dump` come with the NCBI SRA toolkit "
                              f"(https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit), and `pigz` is "
                              f"listed in the anvi'o installation instructions (https://anvio.org/install/).")


    def sanity_check_sra_read_types(self, samples_txt_path):
        """Refuse the kinds of reads this workflow cannot process, before anything is downloaded."""

        single_end = sorted(a for a, e in self.sra_metadata.items() if e['read_type'] == sra.SINGLE_END_SHORT_READS)

        if single_end:
            details = '\n'.join(f"    {a} ({self.sra_metadata[a]['platform']}, "
                                f"{self.sra_metadata[a]['model']})" for a in single_end[:10])
            raise ConfigError(f"{terminal.pluralize('accession', len(single_end))} in '{samples_txt_path}' "
                              f"{'are' if len(single_end) > 1 else 'is'} single-end short reads, and the anvi'o "
                              f"metagenomics workflow does not handle those at this point in time. Here "
                              f"{'they are' if len(single_end) > 1 else 'it is'}:\n\n{details}\n\n"
                              f"To go ahead with the rest of your samples, please remove "
                              f"{'these accessions' if len(single_end) > 1 else 'this accession'} from your "
                              f"samples-txt file. Anvi'o is sorry to be the bearer of this news, and would rather "
                              f"tell you now than after downloading everything else.")


    def sanity_check_samples_do_not_mix_sources(self, samples_txt, samples_txt_path):
        """Refuse a single sample that has both files on disk and accessions of the same kind.

        Different samples may perfectly well come from different places, and often do. But when
        one sample names both, there is no sensible reading of what it wants: anvi'o would have
        to either ignore the files or ignore the accessions, and both would be surprising."""

        offenders = []

        for sample in samples_txt.samples():
            info = samples_txt.get_sample(sample)
            accessions = samples_txt.sra_accessions_for_sample(sample)

            if not accessions:
                continue

            read_types = set(self.sra_metadata[a]['read_type'] for a in accessions)

            if info.get('r1') and sra.PAIRED_END_SHORT_READS in read_types:
                offenders.append((sample, 'short'))
            if info.get('lr') and sra.LONG_READS in read_types:
                offenders.append((sample, 'long'))

        if offenders:
            details = '\n'.join(f"    {s} already has {kind} reads of its own, and also names "
                                f"{kind}-read accessions" for s, kind in offenders)
            raise ConfigError(f"{terminal.pluralize('sample', len(offenders))} in '{samples_txt_path}' "
                              f"{'have' if len(offenders) > 1 else 'has'} both file paths and SRA accessions for "
                              f"the same kind of reads, and anvi'o does not know which of the two you meant:\n\n"
                              f"{details}\n\nPlease use one or the other for any given sample. Having some samples "
                              f"come from files and others from the SRA is fine; it is only mixing the two within "
                              f"one sample that anvi'o cannot make sense of.")


    def get_derived_samples_txt_data(self, samples_txt):
        """Build the samples-txt anvi'o will actually work with.

        Each sample's accessions are sorted by the kind of reads they hold, and every kind gets
        exactly one file: several sequencing runs of the same sample are concatenated into it.
        A sample with both short-read and long-read accessions therefore comes out of here
        looking precisely like a hybrid sample whose files were on disk all along."""

        derived = {}

        # The lr_technology column is all-or-nothing. Anvi'o can either say what every long-read
        # sample was sequenced with, or say nothing and leave the presets to the config file —
        # which is the right answer when every long-read sample came off the same machine, and
        # what `warn_about_long_read_technologies` has already pointed at. Knowing some but not
        # all of them is neither, and is refused rather than papered over.
        technologies = self.get_long_read_technologies(samples_txt)
        self.sanity_check_long_read_technologies(technologies)
        declare_technologies = bool(technologies) and all(technologies.values())

        for sample in samples_txt.samples():
            info = samples_txt.get_sample(sample)
            accessions = samples_txt.sra_accessions_for_sample(sample)

            row = {'group': info.get('group'),
                   'r1': list(info.get('r1') or []),
                   'r2': list(info.get('r2') or []),
                   'lr': list(info.get('lr') or [])}

            if accessions:
                self.sra_accessions_by_sample[sample] = accessions

                if self.get_accessions_for_sample(sample, sra.PAIRED_END_SHORT_READS):
                    row['r1'] = [self.get_reads_path(sample, 'R1')]
                    row['r2'] = [self.get_reads_path(sample, 'R2')]

                if self.get_accessions_for_sample(sample, sra.LONG_READS):
                    row['lr'] = [self.get_reads_path(sample, 'LR')]

            # The key is left out altogether when there is no column to declare. Its mere
            # presence would tell the rest of anvi'o that this samples-txt has an lr_technology
            # column, and the checks that make sure long-read presets are set would then skip
            # themselves, leaving a run to fail later and more obscurely than it needs to.
            if declare_technologies:
                row['lr_technology'] = technologies.get(sample)

            derived[sample] = row

        return derived


    def sanity_check_long_read_technologies(self, technologies):
        """Refuse a run where anvi'o knows some long-read technologies but not all of them.

        The two workable states are knowing all of them, in which case every long-read tool is
        given presets chosen per sample, and knowing none of them, in which case the config file
        decides for everybody. Half an answer is neither: anvi'o would have to either throw away
        what it does know or process one sample's reads with another sample's settings, and both
        of those are worse than saying so."""

        if not technologies or all(technologies.values()) or not any(technologies.values()):
            return

        known = sorted(f"{s} ({t})" for s, t in technologies.items() if t)
        unknown = sorted(s for s, t in technologies.items() if not t)

        raise ConfigError(f"Anvi'o knows which sequencing technology some of your long-read samples were produced "
                          f"with, but not all of them, and it cannot work with half an answer. Long-read mapping "
                          f"and assembly presets are chosen from that technology one sample at a time, so the "
                          f"{terminal.pluralize('sample', len(unknown))} anvi'o is unsure about would have to be "
                          f"processed with some other sample's settings — which is how nanopore reads end up being "
                          f"mapped as though they were PacBio HiFi.\n\n"
                          f"It knows about: {', '.join(known[:10])}{' (and more)' if len(known) > 10 else ''}.\n"
                          f"It does not know about: {', '.join(unknown[:10])}"
                          f"{' (and more)' if len(unknown) > 10 else ''}.\n\n"
                          f"Please fill in the gap. For a sample whose reads are files on your disk, add an "
                          f"`lr_technology` column to your samples-txt. For one anvi'o is downloading, you can do "
                          f"that too, or fill in the `lr_technology` column of "
                          f"{repr(self.sra_metadata_cache_path) if self.sra_metadata_cache_path else 'your SRA metadata file'} "
                          f"for the accessions in question. Valid values are "
                          f"{', '.join(repr(t) for t in sorted(get_valid_lr_technologies()))}.")


    def get_long_read_technologies(self, samples_txt):
        """Return {sample: technology or None} for every sample that will have long reads.

        What the user wrote in the samples-txt wins over what anvi'o worked out from NCBI's
        metadata, since the user is the one who knows."""

        technologies = {}

        for sample in samples_txt.samples():
            info = samples_txt.get_sample(sample)
            accessions = [a for a in samples_txt.sra_accessions_for_sample(sample)
                          if self.sra_metadata.get(a, {}).get('read_type') == sra.LONG_READS]

            if not accessions and not info.get('lr'):
                continue

            declared = info.get('lr_technology')
            inferred = set(self.sra_metadata[a]['lr_technology'] for a in accessions if self.sra_metadata[a]['lr_technology'])

            technologies[sample] = declared or (inferred.pop() if len(inferred) == 1 else None)

        return technologies


    def get_accessions_for_sample(self, sample, read_type):
        """Return a sample's accessions of one kind, in the order the user listed them."""

        return [a for a in self.sra_accessions_by_sample.get(sample, [])
                if self.sra_metadata[a]['read_type'] == read_type]


    def get_reads_path(self, sample, kind):
        """Where a sample's downloaded reads end up. `kind` is one of R1, R2 or LR.

        Everything lands in a single flat directory on purpose: `iu-gen-configs` insists that the
        R1 and R2 of a sample live side by side in a directory that already exists, and it is
        simpler to satisfy that than to work around it."""

        return os.path.join(self.dirs_dict["SRA_DIR"], "reads", f"{sample}_{kind}.fastq.gz")


    def get_accession_fastq_paths(self, accession):
        """The FASTQ files a single SRA run will be turned into.

        Which files these are is decided by what NCBI says the run holds, which is why the
        metadata has to be in hand before the workflow graph can be built."""

        directory = os.path.join(self.dirs_dict["SRA_DIR"], "runs", accession)
        read_type = self.sra_metadata[accession]['read_type']

        if read_type == sra.PAIRED_END_SHORT_READS:
            return [os.path.join(directory, f"{accession}_1.fastq.gz"),
                    os.path.join(directory, f"{accession}_2.fastq.gz")]

        return [os.path.join(directory, f"{accession}.fastq.gz")]


    def init_sra_download_units(self):
        """Work out how much can be downloaded at a time, and in what order.

        This has to happen after assembly groups are known, since co-assembled samples must be on
        disk together and therefore belong to the same release unit."""

        if not self.sra_mode:
            return

        self.sanity_check_sra_workflow_combination()

        # In references mode every sample stands alone: its reads are mapped and then finished
        # with. When samples are assembled, everything that goes into one assembly has to be
        # present at the same time, so the assembly group becomes the unit. And when every sample
        # is mapped against every assembly, there is no way to take them a few at a time at all,
        # so they all end up in a single unit (see `everything_must_be_downloaded_at_once`).
        for sample in self.samples_txt.samples():
            # A sample whose reads were already on disk has nothing to download, and therefore
            # nothing to release either. Release units are about downloaded reads alone, so it
            # does not belong to one. A co-assembly group that mixes the two still becomes a
            # unit, made up of the members that were actually downloaded: the assembly those
            # reads feed is among the unit's consumers either way.
            if not self.sra_accessions_by_sample.get(sample):
                continue

            if self.everything_must_be_downloaded_at_once():
                unit = EVERY_SAMPLE_AT_ONCE
            elif self.references_mode:
                unit = sample
            else:
                unit = self.samples_txt.get_sample(sample).get('group') or sample

            if unit not in self.download_unit_members:
                self.download_units.append(unit)
                self.download_unit_members[unit] = []

            self.download_unit_members[unit].append(sample)

        for unit, members in self.download_unit_members.items():
            self.download_unit_peak_gb[unit] = sum(self.get_predicted_peak_gb(s) for s in members)

        self.sanity_check_disk_budget()
        self.download_unit_gates = self.compute_download_gates()

        self.write_derived_samples_txt()


    def write_derived_samples_txt(self):
        """Put the samples-txt anvi'o actually works with on disk.

        This waits until everything that can refuse the run has had its say, so that a config
        anvi'o turns down leaves nothing behind. It still happens on a dry run: the file is an
        input of the rule that generates `iu-gen-configs` configs, so snakemake has to see it to
        build a graph at all.

        The directory the reads will land in is made here too, because `iu-gen-configs` insists
        that the directory holding a sample's files already exists, though the files need not."""

        os.makedirs(os.path.join(self.dirs_dict["SRA_DIR"], "reads"), exist_ok=True)

        self.derived_samples_txt.write_tsv(self.derived_samples_txt_path, absolute_paths=False,
                                           include_extras=False)


    def warn_about_qc_output_of_local_samples(self):
        """Say that quality-filtered reads go for every sample, not only the downloaded ones.

        Snakemake decides whether a rule's output is temporary when it reads the rule, not when it
        runs a job, so `temp()` cannot be applied to one sample's quality-filtered reads and not
        another's. In a samples-txt that mixes downloaded reads with reads that were already on
        disk, that means the filtered copies of the local samples are deleted too. Nothing the
        user supplied is touched, and the filtering can be redone from files that never went
        anywhere, but it is not what anyone would predict from a setting named `keep_reads`."""

        if not self.sra_mode or self.keep_reads in ('qc', 'both'):
            return

        local = [s for s in self.samples_txt.samples() if s not in self.sra_accessions_by_sample]
        affected = []

        # Reference-based read removal already makes quality-filtered short reads temporary on
        # its own and says so itself, so in a run that does it, short reads are not a surprise
        # that downloading brought with it. It is a short-read step and has no say over what
        # filtlong makes, so local long-read samples are one either way.
        if getattr(self, 'run_qc', False) and not getattr(self, 'remove_short_reads_based_on_references', None):
            affected.extend(s for s in local if self.samples_txt.get_sample(s).get('r1'))

        if getattr(self, 'run_filtlong', False):
            affected.extend(s for s in local if self.samples_txt.get_sample(s).get('lr'))

        local = sorted(set(affected))

        if not local:
            return

        self.run.warning(f"Your samples-txt mixes reads anvi'o is going to download with reads that were "
                         f"already on your disk, and `keep_reads` is set to '{self.keep_reads}'. Quality-"
                         f"filtered reads are deleted once nothing needs them anymore, and snakemake decides "
                         f"that per rule rather than per sample, so it applies to every sample in this run — "
                         f"the {terminal.pluralize('sample', len(local))} whose reads were yours to begin with "
                         f"included: {', '.join(local[:5])}{' (and more)' if len(local) > 5 else ''}.\n\n"
                         f"Nothing you supplied is touched. The FASTQ files you pointed anvi'o at stay exactly "
                         f"where they are, and it is only the quality-filtered copies anvi'o made of them that "
                         f"go. If you would rather keep those, set `keep_reads` to 'qc' or 'both' in the "
                         f"`download_reads` section of your config file — bearing in mind that this keeps the "
                         f"quality-filtered reads of the downloaded samples as well, which is the disk the "
                         f"budget was there to save.",
                         header="FILTERED READS GO FOR YOUR OWN SAMPLES TOO", lc="yellow")


    def everything_must_be_downloaded_at_once(self):
        """Whether this workflow leaves anvi'o no room to download samples a few at a time.

        Mapping every sample against every assembly is such a case. A sample's reads cannot be
        let go of until it has been mapped against the last assembly, and that assembly cannot
        exist until its own reads have been downloaded, so by the time the final assembly is
        built every sample is necessarily still on disk. Reads are still deleted afterwards —
        that part works no matter what — but there is no order in which they could be spread
        out. In references mode the same setting is harmless, because references come from a
        fasta-txt rather than from the reads."""

        return bool(self.get_param_value_from_config('all_against_all')) and not self.references_mode


    def get_predicted_peak_gb(self, sample):
        """The most disk space one sample's reads will take up while they are being dealt with."""

        safety_factor = self.get_param_value_from_config(['download_reads', 'safety_factor']) or 1.3

        return sum(sra.predict_peak_disk_usage_in_gb(self.sra_metadata[a], safety_factor=safety_factor)
                   for a in self.sra_accessions_by_sample.get(sample, []))


    def sanity_check_sra_workflow_combination(self):
        """Refuse the combinations of settings for which a disk budget cannot mean anything."""

        if self.everything_must_be_downloaded_at_once():
            self.run.warning("Your config file asks anvi'o to download reads from the SRA, to assemble them, and to "
                             "map every sample against every assembly (`all_against_all`). This is a perfectly "
                             "reasonable thing to want, but it does mean that every metagenome has to be on disk at "
                             "the same time: a sample cannot be mapped until every assembly exists, and an assembly "
                             "cannot be built until its own reads have been downloaded. So there is no order in "
                             "which anvi'o could download your samples a few at a time.\n\n"
                             "Anvi'o will still delete each sample's reads as soon as nothing needs them anymore, "
                             "which for this kind of run means once the last assembly has been mapped against. What "
                             "it cannot do here is spread the downloads out, so a disk budget can only tell you "
                             "whether you have enough room for all of them — it cannot make do with less.",
                             header="EVERY SAMPLE WILL BE ON DISK AT ONCE", lc="yellow")

        if self.remove_short_reads_based_on_references and not self.run_qc:
            raise ConfigError("Your config file asks anvi'o to download reads from the SRA and to remove short reads "
                              "based on reference sequences, but it also turns off quality filtering "
                              "(`iu_filter_quality_minoche`). In that combination the reference removal step "
                              "consumes the downloaded files themselves rather than a filtered copy of them, which "
                              "would leave anvi'o with nothing to give to the steps that come afterwards, and it "
                              "would have to download everything a second time. Please either leave quality "
                              "filtering on, or drop the reference removal step.")

        if self.additional_params and '--notemp' in self.additional_params and self.keep_reads == 'none':
            raise ConfigError("You passed `--notemp` to snakemake through `--additional-params`, which tells it to "
                              "keep every intermediate file it makes. Since the whole point of downloading reads "
                              "inside this workflow is that they are deleted once they have been used, these two "
                              "things do not go together. Please drop `--notemp`, or set `keep_reads` to something "
                              "other than 'none' in the `download_reads` section of your config file if you meant to "
                              "keep the reads.")


    def get_predicted_sizes_gb(self, unit):
        """The archive and uncompressed-FASTQ sizes a release unit's peak estimate is built from."""

        accessions = [a for sample in self.download_unit_members[unit]
                        for a in self.sra_accessions_by_sample.get(sample, [])]

        archive_gb = sum(sra.predict_archive_bytes(self.sra_metadata[a]) for a in accessions) / (1024 ** 3)
        fastq_gb = sum(sra.predict_fastq_bytes(self.sra_metadata[a]) for a in accessions) / (1024 ** 3)

        return archive_gb, fastq_gb


    def explain_disk_estimate(self, unit, subject):
        """Say where a unit's predicted size came from, in GB anyone can check against NCBI.

        These predictions are several times larger than the download sizes NCBI advertises, and
        the first thing anyone does with a budget error is compare the two and wonder which of
        them is wrong. Neither is: the difference is compression, and this explains that with the
        actual numbers rather than leaving it to be guessed at."""

        archive_gb, fastq_gb = self.get_predicted_sizes_gb(unit)
        safety_factor = self.get_param_value_from_config(['download_reads', 'safety_factor']) or 1.3

        return (f"If these numbers look far larger than the download sizes NCBI advertises, that is because an "
                f"SRA archive is compressed and this workflow works with uncompressed FASTQ. {subject} arrive as "
                f"about {archive_gb:.1f} GB of archives and unpack to about {fastq_gb:.1f} GB of FASTQ. What "
                f"anvi'o counts is the archive, the scratch space `fasterq-dump` needs while it unpacks it, and "
                f"the FASTQ itself, since all three are on disk at the same moment — and then multiplies the "
                f"total by `safety_factor` (currently {safety_factor}). The reads do not stay that large: they "
                f"are compressed as soon as they have been extracted.")


    def sanity_check_disk_budget(self):
        """Make sure the budget can actually fit what it is being asked to fit."""

        if not self.max_disk_gb:
            return

        too_big = {u: gb for u, gb in self.download_unit_peak_gb.items() if gb > self.max_disk_gb}

        if not too_big:
            return

        # When everything has to be downloaded at once there is only ever one unit, and saying
        # that "1 sample or group" does not fit would be a strange way to put it.
        if self.everything_must_be_downloaded_at_once():
            needed = self.download_unit_peak_gb[EVERY_SAMPLE_AT_ONCE]
            raise ConfigError(f"Your disk budget of {self.max_disk_gb} GB (`max_disk_gb` in the `download_reads` "
                              f"section of your config file) is not enough for this workflow. Because you asked for "
                              f"every sample to be mapped against every assembly (`all_against_all`), all "
                              f"{terminal.pluralize('downloaded sample', len(self.download_unit_members[EVERY_SAMPLE_AT_ONCE]))} "
                              f"have to be on disk "
                              f"at the same time, and anvi'o expects that to need about "
                              f"{suggested_disk_budget(needed)} GB. There is no "
                              f"way to do this run with less, so please either raise `max_disk_gb` to at least that "
                              f"much, or reconsider `all_against_all`: with it turned off, anvi'o can work through "
                              f"your samples a few at a time and fit into whatever budget you have.\n\n"
                              f"{self.explain_disk_estimate(EVERY_SAMPLE_AT_ONCE, 'Your samples')}")

        offenders = '\n'.join(f"    {u} needs about {gb:.1f} GB "
                              f"({terminal.pluralize('sample', len(self.download_unit_members[u]))}: "
                              f"{', '.join(self.download_unit_members[u][:5])})"
                              for u, gb in sorted(too_big.items(), key=lambda x: -x[1])[:10])

        biggest = max(too_big, key=lambda u: too_big[u])

        co_assembly_note = ("" if self.references_mode else
                            "\n\nSome of these are co-assembly groups rather than single samples. Everything that "
                            "goes into one assembly has to be on disk at the same time, so a group's samples count "
                            "together against your budget.")

        raise ConfigError(f"Your disk budget of {self.max_disk_gb} GB (`max_disk_gb` in the `download_reads` section "
                          f"of your config file) is smaller than what anvi'o will need to hold at one moment for "
                          f"{terminal.pluralize('sample or group', len(too_big))}, so this workflow could never "
                          f"finish:\n\n{offenders}\n\nPlease raise `max_disk_gb` to at least "
                          f"{suggested_disk_budget(max(too_big.values()))} GB, or leave the biggest of "
                          f"these out."
                          f"{co_assembly_note}\n\n"
                          f"{self.explain_disk_estimate(biggest, 'The reads of ' + biggest)}")


    def compute_download_gates(self):
        """Decide which earlier unit each unit has to wait for.

        Walking backwards from each unit and adding up predicted sizes tells us how far back we
        can go before the budget is used up. Everything from that point on may be in flight at
        the same time, so the unit just before it is the one to wait for. Units that fit inside
        the budget together with everything before them do not wait for anything at all."""

        gates = {}

        if not self.max_disk_gb:
            return gates

        peaks = [self.download_unit_peak_gb[u] for u in self.download_units]

        for k in range(len(self.download_units)):
            total, j = 0.0, k

            while j >= 0:
                total += peaks[j]
                if total > self.max_disk_gb:
                    break
                j -= 1

            # `j` is now the last unit that does not fit alongside this one, so this unit waits
            # for it. `sanity_check_disk_budget` has already ruled out a unit that cannot fit on
            # its own, so `j` can never be `k` and a unit can never end up waiting for itself.
            if j >= 0:
                gates[self.download_units[k]] = self.download_units[j]

        return gates


    def get_gate_flag_for_accession(self, accession):
        """The release flag a run's download has to wait for, if any.

        More than one sample can name the same accession. The run is downloaded once and its
        reads then stay until the last sample that wants them is done, so what governs when the
        download may start is the *earliest* unit that needs it. Waiting on any later one would
        be a promise anvi'o cannot keep: that unit's gate may well be the earlier unit itself,
        whose reads cannot be released until this very download has happened."""

        for unit in self.download_units:
            if any(accession in self.sra_accessions_by_sample.get(sample, [])
                   for sample in self.download_unit_members[unit]):
                gate = self.download_unit_gates.get(unit)

                return [self.get_release_flag_path(gate)] if gate else []

        return []


    def get_gate_flag_for_sample(self, sample):
        """The release flag a sample's download has to wait for, if any."""

        for unit, members in self.download_unit_members.items():
            if sample in members:
                gate = self.download_unit_gates.get(unit)
                return [self.get_release_flag_path(gate)] if gate else []

        return []


    def get_release_flag_path(self, unit):
        """Where the marker that says 'this unit's reads are done with' lives."""

        return os.path.join(self.dirs_dict["SRA_DIR"], "released", f"{unit}.released")


    def get_read_consumer_targets(self, unit):
        """Every file that has to exist before a unit's reads can be considered finished with.

        Anvi'o does not delete the reads itself — snakemake works out from the workflow graph
        which files are still needed and removes them the moment they are not, which is the only
        way to be sure nothing is deleted too early. This list decides something different and
        more forgiving: when the *next* download is allowed to start. If something is missing
        from it the disk budget gets looser than promised, but nothing breaks."""

        targets = []

        for sample in self.download_unit_members[unit]:
            for readset in [rs['id'] for rs in self.readsets if rs['base_sample'] == sample]:
                for group in self.get_groups_a_readset_is_mapped_to(readset):
                    targets.append(os.path.join(self.dirs_dict["PROFILE_DIR"], group, readset, "PROFILE.db"))

                if self.get_param_value_from_config(['import_percent_of_reads_mapped', 'run']) == True:
                    targets.append(os.path.join(self.dirs_dict["QC_DIR"], f"{readset}-total_num_reads.txt"))

                if self.run_krakenuniq and self.readsets_by_id[readset]['type'] == 'SR':
                    targets.append(os.path.join(self.dirs_dict["TAXONOMY_DIR"], f"{readset}-krakenuniq-report.tsv"))

                targets.extend(self.get_qc_report_targets_for_readset(readset))

        if not self.references_mode:
            # The assembly is the reason a co-assembly group's samples had to be on disk together
            # in the first place, so it belongs here too.
            for group in set(self.get_assembly_groups_for_unit(unit)):
                targets.append(os.path.join(self.dirs_dict["FASTA_DIR"], group, "final.contigs.fa"))

        return sorted(set(targets))


    def get_groups_a_readset_is_mapped_to(self, readset):
        """Which assemblies or references a readset is mapped against."""

        return [g for g in self.group_names if readset in self.get_readsets_for_mapping_to_group(g)]


    def get_assembly_groups_for_unit(self, unit):
        """The assembly groups that consume the reads of a release unit."""

        groups = []
        for group, members in self.assembly_members.items():
            for readset in members:
                if self.readsets_by_id[readset]['base_sample'] in self.download_unit_members[unit]:
                    groups.append(group)

        return groups


    def get_qc_report_targets_for_readset(self, readset):
        """The QC reports that read a readset's files, when those are turned on."""

        targets = []

        for parent_dir, readset_ids, stages in self.qc_producers():
            if readset in readset_ids:
                for stage in stages:
                    targets.append(os.path.join(parent_dir, readset, stage))

        return targets


    def get_accessions_of_read_type(self, read_type):
        """Every accession in this workflow that holds one kind of reads."""

        return [a for a, e in self.sra_metadata.items() if e['read_type'] == read_type]


    def get_gather_inputs(self, sample, kind):
        """The per-run FASTQ files that make up one of a sample's read files."""

        if kind == 'LR':
            accessions = self.get_accessions_for_sample(sample, sra.LONG_READS)
            return [self.get_accession_fastq_paths(a)[0] for a in accessions]

        accessions = self.get_accessions_for_sample(sample, sra.PAIRED_END_SHORT_READS)
        mate = 0 if kind == 'R1' else 1

        return [self.get_accession_fastq_paths(a)[mate] for a in accessions]


    def temp_unless_reads_are_kept(self, path, stage):
        """Wrap a path in snakemake's temp() unless the user asked to keep reads of this stage.

        `stage` is 'raw' for the reads as they were downloaded and 'qc' for what quality
        filtering makes of them."""

        from snakemake.io import temp

        if not self.sra_mode:
            return path

        return path if self.keep_reads in (stage, 'both') else temp(path)


    def find_sra_archive(self, directory):
        """Return the archive `prefetch` left in a directory, whichever form it took."""

        import glob

        archives = [f for f in glob.glob(os.path.join(directory, "**", "*"), recursive=True)
                    if f.endswith(('.sra', '.sralite'))]

        return archives[0] if archives else None


    def run_fasterq_dump(self, shell, accession, prefetch_dir, output_dir, threads, log_path, read_type):
        """Extract a run's reads from its archive and compress them.

        Which files should come out of this was decided before the workflow started, from what
        NCBI says the run contains. That metadata is written by whoever submitted the data and is
        occasionally at odds with the data itself, so what actually appeared is checked against
        what was expected and a mismatch is reported as such, rather than as a baffling missing
        output file."""

        import shutil

        archive = self.find_sra_archive(prefetch_dir)
        scratch_dir = os.path.join(output_dir, "scratch")

        if not archive:
            raise ConfigError(f"Anvi'o was about to extract the reads of the SRA run {accession} from its archive, "
                              f"and there is no archive in '{prefetch_dir}' to extract them from. Since `prefetch` "
                              f"checks for one before it reports success, the likeliest explanation is that "
                              f"something removed it afterwards. Deleting that directory and running the workflow "
                              f"again will download it afresh.")

        shell(f"fasterq-dump {archive} --outdir {output_dir} --temp {scratch_dir} --split-3 "
              f"--threads {threads} >> {log_path} 2>&1")

        paired = os.path.join(output_dir, f"{accession}_1.fastq")
        paired_mate = os.path.join(output_dir, f"{accession}_2.fastq")
        unpaired = os.path.join(output_dir, f"{accession}.fastq")

        if read_type == sra.PAIRED_END_SHORT_READS:
            if not (os.path.exists(paired) and os.path.exists(paired_mate)):
                raise ConfigError(f"NCBI describes the SRA run {accession} as paired-end, so anvi'o expected "
                                  f"`fasterq-dump` to produce a pair of FASTQ files for it, and that is not what "
                                  f"happened. Sequencing metadata at NCBI is filled in by whoever submitted the "
                                  f"data, so every now and then it disagrees with the data itself. You can tell "
                                  f"anvi'o what this run really holds by editing its `read_type` in your SRA "
                                  f"metadata file. The log file at '{log_path}' shows what `fasterq-dump` had to "
                                  f"say for itself.")

            # A paired run can also yield a file of reads whose mates did not survive filtering
            # upstream. Anvi'o works with the pairs alone, so those are dropped here rather than
            # quietly finding their way into the analysis.
            if os.path.exists(unpaired):
                with open(log_path, 'a') as log_file:
                    log_file.write(f"Discarding unpaired reads in {unpaired}; anvi'o only uses the pairs.\n")
                os.remove(unpaired)

            to_compress = [paired, paired_mate]
        else:
            if not os.path.exists(unpaired):
                raise ConfigError(f"Anvi'o expected `fasterq-dump` to produce a single FASTQ file for the SRA run "
                                  f"{accession}, and that is not what happened. If NCBI's description of this run "
                                  f"is wrong, you can correct its `read_type` in your SRA metadata file. The log "
                                  f"file at '{log_path}' shows what `fasterq-dump` had to say for itself.")

            to_compress = [unpaired]

        shell(f"pigz --processes {threads} --verbose {' '.join(to_compress)} >> {log_path} 2>&1")

        if os.path.exists(scratch_dir):
            shutil.rmtree(scratch_dir)


    def get_sra_target_files(self):
        """The release flags, which is what pulls the download machinery into the workflow."""

        if not self.sra_mode:
            return []

        return [self.get_release_flag_path(u) for u in self.download_units]
