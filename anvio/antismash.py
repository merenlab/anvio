# pylint: disable=line-too-long
"""Runs antiSMASH on the contigs in a contigs-db and imports the results.

The design mirrors other anvi'o annotation drivers such as `anvio/pfam.py`: the `Antismash`
class does its configuration and sanity checks in `__init__`, the actual work happens
in `process()`. The companion command-line client lives in `anvio/cli/run_antismash.py`.

What we store, and where:
  * Each biosynthetically relevant gene is annotated across four gene functions sources, so every facet
    behaves like any other functional annotation source: 'antiSMASH' (smCOG), 'antiSMASH_ROLE' (role
    within the cluster), 'antiSMASH_DOMAIN' (within gene biosynthetic domains), and 'antiSMASH_REGION'
    (which BGC region it belongs to; accession = an assembly-unique region id, function = the BGCs final
    product).
  * A flat, per-region summary is written to a stand-alone TAB-delimited file, which maps each compact
    region id back to its full contig name, coordinates, and product.
"""

import os
import re
import sys
import glob
import json
import shlex
import shutil
import hashlib
import subprocess

from collections import defaultdict

import anvio
import anvio.tables as t
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.dbops import ContigsDatabase, ContigsSuperclass
from anvio.errors import ConfigError, CommandError
from anvio.genomedescriptions import GenomeDescriptions
from anvio.tables.genefunctions import TableForGeneFunctions


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['jessika-fuessel']


run = terminal.Run()
progress = terminal.Progress()

# the functional annotation source under which per-CDS BGC roles are stored
ANTISMASH_SOURCE = "antiSMASH"

# the functional annotation source under which each gene's role in the cluster is stored
ANTISMASH_ROLE_SOURCE = "antiSMASH_ROLE"

# the functional annotation source under which each gene's biosynthetic domains are stored
ANTISMASH_DOMAIN_SOURCE = "antiSMASH_DOMAIN"

# the functional annotation source under which per-gene BGC-region membership is stored
ANTISMASH_REGION_SOURCE = "antiSMASH_REGION"

# every source anvi-run-antismash writes. They are always added, dropped, and checked for together
ANTISMASH_SOURCES = [ANTISMASH_SOURCE, ANTISMASH_ROLE_SOURCE, ANTISMASH_DOMAIN_SOURCE, ANTISMASH_REGION_SOURCE]

# the name of the file in the report directory that translates the record ids in antiSMASH's own output files
# (HTML, GenBank, JSON) back to anvi'o contig names (see `Antismash._get_record_ids`)
CONTIG_NAME_MAPPING_FILE_NAME = "contig-name-mapping.txt"

# antiSMASH annotation analyses that always run: RRE-Finder (--rre, for RiPP-recognition elements) and the
# Active Site Finder (--asf). Both add BGC-specific information with no anvi'o equivalent. We deliberately do
# NOT add --tigrfam / --pfam2go: antiSMASH's core detection already contributes the relevant TIGRFAM/Pfam
# domains (verified: the imported annotations are identical with or without those two flags), and --pfam2go
# only adds Gene Ontology terms to the HTML report. A user who wants them can pass them via --antismash-add-on.
ANTISMASH_ANNOTATION_ARGS = ["--asf", "--rre"]

# heavy comparison/tree analyses that run only with --include-detailed-output. They are slow and produce
# large outputs (ClusterBlast, MIBiG, smCOG trees) that populate the report's extra tabs, but they do NOT
# change the per-gene annotations imported into the database.
ANTISMASH_DETAILED_ARGS = ["--cb-general", "--cb-knownclusters", "--cb-subclusters", "--cc-mibig", "--smcog-trees"]

# flags `_run_antismash` always passes to antiSMASH on its own, which therefore cannot also come through
# --antismash-add-on
ANTISMASH_FLAGS_SET_BY_ANVIO = ["--cpus", "-c", "--taxon", "--genefinding-tool", "--genefinding-gff3", "--output-dir",
                                "--logfile", "--html-title", "--databases"]

# per-(meta)genome result files antiSMASH always writes, but which we keep only with --include-detailed-output
ANTISMASH_DETAILED_FILE_SUFFIXES = (".gbk", ".json", ".zip")


def _bgc_function_name(smcog_descriptions, domain_names):
    """The smCOG annotation if available, else its domain(s)."""

    return "; ".join(smcog_descriptions) or "; ".join(domain_names)


def get_target_genomes(args):
    """Resolve the target databases from -c (one) or -e (many) into an ordered {name: contigs_db_path} map.

    The external-genomes reader validates each database and, unlike a raw file read, resolves the paths
    relative to the external-genomes file.
    """

    contigs_db_path = args.__dict__.get('contigs_db')
    external_genomes_path = args.__dict__.get('external_genomes')

    if not contigs_db_path and not external_genomes_path:
        raise ConfigError("You didn't tell the script which contig-dbs to work with using either `-c` contig-db or "
                          "`-e` external-genomes file.")
    if contigs_db_path and external_genomes_path:
        raise ConfigError("You need to decide, either we are looking at one assembly or at many, provide either an "
                          "individual contig-db or an external-genomes file, not both.")

    if contigs_db_path:
        utils.is_contigs_db(contigs_db_path)
        return {os.path.splitext(os.path.basename(contigs_db_path))[0]: contigs_db_path}

    genome_descriptions = GenomeDescriptions(args, run=terminal.Run(verbose=False), progress=terminal.Progress(verbose=False))
    genome_descriptions.load_genomes_descriptions(skip_functions=True, init=False)
    return {name: genome_descriptions.genomes[name]['contigs_db_path'] for name in genome_descriptions.genomes}


class Antismash(object):
    """Run antiSMASH on a contigs-db and import its results."""

    def __init__(self, args, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.num_threads = A('num_threads') or 1
        self.taxon = A('taxon') or 'bacteria'
        self.antismash_data_dir = A('antismash_data_dir')
        self.output_dir = A('output_dir')
        self.regions_output_path = A('regions_output')
        self.min_contig_length = A('min_contig_length') or 3000
        self.just_do_it = A('just_do_it') or False
        self.include_detailed_output = A('include_detailed_output') or False
        self.antismash_add_on = A('antismash_add_on') or ''

        # antiSMASH usually lives in its own conda environment (its dependency pins conflict with anvi'o's),
        # so its executable is typically not on anvi'o's PATH. Here we just record the optional installation
        # directory the user gave us; `sanity_check` resolves the actual executable (see `_resolve_antismash_bin`).
        self.antismash_dir = A('antismash_dir')

        # a caller may point us at a directory to keep intermediate files (handy for debugging); otherwise
        # each genome gets its own temporary one, created and removed in `_annotate`.
        self.work_dir_arg = A('work_dir')

        # the target databases: one (with -c) or many (with -e). Resolved once, up front.
        self.genomes = get_target_genomes(args)

        # per-genome state below is (re)set in `_annotate` as we walk the target databases.
        self.contigs_db_path = None
        self.work_dir = None
        self.log_file_path = None
        self.antismash_version = None
        self.regions_by_contig = defaultdict(list)
        self.gene_annotations = []
        self.region_records = []

        self.sanity_check()


    def sanity_check(self):
        # locate the antiSMASH executable (see _resolve_antismash_bin for the search order). This is independent
        # of the target databases, so we do it once here rather than per genome.
        self.antismash_bin = self._resolve_antismash_bin()

        # make sure the add-on flags parse, and that they do not repeat anything anvi'o already sets
        self.antismash_add_on_args = self._get_antismash_add_on_args()

        # the single-database output options make no sense for a batch (each database would clobber the next),
        # so refuse that combination early, with a clear message.
        if len(self.genomes) > 1:
            for flag, value in [("--output-dir", self.output_dir), ("--regions-output", self.regions_output_path),
                                ("--work-dir", self.work_dir_arg)]:
                if value:
                    raise ConfigError(f"Looks like you gave a specific location to store the results in using "
                                      f"`{flag}` and you want to run the script on a set of contig-dbs specified in "
                                      f"your external-genomes file. Each contig-db will need its own output directory "
                                      f"and anvi'o will name it for you using your db's name, just remove the `{flag}` "
                                      f"and you are golden.")

        # antiSMASH refuses to write into an output directory that already has things in it, and leftovers
        # from an earlier run would get mixed with the new ones, so a user-defined work directory must be new
        # (unless the user asks for the existing one to be removed with --force-overwrite).
        if self.work_dir_arg:
            filesnpaths.check_output_directory(self.work_dir_arg)

        # figure out up front which databases are already annotated (so a batch can fail early rather than
        # half-way through, if something is off), and make sure we will not overwrite anything we shouldn't.
        self.genomes_to_skip = set()
        for genome_name, contigs_db_path in self.genomes.items():
            if self._already_annotated(contigs_db_path) and not self.just_do_it:
                # with a single database this is an error; in a batch we skip databases that are already
                # annotated, so a re-run resumes cleanly.
                if len(self.genomes) == 1:
                    raise ConfigError(f"You have already annotated this db with antiSMASH. You want to re-run it? "
                                      f"Just add `--just-do-it` to your command and anvi'o will replace the previous "
                                      f"annotations with new ones (if the report directory from the previous run is "
                                      f"still around, you will also need `--force-overwrite` to replace it, or you "
                                      f"can use `--output-dir` to put the new report elsewhere). If you want to remove previous "
                                      f"annotations yourself, you can run `anvi-delete-functions -c {contigs_db_path} "
                                      f"--annotation-sources {','.join(ANTISMASH_SOURCES)}`. If you can't remember "
                                      f"what you've annotated, you can double-check with `anvi-db-info {contigs_db_path}`.")

                self.genomes_to_skip.add(genome_name)
                continue

            # the report directory must be new, unless the user asks for the existing one to be removed with
            # --force-overwrite (in which case it is removed right here, before anything else happens).
            filesnpaths.check_output_directory(self._report_dir(contigs_db_path))


    def _get_antismash_add_on_args(self):
        """Split the --antismash-add-on string into a list of arguments, refusing flags anvi'o already sets."""

        if not self.antismash_add_on:
            return []

        try:
            add_on_args = shlex.split(self.antismash_add_on)
        except ValueError as e:
            raise ConfigError(f"anvi'o could not make sense of the flags you passed with `--antismash-add-on` "
                              f"('{self.antismash_add_on}'): {e}. Perhaps there is a quote that is not closed?")

        # everything we always pass to antiSMASH ourselves, plus the heavy analyses when they are already on
        flags_set_by_anvio = set(ANTISMASH_FLAGS_SET_BY_ANVIO + ANTISMASH_ANNOTATION_ARGS)
        if self.include_detailed_output:
            flags_set_by_anvio.update(ANTISMASH_DETAILED_ARGS)

        # a flag can also come in the '--flag=value' form
        repeated_flags = sorted(set(arg.split('=', 1)[0] for arg in add_on_args).intersection(flags_set_by_anvio))
        if repeated_flags:
            raise ConfigError(f"Some of the flags you passed with `--antismash-add-on` are flags anvi'o already sets "
                              f"for antiSMASH on its own: {', '.join(repeated_flags)}. Passing them twice would at "
                              f"best be redundant and at worst confuse antiSMASH. If you want to change something "
                              f"anvi'o controls, please use the corresponding anvi-run-antismash flag instead (e.g., "
                              f"`--num-threads` instead of `--cpus`, `--antismash-data-dir` instead of `--databases`, "
                              f"or `--include-detailed-output` for the ClusterBlast/MIBiG analyses).")

        return add_on_args


    def _already_annotated(self, contigs_db_path):
        """Whether a database already carries (any of the) antiSMASH annotations from a previous run."""

        contigs_db = ContigsDatabase(contigs_db_path, run=terminal.Run(verbose=False), progress=terminal.Progress(verbose=False))
        existing_sources = set(contigs_db.meta['gene_function_sources'] or [])
        contigs_db.disconnect()
        return bool(existing_sources.intersection(ANTISMASH_SOURCES))


    def _resolve_antismash_bin(self):
        """Locate the antiSMASH executable and return its absolute path.

        antiSMASH is installed in its own conda environment (its dependencies conflict with
        anvi'o's). We look for it in this order:
          1. the installation directory given via --antismash-dir or the ANTISMASH_DIR environment variable,
          2. an 'antismash' on the PATH (e.g., if it happens to be in the active environment),
          3. an automatic search of the sibling conda environments.
        """

        # 1. an explicit installation directory (the antiSMASH conda environment)
        antismash_dir = self.antismash_dir or os.environ.get('ANTISMASH_DIR')
        if antismash_dir:
            antismash_dir = os.path.expanduser(antismash_dir)
            for candidate in [os.path.join(antismash_dir, "bin", "antismash"), os.path.join(antismash_dir, "antismash")]:
                if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
                    return candidate
            raise ConfigError(f"The antiSMASH program isn't where anvi'o was expecting it :(. anvi'o looked in the "
                              f"directory you provided ('{antismash_dir}'), but there was no 'antismash' executable "
                              f"inside it. Please double check the path you provided with `--antismash-dir` to point "
                              f"anvi'o to the location of antiSMASH.")

        # 2. antismash on the PATH
        on_path = utils.is_program_exists('antismash', dont_raise=True)
        if on_path:
            return on_path

        # 3. auto-discover across the sibling conda environments
        candidates = self._find_antismash_in_conda_envs()
        if len(candidates) == 1:
            self.run.info("antiSMASH found in a conda environment", candidates[0])
            return candidates[0]
        if len(candidates) > 1:
            raise ConfigError(f"Seems like you installed antiSMASH more than once in different environments and now "
                              f"anvi'o is confused. It found these: {', '.join(candidates)}. You can choose one of the "
                              f"installations using the `--antismash-dir` flag.")

        raise ConfigError("anvi'o looked for antiSMASH everywhere, but it doesn't seem to be anywhere :(. antiSMASH is "
                          "not part of the anvi'o installation because both rely on some of the same dependencies but "
                          "different versions of them. Therefore, antiSMASH needs to live in its own environment and be "
                          "installed separately, just follow the installation instructions on the `anvi-run-antismash` "
                          "program webpage. Maybe you also did an excellent job at hiding antiSMASH's location, in that "
                          "case you can let anvi'o know where to find it using the `--antismash-dir` flag.")


    def _find_antismash_in_conda_envs(self):
        """Return any antiSMASH executables found across the conda environments as a sorted list of paths."""

        # the active environment's prefix, from which we can reach the directory that holds all environments
        conda_prefix = os.environ.get('CONDA_PREFIX') or os.path.dirname(os.path.dirname(sys.executable))

        # cover both "we are in a named env, look at sibling envs" and "we are in the base env, look under envs/"
        search_roots = {os.path.dirname(conda_prefix), os.path.join(conda_prefix, "envs")}

        candidates = set()
        for root in search_roots:
            candidates.update(glob.glob(os.path.join(root, "*", "bin", "antismash")))

        return sorted(candidate for candidate in candidates if os.access(candidate, os.X_OK))


    def process(self):
        """Annotate one database (with -c) or every database in an external-genomes file (with -e)."""

        total = len(self.genomes)
        annotated, skipped = [], []
        for i, (genome_name, contigs_db_path) in enumerate(self.genomes.items()):
            if total > 1:
                self.run.warning(None, header=f"[{i + 1}/{total}] {genome_name}", lc="green", nl_before=1)

            if genome_name in self.genomes_to_skip:
                self.run.warning(f"'{genome_name}' already contains antiSMASH annotations, so anvi'o is skipping it. "
                                 f"Use --just-do-it if you want to re-annotate it.")
                skipped.append(genome_name)
                continue

            self._annotate(contigs_db_path)
            annotated.append(genome_name)

        if total > 1:
            message = f"Done. Annotated {len(annotated)} of {total} databases with antiSMASH."
            if skipped:
                message += f" Skipped {len(skipped)} that were already annotated: {', '.join(skipped)}."
            self.run.info_single(message, nl_before=1, nl_after=1, mc="green")


    def _annotate(self, contigs_db_path):
        """Run antiSMASH on a single database and import its results."""

        self.contigs_db_path = contigs_db_path
        self.args.contigs_db = contigs_db_path   # so the ContigsSuperclass below reads this database

        # reset the per-genome result state so nothing leaks between databases in a batch run.
        self.antismash_version = None
        self.regions_by_contig = defaultdict(list)
        self.gene_annotations = []
        self.region_records = []
        self.record_id_of_contig = {}
        self.contig_of_record_id = {}

        # all intermediate files live in a working directory. if the user gave us one (single-database runs
        # only), use it; otherwise make a temporary one and remove it at the end (unless we are in --debug
        # mode). The antiSMASH log is copied out of here into the report directory first, so it always survives.
        made_temporary_work_dir = self.work_dir_arg is None
        self.work_dir = self.work_dir_arg or filesnpaths.get_temp_directory_path()
        if not made_temporary_work_dir:
            filesnpaths.gen_output_directory(self.work_dir)

        # antiSMASH writes its log here; we copy it next to the output at the end, on success or failure.
        self.log_file_path = os.path.join(self.work_dir, "antismash.log")

        # a single ContigsSuperclass is built once and reused throughout, since it is expensive to set up.
        contigs_super = ContigsSuperclass(self.args, r=terminal.Run(verbose=False), p=terminal.Progress(verbose=False))

        try:
            fasta, gff = self._export_inputs(contigs_super)
            antismash_output_dir = self._run_antismash(fasta, gff)
            self._parse_json(antismash_output_dir, fasta)

            # enrich each region with the genes it contains; both the database annotations and the
            # per-region output file are derived from these records, so the membership is computed once.
            self.region_records = self._build_region_records(contigs_super)
            self._store_results()

            # save the antiSMASH report (done before writing the per-region file, which lives inside the
            # report directory). By default we keep a lean, view-only report; with --include-detailed-output
            # we keep the complete report, including the .gbk / JSON files and the comparison data. The report
            # directory does not exist at this point, since `sanity_check` made sure of it.
            report_dir = self._report_dir()
            if os.path.exists(report_dir):
                shutil.rmtree(report_dir)
            shutil.copytree(antismash_output_dir, report_dir)
            if not self.include_detailed_output:
                self._slim_report(report_dir)

            if self.gene_annotations:
                self._write_genes_txt(self._genes_txt_path())
            if self.region_records:
                self._write_regions_txt(self._regions_txt_path())
            self._write_contig_name_mapping(self._contig_name_mapping_path())

            self._report()
        finally:
            # copy the antiSMASH log next to its output (on success or failure) so it is easy to find, then
            # remove the temporary working directory unless we are in --debug mode.
            self._save_antismash_log()
            if made_temporary_work_dir and not anvio.DEBUG:
                shutil.rmtree(self.work_dir, ignore_errors=True)
            elif made_temporary_work_dir:
                self.run.warning(f"Just a heads up, you added the `--debug` flag and therefore anvi'o is not deleting "
                                 f"the temporary files created by anvi-run-antismash. You can find the directory with "
                                 f"the usually temporary files here: {self.work_dir}")


    def _export_inputs(self, contigs_super):
        """Dump a nucleotide FASTA and a GFF3 from the contigs-db, ready for antiSMASH."""

        fasta_raw = os.path.join(self.work_dir, "contigs.raw.fa")
        fasta = os.path.join(self.work_dir, "contigs.fa")
        gff = os.path.join(self.work_dir, "contigs.gff")
        gff_clean = os.path.join(self.work_dir, "contigs.clean.gff")

        self.progress.new("Exporting from the contigs database")

        # antiSMASH cannot meaningfully call clusters on short fragments, so we only export the contigs that
        # are long enough, and below we only keep their gene calls.
        self.progress.update(f"finding contigs >= {self.min_contig_length} bp")
        surviving_contigs = set([contig_name for contig_name, contig_info in contigs_super.contigs_basic_info.items()
                                 if contig_info['length'] >= self.min_contig_length])

        # we must stop here if nothing survives, since asking the export function for an empty list of
        # sequences would give us ALL of them.
        if not surviving_contigs:
            self.progress.end()
            raise ConfigError(f"None of the contigs in '{self.contigs_db_path}' is at least {self.min_contig_length} bp "
                              f"long, so there is nothing for antiSMASH to look at :/ If you really want to search "
                              f"shorter contigs for BGCs, you can lower the cutoff with `--min-contig-length`, but "
                              f"keep in mind that antiSMASH is not very good at finding clusters on short fragments.")

        # antiSMASH quietly renames records with long ids, after which we could no longer match its regions to
        # our contigs. So we give every contig an id antiSMASH will keep, and translate back when parsing.
        # The translation is also written into the report directory (see `_write_contig_name_mapping`).
        self.record_id_of_contig = self._get_record_ids(surviving_contigs)
        self.contig_of_record_id = {record_id: contig_name for contig_name, record_id in self.record_id_of_contig.items()}

        self.progress.update("nucleotide FASTA")
        utils.export_sequences_from_contigs_db(self.contigs_db_path, fasta_raw, seq_names_to_export=sorted(surviving_contigs),
                                               run=terminal.Run(verbose=False))
        with open(fasta_raw) as raw_fh, open(fasta, "w") as fasta_fh:
            for line in raw_fh:
                if line.startswith(">"):
                    line = f">{self.record_id_of_contig[line[1:].strip()]}\n"
                fasta_fh.write(line)

        # we ask for simple headers because we only need each gene's ID and coordinates here (we rewrite
        # the GFF ourselves below); the richer `;Name=...` header carries information antiSMASH ignores.
        self.progress.update("gene calls as GFF3")
        contigs_super.gen_GFF3_file_of_sequences_for_gene_caller_ids(output_file_path=gff, simple_headers=True)

        self.progress.update("cleaning the GFF3")
        # anvi'o exports a flat list of CDS rows with empty source/phase columns, but antiSMASH wants
        # (1) a gene -> CDS hierarchy via a Parent attribute, (2) a non-empty source column, and (3) a
        # valid phase (0/1/2) for every CDS. We synthesise the gene parent on the fly and set phase to 0,
        # which is correct for the single-exon prokaryotic genes in an anvi'o contigs database.
        with open(gff) as source_fh, open(gff_clean, "w") as clean_fh:
            clean_fh.write("##gff-version 3\n")
            for line in source_fh:
                if line.startswith("##gff-version"):
                    continue
                if line.startswith("#"):
                    clean_fh.write(line)
                    continue

                columns = line.rstrip("\n").split("\t")
                if len(columns) < 9 or columns[0] not in surviving_contigs:
                    continue
                if columns[2] != "CDS":
                    clean_fh.write(line)
                    continue

                attributes = dict(kv.split("=", 1) for kv in columns[8].split(";") if "=" in kv)
                anvio_id = attributes.get("ID", "unknown")

                # emit a 'gene' feature for antiSMASH to anchor on, then its 'CDS' with phase and Parent.
                clean_fh.write("\t".join([columns[0], "anvio", "gene", columns[3], columns[4], ".", columns[6], ".",
                                          f"ID=gene_{anvio_id}"]) + "\n")
                clean_fh.write("\t".join([columns[0], "anvio", "CDS", columns[3], columns[4], ".", columns[6], "0",
                                          f"ID=cds_{anvio_id};Parent=gene_{anvio_id};locus_tag={anvio_id}"]) + "\n")

        self.progress.end()
        return fasta, gff_clean


    @staticmethod
    def _get_record_ids(contig_names):
        """Return a {contig name: record id} map with ids antiSMASH will not rename.

        antiSMASH renames records with ids longer than 16 characters (and its rules for what it alters are
        its own business, which may change between versions). Rather than trying to predict which names it
        would keep, we give every contig a short, predictable id, so the naming in antiSMASH's output files
        is consistent and independent of antiSMASH's rules. These ids are at most 16 characters for up to a
        billion contigs.
        """

        return {contig_name: f"anvio_c{i:06d}" for i, contig_name in enumerate(sorted(contig_names))}


    def _run_antismash(self, fasta, gff):
        """Run antiSMASH and return the path to its output directory."""

        output_dir = os.path.join(self.work_dir, "antismash_out")
        log_file_path = self.log_file_path

        cmd = [self.antismash_bin,
               "--cpus", str(self.num_threads),
               "--taxon", self.taxon,
               "--genefinding-tool", "none",
               "--genefinding-gff3", gff,
               *ANTISMASH_ANNOTATION_ARGS,
               "--output-dir", output_dir,
               "--logfile", log_file_path,
               "--html-title", os.path.basename(self.contigs_db_path)]
        # the slow, output-heavy comparison analyses only run when the user asks for detailed output.
        if self.include_detailed_output:
            cmd += ANTISMASH_DETAILED_ARGS
        if self.antismash_data_dir:
            cmd += ["--databases", self.antismash_data_dir]
        cmd += self.antismash_add_on_args
        cmd.append(fasta)

        # put antiSMASH's own environment at the front of PATH so it finds its helper binaries
        # (hmmpfam2, etc.) rather than the ones in the anvi'o environment.
        antismash_env = os.environ.copy()
        antismash_env["PATH"] = os.path.dirname(self.antismash_bin) + os.pathsep + antismash_env.get("PATH", "")

        self.progress.new("Running antiSMASH")
        self.progress.update(f"this can take a while (cpus={self.num_threads})")
        with open(log_file_path, "a") as log_fh:
            return_code = subprocess.call(cmd, env=antismash_env, stdout=log_fh, stderr=subprocess.STDOUT)
        self.progress.end()

        if return_code != 0:
            raise CommandError(f"This didn't work, antiSMASH could not annotate the BGCs in your database :(. One "
                               f"reason could be that you didn't download the antiSMASH databases yet? If that's the "
                               f"case, go to your antiSMASH environment and run `download-antismash-databases` and then "
                               f"come back to your anvi'o environment and try again. If that's not it, we also don't "
                               f"know. Looking at the antismash.log (anvi'o saved it next to your output at "
                               f"{os.path.join(self._report_dir(), 'antismash.log')}) should help you to figure out "
                               f"what went wrong.")

        return output_dir


    def _parse_json(self, antismash_output_dir):
        """Walk the antiSMASH JSON output and collect per-CDS annotations and per-region metadata."""

        json_files = [f for f in os.listdir(antismash_output_dir) if f.endswith(".json")]
        if not json_files:
            raise ConfigError(f"antiSMASH finished without complaining, but then didn't generate the JSON output file "
                              f"anvi'o reads its results from :/ Nothing was added to your contigs database. The "
                              f"antismash.log file (anvi'o saved it at {os.path.join(self._report_dir(), 'antismash.log')}) "
                              f"is the best place to look for what went wrong.")

        with open(os.path.join(antismash_output_dir, json_files[0])) as json_fh:
            payload = json.load(json_fh)

        self.antismash_version = payload.get("version", "unknown")

        for record in payload.get("records", []):
            # translate antiSMASH's record id back to our contig name (see `_get_record_ids`)
            contig = self.contig_of_record_id.get(record.get("id"))
            if contig is None:
                raise ConfigError(f"antiSMASH reported results for a sequence called '{record.get('id')}', but that is "
                                  f"not one of the sequences anvi'o gave it. Most likely antiSMASH renamed it, which "
                                  f"anvi'o tries hard to prevent, so this is not supposed to happen :/ The antismash.log "
                                  f"file may tell you more, and it would be great if you could let the anvi'o "
                                  f"developers know about it.")
            for feature in record.get("features", []):
                feature_type = feature.get("type")
                qualifiers = feature.get("qualifiers", {})

                if feature_type == "region":
                    self.regions_by_contig[contig].append({
                        "region_number": int(qualifiers.get("region_number", ["0"])[0]),
                        "product": ";".join(qualifiers.get("product", [])),
                        "contig_edge": qualifiers.get("contig_edge", ["False"])[0] == "True",
                        "location": feature.get("location", ""),
                    })
                elif feature_type == "CDS":
                    gene_callers_id = self._gene_callers_id_from_locus_tag(qualifiers.get("locus_tag", [None])[0])
                    if gene_callers_id is None:
                        continue

                    # pull the useful, human-meaningful pieces out of antiSMASH's verbose qualifiers:
                    #   - the gene's role in the cluster (biosynthetic / regulatory / transport / ...),
                    #   - the biosynthetic domains it carries (as clean (name, e-value) pairs), and
                    #   - its smCOG classification(s) (id, readable description, e-value), which give the
                    #     clearest functional name when present.
                    role = qualifiers.get("gene_kind", [""])[0]
                    domains = [self._parse_domain(d) for d in qualifiers.get("sec_met_domain", [])]
                    smcogs = [self._parse_smcog(gf) for gf in qualifiers.get("gene_functions", []) if "(smcogs)" in gf]
                    if not (role or domains or smcogs):
                        continue  # a CDS, but with nothing biosynthetically relevant to say about it

                    self.gene_annotations.append({
                        "gene_callers_id": gene_callers_id,
                        "contig": contig,
                        "role": role,
                        "domains": domains,   # [(name, e_value), ...]
                        "smcogs": smcogs,     # [(smcog_id, description, e_value), ...]
                    })


    def _build_region_records(self, contigs_super):
        """Enrich each parsed BGC region with the gene caller ids it contains.

        Returns a list of records, one per region, each a dict with keys: region_id (a compact,
        assembly-unique id, see `_region_id`), contig, product, start, end, contig_edge, and
        gene_callers_ids (the sorted ids of genes fully contained within the region). antiSMASH regions
        do not overlap, so each gene belongs to at most one region. Both the database annotations and the
        per-region output file are derived from these records.
        """

        # group genes by contig once, so each region only scans the genes on its own contig.
        genes_by_contig = defaultdict(list)
        for gene_callers_id, gene_call in contigs_super.genes_in_contigs_dict.items():
            genes_by_contig[gene_call["contig"]].append((gene_callers_id, gene_call["start"], gene_call["stop"]))

        records = []
        for contig, regions in self.regions_by_contig.items():
            for region in regions:
                start, end = self._parse_location(region["location"])
                gene_callers_ids = sorted(gene_callers_id for gene_callers_id, gene_start, gene_stop in genes_by_contig[contig]
                                          if gene_start >= start and gene_stop <= end)
                records.append({
                    "region_id": self._region_id(contig, region["region_number"]),
                    "contig": contig,
                    "product": region["product"] or "unknown",
                    "start": start,
                    "end": end,
                    "contig_edge": region["contig_edge"],
                    "gene_callers_ids": gene_callers_ids,
                })
        return records


    def _store_results(self):
        """Import the parsed antiSMASH results into the contigs database as three function sources.

        The gene functions table is limited to (source, accession, function, e_value) columns, so the
        different facets of an antiSMASH annotation are stored as separate sources:
          * 'antiSMASH': smCOG description, else biosynthetic domain,
          * 'antiSMASH_ROLE': the gene's role in the cluster (biosynthetic / regulatory / transport / ...),
          * 'antiSMASH_DOMAIN': each biosynthetic domain the gene carries (the full domain inventory),
          * 'antiSMASH_REGION': which BGC region the gene belongs to (accession = region id, function = product).
        A gene with several annotations gets several rows (never comma-joined into one cell). Keeping every
        facet in the database is what lets the cross-genome matrix program read them back without the files.
        """

        function_rows, role_rows, domain_rows = [], [], []
        for gene in self.gene_annotations:
            gene_callers_id = gene["gene_callers_id"]

            # the functional name: prefer the smCOG classification(s), else fall back to the domain(s).
            if gene["smcogs"]:
                for smcog_id, description, e_value in gene["smcogs"]:
                    function_rows.append({"gene_callers_id": gene_callers_id, "source": ANTISMASH_SOURCE,
                                          "accession": smcog_id, "function": description, "e_value": e_value})
            elif gene["domains"]:
                for name, e_value in gene["domains"]:
                    function_rows.append({"gene_callers_id": gene_callers_id, "source": ANTISMASH_SOURCE,
                                          "accession": name, "function": name, "e_value": e_value})

            # the biosynthetic domains get their own source, so the full domain inventory is in the database
            # (and thus available to the cross-genome matrix), regardless of whether a smCOG name was found.
            for name, e_value in gene["domains"]:
                domain_rows.append({"gene_callers_id": gene_callers_id, "source": ANTISMASH_DOMAIN_SOURCE,
                                    "accession": name, "function": name, "e_value": e_value})

            if gene["role"]:
                role_rows.append({"gene_callers_id": gene_callers_id, "source": ANTISMASH_ROLE_SOURCE,
                                  "accession": gene["role"], "function": gene["role"], "e_value": 0.0})

        region_rows = [{"gene_callers_id": gene_callers_id, "source": ANTISMASH_REGION_SOURCE,
                        "accession": record["region_id"], "function": record["product"], "e_value": 0.0}
                       for record in self.region_records for gene_callers_id in record["gene_callers_ids"]]

        all_rows = function_rows + role_rows + domain_rows + region_rows

        gene_functions_table = TableForGeneFunctions(self.contigs_db_path,
                                                     run=terminal.Run(verbose=False),
                                                     progress=terminal.Progress(verbose=False))

        # remove every antiSMASH source a previous run left behind (this only happens with --just-do-it). Simply
        # adding the new rows would only replace the sources that show up in them again, so a source with no rows
        # this time (or all of them, if antiSMASH found nothing) would keep its stale annotations.
        contigs_db = ContigsDatabase(self.contigs_db_path, run=terminal.Run(verbose=False), progress=terminal.Progress(verbose=False))
        previous_sources = [source for source in ANTISMASH_SOURCES if source in set(contigs_db.meta['gene_function_sources'] or [])]
        if previous_sources:
            self.run.info("Previous antiSMASH annotations removed", ", ".join(previous_sources))
            gene_functions_table.drop_functions(database=contigs_db.db, sources_to_drop=previous_sources)
        contigs_db.disconnect()

        if all_rows:
            self.progress.new("Storing annotations")
            self.progress.update(f"{len(all_rows):,} rows")
            # `create` wants a dict keyed by a unique entry id. We do NOT pass
            # drop_previous_annotations_first=True: that would delete EVERY functional annotation source in
            # the database (KOfam, COGs, ...), not just ours.
            functions_dict = {i: row for i, row in enumerate(all_rows)}
            try:
                gene_functions_table.create(functions_dict)
            finally:
                self.progress.end()
        else:
            self.run.warning(f"Just to let you know, antiSMASH did not find any BGCs in '{self.contigs_db_path}'. We are "
                             f"not sure if that's what you were expecting, but anvi'o will still register the antiSMASH "
                             f"annotation sources even though there was nothing to annotate, so you will know later that "
                             f"you did look for BGCs in this database.")

        # we always register all antiSMASH sources, including the ones that ended up with no rows, so a database
        # that was searched is always recognizable as such (whatever antiSMASH found), and distinguishable from
        # one that was never processed.
        empty_sources = set(ANTISMASH_SOURCES) - set(row["source"] for row in all_rows)
        if empty_sources:
            gene_functions_table.add_empty_sources_to_functional_sources(empty_sources)

    def _report_dir(self):
        """Where the antiSMASH report is written: the user's --output-dir, or a default next to the db."""

        if self.output_dir:
            return self.output_dir
        return os.path.splitext(self.contigs_db_path)[0] + "-ANTISMASH"


    def _slim_report(self, report_dir):
        """Drop the large per-(meta)genome result files from a report, keeping it a compact, view-only one.

        The .gbk / JSON / zip files are only kept with --include-detailed-output (the heavy ClusterBlast and
        smCOG directories are not generated at all in the default mode). Removing these here does not touch
        the annotations, which live in the contigs database regardless.
        """

        for entry in os.listdir(report_dir):
            path = os.path.join(report_dir, entry)
            if os.path.isfile(path) and entry.endswith(ANTISMASH_DETAILED_FILE_SUFFIXES):
                os.remove(path)


    def _save_antismash_log(self):
        """Copy the antiSMASH log into the report directory, so it sits next to the output it describes.

        Called on both success and failure. On failure the report directory may not exist yet, so we
        create it; if antiSMASH never ran far enough to write a log, there is nothing to copy.
        """

        if not os.path.exists(self.log_file_path):
            return

        report_dir = self._report_dir()
        if not os.path.exists(report_dir):
            os.makedirs(report_dir)
        shutil.copy2(self.log_file_path, os.path.join(report_dir, "antismash.log"))


    def _regions_txt_path(self):
        """Resolve where the per-region TAB-delimited file should be written (inside the report directory)."""

        return self.regions_output_path or os.path.join(self._report_dir(), "antismash_regions.txt")


    def _write_regions_txt(self, txt_path):
        """Write a comprehensive, self-contained per-region summary as a flat TAB-delimited file.

        Each row fully describes one BGC region (its compact id, full contig name, product/type,
        coordinates, and the genes it contains), so the file can be used on its own, without going back
        to the database. The same per-gene membership is also queryable in the 'antiSMASH_REGION' gene
        functions source, linked by region_id.
        """

        region_rows = {}
        for record in self.region_records:
            # the gene ids are joined with ';' rather than ',' on purpose: a comma-separated list of
            # numbers is misread by spreadsheet software (e.g. Excel) as one very large number and shown
            # with trailing zeros, whereas a semicolon-separated list stays intact as text.
            gene_callers_ids = ";".join(str(gene_callers_id) for gene_callers_id in record["gene_callers_ids"])
            region_rows[record["region_id"]] = {
                "contig": record["contig"],
                "product": record["product"],
                "start": record["start"],
                "end": record["end"],
                "length": record["end"] - record["start"],
                "num_genes": len(record["gene_callers_ids"]),
                "gene_callers_ids": gene_callers_ids,
                "contig_edge": record["contig_edge"],
            }

        utils.store_dict_as_TAB_delimited_file(region_rows, txt_path, key_header="region_id",
                                               headers=["region_id", "contig", "product", "start", "end", "length",
                                                        "num_genes", "gene_callers_ids", "contig_edge"])


    def _contig_name_mapping_path(self):
        """Where the record id -> contig name translation is written (inside the report directory)."""

        return os.path.join(self._report_dir(), CONTIG_NAME_MAPPING_FILE_NAME)


    def _write_contig_name_mapping(self, txt_path):
        """Write the TAB-delimited file that translates the record ids in antiSMASH's output back to contig names.

        antiSMASH's own output files (the HTML report and, with --include-detailed-output, the GenBank and JSON
        files) only know the short record ids anvi'o gave the contigs (see `_get_record_ids`), so this file is
        how a user gets from what they see there to the contig in their database.
        """

        mapping = {record_id: {"contig": contig_name} for record_id, contig_name in sorted(self.contig_of_record_id.items())}
        utils.store_dict_as_TAB_delimited_file(mapping, txt_path, key_header="antismash_record_id",
                                               headers=["antismash_record_id", "contig"])


    def _genes_txt_path(self):
        """Where the per-gene TAB-delimited file is written (inside the report directory)."""

        return os.path.join(self._report_dir(), "antismash_genes.txt")


    def _write_genes_txt(self, txt_path):
        """Write a comprehensive, self-contained per-gene summary as a flat TAB-delimited file.

        One row per antiSMASH-annotated gene, with the clear functional name and the specific biosynthetic
        domains as *separate* columns, plus the gene's role and the BGC region it belongs to. The same
        information lives in the database, split across the antiSMASH / antiSMASH_ROLE / antiSMASH_REGION
        sources; this file simply lays it out side by side for direct inspection.
        """

        # map each gene to the region it falls in (regions do not overlap, so at most one per gene).
        gene_to_region = {gene_callers_id: record["region_id"]
                          for record in self.region_records for gene_callers_id in record["gene_callers_ids"]}

        gene_rows = {}
        for gene in self.gene_annotations:
            smcog_ids = [smcog_id for smcog_id, description, e_value in gene["smcogs"]]
            smcog_descriptions = [description for smcog_id, description, e_value in gene["smcogs"]]
            domain_names = [name for name, e_value in gene["domains"]]
            e_values = [e_value for smcog_id, description, e_value in gene["smcogs"]] or [e_value for name, e_value in gene["domains"]]

            gene_rows[str(gene["gene_callers_id"])] = {
                "contig": gene["contig"],
                "region_id": gene_to_region.get(gene["gene_callers_id"], ""),
                "role": gene["role"],
                "function": _bgc_function_name(smcog_descriptions, domain_names),
                "smcog": ";".join(smcog_ids),
                "domains": ";".join(domain_names),
                "e_value": min(e_values) if e_values else "",
            }

        utils.store_dict_as_TAB_delimited_file(gene_rows, txt_path, key_header="gene_callers_id",
                                               headers=["gene_callers_id", "contig", "region_id", "role",
                                                        "function", "smcog", "domains", "e_value"])


    def _report(self):
        """Tell the user what happened."""

        self.run.info("antiSMASH version", self.antismash_version)
        self.run.info("Genes annotated", len(self.gene_annotations))
        self.run.info("BGC regions found", len(self.region_records))
        self.run.info("Genes tagged with a BGC region", sum(len(record["gene_callers_ids"]) for record in self.region_records))
        if self.gene_annotations:
            self.run.info("Genes file", os.path.abspath(self._genes_txt_path()))
        if self.region_records:
            self.run.info("Regions file", os.path.abspath(self._regions_txt_path()))
        report_label = "HTML report (detailed)" if self.include_detailed_output else "HTML report (lean)"
        self.run.info(report_label, os.path.abspath(self._report_dir()))
        self.run.info("Contig names in the report", os.path.abspath(self._contig_name_mapping_path()))
        self.run.info("antiSMASH log", os.path.abspath(os.path.join(self._report_dir(), "antismash.log")))


    @staticmethod
    def _region_id(contig, region_number):
        """Build a compact, assembly-unique identifier for a BGC region.

        antiSMASH numbers regions per contig (restarting at 1 on each contig), so the region number
        alone is not unique across an assembly. We therefore prefix it with a short hash of the contig
        name (the same 12-character sha224 slug anvi'o uses elsewhere) rather than the full contig
        name, which can be long. The full contig name for any gene is always recoverable from its gene
        caller id, and the per-region TAB-delimited file maps each region id back to its contig.
        """

        contig_slug = hashlib.sha224(contig.encode('utf-8')).hexdigest()[0:12]
        return f"{contig_slug}_region_{region_number:03d}"


    @staticmethod
    def _gene_callers_id_from_locus_tag(locus_tag):
        """Recover an anvi'o gene caller id from a locus tag written by our GFF3 export.

        Our GFF3 export uses `ID=<project_name>___<gene_callers_id>`, which we carry over into the
        `locus_tag`. Here we pull the integer tail back off. Returns None for any locus tag that does
        not look like one anvi'o produced (e.g., a gene antiSMASH added on its own)."""

        if locus_tag is None:
            return None

        tail = locus_tag.rsplit("___", 1)[-1]
        return int(tail) if tail.lstrip("-").isdigit() else None


    @staticmethod
    def _parse_location(location):
        """Turn an antiSMASH location string like '[12345:67890](+)' into a (12345, 67890) tuple."""

        body = location.strip("[]").split("]")[0]
        start, end = body.split(":")
        return int(start), int(end)


    @staticmethod
    def _parse_domain(entry):
        """Parse a `sec_met_domain` qualifier into a (name, e_value) pair.

        Example: 'ADH_N (E-value: 1.8e-14, bitscore: 51.6, seeds: 66, tool: rule-based-clusters)'
        becomes ('ADH_N', 1.8e-14). The e-value defaults to 0.0 if it cannot be found.
        """

        name = entry.split(" (", 1)[0].strip()
        match = re.search(r"E-value:\s*([0-9.eE+-]+)", entry)
        return name, (float(match.group(1)) if match else 0.0)


    @staticmethod
    def _parse_smcog(entry):
        """Parse a smCOG `gene_functions` qualifier into a (smcog_id, description, e_value) tuple.

        Example: 'biosynthetic-additional (smcogs) SMCOG1028:crotonyl-CoA reductase / alcohol dehydrogenase
        (Score: 279.1; E-value: 7.4e-85)' becomes ('SMCOG1028', 'crotonyl-CoA reductase / alcohol
        dehydrogenase', 7.4e-85).
        """

        after_tag = entry.split("(smcogs)", 1)[-1].strip()
        smcog_id, _, remainder = after_tag.partition(":")
        description = remainder.split(" (Score", 1)[0].strip().rstrip(",").strip()
        match = re.search(r"E-value:\s*([0-9.eE+-]+)", entry)
        return smcog_id.strip(), description, (float(match.group(1)) if match else 0.0)


class AntismashMatrix(object):
    """Build cross-genome matrices from antiSMASH annotations stored in contigs databases.

    Compares BGCs across several contig-dbs that have been annotated with `anvi-run-antismash`. The script
    accepts external-genomes or internal-genomes file, reads their antiSMASH function sources from the
    databases and pivots them into matrices across the databases. In
    `--genes` mode the rows are distinct gene signatures (product, role, function, smCOG, domains); in
    `--regions` mode the rows are BGC product types and the cells are cluster counts. Two matrices are
    written for each mode: raw frequencies and presence/absence.

    Because a matrix is a cross-genome summary (each cell is a count), it cannot carry the identifiers that
    are meaningful only inside a single database (the gene caller id and the region id). So in `--genes`
    mode we additionally write a combined per-gene 'bridge' file that keeps those identifiers; this is what
    lets a user drop from any matrix cell down to the individual genes and out to their other annotation
    sources. The shared `product` column ties the --genes and --regions matrices together.
    """

    def __init__(self, args, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.external_genomes_path = A('external_genomes')
        self.internal_genomes_path = A('internal_genomes')
        self.output_file_prefix = A('output_file_prefix')
        self.report_genes = A('genes') or False
        self.report_regions = A('regions') or False

        self.sanity_check()


    def sanity_check(self):
        if not self.external_genomes_path and not self.internal_genomes_path:
            raise ConfigError("Looks like you want to compare BGCs across a bunch of contig-dbs but you didn't let "
                              "the script know which contig-dbs to compare. What you need is to point the script to "
                              "your external or internal genomes file? Don't have one yet? Simple, just use "
                              "`anvi-script-gen-genomes-file`.")
        if self.report_genes == self.report_regions:
            raise ConfigError("With this script, you need to decide whether to compare genes of BGCs or entire "
                              "regions, each of which contains several genes. Use either the `--genes` or the "
                              "`--regions` flag, and if you want both, generate one after the other :).")
        if not self.output_file_prefix:
            raise ConfigError("Please use `-O` to define an output prefix. We could generate a generic prefix, "
                              "but it can be surprisingly useful to have meaningful names.")

        filesnpaths.is_output_dir_writable(os.path.dirname(os.path.abspath(self.output_file_prefix)))


    def process(self):
        genomes = self._get_genomes()
        self.run.info("Genomes to include", len(genomes))
        self.run.info("Matrix type", "genes" if self.report_genes else "regions")

        # databases never run through anvi-run-antismash get NA (not 0) in the matrix; warn about them first.
        unannotated = self._unannotated_genomes(genomes)
        if unannotated:
            self._warn_about_unannotated(unannotated)

        if self.report_genes:
            # collect every annotated gene once; both the matrix (by grouping on the signature columns) and
            # the per-gene 'bridge' file (which keeps the gene caller ids) are derived from the same records.
            genes = self._collect_genes(genomes)
            descriptor_headers = ["product", "role", "function", "smcog", "domains"]
            descriptors, counts = self._signatures_from_genes(genes, descriptor_headers)
            self._write_matrices(list(genomes.keys()), descriptor_headers, descriptors, counts, unannotated)
            self._write_genes_bridge_file(genes)
        else:
            descriptor_headers = ["product"]
            descriptors, counts = self._collect_region_products(genomes)
            self._write_matrices(list(genomes.keys()), descriptor_headers, descriptors, counts, unannotated)


    def _get_genomes(self):
        """Resolve the target genomes from -e (external) and/or -i (internal) into an ordered mapping.

        Each value is a dict with `contigs_db_path` and `gene_caller_ids`. For an external genome (a whole
        contigs database) `gene_caller_ids` is None, meaning 'use every gene'. For an internal genome (a bin)
        it is the set of gene caller ids that make up that bin, which is how we restrict the antiSMASH counts
        to a single bin's genes rather than the whole database it happens to share with other bins.
        """

        # a bin's gene list is only worked out when anvi'o initializes the genomes (init=True); for an
        # external-genomes-only run we skip that step, since there the whole database is the genome.
        genome_descriptions = GenomeDescriptions(self.args, run=terminal.Run(verbose=False), progress=terminal.Progress(verbose=False))
        genome_descriptions.load_genomes_descriptions(skip_functions=True, init=bool(self.internal_genomes_path))

        genomes = {}
        for genome_name in genome_descriptions.genomes:
            entry = genome_descriptions.genomes[genome_name]
            is_internal = genome_name in genome_descriptions.internal_genome_names
            genomes[genome_name] = {"contigs_db_path": entry['contigs_db_path'],
                                    "gene_caller_ids": set(entry['gene_caller_ids']) if is_internal else None}
        return genomes


    def _unannotated_genomes(self, genomes):
        """Return, in input order, the names of genomes that were never run through anvi-run-antismash.

        anvi-run-antismash always registers all antiSMASH sources (even when it finds no clusters), so a
        database that has none of them registered was never annotated. Such genomes get NA
        (unknown) in the matrix, as opposed to 0 (annotated, but no cluster of that kind).
        """

        # several internal genomes (bins) can share one contigs database, so we check each database only once
        # and reuse the answer for every genome that lives in it.
        annotated = {}
        unannotated = []
        for genome_name, genome in genomes.items():
            contigs_db_path = genome["contigs_db_path"]
            if contigs_db_path not in annotated:
                contigs_db = ContigsDatabase(contigs_db_path, run=terminal.Run(verbose=False), progress=terminal.Progress(verbose=False))
                annotated[contigs_db_path] = bool(set(contigs_db.meta['gene_function_sources'] or []).intersection(ANTISMASH_SOURCES))
                contigs_db.disconnect()
            if not annotated[contigs_db_path]:
                unannotated.append(genome_name)

        return unannotated


    def _warn_about_unannotated(self, unannotated):
        """Let the user know which databases lack antiSMASH annotations and that they will be scored NA."""

        self.run.warning(f"Seems like not all of your databases contain antiSMASH annotations: "
                         f"({', '.join(unannotated)}). Without those, these contig-dbs will be included in the BGC "
                         f"comparison matrix, but their columns will contain NAs instead of useful numbers. Please "
                         f"run anvi-run-antismash first, and then use this script to generate a matrix comparing a "
                         f"set of contig-dbs.")


    def _antismash_rows(self, contigs_db_path, sources):
        """Read the gene_functions rows for the requested antiSMASH sources from one contigs database.

        This only reads; callers restrict the rows to a single bin's genes themselves, which means the read
        can be cached per database when several bins share one (see `_collect_genes`).
        """

        contigs_db = ContigsDatabase(contigs_db_path, run=terminal.Run(verbose=False), progress=terminal.Progress(verbose=False))
        # smart_get builds the `source IN (...)` clause; error_if_no_data=False lets a database with none
        # of these sources (e.g. one not yet annotated) return no rows instead of raising.
        rows = contigs_db.db.smart_get(t.gene_function_calls_table_name, 'source', sources, error_if_no_data=False)
        contigs_db.disconnect()
        return list(rows.values())


    def _collect_genes(self, genomes):
        """Read every antiSMASH-annotated gene from all the genomes into one flat list of records.

        Each record is a dict with: genome, gene_callers_id, contig, region_id, product, role, function,
        smcog, domains. This is the single source of truth for the --genes output: the matrix is built by
        grouping these records on their descriptor columns, and the per-gene 'bridge' file is just these
        records written out verbatim (it keeps the gene caller ids that the matrix has to collapse away).
        """

        genes = []
        # several internal genomes (bins) can share one contigs database, so read each database once and
        # reuse it; each genome then keeps only the rows for its own genes (all of them, for an external one).
        rows_of_db, contig_of_gene_of_db = {}, {}
        for genome_name, genome in genomes.items():
            contigs_db_path = genome["contigs_db_path"]
            if contigs_db_path not in rows_of_db:
                rows_of_db[contigs_db_path] = self._antismash_rows(contigs_db_path, ANTISMASH_SOURCES)
                contig_of_gene_of_db[contigs_db_path] = self._contig_of_gene(contigs_db_path)

            gene_caller_ids = genome["gene_caller_ids"]
            # assemble the facets of each gene from the four per-gene antiSMASH sources. Region membership
            # (region id + product) comes from antiSMASH_REGION; the rest from the other three sources.
            per_gene = defaultdict(lambda: {"role": "", "smcog_ids": [], "smcog_descriptions": [],
                                            "domains": [], "region_id": "", "product": ""})
            for row in rows_of_db[contigs_db_path]:
                if gene_caller_ids is not None and row['gene_callers_id'] not in gene_caller_ids:
                    continue  # an internal genome (bin) keeps only its own genes
                gene = per_gene[row['gene_callers_id']]
                if row['source'] == ANTISMASH_ROLE_SOURCE:
                    gene["role"] = row['function']
                elif row['source'] == ANTISMASH_DOMAIN_SOURCE:
                    gene["domains"].append(row['accession'])
                elif row['source'] == ANTISMASH_REGION_SOURCE:
                    gene["region_id"] = row['accession']
                    gene["product"] = row['function']
                elif row['source'] == ANTISMASH_SOURCE and row['accession'].startswith("SMCOG"):
                    gene["smcog_ids"].append(row['accession'])
                    gene["smcog_descriptions"].append(row['function'])

            contig_of_gene = contig_of_gene_of_db[contigs_db_path]
            for gene_callers_id, gene in per_gene.items():
                domains = sorted(gene["domains"])
                genes.append({
                    "genome": genome_name,
                    "gene_callers_id": gene_callers_id,
                    "contig": contig_of_gene.get(gene_callers_id, ""),
                    "region_id": gene["region_id"],
                    "product": gene["product"],
                    "role": gene["role"],
                    "function": _bgc_function_name(gene["smcog_descriptions"], domains),
                    "smcog": ";".join(sorted(gene["smcog_ids"])),
                    "domains": ";".join(domains),
                })

        return genes


    def _contig_of_gene(self, contigs_db_path):
        """Return a {gene_callers_id: contig name} mapping for one contigs database."""

        contigs_db = ContigsDatabase(contigs_db_path, run=terminal.Run(verbose=False), progress=terminal.Progress(verbose=False))
        rows = contigs_db.db.get_some_columns_from_table(t.genes_in_contigs_table_name, "gene_callers_id, contig")
        contigs_db.disconnect()
        return dict(rows)


    def _signatures_from_genes(self, genes, descriptor_headers):
        """Group the per-gene records into distinct signatures and count them per genome.

        The signature includes the BGC product type (descriptor_headers[0]), so each --genes row is tied
        to a cluster type, which is what lets it line up with the --regions matrix. A gene signature that
        occurs in more than one product type becomes more than one row.

        Genes that only sit inside a cluster boundary without any functional annotation of their own (every
        facet after `product` is empty) are left out: they are not a 'signature', and would otherwise pile
        up into one meaningless blank-descriptor row per product type. They remain in the per-gene bridge
        file, which is the comprehensive record of every gene in every cluster.
        """

        facet_headers = descriptor_headers[1:]

        descriptors, counts = {}, defaultdict(lambda: defaultdict(int))
        for gene in genes:
            if not any(gene[header] for header in facet_headers):
                continue
            key = tuple(gene[header] for header in descriptor_headers)
            descriptors[key] = {header: gene[header] for header in descriptor_headers}
            counts[key][gene["genome"]] += 1

        return descriptors, counts


    def _collect_region_products(self, genomes):
        """For --regions: build {product: descriptor} and {product: {genome: number of clusters}}."""

        descriptors, counts = {}, defaultdict(lambda: defaultdict(int))
        region_rows_of_db = {}   # read each database once, even when several bins share it
        for genome_name, genome in genomes.items():
            contigs_db_path = genome["contigs_db_path"]
            if contigs_db_path not in region_rows_of_db:
                region_rows_of_db[contigs_db_path] = self._antismash_rows(contigs_db_path, [ANTISMASH_REGION_SOURCE])

            gene_caller_ids = genome["gene_caller_ids"]
            # the antiSMASH_REGION source has one row per gene; collapse to distinct regions, since a region's
            # accession (its id) repeats across all the genes it contains. For an internal genome we keep only
            # the bin's genes, and because a BGC region sits entirely on one contig (hence in one bin), this
            # cleanly counts each region for the bin that owns it.
            product_of_region = {row['accession']: row['function'] for row in region_rows_of_db[contigs_db_path]
                                 if gene_caller_ids is None or row['gene_callers_id'] in gene_caller_ids}
            for product in product_of_region.values():
                descriptors[product] = {"product": product}
                counts[product][genome_name] += 1   # one count per distinct region -> a cluster count

        return descriptors, counts


    def _write_genes_bridge_file(self, genes):
        """Write the combined per-gene table that bridges the matrices back to individual genes.

        One row per antiSMASH-annotated gene across all the genomes. It carries the two identifiers that are
        local to a single database and therefore cannot live in a cross-genome matrix (the gene_callers_id
        and the region_id) next to the same product / role / function / smCOG / domains descriptors the
        --genes matrix groups on. This is the join key for everything downstream: `genome` + `gene_callers_id`
        points back to that genome's other annotation sources (KOfam, COGs, ...); `region_id` groups genes
        into their clusters; and `product` plus the signature columns line up with the --genes and --regions
        matrices.
        """

        headers = ["genome", "gene_callers_id", "contig", "region_id", "product", "role", "function", "smcog", "domains"]
        output_path = f"{os.path.abspath(self.output_file_prefix)}-GENE-DETAILS.txt"
        filesnpaths.is_output_file_writable(output_path)

        with open(output_path, "w") as output:
            output.write("\t".join(headers) + "\n")
            for gene in sorted(genes, key=lambda gene: (gene["genome"], gene["gene_callers_id"])):
                output.write("\t".join(str(gene[header]) for header in headers) + "\n")

        self.run.info("Per-gene details", output_path)


    def _write_matrices(self, genome_names, descriptor_headers, descriptors, counts, unannotated):
        """Write the frequency and presence/absence matrices (descriptor columns first, then one per genome).

        Genomes in `unannotated` (no antiSMASH source stored at all) get NA in every cell rather than a
        count: their biosynthetic content is unknown, not measured-as-zero.
        """

        unannotated = set(unannotated)
        views = [("FREQUENCY", lambda count: str(count)),
                 ("PRESENCE-ABSENCE", lambda count: "1" if count else "0")]

        for view_name, cell in views:
            output_path = f"{os.path.abspath(self.output_file_prefix)}-{view_name}.txt"
            filesnpaths.is_output_file_writable(output_path)
            with open(output_path, "w") as output:
                output.write("\t".join(descriptor_headers + genome_names) + "\n")
                for key in sorted(descriptors):
                    descriptor_columns = [str(descriptors[key][header]) for header in descriptor_headers]
                    value_columns = ["NA" if genome_name in unannotated else cell(counts[key].get(genome_name, 0))
                                     for genome_name in genome_names]
                    output.write("\t".join(descriptor_columns + value_columns) + "\n")

            self.run.info(f"Matrix ({view_name.lower()})", output_path)
