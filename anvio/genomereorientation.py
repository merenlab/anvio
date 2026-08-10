#!/usr/bin/env python
"""Reorient circular genomes to match a reference genome coordinate system."""

import os
import gzip
import shutil

from pathlib import Path

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.drivers.hmmer import HMMer
from anvio.genecalling import GeneCaller
from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2025, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"

P = terminal.pluralize


class PafRecord:
    def __init__(self, qname, qlen, qstart, qend, strand, tname, tlen, tstart, tend, nmatch, alen, mapq, tags):
        self.qname = qname
        self.qlen = qlen
        self.qstart = qstart
        self.qend = qend
        self.strand = strand
        self.tname = tname
        self.tlen = tlen
        self.tstart = tstart
        self.tend = tend
        self.nmatch = nmatch
        self.alen = alen
        self.mapq = mapq
        self.tags = tags

    @property
    def is_primary(self):
        return "tp:A:P" in self.tags

    @property
    def aligned_bases(self):
        return self.alen


class ReorientationResult:
    def __init__(self, name, status, message, output_path=None, trust=None, alignment=None):
        self.name = name
        self.status = status
        self.message = message
        self.output_path = output_path
        self.trust = trust
        self.alignment = alignment  # Store alignment info for post-processing


class GenomeReorienter:
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        # We need these two for this to work.
        filesnpaths.is_program_exists("minimap2", advice_if_not_exists="You should be able to install it in your "
                                      "environment using 'conda install -c bioconda minimap2'.")
        filesnpaths.is_program_exists("seqkit", advice_if_not_exists="You should be able to install it in your "
                                      "environment using 'conda install -c bioconda seqkit'.")

        # Get all the things
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.fasta_txt = os.path.abspath(A('fasta_txt'))
        self.reference_name = A('reference')
        self.output_dir = os.path.abspath(A('output_dir'))
        self.threads = A('threads') or 1
        self.minimap2_preset = A('minimap2_preset') or "asm5"
        self.near_start_bp = A('near_start_bp') or 10000
        self.use_dnaa_for_reference_orientation = A('use_dnaa_for_reference_orientation') or False
        self.use_auto_reference_as_is = A('use_auto_reference_as_is') or False
        self.scaffold_fragmented = A('scaffold_fragmented') or False
        self.min_contig_length = A('min_contig_length') or 1000
        self.keep_query_contigs_intact = A('keep_query_contigs_intact') or False

        # Visualization parameters
        self.skip_visualizing_alignments = A('skip_visualizing_alignments') or False
        self.plot_width = A('plot_width')  # None means auto-detect terminal width
        self.plot_height = A('plot_height') or 20

        # If the user wants extensive logging, we will take care of it through
        # terminal.Run:
        log_file_path = A('log_file_path')
        if log_file_path:
            filesnpaths.is_output_file_writable(log_file_path)
        self.log_run = terminal.Run(log_file_path=log_file_path, verbose=False)

        self.genomes = {}
        self.reference_path = None

        self.sanity_check()


    def process(self):
        self.run.info("Output directory", self.output_dir)

        # Deal with reference selection
        reference_was_user_specified = self.reference_name is not None
        self.select_reference_genome()

        # Decide how to orient the reference
        if self.use_dnaa_for_reference_orientation:
            # Use DnaA gene to orient the reference
            dnaa_position = self._find_dnaa_gene_position()
            if dnaa_position is not None and dnaa_position != 0:
                self.run.warning(None, header="ROTATING REFERENCE TO DnaA POSITION")
                self.run.info("Rotating reference by", f"{dnaa_position:,} nts")
                self.run.info("Reason", "Aligning genome to DnaA gene (origin of replication)")

                # Create a rotated version of the reference
                temp_dir = filesnpaths.get_temp_directory_path()
                rotated_ref = os.path.join(temp_dir, f"{self.reference_name}_rotated_dnaa.fa")
                self._seqkit_rotate(self.reference_path, dnaa_position + 1, rotated_ref)

                # Update reference path to the rotated version
                self.reference_path = rotated_ref
                self.run.info("Rotated reference", rotated_ref, nl_after=1)
            elif dnaa_position == 0:
                self.run.info("Reference rotation", "Not needed (DnaA gene is already at position 0 .. yes, that happened o_O)", nl_after=1)

        # If reference was auto-selected and not using DnaA, find the optimal starting position
        elif not reference_was_user_specified and not self.use_auto_reference_as_is:
            result = self._find_optimal_reference_start()
            if result:
                optimal_position, genomes_covering, total_genomes = result
                coverage_pct = (genomes_covering / total_genomes) * 100 if total_genomes > 0 else 0

                if optimal_position != 0:
                    self.run.warning(None, header="ROTATING REFERENCE TO OPTIMAL START")
                    self.run.info("Rotating reference by", f"{optimal_position:,} nts")
                    self.run.info("Reason", f"This position is conserved across {genomes_covering}/{total_genomes} genomes ({coverage_pct:.1f}%)")

                    # Create a rotated version of the reference
                    temp_dir = filesnpaths.get_temp_directory_path()
                    rotated_ref = os.path.join(temp_dir, f"{self.reference_name}_rotated.fa")
                    self._seqkit_rotate(self.reference_path, optimal_position + 1, rotated_ref)

                    # Update reference path to the rotated version
                    self.reference_path = rotated_ref
                    self.run.info("Rotated reference", rotated_ref, nl_after=1)
                else:
                    self.run.info("Reference rotation", "Not needed (optimal start is already at position 0)", nl_after=1)

        if self.log_run.log_file_path:
            self.run.info("Log file", self.log_run.log_file_path)

        filesnpaths.gen_output_directory(self.output_dir, run=self.run)

        reference_output_path = self._output_path_for(self.reference_name)
        shutil.copyfile(self.reference_path, reference_output_path)

        rotation_msg = "Copied reference genome without reorientation."
        if self.use_dnaa_for_reference_orientation:
            rotation_msg = "Reference genome (rotated to DnaA gene position)."
        elif not reference_was_user_specified and not self.use_auto_reference_as_is:
            rotation_msg = "Reference genome (possibly rotated to optimal start position)."

        results = [ReorientationResult(self.reference_name, "reference",
                                       rotation_msg,
                                       reference_output_path,
                                       trust="REFERENCE")]

        total_to_process = len(self.genomes) - 1
        processed = 0

        if total_to_process > 0:
            self.progress.new("Reorienting genomes", progress_total_items=total_to_process)

        for genome_name, entry in self.genomes.items():
            if genome_name == self.reference_name:
                continue

            processed += 1
            self.progress.clear()
            self.run.warning(None, header=f"REORIENTING {genome_name} ({processed} OF {total_to_process} TOTAL)")

            fasta_path = entry["path"]
            output_path = self._output_path_for(genome_name)
            num_contigs = entry['num_contigs']

            try:
                if num_contigs == 1:
                    # We're working with a ircular genome
                    best_alignment, paf_final, paf_initial, actions = self._process_circular(fasta_path, output_path, genome_name)
                    result = self._report_and_visualize_circular_genome(genome_name, output_path, best_alignment, paf_final, paf_initial, actions, fasta_path)
                    results.append(result)
                else:
                    # We're working with a fragmented genome
                    result_data = self._process_fragmented(fasta_path, output_path, genome_name, num_contigs)
                    result = self._report_and_visualize_fragmented_genome(genome_name, output_path, result_data, fasta_path)
                    results.append(result)
            except (ConfigError, FilesNPathsError, RuntimeError) as e:
                self.progress.clear()
                self.run.info(f"{genome_name} reorientation", f"failed ({e})", mc="red")
                results.append(ReorientationResult(genome_name, "failed", str(e), trust="FAILED"))
            finally:
                if total_to_process > 0:
                    self.progress.increment()

        if total_to_process > 0:
            self.progress.end()

        self._final_report(results)

        failures = [r for r in results if r.status == "failed"]
        if failures:
            raise ConfigError("Reorientation finished with failures. Please review the log output for details.")

        citation_msg = ("Anvi'o used `minimap2` by Li (doi:10.1093/bioinformatics/bty191) to rapidly align genomes "
                       "and `seqkit` by Shen et al (doi:10.1371/journal.pone.0163962) to rotate and reverse/complement "
                       "large sequences.")

        if self.use_dnaa_for_reference_orientation:
            citation_msg += ("Additionally, anvi'o used 'pyrodigal-gv' by Martin Larralde for gene calling. It is an extension of "
                             "'pyrodigal' (doi:10.21105/joss.04296), which builds upon the approach originally implemented by Hyatt et al "
                             "(doi:10.1186/1471-2105-11-119), with additional metagenomics models for giant viruses and viruses with "
                             "alternative genetic codes by Camargo et al (doi:10.1038/s41587-023-01953-y). Anvi'o also used `HMMER` by "
                             "Eddy (doi:10.1371/journal.pcbi.1002195) to identify the DnaA gene for reference orientation.")

        citation_msg += " If you publish your findings, please do not forget to properly credit the authors of these tools."

        self.run.warning(citation_msg, lc='cyan', header="CITATION")

        return results


    def sanity_check(self):
        filesnpaths.is_file_tab_delimited(self.fasta_txt, expected_number_of_fields=2)
        self.genomes = utils.get_TAB_delimited_file_as_dictionary(
            self.fasta_txt,
            expected_fields=['name', 'path'],
            only_expected_fields=True)

        if not len(self.genomes):
            raise ConfigError(f"The fasta-txt file '{self.fasta_txt}' seems to be empty :/")

        if self.reference_name and self.reference_name not in self.genomes:
            raise ConfigError(f"Reference '{self.reference_name}' is not present in fasta-txt '{self.fasta_txt}'.")

        # Get the directory containing fasta.txt for resolving relative paths
        fasta_txt_dir = os.path.dirname(self.fasta_txt)

        for genome_name, entry in self.genomes.items():
            utils.is_this_name_OK_for_database('fasta.txt entry name', genome_name, additional_chars_allowed='.-')

            # Resolve path relative to fasta.txt location if it's not absolute
            if os.path.isabs(entry['path']):
                genome_path = entry['path']
            else:
                genome_path = os.path.join(fasta_txt_dir, entry['path'])

            genome_path = os.path.abspath(genome_path)
            filesnpaths.is_file_exists(genome_path)
            filesnpaths.is_file_fasta_formatted(genome_path)

            num_sequences = utils.get_num_sequences_in_fasta(genome_path)
            self.genomes[genome_name]['num_contigs'] = num_sequences
            self.genomes[genome_name]['path'] = genome_path

        if self.use_auto_reference_as_is and self.reference_name:
            raise ConfigError("9 9 9. You can't use `--reference` with `--use-auto-reference-as-is`. You have to "
                              "choose one of them :/")

        if self.use_auto_reference_as_is and self.use_dnaa_for_reference_orientation:
            raise ConfigError("Well .. `--use-auto-reference-as-is` and `--use-dnaa-for-reference-orientation` are "
                              "not compatible with one another. The former keeps the auto-selected reference untouched. "
                              "The latter rotates it. Incompatible stuff.")

        # please note that whether the reference genome is a single contig or not is checked in
        # `check_reference_is_single_contig`, which is called by `select_reference_genome` so the
        # very same rule applies to user-specified and auto-selected references alike.

        self.output_dir = filesnpaths.check_output_directory(self.output_dir, ok_if_exists=False)

        if self.threads < 1:
            raise ConfigError("Number of threads must be a positive integer.")

        if self.min_contig_length < 0:
            raise ConfigError("--min-contig-length must be a non-negative integer.")


    def _output_path_for(self, genome_name):
        suffix = Path(self.genomes[genome_name]['path']).suffix or ".fa"
        return os.path.join(self.output_dir, f"{genome_name}{suffix}")


    def _report_and_visualize_circular_genome(self, genome_name, output_path, best_alignment, paf_final, paf_initial, actions, fasta_path):
        """Report statistics and visualize circular genome reorientation results."""
        # Quality assessment for circular genomes
        diff = self._get_dv_from_tags(best_alignment.tags)
        if diff is not None:
            approx_ani = (1 - diff) * 100
        elif best_alignment.aligned_bases > 0:
            approx_ani = (best_alignment.nmatch / float(best_alignment.aligned_bases)) * 100
        else:
            approx_ani = 0

        cov_q = 0
        cov_t = 0
        if best_alignment.qlen and best_alignment.tlen:
            cov_q_raw = (best_alignment.aligned_bases / float(best_alignment.qlen)) * 100
            cov_t_raw = (best_alignment.aligned_bases / float(best_alignment.tlen)) * 100
            cov_q = min(cov_q_raw, 100)
            cov_t = min(cov_t_raw, 100)

        avg_cov = (cov_q + cov_t) / 2.0
        if avg_cov < 50:
            trust_label = "NOT TRUSTWORTHY"
            trust_color = "red"
        elif avg_cov < 90:
            trust_label = "SOMEWHAT OK"
            trust_color = "yellow"
        else:
            trust_label = "TRUSTWORTHY"
            trust_color = "green"

        message = (f"Final alignment strand={best_alignment.strand} "
                   f"qstart={best_alignment.qstart} tstart={best_alignment.tstart} "
                   f"alen={best_alignment.aligned_bases} "
                   f"approx_ani={approx_ani:.1f}%")

        self.progress.clear()
        self.run.info("Orientation outcome", trust_label, mc=trust_color)
        self.run.info("Applied actions", "")
        for action in actions:
            self.run.info_single(action, level=2)
        self.run.info("Output FASTA", output_path)
        self.run.info("Final alignment strand", best_alignment.strand)
        self.run.info("Start in query", best_alignment.qstart)
        self.run.info("Start in reference", best_alignment.tstart)
        self.run.info("Query length", best_alignment.qlen)
        self.run.info("Reference length", best_alignment.tlen)
        self.run.info("Aligned length", best_alignment.aligned_bases)
        self.run.info("Query coverage by alignment", f"{cov_q:.1f}%")
        self.run.info("Reference coverage by alignment", f"{cov_t:.1f}%")
        self.run.info("Approx ANI to reference", f"{approx_ani:.1f}%", nl_after=1)

        # Show before and after visualizations
        if not self.skip_visualizing_alignments:
            self.run.info_single("Before reorientation", nl_after=1)
            self._plot_synteny_ribbons(paf_initial, genome_name, label="Before reorientation")
            self.run.info_single("After reorientation", nl_before=1, nl_after=1)
            self._plot_synteny_ribbons(paf_final, genome_name, label="After reorientation")

        return ReorientationResult(genome_name, "success", message, output_path, trust=trust_label, alignment=best_alignment)


    def _report_and_visualize_fragmented_genome(self, genome_name, output_path, result_data, fasta_path):
        """Report statistics and visualize fragmented genome reorientation results."""
        # Quality assessment for fragmented genomes
        ref_cov = result_data['reference_coverage_pct']
        avg_ani = result_data['avg_ani']
        num_aligned = result_data['num_contigs_aligned']
        num_total = result_data['num_contigs_processed']

        # Trust criteria for fragmented genomes
        aligned_pct = (num_aligned / num_total) * 100 if num_total > 0 else 0

        if ref_cov < 50 or aligned_pct < 50:
            trust_label = "NOT TRUSTWORTHY"
            trust_color = "red"
        elif ref_cov < 70 or aligned_pct < 80 or avg_ani < 95:
            trust_label = "SOMEWHAT OK"
            trust_color = "yellow"
        else:
            trust_label = "TRUSTWORTHY"
            trust_color = "green"

        message = (f"Scaffolded {num_aligned}/{num_total} contigs, "
                   f"ref_coverage={ref_cov:.1f}%, avg_ani={avg_ani:.1f}%")

        # Report fragmented genome-specific metrics
        self.progress.clear()
        self.run.info("Scaffolding outcome", trust_label, mc=trust_color)
        self.run.info("Applied actions", "")
        self.run.info_single(result_data['actions_summary'], level=2)
        self.run.info("Output FASTA", output_path)
        self.run.info("Contigs in input", result_data['num_contigs_in_input'])
        if result_data['num_contigs_rotated']:
            self.run.info("Circularly permuted contigs rotated", result_data['num_contigs_rotated'], mc="yellow")
        if result_data['num_contigs_split']:
            self.run.info("Contigs cut at reference boundaries", f"{result_data['num_contigs_split']} (into {result_data['num_fragments']} fragments)", mc="yellow")
        self.run.info("Contigs processed", num_total)
        self.run.info("Contigs aligned", num_aligned)
        self.run.info("Contigs unaligned", result_data['num_contigs_unaligned'])
        self.run.info("Reference coverage", f"{ref_cov:.1f}%")
        self.run.info("Average ANI", f"{avg_ani:.1f}%")
        self.run.info("Total gaps", len(result_data['gaps']))
        self.run.info("Total gap size", f"{result_data['total_gap_size']:,} nts", nl_after=1)

        # Generate alignment plots showing before and after reorientation
        if not self.skip_visualizing_alignments:
            self.progress.update(f"{genome_name}: Generating alignment plots")
            try:
                # Before
                initial_paf = self._minimap2_align(self.reference_path, fasta_path)
                if initial_paf:
                    self.progress.reset()
                    self.run.info_single("Before reorientation", nl_after=1)
                    self._plot_synteny_ribbons(initial_paf, genome_name, label="Before reorientation")

                # After
                reoriented_paf = self._minimap2_align(self.reference_path, output_path)
                if reoriented_paf:
                    self.progress.reset()
                    self.run.info_single("After reorientation", nl_before=1, nl_after=1)
                    self._plot_synteny_ribbons(reoriented_paf, genome_name, label="After reorientation")
            except Exception as e:
                self.log_run.info_single(f"Could not generate alignment plots: {e}", level=2)

        return ReorientationResult(genome_name, "success", message, output_path, trust=trust_label, alignment=None)


    def _process_circular(self, query_fa, output_path, genome_name):
        temp_dir = filesnpaths.get_temp_directory_path()
        self.log_run.info_single(f"Working directory for intermediates: {temp_dir}", level=2)

        try:
            actions = []

            self.progress.update(f"{genome_name}: Initial alignment")
            paf1 = self._minimap2_align(self.reference_path, query_fa)
            anchor1 = self._select_anchor_near_reference_start(paf1, self.near_start_bp)
            cut0 = self._cut0_for_ref0(anchor1)
            self.log_run.info_single(f"Pass-1 anchor strand={anchor1.strand} cut0={cut0}", level=2)
            if anchor1.strand == "-":
                actions.append("reverse-complemented query (pass1 anchor '-')")
            actions.append(f"rotated {cut0 + 1} nts (pass1 snap)")

            self.progress.update(f"{genome_name}: Orienting and rotating (pass 1)")
            oriented_query = os.path.join(temp_dir, "query_oriented.fa")
            if anchor1.strand == "-":
                self._seqkit_reverse_complement(query_fa, oriented_query)
                # After reverse complement, positions transform: i -> qlen-1-i
                # So the cut position must be transformed too
                cut0 = (anchor1.qlen - 1 - cut0) % anchor1.qlen
                self.log_run.info_single(f"Transformed cut0={cut0} after RC", level=2)
            else:
                shutil.copyfile(query_fa, oriented_query)

            # After the initial pass we go straight to verification and correction
            initial_rotated = os.path.join(temp_dir, "query_rot1.fa")
            self._seqkit_rotate(oriented_query, cut0 + 1, initial_rotated)

            self.progress.update(f"{genome_name}: Verification alignment")
            paf_final = self._minimap2_align(self.reference_path, initial_rotated)
            best_final_final = self._best_primary_alignment(paf_final)

            # Iteratively correct until alignment starts at reference position 0
            current_fa = initial_rotated
            max_correction_iterations = 5
            correction_iteration = 0

            self.log_run.info_single(f"Initial: tstart={best_final_final.tstart} qstart={best_final_final.qstart}", level=2)

            # Final correction to ensure query starts at position 0
            # Goal: qstart=0 (tstart offset will be handled by consensus offset fix if needed)
            while best_final_final.qstart != 0 and correction_iteration < max_correction_iterations:
                correction_iteration += 1
                self.progress.update(f"{genome_name}: Correction iteration {correction_iteration}")
                self.log_run.info_single(f"Correction iteration {correction_iteration}: tstart={best_final_final.tstart} qstart={best_final_final.qstart} strand={best_final_final.strand}", level=2)

                # Only fix qstart, not tstart
                # tstart offset will be handled globally by consensus offset fix
                if best_final_final.qstart != 0:
                    if best_final_final.strand == "+":
                        cut_correction = best_final_final.qstart
                    else:
                        cut_correction = self._cut0_for_ref0(best_final_final)
                else:
                    # qstart is already 0, we're done
                    break

                if cut_correction == 0:
                    self.log_run.info_single("Cut correction is 0, alignment already optimal", level=2)
                    break

                corrected_fa = os.path.join(temp_dir, f"query_corrected_{correction_iteration}.fa")
                self._seqkit_rotate(current_fa, cut_correction + 1, corrected_fa)
                actions.append(f"rotated {cut_correction + 1} bp (correction iter {correction_iteration} for qstart={best_final_final.qstart}, tstart={best_final_final.tstart})")

                # Re-align after correction
                paf_final = self._minimap2_align(self.reference_path, corrected_fa)
                best_final_final = self._best_primary_alignment(paf_final)
                current_fa = corrected_fa

                self.log_run.info_single(f"After correction {correction_iteration}: tstart={best_final_final.tstart} qstart={best_final_final.qstart}", level=2)

            if best_final_final.tstart == 0:
                self.log_run.info_single(f"Successfully aligned to tstart=0 after {correction_iteration} correction(s)", level=2)
                shutil.copyfile(current_fa, output_path)
            else:
                # qstart=0 but tstart!=0: query is offset from reference
                # query[0] aligns with ref[tstart], meaning the first tstart bp of ref are at the END of query
                # Rotate BACKWARD to move those bp from end to beginning
                self.log_run.info_single(f"Final tstart={best_final_final.tstart}, applying backward rotation to match reference", level=2)
                final_adjustment_fa = os.path.join(temp_dir, "query_final_adjustment.fa")

                # Rotate backward by tstart (use qlen - tstart for forward rotation equivalent)
                backward_rotation = (best_final_final.qlen - best_final_final.tstart) % best_final_final.qlen
                self._seqkit_rotate(current_fa, backward_rotation + 1, final_adjustment_fa)
                shutil.copyfile(final_adjustment_fa, output_path)
                actions.append(f"rotated backward {best_final_final.tstart} nts (final adjustment to match reference start)")

                # Update alignment info for reporting
                paf_final = self._minimap2_align(self.reference_path, output_path)
                best_final_final = self._best_primary_alignment(paf_final)
            return best_final_final, paf_final, paf1, actions
        finally:
            if not anvio.DEBUG:
                shutil.rmtree(temp_dir)
            else:
                self.run.warning(f"Temp directory kept for debugging: {temp_dir}")


    def _process_fragmented(self, query_fa, output_path, genome_name, num_contigs):
        """
        Reorient contigs in a fragmented genome by aligning them to reference.

        Workflow:
        1. Read and filter contigs by min_contig_length
        2. Align each contig independently to reference
        3. Order contigs by reference position (tstart)
        4. Orient contigs (reverse-complement - strand)
        5. Calculate gaps between consecutive contigs
        6. Write output FASTA (with optional N-padding)
        7. Calculate quality metrics

        Returns dict with reorientation results and metrics.
        """
        temp_dir = filesnpaths.get_temp_directory_path()
        self.log_run.info_single(f"Reorienting {genome_name} ({num_contigs} contigs)", level=2)

        try:
            # Step 1: Read and filter contigs
            contigs = []
            fasta = utils.u.SequenceSource(query_fa)
            while next(fasta):
                contig_len = len(fasta.seq)
                if contig_len >= self.min_contig_length:
                    contigs.append({
                        'id': fasta.id,
                        'seq': fasta.seq,
                        'length': contig_len
                    })
                else:
                    self.log_run.info_single(
                        f"Filtered '{fasta.id}' (length={contig_len} < {self.min_contig_length})",
                        level=2)
            fasta.close()

            num_contigs_after_filter = len(contigs)
            if num_contigs_after_filter == 0:
                raise ConfigError("No contigs remain after length filtering")

            self.log_run.info_single(
                f"Contigs after filter: {num_contigs_after_filter}/{num_contigs}", level=2)

            # Step 2: Align each contig to reference
            self.progress.update(f"{genome_name}: Aligning {num_contigs_after_filter} contigs")
            ref_length = self._get_total_length(self.reference_path)

            alignments = {}
            for idx, contig in enumerate(contigs):
                alignments[contig['id']] = self._align_contig_to_reference(contig, temp_dir, f"contig_{idx}")

            # Step 3: Rotate the contigs an assembler cut out of a cycle in its assembly graph, so
            # that they stop looking scrambled and become co-linear with the reference (see
            # `_find_circular_permutation`). This happens before any cutting, since a contig that
            # has been made co-linear may no longer need to be cut at all
            num_contigs_rotated = self._uncoil_circularly_permuted_contigs(contigs, alignments, ref_length, genome_name, temp_dir)

            # Step 4: Cut the contigs that run past either end of the reference coordinate axis into
            # fragments, since there is no single place on that axis where such a contig could go in
            # one piece (see `_find_reference_boundary_cut_points` for the long story)
            contigs, num_contigs_split, num_fragments = self._split_contigs_at_reference_boundaries(contigs, alignments, ref_length, genome_name, temp_dir)

            # Step 5: Sort every contig -- and every fragment of a contig that had to be cut -- into
            # those that aligned to the reference and those that did not
            contig_alignments = {}
            unaligned_contigs = []

            for contig in contigs:
                primaries = [r for r in (alignments[contig['id']] or []) if r.is_primary]

                if not primaries:
                    self.log_run.info_single(f"'{contig['id']}': no primary alignment to the reference", level=2)
                    unaligned_contigs.append(contig)
                    continue

                # Select best primary alignment
                best = max(primaries, key=lambda r: (r.aligned_bases, r.mapq, r.nmatch))
                contig_alignments[contig['id']] = {
                    'contig_data': contig,
                    'alignment': best
                }
                self.log_run.info_single(
                    f"'{contig['id']}': ref[{best.tstart}:{best.tend}] "
                    f"strand={best.strand} alen={best.aligned_bases}", level=2)

            num_contigs_after_split = len(contigs)
            num_aligned = len(contig_alignments)
            num_unaligned = len(unaligned_contigs)

            self.log_run.info_single(
                f"Alignment summary: {num_aligned}/{num_contigs_after_split} contigs aligned "
                f"({num_contigs_split} contigs were cut into {num_fragments} fragments)",
                level=2)

            if num_unaligned > 0:
                total_unaligned_bp = sum(c['length'] for c in unaligned_contigs)
                self.run.warning(f"{num_unaligned} contig(s) ({total_unaligned_bp:,} nts) in '{genome_name}' had no "
                                 f"alignment to the reference genome and will be appended to the output FASTA as-is "
                                 f"without reorientation.",
                                 header="UNALIGNED CONTIGS DETECTED")

            if num_aligned == 0:
                raise ConfigError("No contigs aligned to reference")

            # Step 6: Order contigs by reference position
            ordered_contig_ids = sorted(
                contig_alignments.keys(),
                key=lambda cid: (
                    contig_alignments[cid]['alignment'].tstart,
                    -contig_alignments[cid]['alignment'].aligned_bases,
                    -contig_alignments[cid]['alignment'].mapq
                )
            )

            # Step 7: Orient contigs and calculate gaps
            gaps = []
            oriented_contigs = []

            for idx, contig_id in enumerate(ordered_contig_ids):
                contig_info = contig_alignments[contig_id]
                contig_data = contig_info['contig_data']
                alignment = contig_info['alignment']

                # Orient based on strand
                if alignment.strand == '-':
                    temp_in = os.path.join(temp_dir, f"orient_in_{idx}.fa")
                    temp_out = os.path.join(temp_dir, f"orient_out_{idx}.fa")
                    with open(temp_in, 'w') as f:
                        f.write(f">{contig_data['id']}\n{contig_data['seq']}\n")

                    self._seqkit_reverse_complement(temp_in, temp_out)

                    rc_fasta = utils.u.SequenceSource(temp_out)
                    next(rc_fasta)
                    oriented_seq = rc_fasta.seq
                    rc_fasta.close()

                    self.log_run.info_single(f"Reverse-complemented '{contig_id}'", level=2)
                else:
                    oriented_seq = contig_data['seq']

                oriented_contigs.append({
                    'id': contig_id,
                    'seq': oriented_seq,
                    'ref_start': alignment.tstart,
                    'ref_end': alignment.tend
                })

                # Calculate gap to next contig
                if idx < len(ordered_contig_ids) - 1:
                    next_contig_id = ordered_contig_ids[idx + 1]
                    next_alignment = contig_alignments[next_contig_id]['alignment']

                    gap_start = alignment.tend
                    gap_end = next_alignment.tstart

                    if gap_end > gap_start:
                        gap_size = gap_end - gap_start
                        gaps.append({
                            'after_contig': contig_id,
                            'gap_size': gap_size,
                            'ref_start': gap_start,
                            'ref_end': gap_end
                        })
                        self.log_run.info_single(
                            f"Gap: {gap_size} nts between '{contig_id}' and '{next_contig_id}'",
                            level=2)
                    elif gap_end < gap_start:
                        overlap_size = gap_start - gap_end
                        self.log_run.info_single(
                            f"WARNING: {overlap_size} nts overlap between '{contig_id}' "
                            f"and '{next_contig_id}'", level=2)
                        gaps.append({
                            'after_contig': contig_id,
                            'gap_size': 0,
                            'overlap': overlap_size
                        })

            # Step 8: Write output FASTA
            self.progress.update(f"{genome_name}: Writing reoriented output")

            # We will do it differently depending on user's wishes
            if self.scaffold_fragmented:
                # Write as a single concatenated sequence with N-padding
                scaffold_sequence = []
                for idx, contig_info in enumerate(oriented_contigs):
                    scaffold_sequence.append(contig_info['seq'])

                    # Add gap padding between contigs
                    if idx < len(gaps):
                        gap_info = gaps[idx]
                        if 'overlap' not in gap_info and gap_info['gap_size'] > 0:
                            num_ns = min(gap_info['gap_size'], 100000)  # Cap at 100kb
                            scaffold_sequence.append('N' * num_ns)

                # Write as single sequence
                full_scaffold = ''.join(scaffold_sequence)
                with open(output_path, 'w') as out_fa:
                    out_fa.write(f">{genome_name}_scaffolded\n")
                    # Write in 80-character lines
                    for i in range(0, len(full_scaffold), 80):
                        out_fa.write(full_scaffold[i:i+80] + '\n')
                    # Append unaligned contigs as separate entries
                    for contig in unaligned_contigs:
                        out_fa.write(f">{contig['id']}\n")
                        seq = contig['seq']
                        for i in range(0, len(seq), 80):
                            out_fa.write(seq[i:i+80] + '\n')
            else:
                # Write as separate contigs (ordered and oriented), then unaligned contigs
                with open(output_path, 'w') as out_fa:
                    for idx, contig_info in enumerate(oriented_contigs):
                        out_fa.write(f">{contig_info['id']}\n")
                        # Write in 80-character lines
                        seq = contig_info['seq']
                        for i in range(0, len(seq), 80):
                            out_fa.write(seq[i:i+80] + '\n')
                    # Append unaligned contigs as-is
                    for contig in unaligned_contigs:
                        out_fa.write(f">{contig['id']}\n")
                        seq = contig['seq']
                        for i in range(0, len(seq), 80):
                            out_fa.write(seq[i:i+80] + '\n')

            # Step 9: Calculate quality metrics
            total_gap_size = sum(g['gap_size'] for g in gaps if 'overlap' not in g)

            # Reference coverage. every primary alignment block of every placed contig counts here,
            # not just the single best block per contig: a contig often aligns to the reference in
            # more than one piece, and counting only the largest of them can understate the
            # coverage of a perfectly good genome by tens of percent (which then drags its trust
            # label down with it for no reason)
            ref_covered_regions = [(record.tstart, record.tend) for contig_id in ordered_contig_ids
                                   for record in alignments[contig_id] if record.is_primary]

            merged_regions = self._merge_intervals(ref_covered_regions)
            total_covered_bases = sum(end - start for start, end in merged_regions)
            reference_coverage_pct = (total_covered_bases / ref_length) * 100

            # Average ANI
            ani_values = []
            for contig_id in ordered_contig_ids:
                alignment = contig_alignments[contig_id]['alignment']
                diff = self._get_dv_from_tags(alignment.tags)
                if diff is not None:
                    ani = (1 - diff) * 100
                elif alignment.aligned_bases > 0:
                    ani = (alignment.nmatch / float(alignment.aligned_bases)) * 100
                else:
                    ani = 0
                ani_values.append(ani)

            avg_ani = sum(ani_values) / len(ani_values) if ani_values else 0.0

            actions_summary = f"ordered {num_aligned} contigs by reference position"
            if num_contigs_rotated > 0:
                actions_summary += f" (after rotating {P('circularly permuted contig', num_contigs_rotated)})"
            if num_contigs_split > 0:
                actions_summary += (f" (after cutting {P('contig', num_contigs_split)} that ran past the ends of the "
                                    f"reference into {P('fragment', num_fragments)})")
            if self.scaffold_fragmented:
                actions_summary += f", inserted {total_gap_size:,} nts of N-padding"
            if unaligned_contigs:
                actions_summary += f", appended {num_unaligned} unaligned contig(s) as-is"

            # Return results
            return {
                'num_contigs_processed': num_contigs_after_split,
                'num_contigs_in_input': num_contigs_after_filter,
                'num_contigs_rotated': num_contigs_rotated,
                'num_contigs_split': num_contigs_split,
                'num_fragments': num_fragments,
                'num_contigs_aligned': num_aligned,
                'num_contigs_unaligned': num_unaligned,
                'reference_coverage_pct': reference_coverage_pct,
                'gaps': gaps,
                'total_gap_size': total_gap_size,
                'avg_ani': avg_ani,
                'actions_summary': actions_summary,
            }

        finally:
            if not anvio.DEBUG:
                shutil.rmtree(temp_dir)


    def _align_contig_to_reference(self, contig, temp_dir, file_tag):
        """Aligns a single contig to the reference and returns its PAF records.

        Parameters
        ==========
        contig : dict
            A contig dictionary with 'id', 'seq', and 'length' keys.
        temp_dir : str
            Directory in which the single-contig FASTA file will be written.
        file_tag : str
            Basename (without extension) for that intermediate FASTA file.

        Returns
        =======
        paf_records : list
            PAF records for this contig against the reference, or None if the contig did not
            align to the reference at all.
        """
        temp_contig_fa = os.path.join(temp_dir, f"{file_tag}.fa")
        with open(temp_contig_fa, 'w') as f:
            f.write(f">{contig['id']}\n{contig['seq']}\n")

        try:
            return self._minimap2_align(self.reference_path, temp_contig_fa)
        except RuntimeError as e:
            self.log_run.info_single(f"'{contig['id']}': no alignment to the reference ({e})", level=2)
            return None


    def _find_circular_permutation(self, contig, paf_records, ref_length):
        """Works out whether a query contig is a circularly permuted version of a reference region.

           Every now and then a query contig aligns to the reference in two big pieces whose order
           is swapped: the first half of the contig matches a *later* part of the reference than
           its second half does. That looks like a large rearrangement, and if it were one, anvi'o
           would have no business touching it. But it usually is not one. It is an artifact of how
           the contig was pulled out of an assembly graph.

           Assemblers report contigs by walking an assembly graph, and that graph contains cycles --
           some of them biological (a whole circular replicon), most of them the result of repeats,
           coverage gaps, and the other ordinary complications of assembling a genome from short
           pieces. When the walk goes around a cycle, the assembler has to break it somewhere to
           report a linear contig, and where it breaks is arbitrary. The result is a contig that is
           perfectly co-linear with a better reference *once you rotate it*, but that looks
           scrambled in the middle as long as you read it from the arbitrary point the assembler
           happened to choose. The give-away is that the two ENDS of such a contig sit right next
           to each other on the reference: the contig closes on itself.

           That is the signature this function looks for, and it demands all of it:

            - the contig's alignment blocks, read along the contig, jump backwards in reference
              coordinates exactly once (a genuine rearrangement would rarely be this tidy, and more
              than one jump is not something a single rotation can explain), and

            - the reference position where the contig's first block begins is immediately after the
              position where its last block ends -- within `--min-contig-length` nucleotides. This
              is the part that matters: it is proof that the contig's two ends are neighbours, so
              rotating it invents no adjacency that is not already in the sequence. It is the same
              standard anvi'o applies to itself before rotating a reference.

           A contig that straddles the position at which the reference begins and ends fails the
           second test by a mile (its ends land at opposite ends of the reference), so it is left
           for `_find_reference_boundary_cut_points` to cut, which is the right treatment for it.

        Parameters
        ==========
        contig : dict
            A contig dictionary with 'id', 'seq', and 'length' keys.
        paf_records : list
            PAF records for this contig against the reference, or None if it did not align.
        ref_length : int
            Length of the reference sequence.

        Returns
        =======
        permutation : dict or None
            None when the contig is not a circular permutation. Otherwise a dictionary with
            'needs_rc' (whether the contig must be reverse-complemented first), 'position' (the
            offset at which to rotate it, in the frame that follows from 'needs_rc'), and
            'end_gap' (how many nucleotides of reference separate the contig's two ends, which is
            the evidence that they are neighbours).
        """
        primaries = [r for r in (paf_records or []) if r.is_primary]

        if len(primaries) < 2:
            return None

        contig_length = contig['length']

        best = max(primaries, key=lambda r: (r.aligned_bases, r.mapq, r.nmatch))
        needs_rc = best.strand == '-'
        oriented = lambda r: (contig_length - r.qend, contig_length - r.qstart) if needs_rc else (r.qstart, r.qend)

        blocks = sorted([oriented(r) + (r.tstart, r.tend) for r in primaries if r.strand == best.strand])

        if len(blocks) < 2:
            return None

        # exactly one place along the contig where the reference coordinate goes backwards
        backward_jumps = [i for i in range(len(blocks) - 1) if blocks[i + 1][2] < blocks[i][2]]

        if len(backward_jumps) != 1:
            return None

        # ... and the contig's two ends have to be neighbours on the reference, or this is a
        # rearrangement after all and rotating the contig would fabricate an adjacency
        end_gap = blocks[0][2] - blocks[-1][3]

        if abs(end_gap) > max(self.min_contig_length, 1):
            return None

        return {'needs_rc': needs_rc, 'position': blocks[backward_jumps[0] + 1][0], 'end_gap': end_gap}


    def _uncoil_circularly_permuted_contigs(self, contigs, alignments, ref_length, genome_name, temp_dir):
        """Rotates query contigs that an assembler cut out of a cycle in its assembly graph.

           See `_find_circular_permutation` for what these contigs are and how anvi'o recognizes
           them. This function rotates each one so that it becomes co-linear with the reference,
           re-aligns it, and reports what it did. The contig keeps its name and every one of its
           nucleotides: only the point at which it is considered to start moves, which is exactly
           the operation anvi'o performs on a circular reference.

           When the user asks for their contigs to be kept intact with `--keep-query-contigs-intact`
           nothing is rotated, since rotating restructures a contig even though it does not break
           it apart. The contigs are still reported, because a permuted contig that is left as it
           is will carry a large chunk of itself in the wrong place.

        Parameters
        ==========
        contigs : list
            List of contig dictionaries with 'id', 'seq', and 'length' keys. Updated in place for
            every contig that is rotated.
        alignments : dict
            Maps each contig id to its PAF records (or None). Entries of rotated contigs are
            replaced with their alignments after rotation.
        ref_length : int
            Length of the reference sequence.
        genome_name : str
            Name of the genome being processed, for reporting.
        temp_dir : str
            Directory for intermediate files.

        Returns
        =======
        num_contigs_rotated : int
            How many contigs were rotated.
        """
        report_lines, num_contigs_rotated = [], 0

        for index, contig in enumerate(contigs):
            permutation = self._find_circular_permutation(contig, alignments[contig['id']], ref_length)

            if not permutation:
                continue

            num_contigs_rotated += 1
            position, end_gap = permutation['position'], permutation['end_gap']
            neighbours = f"{abs(end_gap):,} nts {'apart' if end_gap >= 0 else 'of overlap'}"

            if self.keep_query_contigs_intact:
                report_lines.append((2, f"'{contig['id']}' ({contig['length']:,} nts), which would have become co-linear "
                                        f"with the reference by rotating it {position:,} nts (its two ends are "
                                        f"{neighbours} on the reference)."))
                continue

            sequence = utils.rev_comp(contig['seq']) if permutation['needs_rc'] else contig['seq']
            contig['seq'] = sequence[position:] + sequence[:position]

            alignments[contig['id']] = self._align_contig_to_reference(contig, temp_dir, f"rotated_{index}")

            blocks_before = len([r for r in alignments[contig['id']] or [] if r.is_primary])
            report_lines.append((2, f"'{contig['id']}' ({contig['length']:,} nts) was rotated by {position:,} nts"
                                    f"{' after being reverse-complemented' if permutation['needs_rc'] else ''}, since "
                                    f"its two ends are {neighbours} on the reference and it therefore closes on itself. "
                                    f"It now aligns to the reference in {P('piece', blocks_before)}."))

        if not num_contigs_rotated:
            return 0

        self.progress.clear()

        if self.keep_query_contigs_intact:
            self.run.warning(f"Anvi'o recognized {P('contig', num_contigs_rotated)} in '{genome_name}' that an assembler "
                             f"appears to have cut out of a cycle in its assembly graph. Read from the arbitrary point "
                             f"the assembler chose, such a contig looks scrambled against the reference, but its two "
                             f"ends sit right next to each other on it, which means it closes on itself and would line "
                             f"up perfectly if it were rotated. Anvi'o would normally do exactly that, but you asked "
                             f"for your contigs to be kept intact with `--keep-query-contigs-intact`, so they are left "
                             f"as they are -- which means a large part of each of them is going to sit in the wrong "
                             f"place in the output:",
                             header=f"CIRCULARLY PERMUTED CONTIGS IN {genome_name}")
        else:
            self.run.warning(f"Anvi'o rotated {P('contig', num_contigs_rotated)} in '{genome_name}'. Contigs like these "
                             f"come out of an assembler having been cut out of a cycle in its assembly graph, and the "
                             f"point at which they were cut is arbitrary, so they look scrambled against a better "
                             f"reference even though nothing is actually rearranged. Anvi'o can tell them apart from "
                             f"true rearrangements because their two ends are neighbours on the reference: they close "
                             f"on themselves, so rotating them invents no adjacency that was not already in the "
                             f"sequence. Every nucleotide and every contig name survives this; only the point at which "
                             f"each contig is considered to start has moved:",
                             header=f"CIRCULARLY PERMUTED CONTIGS ROTATED IN {genome_name}")

        for level, line in report_lines:
            self.run.info_single(line, level=level)

        self.run.info_single('', level=0, nl_after=1)

        return 0 if self.keep_query_contigs_intact else num_contigs_rotated


    def _find_reference_boundary_cut_points(self, contig, paf_records, ref_length):
        """Finds the positions at which a query contig must be cut to be placed on the reference.

           A query contig is a linear sequence, but the reference offers a single coordinate axis
           with a beginning and an end, and everything this program does -- ordering contigs,
           orienting them, scaffolding them -- is expressed in terms of that axis. A contig that
           runs past either end of the axis therefore has no single place to go: one part of it
           belongs at the far left of the axis, and another part at the far right. Writing such a
           contig out in one piece would quietly claim an agreement with the reference that the
           contig as a whole does not have, which is why anvi'o cuts it instead.

           Two situations produce this, and both are handled here:

            - The contig straddles the position at which the reference begins and ends. This is
              the common case for circular genomes: the reference was circularized (or rotated) at
              one position and the query genome at another, so a query contig that spans the
              reference's start/end boundary aligns to the very end of the reference with one part
              of itself and to the very beginning of it with another.

            - One end of the contig hangs off the beginning or the end of the reference: it has no
              alignment, AND there is not enough reference left before (or after) the contig's
              alignment for that sequence to fit. Such a stretch has nowhere to go on the axis, and
              it must not travel along inside a contig that is placed, where it would look like it
              had been placed too.

           The second test is the reason unaligned sequence alone is *not* enough to cut a contig.
           A contig that aligns comfortably inside the reference but carries an unaligned stretch
           at one of its ends is simply a genome with a variable region, and cutting it would be
           vandalism: there is plenty of room on the reference axis for that sequence to sit where
           it is. Only when the unaligned stretch is long compared to the reference that remains
           beyond the alignment does it actually run off the axis.

           Between them, these two tests are also why a query genome should not gain more than a
           couple of contigs. The reference axis has exactly one boundary, so at most one contig
           can contain it and be cut in two, and at most one contig on either side of it can hang
           off it. Everything else stays whole.

           Cutting is never done when it would produce a fragment shorter than
           `--min-contig-length`. This keeps the few hundred nucleotides `minimap2` routinely
           soft-clips off the ends of a divergent alignment from being promoted into contigs of
           their own, and it is the same threshold that decides whether a stretch that hangs off
           the reference is long enough to be worth worrying about.

           Rearrangements *within* a contig (inversions, internal translocations) are none of this
           function's business: they describe real differences between two genomes rather than an
           artifact of where a sequence happens to begin, and anvi'o keeps them intact.

        Parameters
        ==========
        contig : dict
            A contig dictionary with 'id', 'seq', and 'length' keys.
        paf_records : list
            PAF records for this contig against the reference, or None if it did not align.
        ref_length : int
            Length of the reference sequence.

        Returns
        =======
        cut_points : list
            Sorted list of (position, reason) tuples, where `position` is a 0-based offset into
            the contig *as it occurs in the input FASTA file* at which it must be cut, and
            `reason` is a short explanation of the cut for reporting purposes.
        """
        primaries = [r for r in (paf_records or []) if r.is_primary]

        # a contig that does not align to the reference at all is not misplaced by being kept in
        # one piece: all of it is equally unplaced, and it is reported as such
        if not primaries:
            return []

        contig_length = contig['length']
        min_fragment_length = max(self.min_contig_length, 1)

        # everything below happens in the 'oriented frame': the contig as it will be written out,
        # running in the same direction as the reference. a position `p` in that frame is the
        # position `contig_length - p` in the input frame when the contig needs reverse
        # complementing, which is how cut points are translated back at the very end
        best = max(primaries, key=lambda r: (r.aligned_bases, r.mapq, r.nmatch))
        needs_rc = best.strand == '-'
        oriented = lambda r: (contig_length - r.qend, contig_length - r.qstart) if needs_rc else (r.qstart, r.qend)

        cut_points = []

        blocks = sorted([oriented(r) + (r.tstart, r.tend) for r in primaries if r.strand == best.strand])
        aligned_stretches = self._merge_intervals([oriented(r) for r in primaries])

        # (1) an end of the contig that hangs off the reference. an end that merely has no
        # alignment is NOT enough: the reference that remains beyond the contig's own alignment has
        # to be too short to accommodate that sequence, since otherwise anvi'o would be cutting up
        # perfectly well-placed contigs simply for having a variable region at one of their ends.
        #
        # a stretch that is covered by an alignment on the other strand is left alone, too: that is
        # an inversion, and inversions are real biology rather than an accident of coordinates
        leading_flank_length, trailing_flank_length = blocks[0][0], contig_length - blocks[-1][1]
        reference_before_contig, reference_after_contig = blocks[0][2], ref_length - blocks[-1][3]

        if leading_flank_length or trailing_flank_length:
            self.log_run.info_single(f"'{contig['id']}': unaligned flanks -- {leading_flank_length:,} nts at the start "
                                     f"(with {reference_before_contig:,} nts of reference before the alignment to hold "
                                     f"it) and {trailing_flank_length:,} nts at the end (with {reference_after_contig:,} "
                                     f"nts of reference after it)", level=2)

        if aligned_stretches[0][0] == blocks[0][0] and leading_flank_length - reference_before_contig >= min_fragment_length:
            cut_points.append((blocks[0][0], "one of its ends hangs off the beginning of the reference"))

        if aligned_stretches[-1][1] == blocks[-1][1] and trailing_flank_length - reference_after_contig >= min_fragment_length:
            cut_points.append((blocks[-1][1], "one of its ends hangs off the end of the reference"))

        # (2) the contig straddles the position at which the reference begins and ends. this shows
        # up as two alignment blocks that follow one another along the contig, the first of which
        # reaches the end of the reference while the second one starts at its beginning. anvi'o
        # requires both blocks to actually reach their respective ends of the reference, so that a
        # contig that merely jumps backwards in reference coordinates -- which is a rearrangement,
        # and thus real biology -- is left alone
        boundary_window = max(min_fragment_length, ref_length // 100)

        for i in range(len(blocks) - 1):
            current_start, current_end, current_tstart, current_tend = blocks[i]
            next_start, next_end, next_tstart, next_tend = blocks[i + 1]

            reaches_end_of_reference = current_tend > ref_length - boundary_window
            starts_at_beginning_of_reference = next_tstart < boundary_window

            if not (reaches_end_of_reference and starts_at_beginning_of_reference and next_tstart < current_tstart):
                continue

            # project the position at which the reference begins and ends onto the contig from
            # each of the two blocks, and cut between them
            projected_from_current = current_end + (ref_length - current_tend)
            projected_from_next = next_start - next_tstart
            lower_bound, upper_bound = current_end, max(current_end, next_start)
            position = min(max((projected_from_current + projected_from_next) // 2, lower_bound), upper_bound)

            cut_points.append((position, "it straddles the position at which the reference begins and ends"))

        # back to the frame in which the contig occurs in the input FASTA file
        if needs_rc:
            cut_points = [(contig_length - position, reason) for position, reason in cut_points]

        # a cut that would leave behind a fragment shorter than the minimum contig length is not
        # worth making
        accepted_cut_points, previous_position = [], 0
        for position, reason in sorted(cut_points):
            if position - previous_position < min_fragment_length or contig_length - position < min_fragment_length:
                self.log_run.info_single(f"'{contig['id']}': not cutting at position {position:,} ({reason}) since it "
                                         f"would leave behind a fragment shorter than {min_fragment_length:,} nts", level=2)
                continue

            accepted_cut_points.append((position, reason))
            previous_position = position

        return accepted_cut_points


    def _split_contigs_at_reference_boundaries(self, contigs, alignments, ref_length, genome_name, temp_dir):
        """Cuts query contigs that run past either end of the reference coordinate axis.

           Anvi'o keeps query contigs intact whenever it can, but sometimes keeping one intact
           would amount to a false claim about how it agrees with the reference. See
           `_find_reference_boundary_cut_points` for what those cases are and why. This function
           replaces every such contig with the fragments it has to be cut into, aligns each
           fragment to the reference on its own merit so that it can be ordered and oriented like
           any other contig, and reports what it did.

           When the user asks for their contigs to be kept intact with
           `--keep-query-contigs-intact`, nothing is cut, but the contigs that *would* have been
           cut are still reported, since the user deserves to know which of the contigs in their
           output file are not really placed by this program.

        Parameters
        ==========
        contigs : list
            List of contig dictionaries with 'id', 'seq', and 'length' keys.
        alignments : dict
            Maps each contig id to its PAF records (or None). Updated in place: the entry of every
            contig that is cut is removed, and one entry per resulting fragment is added.
        ref_length : int
            Length of the reference sequence.
        genome_name : str
            Name of the genome being processed, for reporting.
        temp_dir : str
            Directory for intermediate files.

        Returns
        =======
        contigs : list
            The new list of contigs, in which every contig that had to be cut is replaced by the
            fragments it was cut into.
        num_contigs_split : int
            How many contigs were cut.
        num_fragments : int
            How many fragments those contigs were cut into.
        """
        resulting_contigs, report_lines = [], []
        num_contigs_split, num_fragments = 0, 0

        for contig in contigs:
            cut_points = self._find_reference_boundary_cut_points(contig, alignments[contig['id']], ref_length)

            if not cut_points:
                resulting_contigs.append(contig)
                continue

            num_contigs_split += 1
            reasons = ' and '.join(sorted(set(reason for position, reason in cut_points)))
            positions = ', '.join([f"{position:,}" for position, reason in cut_points])

            if self.keep_query_contigs_intact:
                resulting_contigs.append(contig)
                report_lines.append((2, f"'{contig['id']}' ({contig['length']:,} nts), since {reasons}."))
                continue

            # fragments are numbered along the original contig, from its first nucleotide to its
            # last, so that their names describe how they were related to one another in the input
            # assembly no matter where each of them ends up in the output file
            boundaries = [0] + [position for position, reason in cut_points] + [contig['length']]
            fragments = [{'id': f"{contig['id']}_fragment_{i + 1:02d}",
                          'seq': contig['seq'][boundaries[i]:boundaries[i + 1]],
                          'length': boundaries[i + 1] - boundaries[i]} for i in range(len(boundaries) - 1)]

            del alignments[contig['id']]
            num_fragments += len(fragments)

            report_lines.append((2, f"'{contig['id']}' ({contig['length']:,} nts) was cut into {P('fragment', len(fragments))} at "
                                    f"{P('position', len(cut_points), alt='positions')} {positions}, since {reasons}:"))

            for fragment in fragments:
                alignments[fragment['id']] = self._align_contig_to_reference(fragment, temp_dir, f"fragment_{num_fragments}_{len(resulting_contigs)}")
                resulting_contigs.append(fragment)

                primaries = [r for r in (alignments[fragment['id']] or []) if r.is_primary]
                if primaries:
                    best = max(primaries, key=lambda r: (r.aligned_bases, r.mapq, r.nmatch))
                    report_lines.append((3, f"'{fragment['id']}' ({fragment['length']:,} nts) aligns to the reference at "
                                            f"position {best.tstart:,} on the {best.strand} strand."))
                else:
                    report_lines.append((3, f"'{fragment['id']}' ({fragment['length']:,} nts) does not align to the reference "
                                            f"anywhere, and is reported as an unaligned contig."))

        if not num_contigs_split:
            return resulting_contigs, 0, 0

        self.progress.clear()

        if self.keep_query_contigs_intact:
            self.run.warning(f"Anvi'o found {P('contig', num_contigs_split)} in '{genome_name}' running past the ends of "
                             f"the reference genome, which means there is no single position on the reference where such "
                             f"a contig can be placed in one piece. Normally anvi'o would cut such contigs into fragments "
                             f"and place each fragment independently, but since you asked for your contigs to be kept intact "
                             f"with `--keep-query-contigs-intact`, they will be written out whole. Please keep in mind "
                             f"that parts of these contigs are therefore NOT aligned to the reference in the output file, "
                             f"even though the file will look like they are:",
                             header=f"CONTIGS THAT RUN PAST THE REFERENCE IN {genome_name}")
        else:
            self.run.warning(f"Anvi'o had to cut {P('contig', num_contigs_split)} in '{genome_name}' into "
                             f"{P('fragment', num_fragments)}, which means this genome will have more contigs in the "
                             f"output file than it had in the input file. This is not anvi'o being clumsy: a contig that "
                             f"runs past the beginning or the end of the reference has no single position on the "
                             f"reference to be placed at, since one part of it belongs to the very start of the reference "
                             f"coordinates and another part of it to the very end. Keeping such a contig in one piece "
                             f"would have meant reporting it as aligned when a part of it is not. Each fragment below was "
                             f"aligned to the reference on its own merit, and ordered and oriented accordingly. If you "
                             f"would rather have your contigs kept intact and are willing to accept that some of them "
                             f"will not really be aligned to the reference, please use the flag "
                             f"`--keep-query-contigs-intact`:",
                             header=f"QUERY CONTIGS CUT INTO FRAGMENTS IN {genome_name}")

        for level, line in report_lines:
            self.run.info_single(line, level=level)

        self.run.info_single('', level=0, nl_after=1)

        if self.keep_query_contigs_intact:
            return resulting_contigs, 0, 0

        return resulting_contigs, num_contigs_split, num_fragments


    def _minimap2_align(self, ref_fa, qry_fa, find_all_alignments=False):
        cmd = ["minimap2",
               "-x",
               self.minimap2_preset,
               "-t",
               str(self.threads)]

        # When finding the optimal reference start, we want ALL possible alignments
        # to find truly conserved regions, not just the best alignment
        if find_all_alignments:
            cmd.extend([
                "--secondary=yes",  # Output secondary alignments
                "-N", "100",        # Output up to 100 secondary alignments
                "-p", "0.5"         # Report secondary alignments with score >= 50% of primary
            ])

        cmd.extend([ref_fa, qry_fa])

        self.log_run.info_single(f"Running minimap2: {' '.join(cmd)}", level=2)
        stdout = utils.run_command_and_get_output(cmd, raise_on_error=True)
        return self._read_paf(stdout.splitlines())


    def _seqkit_reverse_complement(self, in_fa, out_fa):
        cmd = ["seqkit", "seq", "-r", "-p", in_fa]
        self.log_run.info_single(f"Running seqkit reverse-complement: {' '.join(cmd)} > {out_fa}", level=2)
        stdout = utils.run_command_and_get_output(cmd, raise_on_error=True)
        with open(out_fa, "w", encoding="utf-8") as f:
            f.write(stdout)


    def _seqkit_rotate(self, in_fa, cut1, out_fa):
        if cut1 < 1:
            raise ConfigError(f"seqkit restart expects 1-based position >=1, got {cut1}")
        cmd = ["seqkit", "restart", "-i", str(cut1), in_fa]
        self.log_run.info_single(f"Running seqkit restart: {' '.join(cmd)} > {out_fa}", level=2)
        stdout = utils.run_command_and_get_output(cmd, raise_on_error=True)
        with open(out_fa, "w", encoding="utf-8") as f:
            f.write(stdout)


    def _read_paf(self, paf_lines):
        recs = []
        for line in paf_lines:
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 12:
                continue
            qname = parts[0]
            qlen = int(parts[1])
            qstart = int(parts[2])
            qend = int(parts[3])
            strand = parts[4]
            tname = parts[5]
            tlen = int(parts[6])
            tstart = int(parts[7])
            tend = int(parts[8])
            nmatch = int(parts[9])
            alen = int(parts[10])
            mapq = int(parts[11])
            tags = "\t".join(parts[12:]) if len(parts) > 12 else ""
            recs.append(
                PafRecord(
                    qname=qname,
                    qlen=qlen,
                    qstart=qstart,
                    qend=qend,
                    strand=strand,
                    tname=tname,
                    tlen=tlen,
                    tstart=tstart,
                    tend=tend,
                    nmatch=nmatch,
                    alen=alen,
                    mapq=mapq,
                    tags=tags,
                )
            )

        if not recs:
            raise RuntimeError("PAF output is empty; minimap2 did not produce any alignments.")

        primaries = [r for r in recs if r.is_primary]
        self.log_run.info_single(f"PAF records: total={len(recs)} primary={len(primaries)}", level=2)

        return recs


    def _select_anchor_near_reference_start(self, recs, near_bp):
        primaries = [r for r in recs if r.is_primary]
        if not primaries:
            raise RuntimeError("No primary alignments (tp:A:P) found in PAF.")

        near = [r for r in primaries if r.tstart < near_bp]
        if near:
            chosen = max(near, key=lambda r: (r.aligned_bases, r.mapq, r.nmatch))
            self.log_run.info_single(
                f"Pass-1: selected anchor near reference start (tstart={chosen.tstart})", level=2)
            return chosen

        chosen = max(primaries, key=lambda r: (r.aligned_bases, r.mapq, r.nmatch))
        self.log_run.info_single("Pass-1: no near-start hit; using best primary.", level=2)
        return chosen


    def _select_anchor_for_ref0(self, recs):
        primaries = [r for r in recs if r.is_primary]
        if not primaries:
            raise RuntimeError("No primary alignments (tp:A:P) found in PAF.")

        def dist_to_zero(r):
            return min(r.tstart, r.tlen - r.tstart)

        return sorted(primaries, key=lambda r: (dist_to_zero(r), -r.aligned_bases, -r.mapq))[0]


    def _select_anchor_smallest_reference_start(self, recs):
        primaries = [r for r in recs if r.is_primary]
        if not primaries:
            raise RuntimeError("No primary alignments (tp:A:P) found in PAF.")
        return sorted(primaries, key=lambda r: (r.tstart, -r.aligned_bases, -r.mapq))[0]


    def _cut0_for_ref0(self, rec):
        """Calculate where to cut the query so that ref[0] aligns with query[0].

        Given an alignment where query[qstart..qend] aligns with ref[tstart..tend],
        we want to find the position in the query that corresponds to ref[0],
        then rotate the query so that position becomes query[0].
        """
        if rec.strand not in {"+", "-"}:
            raise ConfigError(f"Unexpected strand in PAF: {rec.strand}")

        if rec.strand == "+":
            # For + strand: query[qstart] aligns with ref[tstart]
            # So query[qstart - tstart] aligns with ref[0] (may be negative, handled by modulo)
            # We want to rotate query so that this position becomes query[0]
            # seqkit restart cuts at the position and makes it the new start
            # So cut position is: qstart - tstart
            # But when both qstart and tstart have the same offset from 0,
            # we need to rotate by that offset, which is qstart (or tstart)
            # Actually: position in query that aligns with ref[0] is: qstart - tstart
            # That position should become query[0], so we rotate by: qstart - tstart
            cut = rec.qstart - rec.tstart
        else:
            # For - strand alignment (before RC):
            # PAF coordinates are on forward strand. The RC of query[qstart:qend] aligns with ref[tstart:tend]
            # We need to find position in original (pre-RC) query that, after RC of whole sequence, aligns with ref[0]
            # After whole-sequence RC: original position (qend-1) aligns with ref[tstart]
            # Going back by tstart in ref means going forward by tstart in original query before RC
            # Position in original query: (qend - 1) + tstart
            cut = (rec.qend - 1) + rec.tstart

        # Handle negative cuts and wrap around circular genome
        cut %= rec.qlen
        return cut


    def _best_primary_alignment(self, recs):
        primaries = [r for r in recs if r.is_primary]
        if not primaries:
            raise RuntimeError("Final validation: no primary alignments found.")
        best = max(primaries, key=lambda r: (r.aligned_bases, r.mapq, r.nmatch))
        self.log_run.info_single(
            f"Final primary: strand={best.strand} qstart={best.qstart} tstart={best.tstart} "
            f"alen={best.aligned_bases} mapq={best.mapq}",
            level=2)
        return best


    def _final_report(self, results):
        """Show some final information to the user"""

        trust_counts = {"TRUSTWORTHY": [], "SOMEWHAT OK": [], "NOT TRUSTWORTHY": [], "FAILED": [], "REFERENCE": []}
        for r in results:
            label = r.trust or r.status
            if label in trust_counts:
                trust_counts[label].append((r.name, r.output_path))
            else:
                trust_counts["FAILED"].append((r.name, r.output_path))

        num_reference = len(trust_counts["REFERENCE"])
        num_reoriented = len(results) - num_reference
        num_somewhat = len(trust_counts["SOMEWHAT OK"])
        num_not_trustworthy = len(trust_counts["NOT TRUSTWORTHY"])
        num_failed = len(trust_counts["FAILED"])

        if num_not_trustworthy or num_failed:
            review_msg = "the alignment plots above and " if not self.skip_visualizing_alignments else ""
            final_text = (f"Your genome reorientation task considered {num_reoriented} genomes and {num_reference} "
                          f"reference. Some outcomes are not trustworthy (low alignment coverage can lead to unreliable "
                          f"orientation). Please review {review_msg}the summary below to decide which FASTA "
                          f"files to use downstream.")
        elif num_somewhat:
            review_msg = "the alignment plots and " if not self.skip_visualizing_alignments else ""
            final_text = (f"Your genome reorientation task considered {num_reoriented} genomes and {num_reference} "
                          f"reference. Some outcomes are only somewhat OK; check {review_msg}stats below before "
                          f"using them downstream.")
        else:
            review_msg = "A quick glance at the alignment plots above is still recommended before downstream analyses." if not self.skip_visualizing_alignments else ""
            final_text = (f"Your genome reorientation task for {num_reoriented} genomes and {num_reference} reference "
                          f"is complete. All outcomes look trustworthy based on alignment coverage. {review_msg}")

        self.run.warning(final_text, header="FINAL REPORT")

        trustworthy_list = trust_counts["REFERENCE"] + trust_counts["TRUSTWORTHY"]
        self.run.info("Trustworthy", len(trustworthy_list))
        for name, path in trustworthy_list:
            color = "cyan" if (name, path) in trust_counts["REFERENCE"] else "green"
            self.run.info_single(f"{name} -> {path}", mc=color, level=2, cut_after=0)

        self.run.info("Somewhat OK", len(trust_counts["SOMEWHAT OK"]), nl_before=1)
        for name, path in trust_counts["SOMEWHAT OK"]:
            self.run.info_single(f"{name} -> {path}", mc="yellow", level=2, cut_after=0)

        self.run.info("Not trustworthy", len(trust_counts["NOT TRUSTWORTHY"]), nl_before=1)
        for name, path in trust_counts["NOT TRUSTWORTHY"]:
            self.run.info_single(f"{name} -> {path}", mc="red", level=2, cut_after=0)

        self.run.info("Failed", len(trust_counts["FAILED"]), nl_before=1)
        for name, path in trust_counts["FAILED"]:
            self.run.info_single(f"{name} -> {path}", mc="red", level=2, cut_after=0)


    def _get_dv_from_tags(self, tags):
        if not tags:
            return None
        for t in tags.split('\t'):
            if t.startswith('dv:f:'):
                try:
                    return float(t.split(':')[-1])
                except ValueError:
                    return None
        return None


    def _plot_synteny_ribbons(self, recs, genome_name, label):
        """Plot synteny-style ribbons showing alignment blocks between reference and query."""
        try:
            import plotext as plt
        except ImportError:
            self.run.warning(f"plotext is not available; skipping alignment plot for '{genome_name}'.")
            return

        primaries = [r for r in recs if r.is_primary]
        if not primaries:
            self.run.warning(f"No primary alignments to plot for '{genome_name}'.")
            return

        # Sort by alignment size and take top alignments for clarity
        primaries = sorted(primaries, key=lambda r: r.aligned_bases, reverse=True)[:30]

        ref_len = primaries[0].tlen

        # Check if we have multiple contigs (different qnames)
        unique_qnames = list(set([r.qname for r in primaries]))
        has_multiple_contigs = len(unique_qnames) > 1

        # Calculate offsets for multi-contig genomes to avoid overlap
        contig_offsets = {}
        total_query_span = 0

        if has_multiple_contigs:
            # Group records by contig and calculate offsets
            contig_groups = {}
            for r in primaries:
                if r.qname not in contig_groups:
                    contig_groups[r.qname] = []
                contig_groups[r.qname].append(r)

            # Sort contigs by their position on reference (for better visual flow)
            sorted_contigs = sorted(contig_groups.keys(),
                                   key=lambda qn: min(r.tstart for r in contig_groups[qn]))

            current_offset = 0
            gap_size = 50000  # 50kb gap between contigs for visualization

            for qname in sorted_contigs:
                contig_offsets[qname] = current_offset
                contig_len = contig_groups[qname][0].qlen
                current_offset += contig_len + gap_size

            total_query_span = current_offset - gap_size  # Remove last gap
        else:
            # Single contig case
            qry_len = primaries[0].qlen
            total_query_span = qry_len
            contig_offsets[primaries[0].qname] = 0

        # Get plot dimensions - use user-specified or defaults
        plot_width = self.plot_width if self.plot_width else terminal.get_terminal_width()
        plot_height = self.plot_height

        plt.clf()
        plt.plotsize(plot_width, plot_height)
        plt.theme("dark")

        # Define y-coordinates for reference and query tracks
        ref_y = 2.0
        qry_y = 0.0

        # Draw alignment ribbons connecting reference to query first
        for r in primaries:
            ref_start = r.tstart
            ref_end = r.tend

            # Apply offset for this contig
            offset = contig_offsets[r.qname]

            if r.strand == '+':
                # Forward strand: connect in same direction
                qry_start = r.qstart + offset
                qry_end = r.qend + offset
                color = 'green'
            else:
                # Reverse strand: connect in opposite direction (shows inversion)
                qry_start = r.qend + offset
                qry_end = r.qstart + offset
                color = 'red'

            # Draw left edge of ribbon (alignment start) - use braille for vertical lines
            plt.plot([ref_start, qry_start], [ref_y, qry_y], color=color, marker='braille')

            # Draw right edge of ribbon (alignment end) - use braille for vertical lines
            plt.plot([ref_end, qry_end], [ref_y, qry_y], color=color, marker='braille')

            # Draw top and bottom edges of the alignment block - keep solid for horizontal bars
            plt.plot([ref_start, ref_end], [ref_y, ref_y], color=color)
            plt.plot([qry_start, qry_end], [qry_y, qry_y], color=color)

        # Draw genome lines AFTER alignment ribbons so they appear continuous on top
        # Draw reference genome line (top) - always continuous
        plt.plot([0, ref_len], [ref_y, ref_y], color='white')

        # Draw query genome line (bottom) - color-coded by contig type
        if has_multiple_contigs:
            # Draw individual segments for each contig with color coding
            contig_groups = {}
            contig_strands = {}
            for r in primaries:
                if r.qname not in contig_groups:
                    contig_groups[r.qname] = r.qlen
                    contig_strands[r.qname] = r.strand

            for qname, qlen in contig_groups.items():
                offset = contig_offsets[qname]

                # Determine color based on contig type
                if '_wrapPart' in qname:
                    # Wrapped-around contigs are green
                    segment_color = 'green'
                elif contig_strands[qname] == '-':
                    # Reverse-complemented contigs are red
                    segment_color = 'red'
                else:
                    # Forward strand, non-wrapped contigs are white
                    segment_color = 'white'

                plt.plot([offset, offset + qlen], [qry_y, qry_y], color=segment_color)
        else:
            # Single contig - determine color
            strand = primaries[0].strand
            segment_color = 'red' if strand == '-' else 'white'
            plt.plot([0, total_query_span], [qry_y, qry_y], color=segment_color)

        # Format the plot
        contig_info = f" ({len(unique_qnames)} contigs)" if has_multiple_contigs else ""
        plt.title(f"{label}: {genome_name}{contig_info} ↔ {self.reference_name} (green=forward, red=reverse)")
        plt.xlabel("Genome position (nts)")
        plt.ylim(-0.1, 2.1)

        # Set y-axis labels
        plt.yticks([0, 2], [f"{genome_name} [Q]", f"{self.reference_name} [R]"])

        # Set x-axis with human-readable numbers
        try:
            max_len = max(ref_len, total_query_span)
            x_ticks = [0, max_len * 0.25, max_len * 0.5, max_len * 0.75, max_len]
            plt.xticks(x_ticks, [utils.human_readable_number(x, decimals=1) for x in x_ticks])
        except Exception:
            pass

        plt.show()


    def _find_optimal_reference_start(self):
        """Survey all genomes to find the optimal starting position in the reference.

        Returns (optimal_position, num_genomes_covering, total_genomes) or None if reference is user-specified.
        """
        # Get reference length
        ref_len = self._get_total_length(self.reference_path)

        # Use binning to reduce memory and computational overhead
        # Since we use a 1000bp window anyway, we don't need single-base resolution
        bin_size = 1000  # 1kb bins
        num_bins = (ref_len + bin_size - 1) // bin_size  # Ceiling division

        # Initialize coverage bins: for each bin, track SET of genomes that cover it
        coverage_bins = [set() for _ in range(num_bins)]

        self.run.warning(None, header="FINDING OPTIMAL REFERENCE START")
        self.run.info("Reference genome", self.reference_name)
        self.run.info("Reference length", f"{ref_len:,} nts")
        self.run.info("Bin size for coverage", f"{bin_size:,} nts")
        self.run.info("Number of bins", f"{num_bins:,}")

        total_genomes = len([g for g in self.genomes.keys() if g != self.reference_name])
        if total_genomes == 0:
            self.run.warning("No genomes to survey (only reference present)")
            return None

        self.progress.new("Surveying genomes", progress_total_items=total_genomes)
        genome_counter = 0

        for genome_name, entry in self.genomes.items():
            if genome_name == self.reference_name:
                continue

            genome_counter += 1
            query_path = entry['path']

            self.progress.update(f"Aligning {genome_name} ({genome_counter}/{total_genomes})")
            self.log_run.info_single(f"Aligning {genome_name} to find coverage (including secondary alignments)", level=2)
            paf_recs = self._minimap2_align(self.reference_path, query_path, find_all_alignments=True)

            # Mark bins covered by ALL alignments (primary + secondary)
            # This gives us a true picture of all conserved regions, not just the "best" alignment
            self.log_run.info_single(f"Found {len(paf_recs)} alignments (primary + secondary) for {genome_name}", level=2)

            for rec in paf_recs:
                # Calculate which bins this alignment overlaps
                start_bin = rec.tstart // bin_size
                end_bin = min(rec.tend // bin_size, num_bins - 1)

                # Mark all bins in this range as covered by this genome
                for bin_idx in range(start_bin, end_bin + 1):
                    coverage_bins[bin_idx].add(genome_name)

            self.progress.increment()

        self.progress.end()

        # Find the longest contiguous high-coverage region
        # This is better than greedy first-match as it avoids edges of conserved regions
        self.run.info_single("Finding longest conserved region instead of first match", level=2)

        # Determine coverage threshold (at least 80% of genomes, or all if there are few)
        min_coverage_threshold = max(int(total_genomes * 0.8), total_genomes - 1) if total_genomes > 3 else total_genomes

        # Find all contiguous regions with high coverage
        conserved_regions = []
        current_region_start = None
        current_region_end = None

        for bin_idx in range(num_bins):
            bin_coverage = len(coverage_bins[bin_idx])

            if bin_coverage >= min_coverage_threshold:
                if current_region_start is None:
                    # Start a new conserved region
                    current_region_start = bin_idx
                    current_region_end = bin_idx
                else:
                    # Extend current region
                    current_region_end = bin_idx
            else:
                if current_region_start is not None:
                    # End current region and save it
                    region_length = (current_region_end - current_region_start + 1) * bin_size
                    avg_coverage = sum(len(coverage_bins[i]) for i in range(current_region_start, current_region_end + 1)) / (current_region_end - current_region_start + 1)
                    conserved_regions.append({
                        'start_bin': current_region_start,
                        'end_bin': current_region_end,
                        'start_bp': current_region_start * bin_size,
                        'end_bp': min((current_region_end + 1) * bin_size, ref_len),
                        'length': region_length,
                        'avg_coverage': avg_coverage
                    })
                    current_region_start = None
                    current_region_end = None

        # Handle wrap-around for circular genomes
        if current_region_start is not None:
            # Region extends to end of genome, check if it wraps around
            region_length = (current_region_end - current_region_start + 1) * bin_size
            avg_coverage = sum(len(coverage_bins[i]) for i in range(current_region_start, current_region_end + 1)) / (current_region_end - current_region_start + 1)
            conserved_regions.append({
                'start_bin': current_region_start,
                'end_bin': current_region_end,
                'start_bp': current_region_start * bin_size,
                'end_bp': min((current_region_end + 1) * bin_size, ref_len),
                'length': region_length,
                'avg_coverage': avg_coverage
            })

        if not conserved_regions:
            self.run.warning("No highly conserved regions found. Using position 0 as fallback.")
            return 0, 0, total_genomes

        # Find the longest conserved region
        longest_region = max(conserved_regions, key=lambda r: r['length'])

        self.run.info("Conserved regions found", len(conserved_regions))
        self.run.info("Longest conserved region", f"{longest_region['start_bp']:,} - {longest_region['end_bp']:,} bp ({longest_region['length']:,} nts)")
        self.run.info("Average coverage in longest region", f"{longest_region['avg_coverage']:.1f} genomes")

        # Now find the optimal position within this region that doesn't split a gene
        best_position = self._find_position_between_genes(longest_region, ref_len)

        # Get coverage at the selected position
        selected_bin = best_position // bin_size
        genomes_covering = len(coverage_bins[selected_bin])
        coverage_pct = (genomes_covering / total_genomes) * 100 if total_genomes > 0 else 0

        self.run.info("Optimal start position", f"{best_position:,} nts (middle of longest conserved region, between genes)")
        self.run.info("Genomes covering this position", f"{genomes_covering}/{total_genomes} ({coverage_pct:.1f}%)")

        return best_position, genomes_covering, total_genomes


    def _find_position_between_genes(self, conserved_region, ref_len):
        """Find the best position within a conserved region that doesn't split a gene.

        Args:
            conserved_region: dict with 'start_bp' and 'end_bp' keys
            ref_len: total reference length

        Returns:
            Position (bp) in the middle of the conserved region, ideally between genes
        """
        region_start = conserved_region['start_bp']
        region_end = conserved_region['end_bp']
        region_middle = (region_start + region_end) // 2

        self.run.info_single("Calling genes to find intergenic position in conserved region", level=2)

        # Call genes on reference
        temp_dir = filesnpaths.get_temp_directory_path()
        try:
            gene_caller = GeneCaller(self.reference_path, gene_caller='pyrodigal-gv', run=terminal.Run(verbose=False), progress=self.progress)
            gene_calls_dict, _ = gene_caller.process()

            if not gene_calls_dict:
                self.run.warning("No genes found in reference. Using middle of conserved region.")
                return region_middle

            self.run.info_single(f"Found {len(gene_calls_dict)} genes in reference", level=2)

            # Extract gene boundaries
            gene_boundaries = []
            for gene_id, gene_info in gene_calls_dict.items():
                gene_boundaries.append({
                    'id': gene_id,
                    'start': gene_info['start'],
                    'stop': gene_info['stop']
                })

            # Sort by start position
            gene_boundaries.sort(key=lambda g: g['start'])

            # Find intergenic regions within the conserved region
            intergenic_regions = []

            for i in range(len(gene_boundaries)):
                gene1 = gene_boundaries[i]
                gene2 = gene_boundaries[(i + 1) % len(gene_boundaries)]  # Next gene (wrap around)

                # Calculate intergenic region
                if i < len(gene_boundaries) - 1:
                    # Normal case: gap between consecutive genes
                    intergenic_start = gene1['stop']
                    intergenic_end = gene2['start']
                else:
                    # Wrap-around case: from last gene to first gene
                    intergenic_start = gene1['stop']
                    intergenic_end = gene2['start'] + ref_len  # Adjust for circular

                # Check if this intergenic region overlaps with our conserved region
                # For simplicity, use the middle point of the intergenic region
                if intergenic_start < intergenic_end:
                    intergenic_middle = (intergenic_start + intergenic_end) // 2
                    if region_start <= intergenic_middle <= region_end:
                        distance_from_region_middle = abs(intergenic_middle - region_middle)
                        intergenic_regions.append({
                            'position': intergenic_middle,
                            'distance_from_middle': distance_from_region_middle,
                            'length': intergenic_end - intergenic_start
                        })

            if intergenic_regions:
                # Choose intergenic region closest to the middle of conserved region
                best_intergenic = min(intergenic_regions, key=lambda r: r['distance_from_middle'])
                self.run.info_single(f"Selected intergenic position at {best_intergenic['position']:,} nts "
                                     f"({best_intergenic['distance_from_middle']:,} nts from literal region middle, "
                                     f"where the intergenic region is {best_intergenic['length']:,} nts)", level=2, cut_after=None)
                return best_intergenic['position']
            else:
                # No intergenic region found in conserved region, use middle
                self.run.info_single(
                    f"No intergenic regions in conserved region. Using region middle: {region_middle:,} nts",
                    level=2
                )
                return region_middle

        except Exception as e:
            self.run.warning(f"Could not call genes for optimal positioning: {e}. Using middle of conserved region.")
            return region_middle
        finally:
            if not anvio.DEBUG:
                shutil.rmtree(temp_dir)


    def _find_dnaa_gene_position(self):
        """Find the DnaA gene in the reference genome using HMM search.

        Returns the start position of the DnaA gene, or None if not found.
        """
        self.run.warning(None, header="FINDING DnaA GENE IN REFERENCE")
        self.run.info("Reference genome", self.reference_name)

        # Create temporary directory for gene calling and HMM search
        temp_dir = filesnpaths.get_temp_directory_path()

        try:
            # Step 1: Call genes using pyrodigal-gv
            self.run.info("Step 1", "Calling genes")

            gene_caller = GeneCaller(self.reference_path, gene_caller='pyrodigal-gv', run=self.log_run, progress=self.progress)
            gene_calls_dict, aa_sequences_dict = gene_caller.process()

            if not aa_sequences_dict:
                self.run.warning("Gene caller did not find any genes in the reference genome. DnaA orientation will be skipped.")
                return None

            self.run.info("Genes called", f"{len(gene_calls_dict)}")

            # Step 2: Write amino acid sequences to FASTA
            aa_fasta_path = os.path.join(temp_dir, "genes.faa")
            with open(aa_fasta_path, 'w') as f:
                for gene_id, seq in aa_sequences_dict.items():
                    f.write(f">{gene_id}\n{seq}\n")

            # Step 3: Prepare HMM files
            self.run.info("Step 2", "Searching for DnaA gene using HMM")

            # Path to DnaA HMM directory
            dnaa_hmm_dir = os.path.join(os.path.dirname(anvio.__file__),
                                         'data', 'misc', 'HMM', 'ANCHORS', 'Bac_DnaA_C')

            if not os.path.exists(dnaa_hmm_dir):
                raise ConfigError(f"DnaA HMM directory not found at {dnaa_hmm_dir}")

            # Decompress the HMM file if needed
            hmm_gz_path = os.path.join(dnaa_hmm_dir, 'genes.hmm.gz')
            hmm_path = os.path.join(temp_dir, 'Bac_DnaA_C.hmm')

            with gzip.open(hmm_gz_path, 'rb') as f_in:
                with open(hmm_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

            # Press the HMM file
            utils.run_command(['hmmpress', hmm_path], log_file_path=os.path.join(temp_dir, 'hmmpress.log'))

            # Step 4: Run hmmsearch
            hmmer_dir = temp_dir
            target_files = {'AA:GENE': aa_fasta_path}

            hmmer = HMMer(target_files, num_threads_to_use=1, program_to_use='hmmsearch',
                         run=self.log_run, progress=self.progress)

            table_output = hmmer.run_hmmer(
                source='DnaA',
                alphabet='AA',
                context='GENE',
                kind='Bac_DnaA_C',
                domain=None,
                num_genes_in_model=1,
                hmm=hmm_path,
                ref='Pfam',
                noise_cutoff_terms='--cut_ga',
                desired_output='table',
                hmmer_output_dir=hmmer_dir
            )

            # Step 5: Parse results
            if not table_output or not os.path.exists(table_output):
                self.run.warning("No DnaA gene found in the reference genome. DnaA orientation will be skipped.")
                return None

            # Parse the table output
            hits = []
            with open(table_output, 'r') as f:
                for line in f:
                    if line.startswith('#') or not line.strip():
                        continue
                    fields = line.strip().split()
                    if len(fields) >= 19:
                        gene_callers_id = int(fields[0])
                        bit_score = float(fields[5])
                        e_value = float(fields[4])
                        hits.append({
                            'gene_callers_id': gene_callers_id,
                            'bit_score': bit_score,
                            'e_value': e_value
                        })

            if not hits:
                self.run.warning("No DnaA gene found in the reference genome. DnaA orientation will be skipped.")
                return None

            # Get the best hit (highest bit score)
            best_hit = max(hits, key=lambda x: x['bit_score'])
            gene_id = best_hit['gene_callers_id']

            # Get the gene position from the gene calls
            if gene_id not in gene_calls_dict:
                self.run.warning(f"Gene {gene_id} not found in gene calls. DnaA orientation will be skipped.")
                return None

            gene_info = gene_calls_dict[gene_id]
            dnaa_start = gene_info['start']
            dnaa_stop = gene_info['stop']

            self.run.info("DnaA gene found", f"Gene {gene_id}")
            self.run.info("DnaA location", f"{dnaa_start:,} - {dnaa_stop:,} nts")
            self.run.info("DnaA bit score", f"{best_hit['bit_score']:.1f}")
            self.run.info("DnaA e-value", f"{best_hit['e_value']:.2e}", nl_after=1)

            return dnaa_start

        finally:
            # Clean up temporary directory
            if not anvio.DEBUG:
                shutil.rmtree(temp_dir)
            else:
                self.run.warning(f"Temp directory kept for debugging: {temp_dir}")


    def select_reference_genome(self):
        """Determine the reference genome whether it is from the user or de novo"""

        user_specified = self.reference_name is not None

        if user_specified:
            self.run.info("Reference genome", f"'{self.reference_name}' (specified by the user)")
            self.reference_path = self.genomes[self.reference_name]['path']
        else:
            candidate_stats = []
            for genome_name, entry in self.genomes.items():
                # `num_contigs` is already in `self.genomes` thanks to `sanity_check`, so there is no
                # need to go through every FASTA file one more time here just to count sequences.
                candidate_stats.append((entry['num_contigs'], -self._get_total_length(entry['path']), genome_name, entry['path']))

            candidate_stats.sort()
            chosen = candidate_stats[0]
            self.reference_name = chosen[2]
            self.reference_path = chosen[3]
            chosen_contigs = chosen[0]
            chosen_len = -chosen[1]
            self.run.info("Reference genome", f"'{self.reference_name}' (selected by anvi'o)")
            self.run.info("Reference contigs", chosen_contigs)
            self.run.info("Reference length", chosen_len)

        # everything this program does to the reference that involves 'rotations' assumes
        # that the reference is a complete circular sequence. we can't really be 100% sure
        # it the sequence is circular, but we can make sure at least some part of it (i.e.,
        # being a single contig actually holds:
        self.check_reference_is_single_contig(user_specified)


    def check_reference_is_single_contig(self, user_specified):
        """Make sure the reference genome is a single contig, and complain loudly if it is not.

           This program needs the reference to be a single contig, since everything it does for the
           user, such as orienting contigs, ordering them, scaffolding them, requires a single continuous
           coordinate system to work against. Coordinates that come from different contigs of a
           fragmented reference simply do not form one axis.

           Please note that this is a separate concern from *circularity*. A reference must also be
           a complete, circular sequence if (and only if) anvi'o is going to ROTATE it, which happens
           when the reference is auto-selected, or when `--use-dnaa-for-reference-orientation` is
           used. Rotating a linear sequence or a fragment of a chromosome is obviously meaningless.
           But a single-contig genomic locus that is not circular at all is a perfectly good reference
           for someone who only wants to orient/order other contigs against it, and in that case anvi'o
           will not rotate anything (which is why `--reference` alone never triggers a rotation).

           The only way out of the single-contig requirement is `--use-auto-reference-as-is`, which
           also skips every rotation step. Even then the user gets a warning, since a multi-contig
           reference still does not offer that single continuous coordinate system.

        Parameters
        ==========
        user_specified : bool
            Whether the reference genome was set by the user through `--reference`, or was chosen
            by anvi'o automatically. This only influences which error message the user gets.
        """

        num_contigs = self.genomes[self.reference_name]['num_contigs']

        if num_contigs == 1:
            return

        if self.use_auto_reference_as_is:
            # the user explicitly asked anvi'o to not tinker with the reference, so there will be no
            # rotation, and thus nothing biologically meaningless will happen to the reference. but
            # they still deserve to know what they are getting (and what they are not getting):
            self.run.warning(f"The reference genome anvi'o selected for you, '{self.reference_name}', is not a single "
                             f"contig :/ It is composed of {P('contig', num_contigs)}, and normally anvi'o would have "
                             f"refused to work with it. But since you used the flag `--use-auto-reference-as-is`, anvi'o "
                             f"will use it but NOT rotate it to any start position (which is exactly what that "
                             f"flag asks for), so nothing biologically meaningless will be done to it. That said, "
                             f"please keep the following in mind while you are interpreting the results: (1) the only "
                             f"thing anvi'o is really doing for you here is putting all your sequences on a consistent "
                             f"strand, (2) alignment coordinates that come from different contigs of the reference do "
                             f"NOT form a single continuous axis, so the contig ordering, the output of "
                             f"`--scaffold-fragmented`, and the 'Start in reference' / 'Start in query' numbers you "
                             f"will see below should not be trusted, and (3) the trust labels in the final report are "
                             f"computed on that same shaky basis, so they are much weaker statements than they would "
                             f"be with a single-contig reference. In an ideal world you would add a single-contig "
                             f"reference to your fasta-txt file and point anvi'o to it with `--reference`. That "
                             f"reference does not have to be a complete circular genome, by the way. If you are "
                             f"working with a genomic region rather than whole chromosomes, a single contig that "
                             f"covers the locus of interest will do just fine.",
                             header="YOUR REFERENCE IS NOT A SINGLE CONTIG", lc='yellow')
            return

        if user_specified:
            raise ConfigError(f"The reference must be a single contig, but the one you chose with `--reference`, "
                              f"'{self.reference_name}', is composed of {P('contig', num_contigs)} :/ This will "
                              f"not work as everything this program does for you needs a single, continuous "
                              f"coordinate system to work against, and coordinates that come from different contigs "
                              f"of a fragmented reference can't offer that. IS THERE A SOLUTION? There is no "
                              f"solution that will maek a multi-contig reference work for most applications, but "
                              f"please look at the online help for this program anyway.")

        genomes_with_a_single_contig = [g for g in self.genomes if self.genomes[g]['num_contigs'] == 1]

        if not len(genomes_with_a_single_contig):
            raise ConfigError(f"Anvi'o auto-selected '{self.reference_name}' as your reference since it had the fewest "
                              f"contigs of the bunch, but it still is composed of {P('contig', num_contigs)}. In fact, "
                              f"not a single one of the {P('entr', len(self.genomes), sfp='ies', sfs='y')} in your "
                              f"fasta-txt file is a single contig, so there is nothing anvi'o can use here as a "
                              f"reference :( Two things are going wrong at once here. First, anvi'o needs the reference to be a single "
                              f"contig, since orienting and ordering contigs requires a single, continuous coordinate "
                              f"system to work against. Second, when anvi'o auto-selects a reference it also ROTATES it "
                              f"to a start position that is conserved across your dataset, and rotating a fragment of a "
                              f"chromosome is a biologically meaningless operation. Your options: (1) add a "
                              f"single-contig reference to your fasta-txt file and point anvi'o to it with "
                              f"`--reference` (this can be a complete circular genome if you are working with whole "
                              f"chromosomes, but it can just as well be a single contig covering a genomic locus of "
                              f"interest, since anvi'o never rotates a reference you name explicitly), or (2) if all you "
                              f"want is for everything to end up on the same strand with no rotation whatsoever, use "
                              f"the flag `--use-auto-reference-as-is` (and please do read what it does first, since the "
                              f"results will be much more limited, and for that kind of reading, the online help is likely "
                              f"much more useful here than the program help menu).")

        # if we are here, anvi'o somehow managed to pick a multi-contig reference even though there were
        # single-contig genomes to choose from. so, if this is happening, someone changed something
        # somehwere they shouldn't have and we need to take a look:
        raise ConfigError("Something unexpected happened, and anvi'o needs a visit from its programmers :(")


    def _get_total_length(self, fasta_path):
        fasta = utils.u.SequenceSource(fasta_path)
        total = 0
        while next(fasta):
            total += len(fasta.seq)
        return total
