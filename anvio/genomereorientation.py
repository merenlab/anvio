#!/usr/bin/env python
# coding: utf-8
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
    def __init__(self, name, status, message, output_path=None, trust=None):
        self.name = name
        self.status = status
        self.message = message
        self.output_path = output_path
        self.trust = trust


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
        self.scaffold_fragmented = A('scaffold_fragmented') or False
        self.min_contig_length = A('min_contig_length') or 1000

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
        reference_was_user_specified = self.reference_name is not None
        self.select_reference_genome()

        # Decide how to orient the reference
        if self.use_dnaa_for_reference_orientation:
            # Use DnaA gene to orient the reference
            dnaa_position = self._find_dnaa_gene_position()
            if dnaa_position is not None and dnaa_position != 0:
                self.run.warning(None, header="ROTATING REFERENCE TO DnaA POSITION")
                self.run.info("Rotating reference by", f"{dnaa_position:,} bp")
                self.run.info("Reason", "Aligning genome to DnaA gene (origin of replication)")

                # Create a rotated version of the reference
                temp_dir = filesnpaths.get_temp_directory_path()
                rotated_ref = os.path.join(temp_dir, f"{self.reference_name}_rotated_dnaa.fa")
                self._seqkit_rotate(self.reference_path, dnaa_position + 1, rotated_ref)

                # Update reference path to the rotated version
                self.reference_path = rotated_ref
                self.run.info("Rotated reference", rotated_ref, nl_after=1)
            elif dnaa_position == 0:
                self.run.info("Reference rotation", "Not needed (DnaA gene is already at position 0)", nl_after=1)
            else:
                self.run.warning("DnaA gene not found. Will proceed without rotating the reference.", nl_after=1)

        # If reference was auto-selected and not using DnaA, find the optimal starting position
        elif not reference_was_user_specified:
            result = self._find_optimal_reference_start()
            if result:
                optimal_position, genomes_covering, total_genomes = result
                coverage_pct = (genomes_covering / total_genomes) * 100 if total_genomes > 0 else 0

                if optimal_position != 0:
                    self.run.warning(None, header="ROTATING REFERENCE TO OPTIMAL START")
                    self.run.info("Rotating reference by", f"{optimal_position:,} bp")
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

        self.run.info("Reference genome", self.reference_name)
        self.run.info("Output directory", self.output_dir)
        if self.log_run.log_file_path:
            self.run.info("Log file", self.log_run.log_file_path)

        filesnpaths.gen_output_directory(self.output_dir, run=self.run)

        reference_output_path = self._output_path_for(self.reference_name)
        shutil.copyfile(self.reference_path, reference_output_path)

        rotation_msg = "Copied reference genome without reorientation."
        if self.use_dnaa_for_reference_orientation:
            rotation_msg = "Reference genome (rotated to DnaA gene position)."
        elif not reference_was_user_specified:
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

        # Validate reference is single-contig if user-specified
        if self.reference_name:
            ref_contigs = self.genomes[self.reference_name]['num_contigs']
            if ref_contigs != 1:
                raise ConfigError(f"Reference genome '{self.reference_name}' must be a single-contig "
                                  f"circular genome, but it has {ref_contigs} contigs. Please choose a "
                                  f"different reference or let anvi'o auto-select one.")

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

        return ReorientationResult(genome_name, "ok", message, output_path, trust=trust_label)


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
        self.run.info("Contigs processed", num_total)
        self.run.info("Contigs aligned", num_aligned)
        self.run.info("Contigs unaligned", result_data['num_contigs_unaligned'])
        self.run.info("Reference coverage", f"{ref_cov:.1f}%")
        self.run.info("Average ANI", f"{avg_ani:.1f}%")
        self.run.info("Total gaps", len(result_data['gaps']))
        self.run.info("Total gap size", f"{result_data['total_gap_size']:,} bp", nl_after=1)

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

        return ReorientationResult(genome_name, "ok", message, output_path, trust=trust_label)


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
            actions.append(f"rotated {cut0 + 1} bp (pass1 snap)")

            self.progress.update(f"{genome_name}: Orienting and rotating (pass 1)")
            oriented_query = os.path.join(temp_dir, "query_oriented.fa")
            if anchor1.strand == "-":
                self._seqkit_reverse_complement(query_fa, oriented_query)
            else:
                shutil.copyfile(query_fa, oriented_query)

            rot1 = os.path.join(temp_dir, "query_rot1.fa")
            self._seqkit_rotate(oriented_query, cut0 + 1, rot1)

            self.progress.update(f"{genome_name}: Re-aligning (pass 2)")
            paf2 = self._minimap2_align(self.reference_path, rot1)
            anchor2 = self._select_anchor_smallest_reference_start(paf2)
            cut1 = self._cut0_for_ref0(anchor2)
            self.log_run.info_single(f"Pass-2 anchor strand={anchor2.strand} cut1={cut1}", level=2)
            actions.append(f"rotated {cut1 + 1} bp (pass2 snap)")

            self.progress.update(f"{genome_name}: Rotating (pass 2)")
            rot2 = os.path.join(temp_dir, "query_rot2.fa")
            self._seqkit_rotate(rot1, cut1 + 1, rot2)

            self.progress.update(f"{genome_name}: Re-aligning (pass 3)")
            paf3 = self._minimap2_align(self.reference_path, rot2)
            anchor_snap = self._select_anchor_for_ref0(paf3)
            cut_snap = self._cut0_for_ref0(anchor_snap)
            actions.append(f"rotated {cut_snap + 1} bp (snap-to-ref0 anchor)")

            self.progress.update(f"{genome_name}: Rotating (pass 3)")
            rot3 = os.path.join(temp_dir, "query_rot3.fa")
            self._seqkit_rotate(rot2, cut_snap + 1, rot3)

            self.progress.update(f"{genome_name}: Final alignment (pass 4)")
            paf4 = self._minimap2_align(self.reference_path, rot3)
            best_final = self._best_primary_alignment(paf4)

            # final snap: force reference position 0 to align to query position 0
            self.progress.update(f"{genome_name}: Final rotation snap")
            cut_final = self._cut0_for_ref0(best_final)
            final_fa = os.path.join(temp_dir, "query_rot_final.fa")
            if cut_final == 0:
                shutil.copyfile(rot3, final_fa)
            else:
                self._seqkit_rotate(rot3, cut_final + 1, final_fa)
            actions.append(f"rotated {cut_final + 1} bp (final snap to ref0)")

            self.progress.update(f"{genome_name}: Verification alignment")
            paf_final = self._minimap2_align(self.reference_path, final_fa)
            best_final_final = self._best_primary_alignment(paf_final)

            # Iteratively correct until alignment starts at reference position 0
            current_fa = final_fa
            max_correction_iterations = 5
            correction_iteration = 0

            self.log_run.info_single(f"Initial: tstart={best_final_final.tstart} qstart={best_final_final.qstart}", level=2)

            while best_final_final.tstart != 0 and correction_iteration < max_correction_iterations:
                correction_iteration += 1
                self.progress.update(f"{genome_name}: Correction iteration {correction_iteration}")
                self.log_run.info_single(f"Correction iteration {correction_iteration}: tstart={best_final_final.tstart} qstart={best_final_final.qstart} strand={best_final_final.strand}", level=2)

                # For circular genomes, we need to rotate the query so that
                # the query position that aligns with ref[0] becomes query[0]
                #
                # Current state: query[qstart] aligns with ref[tstart]
                # Goal: make ref[0] align with query[0]
                #
                # For + strand:
                # If query[qstart] aligns with ref[tstart], then:
                # - ref[0] aligns with query[qstart - tstart] (handle negative with modulo)
                if best_final_final.strand == "+":
                    # Position in query that aligns with ref[0]
                    cut_correction = (best_final_final.qstart - best_final_final.tstart) % best_final_final.qlen
                else:
                    # For - strand, use the existing logic which handles reverse complement
                    cut_correction = self._cut0_for_ref0(best_final_final)

                self.log_run.info_single(f"Calculated cut_correction={cut_correction} (qstart={best_final_final.qstart}, qend={best_final_final.qend})", level=2)
                corrected_fa = os.path.join(temp_dir, f"query_corrected_{correction_iteration}.fa")

                if cut_correction == 0:
                    self.log_run.info_single("Cut correction is 0, alignment already starts at position 0", level=2)
                    break

                self._seqkit_rotate(current_fa, cut_correction + 1, corrected_fa)
                actions.append(f"rotated {cut_correction + 1} bp (correction iter {correction_iteration} for tstart={best_final_final.tstart})")

                # Re-align after correction
                paf_final = self._minimap2_align(self.reference_path, corrected_fa)
                best_final_final = self._best_primary_alignment(paf_final)
                current_fa = corrected_fa

                self.log_run.info_single(f"After correction {correction_iteration}: tstart={best_final_final.tstart} qstart={best_final_final.qstart}", level=2)

            if best_final_final.tstart == 0:
                self.log_run.info_single(f"Successfully aligned to tstart=0 after {correction_iteration} correction(s)", level=2)
            else:
                self.log_run.info_single(f"Warning: Could not achieve tstart=0 after {correction_iteration} attempts (final tstart={best_final_final.tstart})", level=2)

            shutil.copyfile(current_fa, output_path)
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
            contig_alignments = {}
            ref_length = self._get_total_length(self.reference_path)
            num_wrap_around_contigs = 0

            for idx, contig in enumerate(contigs):
                # Write contig to temp file
                temp_contig_fa = os.path.join(temp_dir, f"contig_{idx}.fa")
                with open(temp_contig_fa, 'w') as f:
                    f.write(f">{contig['id']}\n{contig['seq']}\n")

                # Align to reference
                try:
                    paf_records = self._minimap2_align(self.reference_path, temp_contig_fa)

                    if not paf_records:
                        self.log_run.info_single(f"'{contig['id']}': no alignment", level=2)
                        continue

                    # Check for wrap-around: contig spans the circularization point
                    # This happens when we have alignments at both ends of the reference
                    wrap_around_detected = False

                    if len(paf_records) >= 2:
                        # Sort alignments by reference start position
                        sorted_aligns = sorted(paf_records, key=lambda r: r.tstart)

                        # Check if we have alignment near start AND near end of reference
                        # Near start: tstart < 10% of ref_length
                        # Near end: tstart > 90% of ref_length
                        near_start = [a for a in sorted_aligns if a.tstart < ref_length * 0.1]
                        near_end = [a for a in sorted_aligns if a.tstart > ref_length * 0.9]

                        if near_start and near_end:
                            # Likely wrap-around: check if alignments together cover most of contig
                            total_aligned = sum(a.aligned_bases for a in paf_records)
                            if total_aligned > contig['length'] * 0.7:  # At least 70% of contig aligned
                                wrap_around_detected = True
                                num_wrap_around_contigs += 1

                                # Get the two main alignments (one near end, one near start)
                                end_align = max(near_end, key=lambda r: r.aligned_bases)
                                start_align = max(near_start, key=lambda r: r.aligned_bases)

                                self.log_run.info_single(
                                    f"'{contig['id']}': WRAP-AROUND detected! Splitting into 2 parts", level=2)
                                self.log_run.info_single(
                                    f"  Part 1 (end): ref[{end_align.tstart}:{end_align.tend}] "
                                    f"qry[{end_align.qstart}:{end_align.qend}] strand={end_align.strand}", level=2)
                                self.log_run.info_single(
                                    f"  Part 2 (start): ref[{start_align.tstart}:{start_align.tend}] "
                                    f"qry[{start_align.qstart}:{start_align.qend}] strand={start_align.strand}", level=2)

                                # Extract the two parts of the contig based on query positions
                                # Note: Need to handle strand orientation
                                if end_align.strand == '+':
                                    # Forward strand: extract by query positions directly
                                    part1_seq = contig['seq'][end_align.qstart:end_align.qend]
                                else:
                                    # Reverse strand: query positions are on RC, need to extract and RC
                                    part1_seq = contig['seq'][end_align.qstart:end_align.qend]

                                if start_align.strand == '+':
                                    part2_seq = contig['seq'][start_align.qstart:start_align.qend]
                                else:
                                    part2_seq = contig['seq'][start_align.qstart:start_align.qend]

                                # Create two contig entries
                                part1_id = f"{contig['id']}_wrapPart1"
                                part2_id = f"{contig['id']}_wrapPart2"

                                contig_alignments[part1_id] = {
                                    'contig_data': {
                                        'id': part1_id,
                                        'seq': part1_seq,
                                        'length': len(part1_seq),
                                        'original_id': contig['id']
                                    },
                                    'alignment': end_align
                                }

                                contig_alignments[part2_id] = {
                                    'contig_data': {
                                        'id': part2_id,
                                        'seq': part2_seq,
                                        'length': len(part2_seq),
                                        'original_id': contig['id']
                                    },
                                    'alignment': start_align
                                }

                    # If not wrap-around, process normally
                    if not wrap_around_detected:
                        primaries = [r for r in paf_records if r.is_primary]

                        if primaries:
                            # Select best primary alignment
                            best = max(primaries, key=lambda r: (r.aligned_bases, r.mapq, r.nmatch))
                            contig_alignments[contig['id']] = {
                                'contig_data': contig,
                                'alignment': best
                            }
                            self.log_run.info_single(
                                f"'{contig['id']}': ref[{best.tstart}:{best.tend}] "
                                f"strand={best.strand} alen={best.aligned_bases}", level=2)
                        else:
                            self.log_run.info_single(f"'{contig['id']}': no primary alignment", level=2)

                except RuntimeError as e:
                    self.log_run.info_single(f"'{contig['id']}': alignment failed ({e})", level=2)

            # Count unique original contigs that aligned
            # (accounting for split wrap-around contigs)
            aligned_original_contigs = set()
            for contig_id, info in contig_alignments.items():
                original_id = info['contig_data'].get('original_id', info['contig_data']['id'])
                aligned_original_contigs.add(original_id)

            num_aligned_original = len(aligned_original_contigs)
            num_unaligned = num_contigs_after_filter - num_aligned_original
            num_alignment_pieces = len(contig_alignments)  # Total pieces (including split contigs)

            self.log_run.info_single(
                f"Alignment summary: {num_aligned_original} contigs aligned "
                f"({num_alignment_pieces} pieces total, {num_wrap_around_contigs} split for wrap-around)",
                level=2)

            if num_aligned_original == 0:
                raise ConfigError("No contigs aligned to reference")

            # Step 3: Order contigs by reference position
            ordered_contig_ids = sorted(
                contig_alignments.keys(),
                key=lambda cid: (
                    contig_alignments[cid]['alignment'].tstart,
                    -contig_alignments[cid]['alignment'].aligned_bases,
                    -contig_alignments[cid]['alignment'].mapq
                )
            )

            # Step 4: Orient contigs and calculate gaps
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
                            f"Gap: {gap_size} bp between '{contig_id}' and '{next_contig_id}'",
                            level=2)
                    elif gap_end < gap_start:
                        overlap_size = gap_start - gap_end
                        self.log_run.info_single(
                            f"WARNING: {overlap_size} bp overlap between '{contig_id}' "
                            f"and '{next_contig_id}'", level=2)
                        gaps.append({
                            'after_contig': contig_id,
                            'gap_size': 0,
                            'overlap': overlap_size
                        })

            # Step 5: Write output FASTA
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
            else:
                # Write as separate contigs (ordered and oriented)
                with open(output_path, 'w') as out_fa:
                    for idx, contig_info in enumerate(oriented_contigs):
                        out_fa.write(f">{contig_info['id']}\n")
                        # Write in 80-character lines
                        seq = contig_info['seq']
                        for i in range(0, len(seq), 80):
                            out_fa.write(seq[i:i+80] + '\n')

            # Step 6: Calculate quality metrics
            total_gap_size = sum(g['gap_size'] for g in gaps if 'overlap' not in g)

            # Reference coverage
            ref_covered_regions = [(contig_alignments[cid]['alignment'].tstart,
                                    contig_alignments[cid]['alignment'].tend)
                                   for cid in ordered_contig_ids]
            ref_covered_regions.sort()

            # Merge overlapping intervals
            merged_regions = [ref_covered_regions[0]]
            for start, end in ref_covered_regions[1:]:
                if start <= merged_regions[-1][1]:
                    merged_regions[-1] = (merged_regions[-1][0], max(merged_regions[-1][1], end))
                else:
                    merged_regions.append((start, end))

            total_covered_bases = sum(end - start for start, end in merged_regions)
            ref_length = self._get_total_length(self.reference_path)
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

            actions_summary = f"ordered {num_aligned_original} contigs by reference position"
            if num_wrap_around_contigs > 0:
                actions_summary += f" ({num_wrap_around_contigs} split for wrap-around)"
            if self.scaffold_fragmented:
                actions_summary += f", inserted {total_gap_size:,} bp of N-padding"

            # Return results
            return {
                'num_contigs_processed': num_contigs_after_filter,
                'num_contigs_aligned': num_aligned_original,
                'num_contigs_unaligned': num_unaligned,
                'reference_coverage_pct': reference_coverage_pct,
                'gaps': gaps,
                'total_gap_size': total_gap_size,
                'avg_ani': avg_ani,
                'actions_summary': actions_summary,
                'num_wrap_around_contigs': num_wrap_around_contigs
            }

        finally:
            if not anvio.DEBUG:
                shutil.rmtree(temp_dir)


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
            # For - strand: query is reverse complemented
            # query[qend-1] aligns with ref[tstart] (coordinates are on forward strand)
            # Position that aligns with ref[0]: qend-1 - tstart (but on reverse)
            # After reverse complement rotation correction:
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
        plt.xlabel("Genome position (bp)")
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
        self.run.info("Reference length", f"{ref_len:,} bp")
        self.run.info("Bin size for coverage", f"{bin_size:,} bp")
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

        # Find the position with maximum coverage (in a greedy fashion)
        # Use a sliding window to smooth out noise and find a good anchor region
        window_size = min(1000, ref_len // 100)  # 1kb or 1% of genome, whichever is smaller
        if window_size < 100:
            window_size = min(100, ref_len)

        # Convert window size to number of bins
        window_bins = max(1, window_size // bin_size)

        self.progress.new("Finding optimal start position", progress_total_items=num_bins)
        best_score = -1
        best_position = 0
        best_window_genomes = set()

        # Search at bin resolution rather than base-pair resolution
        for bin_idx in range(num_bins):
            if bin_idx % 10 == 0:  # Update every 10 bins
                self.progress.update(f"Scanning bin {bin_idx:,} / {num_bins:,}")
                self.progress.increment(increment_to=bin_idx)

            # Calculate coverage in a window of bins around this bin
            # Coverage = number of unique genomes covering the bins in the window
            window_genomes = set()
            for offset in range(-window_bins // 2, window_bins // 2):
                bin_pos = (bin_idx + offset) % num_bins  # Circular genome
                window_genomes.update(coverage_bins[bin_pos])

            avg_coverage = len(window_genomes)

            if avg_coverage > best_score:
                best_score = avg_coverage
                # Convert bin index back to position (use middle of bin)
                best_position = bin_idx * bin_size + bin_size // 2
                # Make sure we don't go past the reference length
                if best_position >= ref_len:
                    best_position = bin_idx * bin_size
                best_window_genomes = window_genomes.copy()

                # Early termination: if all genomes cover this position, we can't do better. Even though
                # we may have to think about this more carefully. As in not stopping immediately after
                # the first match, but stopping after three consecutive matches or something. Difficult
                # to explain why, but I know why. TRUST MEH.
                if avg_coverage == total_genomes:
                    self.progress.update(f"Found optimal position at {best_position:,} (100% coverage)")
                    self.log_run.info_single(f"Early termination: found position with 100% coverage at {best_position:,}", level=2)
                    break

        self.progress.end()

        # Report the window coverage that was actually used for selection
        genomes_covering = len(best_window_genomes)
        coverage_pct = (genomes_covering / total_genomes) * 100 if total_genomes > 0 else 0

        self.run.info("Optimal start position", f"{best_position:,} bp")
        self.run.info("Window size used for selection", f"{window_size:,} bp")
        self.run.info("Genomes covering this position", f"{genomes_covering}/{total_genomes} ({coverage_pct:.1f}%)")

        return best_position, genomes_covering, total_genomes


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
            self.run.info("DnaA location", f"{dnaa_start:,} - {dnaa_stop:,} bp")
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

        if self.reference_name:
            self.run.info("Reference genome", f"'{self.reference_name}' (specified by the user)")
            self.reference_path = self.genomes[self.reference_name]['path']
            return

        candidate_stats = []
        for genome_name, entry in self.genomes.items():
            path = entry['path']
            num_sequences = utils.get_num_sequences_in_fasta(path)
            total_len = self._get_total_length(path)
            candidate_stats.append((num_sequences, -total_len, genome_name, path))

        candidate_stats.sort()
        chosen = candidate_stats[0]
        self.reference_name = chosen[2]
        self.reference_path = chosen[3]
        chosen_contigs = chosen[0]
        chosen_len = -chosen[1]
        self.run.info("Reference genome", f"'{self.reference_name}' (selected by anvi'o)")
        self.run.info("Reference contigs", chosen_contigs)
        self.run.info("Reference length", chosen_len)


    def _get_total_length(self, fasta_path):
        fasta = utils.u.SequenceSource(fasta_path)
        total = 0
        while next(fasta):
            total += len(fasta.seq)
        return total
