#!/usr/bin/env python
# coding: utf-8
"""Reorient circular genomes to match a reference genome coordinate system."""

import os
import shutil

from pathlib import Path

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

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

        # If reference was auto-selected, find the optimal starting position
        if not reference_was_user_specified:
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
        if not reference_was_user_specified:
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

            try:
                best_alignment, paf_final, _, actions = self._reorient_single(fasta_path, output_path, genome_name)
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
                self.run.info("Applied actions", ", ".join(actions))
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
                results.append(ReorientationResult(genome_name, "ok", message, output_path, trust=trust_label))
                self._plot_dotplot(paf_final, genome_name, label="After reorientation")
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

        self.run.warning("Anvi'o used `minimap2` by Li (doi:10.1093/bioinformatics/bty191) to rapidly align genomes "
                         "and `seqkit` by Shen et al (doi:10.1371/journal.pone.0163962) to rotate and reverse/complement "
                         "large sequences. If you publish your findings, please do not forget to properly credit the "
                         "authors of these tools.", lc='cyan', header="CITATION")

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

        for genome_name, entry in self.genomes.items():
            utils.is_this_name_OK_for_database('fasta.txt entry name', genome_name, additional_chars_allowed='.-')

            genome_path = os.path.abspath(entry['path'])
            filesnpaths.is_file_exists(genome_path)
            filesnpaths.is_file_fasta_formatted(genome_path)

            num_sequences = utils.get_num_sequences_in_fasta(genome_path)
            if num_sequences != 1:
                raise ConfigError(f"The FASTA for '{genome_name}' must contain exactly one sequence (found {num_sequences}). "
                                  "Anvi'o currently expects circular single-contig genomes here.")

            self.genomes[genome_name]['path'] = genome_path

        self.output_dir = filesnpaths.check_output_directory(self.output_dir, ok_if_exists=False)

        if self.threads < 1:
            raise ConfigError("Number of threads must be a positive integer.")


    def _output_path_for(self, genome_name):
        suffix = Path(self.genomes[genome_name]['path']).suffix or ".fa"
        return os.path.join(self.output_dir, f"{genome_name}{suffix}")


    def _reorient_single(self, query_fa, output_path, genome_name):
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


    def _minimap2_align(self, ref_fa, qry_fa):
        cmd = ["minimap2",
               "-x",
               self.minimap2_preset,
               "-t",
               str(self.threads),
               ref_fa,
               qry_fa]
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
            final_text = (f"Your genome reorientation task considered {num_reoriented} genomes and {num_reference} "
                          f"reference. Some outcomes are not trustworthy (low alignment coverage can lead to unreliable "
                          f"orientation). Please review the dotplots above and the summary below to decide which FASTA "
                          f"files to use downstream.")
        elif num_somewhat:
            final_text = (f"Your genome reorientation task considered {num_reoriented} genomes and {num_reference} "
                          f"reference. Some outcomes are only somewhat OK; check the dotplots and stats below before "
                          f"using them downstream.")
        else:
            final_text = (f"Your genome reorientation task for {num_reoriented} genomes and {num_reference} reference "
                          f"is complete. All outcomes look trustworthy based on alignment coverage. A quick glance at "
                          f"the dotplots above is still recommended before downstream analyses.")

        self.run.warning(final_text, header="FINAL REPORT")

        trustworthy_list = trust_counts["REFERENCE"] + trust_counts["TRUSTWORTHY"]
        self.run.info("Trustworthy", len(trustworthy_list))
        for name, path in trustworthy_list:
            color = "cyan" if (name, path) in trust_counts["REFERENCE"] else "green"
            self.run.info_single(f"{name} -> {path}", mc=color, level=2)

        self.run.info("Somewhat OK", len(trust_counts["SOMEWHAT OK"]), nl_before=1)
        for name, path in trust_counts["SOMEWHAT OK"]:
            self.run.info_single(f"{name} -> {path}", mc="yellow", level=2)

        self.run.info("Not trustworthy", len(trust_counts["NOT TRUSTWORTHY"]), nl_before=1)
        for name, path in trust_counts["NOT TRUSTWORTHY"]:
            self.run.info_single(f"{name} -> {path}", mc="red", level=2)

        self.run.info("Failed", len(trust_counts["FAILED"]), nl_before=1)
        for name, path in trust_counts["FAILED"]:
            self.run.info_single(f"{name} -> {path}", mc="red", level=2)


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


    def _plot_dotplot(self, recs, genome_name, label):
        try:
            import plotext as plt
        except ImportError:
            self.run.warning(f"plotext is not available; skipping dotplot for '{genome_name}'.")
            return

        primaries = [r for r in recs if r.is_primary]
        if not primaries:
            self.run.warning(f"No primary alignments to plot for '{genome_name}'.")
            return

        primaries = sorted(primaries, key=lambda r: r.aligned_bases, reverse=True)[:500]

        plus_x, plus_y, minus_x, minus_y = [], [], [], []
        for r in primaries:
            ref_mid = (r.tstart + r.tend) / 2.0
            qry_mid = (r.qstart + r.qend) / 2.0
            if r.strand == '+':
                plus_x.append(ref_mid)
                plus_y.append(qry_mid)
            else:
                minus_x.append(ref_mid)
                minus_y.append(r.qlen - qry_mid)

        plt.clf()
        plt.plotsize(80, 20)
        plt.theme("dark")
        if plus_x:
            plt.scatter(plus_x, plus_y, label='+ strand', color='lightgreen')
        if minus_x:
            plt.scatter(minus_x, minus_y, label='- strand', color='orange')
        # draw reference diagonal for visual guide
        ref_len = primaries[0].tlen
        qry_len = primaries[0].qlen
        diag_len = min(ref_len, qry_len)
        plt.plot([0, diag_len], [0, diag_len], color='gray')
        plt.xlabel("Reference start")
        plt.ylabel("Query start")
        plt.title(f"{label}: {genome_name} vs {self.reference_name}")
        try:
            plt.xlim(0, ref_len)
            plt.ylim(0, qry_len)
            x_ticks = [0, ref_len * 0.25, ref_len * 0.5, ref_len * 0.75, ref_len]
            y_ticks = [0, qry_len * 0.25, qry_len * 0.5, qry_len * 0.75, qry_len]
            plt.xticks(x_ticks, [utils.human_readable_number(x, decimals=1) for x in x_ticks])
            plt.yticks(y_ticks, [utils.human_readable_number(y, decimals=1) for y in y_ticks])
        except Exception:
            pass
        try:
            plt.legend(True)
        except Exception:
            pass
        plt.show()


    def _find_optimal_reference_start(self):
        """Survey all genomes to find the optimal starting position in the reference.

        Returns (optimal_position, num_genomes_covering, total_genomes) or None if reference is user-specified.
        """
        # Get reference length
        ref_len = self._get_total_length(self.reference_path)

        # Initialize coverage array: for each position, track SET of genomes that cover it
        coverage = [set() for _ in range(ref_len)]

        self.run.warning(None, header="FINDING OPTIMAL REFERENCE START")
        self.run.info("Reference genome", self.reference_name)
        self.run.info("Reference length", f"{ref_len:,} bp")

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
            self.log_run.info_single(f"Aligning {genome_name} to find coverage", level=2)
            paf_recs = self._minimap2_align(self.reference_path, query_path)

            # Mark positions covered by primary alignments
            primaries = [r for r in paf_recs if r.is_primary]
            for rec in primaries:
                for pos in range(rec.tstart, rec.tend):
                    if pos < ref_len:
                        coverage[pos].add(genome_name)

            self.progress.increment()

        self.progress.end()

        # Find the position with maximum coverage (in a greedy fashion)
        # Use a sliding window to smooth out noise and find a good anchor region
        window_size = min(1000, ref_len // 100)  # 1kb or 1% of genome, whichever is smaller
        if window_size < 100:
            window_size = min(100, ref_len)

        self.progress.new("Finding optimal start position", progress_total_items=ref_len)
        best_score = -1
        best_position = 0
        best_window_genomes = set()

        for i in range(ref_len):
            if i % 10000 == 0:  # Update every 10k positions to avoid too many updates
                self.progress.update(f"Scanning position {i:,} / {ref_len:,}")
                self.progress.increment(increment_to=i)

            # Calculate average coverage in a window around this position
            # Coverage = number of unique genomes covering each position
            window_genomes = set()
            for offset in range(-window_size // 2, window_size // 2):
                pos = (i + offset) % ref_len  # Circular genome
                window_genomes.update(coverage[pos])

            avg_coverage = len(window_genomes)

            if avg_coverage > best_score:
                best_score = avg_coverage
                best_position = i
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
