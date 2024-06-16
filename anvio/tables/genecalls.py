# -*- coding: utf-8
# pylint: disable=line-too-long

import os
import numpy

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.fastalib as u
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.filesnpaths as filesnpaths
import anvio.genecalling as genecalling

from anvio.tables.tableops import Table
from anvio.errors import ConfigError


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


pp = terminal.pretty_print


class TablesForGeneCalls(Table):
    def __init__(self, db_path, contigs_fasta=None, args=None, run=terminal.Run(), progress=terminal.Progress(), debug=False):
        self.run = run
        self.progress = progress

        self.db_path = db_path
        self.contigs_fasta = contigs_fasta
        self.args = args
        self.debug = debug

        utils.is_contigs_db(self.db_path)

        if self.contigs_fasta:
            filesnpaths.is_file_exists(self.contigs_fasta)
            filesnpaths.is_file_fasta_formatted(self.contigs_fasta)


    def check_gene_calls_dict(self, gene_calls_dict):
        if not isinstance(gene_calls_dict, type({})):
            raise ConfigError("Gene calls dict must be a dict instance :/")

        try:
            [int(g) for g in list(gene_calls_dict.keys())]
        except ValueError:
            raise ConfigError("Keys of a gene calls dict must be integers!")

        if False in [x['direction'] in ['f', 'r'] for x in list(gene_calls_dict.values())]:
            raise ConfigError("The values in 'direction' column can't be anything but 'f' (for forward) "
                               "or 'r' (for reverse). You have other stuff, and it is not cool.")

        if len([True for x in list(gene_calls_dict.values()) if 'call_type' not in x]):
            raise ConfigError("Gene calls must have key 'call_type'.")

        call_types_known = set(list(constants.gene_call_types.values()))
        if len([True for x in list(gene_calls_dict.values()) if x['call_type'] not in call_types_known]):
            raise ConfigError("Call types must be one of these (note, they're all integers): %s." % (', '.join([str(c) for c in call_types_known])))

        if False in [x['stop'] > x['start'] for x in list(gene_calls_dict.values())]:
            raise ConfigError("For each gene call, the stop position must be bigger than the start position. "
                               "Your gene calls dict does not conform to that. If you have reverse gene calls "
                               "you must use the 'direction' column to declare that.")

        if False in [(x['stop'] - float(x['start'])) % 3.0 == 0 for x in list(gene_calls_dict.values())]:
            raise ConfigError("Something is wrong with your gene calls. For every gene call, the (stop - start) "
                               "should be multiply of 3. It is not the case for all, which is a deal breaker.")


    def use_external_gene_calls_to_populate_genes_in_contigs_table(self, input_file_path, gene_calls_dict=None, ignore_internal_stop_codons=False,
                                                                   skip_predict_frame=False, skip_amino_acid_sequences=False):
        """Add genes to the contigs database.

        Primary input is either an `input_file_path` for external gene calls, or an
        external `gene_calls_dict` dictionary object.

        Parameters
        ==========
        input_file_path : str
            Path to file with one of the following structures.

            Option 1:
                gene_callers_id contig          start  stop  direction  partial  call_type  source    version
                0               CACHJY01_00016  0      693   r          1        1          prodigal  v2.6.3
                1               CACHJY01_00016  711    1140  r          0        1          prodigal  v2.6.3

            Option 2:
                gene_callers_id contig          start  stop  direction  partial  call_type  source    version  aa_sequence
                0               CACHJY01_00016  0      693   r          1        1          prodigal  v2.6.3   MSKKIYFTEYSKVNRLQTISNFTGSA
                1               CACHJY01_00016  711    1140  r          0        1          prodigal  v2.6.3   MVNVDYHGLIAGAGSGKTKVLTSRIAHIIK

        gene_calls_dict : dict, None
            Alternative to `input_file_path`. If provided, entries will be APPENDED to the database.
            So you need to make sure gene caller ids in your dict does not overlap with the ones in
            the database. Should look like:

                {
                    "1": {
                        "contig": "contig_name",
                        "start": 20,
                        "stop": 1544,
                        "direction": "f",
                        "partial": 0,
                        "call_type": 1,
                        "source": "source_name",
                        "version": "unknown",
                        "aa_sequence": "MSKKIYFTEYSKVNRLQTISNFTGSA"
                    },

                    "2": {
                      (...)
                    },

                    (...)
                }

            All entries are required except "aa_sequence", which is optional. If provided, it should
            be present for ALL entries, even if it is an empty string. It's presence will be used to
            populate `gene_amino_acid_sequences`.

        ignore_internal_stop_codons : bool, False
            If False, ConfigError will be raised if a stop codon is found inside any gene. If True,
            this is suppressed and the stop codon is replaced with the character `X`.

        skip_predict_frame : bool, False
            If True, ConfigError will be raised if a gene is not divisible by 3. If False, anvi'o predicts
            the most likley open reading frame and trims the start/stop of the gene call to reflect this
            change so that the gene *is* divisible by 3. This flag allows the retention of amino acid
            sequences even if genes are not divisible by 3, or when it is flagged as partial.

        skip_amino_acid_sequences : bool, False
            Should the gene_amino_acid_sequences table be populated? This may be useful if genes
            that are not translated are being added, such as ribosomal RNA genes, etc.
        """

        # by default we assume that this is a pristine run. but if the user sends a dictionary
        append_to_the_db = False

        gene_calls_found = False
        # let's do a rigorous check whether the user provided a gene_calls_dict.
        if (gene_calls_dict is not None and gene_calls_dict is not False):
            if not isinstance(gene_calls_dict, dict):
                raise ConfigError("'Use external gene calls' function received a non-empty gene_calls_dict object,\
                                    but it is of type '%s', and not '%s'" % (type(gene_calls_dict), type({})))

            # congrats, we have a dict.
            gene_calls_found = True

            has_aa_seq =  lambda x: True if 'aa_sequence' in x else False
            num_with_aa_seqs = sum([has_aa_seq(gene_call) for gene_call in gene_calls_dict.values()])
            num_gene_calls = len(gene_calls_dict)
            if num_with_aa_seqs != 0 and num_with_aa_seqs != num_gene_calls:
                raise ConfigError("The gene_calls_dict passed to use_external_gene_calls_to_populate_genes_in_contigs_table "
                                  "has %d entries with 'aa_sequence' and %d without. Either 0 or all (%d) should have "
                                  "'aa_sequence'" % (num_with_aa_seqs, num_gene_calls-num_with_aa_seqs, num_gene_calls))

            if not len(gene_calls_dict):
                # but it is empty ... silly user.
                self.run.info_single("'Use external gene calls' function found an empty gene calls dict, returning "
                                     "prematurely and assuming you know what's up. If you don't, stop here and try to "
                                     "identify what decisions you've made might have led you to this weird point your "
                                     "workflow (or 'life', totally up to you and your mood, but anvi'o thinks you've "
                                     "done great so far.", nl_before=1, nl_after=1)
                return


        if (not input_file_path and not gene_calls_found) or (input_file_path and gene_calls_found):
            raise ConfigError("You must provide either an input file, or an gene calls dict to process external "
                              "gene calls. You called `use_external_gene_calls_to_populate_genes_in_contigs_table` "
                              "with wrong parameters.")

        Table.__init__(self, self.db_path, anvio.__contigs__version__, self.run, self.progress, simple=True)

        # take care of gene calls dict
        if not gene_calls_found:
            expected_fields = t.genes_in_contigs_table_structure
            column_mapping = [int, str, int, int, str, int, int, str, str]

            if 'aa_sequence' in utils.get_columns_of_TAB_delim_file(input_file_path):
                expected_fields = t.genes_in_contigs_table_structure + ['aa_sequence']
                column_mapping.append(lambda x: '' if x is None else str(x)) # str(None) is 'None', amazingly

            gene_calls_dict = utils.get_TAB_delimited_file_as_dictionary(input_file_path,
                                                                         expected_fields=expected_fields,
                                                                         only_expected_fields=True,
                                                                         column_mapping=column_mapping)

            if not len(gene_calls_dict):
                raise ConfigError("You provided an external gene calls file, but it returned zero gene calls. Assuming that "
                                  "this is an error, anvi'o will stop here and complain. If this is not an error and you "
                                  "in fact expected this, the proper way of doing this is to use `--skip-gene-calls` flag, "
                                  "instead of providing an emtpy external gene calls file. You don't agree? You need this "
                                  "for some weird step for you weird pipeline? Let us know, and we will consider changing "
                                  "this.")

            self.run.info("External gene calls", "%d gene calls recovered and will be processed." % len(gene_calls_dict))
        else:
            # FIXME: we need to make sure the gene caller ids in the incoming directory is not going to
            #        overwrite an existing gene call. Something like this would have returned the
            #        current max, which could be cross-checked with what's in the dict:
            #
            #            contigs_db = ContigsDatabase(self.db_path)
            #            next_id = contigs_db.db.get_max_value_in_column('genes_in_contigs', 'gene_callers_id') + 1
            #            contigs_db.disconnect()
            append_to_the_db = True

        # recover amino acid sequences or create a blank dictionary
        if skip_amino_acid_sequences:
            amino_acid_sequences = dict([(g, '') for g in gene_calls_dict])
        else:
            gene_calls_dict, amino_acid_sequences = self.get_amino_acid_sequences_for_genes_in_gene_calls_dict(
                gene_calls_dict,
                ignore_internal_stop_codons=ignore_internal_stop_codons,
                skip_predict_frame=skip_predict_frame,
            )

        # populate genes_in_contigs, and gene_amino_acid_sequences table in contigs db.
        self.populate_genes_in_contigs_table(gene_calls_dict, amino_acid_sequences, append_to_the_db=append_to_the_db)


    def get_amino_acid_sequences_for_genes_in_gene_calls_dict(self, gene_calls_dict, ignore_internal_stop_codons=False, skip_predict_frame=False):
        """Recover amino acid sequences for gene calls in a gene_calls_dict.

        If 'aa_sequence' exists as keys in the gene_calls_dict[<key>] objects, this function will take
        those seqeunces into consideration and use them without trying to predict frames even if the gene
        call is partial. So user-defined aa sequences in `aa_sequence` column will have priority.

        Please note this FIXME: By reading all contig sequences into memory, anvi'o does a pretty bad job
        at memory management throughout this function :(

        Parameters
        ==========
        ignore_internal_stop_codons : bool, False
            If False, ConfigError will be raised if a stop codon is found inside any gene. If True,
            this is suppressed and the stop codon is replaced with the character `X`.

        skip_predict_frame : bool, False
            If True, ConfigError will be raised if a gene is not divisible by 3. If False, anvi'o predicts
            the most likley open reading frame and trims the start/stop of the gene call to reflect this
            change so that the gene *is* divisible by 3. This flag allows the retention of amino acid
            sequences even if genes are not divisible by 3, or when it is flagged as partial.
        """

        predict_frame = (not skip_predict_frame)

        if predict_frame:
            # Preload the markov model to predict frames and assign null codon and stop codon transition probabilities
            model_path = os.path.join(os.path.dirname(anvio.__file__), 'data/seq_transition_models/AA/MM_GC_0-39.npy')
            if not filesnpaths.is_file_exists(model_path, dont_raise=True):
                raise ConfigError("The task at hand calls for the use of the anvi'o Markov model to predict proper open reading "
                                  "frames for external gene calls when necessary, but the model does not seem to be in the right "
                                  "place in the anvi'o codebase. FAILING BIG HERE.")

            model = numpy.load(model_path)
            null_prob = numpy.median(model)
            stop_prob = model.min()/1e6

        gene_caller_ids_with_user_provided_amino_acid_sequences = set([])

        # get all the amino acids sorted out. either we will start with an empty dict, or take user defined
        # aa seqs as starting material
        if 'aa_sequence' in gene_calls_dict[list(gene_calls_dict.keys())[0]]:
            # the external gene calls file include amino acid sequences for at least some of the gene calls
            # here we will learn about them, and then use them we already have AA sequences
            amino_acid_sequences = {gene_caller_id: info['aa_sequence'].strip() for gene_caller_id, info in gene_calls_dict.items()}

            for gene_callers_id in amino_acid_sequences:
                if len(amino_acid_sequences[gene_callers_id]):
                    gene_caller_ids_with_user_provided_amino_acid_sequences.add(gene_callers_id)

                    gene_length = gene_calls_dict[gene_callers_id]['stop'] - gene_calls_dict[gene_callers_id]['start']
                    estimated_length_for_aa = gene_length / 3
                    user_provided_aa_length = len(amino_acid_sequences[gene_callers_id])

                    # there is already a sanity check for htis, but one can't be too careful
                    if gene_calls_dict[gene_callers_id]['call_type'] != constants.gene_call_types['CODING'] and user_provided_aa_length:
                        raise ConfigError("You have provided an amino acid sequence for at least one gene call in your external gene calls "
                                           "(%d) file that was not marked as CODING type :(" % gene_callers_id)

                    if user_provided_aa_length > estimated_length_for_aa:
                        raise ConfigError("Bad news :( There seems to be at least one gene call in your external gene calls file "
                                          "that has an aminio acid sequence that is longer than the expected length of it given the "
                                          "start/stop positions of the gene call. This is certainly true for gene call number %d "
                                          "but anvi'o doesn't know if there are more of these in your file or not :/" % gene_callers_id)

            self.run.warning("Anvi'o found amino acid sequences in your external gene calls file that match to %d of %d gene "
                             "in it and will use these amino acid sequences for everything." % (len(amino_acid_sequences), len(gene_calls_dict)))
        else:
            amino_acid_sequences = {}

        # FIXME: this is a very poor practice for memory management:
        contig_sequences = {}
        if self.contigs_fasta:
            fasta = u.SequenceSource(self.contigs_fasta)
            while next(fasta):
                contig_sequences[fasta.id] = {'sequence': fasta.seq}
            fasta.close()
        else:
            database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))
            contig_sequences = database.get_table_as_dict(t.contig_sequences_table_name)

        # keep track of things to report later
        P = lambda x: x % ("partial_" if partial else "")
        report = {"num_non_coding_gene_calls": 0,
                  "num_partial_gene_calls": 0,
                  "num_partial_gene_calls_with_user_provided_amino_acid_sequences": 0,
                  "num_gene_calls_with_user_provided_amino_acid_sequences": 0,
                  "num_partial_gene_calls_with_predicted_frame": 0,
                  "num_gene_calls_with_predicted_frame": 0,
                  "num_partial_gene_calls_with_no_predicted_frame": 0,
                  "num_gene_calls_with_no_predicted_frame": 0,
                  "num_partial_genes_with_internal_stops": 0,
                  "num_genes_with_internal_stops": 0}

        # the main loop to go through all the gene calls.
        for gene_callers_id in gene_calls_dict:
            gene_call = gene_calls_dict[gene_callers_id]
            partial = gene_call['partial']
            contig_name = gene_call['contig']

            if contig_name not in contig_sequences:
                # remove the partial contigs database so things don't get screwed later
                os.remove(self.db_path)
                raise ConfigError("You are in big trouble :( The contig name '%s' in your external gene calls file "
                                  "does not appear to be among your contigs. Rhetorical question time: "
                                  "HOW DID THIS HAPPEN?" % contig_name)

            # if this is a gene call that is not CODING, we have no interest in trying to get amino acid seqeunces for it
            if gene_calls_dict[gene_callers_id]['call_type'] != constants.gene_call_types['CODING']:
                report["num_non_coding_gene_calls"] += 1
                continue

            sequence = contig_sequences[contig_name]['sequence'][gene_call['start']:gene_call['stop']]
            if gene_call['direction'] == 'r':
                sequence = utils.rev_comp(sequence)

            # a let's keep track of partial gene calls
            if partial:
                report["num_partial_gene_calls"] += 1

            if gene_callers_id in gene_caller_ids_with_user_provided_amino_acid_sequences:
                # FIXME / NOTE: Here we actually move on with the assumption that the start/stop positions
                #               for the gene call are appropriate, and the user-provided amino acid sequence actually
                #               matches to those essential information. It may have been a better strategy to
                #               spend just a little more here test the frame, and start/stop positions. This probably
                #               will explode at some point due to some user error, and some poor soul will spend hours
                #               in the codebase to figure out how the hell did it happen.
                report[P("num_%sgene_calls_with_user_provided_amino_acid_sequences")] += 1
                amino_acid_sequence = amino_acid_sequences[gene_callers_id]
            elif predict_frame:
                # no amino acid sequence is provided, BUT USER WANTS FRAME TO BE PREDICTED
                # we may be good, if we can try to predict one for it.
                frame, amino_acid_sequence = utils.get_most_likely_translation_frame(sequence, model=model, stop_prob=stop_prob, null_prob=null_prob)

                if frame is None:
                    # we not good because we couldn't find a frame for it. because this gene call has no predicted frame,
                    # and no user-provided amino acid sequence, we will mark this one as noncoding. BAM:
                    gene_calls_dict[gene_callers_id]['call_type'] = constants.gene_call_types['NONCODING']
                    report[P("num_%sgene_calls_with_no_predicted_frame")] += 1
                    continue
                else:
                    # we good. found the amino acid sequence. we will update the gene call so start/stop
                    # matches to the frame, and report the amino acid sequence
                    report[P("num_%sgene_calls_with_predicted_frame")] += 1
                    gene_calls_dict[gene_callers_id] = self.update_gene_call(gene_call, frame)
                    amino_acid_sequences[gene_callers_id] = amino_acid_sequence
            elif not predict_frame:
                # no amino acid sequence is provided, AND USER DOES NOW WANTS FRAME TO BE PREDICTED (what an a-hole)
                # we will do the dumb thing, and try to translate the DNA sequence directly
                try:
                    amino_acid_sequence = utils.get_translated_sequence_for_gene_call(sequence, gene_callers_id)
                except ConfigError as non_divisible_by_3_error:
                    raise ConfigError(non_divisible_by_3_error.e + ". Since you are creating a contigs database, "
                                      "anvi'o is willing to strike you a deal, but it will require you to trust her a bit more and give her "
                                      "the power to modify the external gene calls you provided. In your external gene calls file you "
                                      "have at least one gene call for which you did not provide an amino acid sequence, and marked it as a "
                                      "CODING type gene call. But becasue YOU ELECTED TO SKIP anvi'o frame prediction to estimate a proper amino "
                                      "acid sequence, anvi'o simply tried to translate your DNA sequence from the start to the end. But as you "
                                      "can tell, it didn't go well. This may be happening because you simply didn't follow the instructions for "
                                      "external gene calls file format. We hope it is not the case as we have described the format of this file "
                                      "here: http://merenlab.org/software/anvio/help/artifacts/external-gene-calls/. But this may also happen "
                                      "even if you have follwed it carefully, but your amino acid sequences are simply not translatable because "
                                      "they are partial. In these cases you have two options. Either you provide amino acid sequences for these "
                                      "gene calls explicitly, or do not use the `--skip-predict-frame` so anvi'o can do the following whenever "
                                      "this problem arises using a Markov model: (1) translate all 3 possible amino acid sequences for the "
                                      "gene (one for each frame), (2) determine which is the most likely based on the tendency that amino acids "
                                      "tend to co-occur as neighbors [nerd speak: a 4th order markov state model trained on the uniprot50 dataset], "
                                      "and finally (3) trim the start and/or stop of your gene to match the most likley frame. The trimming of your "
                                      "start/stop positions will be reflected in the anvi'o contigs database, but will *not* be changed in the external "
                                      "gene calls file you've provided (if you want the modified gene calls as they will appear in your contigs database, "
                                      "you can use `anvi-export-gene-calls` after your contigs database has been created). If all this sounds good to you, "
                                      "go ahead and remove the --skip-predict-frame flag. If not, well then you are on your own :( Find more info here: "
                                      "http://merenlab.org/software/anvio/help/programs/anvi-gen-contigs-database/")

            else:
                raise ConfigError("You broke anvi'o and ended up somewhere no one should ever end up in its codebase. Not nice.")

            # when we are here, we one way or another recovered amino acid sequences either by predicting them or by relying
            # upon user provided data. we have one last control before moving on with our lives:
            if amino_acid_sequence.find('*') > -1:
                if ignore_internal_stop_codons:
                    amino_acid_sequence = amino_acid_sequence.replace('*', 'X')
                    report[P("num_%sgenes_with_internal_stops")] += 1
                else:
                    os.remove(self.db_path)
                    raise ConfigError("Oops. Anvi'o run into an amino acid sequence (that corresponds to the gene callers id '%s') "
                                      "which had an internal stop codon :/ This is sometimes due to errors in the external gene "
                                      "calls file, but more often it is due to the non-standard genetic code. You still can continue "
                                      "by asking anvi'o to ignore internal stop codons via the flag `--ignore-internal-stop-codons`. "
                                      "It will probably look very ugly on your screen, but here is the "
                                      "DNA sequence for that gene in case you don't trust anvi'o (which only would be fair since "
                                      "anvi'o does not trust you either): '%s'. And here is the amino acid sequence of the "
                                      "same gene call if you would like to BLAST it around and see whether it actually makes "
                                      "sense as a gene call: '%s'." % (str(gene_callers_id), sequence, amino_acid_sequence))

            amino_acid_sequences[gene_callers_id] = amino_acid_sequence

        # reporting time
        self.run.warning(None, header="EXTERNAL GENE CALLS PARSER REPORT", lc="cyan")
        self.run.info("Num gene calls in file", len(gene_calls_dict))
        self.run.info("Non-coding gene calls", report["num_non_coding_gene_calls"])
        self.run.info("Partial gene calls", report["num_partial_gene_calls"])
        self.run.info("Num amino acid sequences provided", report["num_partial_gene_calls_with_user_provided_amino_acid_sequences"] + report["num_gene_calls_with_user_provided_amino_acid_sequences"], mc="green")
        self.run.info("  - For complete gene calls", report["num_gene_calls_with_user_provided_amino_acid_sequences"])
        self.run.info("  - For partial gene calls", report["num_partial_gene_calls_with_user_provided_amino_acid_sequences"])
        self.run.info("Frames predicted", report["num_partial_gene_calls_with_predicted_frame"] + report["num_gene_calls_with_predicted_frame"])
        self.run.info("  - For complete gene calls", report["num_gene_calls_with_predicted_frame"])
        self.run.info("  - For partial gene calls", report["num_partial_gene_calls_with_predicted_frame"])
        self.run.info("Gene calls marked as NONCODING", report["num_partial_gene_calls_with_no_predicted_frame"] + report["num_gene_calls_with_no_predicted_frame"], mc="red")
        self.run.info("  - For complete gene calls", report["num_gene_calls_with_no_predicted_frame"], mc="red")
        self.run.info("  - For partial gene calls", report["num_partial_gene_calls_with_no_predicted_frame"], mc="red")
        self.run.info("Gene calls with internal stops", report["num_genes_with_internal_stops"] + report["num_partial_genes_with_internal_stops"])
        self.run.info("  - For complete gene calls", report["num_genes_with_internal_stops"])
        self.run.info("  - For partial gene calls", report["num_partial_genes_with_internal_stops"], nl_after=1)

        return gene_calls_dict, amino_acid_sequences


    def update_gene_call(self, gene_call, frame):
        """Updates gene call that has had its codon frame predicted"""

        reverse = True if gene_call['direction'] == 'r' else False

        trim_5p = frame
        trim_3p = (gene_call['stop'] - gene_call['start'] - frame) % 3

        if reverse:
            gene_call['start'] += trim_3p
            gene_call['stop'] -= trim_5p
        else:
            gene_call['start'] += trim_5p
            gene_call['stop'] -= trim_3p

        return gene_call


    def call_genes_and_populate_genes_in_contigs_table(self, gene_caller='prodigal'):
        Table.__init__(self, self.db_path, anvio.__contigs__version__, self.run, self.progress, simple=True)

        # get gene calls and amino acid sequences
        gene_calls_dict, amino_acid_sequences = self.run_gene_caller(gene_caller)

        # make sure the returning gene calls dict is proper
        self.check_gene_calls_dict(gene_calls_dict)

        # populate genes_in_contigs, and gene_amino_acid_sequences table in contigs db.
        self.populate_genes_in_contigs_table(gene_calls_dict, amino_acid_sequences)


    def run_gene_caller(self, gene_caller='prodigal'):
        """Runs gene caller, and returns gene_calls_dict, and amino acid sequences."""
        remove_fasta_after_processing = False

        if not self.contigs_fasta:
            self.contigs_fasta = filesnpaths.get_temp_file_path()
            utils.export_sequences_from_contigs_db(self.db_path,
                                                   output_file_path=self.contigs_fasta,
                                                   run=self.run)
            remove_fasta_after_processing = True

        if self.debug:
            self.run.info_single('--debug flag is [ON], which means temporary directories generated by '
                                'this run will not be removed', nl_after=2)

        gene_caller = genecalling.GeneCaller(self.contigs_fasta, gene_caller=gene_caller, args=self.args, debug=self.debug)

        gene_calls_dict, amino_acid_sequences = gene_caller.process()

        if not self.debug and remove_fasta_after_processing:
            os.remove(self.contigs_fasta)

        return gene_calls_dict, amino_acid_sequences


    def populate_genes_in_contigs_table(self, gene_calls_dict, amino_acid_sequences, append_to_the_db=False):
        utils.is_contigs_db(self.db_path)
        database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))

        if not append_to_the_db:
            database._exec('''DELETE FROM %s''' % (t.genes_in_contigs_table_name))
            database._exec('''DELETE FROM %s''' % (t.gene_amino_acid_sequences_table_name))
        else:
            # so we are in the append mode. We must remove all the previous entries from genes in contigs
            # that matches to the incoming sources. otherwise we may end up with many duplicates in the db.
            sources = set([v['source'] for v in gene_calls_dict.values()])

            # basically here we will go through those sources, find gene caller ids associated with them in
            # the genes in contigs table, and then remove entries for those gene caller ids both from the
            # genes in contigs and genes in splits tables.
            for source in sources:
                gene_caller_ids_for_source = database.get_single_column_from_table(t.genes_in_contigs_table_name,
                                                                                   'gene_callers_id',
                                                                                   where_clause="""source='%s'""" % source)

                if gene_caller_ids_for_source:
                    for table_name in [t.genes_in_contigs_table_name, t.genes_in_splits_table_name]:
                        database._exec('''DELETE FROM %s WHERE gene_callers_id IN (%s)''' % \
                                                    (table_name, ','.join([str(g) for g in gene_caller_ids_for_source])))

        self.progress.new('Processing')
        self.progress.update('Entering %d gene calls into the db ...' % (len(gene_calls_dict)))

        db_entries = [tuple([gene_callers_id] + [gene_calls_dict[gene_callers_id][h] for h in t.genes_in_contigs_table_structure[1:]]) for gene_callers_id in gene_calls_dict]
        database._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?,?,?,?,?)''' % t.genes_in_contigs_table_name, db_entries)

        db_entries = [tuple([gene_callers_id, amino_acid_sequences[gene_callers_id] if gene_callers_id in amino_acid_sequences else '']) for gene_callers_id in gene_calls_dict]
        database._exec_many('''INSERT INTO %s VALUES (?,?)''' % t.gene_amino_acid_sequences_table_name, db_entries)

        self.progress.end()

        database.disconnect()


    def populate_genes_in_splits_tables(self, gene_calls_dict=None):
        utils.is_contigs_db(self.db_path)
        Table.__init__(self, self.db_path, anvio.__contigs__version__, self.run, self.progress)
        self.init_gene_calls_dict()

        if not gene_calls_dict:
            gene_calls_dict = self.gene_calls_dict

        genes_in_splits = GenesInSplits()
        # build a dictionary for fast access to all genes identified within a contig
        gene_calls_in_contigs_dict = {}
        for gene_callers_id in gene_calls_dict:
            contig = gene_calls_dict[gene_callers_id]['contig']
            if contig in gene_calls_in_contigs_dict:
                gene_calls_in_contigs_dict[contig].add(gene_callers_id)
            else:
                gene_calls_in_contigs_dict[contig] = set([gene_callers_id])

        contigs_without_any_gene_calls = list(set(self.contigs_info.keys()) - set(gene_calls_in_contigs_dict.keys()))
        self.run.info('Contigs with at least one gene call', '%d of %d (%.1f%%)' % (len(gene_calls_in_contigs_dict),
                                                                                    len(self.contigs_info),
                                                                                    len(gene_calls_in_contigs_dict) * 100.0 / len(self.contigs_info)))

        for contig in contigs_without_any_gene_calls:
            gene_calls_in_contigs_dict[contig] = set([])

        splits_dict = {}
        for contig in self.contigs_info:
            for split_name in self.contig_name_to_splits[contig]:
                start = self.splits_info[split_name]['start']
                stop = self.splits_info[split_name]['end']

                gene_start_stops = []
                # here we go through all genes in the contig and identify the all the ones that happen to be in
                # this particular split to generate summarized info for each split. BUT one important that is done
                # in the following loop is genes_in_splits.add call, which populates GenesInSplits class.
                for gene_callers_id in gene_calls_in_contigs_dict[contig]:
                    if gene_calls_dict[gene_callers_id]['stop'] > start and gene_calls_dict[gene_callers_id]['start'] < stop:
                        gene_start_stops.append((gene_calls_dict[gene_callers_id]['start'], gene_calls_dict[gene_callers_id]['stop']), )
                        genes_in_splits.add(split_name, start, stop, gene_callers_id, gene_calls_dict[gene_callers_id]['start'], gene_calls_dict[gene_callers_id]['stop'])

                # here we identify genes that are associated with a split even if one base of the gene spills into
                # the defined start or stop of a split, which means, split N, will include genes A, B and C in this
                # scenario:
                #
                # contig: (...)------[ gene A ]--------[     gene B    ]----[gene C]---------[    gene D    ]-----(...)
                #         (...)----------x---------------------------------------x--------------------------------(...)
                #                        ^ (split N start)                       ^ (split N stop)
                #                        |                                       |
                #                        |<-              split N              ->|
                #
                # however, when looking at the coding versus non-coding nucleotide ratios in a split, we have to make
                # sure that only the relevant portion of gene A and gene C is counted:
                total_coding_nts = 0
                for gene_start, gene_stop in gene_start_stops:
                    total_coding_nts += (gene_stop if gene_stop < stop else stop) - (gene_start if gene_start > start else start)

                splits_dict[split_name] = {'num_genes': len(gene_start_stops),
                                           'avg_gene_length': numpy.mean([(l[1] - l[0]) for l in gene_start_stops]) if len(gene_start_stops) else 0.0,
                                           'ratio_coding': total_coding_nts * 1.0 / (stop - start),
                                           }

        # open connection
        database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))

        # push entries for genes in splits table
        db_entries = [[d[h] for h in t.genes_in_splits_table_structure] for d in genes_in_splits.d]
        database._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?)''' % t.genes_in_splits_table_name, db_entries)

        # disconnect
        database.disconnect()


class GenesInSplits:
    def __init__(self):
        self.d = []


    def add(self, split_name, split_start, split_end, gene_callers_id, prot_start, prot_end):
        gene_length = prot_end - prot_start

        if gene_length <= 0:
            raise ConfigError("tables/genecalls/GeneInSplits: OK. There is something wrong. We have this gene, '%s',\
                                which starts at position %d and ends at position %d. Well, it doesn't look right,\
                                does it?" % (gene_callers_id, prot_start, prot_end))

        # if only a part of the gene is in the split:
        start_in_split = (split_start if prot_start < split_start else prot_start) - split_start
        stop_in_split = (split_end if prot_end > split_end else prot_end) - split_start
        percentage_in_split = (stop_in_split - start_in_split) * 100.0 / gene_length

        self.d.append({'split': split_name,
                       'gene_callers_id': gene_callers_id,
                       'start_in_split': start_in_split,
                       'stop_in_split': stop_in_split,
                       'percentage_in_split': percentage_in_split})
