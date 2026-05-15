from collections import defaultdict
from glob import glob
import hashlib

import pandas as pd

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.fastalib as u
import anvio.dbops as dbops
import anvio.utils as utils
from anvio.genomestorage import GenomeStorage
from anvio.tables.genefunctions import TableForGeneFunctions


class PanRepresenter:
    def __init__(self, args, tmpdir):
        self.pan = dbops.PanSuperclass(args)
        self.pan.init_gene_clusters()
        self.storage = GenomeStorage(args.genomes_storage)
        self._external_genomes = args.external_genomes
        self.genome_to_clusters = self.pan.get_gene_clusters_in_genomes_dict(
            self.pan.gene_clusters
        )
        self.flat_lookup = {
            (genome, gene_id): cluster
            for genome, cluster_subset in self.genome_to_clusters.items()
            for cluster in cluster_subset
            for gene_id in self.pan.gene_clusters[cluster][genome]
        }
        self.first_iteration = True
        self.current_id = 0
        self.current_pos = 0
        self.all_clusters = set(self.pan.gene_clusters)
        self.seen_clusters = set()
        self.functions = dict()
        self.gene_calls = dict()
        self.best_sequences = dict()
        self.accessory_contig = list()
        self.args = args
        self.tmpdir = tmpdir
        self.suplement_contig_name = ""

    @property
    def external_genomes(self):
        # Read the external-genomes.txt file
        df = pd.read_csv(self._external_genomes, sep="\t")
        # Turn into a dictionary: {genome_name: contigs_db_path}
        return dict(zip(df["name"], df["contigs_db_path"]))

    def get_most_common_genome(self):
        """
        Calculates the number of gene_clusters each genome is missing and returns the name of the genome missing the least.

        :return: The name of the most common genome
        """

        num_gene_clusters_missing_per_genome = dict(
            [(genome_name, 0) for genome_name in self.storage.genome_names]
        )

        for genome_name in self.storage.genome_names:
            for gene_cluster_name in self.all_clusters:
                if gene_cluster_name not in self.genome_to_clusters[genome_name]:
                    num_gene_clusters_missing_per_genome[genome_name] += 1
        sorted_dict = sorted(
            num_gene_clusters_missing_per_genome.items(),
            key=lambda x: x[1],
            reverse=False,
        )
        return sorted_dict[0][0]

    def get_best_genome(self, alpha=0.8):
        """
        Returns the name of the best genome based on completeness, redundancy and number of contigs.
        It maximizes completeness and minimizes redundancy and number of contigs.

        :param alpha: The weigth given to completeness
        :return: The name of the best genome
        """
        if self.args.best_genome:
            return self.args.best_genome

        genomes_dict = self.storage.get_genomes_dict()
        genomes = [
            {
                "name": name,
                "percent_completion": data["percent_completion"],
                "percent_redundancy": data["percent_redundancy"],
                "num_contigs": data["num_contigs"],
                "num_genes": data.get("num_genes"),
            }
            for name, data in genomes_dict.items()
            if data.get("num_contigs") < 5
        ]
        # Extract values for normalization
        comp_values = [
            g["percent_completion"] - g["percent_redundancy"] for g in genomes
        ]
        contig_values = [g["num_contigs"] for g in genomes]

        comp_min, comp_max = min(comp_values), max(comp_values)
        contig_min, contig_max = min(contig_values), max(contig_values)

        # Unified normalization
        def normalize(x, min_val, max_val, maximize=True):
            if max_val == min_val:
                return 0.5
            norm = (x - min_val) / (max_val - min_val)
            return norm if maximize else 1 - norm

        # Compute scores
        for g in genomes:
            comp_norm = normalize(
                g["percent_completion"] - g["percent_redundancy"],
                comp_min,
                comp_max,
                maximize=True,
            )
            contig_norm = normalize(
                g["num_contigs"], contig_min, contig_max, maximize=False
            )
            g["score"] = alpha * comp_norm + (1 - alpha) * contig_norm

        # Return best genome using a tiebreaker
        return max(genomes, key=lambda g: (g["score"], g["num_genes"]))["name"]

    def get_contigs_db(self, genome):
        """Create ContigsSuperclass and DB insentences

        Args:
            genome (str): the name of the genome

        Returns:
            contigs: An instence of the ContigsSuperclass
            contigs_db: An instence of the DB class
        """
        # Create the contigs object for the genome
        genome_db_path = self.external_genomes[genome]

        # Instantiate a contigs database
        self.args.contigs_db = genome_db_path  # I have to do this
        contigs = dbops.ContigsSuperclass(self.args)
        contigs.init_contig_sequences()  # This line was added to get rid of the annoying warning message, remove it you will see what i'm talking about
        # This one allows SQL queries on any database
        contigs_db = db.DB(genome_db_path, anvio.__contigs__version__)
        return contigs, contigs_db

    def build_contig_to_genes_dict(self, genome, contigs):
        """Maps the contigs of the genome to genes

        Args:
            genome (str): the name of the genome
            contigs (Contig): An instence of the ContigsSuperclass

        Returns:
            dict: A dictionary where keys are contigs and values are lists of gene_ids sorted.
        """
        contig_to_genes = defaultdict(list)
        for gene_id, gene_info in contigs.genes_in_contigs_dict.items():
            if gene_id in self.storage.get_gene_caller_ids(genome):
                contig_name = gene_info.get("contig")
                contig_to_genes[contig_name].append(gene_id)
        # Sort gene IDs in genomic order if start positions are available
        for contig, genes in contig_to_genes.items():
            contig_to_genes[contig] = sorted(
                genes,
                key=lambda gid: contigs.genes_in_contigs_dict[gid].get("start", 0),
            )
        return dict(contig_to_genes)

    def build_functions_lookup(self, contigs_db):
        # Create functions Lookup table
        functions_dict = contigs_db.get_table_as_dict(
            t.gene_function_calls_table_name, string_the_key=True
        )
        filtered_funcs = defaultdict(list)
        for fun in functions_dict.values():
            if (
                gene_id := fun.get("gene_callers_id")
            ) is not None:  # Walrus operator (Python 3.8+)
                filtered_funcs[gene_id].append(fun)
        return filtered_funcs

    def set_gene_calls(self, raw_data):
        raw_df = pd.DataFrame(raw_data).T
        raw_df = raw_df.drop(columns=["length", "rev_compd", "sequence", "header"])
        raw_df["gene_callers_id"] = raw_df.index
        raw_df = raw_df[[raw_df.columns[-1]] + list(raw_df.columns[:-1])]
        self.gene_calls = raw_df.T.to_dict()
        return True

    def get_gene_calls_data(self, contigs, target=None):
        result_2 = []
        if not self.first_iteration and target:
            for gene_ids in target.values():
                result_2.extend(gene_ids)

        return contigs.get_sequences_for_gene_callers_ids(
            result_2, report_aa_sequences=True, include_aa_sequences=True
        )[1]

    def get_filtered_genes(self, genome, contig_name_to_genes):
        result = {}
        for contig, gene_ids in contig_name_to_genes.items():
            kept_genes = []
            for gene_id in gene_ids:
                cluster = self.flat_lookup.get((genome, gene_id))
                if cluster and cluster not in self.seen_clusters:
                    kept_genes.append(gene_id)
                    self.seen_clusters.add(cluster)
            result[contig] = kept_genes
        return result

    def add_function(self, genome, filtered_funcs, source, id, contig=None):
        if self.first_iteration:
            contig = source.get("contig")

        if not source.get("version") == "unknown":
            current_funcs = filtered_funcs.get(id, [])

            for func in current_funcs:
                if not self.first_iteration:
                    func["gene_callers_id"] = self.current_id

                self.functions[len(self.functions) + 1] = func

            self.functions[len(self.functions) + 1] = {
                "gene_callers_id": id if self.first_iteration else self.current_id,
                "source": "MedBel",
                "accession": self.generate_accession(source.get("aa_sequence"), 7),
                "function": f"{genome}__{contig}__{id}",
                "e_value": 0,
            }

    def add_gene_call(self, source, start=None, stop=None):
        gene_info = {
            "gene_callers_id": self.current_id,
            "contig": self.suplement_contig_name,
            "start": start if start else self.current_pos,
            "stop": stop if stop else self.current_pos + source.get("length"),
            "direction": source.get("direction"),
            "partial": source.get("partial"),
            "call_type": source.get("call_type"),
            "source": source.get("source"),
            "version": source.get("version"),
            "aa_sequence": source.get("aa_sequence"),
        }
        self.gene_calls[self.current_id] = gene_info
        return True

    def get_stretches(self, source):
        """Takes a array of numbers, and turns reports back stretches

        For example, for a `numbers_list` that looks like this:

        >>> [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 15, 16, 17, 18, 99]

        This function will return the following:

        >>> [(2, 12), (15, 18), (99, 99)]
        """
        streches = dict()
        for contig, gene_ids in source.items():
            streches[contig] = utils.get_stretches_for_numbers_list(
                gene_ids, discard_singletons=False
            )
        return streches

    def generate_accession(self, sequence, size):
        """Generate a hash from an aa-sequence

        Args:
            sequence (str): aa-sequence
            size (int): the size to subset from the hash string

        Returns:
            str: a substring of the hash string of length `size`
        """
        return hashlib.sha256(sequence.encode()).hexdigest()[:size].upper()

    def export_gene_calls(self):
        """Exports the gene calls table as a tab delimited file"""

        gene_calls = pd.DataFrame(self.gene_calls).T
        utils.store_dataframe_as_TAB_delimited_file(
            gene_calls, self.args.external_gene_calls
        )

    # This function is not actually need, but I might keep it for debuging purposes
    def export_functions(self):
        """Exports the functions table as a tab delimited file"""

        functions = pd.DataFrame(self.functions).T
        utils.store_dataframe_as_TAB_delimited_file(functions, "functions.tsv")

    def export_fasta(self):
        """Writes the DNA sequnces to a fasta file `output_file_path`

        Args:
            output_file_path (str): the path to the output file

        Returns:
            bool: retunrs true if successful
        """
        output_fasta = u.FastaOutput(self.args.contigs_fasta)
        for name, seq in self.best_sequences.items():
            output_fasta.write_id(name)
            output_fasta.write_seq(seq)
        return True

    def process_best_genome(self):
        best_genome = self.get_best_genome()
        contigs, contigs_db = self.get_contigs_db(best_genome)
        contig_sequences = contigs_db.get_table_as_dict(
            t.contig_sequences_table_name, string_the_key=True
        )
        filtered_funcs = self.build_functions_lookup(contigs_db)
        raw_data = self.get_gene_calls_data(contigs)

        for k, v in raw_data.items():
            self.add_function(best_genome, filtered_funcs, v, k)

        self.set_gene_calls(raw_data)
        self.best_sequences = {
            k: v.get("sequence") for k, v in contig_sequences.items()
        }
        self.seen_clusters.update(self.genome_to_clusters[best_genome])
        self.current_id = max(contigs.genes_in_contigs_dict) + 1
        self.all_clusters -= self.seen_clusters
        self.first_iteration = False
        # Set the suplement contig's name
        self.suplement_contig_name = f"{best_genome}_pangenome_supplement"

    def process_additional_genomes(self):
        most_common_genome = self.get_most_common_genome()
        contigs, contigs_db = self.get_contigs_db(most_common_genome)
        contig_name_to_genes = self.build_contig_to_genes_dict(
            most_common_genome, contigs
        )
        filtered_funcs = self.build_functions_lookup(contigs_db)
        filtered_genes = self.get_filtered_genes(
            most_common_genome, contig_name_to_genes
        )
        if self.args.keep_promoter:
            raw_data = self.get_gene_calls_data(contigs)
        else:
            raw_data = self.get_gene_calls_data(contigs, filtered_genes)

        streches = self.get_stretches(filtered_genes)
        contig_to_seq = contigs_db.get_table_as_dict(
            t.contig_sequences_table_name, string_the_key=True
        )

        if not self.args.keep_senteny:
            for k, v in raw_data.items():
                self.accessory_contig.append(v.get("sequence"))
                self.add_gene_call(v)
                self.add_function(most_common_genome, filtered_funcs, v, k)
                self.current_id += 1
                self.current_pos += v.get("length") + 20

            self.best_sequences[self.suplement_contig_name] = (
                "NNNNNNNNNNNNNNNNNNNN".join(self.accessory_contig)
            )
            self.all_clusters -= self.seen_clusters
        else:
            for contig, limits in streches.items():
                if not limits:
                    continue

                contig_seq = contig_to_seq.get(contig).get("sequence")
                for start, stop in limits:
                    first = raw_data.get(start)
                    last = raw_data.get(stop)

                    begin = first.get("start")
                    end = last.get("stop")
                    if self.args.keep_promoter:
                        valid_ids = contig_name_to_genes.get(contig)

                        before = (
                            raw_data.get(start - 1) if start > min(valid_ids) else None
                        )
                        after = (
                            raw_data.get(stop + 1) if stop < max(valid_ids) else None
                        )

                        begin = before.get("stop") if before else 0
                        end = after.get("start") if after else len(contig_seq)

                    extracted = contig_seq[begin:end]
                    self.accessory_contig.append(extracted)

                    for i in range(start, stop + 1):
                        current = raw_data.get(i)
                        new_start = current.get("start") - begin + self.current_pos
                        new_stop = current.get("stop") - begin + self.current_pos

                        self.add_gene_call(current, new_start, new_stop)
                        self.add_function(
                            most_common_genome, filtered_funcs, current, i, contig
                        )
                        self.current_id += 1

                    self.current_pos += end - begin + 20

    def assemble_sumplementary_contig(self, size=20):
        gap = "N" * size
        self.best_sequences[self.suplement_contig_name] = gap.join(
            self.accessory_contig
        )

    def write_outputs(self):
        self.export_gene_calls()
        #  self.export_functions()

        # self.export_fasta(f"{dir_path}/My_Perfect_FASTA.fasta")
        self.export_fasta()

    def build_contigs_db(self):
        # Create the contigs_db
        a = dbops.ContigsDatabase(self.args.output_file, quiet=False, skip_init=True)
        a.create(self.args)

        gene_function_calls_table = TableForGeneFunctions(self.args.output_file)
        gene_function_calls_table.create(
            self.functions, drop_previous_annotations_first=False
        )

    def main(self):
        while self.all_clusters:
            if self.first_iteration:
                self.process_best_genome()

            self.process_additional_genomes()
            self.all_clusters -= self.seen_clusters
        # end while

        self.assemble_sumplementary_contig()
        self.write_outputs()
        self.build_contigs_db()
