import shutil
import hashlib
import argparse
from collections import defaultdict

import pandas as pd

import anvio
import anvio.tables as t
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.fastalib as fastalib
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.genomedescriptions import GenomeDescriptions
from anvio.tables.genefunctions import TableForGeneFunctions


class PanRepresenter:
    def __init__(self, args, r, p):
        self.args = args
        self.run = r
        self.progress = p

        self.rq = terminal.Run(verbose=False)
        self.pq = terminal.Progress(verbose=False)

        # figuring out the args
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.output_file = A('output_file')
        self.keep_promoter = A('keep_promoter')
        self.keep_synteny = A('keep_synteny') or A('keep_promoter')
        self.gap_size = A('gap_size')
        self.alpha = A('alpha')
        self.max_num_contigs = A('max_num_contigs')
        self.representative = A('representative')

        # Check if output file is writeable .. we do it here rather than in sanity_check
        # because we have to initialize the PanSuperclass before sanity check, which takes
        # a long time and we don't want the user to wait for that before learning their
        # output path is bad or output file already exists.
        filesnpaths.is_output_file_writable(self.output_file, ok_if_exists=False)

        # initialize stuff!
        self.progress.new("Bleep bloop")
        self.progress.update('Initializing the pangenome, genome descriptions, and gene clusters ..')
        self.pan = dbops.PanSuperclass(args, r=self.rq, p=self.pq)
        self.pan.init_gene_clusters()
        self.project_name = A('project_name') or self.pan.p_meta.get('project_name')

        self.progress.end()

        # SANITY CHECK
        self.sanity_check()

        self.external_genomes = GenomeDescriptions(self.args)
        self.external_genomes.load_genomes_descriptions()
        self.genome_to_clusters = self.pan.get_gene_clusters_in_genomes_dict(self.pan.gene_clusters)

        # some default variables that will be filled in by the class methods
        self.first_iteration = True
        self.current_id = 0
        self.current_pos = 0
        self.all_clusters = set(self.pan.gene_clusters)
        self.seen_clusters = set()
        self.functions = {}
        self.gene_calls = {}
        self.representative_contigs = {}
        self.supplementary_contig = []
        self.supplement_contig_name = ""
        self.flat_lookup = {}
        for genome, cluster_subset in self.genome_to_clusters.items():
            for cluster in cluster_subset:
                for gene_id in self.pan.gene_clusters[cluster][genome]:
                    self.flat_lookup[(genome, gene_id)] = cluster

    def sanity_check(self):
        if self.representative and self.representative not in self.pan.genomes_storage.genome_names:
            raise ConfigError(f"The genome '{self.representative}' was provided as the representative "
                            f"but does not exist in the genomes storage.")

        if self.keep_promoter:
            self.run.warning("Since you chose to keep the promoter region that means you also keep the synteny by definition")

        if self.max_num_contigs <= 0:
            raise ConfigError(f"--max-num-contigs must be a positive integer, but you provided {self.max_num_contigs}. "
                            f"Please try again with a value greater than 0.")

        if self.gap_size <= 0:
            raise ConfigError(f"--gap-size must be a positive integer, but you provided {self.gap_size}. "
                            f"Please try again with a value greater than 0.")
        if self.alpha < 0 or self.alpha > 1:
            raise ConfigError(f"--alpha must be a value between 0.0 and 1.0 inclusive, but you provided {self.alpha}.")


    def get_most_common_genome(self):
        """Calculates the number of gene_clusters each genome is missing and returns the name of the genome missing the least.

        Returns
        =======
        str
            The name of the most common genome.
        """

        num_gene_clusters_missing_per_genome = {genome_name: 0 for genome_name in self.pan.genomes_storage.genome_names}

        for genome_name in self.pan.genomes_storage.genome_names:
            for gene_cluster_name in self.all_clusters:
                if gene_cluster_name not in self.genome_to_clusters[genome_name]:
                    num_gene_clusters_missing_per_genome[genome_name] += 1
        sorted_genomes = sorted(
            num_gene_clusters_missing_per_genome.items(),
            key=lambda x: x[1],
        )
        return sorted_genomes[0][0]


    def get_representative(self):
        """Returns the name of the representative genome based on completeness, redundancy, and number of contigs.
        It maximizes completeness and minimizes redundancy and number of contigs.
        In case of a tie, it returns the genome with the highest number of genes.

        Returns
        =======
        str
            The name of the representative genome.
        """

        if self.representative:
            return self.representative

        genomes_dict = self.pan.genomes_storage.get_genomes_dict()
        genomes = [
            {
                "name": name,
                "percent_completion": data["percent_completion"],
                "percent_redundancy": data["percent_redundancy"],
                "num_contigs": data["num_contigs"],
                "num_genes": data.get("num_genes"),
            }
            for name, data in genomes_dict.items()
            if data.get("num_contigs") < self.max_num_contigs
        ]
        if not genomes:
            raise ConfigError(f"I'm sorry to tell you, but no genome has less than {self.max_num_contigs} contigs, maybe try to increase the limit")

        comp_values = [
            g["percent_completion"] - g["percent_redundancy"] for g in genomes
        ]
        num_contigs = [g["num_contigs"] for g in genomes]

        comp_min, comp_max = min(comp_values), max(comp_values)
        contig_min, contig_max = min(num_contigs), max(num_contigs)

        def normalize(x, min_val, max_val, maximize=True):
            if max_val == min_val:
                return 0.5
            norm = (x - min_val) / (max_val - min_val)
            return norm if maximize else 1 - norm

        for g in genomes:
            comp_norm = normalize(
                g["percent_completion"] - g["percent_redundancy"],
                comp_min,
                comp_max,
            )
            contig_norm = normalize(g["num_contigs"], contig_min, contig_max, maximize=False)
            g["score"] = self.alpha * comp_norm + (1 - self.alpha) * contig_norm

        return max(genomes, key=lambda g: (g["score"], g["num_genes"]))["name"]


    def get_contigs_db(self, genome):
        """Create ContigsSuperclass and ContigsDatabase instances

        Parameters
        ==========
        genome : str
            the name of the genome

        Returns
        =======
        contigs : ContigsSuperclass
            An instance of the ContigsSuperclass
        contigs_db : ContigsDatabase
            An instance of the ContigsDatabase class
        """
        genome_db_path = self.external_genomes.genomes.get(genome).get("contigs_db_path")
        args = argparse.Namespace(**vars(self.args), contigs_db=genome_db_path)
        contigs = dbops.ContigsSuperclass(args, r=self.rq, p=self.pq)
        contigs.init_contig_sequences()
        contigs_db = dbops.ContigsDatabase(genome_db_path, run=self.rq, progress=self.pq, quiet=True)
        return contigs, contigs_db


    def build_contig_to_genes_dict(self, genome, contigs):
        """Maps the contigs of the genome to genes

        Parameters
        ==========
        genome : str
            the name of the genome
        contigs : ContigsSuperclass
            An instance of the ContigsSuperclass

        Returns
        =======
        dict
            A dictionary where keys are contig names and values are lists of gene_ids sorted by genomic order.
        """
        contig_to_genes = defaultdict(list)
        valid_gene_ids = set(self.pan.genomes_storage.get_gene_caller_ids(genome))
        for gene_id, gene_info in contigs.genes_in_contigs_dict.items():
            if gene_id in valid_gene_ids:
                contig_name = gene_info["contig"]
                contig_to_genes[contig_name].append(gene_id)
        for contig, genes in contig_to_genes.items():
            contig_to_genes[contig] = sorted(
                genes,
                key=lambda gid: contigs.genes_in_contigs_dict[gid]["start"],
            )
        return dict(contig_to_genes)


    def build_functions_lookup(self, contigs_db):
        functions_dict = contigs_db.db.get_table_as_dict(
            t.gene_function_calls_table_name, string_the_key=True
        )
        filtered_funcs = defaultdict(list)
        for fun in functions_dict.values():
            if (gene_id := fun.get("gene_callers_id")) is not None:
                filtered_funcs[gene_id].append(fun)
        return filtered_funcs


    def set_gene_calls(self, gene_calls_data):
        raw_df = pd.DataFrame(gene_calls_data).T
        raw_df = raw_df.drop(columns=["length", "rev_compd", "sequence", "header"])
        raw_df["gene_callers_id"] = raw_df.index
        raw_df = raw_df[[raw_df.columns[-1]] + list(raw_df.columns[:-1])]
        self.gene_calls = raw_df.T.to_dict()


    def get_gene_calls_data(self, contigs, target=None):
        if target is None:
            gene_ids_of_interest = list(contigs.genes_in_contigs_dict.keys())
        else:
            gene_ids_of_interest = []
            for gene_ids in target.values():
                gene_ids_of_interest.extend(gene_ids)

        return contigs.get_sequences_for_gene_callers_ids(
            gene_ids_of_interest, report_aa_sequences=True, include_aa_sequences=True
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


    def add_function(self, genome, filtered_funcs, source, gene_id, contig=None):
        if self.first_iteration:
            contig = source.get("contig")

        if source["version"] != "unknown":
            current_funcs = filtered_funcs.get(gene_id, [])

            for func in current_funcs:
                if not self.first_iteration:
                    func["gene_callers_id"] = self.current_id

                self.functions[len(self.functions) + 1] = func

            self.functions[len(self.functions) + 1] = {
                "gene_callers_id": gene_id if self.first_iteration else self.current_id,
                "source": "PanRepresentative",
                "accession": self.generate_accession(source.get("aa_sequence"), 7),
                "function": f"Source genome: {genome}; Source contig: {contig}; Gene id in source genome: {gene_id}{cluster_part}",
                "e_value": 0,
            }


    def add_gene_call(self, source, start=None, stop=None):
        gene_info = {
            "gene_callers_id": self.current_id,
            "contig": self.supplement_contig_name,
            "start": start if start is not None else self.current_pos,
            "stop": stop if stop is not None else self.current_pos + source["length"],
            "direction": source["direction"],
            "partial": source["partial"],
            "call_type": source["call_type"],
            "source": source["source"],
            "version": source["version"],
            "aa_sequence": source.get("aa_sequence"),
        }
        self.gene_calls[self.current_id] = gene_info


    def get_stretches(self, source):
        """For each contig, returns the consecutive stretches of gene IDs.

        Parameters
        ==========
        source : dict
            A dictionary mapping contig names to lists of gene IDs.

        Returns
        =======
        stretches : dict
            A dictionary mapping contig names to lists of (first_id, last_id) tuples,
            where each tuple represents a consecutive stretch of gene IDs.
        """
        stretches = {}
        for contig, gene_ids in source.items():
            stretches[contig] = utils.get_stretches_for_numbers_list(gene_ids, discard_singletons=False)
        return stretches


    def generate_accession(self, sequence, size):
        """Generate a hash from an aa-sequence

        Parameters
        ==========
        sequence : str
            aa-sequence
        size : int
            the size to subset from the hash string

        Returns
        =======
        str
            A substring of the hash string of length `size`.
        """
        if not sequence:
            return ""
        return hashlib.sha256(sequence.encode()).hexdigest()[:size].upper()


    def export_gene_calls(self, tmp_dir):
        """Exports the gene calls table as a tab delimited file"""

        external_gene_calls = f"{tmp_dir}/{self.project_name}_gene_calls.txt"
        gene_calls = pd.DataFrame(self.gene_calls).T
        utils.store_dataframe_as_TAB_delimited_file(gene_calls, external_gene_calls)
        self.args.external_gene_calls = external_gene_calls


    def export_functions(self, tmp_dir):
        """Exports the functions table as a tab delimited file"""

        functions = pd.DataFrame(self.functions).T
        utils.store_dataframe_as_TAB_delimited_file(functions, f"{tmp_dir}/{self.project_name}_functions.tsv")


    def export_fasta(self, tmp_dir):
        self.args.contigs_fasta = f"{tmp_dir}/{self.project_name}.fasta"
        output_fasta = fastalib.FastaOutput(self.args.contigs_fasta)
        for name, seq in self.representative_contigs.items():
            output_fasta.write_id(name)
            output_fasta.write_seq(seq)


    def process_representative_genome(self):
        self.representative = self.get_representative()
        self.run.info('Representative genome', self.representative, nl_before=1, mc='green')
        contigs, contigs_db = self.get_contigs_db(self.representative)
        filtered_funcs = self.build_functions_lookup(contigs_db)
        gene_calls_data = self.get_gene_calls_data(contigs)

        for gene_id, gene_data in gene_calls_data.items():
            self.add_function(self.representative, filtered_funcs, gene_data, gene_id)

        self.set_gene_calls(gene_calls_data)
        self.representative_contigs = {contig_name: data["sequence"] for contig_name, data in contigs.contig_sequences.items()}
        self.seen_clusters.update(self.genome_to_clusters[self.representative])
        self.current_id = max(contigs.genes_in_contigs_dict) + 1
        self.all_clusters -= self.seen_clusters
        self.first_iteration = False
        self.supplement_contig_name = f"{self.representative}_pangenome_supplement"
        self.run.info('Gene clusters covered', len(self.seen_clusters))


    def process_additional_genomes(self):
        most_common_genome = self.get_most_common_genome()
        clusters_before = len(self.seen_clusters)
        self.run.info('Supplementing from', most_common_genome, nl_before=1)
        contigs, contigs_db = self.get_contigs_db(most_common_genome)
        contig_name_to_genes = self.build_contig_to_genes_dict(most_common_genome, contigs)
        filtered_funcs = self.build_functions_lookup(contigs_db)
        filtered_genes = self.get_filtered_genes(most_common_genome, contig_name_to_genes)

        if self.keep_promoter:
            gene_calls_data = self.get_gene_calls_data(contigs)
        else:
            gene_calls_data = self.get_gene_calls_data(contigs, filtered_genes)

        stretches = self.get_stretches(filtered_genes)

        if not self.keep_synteny:
            for gene_id, gene_data in gene_calls_data.items():
                self.supplementary_contig.append(gene_data["sequence"])
                self.add_gene_call(gene_data)
                self.add_function(most_common_genome, filtered_funcs, gene_data, gene_id)
                self.current_id += 1
                self.current_pos += gene_data["length"] + self.gap_size

        else:
            for contig, limits in stretches.items():
                if not limits:
                    continue

                contig_seq = contigs.contig_sequences[contig]["sequence"]
                valid_ids = contig_name_to_genes.get(contig)
                min_valid_id, max_valid_id = min(valid_ids), max(valid_ids)
                for first_id, last_id in limits:
                    first_gene_call = gene_calls_data.get(first_id)
                    last_gene_call = gene_calls_data.get(last_id)

                    start = first_gene_call["start"]
                    stop = last_gene_call["stop"]
                    if self.keep_promoter:
                        previous_gene_id = (gene_calls_data.get(first_id - 1) if first_id > min_valid_id else None)
                        next_gene_id = (gene_calls_data.get(last_id + 1) if last_id < max_valid_id else None)

                        start = previous_gene_id["stop"] if previous_gene_id else 0
                        stop = (max(last_gene_call["stop"], next_gene_id["start"]) if next_gene_id else len(contig_seq))

                    extracted = contig_seq[start:stop]
                    self.supplementary_contig.append(extracted)

                    for i in range(first_id, last_id + 1):
                        current = gene_calls_data.get(i)
                        new_start = current["start"] - start + self.current_pos
                        new_stop = current["stop"] - start + self.current_pos

                        self.add_gene_call(current, new_start, new_stop)
                        self.add_function(most_common_genome, filtered_funcs, current, i, contig)
                        self.current_id += 1

                    self.current_pos += stop - start + self.gap_size
        self.all_clusters -= self.seen_clusters
        self.run.info('New clusters added', len(self.seen_clusters) - clusters_before)
        self.run.info('Clusters still missing', len(self.all_clusters))


    def assemble_supplementary_contig(self):
        gap = "N" * self.gap_size
        self.representative_contigs[self.supplement_contig_name] = gap.join(self.supplementary_contig)


    def build_contigs_db(self):
        self.args.db_variant = "pan-genome"

        panrep_db = dbops.ContigsDatabase(self.output_file, quiet=True, skip_init=True)
        panrep_db.create(self.args)

        gene_function_calls_table = TableForGeneFunctions(self.output_file)
        gene_function_calls_table.create(self.functions, drop_previous_annotations_first=False)
        self.run.info('Pan-representative contigs DB', self.output_file, mc='green',nl_before=1, nl_after=1)


    def cleanup(self, tmp_dir):
        if anvio.DEBUG:
            self.run.warning(f"The temp directory, {tmp_dir}, is kept. Please don't forget to clean it up later", header="Debug")
        else:
            self.run.info_single("Cleaning up the temp directory (you can use `--debug` if you would like to keep it for testing purposes)", nl_before=1, nl_after=1)
            shutil.rmtree(tmp_dir)


    def process(self):
        tmp_dir = filesnpaths.get_temp_directory_path()

        while self.all_clusters:
            if self.first_iteration:
                self.process_representative_genome()

            self.process_additional_genomes()

        self.assemble_supplementary_contig()
        self.export_gene_calls(tmp_dir)
        self.export_functions(tmp_dir)
        self.export_fasta(tmp_dir)
        self.build_contigs_db()
        self.cleanup(tmp_dir)
