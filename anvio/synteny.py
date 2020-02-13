import anvio
import anvio.filesnpaths as filesnpaths
import anvio.terminal as terminal
import anvio.utils as utils
import anvio.dbops as dbops
import anvio.tables as t


class NGram(object):

    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.external_genomes = A('external_genomes')
        self.genes = {}
        self.populate_genes()


    def populate_genes(self):
        filepaths = utils.get_TAB_delimited_file_as_dictionary(self.external_genomes)

        for locus_name in filepaths:
            path = filepaths[locus_name]["contigs_db_path"]

            self.genes[locus_name] = self.get_genes_from_contigs_db_path(path)


    def get_genes_from_contigs_db_path(self, path):
        contigs_db = dbops.ContigsDatabase(path)
        annotations_dict = contigs_db.db.get_table_as_dict(t.gene_function_calls_table_name)
        print(annotations_dict)

    def countSynteny(self, k):

        kFreq = {}
        # Make sliding window of k length across all genes
        for i in range(0, len(genes) - k + 1):
            # extract window
            window = sorted(genes[i:i + k])
            ngram = "_".join(map(str, list(window)))
            # if ngram is not in dictionary add it
            # if it is add + 1
            if ngram in kFreq:
                kFreq[ngram] +=  1
            else:
                kFreq[ngram] = 1

        return kFreq
