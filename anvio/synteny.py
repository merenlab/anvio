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
        self.annotation_sources = A('annotation_sources')
        self.window_size = A('window_size')
        self.genes = {}
        self.populate_genes()

    def populate_genes(self):
        filepaths = utils.get_TAB_delimited_file_as_dictionary(self.external_genomes)

        for locus_name in filepaths:
            path = filepaths[locus_name]["contigs_db_path"]
            self.genes[locus_name] = self.get_genes_from_contigs_db_path(path,self.annotation_sources)
        l = []
        l_dict = {}
        for path, d in self.genes.items():
            for k, e in d.items():
                if path not in l_dict:
                    l_dict[path] = (k, e['function'])
                else:
                    pass
                l.append((path, k, e['function']))
        print(l_dict)

    def get_genes_from_contigs_db_path(self, path, annotation_sources):
        contigs_db = dbops.ContigsDatabase(path)
        annotations_dict = contigs_db.db.get_table_as_dict(t.gene_function_calls_table_name)
        requested_sources = [s.strip() for s in annotation_sources.split(',')]
        missing_sources = [s for s in requested_sources if s not in annotation_sources]
        if len(missing_sources):
            raise ConfigError("One or more of the annotation sources you requested does not appear to be in the\
                                contigs database :/ Here is the list: %s." % (', '.join(missing_sources)))

        annotations_dict = utils.get_filtered_dict(annotations_dict, 'source', set(requested_sources))
        return annotations_dict

    k = 3
    def countSynteny(self, k, genes):
        kFreq = {}
        # Make sliding window of k length across all gene
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
    

    def driveSynteny(self):
        for key,value in self.genes:
            self.genes[key] = countSynteny(value)
