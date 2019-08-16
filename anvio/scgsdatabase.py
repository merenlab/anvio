#!/usr/bin/env python3
# -*- coding: utf-8
"""
    This file contains SCGsSetup, SCGsDataBase, LowestIdentity classes.

"""
import os
import sys
import gzip
import shutil
import requests
import pickle
import re
import shutil
import subprocess
import tarfile
import multiprocessing

from io import BytesIO

import anvio
import anvio.fastalib as u
import anvio.terminal as terminal
import anvio.pfam as pfam
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.errors import FilesNPathsError

from anvio.drivers.diamond import Diamond


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print
run_quiet = terminal.Run(log_file_path=None,verbose=False)
progress_quiet = terminal.Progress(verbose=False)


_author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Quentin Clayssen"
__email__ = "quentin.clayssen@gmail.com"

def read_remote_file(url, is_gzip=True):
    remote_file = requests.get(url)

    if remote_file.status_code == 404:
        raise Exception("'%s' returned 404 Not Found. " % url)

    if is_gzip:
        buf = BytesIO(remote_file.content)
        fg = gzip.GzipFile(fileobj=buf)
        return fg.read().decode('utf-8')

    return remote_file.content.decode('utf-8')


class SCGsSetup(object):
    def __init__(self, args, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress
        self.SCG_data_dir = args.scgs_data_dir

        if not self.SCG_data_dir:
            self.SCG_data_dir = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/SCG')

        if not args.reset:
            self.is_database_exists()

        filesnpaths.gen_output_directory(self.SCG_data_dir, delete_if_exists=args.reset)

        self.database_url = "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/"
        self.files = ['VERSION', 'ar122_msa_individual_genes.tar.gz', 'ar122_taxonomy.tsv','bac120_msa_individual_genes.tar.gz','bac120_taxonomy.tsv']


    def is_database_exists(self):
        if os.path.exists(os.path.join(self.SCG_data_dir, 'SCG-A.hmm.gz')):
            raise ConfigError("It seems you already have SCG database installed in '%s', please use --reset flag if you want to re-download it." % self.SCG_data_dir)


    def get_remote_version(self):
        content = read_remote_file(self.database_url + 'VERSION',is_gzip=False)

        # below we are parsing this, not so elegant.
        #v89
        #
        #Released June 17th, 2019
        version = content.strip().split('\n')[0].strip()
        release_date = content.strip().split('\n')[2].strip()

        self.run.info("Current gtdb release", "%s (%s)" % (version, release_date))


    def download(self):
        self.run.info("Database URL", self.database_url)

        for file_name in self.files:
            utils.download_file(self.database_url + '/' + file_name,
                os.path.join(self.SCG_data_dir, file_name), progress=self.progress, run=self.run)

        self.confirm_downloaded_files()
        self.decompress_files()
        self.merge_tsv_files()


    def confirm_downloaded_files(self):

        for file_name in self.files:
            if not filesnpaths.is_file_exists(os.path.join(self.SCG_data_dir, file_name), dont_raise=True):
                 # TO DO: Fix messages :(
                raise ConfigError("Have missing file %s, please run --reset" % file_name)

            hash_on_disk = utils.get_file_md5(os.path.join(self.SCG_data_dir, file_name))


    def decompress_files(self):

        file_to_dextract=[file_name for file_name in self.files if file_name.endswith(".tar.gz")]
        for file_name in file_to_dextract:
            full_path = os.path.join(self.SCG_data_dir, file_name)
            tar=tarfile.open(full_path)
            extractpaht=os.path.join(self.SCG_data_dir,"msa_individual_genes")
            self.run.info("Extracting %s" % (file_name), "%s" % (extractpaht))
            tar.extractall(path=extractpaht)
            tar.close()
            #os.remove(full_path)

    def merge_tsv_files(self):
        file_to_dextract=[file_name for file_name in self.files if file_name.endswith(".tsv")]
        full_path = os.path.join(self.SCG_data_dir, "merge_taxonomy.tsv")
        self.run.info("Mergin tsv file ", ','.join(file_to_dextract))
        with open(full_path, 'w') as outfile:
            for fname in file_to_dextract:
                with open(self.SCG_data_dir+"/"+fname) as infile:
                    for line in infile:
                        outfile.write(line)



class SCGsDataBase():

    """Generate diamond databases of GTDB or import's Single Copy core genes Fasta whom matching with anvi'o HMM's"""
    def __init__(self, args, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None

        self.genes_files_directory = A('genes_files_directory')
        self.scgs_directory =A('scgs_directory')
        self.hmms=A('hmms')
        self.path_tsv_taxonomy=A('taxonomy_files')
        self.output_directory=A('out_put_directory')
        self.num_threads=A('num_threads')
        self.taxonomy_level=A('taxonomy_level')


        self.classic_input_directory = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/SCG/')

        if not self.path_tsv_taxonomy:
            self.path_tsv_taxonomy = os.path.join(self.classic_input_directory,'merge_taxonomy.tsv')

        self.num_threads=args.num_threads

        if not self.num_threads:
            self.num_threads=2


        if not self.genes_files_directory:
            self.genes_files_directory = os.path.join(self.classic_input_directory,'msa_individual_genes')


        if not self.output_directory:
             self.output_directory=os.path.join(os.path.dirname(anvio.__file__), 'data/misc/SCG/mergedb')


        self.outtsv = os.path.join(self.output_directory, 'matching_taxonomy.tsv')


        self.SCGs_fasta = [files for files in os.listdir(
            self.scgs_directory) if not files.endswith(".dmnd")] #FIXME useless?

        self.tanomie_index_tsv={"domain": 1,
        "phylum": 2,
        "class": 3,
        "order": 4,
        "family": 5,
        "genus": 6,
        "species": 7}

        self.sanity_check()

        if self.classic_input_directory:
            self.check_latest_version()



    def sanity_check(self):



        if not filesnpaths.is_file_exists(self.genes_files_directory, dont_raise=True):
            raise ConfigError("Anvi'o could not find gene list file '%s'. If you did not provided any as a parameter \
                               anvi'o looks in '%s'. You can download file by using the commande 'anvi-setup-scgs'."\
                               % (self.genes_files_directory, self.classic_input_directory))

        self.genes_files = os.listdir(self.genes_files_directory)

        if not filesnpaths.is_file_exists(self.path_tsv_taxonomy, dont_raise=True):
            raise ConfigError("Anvi'o could not find gene list file '%s'. If you did not provided any as a parameter \
                               anvi'o looks in '%s'. You can download file by using the commande 'anvi-setup-scgs'."\
                               % (self.path_tsv_taxonomy, self.classic_input_directory))

        if not filesnpaths.is_file_exists(self.hmms, dont_raise=True):
            raise ConfigError("Anvi'o could not find gene list file '%s'. You must declare one before continue."\
                               % self.hmms)


        if filesnpaths.is_output_dir_writable(self.output_directory):
                filesnpaths.gen_output_directory(self.output_directory)

    def check_latest_version(self):
        self.database_url = "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/"

        with open(os.path.join(self.classic_input_directory, 'VERSION')) as f:
            content_local = f.read()

            version_local = content_local.strip().split('\n')[0].strip()
            release_date_local = content_local.strip().split('\n')[2].strip()


        content_online = read_remote_file(self.database_url + 'VERSION',is_gzip=False)
        version_online = content_online.strip().split('\n')[0].strip()
        release_date_online= content_online.strip().split('\n')[2].strip()

        self.run.info("Your current gtdb release", "%s (%s)" % (version_local, release_date_local))

        if version_local != version_online or release_date_local != release_date_online :
            run.warning("There is a more recent version online %s (%s).\
                         You can download it by using the commande 'anvi-setup-scgs'" % (version_online, release_date_online))
        else:
            self.run.info_single("and this is the lates version")


    def make_scg_db(self):

        self.run.info('SCG genes directory', self.genes_files_directory)
        self.run.info('Taxonomy file', self.path_tsv_taxonomy)
        self.run.info('Fasta for HMM reference', self.hmms)
        self.run.info('Output directory', self.output_directory)

        self.dictionary_correspondance_SCGs={}



        self.matrix_taxonomy=self.do_taxonomy_dictonnrary_with_tsv(self.path_tsv_taxonomy)

        self.pathpickle_dictionary_correspondance_SCGs = os.path.join(
            self.output_directory, 'dictionary_correspondance_SCGs.pickle')

        if not os.path.exists(self.pathpickle_dictionary_correspondance_SCGs):
            self.make_correpondance_dictonnary_SCG()

        else:
            with open(self.pathpickle_dictionary_correspondance_SCGs, 'rb') as handle:
                self.dictionary_correspondance_SCGs = pickle.load(handle)


        self.dictionnary_pickel_taxo={}
        if not self.taxonomy_level:
            name_out_put_final_SCG_db="SCG_diamond_db"
        else:
            name_out_put_final_SCG_db=self.taxonomy_level+"SCG_diamond_db"
        for keycorres in self.dictionary_correspondance_SCGs.keys():

            self.output_directory_SCG=os.path.join(self.output_directory, name_out_put_final_SCG_db)

            self.path_new_fasta_SCG=os.path.join(self.output_directory_SCG, keycorres)

            self.path_original_genes=os.path.join(self.genes_files_directory, str(keycorres))

            filesnpaths.gen_output_directory(self.output_directory_SCG)

            self.create_SCG_db(keycorres,self.taxonomy_level)
            path_dictionnary_pickel_taxo_anvi = os.path.join(
                self.output_directory, 'dictionnary_msa_taxonomy_anvio.pickle')

            with open(path_dictionnary_pickel_taxo_anvi, 'wb') as handle:
                pickle.dump(self.dictionnary_pickel_taxo, handle, protocol=pickle.HIGHEST_PROTOCOL)


    def create_SCG_db(self, keycorres, dictionnary_level=None):

        self.listtaxo=[]
        unmatch=[]
        newfasta=""
        with open(self.outtsv,'a') as tsv:
            for refence in self.dictionary_correspondance_SCGs[keycorres]:
                self.path_new_fasta_SCG=os.path.join(self.output_directory_SCG, keycorres+".faa")
                fasta=u.ReadFasta(os.path.join(self.genes_files_directory, str(refence)+".faa"), quiet=True)
                if dictionnary_level!=None:
                    for matrix_code, matrix_taxon in self.matrix_taxonomy.items():
                        index,name=self.match(fasta, matrix_code, matrix_taxon, dictionnary_level)
                        if not index:
                            unmatch.append(name)
                            continue
                        else:
                            self.listtaxo.append(name)
                            print(index,fasta.ids[index],self.matrix_taxonomy[fasta.ids[index]])
                            #sys.exit(status=None)
                            newfasta=newfasta+">"+fasta.ids[index]+"\n"+fasta.sequences[index]+"\n"
                            tsv.write(fasta.ids[index]+"\t"+';'.join(self.matrix_taxonomy[fasta.ids[index]])+"\n")
                            self.dictionnary_pickel_taxo[fasta.ids[index]]=self.matrix_taxonomy[fasta.ids[index]]
                else:
                    i=0
                    for fasta_id in fasta.ids:
                        if fasta_id in self.matrix_taxonomy:
                            newfasta=newfasta+">"+fasta_id+"\n"+fasta.sequences[i]+"\n"
                            if fasta_id not in self.dictionnary_pickel_taxo:
                                self.dictionnary_pickel_taxo[fasta_id]=self.matrix_taxonomy[fasta_id]
                                tsv.write(fasta_id+"\t"+';'.join(self.matrix_taxonomy[fasta_id])+"\n")
                        else:
                            unmatch.append(fasta_id)
                        i+=1
            with open(self.path_new_fasta_SCG,'a') as file:
                file.write(newfasta)
            self.diamonddb_stdin(newfasta,self.path_new_fasta_SCG)


    def match(self, fasta, matrix_code, matrix_taxon, dictionnary_level):
        if dictionnary_level:
            name=str(matrix_taxon[self.taxonomy_level]).rstrip()
        else:
            name=str(matrix_taxon[-1]).rstrip()
        if name not in self.listtaxo and matrix_code in fasta.ids:
            index=fasta.ids.index(matrix_code)
            return(index,name)
        else:
            index=False
            name=(name+"\t"+matrix_taxon[0]+"\n")
            return(index,name)


    def make_correpondance_dictonnary_SCG(self):

        self.diamonddb(self.hmms)
        number_genes_files=len(self.genes_files)
        self.progress.new('Searching for corresponding SCG with Anvio reference', progress_total_items=number_genes_files)
        self.progress.update('...')
        genes_files_number=0
        for gene_file in self.genes_files:
            genes_files_number+=1

            self.name=gene_file.replace('.faa','')
            self.path_diamondb=os.path.join(self.output_directory, "diamonddb")
            self.path_refund_genes=os.path.join(self.path_diamondb, gene_file)
            self.path_original_genes=os.path.join(self.genes_files_directory,gene_file)
            filesnpaths.gen_output_directory(self.path_diamondb)
            self.sequence_to_blast=self.refund_genes()

            #self.progress.update("Reference %s : %d / %d"% self.name, genes_files_number, number_genes_files)
            self.diamond_output=self.diamondblast_stdin(self.sequence_to_blast,self.hmms)

            if not self.diamond_output:
                #run.info("no match", self.name)
                continue

            self.hmm_gene=str(self.diamond_output).split('\n')[0].split('\t')[1]
            if self.hmm_gene not in self.dictionary_correspondance_SCGs:
                self.dictionary_correspondance_SCGs[self.hmm_gene]=[self.name]

            else:
                self.dictionary_correspondance_SCGs[self.hmm_gene]+=[self.name]

            os.remove(self.path_refund_genes)
            os.remove(self.pathrefpath_refund_genesundgenes+'.dmnd')

        self.progress.end()
        if not len(self.dictionary_correspondance_SCGs):
            raise ConfigError("There is no match with SCGs fasta provide and anvi'o ressource")

        self.run.info("Path correspondance dictonnary",path_tsv_taxonomy)

        with open(self.pathpickle_dictionary_correspondance_SCGs, 'wb') as handle:
            pickle.dump(self.dictionary_correspondance_SCGs, handle, protocol=pickle.HIGHEST_PROTOCOL)


    def do_taxonomy_dictonnrary_with_tsv(self, path_tsv_taxonomy, filter=False):
        """format tsv in dictonnary with code as key and liste for taxonomy for value"""
        self.run.info("Correspondance taxonomy file",path_tsv_taxonomy)

        self.matrix={}
        individues_filterd=[]
        number_individues_filterd=0
        number_individues_selected=0
        linestaxo= [line for line in open(path_tsv_taxonomy).readlines()]
        self.run.info("Total Indivues"," %d"% len(linestaxo))
        self.progress.new('Parse tsv file')
        self.progress.update('...')
        for line in linestaxo:
            taxo=[]
            for taxonomy in line.split("\t")[1].split(";"):
                if filter:
                    level_taxo=taxonomy.rstrip()
                    level_taxo=re.sub(r'(\\*[a-z])\_[A-Z]', r'\1',level_taxo)
                    level_taxo=re.sub(r"sp[0-9]*","",level_taxo)
                    level_taxo=re.sub(" ","_",level_taxo)
                    if re.match(r'[a-z]{3,}__[A-Z][a-z]*_[a-z]{3,}',level_taxo) or re.match(r'[a-z]__[A-Z][a-z]{3,}',level_taxo):
                        taxo.append(level_taxo)
                    if not len(taxo):
                        number_individues_filterd+=1
                        continue

                else:
                    level_taxo=re.sub(r'[a-z]__', '',taxonomy).rstrip()
                    taxo.append(level_taxo)

            number_individues_selected+=1
            self.progress.update('Number of indivues selected %s' % number_individues_selected)

            self.matrix[str(line.split("\t")[0])]=taxo

            #.replace(r'_[A-Z]*$','').replace(r"_[A-Z] ","")
        self.progress.end()
        if number_individues_selected==0:
            raise ConfigError("There's a problem with the file '%s'. Anvi'o couldn't use it. \
                              If you are sure that the file is not empty and in the right format try the option '--no-filter'."\
                               % self.path_tsv_taxonomy)

        self.run.info("Indivues Selected"," %d"% int(number_individues_selected))
        return(self.matrix)


    def diamonddb(self, pathgdb):

        diamond = Diamond(run=run_quiet, progress=progress_quiet, num_threads=self.num_threads)
        diamond.query_fasta=pathgdb
        diamond.makedb(pathgdb)


    def diamondblast_stdin(self, sequencetoblast, path_reference):
        diamond = Diamond(run=run_quiet, progress=progress_quiet, num_threads=self.num_threads)
        diamond.target_fasta = path_reference+".dmnd"
        diamond.blastp_stdin(sequencetoblast)

    def diamonddb_stdin(self,sequencetoblast, path_reference):
        diamond = Diamond(run=run_quiet, progress=progress_quiet, num_threads=self.num_threads)
        diamond.makedb_stdin(sequencetoblast, output_file_path=path_reference)



    def diamondblast_stdou(self, path_fasta_db, max_target_seqs=20, evalue=1e-05, min_pct_id=90):
        pathdb=path_fasta_db+".dmnd"
        diamond = Diamond(pathdb, run=run_quiet, progress=progress_quiet, num_threads=self.num_threads)
        diamond_output = diamond.blastp_stdout(query_fasta=self.path_refund_genes, max_target_seqs = max_target_seqs,min_pct_id=min_pct_id)
        return(diamond_output)


    def refund_genes(self, percent_gap=0.5):
        """Use only the sequence with less than percent_gap"""
        fasta=u.ReadFasta(self.path_original_genes, quiet=True)

        i=0
        sequence_to_blast=""
        while i < len(fasta.ids):
            if fasta.sequences[i].count('-') < len(fasta.sequences[i])*percent_gap:
                sequence_to_blast=sequence_to_blast+">"+fasta.ids[i]+"\n"+fasta.sequences[i]+"\n"
            i+=1
        return(sequence_to_blast)


    def trie_blast(self,out_put_path,dictionary_correspondance_SCGs):#FIXME useless
        """ return dictonnary with the name of Genes from given sequence corresponding with the one from anvi'o"""

        pathblastfile=os.path.join(outpath,blastfile)
        if os.stat(pathblastfile).st_size == 0:
            os.remove(pathblastfile)
        hmm_gene, bacgene=self.corre(pathblastfile)
        dictionary_correspondance_SCGs[bacgene]=hmm_gene
        return dictionary_correspondance_SCGs;





class lowident():
    def __init__(self, args, run=run, progress=progress):

        self.args = args
        self.run = run
        self.progress = progress


        A = lambda x: args.__dict__[x] if x in args.__dict__ else None

        self.scgs_directory =A('scgs_directory')
        self.path_tsv_taxonomy=A('taxonomy_files')
        self.output_file=A('output_file')
        self.num_threads=A('num_threads')
        self.num_process=A('num_process')


        self.classic_input_directory = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/SCG/')

        if not self.path_tsv_taxonomy:
            self.path_tsv_taxonomy = os.path.join(self.classic_input_directory,'merge_taxonomy.tsv')

        self.num_threads=args.num_threads

        if not self.num_threads:
            self.num_threads=1

        if not self.num_process:
            self.num_process=1


        if not self.scgs_directory:
            self.scgs_directory = os.path.join(self.classic_input_directory,'msa_individual_genes')


        if not self.output_file:
            self.output_file = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/SCG/mergedb/dico_low_ident.pickle')



        self.genes_files = [files for files in os.listdir(
            self.scgs_directory) if not files.endswith(".dmnd")]

        self.dicolevel = {}


    def process(self):




        self.output_directory = filesnpaths.get_temp_directory_path()
        self.pathpickle_dico_taxo = os.path.join(
            self.output_directory, 'dico_taxo_code_species.pickle')

        if not os.path.exists(self.pathpickle_dico_taxo):
            self.make_dicolevel(self.path_tsv_taxonomy)
            with open(self.pathpickle_dico_taxo, 'wb') as handle:
                pickle.dump(self.dicolevel, handle, protocol=pickle.HIGHEST_PROTOCOL)
        else:
            with open(self.pathpickle_dico_taxo, 'rb') as handle:
                self.dicolevel = pickle.load(handle)


        dico_low_ident = {}
        for genes in self.genes_files:
            dico_low_ident_genes={}
            path_fasta = os.path.join(self.scgs_directory, genes)
            dico_low_ident_genes=self.creatsubfa_ident( path_fasta, genes,dico_low_ident_genes)
            pathpickle_dico_ident = os.path.join(
                self.output_directory, genes+'_dico_low_ident.pickle')
            with open(pathpickle_dico_ident, 'wb') as handle:
                pickle.dump(dico_low_ident_genes, handle, protocol=pickle.HIGHEST_PROTOCOL)
            if genes not in dico_low_ident:
                dico_low_ident[genes]=dico_low_ident_genes
            else:
                dico_low_ident[genes]=dico_low_ident[genes].update(dico_low_ident_genes)

        pathpickle_dico_ident = os.path.join(
            self.output_file)
        with open(self.output_file, 'wb') as handle:
            pickle.dump(dico_low_ident, handle, protocol=pickle.HIGHEST_PROTOCOL)


    def make_dicolevel(self,path_tsv_taxonomy):
        with open(path_tsv_taxonomy, 'r') as taxotsv:
            linestaxo = taxotsv.readlines()
            for line in linestaxo:
                names = line.split('\t')
                code = names[0]
                leveltaxos = names[1].split(';')
                for leveltaxo in leveltaxos:
                    leveltaxo=leveltaxo.rstrip()
                    if leveltaxo not in self.dicolevel:
                        self.dicolevel[leveltaxo] = [code]
                        continue
                    if code in self.dicolevel[leveltaxo]:
                        continue
                    else:
                        self.dicolevel[leveltaxo] = self.dicolevel[leveltaxo] + [code]


    def match_ident(self, fasta, codes, listindex):
        for code in codes:
            if code in fasta.ids:
                index = fasta.ids.index(code)
                if index:
                    listindex.append(index)
        return(listindex)


    def multidiamond(self, listeprocess,dico_low_ident):
        manager = multiprocessing.Manager()
        input_queue = manager.Queue()
        output_queue = manager.Queue()

        num_listeprocess = len(listeprocess)

        progress.new('Aligning amino acid sequences for genes in gene clusters', progress_total_items=num_listeprocess)
        progress.update('...')

        for pathquery in listeprocess:
            input_queue.put(pathquery)

        workers = []
        for i in range(0, self.num_process):
            worker = multiprocessing.Process(target=self.diamond,
                args=(input_queue,output_queue))

            workers.append(worker)
            worker.start()


        finish_process = 0
        while finish_process < num_listeprocess:
            try:
                taxo_ident_item = output_queue.get()

                if taxo_ident_item:
                    dico_low_ident[taxo_ident_item['taxonomy']]= taxo_ident_item['cutoff']


                finish_process += 1
                progress.increment()
                progress.update("Processed %d of %d non-singleton GCs in 10 threads." %
                    (finish_process, num_listeprocess ))

            except KeyboardInterrupt:
                print("Anvi'o profiler recieved SIGINT, terminating all processes...")
                break

        for worker in workers:
            worker.terminate()

        progress.end()

        return dico_low_ident


    def creatsubfa_ident(self, path_fasta, genes,dico_low_ident):

        try:
            fasta=u.ReadFasta(path_fasta, quiet=True)
        except:
            self.run.warning('%s does not seem to be a FASTA file.' % path_fasta)

        listeprocess=[]
        for taxonomy, codes in self.dicolevel.items():
            if taxonomy.startswith("Archaea") or taxonomy.startswith("Bacteria"):
                continue
            listindex = []
            riboname = genes.replace(".faa", "")
            path_new_fasta_SCG = os.path.join(self.output_directory,taxonomy)
            pathpickle_dico_ident = path_new_fasta_SCG + "_dico_low_ident.pickle"

            if not os.path.exists(path_new_fasta_SCG):
                listindex = self.match_ident(fasta, codes, listindex)
                if len(listindex) > 1:
                    listeprocess.append(path_new_fasta_SCG)
                    with open(path_new_fasta_SCG, 'w+') as file:
                        for index in listindex:
                            file.write(">" + fasta.ids[index] +
                                       "\n" + fasta.sequences[index] + "\n")
                else:
                    if genes not in dico_low_ident:
                        dico_low_ident[riboname]={}
                        dico_low_ident[riboname][taxonomy]=100
                    else:
                        dico_low_ident[riboname][taxonomy]=100

        if listeprocess:
            dico_low_ident=self.multidiamond(listeprocess,dico_low_ident)
            return(dico_low_ident)
        else:
            dico_low_ident={}
            return(dico_low_ident)



    def diamond(self, input_queue,output_queue):
        while True:
            pathquery = input_queue.get(True)
            pathdb=pathquery+".dmnd"
            path_diamond=pathquery+'.txt'
            taxonomy=os.path.basename(pathquery)

            if not os.path.exists(path_diamond):
                self.diamonddb(pathquery,pathdb)
                ouputdiamond=self.run_diamond(pathquery,pathdb)

                os.remove(pathdb)
                os.remove(pathquery+'log_file')
            low_ident = self.select_low_ident(ouputdiamond)
            os.remove(pathquery)
            output = {'taxonomy': taxonomy, 'cutoff': low_ident}
            output_queue.put(output)


    def diamonddb(self, pathquery,pathdb):

        diamond = Diamond(query_fasta=pathquery,run=run_quiet, progress=progress_quiet, num_threads=self.num_threads)
        diamond.query_fasta=pathquery
        diamond.run.log_file_path=pathquery+'log_file'
        diamond.target_fasta = pathquery
        diamond.num_threads = self.num_threads
        diamond.makedb()


    def run_diamond(self, pathquery,pathdb):
        diamond = Diamond(run=run_quiet, progress=progress_quiet, num_threads=self.num_threads)

        diamond.evalue=None
        diamond.run.log_file_path=pathquery+'log_file'
        diamond.target_fasta = pathdb
        diamond.query_fasta=pathquery
        diamond.max_target_seqs=None
        diamond.search_output_path = pathquery
        diamond.tabular_output_path = pathquery + '.txt'
        output=diamond.blastp_stdout()

        return output


    def select_low_ident(self,str_diamond_output,lowest_ident=100):
        "Select the lowest percent identity on aligment output"
        for line in str_diamond_output.split('\n'):
            if line:
                ident = line.strip().split('\t')[2]
                if float(ident) < float(lowest_ident):
                    lowest_ident = ident
        return(lowest_ident)
