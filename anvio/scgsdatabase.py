#!/usr/bin/env python3

import os
import gzip
import shutil
import requests
from io import BytesIO
import anvio
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
import os
import pickle
import re
import shutil
import subprocess
import sys
import anvio.fastalib as u
from anvio.drivers.diamond import Diamond
import anvio.terminal as terminal
import anvio.pfam as pfam
import tarfile


#run_quiet = terminal.Run(verbose=False)
#progress_quiet = terminal.Progress(verbose=False)

_author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Quentin Clayssen"
__email__ = "quentin.clayssen@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print
run_quiet = terminal.Run(verbose=False)
progress_quiet = terminal.Progress(verbose=False)



class scgsdatabase():
    def __init__(self, args, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None

        self.genes_files_directory = args.genes_files_directory

        self.scgs_directory = args.scgs_directory

        #self.num_threads=args.num_threads



        self.dictionnary_level=self.dictionnary_level()



        if not args.taxofiles:
            self.taxofiles = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/SCG/merge_taxonomy.tsv')

        #if not args.outtsv:
        self.outtsv = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/SCG/mergedb/matching_taxonomy.tsv')

        #if not args.num_threads:
        self.num_threads=args.num_threads

        if not args.num_threads:
            self.num_threads=2


        self.hmms=args.hmms

        if not args.genes_files_directory:
            self.genes_files_directory = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/SCG/msa_individual_genes')

        self.genesfiles=os.listdir(self.genes_files_directory)


        if not args.outputdirectory:
            outputdirectory = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/SCG/mergedb')
            if not os.path.exists(outputdirectory):
                os.mkdir(os.path.dirname(anvio.__file__)+'/data/misc/SCG/mergedb')
            self.outputdirectory = outputdirectory



        self.scgsfasa = [files for files in os.listdir(
            self.scgs_directory) if not files.endswith(".dmnd")]




        dictionnary_corres={}

        matrix_taxonomy=self.domatrix(self.taxofiles)

        pathpickle_dictionnary__corres = os.path.join(
            self.outputdirectory, 'dictionnary__corres.pickle')
        if not os.path.exists(pathpickle_dictionnary__corres):
            for gene in self.genesfiles:
                #for dbdiamond
                name=gene.replace('.faa','')
                path_diamondb=os.path.join(self.outputdirectory, "diamonddb")
                pathdb=os.path.join(self.genes_files_directory, name)
                pathgenes=os.path.join(self.genes_files_directory, gene)
                outpath=os.path.join(self.outputdirectory, "diamondblast")
                pathrefundgenes=os.path.join(path_diamondb, gene)
                pathblast=os.path.join(outpath,name)
                if not os.path.exists(path_diamondb):
                    os.mkdir(path_diamondb)
                    #pathdb=self.diamonddb(pathdb,pathgenes)
                    #self.diamondblast(pathdb,pathblast)


                if not os.path.exists(outpath):
                    os.mkdir(outpath)
                    outpath=os.path.join(self.outputdirectory, "diamondblast")


                #sequencetoblast=self.refundgenes(pathgenes,pathrefundgenes,pathdb)
                #self.diamonddb_stdin(sequencetoblast,pathrefundgenes)
                if not os.path.exists(pathrefundgenes):
                    continue
                diamond_output=self.diamondblast_stdou(self.hmms,pathrefundgenes)
                #print(diamond_output)
                if not diamond_output:
                    run.info("no match", gene)
                    os.remove(pathrefundgenes)
                    os.remove(pathrefundgenes+'.dmnd')
                    continue

                #self.trie_blast(diamond_output,dictionnary_corres)
                hmm_gene=str(diamond_output).split('\n')[0].split('\t')[0]
                #hmm_gene, bacgene=self.make_dictionnary__correspondance(diamond_output)
                if hmm_gene not in dictionnary_corres:
                    dictionnary_corres[hmm_gene]=[name]

                else:
                    dictionnary_corres[hmm_gene]=dictionnary_corres[hmm_gene]+[name]

                with open(pathpickle_dictionnary__corres, 'wb') as handle:
                    pickle.dump(dictionnary_corres, handle, protocol=pickle.HIGHEST_PROTOCOL)

        else:
            with open(pathpickle_dictionnary__corres, 'rb') as handle:
                dictionnary_corres = pickle.load(handle)




        keylevel = "species"

        dictionnary__pickel_taxo={}
        for keycorres in dictionnary_corres.keys():

            outputsub=os.path.join(self.outputdirectory, keylevel)

            pathfile=os.path.join(outputsub, keycorres)

            pathfa=os.path.join(self.genes_files_directory, str(keycorres))

            if not os.path.exists(outputsub):
                os.mkdir(outputsub)

            dictionnary__pickel_taxo=self.creatsubfa(matrix_taxonomy,dictionnary_corres[keycorres],pathfile,keylevel,dictionnary__pickel_taxo)
            pathpickle_dictionnary__taxo_anvi = os.path.join(
                self.outputdirectory, 'dictionnary__msa_taxonomy_anvio.pickle')

            with open(pathpickle_dictionnary__taxo_anvi, 'wb') as handle:
                pickle.dump(dictionnary__pickel_taxo, handle, protocol=pickle.HIGHEST_PROTOCOL)

        #self.cleandir(outpath)
        #self.cleandir(path_diamondb)
        sys.exit()
        if True:
            pathpickle_dictionnary__taxo = os.path.join(
                self.outputdirectory, 'dictionnary__taxo_code_species.pickle')

            if not os.path.exists(pathpickle_dictionnary__taxo):
                dictionnary_level = self.make_dictionnary_level(self.taxofiles)
                with open(pathpickle_dictionnary__taxo, 'wb') as handle:
                    pickle.dump(dictionnary_level, handle, protocol=pickle.HIGHEST_PROTOCOL)
            else:
                with open(pathpickle_dictionnary__taxo, 'rb') as handle:
                    dictionnary_level = pickle.load(handle)

            dictionnary__low_ident = {}
            for scgsfasta in os.listdir(outputsub):
                pathfa = os.path.join(self.genes_files_directory, scgsfasta )
                self.creatsubfa_ident(dictionnary_level, pathfa, self.outputdirectory, scgsfasta )

            pathpickle_dictionnary__ident = os.path.join(
                self.outputdirectory, 'dictionnary__low_ident.pickle')
            with open(pathpickle_dictionnary__ident, 'wb') as handle:
                pickle.dump(dictionnary__low_ident, handle, protocol=pickle.HIGHEST_PROTOCOL)


    def is_database_exists(self):
        if os.path.exists(os.path.join(self.SCG_data_dir, 'SCG-A.hmm.gz')):
            raise ConfigError("It seems you already have SCG database installed in '%s', please use --reset flag if you want to re-download it." % self.SCG_data_dir)


    def get_remote_version(self):
        content = read_remote_file(self.database_url + '/VERSION')


    def getfasta(self,pathfa):
        fasta = u.ReadFasta(pathfa)
        return(fasta)

    def domatrix(self,taxofile):
        """format tsv in dictonnary with code as key and liste for taxonomy for value"""
        with open(taxofile,'r') as taxo:
            linestaxo=taxo.readlines()
            listspecies=[]
            matrix={}
            for line in linestaxo:
                names=line.split("\t")
                code=str(names[0])
                taxos=names[1].split(";")
                """for name in taxos[0:-2]:
                    if not re.match('[a-z]__[A-Z]{1}[a-z]{1,}',name):"""

                taxo= re.sub(r'_[A-Z]*$', '', taxos[-1]).rstrip()
                leveltaxo3 = re.sub(r'(\\*[a-z])\_[A-Z]', r'\1', taxo)
                last_taxo1 = re.sub(r'sp[0-9]{1,}$', '', leveltaxo3)
                last_taxo=last_taxo1.rstrip(" ").replace(" ","_")
                if not re.match('s__[A-Z]{1}[a-z]{1,}',last_taxo):
                    continue
                taxos[-1]=last_taxo

                matrix[code]=taxos
            return(matrix)



    def creatsubfa(self,matrix_taxonomy,listerefence,pathfile,keylevel,dictionnary__pickel_taxo,dictionnary_level=None):

        listtaxo=[]
        unmatch=[]
        newfasta=""
        with open(self.outtsv,'a') as tsv:
            with open(pathfile,'a') as file:
                for refence in listerefence:
                    pathfa=os.path.join(self.outputdirectory, "diamonddb",refence+".faa")
                    run.info("pathfa",pathfa)
                    fasta=u.ReadFasta(pathfa)
                    if dictionnary_level!=None:
                        for code, taxonomy in matrix_taxonomy.items():
                            index,name=self.match(fasta,keylevel,taxonomy,code,listtaxo)
                            if not index:
                                unmatch.append(name)
                                continue
                            else:
                                listtaxo.append(name)
                                newfasta=newfasta+">"+fasta.ids[index]+"\n"+fasta.sequences[index]+"\n"
                                tsv.write(fasta.ids[index]+"\t"+';'.join([matrix_taxonomy[fasta.ids[index]]])+"\n")
                                dictionnary__pickel_taxo[fasta.ids[index]]=matrix_taxonomy[fasta.ids[index]]
                    else:
                        i=0
                        while i<len(fasta.ids):
                            if fasta.ids[i] in matrix_taxonomy:
                                newfasta=newfasta+">"+fasta.ids[i]+"\n"+fasta.sequences[i]+"\n"
                                if fasta.ids[i] not in dictionnary__pickel_taxo:
                                    dictionnary__pickel_taxo[fasta.ids[i]]=matrix_taxonomy[fasta.ids[i]]
                                    tsv.write(fasta.ids[i]+"\t"+';'.join(matrix_taxonomy[fasta.ids[i]])+"\n")

                            else:
                                unmatch.append(fasta.ids[i])
                            i+=1
                file.write(newfasta)
                self.diamonddb_stdin(newfasta,pathfile)
                #print(unmatch)
        return(dictionnary__pickel_taxo)

    def match(self,fasta,keylevel,taxonomy,code,listtaxo,dictionnary_level=False):
        if dictionnary_level:
            name=str(taxonomy[keylevel]).rstrip()
        else:
            name=str(taxonomy[-1]).rstrip()
        if name not in listtaxo and code in fasta.ids:
            index=fasta.ids.index(code)
            return(index,name)
        else:
            index=False
            name=(name+"\t"+taxonomy[0]+"\n")
            return(index,name)

    """def db_blast(self,genesfiles,outputdirectory,pathdb,pathgenes):
        for genes in genesfiles:
            pathrefundgenes=self.refundgenes(genes)
            pathdb, path_diamondb=self.diamonddb(genes,pathdb,pathgenes)
            outpath=self.diamondblast(pathdb,genes,self.hmms,self.outputdirectory)
        return(outpath,path_diamondb)

    def diamonddb(self,pathdb,pathgenes):
        cmddb = "diamond makedb --in "+pathgenes+" -d "+pathdb
        os.system(cmddb)
        return(pathdb)

    def diamondblast(self,pathdb,pathblast):
        #with open(self.hmms,'r') as hmmfile:
        #print(pathdb, self.hmms)
        cmdblast = "diamond blastp -d "+pathdb+".dmnd"+" -q "+self.hmms+" -o "+pathblast
        print(cmdblast)
        os.system(cmdblast)
        #return(outpath)"""
    """def get_hit_diamond(self,sequence,pathdb, max_target_seqs=20, evalue=1e-05, min_pct_id=90):
        pathdb=pathdb+".dmnd"
        #diamond = Diamond(pathdb)#, run=run_quiet, progress=progress_quiet
        diamond.max_target_seqs = max_target_seqs
        diamond.evalue = evalue
        diamond.min_pct_id = min_pct_id

        diamond_output = diamond.blastp_stdin(sequence,pathdb)

    def run_diamond(self,query_fasta, target_fasta,name):
        target_fasta=target_fasta+".dmnd"
        diamond = Diamond(query_fasta, target_fasta)

        diamond.max_target_seqs=20
        diamond.search_output_path = ("search_"+name)
        diamond.tabular_output_path = ("tabular_"+name)

        diamond.blastp()
        diamond_output = diamond.view()
        #sys.exit(status=None)
        return diamond_output"""


    def diamonddb_stdin(self,sequencetoblast, pathgene):

        diamond = Diamond()
        diamond.target_fasta = pathgene
        run.info("diamonddb_stdin", pathgene)
        diamond.makedb_stdin(sequencetoblast)

    def diamondblast_stdin(self,sequencetoblast, pathgene):

        diamond = Diamond()
        diamond.target_fasta = pathgene
        run.info("diamonddb_stdin", pathgene)
        diamond.makedb_stdin(sequencetoblast)



    def diamondblast_stdou(self,pathquery,pathdb, max_target_seqs=20, evalue=1e-05, min_pct_id=90):
        pathdb=pathdb+".dmnd"
        diamond = Diamond(pathdb,run=self.run, progress=self.progress, num_threads=self.num_threads)
        #diamond = Diamond(pathdb)#, run=run_quiet, progress=progress_quiet
        diamond.max_target_seqs = max_target_seqs
        diamond.query_fasta=pathquery
        run.info("diamondblast_stdou", pathdb)
        diamond_output = diamond.blastp_stdout()
        str_diamond_output=diamond_output.decode("utf-8")
        return(str_diamond_output)



    def make_dictionnary__correspondance(self,diamond_output):
        listeline=str(diamond_output).split('\n')[0].split('\t')[0]
        bestscore = 0
        i=0
        for line in listeline:
            score = line[11].rstrip()
            if i==3:
                bacgene=os.path.basename(pathblast)
                hmm_gene=line
                return(hmm_gene,bacgene)
                break
            if float(score) > float(bestscore):
                bestscore=score
                bestname=line[0]
            if bestname==line[0]:
                i+=1



    def refundgenes(self,pathgenes,pathrefundgenes,path_diamondb):
        print(pathgenes)
        fasta=u.ReadFasta(pathgenes)
        i=0
        sequencetoblast=""
        while i < len(fasta.ids):
            if fasta.sequences[i].count('-') < len(fasta.sequences[i])*0.5:
                sequencetoblast=sequencetoblast+">"+fasta.ids[i]+"\n"+fasta.sequences[i]+"\n"
            i+=1
        with open(pathrefundgenes,'w') as newfastadb:
            newfastadb.write(sequencetoblast)
        return(sequencetoblast)




    def trie_blast(self,outpath,dictionnary_corres):
        """ return dictonnary with the name of Genes from given sequence corresponding with the one from anvi'o"""

        pathblastfile=os.path.join(outpath,blastfile)
        if os.stat(pathblastfile).st_size == 0:
            os.remove(pathblastfile)
        hmm_gene, bacgene=self.corre(pathblastfile)
        dictionnary_corres[bacgene]=hmm_gene
        return dictionnary_corres;



    def cleandir(self,dir):
        shutil.rmtree(dir)

    def dictionnary_level(self):
        dictionnary_level= {}
        dictionnary_level["domain"]=1
        dictionnary_level["phylum"]=2
        dictionnary_level["class"]=3
        dictionnary_level["order"]=4
        dictionnary_level["family"]=5
        dictionnary_level["genus"]=6
        dictionnary_level["species"]=7
        return(dictionnary_level)

    def make_dictionnary_level(self,taxofiles):
        with open(taxofiles, 'r') as taxotsv:
            linestaxo = taxotsv.readlines()
            dictionnary_level = {}
            for line in linestaxo:
                names = line.split('\t')
                code = names[0]
                leveltaxos = names[1].split(';')
                for leveltaxo in leveltaxos[::-1]:
                    leveltaxo1 = re.sub(r'sp[1-9]*', '', leveltaxo)
                    leveltaxo3 = re.sub(r'(\\*[a-z])\_[A-Z]', r'\1', leveltaxo1)
                    leveltaxo2 = leveltaxo3.rstrip()
                    leveltaxo = re.sub(r'\s+', r'_', leveltaxo2)
                    if leveltaxo == 's__':
                        continue
                    if leveltaxo == 'd__Bacteria':
                        continue
                    if leveltaxo not in dictionnary_level:
                        dictionnary_level[leveltaxo] = [code]
                        continue
                    if code in dictionnary_level[leveltaxo]:
                        continue
                    else:
                        dictionnary_level[leveltaxo] = dictionnary_level[leveltaxo] + [code]
        return(dictionnary_level)


    def match_ident(self,fasta, taxonomy, codes, listindex, dictionnary_level=False):
        for code in codes:
            if code in fasta.ids:
                index = fasta.ids.index(code)
                if index:
                    listindex.append(index)
        return(listindex)


    def creatsubfa_ident(self,dictionnary_level, pathfa, outputdirectory, genes):
        fasta = self.getfasta(pathfa)

        unmatch = []

        dictionnary__low_ident = {}
        for taxonomy, codes in dictionnary_level.items():
            listindex = []

            if taxonomy == "d__Bacteria":
                continue
            listindex = self.match_ident(fasta, taxonomy, codes, listindex)
            if listindex:
                riboname = genes.replace(".faa", "")
                pathfile = os.path.join(self.outputdirectory, riboname + "_" + taxonomy)
                with open(pathfile, 'a') as file:
                    for index in listindex:
                        file.write(">" + fasta.ids[index] +
                                   "\n" + fasta.sequences[index] + "\n")
                if not os.path.exists(pathfile + ".diamond"):
                    cmddb = "diamond makedb --in " + pathfile + " -d " + pathfile
                    os.system(cmddb)
                    cmddb = "diamond blastp -d " + pathfile + ".dmnd -q " + \
                        pathfile + " -o " + pathfile + ".diamond"
                    os.system(cmddb)
                #result = subprocess.run(['diamond', 'blastp', '-d', pathfile+'.dmnd', '-q', pathfile], stdout=subprocess.PIPE)
                # diamond=result.stdout
                # print(diamond)
                dictionnary__low_ident = select_low_ident(pathfile, genes, dictionnary__low_ident)
                pathpickle_dictionnary__ident = pathfile + "_dictionnary__low_ident.pickle"
                with open(pathpickle_dictionnary__ident, 'wb') as handle:
                    pickle.dump(dictionnary__low_ident, handle,
                                protocol=pickle.HIGHEST_PROTOCOL)
        return(dictionnary__low_ident)


    def select_low_ident(self,pathfile, genes, dictionnary__low_ident):
        low_ident = 101
        genes_taxo = os.path.basename(pathfile)
        with open(pathfile + ".diamond", 'r') as file:
            lines = file.readlines()
            for line in lines:
                info = line.split("\t")
                ident = info[3].rstrip()
                if float(ident) < float(low_ident):
                    low_ident = ident
            dictionnary__low_ident[genes_taxo] = low_ident
        return(dictionnary__low_ident)


def read_remote_file(url, is_gzip=True):
    remote_file = requests.get(url)

    if remote_file.status_code == 404:
        raise Exception("'%s' returned 404 Not Found. " % url)

    return remote_file.content.decode('utf-8')


class SCGsSetup(object):
    def __init__(self, args, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress
        self.SCG_data_dir = args.scgs_data_dir

        filesnpaths.is_program_exists('hmmpress')

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
        content = read_remote_file(self.database_url + '/VERSION')

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
