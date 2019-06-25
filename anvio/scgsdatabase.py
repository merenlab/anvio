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

        self.genesfilesdirectory=args.genesfilesdirectory

        self.scgsdirectory=args.scgsdirectory

        #self.num_threads=args.num_threads



        self.dicolevel=self.diconlevel()



        if not args.taxofiles:
            self.taxofiles = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/SCG/merge_taxonomy.tsv')

        #if not args.outtsv:
        self.outtsv = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/SCG/mergedb/matching_taxonomy.tsv')

        #if not args.num_threads:
        self.num_threads=1


        self.hmms=args.hmms

        if not args.genesfilesdirectory:
            self.genesfilesdirectory = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/SCG/msa_individual_genes')

        self.genesfiles=os.listdir(self.genesfilesdirectory)


        if not args.outputdirectory:
            outputdirectory = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/SCG/mergedb')
            if not os.path.exists(outputdirectory):
                os.mkdir(os.path.dirname(anvio.__file__)+'/data/misc/SCG/mergedb')
            self.outputdirectory = outputdirectory



        self.scgsfasa = [files for files in os.listdir(
            self.scgsdirectory) if not files.endswith(".dmnd")]




        dicocorres={}
        taxomatrix=self.domatrix(self.taxofiles)
        dico_pickel_taxo={}
        for gene in self.genesfiles:
            #for dbdiamond
            name=gene.replace('.faa','')
            outpathdb=os.path.join(self.outputdirectory, "diamonddb")
            pathdb=os.path.join(self.genesfilesdirectory, name)
            pathgenes=os.path.join(self.genesfilesdirectory, gene)
            outpath=os.path.join(self.outputdirectory, "diamondblast")
            pathrefundgenes=os.path.join(outpathdb, gene)
            pathblast=os.path.join(outpath,name)
            if not os.path.exists(outpathdb):
                os.mkdir(outpathdb)
                #pathdb=self.diamonddb(pathdb,pathgenes)
                #self.diamondblast(pathdb,pathblast)


            if not os.path.exists(outpath):
                os.mkdir(outpath)
                outpath=os.path.join(self.outputdirectory, "diamondblast")


            sequencetoblast=self.refundgenes(pathgenes,pathrefundgenes,pathdb)
            self.diamonddb_stdin(sequencetoblast,pathrefundgenes)
            diamond_output=self.diamondblast_stdou(self.hmms,pathrefundgenes)
            #print(diamond_output)
            if not diamond_output:
                run.info("no match", gene)
                os.remove(pathrefundgenes)
                os.remove(pathrefundgenes+'.dmnd')
                continue

            #self.trie_blast(diamond_output,dicocorres)
            hmmgene=str(diamond_output).split('\n')[0].split('\t')[0]
            #hmmgene, bacgene=self.make_dico_correspondance(diamond_output)
            if hmmgene not in dicocorres:
                dicocorres[hmmgene]=[name]

            else:
                dicocorres[hmmgene]=dicocorres[hmmgene]+[name]





        keylevel="species"
        print(dicocorres)
        for keycorres in dicocorres.keys():
            outputsub=os.path.join(self.outputdirectory, keylevel)
            pathfile=os.path.join(outputsub, keycorres)
            pathfa=os.path.join(self.genesfilesdirectory, str(keycorres))
            if not os.path.exists(outputsub):
                os.mkdir(outputsub)
            dico_pickel_taxo=self.creatsubfa(taxomatrix,dicocorres[keycorres],pathfile,keylevel,dico_pickel_taxo)

        #self.cleandir(outpath)
        #self.cleandir(outpathdb)
        sys.exit()
        if True:
            pathpickle_dico_taxo = os.path.join(
                self.outputdirectory, 'dico_taxo_code_species.pickle')

            if not os.path.exists(pathpickle_dico_taxo):
                dicolevel = self.make_dicolevel(self.taxofiles)
                with open(pathpickle_dico_taxo, 'wb') as handle:
                    pickle.dump(dicolevel, handle, protocol=pickle.HIGHEST_PROTOCOL)
            else:
                with open(pathpickle_dico_taxo, 'rb') as handle:
                    dicolevel = pickle.load(handle)

            dico_low_ident = {}
            for scgsfasta in os.listdir(outputsub):
                pathfa = os.path.join(self.genesfilesdirectory, scgsfasta )
                self.creatsubfa_ident(dicolevel, pathfa, self.outputdirectory, scgsfasta )

            pathpickle_dico_ident = os.path.join(
                self.outputdirectory, 'dico_low_ident.pickle')
            with open(pathpickle_dico_ident, 'wb') as handle:
                pickle.dump(dico_low_ident, handle, protocol=pickle.HIGHEST_PROTOCOL)


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
                for name in taxos:
                    taxo= re.sub(r'_[A-Z]*$', '', name).rstrip()
                    leveltaxo1 = re.sub(r'sp[1-9]*', '', taxo)
                    leveltaxo3 = re.sub(r'(\\*[a-z])\_[A-Z]', r'\1', leveltaxo1)
                    leveltaxo2 = leveltaxo3.rstrip()
                    leveltaxo = re.sub(r'\s+', r'_', leveltaxo2)
                    last_taxo=leveltaxo.rstrip().replace(" ","_")
                matrix[code]=last_taxo
            return(matrix)


    def creatsubfa(self,taxomatrix,listerefence,pathfile,keylevel,dico_pickel_taxo,dicolevel=None):

        listtaxo=[]
        unmatch=[]
        newfasta=""
        with open(self.outtsv,'a') as tsv:
            with open(pathfile,'a') as file:
                print(listerefence)
                for refence in listerefence:
                    print(refence)
                    pathfa=os.path.join(self.outputdirectory, "diamonddb",refence+".faa")
                    run.info("pathfa",pathfa)
                    fasta=u.ReadFasta(pathfa)
                    if dicolevel!=None:
                        for code, taxonomy in taxomatrix.items():
                            index,name=self.match(fasta,keylevel,taxonomy,code,listtaxo)
                            if not index:
                                unmatch.append(name)
                                continue
                            else:
                                listtaxo.append(name)
                                newfasta=newfasta+">"+fasta.ids[index]+"\n"+fasta.sequences[index]+"\n"
                                tsv.write(fasta.ids[index]+"\t"+';'.join(taxomatrix[fasta.ids[index]])+"\n")
                                dico_pickel_taxo[fasta.ids[index]]=taxomatrix[fasta.ids[index]]
                    else:
                        i=0
                        while i<len(fasta.ids):
                            """some fasta id are not relate to any taxonomy ?!"""
                            try:
                                tsv.write(fasta.ids[i]+"\t"+';'.join(taxomatrix[fasta.ids[i]])+"\n")
                                dico_pickel_taxo[fasta.ids[i]]=taxomatrix[fasta.ids[i]]
                                newfasta=newfasta+">"+fasta.ids[i]+"\n"+fasta.sequences[i]+"\n"
                            except:
                                unmatch.append(fasta.ids[i])

                            i+=1
                file.write(newfasta)
                self.diamonddb_stdin(newfasta,pathfile)
                #print(unmatch)
        return(dico_pickel_taxo)

    def match(self,fasta,keylevel,taxonomy,code,listtaxo,dicolevel=False):
        if dicolevel:
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
            pathdb, outpathdb=self.diamonddb(genes,pathdb,pathgenes)
            outpath=self.diamondblast(pathdb,genes,self.hmms,self.outputdirectory)
        return(outpath,outpathdb)

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



    def make_dico_correspondance(self,diamond_output):
        listeline=str(diamond_output).split('\n')[0].split('\t')[0]
        bestscore = 0
        i=0
        for line in listeline:
            score = line[11].rstrip()
            if i==3:
                bacgene=os.path.basename(pathblast)
                hmmgene=line
                return(hmmgene,bacgene)
                break
            if float(score) > float(bestscore):
                bestscore=score
                bestname=line[0]
            if bestname==line[0]:
                i+=1



    def refundgenes(self,pathgenes,pathrefundgenes,outpathdb):
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




    def trie_blast(self,outpath,dicocorres):
        """ return dictonnary with the name of Genes from given sequence corresponding with the one from anvi'o"""

        pathblastfile=os.path.join(outpath,blastfile)
        if os.stat(pathblastfile).st_size == 0:
            os.remove(pathblastfile)
        hmmgene, bacgene=self.corre(pathblastfile)
        dicocorres[bacgene]=hmmgene
        return dicocorres;



    def cleandir(self,dir):
        shutil.rmtree(dir)

    def diconlevel(self):
        dicolevel= {}
        dicolevel["domain"]=1
        dicolevel["phylum"]=2
        dicolevel["class"]=3
        dicolevel["order"]=4
        dicolevel["family"]=5
        dicolevel["genus"]=6
        dicolevel["species"]=7
        return(dicolevel)

    def make_dicolevel(self,taxofiles):
        with open(taxofiles, 'r') as taxotsv:
            linestaxo = taxotsv.readlines()
            dicolevel = {}
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
                    if leveltaxo not in dicolevel:
                        dicolevel[leveltaxo] = [code]
                        continue
                    if code in dicolevel[leveltaxo]:
                        continue
                    else:
                        dicolevel[leveltaxo] = dicolevel[leveltaxo] + [code]
        return(dicolevel)


    def match_ident(self,fasta, taxonomy, codes, listindex, dicolevel=False):
        for code in codes:
            if code in fasta.ids:
                index = fasta.ids.index(code)
                if index:
                    listindex.append(index)
        return(listindex)


    def creatsubfa_ident(self,dicolevel, pathfa, outputdirectory, genes):
        fasta = self.getfasta(pathfa)

        unmatch = []

        dico_low_ident = {}
        for taxonomy, codes in dicolevel.items():
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
                dico_low_ident = select_low_ident(pathfile, genes, dico_low_ident)
                pathpickle_dico_ident = pathfile + "_dico_low_ident.pickle"
                with open(pathpickle_dico_ident, 'wb') as handle:
                    pickle.dump(dico_low_ident, handle,
                                protocol=pickle.HIGHEST_PROTOCOL)
        return(dico_low_ident)


    def select_low_ident(self,pathfile, genes, dico_low_ident):
        low_ident = 101
        genes_taxo = os.path.basename(pathfile)
        with open(pathfile + ".diamond", 'r') as file:
            lines = file.readlines()
            for line in lines:
                info = line.split("\t")
                ident = info[3].rstrip()
                if float(ident) < float(low_ident):
                    low_ident = ident
            dico_low_ident[genes_taxo] = low_ident
        return(dico_low_ident)


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
