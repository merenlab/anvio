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



class scgsdatabase():
    def __init__(self, args, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None

        self.taxofiles=args.taxofiles

        self.hmms=args.hmms


        if not args.genesfilesdirectory:
            self.genesfilesdirectory = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/SCG')
        else:
            self.genesfilesdirectory=args.genesfilesdirectory

        if not args.outputdirectory:
            self.outputdirectory = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/SCG/mergedb')
        else:
            self.genesfilesdirectory=args.genesfilesdirectory

        self.scgsdirectory=args.scgsdirectory

        self.genesfiles=os.listdir(self.genesfilesdirectory)

        self.dicolevel=self.diconlevel()

        self.scgsfasa = [files for files in os.listdir(
            self.scgsdirectory) if not files.endswith(".dmnd")]


        def is_database_exists(self):
            if not os.path.exists(os.path.join(self.pfam_data_dir, 'bac120_msa_individual_genes.tar.gz')):
                raise ConfigError("It seems you do not have SCGS database installed, please run 'anvi-setup-scgs' to download it.")


        def get_version(self):
            with open(os.path.join(self.pfam_data_dir, 'Pfam.version')) as f:
                content = f.read()


        for genes in self.genesfiles:
            #for dbdiamond
            name=genes.replace('.faa','')
            outpathdb=os.path.join(self.outputdirectory, "diamonddb")
            pathdb=os.path.join(str(outpathdb), str(name))
            pathgenes=os.path.join(self.genesfilesdirectory, genes)
            if not os.path.exists(outpathdb):
                os.mkdir(outpathdb)

            #for diamond sequences
            outpath=os.path.join(self.outputdirectory, "diamondblast")
            pathblast=os.path.join(outpath,genes)
            if not os.path.exists(outpath):
                os.mkdir(outpath)
            #print(pathrefundgenes)

            #sequencetoblast=refundgenes(pathrefundgenes,pathdb)
            pathdb=self.diamonddb(genes,pathdb,pathgenes)
            tabular_output_path=self.diamondblast(pathdb,genes,self.hmms,pathblast)
            #diamondblast(pathdb,genes,hmms,self.outputdirectory)

            #outpath, outpathdb=db_blast(genesfiles,hmms,self.outputdirectory,pathdb,pathgenes)
            outpath=os.path.join(self.outputdirectory, "diamondblast")
            dicocorres=self.trie_blast(outpath)
            taxomatrix=self.domatrix(self.taxofiles)
            keylevel="species"
            for keycorres in dicocorres:
                outputsub=os.path.join(self.outputdirectory, keylevel)
                pathfile=os.path.join(outputsub, dicocorres[keycorres])
                pathfa=os.path.join(self.genesfilesdirectory, str(keycorres))
                if not os.path.exists(outputsub):
                    os.mkdir(outputsub)
                self.creatsubfa(taxomatrix,pathfa,pathfile,keylevel,dicocorres[keycorres])

        #self.cleandir(outpath)
        #self.cleandir(outpathdb)

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
            for scgsfasta in os.listdir(pathfile):
                pathfa = os.path.join(genesfilesdir, scgsfasta )
                self.creatsubfa_ident(dicolevel, pathfa, self.outputdirectory, scgsfasta )

            pathpickle_dico_ident = os.path.join(
                self.outputdirectory, 'dico_low_ident.pickle')
            with open(pathpickle_dico_ident, 'wb') as handle:
                pickle.dump(dico_low_ident, handle, protocol=pickle.HIGHEST_PROTOCOL)



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
                taxo=names[1].split(";")
                taxo[-1]=taxo[-1].rstrip().replace(" ","_")

                for name in names:
                    name = re.sub(r'_[A-Z]*$', '', name)
                matrix[code]=taxo
            return(matrix)


    def creatsubfa(self,taxomatrix,pathfa,pathfile,keylevel,dicolevel=None):
        fasta=self.getfasta(pathfa)
        listtaxo=[]
        unmatch=[]
        with open(pathfile,'a') as file:
            for code, taxonomy in taxomatrix.items():
                index,name=self.match(fasta,keylevel,taxonomy,code,listtaxo)
                if not index:
                    unmatch.append(name)
                    continue
                else:
                    listtaxo.append(name)
                    file.write(">"+fasta.ids[index]+"\n"+fasta.sequences[index]+"\n")
        cmddb = "diamond makedb --in "+pathfile+" -d "+pathfile
        os.system(cmddb)
        return(pathfile)

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

    def db_blast(self,genesfiles,hmms,outputdirectory,pathdb,pathgenes):
        for genes in genesfiles:
            pathrefundgenes=self.refundgenes(genes)
            pathdb, outpathdb=self.diamonddb(genes,pathdb,pathgenes)
            outpath=self.diamondblast(pathdb,genes,self.hmms,self.outputdirectory)
        return(outpath,outpathdb)

    def diamonddb(self,genes,pathdb,pathgenes):
        cmddb = "diamond makedb --in "+pathgenes+" -d "+pathdb
        os.system(cmddb)
        return(pathdb)

    def diamondblast(self,pathdb,genes,hmms,pathblast):
        cmdblast = "diamond blastp -d "+pathdb+".dmnd"+" -q "+self.hmms+" -o "+pathblast
        os.system(cmdblast)
        #return(outpath)


    def diamonddb_stdin(self,sequences, db_path,  max_target_seqs=20, evalue=1e-05, min_pct_id=90):

        diamond = Diamond(db_path, run=run_quiet, progress=progress_quiet)
        diamond.max_target_seqs = max_target_seqs
        diamond.evalue = evalue
        diamond.min_pct_id = min_pct_id

        expected_output = diamond.makedb_stdin(sequence)
        return expected_output

    def get_hit_diamond(self,sequence,pathdb, max_target_seqs=20, evalue=1e-05, min_pct_id=90):
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
        return diamond_output

    def corre(self,pathblast):
            with open(pathblast, mode='r') as file:
                lines=file.readlines()
                for line in lines:
                    listeline=line.split("\t")
                    bestscore = 0
                    i=0
                    for line in lines:
                        listeline=line.split("\t")
                        score = listeline[11].rstrip()
                        if i==3:
                            bacgene=os.path.basename(pathblast)
                            hmmgene=listeline[0]
                            return(hmmgene,bacgene)
                            break
                        if float(score) > float(bestscore):
                            bestscore=score
                            bestname=listeline[0]
                        if bestname==listeline[0]:
                            i+=1


    def refundgenes(self,pathgenes,pathdb):
        fasta=getfasta(pathgenes)
        i=0
        sequencetoblast=""
        while i < len(fasta.ids):
            if fasta.sequences[i].count('-') < len(fasta.sequences[i])*0.5:
                self.diamonddb_stdin(sequencetoblast+">"+fasta.ids[i]+"\n"+fasta.sequences[i],pathdb)

                sequencetoblast=sequencetoblast+">"+fasta.ids[i]+"\n"+fasta.sequences[i]+"\n"
            i+=1
        return(sequencetoblast)


    def trie_blast(self,outpath):
        """ return dictonnary with the name of Genes from given sequence corresponding with the one from anvi'o"""
        blastfiles=os.listdir(outpath)
        dicocorres= {}
        for blastfile in blastfiles:
            pathblastfile=os.path.join(outpath,blastfile)
            if os.stat(pathblastfile).st_size == 0:
                os.remove(pathblastfile)
                continue
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
                    leveltaxo1 = re.sub(r'sp[1-9]', '', leveltaxo)
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
        fasta = lvl.getfasta(pathfa)

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
        self.files = ['VERSION', 'ar122_msa_individual_genes.tar.gz', 'ar122_msa_marker_info.tsv','bac120_msa_individual_genes.tar.gz','bac120_msa_marker_info.tsv', 'VERSION']


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
