"""Interface to mmseqs2"""

import anvio
import anvio.filesnpaths as filesnpaths
import anvio.terminal as terminal
import anvio.utils as utils

from anvio.errors import ConfigError

import os

from multiprocessing import cpu_count

__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2020, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Samuel Miller"
__email__ = "samuelmiller10@gmail.com"

class MMseqs2:
    """ Parent class for MMseqs2 usage

    Parameters
    ==========
    program_name: str, mmseqs
        The program name that gives access to MMseqs2 functionality
    num_threads: int, multiprocessing.cpu_count()
        Used to multithread various mmseqs subcommands, by default equal to available cores
    verbosity_level: int, 3
        Level of MMseqs2 reporting: 0=nothing, 1: +errors, 2: +warnings, 3: +info
    top_output_dir: str, anvio.filesnpaths.get_temp_directory_path()
        Top directory of MMseqs2 output, containing subdirectories for databases
    tmp_dir: str, None
        Path to temp file directory, by default, <top_output_dir>/MMSEQS_TMP
    """
    def __init__(
        self,
        program_name='mmseqs',
        num_threads=1,
        verbosity_level=3,
        top_output_dir=filesnpaths.get_temp_directory_path(),
        tmp_dir=None,
        log_file=None,
        quiet_anvio=False,
        run=None,
        progress=None,
        skip_sanity_check=False):

        self.program_name = program_name
        self.tested_versions = ['e1a1c1226ef22ac3d0da8e8f71adb8fd2388a249']

        self.num_threads = num_threads
        self.verbosity_level = verbosity_level
        self.top_output_dir = top_output_dir
        self.tmp_dir = os.path.join(
            self.top_output_dir, "MMSEQS_TMP") if not tmp_dir else tmp_dir
        if not os.path.exists(self.tmp_dir):
            os.mkdir(self.tmp_dir)
        self.log_file_path = os.path.join(
            self.top_output_dir, '00_log') if not log_file else log_file

        self.quiet_anvio = quiet_anvio
        self.run = run or terminal.Run(verbose=(not self.quiet_anvio))
        self.progress = progress or terminal.Progress(not self.quiet_anvio)

        if not skip_sanity_check:
            self.sanity_check()


        # self.run.warning("Anvi'o will now use 'MMseqs2' "
        #                  "by Steinegger et al. (DOI: 10.1038/nbt.3988). "
        #                  "If you publish your findings, "
        #                  "please do not forget to properly credit their work.",
        #                  lc='green', header="CITATION")

        # UNCOMMENT
        # utils.is_program_exists('mmseqs')


    def sanity_check(self):
        """ Executes rudimentary checks

        Parameters
        ==========
        N\A

        Returns
        =======
        N\A
        """

        self.check_programs()

        if self.num_threads not in range(1, cpu_count() + 1):
            raise ConfigError("The number of threads cannot equal %d and must lie between 1 and %d"
                              % self.num_threads, cpu_count())

        if self.verbosity_level not in (0, 1, 2, 3):
            raise ConfigError("Verbosity level must be 0, 1, 2, or 3: "
                              "0=nothing, 1: +errors, 2: +warnings, 3: +info")

        filesnpaths.is_output_file_writable(self.log_file_path)


    def dereplicate(self, fasta_file_paths, output_name=None, tmp_dir=None):
        """ Wrapper for dereplication using easy-linclust

        Parameters
        ==========
        fasta_file_paths: str or list of str
            All sequences in one or more input files will be clustered
        output_name: str, None
            Filename prefix of output files, by default, first FASTA file prefix + "_DEREP"

        Returns
        =======
        fasta_all_seqs_path: str
            Path to FASTA file of all sequences from each cluster
            with a header line for the representative sequence
        tsv_cluster: str
            Path to TSV file of cluster membership
        fasta_rep_seq_path: str
            Path to FASTA file of representative sequence from each cluster
        """

        if type(fasta_file_paths) == str:
            fasta_file_paths = [fasta_file_paths]

        if output_name is None:
            output_name = os.path.basename(os.path.splitext(fasta_file_paths[0])[0]) + "_DEREP"
        output_path_prefix = os.path.join(self.top_output_dir, output_name)

        cmd_line = [self.program_name,
                    'easy-cluster',
                    *fasta_file_paths,
                    output_path_prefix,
                    os.path.join(self.tmp_dir, 'TMP-' + output_name),
                    '-c', '1',
                    '--cov-mode', '2', # Threshold based on coverage of QUERY SEQUENCE
                    '--cluster-mode', '2', # Greedy clustering by sequence length
                    '--min-seq-id', '1',
                    '--threads', self.num_threads]

        self.progress.new('Processing')
        self.progress.update('Dereplicating sequences...')
        self.check_return_value(
            utils.run_command(cmd_line, self.log_file_path, remove_log_file_if_exists=False))
        self.progress.end()

        return (output_path_prefix + "_all_seqs.fasta",
                output_path_prefix + "_cluster.tsv",
                output_path_prefix + "_rep_seq.fasta")


    def run_cluster_pipeline(self, fasta_file_paths, cov_thresh=0.8, min_seq_id=0.98):
        """ Procedure for clustering and recovery of seed sequences and alignment information

        Parameters
        ==========
        fasta_file_paths: str or list of str
            All sequences in one or more input files will be clustered
        cov_thresh: float, 0.8
            Minimum proportion of aligned residues needed for clustering
        min_seq_id: float, 0.98
            Minimum identity of alignment needed for clustering

        Returns
        =======
        fasta_file_path: str
            Cluster seed sequences, the longest sequence in each cluster
        align_file_path: str
            Tabular m8 file of alignment information
        """

        # print(self.top_output_dir)
        seq_db_path = self.run_createdb(fasta_file_paths)
        cluster_db_path = self.run_cluster(seq_db_path)
        align_db_path = self.run_align(seq_db_path, seq_db_path, cluster_db_path)
        align_file_path = self.run_convertalis(seq_db_path, seq_db_path, align_db_path)
        cluster_tsv_path = self.run_createtsv(seq_db_path, seq_db_path, cluster_db_path)
        subset_db_path = self.run_createsubdb(cluster_db_path, seq_db_path, subset_db_name='REPSEQDB')
        fasta_file_path = self.run_result2flat(seq_db_path, seq_db_path, subset_db_path)
        return fasta_file_path, align_file_path


    def check_programs(self, quiet=False):
        utils.is_program_exists(self.program_name)

        output, ret_code = utils.get_command_output_from_shell('%s -h' % self.program_name)
        try:
            version_found = output.split(b'\n')[5].split()[2].decode()
            if not quiet:
                self.run.info(
                    "%s version found" % self.program_name, version_found, mc='green', nl_after=1)
        except:
            version_found = 'Unknown'
            self.run.warning(
                "Anvi'o failed to learn the version of %s installed on this system :/"
                % self.program_name)

        if version_found not in self.tested_versions:
            self.run.warning(
                "The version of %s installed on your system ('%s') "
                "is not one of those that we tested its anvi'o driver with. "
                "Anvi'o will continue to try to run everything as if this didn't happen. "
                "If you see this warning but everything works fine, "
                "let us know so we can include this version number "
                "in the list of 'tested' version numbers. "
                "If you see an unexpected error, please consider installing one of these versions "
                "of MMseqs2 (and again please let us know anyway so we can address it for later): "
                "'%s'" % (self.program_name, version_found, ', '.join(list(self.tested_versions))))

        self.installed_version = version_found


    def check_return_value(self, ret_value):
        if ret_value:
            raise ConfigError(
                "The last call didn't work."
                "Please look at the log file (%s) to see what went wrong."
                "The error could be due to the version of MMseqs2 that you are using."
                "We have tested the following versions: %s."
                "You can learn which version you have on your system by typing '%s -h'."
                "You can try changing the version of MMseqs2 by downloading from %s"
                % (self.log_file_path,
                   self.tested_versions,
                   self.program_name,
                   'https://github.com/soedinglab/MMseqs2/releases'))


    def run_createdb(self, fasta_file_paths, seq_db_name='INPUTDB', seq_db_subdir='INPUTDB'):
        """ Runs createdb subcommand, creating a seqDB

        Parameters
        ==========
        fasta_file_paths: str or list of str
            FASTA file input is transformed into a database
        seq_db_name: str, 'INPUTDB'
            Filename of the output database
        seq_db_subdir: str, 'INPUTDB'
            Subdirectory relative to top mmseqs output directory

        Returns
        =======
        seq_db_path: str
        """

        if type(fasta_file_paths) == str:
            fasta_file_paths = [fasta_file_paths]

        seq_db_dir = os.path.join(self.top_output_dir, seq_db_subdir)
        filesnpaths.check_output_directory(seq_db_dir)
        filesnpaths.gen_output_directory(seq_db_dir)

        seq_db_path = os.path.join(seq_db_dir, seq_db_name)

        self.run.info('[MMseqs2 createdb] Input FASTA file paths', ' ; '.join(fasta_file_paths))
        self.run.info('[MMseqs2 createdb] Output sequence database', seq_db_path)
        self.run.info('[MMSeqs2 createdb] Log file path', self.log_file_path, nl_after=True)

        cmd_line = [self.program_name,
                    'createdb',
                    *fasta_file_paths,
                    seq_db_path,
                    '-v', self.verbosity_level]

        self.progress.new('Processing')
        self.progress.update('Creating sequence database from FASTA input...')
        self.check_return_value(
            utils.run_command(cmd_line, self.log_file_path, remove_log_file_if_exists=False))
        self.progress.end()

        return seq_db_path


    def run_cluster(
        self,
        seq_db_path,
        cluster_db_name='CLUSTERDB',
        cluster_db_subdir='CLUSTERDB',
        cov_thresh=0.8,
        min_seq_id=0.98):
        """ Runs cluster subcommand, creating a clusterDB

        Parameters
        ==========
        seq_db_path: str
            mmseqs database of sequences to be clustered
        cluster_db_name: str, 'CLUSTERDB'
            Name of cluster database output
        cluster_db_subdir: str, 'CLUSTERDB'
            Output subdirectory relative to top mmseqs output directory
        cov_thresh: float, 0.8
            Minimum proportion of aligned residues needed for clustering
        min_seq_id: float, 0.98
            Minimum identity of alignment needed for clustering

        Returns
        =======
        cluster_db_dir: str
            Path to the cluster database, which is actually a set of files, e.g. CLUSTERDB.0
        """

        cluster_db_dir = os.path.join(self.top_output_dir, cluster_db_subdir)
        filesnpaths.gen_output_directory(cluster_db_dir)
        cluster_tmp_dir = os.path.join(cluster_db_dir, 'TMP')

        cluster_db_path = os.path.join(cluster_db_dir, cluster_db_name)

        self.run.info('[MMseqs2 cluster] Input sequence database', seq_db_path)
        self.run.info('[MMseqs2 cluster] Output cluster database', cluster_db_path)
        self.run.info('[MMseqs2 cluster] Coverage threshold', cov_thresh)
        self.run.info('[MMseqs2 cluster] Sequence identity threshold', min_seq_id)
        self.run.info('[MMseqs2 cluster] Number of threads', self.num_threads)
        self.run.info('[MMseqs2 cluster] Log file path', self.log_file_path, nl_after=True)

        cmd_line = [self.program_name,
                    'cluster',
                    seq_db_path,
                    cluster_db_path,
                    cluster_tmp_dir,
                    '--cluster-mode', 2, # Greedy clustering by sequence length
                    '--cov-mode', 2, # Threshold based on coverage of QUERY SEQUENCE
                    '-c', cov_thresh,
                    '--min-seq-id', min_seq_id,
                    '--threads', self.num_threads,
                    '-v', self.verbosity_level]

        self.progress.new('Processing')
        self.progress.update('Creating cluster database from sequence input...')
        self.check_return_value(
            utils.run_command(cmd_line, self.log_file_path, remove_log_file_if_exists=False))
        self.progress.end()

        return cluster_db_path


    def run_align(
        self,
        query_db_path,
        target_db_path,
        result_db_path,
        align_db_name='ALIGNDB',
        align_db_subdir=None):
        """ Runs align subcommand, creating an alignDB

        Parameters
        ==========
        query_db_path: str
            mmseqs database of query sequences
        target_db_path: str
            mmseqs database of target sequences
        result_db_path: str
            mmseqs "result" database of "prefiltered" sequences, e.g., a cluster database
        align_db_name: str, 'ALIGNDB'
            Name of database of aligned sequence database output
        align_db_subdir: str, None
            Output subdirectory relative to top mmseqs output directory

        Returns
        =======
        align_db_path: str
            Path to the aligned sequence database, which is actually a set of files, e.g. ALIGNDB.0
        """

        if align_db_subdir is None:
            align_db_subdir = os.path.dirname(result_db_path)
        align_db_dir = os.path.join(self.top_output_dir, align_db_subdir)
        filesnpaths.gen_output_directory(align_db_dir)

        align_db_path = os.path.join(align_db_dir, align_db_name)

        self.run.info('[MMseqs2 align] Query sequence database', query_db_path)
        self.run.info('[MMseqs2 align] Target sequence database', target_db_path)
        self.run.info('[MMseqs2 align] Result database', result_db_path)
        self.run.info('[MMseqs2 align] Output alignment database', align_db_path)
        self.run.info('[MMseqs2 align] Number of threads', self.num_threads)
        self.run.info('[MMseqs2 align] Log file path', self.log_file_path, nl_after=True)

        cmd_line = [self.program_name,
                    'align',
                    query_db_path,
                    target_db_path,
                    result_db_path,
                    align_db_path,
                    '--threads', self.num_threads,
                    '-v', self.verbosity_level]

        self.progress.new('Processing')
        self.progress.update('Creating alignment database from "prefiltered" sequence result...')
        self.check_return_value(
            utils.run_command(cmd_line, self.log_file_path, remove_log_file_if_exists=False))
        self.progress.end()

        return align_db_path


    def run_convertalis(
        self,
        query_db_path,
        target_db_path,
        align_db_path,
        align_file_name=None):
        """ Runs convertalis subcommand, creating a tabular m8 file of alignment information

        Parameters
        ==========
        query_db_path: str
            mmseqs database of query sequences
        target_db_path: str
            mmseqs database of target sequences
        align_db_path: str
            mmseqs database of sequence alignments
        align_file_name: str, None
            Name not including extension of output file, by default the alignment database name

        Returns
        =======
        align_file_path: str
            Placed in the directory containing the alignment database
        """

        if align_file_name is None:
            align_file_path = os.path.join(
                os.path.dirname(align_db_path), os.path.basename(align_db_path) + '.m8')
        else:
            align_file_path = os.path.join(os.path.dirname(align_db_path), align_file_name + '.m8')

        self.run.info('[MMseqs2 convertalis] Query sequence database', query_db_path)
        self.run.info('[MMseqs2 convertalis] Target sequence database', target_db_path)
        self.run.info('[MMseqs2 convertalis] Alignment database', align_db_path)
        self.run.info('[MMseqs2 convertalis] Output alignment file', align_file_path)
        self.run.info('[MMseqs2 convertalis] Number of threads', self.num_threads)
        self.run.info('[MMseqs2 convertalis] Log file path', self.log_file_path, nl_after=True)

        cmd_line = [self.program_name,
                    'convertalis',
                    query_db_path,
                    target_db_path,
                    align_db_path,
                    align_file_path,
                    '--format-output', 'query,target,qseq,pident,alnlen,mismatch,gapopen,'
                                       'qstart,qend,tstart,tend,evalue,bits',
                    '--search-type', 3]

        self.progress.new('Processing')
        self.progress.update('Recovering alignment information from database...')
        self.check_return_value(
            utils.run_command(cmd_line, self.log_file_path, remove_log_file_if_exists=False))
        self.progress.end()

        return align_file_path


    def run_createtsv(self, query_db_path, target_db_path, result_db_path, tsv_file_name=None):
        """ Run createtsv subcommand

        Parameters
        ==========
        query_db_path: str
            mmseqs database of query sequences
        target_db_path: str
            mmseqs database of target sequences
        result_db_path: str
            mmseqs "result" database of "prefiltered" sequences, e.g., a cluster database
        tsv_file_name: str, None
            Name not including extension of output file, by default the result database name

        Returns
        =======
        tsv_file_path: str
            Columns in output are
            query ID, target ID, ungapped score, diagonal (match position in query - target seq)
        """

        if tsv_file_name is None:
            tsv_file_path = os.path.join(
                os.path.dirname(result_db_path), os.path.basename(result_db_path) + '.tsv')
        else:
            tsv_file_path = os.path.join(os.path.dirname(result_db_path), tsv_file_name + '.tsv')

        self.run.info('[MMseqs2 createtsv] Query sequence database', query_db_path)
        self.run.info('[MMseqs2 createtsv] Target sequence database', target_db_path)
        self.run.info('[MMseqs2 createtsv] Result database directory', result_db_path)
        self.run.info('[MMseqs2 createtsv] Output tsv file', tsv_file_path)
        self.run.info('[MMseqs2 createtsv] Number of threads', self.num_threads)
        self.run.info('[MMseqs2 createtsv] Log file path', self.log_file_path, nl_after=True)

        cmd_line = [self.program_name,
                    'createtsv',
                    query_db_path,
                    target_db_path,
                    result_db_path,
                    tsv_file_path,
                    '--threads', self.num_threads,
                    '-v', self.verbosity_level]

        self.progress.new('Processing')
        self.progress.update('Producing tsv representation of "result" database...')
        self.check_return_value(
            utils.run_command(cmd_line, self.log_file_path, remove_log_file_if_exists=False))
        self.progress.end()

        return tsv_file_path


    def run_createsubdb(
        self,
        subset_source_db_path,
        result_db_path,
        subset_db_name='SUBSETDB'):
        """ Run createsubdb subcommand

        Parameters
        ==========
        subset_source_db_path: str
            mmseqs database used to subset sequences in other input database
        result_db_path: str
            Sequences in this database are subsetted
        subset_db_name: str, 'SUBSETDB'
            Name of output database

        Returns
        =======
        subset_db_path: str
            Path to subsetted database, placed in the directory of the database used to subset
        """

        subset_db_path = os.path.join(os.path.dirname(subset_source_db_path), subset_db_name)

        self.run.info('[MMseqs2 createsubdb] Input subset source database', subset_source_db_path)
        self.run.info('[MMseqs2 createsubdb] Input result database', result_db_path)
        self.run.info('[MMseqs2 createsubdb] Output subset database', subset_db_path)
        self.run.info('[MMseqs2 createsubdb] Log file path', self.log_file_path, nl_after=True)

        cmd_line = ['mmseqs',
                    'createsubdb',
                    subset_source_db_path,
                    result_db_path,
                    subset_db_path,
                    '-v', self.verbosity_level]

        self.progress.new('Processing')
        self.progress.update('Subsetting using %s ...' % os.path.basename(subset_source_db_path))
        self.check_return_value(
            utils.run_command(cmd_line, self.log_file_path, remove_log_file_if_exists=False))
        self.progress.end()

        return subset_db_path


    def run_result2flat(self, query_db_path, target_db_path, result_db_path, fasta_file_name=None):
        """ Run result2flat subcommand

        Parameters
        ==========
        query_db_path: str
            mmseqs database of query sequences
        target_db_path: str
            mmseqs database of target sequences
        result_db_path: str
            mmseqs "result" database to be converted to FASTA output
        fasta_file_name: str, None
            Name of FASTA output without extension, by default the name of the "result" database

        Returns
        =======
        fasta_file_path: str
            FASTA output is placed in the "result" database subdirectory
        """

        if fasta_file_name is None:
            fasta_file_name = os.path.basename(result_db_path)
        fasta_file_path = os.path.join(
            os.path.dirname(result_db_path), fasta_file_name + '.fasta')

        self.run.info('[MMseqs2 result2flat] Input query sequence database', query_db_path)
        self.run.info('[MMseqs2 result2flat] Input target sequence database', target_db_path)
        self.run.info('[MMseqs2 result2flat] Input result database', result_db_path)
        self.run.info('[MMseqs2 result2flat] Output fasta file path', fasta_file_path)
        self.run.info('[MMseqs2 result2flat] Log file path', self.log_file_path, nl_after=True)

        cmd_line = [self.program_name,
                    'result2flat',
                    query_db_path,
                    target_db_path,
                    result_db_path,
                    fasta_file_path,
                    '-v', self.verbosity_level]

        self.progress.new('Processing')
        self.progress.update('Generating FASTA file from database, %s ...'
                             % os.path.basename(result_db_path))
        self.check_return_value(
            utils.run_command(cmd_line, self.log_file_path, remove_log_file_if_exists=False))
        self.progress.end()

        return fasta_file_path

# mmseqs = MMseqs2()
# mmseqs.run_cluster_pipeline('/Users/sammiller/Documents/work/mmseqs2_test/FASTA/test.fasta')