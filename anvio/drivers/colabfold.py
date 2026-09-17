"""Interface to ColabFold (AlphaFold2 protein structure prediction).

ColabFold exposes two command-line programs:

    - colabfold_search : CPU-heavy multiple sequence alignment (MSA) generation against a
                         local database. Produces `.a3m` files.
    - colabfold_batch  : the predictor. In online mode it generates the MSA using the public
                         MMseqs2 server AND predicts the structure in one call. Given a
                         directory of `.a3m` files it skips the MSA and only predicts.

So the MSA step runs either online (public server) or locally (colabfold_search against a
downloaded database). This driver handles both. ColabFold is not assumed to be on the $PATH:
if the user provides a conda environment name, every command is prefixed with
`conda run -n <NAME>`; otherwise the colabfold programs are expected to be directly callable
(e.g. on the $PATH).

This driver is a thin wrapper around the two programs. Locating and parsing the per-gene output
files (best model PDB and confidence scores) happens in anvio.structureops, which orchestrates
a single batched run over all genes of interest.
"""

import os
import glob
import shlex

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['Florian_Trigodet']


run = terminal.Run()
progress = terminal.Progress()


class ColabFold:
    def __init__(self, args, run=run, progress=progress, skip_sanity_check=False):
        self.run = run
        self.progress = progress
        self.args = args

        A = lambda x, t: t(args.__dict__[x]) if x in args.__dict__ and args.__dict__[x] is not None else None
        null = lambda x: x

        self.conda_env = A('colabfold_conda_env', str)
        self.use_msa_server = A('colabfold_msa_server', bool)
        self.local_db = A('colabfold_db', null)
        self.num_models = A('num_models', int)
        self.num_recycle = A('num_recycle', int)
        self.amber = A('amber', bool)
        self.additional_params = A('colabfold_additional_parameters', str)

        # checkpoint flags that split the run into its MSA (colabfold_search) and prediction
        # (colabfold_batch) halves. Only one may be set (see sanity_check).
        self.only_msa = A('only_msa', bool)
        self.only_predict = A('only_predict', bool)

        # anvi'o's --num-threads is authoritative: it is always forwarded to the colabfold programs,
        # overriding their own defaults (e.g. colabfold_search defaults to 64 threads). So if the user
        # leaves anvi'o's default of 1, colabfold_search is run with a single thread too.
        self.num_threads = A('num_threads', int)

        self.citation = ("ColabFold: Mirdita et al. 2022 (doi:10.1038/s41592-022-01488-1); "
                         "AlphaFold2: Jumper et al. 2021 (doi:10.1038/s41586-021-03819-2)")
        self.web = "https://github.com/sokrypton/ColabFold"

        # every colabfold command is run inside the user's conda environment. when no environment
        # name is given, the programs are expected to be directly available (e.g. on the $PATH).
        self.command_prefix = ['conda', 'run', '--no-capture-output', '-n', self.conda_env] if self.conda_env else []

        # the MSA is generated either locally (colabfold_search against self.local_db) or online
        # (colabfold_batch queries the public MMseqs2 server)
        self.msa_source = 'local' if self.local_db else 'server'

        # --only-predict resumes from a directory of .a3m files an earlier --only-msa run produced, so
        # its prediction is always fed a local MSA directory regardless of whether an MSA source was
        # given on the command line
        if self.only_predict:
            self.msa_source = 'local'

        self.search_program = 'colabfold_search'
        self.batch_program = 'colabfold_batch'

        if not skip_sanity_check:
            self.sanity_check()


    def sanity_check(self):
        if self.num_models is not None and not (1 <= self.num_models <= 5):
            raise ConfigError("ColabFold can predict between 1 and 5 models per protein, but you asked for %d with "
                              "--num-models. Please pick a value in that range." % self.num_models)

        if self.only_msa and self.only_predict:
            raise ConfigError("You asked for both --only-msa and --only-predict, but these are the two halves of a "
                              "single split ColabFold run and are mutually exclusive. Run --only-msa first (it "
                              "generates the MSAs on, say, a CPU node), then --only-predict (it predicts structures "
                              "from those MSAs on a GPU node).")

        # --only-predict resumes from a local MSA checkpoint on disk: it needs only the predictor (and a
        # GPU), not an MSA source, since the MSAs already exist.
        if self.only_predict:
            self.check_program(self.batch_program)
            self.check_gpu()
            return

        # from here on, the MSA step will run (either as part of a full run, or standalone with
        # --only-msa), so an MSA source is required
        if self.local_db and self.use_msa_server:
            raise ConfigError("You asked ColabFold to use both a local database (--colabfold-db) and the "
                              "public MSA server (--colabfold-msa-server) for the MSA step. Please pick only one.")

        if not self.local_db and not self.use_msa_server:
            raise ConfigError("ColabFold needs to know how to generate the multiple sequence alignment (MSA). "
                              "Please use either --colabfold-msa-server to query the public MMseqs2 server (fine "
                              "for a handful of sequences), or --colabfold-db to point to a local ColabFold "
                              "database directory (recommended for many sequences).")

        if self.only_msa and self.msa_source != 'local':
            raise ConfigError("--only-msa splits the ColabFold run at the MSA step, which anvi'o only runs locally "
                              "(with colabfold_search against --colabfold-db). It cannot be combined with "
                              "--colabfold-msa-server, because the public server generates the MSA and predicts the "
                              "structure in a single step that cannot be split. Please provide a local ColabFold "
                              "database with --colabfold-db instead.")

        if self.local_db:
            if not filesnpaths.is_file_exists(self.local_db, dont_raise=True) or not os.path.isdir(self.local_db):
                raise ConfigError("The ColabFold database path you provided with --colabfold-db does not seem to be a "
                                  "directory that exists: '%s'. This should be the directory you set up with ColabFold's "
                                  "`setup_databases.sh` script." % self.local_db)
            self.check_program(self.search_program)

        # --only-msa stops after the MSA step, so it needs neither the predictor nor a GPU
        if not self.only_msa:
            self.check_program(self.batch_program)
            self.check_gpu()


    def check_program(self, program):
        """Confirm a colabfold program is callable (inside the conda env, if one was given)"""

        log_file_path = filesnpaths.get_temp_file_path()

        try:
            ret_val = utils.run_command(self.command_prefix + [program, '--help'], log_file_path)
        except ConfigError:
            ret_val = -1

        if ret_val != 0:
            env_msg = (" inside the conda environment '%s'" % self.conda_env) if self.conda_env else \
                      (". Since you did not provide a --colabfold-conda-env, anvi'o looked for it on your $PATH")
            raise ConfigError("Anvi'o could not run ColabFold's '%s' program%s. If ColabFold is installed in a conda "
                              "environment, provide its name with --colabfold-conda-env so anvi'o can run it via "
                              "`conda run -n <NAME>`. The command anvi'o tried to run is in this log file, in case it "
                              "helps: %s" % (program, env_msg, log_file_path))

        os.remove(log_file_path) if filesnpaths.is_file_exists(log_file_path, dont_raise=True) else None


    def check_gpu(self):
        """Best-effort, pre-run check of whether ColabFold will see a GPU, so the user is not surprised
        by a slow CPU-only run.

        We ask the JAX runtime (the same one colabfold_batch uses) which devices it can see. This is a
        best-effort check: if we can't run it (e.g. JAX is not importable from the python we reach), we
        say so and rely on the post-run log check (see report_gpu_usage_from_log) instead. It never
        raises.
        """

        # a one-liner that asks JAX what devices it sees and prints a token we can parse
        probe = ("import jax; "
                 "devices = jax.devices(); "
                 "gpus = [d for d in devices if d.platform != 'cpu']; "
                 "print('ANVIO_GPU_CHECK|%d|%s' % (len(gpus), ','.join(sorted({d.platform for d in devices}))))")

        try:
            output = utils.run_command_and_get_output(self.command_prefix + ['python', '-c', probe], raise_on_error=False)
        except Exception:
            output = ''

        token = None
        for line in (output or '').splitlines():
            if line.startswith('ANVIO_GPU_CHECK|'):
                token = line.strip()
                break

        if token is None:
            self.run.warning("Anvi'o could not check ahead of time whether ColabFold will use a GPU (it was not able "
                             "to query the JAX runtime that ColabFold uses). No worries: the structure prediction will "
                             "still run, and anvi'o will tell you at the end whether a GPU or the CPU was actually "
                             "used.", header="GPU CHECK: COULD NOT DETERMINE", lc="yellow")
            return

        _, num_gpus, platforms = token.split('|')

        if int(num_gpus) > 0:
            self.run.info_single("ColabFold's runtime sees %s GPU device(s) (%s), so the structure prediction should "
                                 "run on the GPU. Nice." % (num_gpus, platforms), mc='green', nl_before=1, nl_after=1)
        else:
            self.run.warning("ColabFold's runtime does not see a GPU, so the structure prediction will run on the CPU. "
                             "This works, but it is *much* slower -- a single protein can take many minutes. If you do "
                             "have a GPU, make sure ColabFold (and its JAX installation) can see it before running. "
                             "Press CTRL+C now if you would like to stop.", header="GPU CHECK: NO GPU DETECTED")


    def report_gpu_usage_from_log(self, log_file_path):
        """After a run, report from ColabFold's own log whether it used a GPU or fell back to the CPU.

        ColabFold logs `no GPU detected, will be using CPU` when it falls back to the CPU; the absence
        of that message means it used the GPU. This is best-effort and never raises.
        """

        if not filesnpaths.is_file_exists(log_file_path, dont_raise=True):
            return

        try:
            with open(log_file_path) as log_file:
                log_content = log_file.read().lower()
        except Exception:
            return

        if 'no gpu detected' in log_content:
            self.run.warning("Heads up: ColabFold reported that no GPU was detected, so this run used the CPU. That is "
                             "why it may have been slow. If you expected a GPU to be used, check that ColabFold and its "
                             "JAX installation can see your GPU.", header="COLABFOLD USED THE CPU")
        else:
            self.run.info_single("ColabFold did not report falling back to the CPU, so the GPU was used for this run.",
                                 mc='green', nl_before=1, nl_after=1)


    def run_search(self, fasta_path, msa_dir, log_file_path):
        """Generate MSAs locally with colabfold_search. Produces .a3m files in msa_dir.

        Kept as a distinct step so a future PR can expose an --only-msa / --only-predict checkpoint.
        """

        cmd_line = self.command_prefix + [self.search_program, fasta_path, self.local_db, msa_dir]

        if self.num_threads is not None:
            cmd_line += ['--threads', str(self.num_threads)]

        self.run.info('ColabFold MSA source', 'local database (%s)' % self.local_db)
        self.run.info('ColabFold search cmd', ' '.join([str(x) for x in cmd_line]), quiet=(not anvio.DEBUG))

        self.progress.new('ColabFold')
        self.progress.update('Generating MSAs locally with %s (this is CPU-heavy) ...' % self.search_program)
        # append so we don't clobber a log the prediction step will also write to (see process())
        ret_val = utils.run_command(cmd_line, log_file_path, remove_log_file_if_exists=False)
        self.progress.end()

        if ret_val != 0 or not len(glob.glob(os.path.join(msa_dir, '*.a3m'))):
            raise ConfigError("Something went wrong while ColabFold was generating MSAs with '%s'. Please take a look at "
                              "the log file to find out what happened: %s" % (self.search_program, log_file_path))

        return msa_dir


    def run_batch(self, input_path, out_dir, log_file_path):
        """Run colabfold_batch to predict structures.

        When self.msa_source is 'server', input_path is the query FASTA and the MSA is generated on the
        fly via the public server. When it is 'local', input_path is the directory of .a3m files
        produced by run_search.
        """

        # colabfold_batch has no threads flag: the prediction is GPU-bound and manages its own device
        # parallelism. The only case where a CPU thread count matters is a CPU-only run, which respects
        # OMP_NUM_THREADS, so we set it via `env` (it has no effect when a GPU is used).
        thread_env = ['env', 'OMP_NUM_THREADS=%d' % self.num_threads] if self.num_threads is not None else []

        cmd_line = self.command_prefix + thread_env + [self.batch_program, input_path, out_dir]

        if self.msa_source == 'server':
            cmd_line += ['--msa-mode', 'mmseqs2_uniref_env']

        if self.num_models is not None:
            cmd_line += ['--num-models', str(self.num_models)]

        if self.num_recycle is not None:
            cmd_line += ['--num-recycle', str(self.num_recycle)]

        if self.amber:
            cmd_line += ['--amber', '--num-relax', '1']

        if self.additional_params:
            # shlex.split (not str.split) so quoted values with spaces, e.g. paths, survive intact
            cmd_line += shlex.split(self.additional_params)

        self.run.info('ColabFold predict cmd', ' '.join([str(x) for x in cmd_line]), quiet=(not anvio.DEBUG))

        self.progress.new('ColabFold')
        self.progress.update('Predicting structures with %s (this is GPU-heavy) ...' % self.batch_program)
        # append so we don't clobber the MSA step's log (both steps share one log; see process())
        ret_val = utils.run_command(cmd_line, log_file_path, remove_log_file_if_exists=False)
        self.progress.end()

        if ret_val != 0:
            raise ConfigError("Something went wrong while ColabFold was predicting structures with '%s'. Please take a "
                              "look at the log file to find out what happened: %s" % (self.batch_program, log_file_path))

        return out_dir


    def process(self, fasta_path, out_dir):
        """Run ColabFold end-to-end over a (multi-sequence) FASTA and return the results directory.

        Per-gene output parsing is done by the caller (anvio.structureops) via get_output_paths.
        """

        filesnpaths.gen_output_directory(out_dir, delete_if_exists=False, dont_warn=True)
        log_file_path = os.path.join(out_dir, '00_log.txt')
        msa_dir = os.path.join(out_dir, 'msas')

        # the MSA and prediction steps both append to this single log (each writing its own command
        # header + output), so we start it fresh here in case a previous run left one behind (e.g.
        # when reusing a --dump-dir). --only-predict is the exception: it resumes an --only-msa run, so
        # we keep that run's MSA log and append the prediction log to it for a complete record.
        if not self.only_predict and filesnpaths.is_file_exists(log_file_path, dont_raise=True):
            os.remove(log_file_path)

        self.run.info('Log file path', log_file_path)

        if self.only_msa:
            # generate the MSAs locally and stop; --only-predict will resume from this directory
            self.run_search(fasta_path, msa_dir, log_file_path)
            return out_dir

        if self.only_predict:
            # the MSAs were produced by an earlier --only-msa run; skip straight to prediction
            self.run_batch(msa_dir, out_dir, log_file_path)
        elif self.msa_source == 'local':
            self.run_search(fasta_path, msa_dir, log_file_path)
            self.run_batch(msa_dir, out_dir, log_file_path)
        else:
            self.run_batch(fasta_path, out_dir, log_file_path)

        # report from ColabFold's own log whether it used the GPU or fell back to the CPU
        self.report_gpu_usage_from_log(log_file_path)

        return out_dir


    def get_output_paths(self, out_dir, jobname):
        """Locate the best-model PDB and its scores JSON for a given jobname in a results directory.

        ColabFold names the best model `rank_001` and writes a relaxed variant when --amber is used.
        Returns a (pdb_path, scores_path) tuple, or (None, None) if the prediction is not found.
        """

        tag = 'relaxed' if self.amber else 'unrelaxed'

        pdb_hits = sorted(glob.glob(os.path.join(out_dir, '%s_%s_rank_001_*.pdb' % (jobname, tag))))
        if not pdb_hits:
            # fall back to the unrelaxed model in case relaxation was skipped for this query
            pdb_hits = sorted(glob.glob(os.path.join(out_dir, '%s_unrelaxed_rank_001_*.pdb' % jobname)))

        scores_hits = sorted(glob.glob(os.path.join(out_dir, '%s_scores_rank_001_*.json' % jobname)))

        pdb_path = pdb_hits[0] if pdb_hits else None
        scores_path = scores_hits[0] if scores_hits else None

        return pdb_path, scores_path
