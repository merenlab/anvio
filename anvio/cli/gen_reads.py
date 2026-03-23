#!/usr/bin/env python
# -*- coding: utf-8

import os
import sys
import math
import random
import bisect

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.fastalib as fastalib
import anvio.filesnpaths as filesnpaths

from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2026, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['FlorianTrigodet']
__requires__ = ['fasta']
__provides__ = ['paired-end-fastq', 'single-end-fastq']
__description__ = ("Generate synthetic sequencing reads (Illumina, PacBio HiFi, ONT) from reference "
                   "FASTA files with optional SNV injection at controlled multi-allele frequencies")


PRESETS = {
    'illumina-paired': {
        'read_type': 'paired-end',
        'read_length': 150,
        'insert_size': 450,
        'insert_size_std': 50,
        'error_rate': 0.005,
        'quality_score': '?',  # Q30
    },
    'illumina-single': {
        'read_type': 'single-end',
        'read_length': 150,
        'error_rate': 0.005,
        'quality_score': '?',  # Q30
    },
    'pacbio-hifi': {
        'read_type': 'long-distributed',
        'read_length': 15000,
        'read_length_std': 3500,
        'length_distribution': 'normal',
        'min_read_length': 5000,
        'error_rate': 0.001,
        'quality_score': 'F',  # Q37
    },
    'pacbio-clr': {
        'read_type': 'long-distributed',
        'read_length': 15000,
        'read_length_std': 8000,
        'length_distribution': 'normal',
        'min_read_length': 1000,
        'error_rate': 0.12,
        'quality_score': '.',  # Q13
    },
    'ont-r9': {
        'read_type': 'long-distributed',
        'read_length': 5000,
        'read_length_std': 4000,
        'length_distribution': 'lognormal',
        'min_read_length': 200,
        'error_rate': 0.06,
        'quality_score': '3',  # Q18
    },
    'ont-r10': {
        'read_type': 'long-distributed',
        'read_length': 8000,
        'read_length_std': 5000,
        'length_distribution': 'lognormal',
        'min_read_length': 200,
        'error_rate': 0.01,
        'quality_score': '=',  # Q28
    },
    'ont-ultralong': {
        'read_type': 'long-distributed',
        'read_length': 50000,
        'read_length_std': 40000,
        'length_distribution': 'lognormal',
        'min_read_length': 1000,
        'error_rate': 0.02,
        'quality_score': ':',  # Q25
    },
}

BASES = ['A', 'T', 'C', 'G']


def main():
    try:
        run_program()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def resolve_parameters(args):
    """Merge preset defaults with explicit CLI overrides.

    Any CLI argument that is not None takes precedence over the preset value.
    """

    params = {}

    if args.preset:
        if args.preset not in PRESETS:
            raise ConfigError(f"Unknown preset '{args.preset}'. Available presets: {', '.join(PRESETS.keys())}")
        params = dict(PRESETS[args.preset])

    # overlay explicit CLI flags (anything not None overrides preset)
    cli_overrides = {
        'read_type': args.read_type,
        'read_length': args.read_length,
        'read_length_std': args.read_length_std,
        'min_read_length': args.min_read_length,
        'coverage': args.coverage,
        'insert_size': args.insert_size,
        'insert_size_std': args.insert_size_std,
        'error_rate': args.error_rate,
        'quality_score': args.quality_score,
    }

    for key, value in cli_overrides.items():
        if value is not None:
            params[key] = value

    # length_distribution can be overridden via CLI
    if args.length_distribution is not None:
        params['length_distribution'] = args.length_distribution

    # apply remaining defaults for things not set by preset or CLI
    params.setdefault('coverage', 50)
    params.setdefault('error_rate', 0.0)
    params.setdefault('quality_score', '?')
    params.setdefault('min_read_length', 200)
    params.setdefault('length_distribution', 'normal')

    # validation
    if 'read_type' not in params:
        raise ConfigError("You need to specify either --preset or --read-type so anvi'o knows what "
                          f"kind of reads to generate. Available presets: {', '.join(PRESETS.keys())}. "
                          f"Available read types: paired-end, single-end, long-fixed, long-distributed.")

    if 'read_length' not in params:
        raise ConfigError("Read length is required. Either use a --preset or provide --read-length explicitly.")

    if params['read_type'] == 'paired-end':
        if 'insert_size' not in params:
            raise ConfigError("Paired-end mode requires --insert-size (or use a preset like 'illumina-paired').")
        params.setdefault('insert_size_std', 0)

    if params['read_type'] == 'long-distributed':
        if 'read_length_std' not in params:
            raise ConfigError("Long-distributed mode requires --read-length-std (or use a long-read preset "
                              "like 'ont-r10', 'pacbio-hifi', etc.).")

    if params['coverage'] <= 0:
        raise ConfigError("Coverage must be greater than 0.")

    if not (0 <= params['error_rate'] < 1):
        raise ConfigError("Error rate must be between 0 and 1 (exclusive).")

    if params['read_length'] <= 0:
        raise ConfigError("Read length must be greater than 0.")

    return params


def load_mutations_file(filepath, contigs):
    """Parse a mutations TSV file and return a mutation lookup structure.

    Returns
    =======
    dict : {contig_name: [(position, cumulative_freqs), ...]} sorted by position
        cumulative_freqs is [cum_A, cum_T, cum_C, cum_G] for use with random.random()
    """

    filesnpaths.is_file_exists(filepath)
    filesnpaths.is_file_tab_delimited(filepath)

    mutations = {}
    line_num = 0

    with open(filepath) as f:
        header = f.readline().strip().split('\t')
        line_num = 1

        expected_header = ['contig_name', 'position', 'freq_A', 'freq_T', 'freq_C', 'freq_G']
        if header != expected_header:
            raise ConfigError(f"The mutations file header doesn't look right. Expected columns: "
                              f"{', '.join(expected_header)}, but got: {', '.join(header)}")

        for line in f:
            line_num += 1
            fields = line.strip().split('\t')

            if len(fields) != 6:
                raise ConfigError(f"Line {line_num} in the mutations file has {len(fields)} columns "
                                  f"instead of the expected 6.")

            contig_name = fields[0]
            try:
                position = int(fields[1])
                freqs = [float(fields[i]) for i in range(2, 6)]
            except ValueError:
                raise ConfigError(f"Line {line_num} in the mutations file has non-numeric values "
                                  f"where numbers were expected.")

            if contig_name not in contigs:
                raise ConfigError(f"Line {line_num} in the mutations file references contig "
                                  f"'{contig_name}', which is not in the input FASTA.")

            if position < 0 or position >= len(contigs[contig_name]):
                raise ConfigError(f"Line {line_num}: position {position} is out of bounds for "
                                  f"contig '{contig_name}' (length {len(contigs[contig_name])}).")

            freq_sum = sum(freqs)
            if abs(freq_sum - 1.0) > 0.01:
                raise ConfigError(f"Line {line_num}: base frequencies sum to {freq_sum:.4f} instead "
                                  f"of 1.0. Please make sure freq_A + freq_T + freq_C + freq_G = 1.0.")

            # normalize to exactly 1.0
            freqs = [f / freq_sum for f in freqs]

            # build cumulative frequencies for sampling
            cumulative = []
            running = 0.0
            for freq in freqs:
                running += freq
                cumulative.append(running)
            cumulative[-1] = 1.0  # avoid floating point edge case

            if contig_name not in mutations:
                mutations[contig_name] = []
            mutations[contig_name].append((position, cumulative))

    # sort by position and pre-build positions index for bisect lookups
    for contig_name in mutations:
        mutations[contig_name].sort(key=lambda x: x[0])

    # convert to {contig: (positions_list, entries_list)} for efficient bisect
    return {contig: ([m[0] for m in entries], entries) for contig, entries in mutations.items()}


def generate_random_mutations(contigs, density, num_alleles, rng):
    """Generate random multi-allele mutations at a given density.

    Returns the same structure as load_mutations_file.
    """

    mutations = {}

    for contig_name, seq in contigs.items():
        num_snvs = int(len(seq) * density)
        if num_snvs == 0:
            continue

        positions = sorted(rng.sample(range(len(seq)), min(num_snvs, len(seq))))
        contig_mutations = []

        for pos in positions:
            ref_base = seq[pos]
            if ref_base not in BASES:
                continue  # skip N or ambiguous bases

            # pick which bases will be present (always include reference)
            alt_bases = [b for b in BASES if b != ref_base]
            chosen_alts = rng.sample(alt_bases, num_alleles - 1)
            allele_bases = [ref_base] + chosen_alts

            # generate random frequencies for the chosen alleles
            weights = [rng.random() for _ in allele_bases]
            total_weight = sum(weights)
            allele_freqs = {b: w / total_weight for b, w in zip(allele_bases, weights)}

            # build frequency array in A, T, C, G order
            freqs = [allele_freqs.get(b, 0.0) for b in BASES]

            # cumulative
            cumulative = []
            running = 0.0
            for f in freqs:
                running += f
                cumulative.append(running)
            cumulative[-1] = 1.0

            contig_mutations.append((pos, cumulative))

        if contig_mutations:
            mutations[contig_name] = contig_mutations

    # convert to {contig: (positions_list, entries_list)} for efficient bisect
    return {contig: ([m[0] for m in entries], entries) for contig, entries in mutations.items()}


def apply_mutations(seq_list, start_pos, contig_name, mutations, rng):
    """Apply mutations from the mutation table to a read (as a list of chars).

    Modifies seq_list in place. Uses bisect for efficient lookup.
    """

    if contig_name not in mutations:
        return

    positions, entries = mutations[contig_name]

    left = bisect.bisect_left(positions, start_pos)
    right = bisect.bisect_right(positions, start_pos + len(seq_list) - 1)

    for i in range(left, right):
        pos, cumulative = entries[i]
        idx_in_read = pos - start_pos
        r = rng.random()
        for base_idx, cum_freq in enumerate(cumulative):
            if r < cum_freq:
                seq_list[idx_in_read] = BASES[base_idx]
                break


def apply_errors(seq_list, error_rate, rng):
    """Apply random base substitution errors to a read (as a list of chars)."""

    if error_rate <= 0:
        return

    for i in range(len(seq_list)):
        if rng.random() < error_rate:
            current = seq_list[i]
            if current in BASES:
                alternatives = [b for b in BASES if b != current]
                seq_list[i] = rng.choice(alternatives)


def write_fastq_record(fh, read_id, sequence, quality_char):
    """Write a single FASTQ record."""

    fh.write(f"@{read_id}\n{sequence}\n+\n{quality_char * len(sequence)}\n")


def generate_paired_end_reads(contig_name, contig_seq, params, mutations, rng, read_counter):
    """Generator yielding (r1_id, r1_seq, r2_id, r2_seq, new_counter) tuples."""

    read_length = params['read_length']
    insert_size = params['insert_size']
    insert_size_std = params['insert_size_std']
    error_rate = params['error_rate']
    contig_length = len(contig_seq)

    min_fragment = 2 * read_length
    if contig_length < min_fragment:
        return

    num_pairs = int(contig_length * params['coverage'] / (2 * read_length + insert_size))

    for _ in range(num_pairs):
        I = max(0, int(round(rng.gauss(insert_size, insert_size_std))))
        fragment = 2 * read_length + I

        if fragment > contig_length:
            I = 0
            fragment = 2 * read_length
            if fragment > contig_length:
                continue

        start = rng.randint(0, contig_length - fragment)

        # R1: forward
        r1_list = list(contig_seq[start:start + read_length])
        apply_mutations(r1_list, start, contig_name, mutations, rng)
        apply_errors(r1_list, error_rate, rng)

        # R2: reverse complement
        r2_start = start + read_length + I
        r2_list = list(contig_seq[r2_start:r2_start + read_length])
        apply_mutations(r2_list, r2_start, contig_name, mutations, rng)
        apply_errors(r2_list, error_rate, rng)
        r2_seq = utils.rev_comp(''.join(r2_list))

        read_counter += 1
        r1_id = f"{contig_name}_read_{read_counter}/1"
        r2_id = f"{contig_name}_read_{read_counter}/2"

        yield r1_id, ''.join(r1_list), r2_id, r2_seq, read_counter


def generate_single_end_reads(contig_name, contig_seq, params, mutations, rng, read_counter):
    """Generator yielding (read_id, read_seq, new_counter) tuples."""

    read_length = params['read_length']
    error_rate = params['error_rate']
    contig_length = len(contig_seq)

    if contig_length < read_length:
        return

    num_reads = int(contig_length * params['coverage'] / read_length)

    for _ in range(num_reads):
        start = rng.randint(0, contig_length - read_length)
        seq_list = list(contig_seq[start:start + read_length])
        apply_mutations(seq_list, start, contig_name, mutations, rng)
        apply_errors(seq_list, error_rate, rng)

        read_counter += 1
        seq = ''.join(seq_list)

        # randomly choose strand
        if rng.random() < 0.5:
            seq = utils.rev_comp(seq)

        yield f"{contig_name}_read_{read_counter}", seq, read_counter


def sample_lognormal_length(desired_mean, desired_std, rng):
    """Sample a read length from a lognormal distribution.

    Converts desired mean and std (in bp) to the mu and sigma parameters of the
    underlying normal distribution that the lognormal is based on.
    """

    # lognormal parameterization: if X ~ LogNormal(mu, sigma), then
    # E[X] = exp(mu + sigma^2/2) and Var[X] = (exp(sigma^2) - 1) * exp(2*mu + sigma^2)
    # solving for mu and sigma given desired mean M and std S:
    sigma_sq = math.log(1 + (desired_std / desired_mean) ** 2)
    mu = math.log(desired_mean) - sigma_sq / 2

    return int(rng.lognormvariate(mu, math.sqrt(sigma_sq)))


def generate_long_reads(contig_name, contig_seq, params, mutations, rng, read_counter):
    """Generator yielding (read_id, read_seq, new_counter) tuples.

    Supports both fixed-length and distributed-length modes, with normal or
    lognormal length distributions.
    """

    mean_length = params['read_length']
    error_rate = params['error_rate']
    contig_length = len(contig_seq)
    distributed = params['read_type'] == 'long-distributed'
    use_lognormal = params.get('length_distribution') == 'lognormal'

    min_length = params.get('min_read_length', 200) if distributed else mean_length

    if contig_length < min_length:
        return

    num_reads = int(contig_length * params['coverage'] / mean_length)

    for _ in range(num_reads):
        if distributed:
            if use_lognormal:
                read_length = sample_lognormal_length(mean_length, params['read_length_std'], rng)
                read_length = max(min_length, read_length)
            else:
                read_length = max(min_length, int(rng.gauss(mean_length, params['read_length_std'])))
            read_length = min(read_length, contig_length)
        else:
            read_length = min(mean_length, contig_length)

        start = rng.randint(0, contig_length - read_length)
        seq_list = list(contig_seq[start:start + read_length])
        apply_mutations(seq_list, start, contig_name, mutations, rng)
        apply_errors(seq_list, error_rate, rng)

        read_counter += 1
        seq = ''.join(seq_list)

        if rng.random() < 0.5:
            seq = utils.rev_comp(seq)

        yield f"{contig_name}_read_{read_counter}", seq, read_counter


def run_program():
    args = get_args()
    run = terminal.Run()
    progress = terminal.Progress()
    pp = terminal.pretty_print

    # resolve parameters from preset + CLI
    params = resolve_parameters(args)
    read_type = params['read_type']

    # set random seed
    seed = args.seed
    rng = random.Random(seed)

    # validate SNV injection args
    if args.mutations_file and args.snv_density:
        raise ConfigError("You can't use --mutations-file and --snv-density at the same time. "
                          "Please pick one or the other for SNV injection.")

    if args.num_alleles is not None and not args.snv_density:
        raise ConfigError("--num-alleles only makes sense with --snv-density. If you're using "
                          "--mutations-file, allele frequencies are specified per position in the file.")

    num_alleles = args.num_alleles or 2
    if num_alleles < 2 or num_alleles > 4:
        raise ConfigError("--num-alleles must be between 2 and 4.")

    # validate output
    output_prefix = args.output_file_prefix
    if not output_prefix:
        raise ConfigError("Please provide an output file prefix with -o / --output-file-prefix.")

    output_dir = os.path.dirname(os.path.abspath(output_prefix))
    filesnpaths.is_output_dir_writable(output_dir)

    # load contigs
    filesnpaths.is_file_exists(args.fasta_file)

    run.info('Input FASTA', args.fasta_file)

    contigs = {}
    fasta = fastalib.SequenceSource(args.fasta_file)
    while next(fasta):
        contigs[fasta.id] = fasta.seq.upper()
    fasta.close()

    if not contigs:
        raise ConfigError("The input FASTA file appears to be empty. Nothing to do here :(")

    run.info('Contigs loaded', f"{len(contigs)} ({pp(sum(len(s) for s in contigs.values()))} bp total)")

    # filter short contigs
    min_contig_length = params['read_length']
    if read_type == 'paired-end':
        min_contig_length = 2 * params['read_length'] + params.get('insert_size', 0)

    short_contigs = [name for name, seq in contigs.items() if len(seq) < min_contig_length]
    if short_contigs:
        for name in short_contigs:
            del contigs[name]
        run.warning(f"{len(short_contigs)} contig(s) shorter than {min_contig_length} bp were skipped: "
                    f"{', '.join(short_contigs[:5])}{'...' if len(short_contigs) > 5 else ''}",
                    header="SHORT CONTIGS SKIPPED")

    if not contigs:
        raise ConfigError("All contigs were too short for the requested read type and length. "
                          "Nothing to generate :(")

    # load or generate mutations
    mutations = {}
    if args.mutations_file:
        mutations = load_mutations_file(args.mutations_file, contigs)
        total_snv_positions = sum(len(positions) for positions, _ in mutations.values())
        run.info('Mutations file', args.mutations_file)
        run.info('SNV positions loaded', pp(total_snv_positions))
    elif args.snv_density:
        if args.snv_density <= 0 or args.snv_density >= 1:
            raise ConfigError("--snv-density must be between 0 and 1 (exclusive).")
        mutations = generate_random_mutations(contigs, args.snv_density, num_alleles, rng)
        total_snv_positions = sum(len(positions) for positions, _ in mutations.values())
        run.info('SNV density', args.snv_density)
        run.info('Num alleles per SNV', num_alleles)
        run.info('SNV positions generated', pp(total_snv_positions))
    else:
        run.info('SNV injection', 'None')

    # report parameters
    run.info('Preset', args.preset or 'None')
    run.info('Read type', read_type)
    run.info('Read length', params['read_length'])
    if read_type == 'long-distributed':
        run.info('Read length std', params['read_length_std'])
        run.info('Length distribution', params['length_distribution'])
        run.info('Min read length', params['min_read_length'])
    if read_type == 'paired-end':
        run.info('Insert size', params['insert_size'])
        run.info('Insert size std', params['insert_size_std'])
    run.info('Coverage', params['coverage'])
    run.info('Error rate', params['error_rate'])
    run.info('Quality score', params['quality_score'])
    run.info('Random seed', seed)

    # estimate total reads for progress bar
    total_bp = sum(len(s) for s in contigs.values())
    if read_type == 'paired-end':
        est_total = int(total_bp * params['coverage'] / (2 * params['read_length'] + params['insert_size'])) * 2
    else:
        est_total = int(total_bp * params['coverage'] / params['read_length'])

    # open output files
    quality_char = params['quality_score']

    file_handles = []
    if read_type == 'paired-end':
        r1_path = f"{output_prefix}-R1.fastq"
        r2_path = f"{output_prefix}-R2.fastq"
        r1_fh = open(r1_path, 'w')
        file_handles.append(r1_fh)
        r2_fh = open(r2_path, 'w')
        file_handles.append(r2_fh)
        output_files = [r1_path, r2_path]
    else:
        out_path = f"{output_prefix}.fastq"
        out_fh = open(out_path, 'w')
        file_handles.append(out_fh)
        output_files = [out_path]

    # generate reads
    read_counter = 0
    total_reads = 0

    progress.new('Generating reads', progress_total_items=est_total)

    try:
        for contig_name, contig_seq in contigs.items():
            if read_type == 'paired-end':
                for r1_id, r1_seq, r2_id, r2_seq, read_counter in generate_paired_end_reads(contig_name, contig_seq, params, mutations, rng, read_counter):
                    write_fastq_record(r1_fh, r1_id, r1_seq, quality_char)
                    write_fastq_record(r2_fh, r2_id, r2_seq, quality_char)
                    total_reads += 2
                    if total_reads % 10000 == 0:
                        progress.increment(10000)
                        progress.update(f"{pp(total_reads)} reads written")

            elif read_type in ('single-end', 'long-fixed', 'long-distributed'):
                generator = generate_long_reads if read_type in ('long-fixed', 'long-distributed') else generate_single_end_reads
                for read_id, read_seq, read_counter in generator(contig_name, contig_seq, params, mutations, rng, read_counter):
                    write_fastq_record(out_fh, read_id, read_seq, quality_char)
                    total_reads += 1
                    if total_reads % 10000 == 0:
                        progress.increment(10000)
                        progress.update(f"{pp(total_reads)} reads written")
    finally:
        for fh in file_handles:
            fh.close()

    progress.end()

    run.info('Total reads generated', pp(total_reads))
    for path in output_files:
        run.info('Output', path)


def get_args():
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('INPUT', "The reference FASTA file from which reads will be generated.")
    groupA.add_argument('-f', '--fasta-file', required=True, metavar='FASTA', help="Input reference FASTA file (can contain multiple contigs).")
    groupA.add_argument('-o', '--output-file-prefix', required=True, metavar='PREFIX', help="Output file prefix. Paired-end mode produces PREFIX-R1.fastq and PREFIX-R2.fastq, other modes produce PREFIX.fastq.")

    groupB = parser.add_argument_group('READ TYPE', "Choose a preset or specify the read type explicitly. "
                        "Presets set sensible defaults for all read parameters, which you can override individually.")
    groupB.add_argument('--preset', type=str, default=None, choices=list(PRESETS.keys()), help="Use a preset configuration. Available: " + ', '.join(PRESETS.keys()) + ".")
    groupB.add_argument('--read-type', type=str, default=None, choices=['paired-end', 'single-end', 'long-fixed', 'long-distributed'], help="Read type to generate. Inferred from preset if not given.")

    groupC = parser.add_argument_group('READ PARAMETERS', "Fine-tune read generation. These override preset values when both are given.")
    groupC.add_argument('--read-length', type=int, default=None, metavar='INT', help="Read length in bp. For long-distributed mode, this is the mean length.")
    groupC.add_argument('--read-length-std', type=float, default=None, metavar='FLOAT', help="Standard deviation of read length (only for long-distributed / ONT mode).")
    groupC.add_argument('--min-read-length', type=int, default=None, metavar='INT', help="Minimum read length for long-distributed mode. Default: 200.")
    groupC.add_argument('--coverage', type=float, default=None, metavar='FLOAT', help="Target coverage depth. Default: 50.")
    groupC.add_argument('--insert-size', type=int, default=None, metavar='INT', help="Mean insert size between R1 and R2 (paired-end only).")
    groupC.add_argument('--insert-size-std', type=int, default=None, metavar='INT', help="Standard deviation of insert size (paired-end only).")
    groupC.add_argument('--error-rate', type=float, default=None, metavar='FLOAT', help="Per-base substitution error probability. Default: 0 (no errors).")
    groupC.add_argument('--quality-score', type=str, default=None, metavar='CHAR', help="Fixed FASTQ quality score character. Default depends on preset (e.g., '?' for Q30 Illumina, 'F' for Q37 HiFi).")
    groupC.add_argument('--length-distribution', type=str, default=None, choices=['normal', 'lognormal'], help="Length distribution for long-distributed mode. 'lognormal' produces a right-skewed distribution realistic for ONT data. Default: normal (or set by preset).")

    groupD = parser.add_argument_group('SNV INJECTION', "Introduce SNVs at specific positions with controlled multi-allele frequencies. "
                        "Use --mutations-file for precise control, or --snv-density for random placement.")
    groupD.add_argument('--mutations-file', type=str, default=None, metavar='FILE', help="Path to a TAB-delimited file with columns: contig_name, position (0-indexed), freq_A, freq_T, freq_C, freq_G. Frequencies must sum to 1.0.")
    groupD.add_argument('--snv-density', type=float, default=None, metavar='FLOAT', help="Randomly place SNVs at this density (e.g., 0.01 = 1 SNV per 100 bp).")
    groupD.add_argument('--num-alleles', type=int, default=None, metavar='INT', help="Number of alleles (2-4) at each random SNV position. Only used with --snv-density. Default: 2.")

    groupE = parser.add_argument_group('REPRODUCIBILITY')
    groupE.add_argument('--seed', type=int, default=42, metavar='INT', help="Random seed for reproducibility. Default: 42.")

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
