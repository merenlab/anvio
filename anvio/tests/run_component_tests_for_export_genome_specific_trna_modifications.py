"""
Runs anvi-export-genome-specific-trna-modifications against pre-built mock data and
checks key properties of the output.

Mock data lives in:
  anvio/tests/sandbox/mock-data-export-genome-specific-tRNA-modification-profiles/

Run from any directory:
    python create_mock_data_and_test.py

To regenerate the mock data from scratch, see the git history for the original
create_mock_data_and_test.py before it was split.
"""

import os
import sys
import subprocess
import pandas as pd

# When this script lives in anvio/tests/, sandbox/ is a sibling directory.
# Pass an explicit path as argv[1] to override (useful before the script is moved).
_default_sandbox = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    'sandbox', 'mock-data-export-genome-specific-tRNA-modification-profiles'
)
SANDBOX = os.path.realpath(sys.argv[1] if len(sys.argv) > 1 else _default_sandbox)

TRNASEQ_DB_PATH = os.path.join(SANDBOX, 'TRNASEQ_CONTIGS.db')
GENOME_A_PATH   = os.path.join(SANDBOX, 'GENOME_A.db')
GENOME_B_PATH   = os.path.join(SANDBOX, 'GENOME_B.db')
EXT_GENOMES_PATH = os.path.join(SANDBOX, 'external_genomes.txt')
ENZYME_LIST_PATH = os.path.join(SANDBOX, 'enzyme_list.txt')
MODS_PATH       = os.path.join(SANDBOX, 'MODIFICATIONS.txt')
SEEDS_PATH      = os.path.join(SANDBOX, 'SEEDS_SPECIFIC.txt')
OUTPUT_PATH     = os.path.join(SANDBOX, 'output.txt')
ENZDIST_PATH    = os.path.join(SANDBOX, 'enzyme_distribution.txt')

# Scenario constants (must match what was used to create the mock data)
CONTIG_1 = 'c_001'
CONTIG_2 = 'c_002'

PASS = '\033[92mPASS\033[0m'
FAIL = '\033[91mFAIL\033[0m'

results = []

def check(label, condition, detail=''):
    status = PASS if condition else FAIL
    msg = f'  [{status}] {label}'
    if not condition and detail:
        msg += f'\n         {detail}'
    print(msg)
    results.append(condition)


def fix_external_genomes():
    """Rewrite external_genomes.txt with paths correct for this machine."""
    with open(EXT_GENOMES_PATH, 'w') as f:
        f.write('name\tcontigs_db_path\n')
        f.write(f'genome_A\t{GENOME_A_PATH}\n')
        f.write(f'genome_B\t{GENOME_B_PATH}\n')


def run_program():
    cmd = [
        'anvi-export-genome-specific-trna-modifications',
        '-t', TRNASEQ_DB_PATH,
        '-m', MODS_PATH,
        '-s', SEEDS_PATH,
        '-e', EXT_GENOMES_PATH,
        '--modification-enzyme-list', ENZYME_LIST_PATH,
        '-o', OUTPUT_PATH,
        '--enzyme-distribution-output', ENZDIST_PATH,
        '--min-coverage-for-detection', '20',
        '-W',
    ]
    print(f'\n  Running: {" ".join(cmd)}\n')
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print('  STDERR:', result.stderr[-2000:])
    return result.returncode == 0


def check_output():
    df = pd.read_csv(OUTPUT_PATH, sep='\t')
    print(f'\n  Output has {len(df)} rows and {len(df.columns)} columns')
    print(f'  detection_status counts:\n{df["detection_status"].value_counts().to_string()}')
    print()

    check('Output is non-empty', len(df) > 0)
    check('All expected columns present',
          set(df.columns) >= {'modifying_enzyme_name', 'detection_status', 'coverage_at_position',
                              'modified.fraction.reference_genome', 'A', 'C', 'G', 'T'})

    # detection status for cmoA / seed_1 / pos 34
    s1_34 = df[(df['contig_name']==CONTIG_1) & (df['canonical_position']==34) &
               (df['sample_name']=='S1') & (df['modifying_enzyme_name']=='cmoA')]
    check('cmoA seed_1 S1 pos34 — detection_status = detected',
          len(s1_34)==1 and s1_34.iloc[0]['detection_status']=='detected')

    s2_34 = df[(df['contig_name']==CONTIG_1) & (df['canonical_position']==34) &
               (df['sample_name']=='S2') & (df['modifying_enzyme_name']=='cmoA')]
    check('cmoA seed_1 S2 pos34 — detection_status = not_detected',
          len(s2_34)==1 and s2_34.iloc[0]['detection_status']=='not_detected')

    s3_34 = df[(df['contig_name']==CONTIG_1) & (df['canonical_position']==34) &
               (df['sample_name']=='S3') & (df['modifying_enzyme_name']=='cmoA')]
    check('cmoA seed_1 S3 pos34 — detection_status = insufficient_coverage',
          len(s3_34)==1 and s3_34.iloc[0]['detection_status']=='insufficient_coverage')

    # not_detected row has full coverage assigned to ref_genome nucleotide (T)
    if len(s2_34) == 1:
        r = s2_34.iloc[0]
        check('not_detected row: T = SEEDS_SPECIFIC coverage (30)', r['T'] == 30, f'got T={r["T"]}')
        check('not_detected row: A = 0', r['A'] == 0, f'got A={r["A"]}')
        check('not_detected row: C = 0', r['C'] == 0, f'got C={r["C"]}')
        check('not_detected row: G = 0', r['G'] == 0, f'got G={r["G"]}')
        check('not_detected row: coverage_at_position = 30',
              r['coverage_at_position'] == 30, f'got {r["coverage_at_position"]}')
        check('not_detected row: modified.fraction.reference_genome = 0.0',
              r['modified.fraction.reference_genome'] == 0.0,
              f'got {r["modified.fraction.reference_genome"]}')

    # insufficient_coverage row still reports position coverage
    if len(s3_34) == 1:
        r = s3_34.iloc[0]
        check('insufficient_coverage row: coverage_at_position still reported (5)',
              r['coverage_at_position'] == 5, f'got {r["coverage_at_position"]}')
        check('insufficient_coverage row: fractions are NA',
              pd.isna(r['modified.fraction.reference_genome']) or r['modified.fraction.reference_genome'] == 'NA',
              f'got {r["modified.fraction.reference_genome"]}')

    # detected row has correct fractions (A=0, C=80, G=0, T=20; ref=T)
    if len(s1_34) == 1:
        r = s1_34.iloc[0]
        check('detected row: modified.fraction.reference_genome = 0.80',
              abs(float(r['modified.fraction.reference_genome']) - 0.80) < 1e-6,
              f'got {r["modified.fraction.reference_genome"]}')
        check('detected row: modified.fraction.reference_abundant_nucleotide = 0.20',
              abs(float(r['modified.fraction.reference_abundant_nucleotide']) - 0.20) < 1e-6,
              f'got {r["modified.fraction.reference_abundant_nucleotide"]}')
        check('detected row: reference_abundant_nucleotide = C',
              r['reference_abundant_nucleotide'] == 'C', f'got {r["reference_abundant_nucleotide"]}')

    # isoacceptor specificity: cmoA (Specific: Leu-TAG) must not appear for Phe-GAA seed
    cmoa_seed2 = df[(df['contig_name']==CONTIG_2) & (df['modifying_enzyme_name']=='cmoA')]
    check('cmoA (Specific: Leu-TAG) does NOT appear for Phe-GAA seed', len(cmoa_seed2) == 0)

    # miaA (Non-specific) appears for both seeds
    miaa_seed1 = df[(df['contig_name']==CONTIG_1) & (df['modifying_enzyme_name']=='miaA')]
    miaa_seed2 = df[(df['contig_name']==CONTIG_2) & (df['modifying_enzyme_name']=='miaA')]
    check('miaA (Non-specific) appears for Leu-TAG seed', len(miaa_seed1) > 0)
    check('miaA (Non-specific) appears for Phe-GAA seed', len(miaa_seed2) > 0)

    # genome_B has no cmoA gene — no cmoA rows from it
    cmoa_gB = df[(df['genome_name']=='genome_B') & (df['modifying_enzyme_name']=='cmoA')]
    check('genome_B (no cmoA gene) has no cmoA rows', len(cmoa_gB) == 0)

    # seed_2 appears for both genome_A and genome_B via miaA
    miaa_gA_s2 = df[(df['genome_name']=='genome_A') & (df['contig_name']==CONTIG_2) & (df['modifying_enzyme_name']=='miaA')]
    miaa_gB_s2 = df[(df['genome_name']=='genome_B') & (df['contig_name']==CONTIG_2) & (df['modifying_enzyme_name']=='miaA')]
    check('miaA seed_2 appears in genome_A output', len(miaa_gA_s2) > 0)
    check('miaA seed_2 appears in genome_B output', len(miaa_gB_s2) > 0)

    # all detected fractions in [0, 1]
    det = df[df['detection_status']=='detected']
    fracs = pd.to_numeric(det['modified.fraction.reference_genome'], errors='coerce')
    check('All detected fractions are in [0, 1]',
          fracs.notna().all() and (fracs >= 0).all() and (fracs <= 1).all(),
          f'out-of-range: {fracs[(fracs<0)|(fracs>1)].tolist()}')

    # coverage_at_position == A+C+G+T for detected rows
    det2 = det.copy()
    for nt in ['A', 'C', 'G', 'T']:
        det2[nt] = pd.to_numeric(det2[nt], errors='coerce')
    det2['total'] = det2['A'] + det2['C'] + det2['G'] + det2['T']
    det2['cov'] = pd.to_numeric(det2['coverage_at_position'], errors='coerce')
    mismatch = det2[det2['total'] != det2['cov']]
    check('coverage_at_position == A+C+G+T for all detected rows',
          len(mismatch) == 0, f'{len(mismatch)} mismatches')

    # enzyme distribution table
    edist = pd.read_csv(ENZDIST_PATH, sep='\t')
    cmoa_gA_dist = edist[(edist['modifying_enzyme_name']=='cmoA') & (edist['genome_name']=='genome_A')]
    cmoa_gB_dist = edist[(edist['modifying_enzyme_name']=='cmoA') & (edist['genome_name']=='genome_B')]
    check('Enzyme distribution: cmoA gene_count=1 in genome_A',
          len(cmoa_gA_dist)==1 and cmoa_gA_dist.iloc[0]['gene_count']==1,
          f'got {cmoa_gA_dist["gene_count"].tolist()}')
    check('Enzyme distribution: cmoA gene_count=0 in genome_B',
          len(cmoa_gB_dist)==1 and cmoa_gB_dist.iloc[0]['gene_count']==0,
          f'got {cmoa_gB_dist["gene_count"].tolist()}')


if __name__ == '__main__':
    if not os.path.isdir(SANDBOX):
        print(f'Mock data directory not found: {SANDBOX}')
        print('Copy the mock data to anvio/tests/sandbox/mock-data-export-genome-specific-tRNA-modification-profiles/')
        sys.exit(1)

    fix_external_genomes()

    print('\n=== Running program ===')
    ok = run_program()
    if not ok:
        print(f'  [{FAIL}] Program exited with non-zero status — cannot check output')
        sys.exit(1)

    print('=== Checking output ===')
    check_output()

    passed = sum(results)
    total  = len(results)
    print(f'\n=== {passed}/{total} checks passed ===\n')
    sys.exit(0 if passed == total else 1)
