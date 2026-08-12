"""Asking NCBI about SRA runs, and making sense of what comes back.

Anvi'o workflows that download sequencing reads from the SRA need to know a few things about
each run *before* they can do anything with it: whether it is paired-end or single-end, whether
it came off a short-read or a long-read instrument, and roughly how much disk space it is going
to take up. None of that can be learned from an accession alone, and all of it has to be known
up front, since it determines which files a workflow will produce and how many runs it can
afford to keep on disk at once.

This module answers those questions by asking NCBI's Entrez API for the 'runinfo' table, turns
the answers into anvi'o's own vocabulary, and keeps them in a plain TAB-delimited cache file so
the question is asked exactly once. The cache is meant to be read by humans and edited by them:
when NCBI's metadata is wrong or incomplete (which happens, since much of it is whatever the
original submitter typed in), correcting the file is the way to set anvi'o straight.
"""

import os
import csv
import time

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['FlorianTrigodet']


run = terminal.Run()
progress = terminal.Progress()


# The read types anvi'o distinguishes. Everything downstream branches on these.
PAIRED_END_SHORT_READS = 'PE-SR'
SINGLE_END_SHORT_READS = 'SE-SR'
LONG_READS = 'LR'

# NCBI's `Platform` values, sorted into the two worlds anvi'o cares about. Anything that is not
# in either set is refused rather than guessed at, because guessing wrong here means running a
# short-read assembler on long reads (or the other way around).
SHORT_READ_PLATFORMS = {'ILLUMINA', 'BGISEQ', 'DNBSEQ', 'ELEMENT', 'ULTIMA',
                        'ION_TORRENT', 'LS454', 'ABI_SOLID', 'CAPILLARY'}
LONG_READ_PLATFORMS = {'OXFORD_NANOPORE', 'PACBIO_SMRT'}

# Which anvi'o long-read technology token to use for a given instrument. Nanopore is a single
# token regardless of the instrument, but PacBio chemistry is not always recoverable: the RS II
# and the original Sequel are CLR machines, the Revio only produces HiFi, and the Sequel II does
# both, so anvi'o refuses to pick for it and says so instead.
PACBIO_MODELS_TO_TECHNOLOGY = {'PACBIO RS': 'pb-clr',
                               'SEQUEL II': None,
                               'SEQUEL IIE': None,
                               'SEQUEL': 'pb-clr',
                               'REVIO': 'pb-hifi'}

# The columns of the cache file, in the order they are written.
CACHE_COLUMNS = ['accession', 'read_type', 'lr_technology', 'platform', 'model',
                 'library_layout', 'spots', 'bases', 'size_mb', 'source']

# Where the values in a cache row came from: straight from NCBI, or edited by a human. Anvi'o
# never overwrites a row a human has touched, and never nags about one either.
SOURCE_NCBI = 'ncbi'
SOURCE_USER = 'user'

# The columns of NCBI's runinfo table, in the order NCBI emits them. The table comes back as CSV
# with no header when it is fetched through a history server query, so the order is what
# identifies the fields.
RUNINFO_COLUMNS = ['Run', 'ReleaseDate', 'LoadDate', 'spots', 'bases', 'spots_with_mates',
                   'avgLength', 'size_MB', 'AssemblyName', 'download_path', 'Experiment',
                   'LibraryName', 'LibraryStrategy', 'LibrarySelection', 'LibrarySource',
                   'LibraryLayout', 'InsertSize', 'InsertDev', 'Platform', 'Model', 'SRAStudy',
                   'BioProject', 'Study_Pubmed_id', 'ProjectID', 'Sample', 'BioSample',
                   'SampleType', 'TaxID', 'ScientificName']

EUTILS_URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'

# How many accessions to ask about at once, and how long to wait between requests. NCBI asks for
# no more than three requests per second from clients without an API key.
ACCESSIONS_PER_REQUEST = 100
SECONDS_BETWEEN_REQUESTS = 0.4

# A FASTQ file holds each base twice (once as a nucleotide, once as a quality score) plus the
# read names, so it lands at roughly this many bytes per base once uncompressed.
BYTES_PER_BASE_IN_FASTQ = 2.2

# When NCBI reports no base count for a run we fall back on the size of the archive, which is
# compressed, and assume it expands by about this much.
ARCHIVE_TO_FASTQ_EXPANSION = 4.0


def is_valid_accession(accession):
    """Return True if this looks like an SRA run accession."""

    return accession.startswith(('SRR', 'ERR', 'DRR'))


def sanity_check_accessions(accessions, source_of_accessions):
    """Complain about anything that is not an SRA *run* accession.

    The most common mistake by far is reaching for a BioSample or a BioProject accession, since
    those are what papers tend to print, so those get an error of their own."""

    biosamples = [a for a in accessions if a.startswith(('SAMEA', 'SAMN', 'SAMD'))]
    bioprojects = [a for a in accessions if a.startswith('PRJ')]
    others = [a for a in accessions if not is_valid_accession(a) and a not in biosamples and a not in bioprojects]

    if biosamples or bioprojects:
        wrong_ones = biosamples + bioprojects
        raise ConfigError(f"Anvi'o found {terminal.pluralize('accession', len(wrong_ones))} in {source_of_accessions} "
                          f"that identify samples or projects rather than sequencing runs: "
                          f"{', '.join(sorted(wrong_ones)[:5])}. Only run accessions (the ones starting with SRR, "
                          f"ERR, or DRR) point at actual sequencing data. If you search for one of these on the NCBI "
                          f"SRA website (https://www.ncbi.nlm.nih.gov/sra) you will find the run accessions that "
                          f"belong to it, and those are what anvi'o needs here.")

    if others:
        raise ConfigError(f"Anvi'o does not recognize {terminal.pluralize('accession', len(others))} in "
                          f"{source_of_accessions} as an SRA run accession: {', '.join(sorted(others)[:5])}. "
                          f"SRA run accessions start with SRR, ERR, or DRR.")


def infer_read_type(platform, library_layout, model):
    """Work out what kind of reads an SRA run holds, from NCBI's own description of it.

    Returns a (read_type, lr_technology) tuple. `lr_technology` is one of anvi'o's long-read
    technology tokens when it can be determined with confidence, and None otherwise (either
    because these are short reads, or because the instrument model does not say which PacBio
    chemistry was used)."""

    platform = (platform or '').strip().upper()
    library_layout = (library_layout or '').strip().upper()
    model = (model or '').strip().upper()

    if platform == 'OXFORD_NANOPORE':
        return LONG_READS, 'ont'

    if platform == 'PACBIO_SMRT':
        # Longest match first, so that 'Sequel II' is not mistaken for 'Sequel'.
        for known_model in sorted(PACBIO_MODELS_TO_TECHNOLOGY, key=len, reverse=True):
            if known_model in model:
                return LONG_READS, PACBIO_MODELS_TO_TECHNOLOGY[known_model]
        return LONG_READS, None

    if platform in SHORT_READ_PLATFORMS:
        return (PAIRED_END_SHORT_READS if library_layout == 'PAIRED' else SINGLE_END_SHORT_READS), None

    raise ConfigError(f"Anvi'o has no idea what to make of the sequencing platform '{platform}', which is what NCBI "
                      f"reports for one of your accessions. It could be a platform anvi'o has not learned about yet. "
                      f"If you know what kind of reads these are, you can say so yourself by editing the "
                      f"`read_type` column of your SRA metadata file (use '{PAIRED_END_SHORT_READS}' for paired-end "
                      f"short reads or '{LONG_READS}' for long reads), and anvi'o will take your word for it.")


def predict_fastq_bytes(entry):
    """Estimate how big the uncompressed FASTQ file(s) of a run will be.

    NCBI reports the number of bases for most runs, which makes for a decent estimate. For the
    runs where it reports nothing (or zero, which happens more often than you would hope) we fall
    back on the size of the compressed archive and assume a generous expansion factor."""

    bases = entry.get('bases') or 0
    size_mb = entry.get('size_mb') or 0

    if bases > 0:
        return bases * BYTES_PER_BASE_IN_FASTQ

    return size_mb * 1024 * 1024 * ARCHIVE_TO_FASTQ_EXPANSION


def predict_peak_disk_usage_in_gb(entry, safety_factor=1.3):
    """Estimate the most disk space a single run will occupy at any one moment.

    A run does not simply appear as a FASTQ file. It arrives as an `.sra` archive, is unpacked
    into FASTQ (while `fasterq-dump` keeps scratch space of its own roughly the size of the
    output), gets compressed, and is then quality filtered into a second copy. The peak is what
    matters for a disk budget, so this counts the archive, the scratch space, the uncompressed
    FASTQ and one filtered copy of it."""

    fastq_bytes = predict_fastq_bytes(entry)
    archive_bytes = (entry.get('size_mb') or 0) * 1024 * 1024

    # archive + fasterq-dump scratch + the FASTQ itself + a quality-filtered copy
    peak_bytes = archive_bytes + (fastq_bytes * 3)

    return (peak_bytes / (1024 ** 3)) * safety_factor


def fetch_runinfo(accessions, run=run, progress=progress):
    """Ask NCBI for what it knows about a list of SRA run accessions.

    Returns {accession: {...}} using anvi'o's own vocabulary. Accessions NCBI has nothing to say
    about are simply absent from the result: it is the caller's job to notice and complain, since
    what to do about a missing accession depends on why it was asked for."""

    accessions = list(dict.fromkeys(accessions))
    chunks = [accessions[i:i + ACCESSIONS_PER_REQUEST] for i in range(0, len(accessions), ACCESSIONS_PER_REQUEST)]

    entries = {}

    progress.new('Asking NCBI about SRA runs', progress_total_items=len(chunks))
    for i, chunk in enumerate(chunks):
        progress.update(f"Batch {i + 1} of {len(chunks)} ({len(chunk)} accessions) ...")
        progress.increment(increment_to=i + 1)

        for accession, entry in _fetch_runinfo_for_one_batch(chunk).items():
            entries[accession] = entry

        if i < len(chunks) - 1:
            time.sleep(SECONDS_BETWEEN_REQUESTS)
    progress.end()

    run.info('SRA runs described by NCBI', f"{len(entries)} of {len(accessions)}")

    return entries


def _fetch_runinfo_for_one_batch(accessions):
    """Run a single esearch/efetch round trip for a batch of accessions."""

    search_url = f"{EUTILS_URL}esearch.fcgi?db=sra&usehistory=y&retmax={len(accessions)}&term={'+OR+'.join(accessions)}"

    try:
        search_response = utils.get_remote_file_content(search_url, timeout=60)
    except Exception as e:
        raise ConfigError(f"Anvi'o could not reach NCBI to look up the metadata of your SRA accessions. If you are "
                          f"on a computer without internet access, you can generate the metadata file elsewhere with "
                          f"`anvi-script-get-sra-metadata` and bring it along. This is what we know: {e}")

    query_key = _read_xml_value(search_response, 'QueryKey')
    web_env = _read_xml_value(search_response, 'WebEnv')

    if not query_key or not web_env:
        raise ConfigError("NCBI's answer to anvi'o's search for your SRA accessions did not include the bits anvi'o "
                          "needs to ask its follow-up question (namely a query key and a web environment). This "
                          "usually means NCBI is having a moment; trying again in a few minutes tends to work.")

    fetch_url = f"{EUTILS_URL}efetch.fcgi?db=sra&query_key={query_key}&WebEnv={web_env}&rettype=runinfo&retmode=csv"

    try:
        runinfo = utils.get_remote_file_content(fetch_url, timeout=120)
    except Exception as e:
        raise ConfigError(f"Anvi'o found your SRA accessions at NCBI, but then failed to download the table that "
                          f"describes them. This is what we know: {e}")

    return _parse_runinfo(runinfo)


def _read_xml_value(text, tag):
    """Pull a single value out of an XML document without pretending to be an XML parser."""

    opening, closing = f"<{tag}>", f"</{tag}>"

    if opening not in text or closing not in text:
        return None

    return text.split(opening, 1)[1].split(closing, 1)[0].strip()


def _parse_runinfo(runinfo_text):
    """Turn NCBI's runinfo CSV into anvi'o's vocabulary, keyed by accession."""

    entries = {}

    for fields in csv.reader(runinfo_text.splitlines()):
        if not fields or not fields[0].strip():
            continue

        # Depending on how the table was requested it may or may not carry a header line.
        if fields[0] == 'Run':
            continue

        if len(fields) < len(RUNINFO_COLUMNS):
            continue

        row = dict(zip(RUNINFO_COLUMNS, fields))
        accession = row['Run'].strip()

        if not is_valid_accession(accession):
            continue

        read_type, lr_technology = infer_read_type(row['Platform'], row['LibraryLayout'], row['Model'])

        entries[accession] = {'accession': accession,
                              'read_type': read_type,
                              'lr_technology': lr_technology,
                              'platform': row['Platform'].strip(),
                              'model': row['Model'].strip(),
                              'library_layout': row['LibraryLayout'].strip(),
                              'spots': _as_int(row['spots']),
                              'bases': _as_int(row['bases']),
                              'size_mb': _as_int(row['size_MB']),
                              'source': SOURCE_NCBI}

    return entries


def _as_int(value):
    """Read a number NCBI may well have left blank."""

    try:
        return int(float(str(value).strip()))
    except (TypeError, ValueError):
        return 0


def read_cache(cache_path):
    """Read a metadata cache file into {accession: {...}}."""

    if not os.path.exists(cache_path):
        return {}

    filesnpaths.is_file_tab_delimited(cache_path)

    columns_found = utils.get_columns_of_TAB_delim_file(cache_path, include_first_column=True)
    missing_columns = [c for c in CACHE_COLUMNS if c not in columns_found]

    if missing_columns:
        raise ConfigError(f"The SRA metadata file at '{cache_path}' is missing "
                          f"{terminal.pluralize('column', len(missing_columns))} anvi'o needs: "
                          f"{', '.join(missing_columns)}. If you have been editing this file by hand and something "
                          f"went sideways, the simplest fix is to delete it and let anvi'o build a new one.")

    entries = {}

    for accession, row in utils.get_TAB_delimited_file_as_dictionary(cache_path).items():
        entry = {'accession': accession,
                 'read_type': (row.get('read_type') or '').strip(),
                 'lr_technology': (row.get('lr_technology') or '').strip() or None,
                 'platform': (row.get('platform') or '').strip(),
                 'model': (row.get('model') or '').strip(),
                 'library_layout': (row.get('library_layout') or '').strip(),
                 'spots': _as_int(row.get('spots')),
                 'bases': _as_int(row.get('bases')),
                 'size_mb': _as_int(row.get('size_mb')),
                 'source': (row.get('source') or SOURCE_NCBI).strip()}

        if entry['read_type'] not in (PAIRED_END_SHORT_READS, SINGLE_END_SHORT_READS, LONG_READS):
            raise ConfigError(f"The accession {accession} in your SRA metadata file at '{cache_path}' has a "
                              f"`read_type` of '{entry['read_type']}', which anvi'o does not recognize. It must be "
                              f"one of '{PAIRED_END_SHORT_READS}' (paired-end short reads), "
                              f"'{SINGLE_END_SHORT_READS}' (single-end short reads), or '{LONG_READS}' (long reads).")

        entries[accession] = entry

    return entries


def write_cache(cache_path, entries):
    """Write {accession: {...}} out as a TAB-delimited file, sorted by accession."""

    filesnpaths.is_output_file_writable(cache_path)

    with open(cache_path, 'w') as output:
        output.write('\t'.join(CACHE_COLUMNS) + '\n')
        for accession in sorted(entries):
            entry = entries[accession]
            output.write('\t'.join([str(entry.get(c) if entry.get(c) is not None else '') for c in CACHE_COLUMNS]) + '\n')


def get_metadata_for_accessions(accessions, cache_path, source_of_accessions, fetch_missing=True,
                                run=run, progress=progress):
    """Return everything anvi'o knows about a set of accessions, asking NCBI only if it must.

    Whatever is already in the cache file is used as-is, so a row someone has corrected by hand
    stays corrected. Only the accessions that are missing from it are looked up, and the file is
    then rewritten with the new rows added to the old ones."""

    sanity_check_accessions(accessions, source_of_accessions)

    entries = read_cache(cache_path)
    missing = [a for a in accessions if a not in entries]

    if missing and not fetch_missing:
        raise ConfigError(f"Anvi'o needs to know a few things about {terminal.pluralize('SRA run', len(missing))} "
                          f"before it can do anything with them, and the metadata file at '{cache_path}' has nothing "
                          f"to say about {'them' if len(missing) > 1 else 'it'}: "
                          f"{', '.join(sorted(missing)[:5])}. You can fill in the gaps by running "
                          f"`anvi-script-get-sra-metadata` on a computer with internet access.")

    if missing:
        run.warning(f"Anvi'o is about to ask NCBI about {terminal.pluralize('SRA run', len(missing))} it has not "
                    f"seen before. The answers will be kept in '{cache_path}', so this only happens once.",
                    header="REACHING OUT TO NCBI", lc="green")

        entries.update(fetch_runinfo(missing, run=run, progress=progress))
        write_cache(cache_path, entries)

    still_missing = [a for a in accessions if a not in entries]
    if still_missing:
        raise ConfigError(f"NCBI had nothing to say about {terminal.pluralize('accession', len(still_missing))} from "
                          f"{source_of_accessions}: {', '.join(sorted(still_missing)[:5])}. Accessions that have been "
                          f"suppressed or withdrawn behave exactly like this, and so do typos. Please double check "
                          f"{'them' if len(still_missing) > 1 else 'it'} at "
                          f"https://www.ncbi.nlm.nih.gov/sra before trying again.")

    return {a: entries[a] for a in accessions}


def warn_about_long_read_technologies(entries, run=run, cache_path=None):
    """Point out long-read runs whose sequencing chemistry anvi'o could not pin down.

    Anvi'o does not stop for these. A workflow where every long-read sample comes off the same
    instrument does not need per-sample technologies at all, since the tool presets can simply be
    set in the config file. But a mixed set of samples does need them, and nobody wants to find
    that out after the fact, so we say something either way."""

    long_read_accessions = [a for a, e in entries.items() if e['read_type'] == LONG_READS]

    if not long_read_accessions:
        return

    undetermined = sorted(a for a in long_read_accessions
                          if not entries[a]['lr_technology'] and entries[a]['source'] != SOURCE_USER)

    if not undetermined:
        return

    models = ', '.join(sorted(set(entries[a]['model'] for a in undetermined if entries[a]['model'])))

    run.warning(f"Anvi'o could tell that {terminal.pluralize('accession', len(undetermined))} in your samples-txt "
                f"holds long reads, but not which sequencing chemistry produced them, because the instrument "
                f"{'models' if len(undetermined) > 1 else 'model'} NCBI reports ({models}) can be used for more than "
                f"one. Here {'they are' if len(undetermined) > 1 else 'it is'}: {', '.join(undetermined[:10])}"
                f"{' (and more)' if len(undetermined) > 10 else ''}.\n\n"
                f"If all of your long-read samples came off the same kind of machine, you can ignore this and set "
                f"the long-read presets in your config file as usual. If they did not, anvi'o needs to be told which "
                f"is which: either add an `lr_technology` column to your samples-txt, or fill in the `lr_technology` "
                f"column of "
                f"{'the metadata file at ' + repr(cache_path) if cache_path else 'your SRA metadata file'} "
                f"(anvi'o will not touch rows you have edited). Valid values are 'ont', 'pb-clr', and 'pb-hifi'.",
                header="A WORD ABOUT YOUR LONG READS", lc="yellow")
