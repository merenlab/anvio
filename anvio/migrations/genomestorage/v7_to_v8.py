#!/usr/bin/env python
# -*- coding: utf-8

"""Migration script that detects and repairs bytes-string corruption in genomes
storage databases that were converted from HDF5 to SQLite via the v4-to-v5
migration. In Python 3, h5py returns bytes for string data. The v4-to-v5 script
called str() on those bytes values, producing literal 'b"..."' or "b'...'" text
in the database instead of clean strings."""

import sys
import argparse

import anvio.db as db
import anvio.terminal as terminal

from anvio.errors import ConfigError

current_version, next_version = [x[1:] for x in __name__.split('_to_')]

run = terminal.Run()
progress = terminal.Progress()


def strip_bytes_wrapper(value):
    """If a string value looks like a Python bytes repr (b'...' or b"..."),
    return the inner string. Otherwise return None."""

    if not isinstance(value, str):
        return None

    if (value.startswith("b'") and value.endswith("'")) or \
       (value.startswith('b"') and value.endswith('"')):
        return value[2:-1]

    return None


def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    genomes_db = db.DB(db_path, None, ignore_version=True)

    if str(genomes_db.get_version()) != current_version:
        genomes_db.disconnect()
        raise ConfigError("Version of this genome storage is not %s (hence, this script cannot "
                          "really do anything)." % current_version)

    total_fixes = 0

    progress.new("Checking for bytes-string corruption")

    # -----------------------------------------------------------------------
    # A. Fix meta 'hash' value
    # -----------------------------------------------------------------------
    progress.update("Checking meta hash value")

    hash_value = genomes_db.get_meta_value('hash')
    clean_hash = strip_bytes_wrapper(str(hash_value))
    if clean_hash is not None:
        genomes_db.set_meta_value('hash', clean_hash)
        run.info_single("Fixed corrupted hash meta value: '%s' -> '%s'" % (hash_value, clean_hash),
                        nl_before=1)
        total_fixes += 1

    # -----------------------------------------------------------------------
    # B. Fix gene_function_calls accession + function (correlated)
    #
    # The v4-to-v5 bug: str(b'accession|||function').split('|||') produced
    # accession = "b'real_accession" and function = "real_function'"
    # (or b"/..." variant when the value contained single quotes).
    # -----------------------------------------------------------------------
    progress.update("Scanning gene_function_calls for corrupted accession/function values")

    # Count single-quote variant
    response = genomes_db._exec("""SELECT COUNT(*) FROM gene_function_calls WHERE accession LIKE 'b''%'""")
    count_sq = response.fetchone()[0]

    # Count double-quote variant
    response = genomes_db._exec('''SELECT COUNT(*) FROM gene_function_calls WHERE accession LIKE 'b"%' ''')
    count_dq = response.fetchone()[0]

    func_corrupted = count_sq + count_dq

    if func_corrupted > 0:
        progress.update("Fixing %d corrupted accession/function entries" % func_corrupted)

        # Fix single-quote variant: accession starts with b', function ends with '
        if count_sq > 0:
            genomes_db._exec("""UPDATE gene_function_calls
                                   SET accession = SUBSTR(accession, 3),
                                       function = SUBSTR(function, 1, LENGTH(function) - 1)
                                 WHERE accession LIKE 'b''%'""")

        # Fix double-quote variant: accession starts with b", function ends with "
        if count_dq > 0:
            genomes_db._exec('''UPDATE gene_function_calls
                                   SET accession = SUBSTR(accession, 3),
                                       function = SUBSTR(function, 1, LENGTH(function) - 1)
                                 WHERE accession LIKE 'b"%' ''')

        total_fixes += func_corrupted
        run.info_single("Fixed %d corrupted accession/function pairs in gene_function_calls "
                        "(%d single-quote, %d double-quote variant)." % (func_corrupted, count_sq, count_dq))

    # -----------------------------------------------------------------------
    # C. Fix gene_function_calls source column (defensive — only affected if
    #    migration ran with h5py < 3.0 where group keys were bytes)
    # -----------------------------------------------------------------------
    progress.update("Checking gene_function_calls source column")

    response = genomes_db._exec("""SELECT COUNT(*) FROM gene_function_calls
                                    WHERE (source LIKE 'b''%''' OR source LIKE 'b"%"')""")
    count_source = response.fetchone()[0]

    if count_source > 0:
        genomes_db._exec("""UPDATE gene_function_calls
                               SET source = SUBSTR(source, 3, LENGTH(source) - 3)
                             WHERE source LIKE 'b''%'''""")
        genomes_db._exec('''UPDATE gene_function_calls
                               SET source = SUBSTR(source, 3, LENGTH(source) - 3)
                             WHERE source LIKE 'b"%"' ''')
        total_fixes += count_source
        run.info_single("Fixed %d corrupted source values in gene_function_calls." % count_source)

    # -----------------------------------------------------------------------
    # D. Fix genome_info genome_hash column
    # -----------------------------------------------------------------------
    progress.update("Checking genome_info genome_hash column")

    response = genomes_db._exec("""SELECT COUNT(*) FROM genome_info
                                    WHERE (genome_hash LIKE 'b''%''' OR genome_hash LIKE 'b"%"')""")
    count_ghash = response.fetchone()[0]

    if count_ghash > 0:
        genomes_db._exec("""UPDATE genome_info
                               SET genome_hash = SUBSTR(genome_hash, 3, LENGTH(genome_hash) - 3)
                             WHERE genome_hash LIKE 'b''%'''""")
        genomes_db._exec('''UPDATE genome_info
                               SET genome_hash = SUBSTR(genome_hash, 3, LENGTH(genome_hash) - 3)
                             WHERE genome_hash LIKE 'b"%"' ''')
        total_fixes += count_ghash
        run.info_single("Fixed %d corrupted genome_hash values in genome_info." % count_ghash)

    # -----------------------------------------------------------------------
    # E. Fix genome_name across all three tables (defensive)
    # -----------------------------------------------------------------------
    progress.update("Checking genome_name columns")

    for table_name in ['genome_info', 'gene_info', 'gene_function_calls']:
        response = genomes_db._exec("""SELECT COUNT(*) FROM %s
                                        WHERE (genome_name LIKE 'b''%%''' OR genome_name LIKE 'b"%%"')""" % table_name)
        count_gn = response.fetchone()[0]

        if count_gn > 0:
            genomes_db._exec("""UPDATE %s
                                   SET genome_name = SUBSTR(genome_name, 3, LENGTH(genome_name) - 3)
                                 WHERE genome_name LIKE 'b''%%'''""" % table_name)
            genomes_db._exec('''UPDATE %s
                                   SET genome_name = SUBSTR(genome_name, 3, LENGTH(genome_name) - 3)
                                 WHERE genome_name LIKE 'b"%%"' ''' % table_name)
            total_fixes += count_gn
            run.info_single("Fixed %d corrupted genome_name values in %s." % (count_gn, table_name))

    # -----------------------------------------------------------------------
    # F. Fix BLOB sequences in gene_info (defensive — bytes inserted directly
    #    into sqlite3 are stored as BLOBs rather than TEXT)
    # -----------------------------------------------------------------------
    progress.update("Checking for BLOB-typed sequences in gene_info")

    response = genomes_db._exec("""SELECT COUNT(*) FROM gene_info WHERE typeof(aa_sequence) = 'blob'""")
    blob_aa = response.fetchone()[0]

    response = genomes_db._exec("""SELECT COUNT(*) FROM gene_info WHERE typeof(dna_sequence) = 'blob'""")
    blob_dna = response.fetchone()[0]

    if blob_aa > 0 or blob_dna > 0:
        if blob_aa > 0:
            genomes_db._exec("""UPDATE gene_info SET aa_sequence = CAST(aa_sequence AS TEXT)
                                 WHERE typeof(aa_sequence) = 'blob'""")
        if blob_dna > 0:
            genomes_db._exec("""UPDATE gene_info SET dna_sequence = CAST(dna_sequence AS TEXT)
                                 WHERE typeof(dna_sequence) = 'blob'""")
        total_fixes += blob_aa + blob_dna
        run.info_single("Converted %d BLOB aa_sequences and %d BLOB dna_sequences to TEXT in "
                        "gene_info." % (blob_aa, blob_dna))

    # -----------------------------------------------------------------------
    # G. Rebuild gene_function_sources meta value from now-clean data
    # -----------------------------------------------------------------------
    progress.update("Rebuilding gene_function_sources meta value")

    sources_response = genomes_db._exec("""SELECT DISTINCT source FROM gene_function_calls""")
    sources = sorted([row[0] for row in sources_response.fetchall()])

    if sources:
        genomes_db.set_meta_value('gene_function_sources', ','.join(sources))

    # -----------------------------------------------------------------------
    # H. Version bump
    # -----------------------------------------------------------------------
    genomes_db.remove_meta_key_value_pair('version')
    genomes_db.set_version(next_version)
    genomes_db.disconnect()

    progress.end()

    if total_fixes > 0:
        run.info_single("Your genomes storage is now version %s. This migration fixed bytes-string "
                        "corruption that was introduced during a previous HDF5-to-SQLite conversion "
                        "(the v4-to-v5 migration). See above for details on what was fixed." % next_version,
                        nl_after=1, nl_before=1, mc='green')
    else:
        run.info_single("Your genomes storage is now version %s. No bytes-string corruption was "
                        "detected -- your database was either not affected by the v4-to-v5 migration "
                        "bug, or was created natively as SQLite." % next_version,
                        nl_after=1, nl_before=1, mc='green')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A simple script to upgrade a genome storage '
                                     'from version %s to version %s' % (current_version, next_version))
    parser.add_argument('genomes_storage', metavar='GENOMES_STORAGE',
                        help="An anvi'o genomes storage of version %s" % current_version)
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.genomes_storage)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
