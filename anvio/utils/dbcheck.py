from anvio.errors import ConfigError
from anvio.dbinfo import DBInfo as dbi


def is_contigs_db(db_path, dont_raise=False):
    dbi(db_path, expecting='contigs', dont_raise=dont_raise)
    return True



def is_trnaseq_db(db_path):
    dbi(db_path, expecting='trnaseq')
    return True



def is_pan_or_profile_db(db_path, genes_db_is_also_accepted=False):
    ok_db_types = ['pan', 'profile'] + (['genes'] if genes_db_is_also_accepted else [])
    dbi(db_path, expecting=ok_db_types)
    return True



def is_profile_db(db_path):
    dbi(db_path, expecting='profile')
    return True



def is_structure_db(db_path):
    dbi(db_path, expecting='structure')
    return True



def is_blank_profile(db_path):
    database = dbi(db_path, dont_raise=True)

    if database.db_type != 'profile':
        return False

    return database.blank



def is_pan_db(db_path):
    dbi(db_path, expecting='pan')
    return True



def is_genome_storage(db_path):
    dbi(db_path, expecting='genomestorage')
    return True



def is_genes_db(db_path):
    dbi(db_path, expecting='genes')
    return True



def is_kegg_modules_db(db_path):
    dbi(db_path, expecting='modules')
    return True



def is_profile_db_merged(profile_db_path):
    return dbi(profile_db_path, expecting='profile').merged



def is_profile_db_and_contigs_db_compatible(profile_db_path, contigs_db_path):
    # let's make sure we did get some paths
    if not profile_db_path or not contigs_db_path:
        if not profile_db_path and not contigs_db_path:
            missing = 'any paths for profile-db or contigs-db'
        else:
            missing = 'a profile-db path' if not profile_db_path else 'a contigs-db path'

        raise ConfigError(f"A function in anvi'o was about to see if the profile-db and the contigs-db someone "
                          f"wanted to work with were compatible with one another, but it was called without "
                          f"{missing} :/")

    pdb = dbi(profile_db_path)
    cdb = dbi(contigs_db_path)

    if cdb.hash != pdb.hash:
        raise ConfigError(f"The contigs database and the profile database at '{profile_db_path}' "
                          f"does not seem to be compatible. More specifically, this contigs "
                          f"database is not the one that was used when %s generated this profile "
                          f"database (%s != %s)." % ('anvi-merge' if pdb.merged else 'anvi-profile', cdb.hash, pdb.hash))
    return True



def is_structure_db_and_contigs_db_compatible(structure_db_path, contigs_db_path):
    sdb = dbi(structure_db_path)
    cdb = dbi(contigs_db_path)

    if cdb.hash != sdb.hash:
        raise ConfigError('The contigs and structure databases do not seem compatible. '
                          'More specifically, the contigs database is not the one that '
                          'was used when the structure database was created (%s != %s).'\
                               % (cdb.hash, sdb.hash))

    return True



def is_pan_db_and_genomes_storage_db_compatible(pan_db_path, genomes_storage_path):
    pdb = dbi(pan_db_path)
    gdb = dbi(genomes_storage_path)

    if pdb.hash != gdb.hash:
        raise ConfigError(f"The pan database and the genomes storage database do not seem to "
                          f"be compatible. More specifically, the genomes storage database is "
                          f"not the one that was used when the pangenome was created. "
                          f"({pdb.hash} != {gdb.hash})")

    return True
