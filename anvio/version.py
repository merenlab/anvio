"""
Version information for anvi'o and its database artifacts.

This module centralizes all version information to avoid circular imports
during package installation and to provide a single source of truth for
version-related data.
"""

# Core anvi'o version information
anvio_version = '8-dev'
anvio_codename = 'marie'  # after Marie Tharp -- see the release notes for details: https:ithub.com/merenlab/anvio/releases/tag/v8

# Mandatory Python version
major_python_version_required = 3
minor_python_version_required = 10

# Version numbers for anvi'o artifacts -- if you need to change any
# of these numbers, you better have a well-tested migration script
# ready to go along with that change!
contigs_db_version = "24"
profile_db_version = "40"
genes_db_version = "6"
pan_db_version = "21"
pangraph_db_version = "2"
auxiliary_data_version = "2"
structure_db_version = "2"
genomes_storage_version = "7"
trnaseq_db_version = "2"
workflow_config_version = "3"
metabolic_modules_db_version = "4"

versions_for_db_types = {'contigs': contigs_db_version,
                         'profile': profile_db_version,
                         'genes': genes_db_version,
                         'structure': structure_db_version,
                         'pan': pan_db_version,
                         'pan-graph': pangraph_db_version,
                         'genomestorage': genomes_storage_version,
                         'auxiliary data for coverages': auxiliary_data_version,
                         'trnaseq': trnaseq_db_version,
                         'config': workflow_config_version,
                         'modules': metabolic_modules_db_version}

def get_versions():
    """
    Return version tuple for all anvi'o components.

    This function mirrors the original set_version() function from __init__.py
    but ensures database versions are loaded first.
    """

    # If you change anything here, it is extremely important to also
    # change anvio/__init__.py to make sure versions are read and
    # assigned to key variables in the same order to those that are here.
    return (anvio_version,
            anvio_codename,
            contigs_db_version,
            profile_db_version,
            genes_db_version,
            pan_db_version,
            pangraph_db_version,
            auxiliary_data_version,
            structure_db_version,
            genomes_storage_version,
            trnaseq_db_version,
            workflow_config_version,
            metabolic_modules_db_version)
