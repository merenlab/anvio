#!/usr/bin/env python
# -*- coding: utf-8

import sys
import h5py
import argparse
import anvio.terminal as terminal

from anvio.errors import ConfigError

current_version = '3'
next_version    = '4'

run = terminal.Run()
progress = terminal.Progress()

def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    fp = h5py.File(db_path, 'a')

    if fp.attrs['version'] != current_version:
      fp.close()
      raise ConfigError("Genome storage version is not %s." % current_version)

    progress.new('Upgrading genome storage')
    genomes = fp['/data/genomes/'].keys()
    for genome in genomes:
      gene_caller_ids = fp['/data/genomes/%s/' % genome].keys()
      for genome_caller_id in gene_caller_ids:
        fp.move('/data/genomes/%s/%s/sequence' % (genome, genome_caller_id), '/data/genomes/%s/%s/aa_sequence' % (genome, genome_caller_id))
        fp['/data/genomes/%s/%s/dna_sequence' % (genome, genome_caller_id)] = ''
        progress.update('Upgrading genome "%s" and gene caller id "%s"' % (genome, genome_caller_id))
    progress.end()

    fp.attrs['version'] = next_version
    fp.close()

    run.info_single('The update is done! Your old genome storage file is saved as a backup at "%s". Feel free to\
                     delete it if everything seems to be working (if you are seeing this message and yet you still\
                     have not your command line back, you can kill this process by pressing CTRL + C and it will not\
                     affect anything --for some reason in some cases the process just hangs, and we have not been\
                     able to identify the problem).' % (backup_db_path), nl_after=1, nl_before=1)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A simple script to upgrade genomes storage from version %s to version %s' % (current_version, next_version))
    parser.add_argument('genomes_storage', metavar = 'GENOMES_STORAGE', help = "An anvi'o genomes storage of version %s" % current_version)
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.genomes_storage)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
