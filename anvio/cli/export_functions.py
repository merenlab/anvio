#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio
import anvio.tables as t
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ['contigs-db', 'functions']
__provides__ = ['functions-txt']
__description__ = "Export functions of genes from an anvi'o contigs database for a given annotation source"



def main():
    args = get_args()
    run = terminal.Run()
    progress = terminal.Progress()

    try:
        utils.is_contigs_db(args.contigs_db)

        progress.new('Initializing')
        progress.update('...')
        contigs_db = dbops.ContigsDatabase(args.contigs_db)
        annotation_sources = contigs_db.meta['gene_function_sources']
        functions_dict = contigs_db.db.get_table_as_dict(t.gene_function_calls_table_name)
        contigs_db.disconnect()
        progress.end()

        if not annotation_sources:
            raise ConfigError("This contigs database does not have any functional functions :/")

        if args.list_annotation_sources:
            run.warning('', 'FUNCTIONAL ANNOTATION SOURCE%s FOUND' % 'S' if len(annotation_sources) > 1 else '', lc='yellow')
            for annotation_source in annotation_sources:
                if annotation_source == annotation_sources[-1]:
                    run.info_single('%s' % annotation_source, nl_after=1)
                else:
                    run.info_single('%s' % annotation_source)
            sys.exit(0)

        if not args.output_file:
            raise ConfigError("You should provide an output file name.")

        filesnpaths.is_output_file_writable(args.output_file)

        if args.genes_of_interest:
            requested_genes = [x.strip() for x in open(args.genes_of_interest).readlines()]
            run.info('Gene calls requested', f"{', '.join(requested_genes)}")

            for g in requested_genes:
                if not g.isnumeric():
                    raise ConfigError(f"Anvi'o found a gene caller ID that was not a number. Specifically, this one: {g}. "
                                      f"This is a no-no, so please make sure all your gene caller IDs are integers.")

            progress.new("Filtering the functions dict for requested gene calls")
            functions_dict = utils.get_filtered_dict(functions_dict, 'gene_callers_id', set([int(x) for x in requested_genes]))
            progress.end()

        run.info('Annotation sources', '%s.' % (', '.join(annotation_sources)))

        if args.annotation_sources:
            requested_sources = [s.strip() for s in args.annotation_sources.split(',')]

            missing_sources = [s for s in requested_sources if s not in annotation_sources]
            if len(missing_sources):
                raise ConfigError("One or more of the annotation sources you requested does not appear to be in the "
                                   "contigs database :/ Here is the list: %s." % (', '.join(missing_sources)))

            run.info('Requested sources', '%s.' % (', '.join(requested_sources)))

            progress.new('Initializing')
            progress.update('Filtering the functions dict ...')
            functions_dict = utils.get_filtered_dict(functions_dict, 'source', set(requested_sources))
            progress.end()


        progress.new('Exporting functions')
        progress.update('...')
        output = open(args.output_file, 'w')
        output.write('\t'.join(t.gene_function_calls_table_structure) + '\n')
        num_entries_reported = 0
        for entry in list(functions_dict.values()):
            if entry['e_value'] == None:
                entry['e_value'] = 0.0

            output.write('\t'.join([str(entry[key]) for key in t.gene_function_calls_table_structure]) + '\n')
            num_entries_reported += 1
        output.close()
        progress.end()

        run.info('Number of functions reported', num_entries_reported, mc='green')
        run.info('Output file', args.output_file)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    parser.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))
    parser.add_argument(*anvio.A('output-file'), **anvio.K('output-file'))
    parser.add_argument(*anvio.A('annotation-sources'), **anvio.K('annotation-sources'))
    parser.add_argument(*anvio.A('list-annotation-sources'), **anvio.K('list-annotation-sources'))
    parser.add_argument(*anvio.A('genes-of-interest'), **anvio.K('genes-of-interest'))

    return parser.get_args(parser)

if __name__ == '__main__':
    main()
