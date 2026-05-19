#!/usr/bin/env python

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
__requires__ = ['contigs-db', 'functions', 'genes-of-interest-txt']
__provides__ = ['functions-txt']
__description__ = "Export functions of genes from an anvi'o contigs database for a given annotation source"


def export_functions_in_matrix_format(functions_dict, sources_to_include, output_file, run, progress):
    """Export functions in wide format: one row per gene, sources as columns."""

    progress.new('Exporting functions in matrix format')
    progress.update('Building wide-format dictionary ...')

    matrix_dict = {}
    for entry in functions_dict.values():
        gene_id = entry['gene_callers_id']
        source = entry['source']

        if gene_id not in matrix_dict:
            matrix_dict[gene_id] = {}

        matrix_dict[gene_id][source] = entry['function'] if entry['function'] else ''
        matrix_dict[gene_id][f"{source} (ACCESSION)"] = entry['accession'] if entry['accession'] else ''

    headers = ['gene_callers_id']
    for source in sources_to_include:
        headers.append(f"{source} (ACCESSION)")
        headers.append(source)

    for gene_id in matrix_dict:
        for source in sources_to_include:
            if source not in matrix_dict[gene_id]:
                matrix_dict[gene_id][source] = ''
            if f"{source} (ACCESSION)" not in matrix_dict[gene_id]:
                matrix_dict[gene_id][f"{source} (ACCESSION)"] = ''

    progress.update('Writing output file ...')
    utils.store_dict_as_TAB_delimited_file(matrix_dict, output_file, headers=headers, key_header='gene_callers_id')
    progress.end()

    run.info('Number of genes reported', len(matrix_dict), mc='green')
    run.info('Output file', output_file)


CONTIG_INFO_FIELDS = ['contig_name', 'start', 'stop', 'direction', 'partial', 'contig_length', 'contig_gc_content']


def export_functions_in_long_format(functions_dict, output_file, run, progress, include_contig_info=False):
    """Export functions in long format: one row per annotation."""

    progress.new('Exporting functions')
    progress.update('...')

    headers = t.gene_function_calls_table_structure + (CONTIG_INFO_FIELDS if include_contig_info else [])

    output = open(output_file, 'w')
    output.write('\t'.join(headers) + '\n')

    sorted_entries = sorted(functions_dict.values(), key=lambda x: x['gene_callers_id'])

    num_entries_reported = 0
    for entry in sorted_entries:
        if entry['e_value'] is None:
            entry['e_value'] = 0.0

        output.write('\t'.join([str(entry[key]) for key in headers]) + '\n')
        num_entries_reported += 1

    output.close()
    progress.end()

    run.info('Number of functions reported', num_entries_reported, mc='green')
    run.info('Output file', output_file)


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
        if args.include_contig_info:
            genes_in_contigs_dict = contigs_db.db.get_table_as_dict(t.genes_in_contigs_table_name)
            contigs_info_dict = contigs_db.db.get_table_as_dict(t.contigs_info_table_name)
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

        if args.include_contig_info and args.matrix_format:
            raise ConfigError("The flag --include-contig-info has no effect with --matrix-format. Please "
                              "remove one of these flags and try again.")

        if args.include_contig_info and not args.matrix_format:
            progress.new('Enriching output with contig info')
            progress.update('...')
            for entry in functions_dict.values():
                gene_call = genes_in_contigs_dict.get(entry['gene_callers_id'], {})
                contig_name = gene_call.get('contig', '')
                contig_info = contigs_info_dict.get(contig_name, {})
                entry['contig_name'] = contig_name
                entry['start'] = gene_call.get('start', '')
                entry['stop'] = gene_call.get('stop', '')
                entry['direction'] = gene_call.get('direction', '')
                entry['partial'] = gene_call.get('partial', '')
                entry['contig_length'] = contig_info.get('length', '')
                entry['contig_gc_content'] = contig_info.get('gc_content', '')
            progress.end()

        if args.matrix_format:
            if args.annotation_sources:
                sources_to_include = [s.strip() for s in args.annotation_sources.split(',')]
            else:
                sources_to_include = sorted(annotation_sources)

            export_functions_in_matrix_format(functions_dict, sources_to_include, args.output_file, run, progress)
        else:
            export_functions_in_long_format(functions_dict, args.output_file, run, progress, include_contig_info=args.include_contig_info)
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
    parser.add_argument(*anvio.A('matrix-format'), **anvio.K('matrix-format', {'help': "Export functions in wide format: "
                    "one row per gene, sources as columns."}))
    parser.add_argument('--include-contig-info', action='store_true', default=False,
                        help="Extend the long-format output with positional and contig-level columns: "
                             "contig name, gene start/stop coordinates, strand direction, partial call "
                             "flag, contig length, and contig GC content. Has no effect with --matrix-format.")

    return parser.get_args(parser)

if __name__ == '__main__':
    main()
