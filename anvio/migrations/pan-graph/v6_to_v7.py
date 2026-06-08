#!/usr/bin/env python

import sys
import argparse

import anvio.dbinfo as dbinfo
import anvio.tables as t
import anvio.terminal as terminal

from anvio.errors import ConfigError

current_version, next_version = [x[1:] for x in __name__.split('_to_')]

run = terminal.Run()
progress = terminal.Progress()


def _column_names(pan_graph_db, table_name):
    """Return the column names of an existing SQLite table without reading rows."""
    return [d[0] for d in pan_graph_db._exec(
        f"SELECT * FROM {table_name} LIMIT 0").description]


def _migrate_regions(pan_graph_db):
    """v6 stored region_id as ``"<component_id>_<rid>"`` with no separate
    ``component_id`` column. v7 splits these into two columns and strips
    the prefix so per-component region_ids start from 0."""
    old_columns = _column_names(pan_graph_db, "pan_graph_regions")
    if 'component_id' in old_columns:
        # Already migrated (manual repair) -- skip silently.
        return

    region_id_idx = old_columns.index('region_id')
    rows = pan_graph_db._exec("SELECT * FROM pan_graph_regions").fetchall()

    new_rows = []
    for row in rows:
        prefixed = row[region_id_idx]
        if prefixed is None:
            cid, rid = 0, ''
        else:
            head, _, tail = str(prefixed).partition('_')
            if tail == '':
                # Shouldn't happen in a valid v6 db, but fall back gracefully.
                cid, rid = 0, head
            else:
                try:
                    cid = int(head)
                except ValueError:
                    cid = 0
                rid = tail
        rest = [row[i] for i, c in enumerate(old_columns) if c != 'region_id']
        # Map old column order -> new (without region_id), then prepend the
        # split (component_id, region_id) pair before reordering.
        carryover = {c: v for c, v in zip(
            [c for c in old_columns if c != 'region_id'], rest)}
        carryover['component_id'] = cid
        carryover['region_id'] = rid
        new_rows.append([carryover.get(col) for col in t.pan_graph_regions_table_structure])

    pan_graph_db._exec("DROP TABLE pan_graph_regions")
    pan_graph_db.create_table(
        t.pan_graph_regions_table_name,
        t.pan_graph_regions_table_structure,
        t.pan_graph_regions_table_types,
    )
    placeholders = ','.join(['?'] * len(t.pan_graph_regions_table_structure))
    pan_graph_db._exec_many(
        f"INSERT INTO {t.pan_graph_regions_table_name} VALUES ({placeholders})",
        new_rows,
    )


def _migrate_edges(pan_graph_db):
    """v7 retires the pre-``genomes_json`` ``directions`` column on
    pan_graph_edges: if a db still carries it, the column is dropped and
    ``genomes_json`` is filled with an empty list (data semantics differ;
    the original ``directions`` payload is not convertible). A per-edge
    ``component_id`` is intentionally NOT added -- the component can be
    joined from either endpoint in pan_graph_nodes."""
    old_columns = _column_names(pan_graph_db, "pan_graph_edges")
    if 'genomes_json' in old_columns and 'directions' not in old_columns:
        # Already on the target shape.
        return

    has_genomes_json = 'genomes_json' in old_columns
    genomes_idx = old_columns.index('genomes_json') if has_genomes_json else None

    rows = pan_graph_db._exec("SELECT * FROM pan_graph_edges").fetchall()
    new_rows = []
    for row in rows:
        genomes_json = row[genomes_idx] if has_genomes_json else '[]'

        carryover = {c: v for c, v in zip(old_columns, row)}
        carryover['genomes_json'] = genomes_json
        # Drop any columns no longer in v7 (e.g. legacy 'directions') by
        # only emitting the v7 structure.
        new_rows.append([carryover.get(col) for col in t.pan_graph_edges_table_structure])

    pan_graph_db._exec("DROP TABLE pan_graph_edges")
    pan_graph_db.create_table(
        t.pan_graph_edges_table_name,
        t.pan_graph_edges_table_structure,
        t.pan_graph_edges_table_types,
    )
    placeholders = ','.join(['?'] * len(t.pan_graph_edges_table_structure))
    pan_graph_db._exec_many(
        f"INSERT INTO {t.pan_graph_edges_table_name} VALUES ({placeholders})",
        new_rows,
    )


def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    pan_graph_db_info = dbinfo.PanGraphDBInfo(db_path)
    if str(pan_graph_db_info.version) != current_version:
        raise ConfigError(
            f"The version of the provided pan-graph database is {pan_graph_db_info.version}, "
            f"not the required version, {current_version}, so this script cannot upgrade the database.")

    pan_graph_db = pan_graph_db_info.load_db()

    progress.new("Migrating")

    progress.update("Rewriting pan_graph_regions (split component_id + plain region_id) ...")
    _migrate_regions(pan_graph_db)

    progress.update("Rewriting pan_graph_edges (drop legacy 'directions' if present) ...")
    _migrate_edges(pan_graph_db)

    progress.update("Updating version")
    pan_graph_db.remove_meta_key_value_pair('version')
    pan_graph_db.set_version(next_version)

    progress.update("Committing changes")
    pan_graph_db.disconnect()

    progress.end()

    run.info_single(
        f"Your pan-graph database is now version {next_version}. pan_graph_regions now carries a "
        f"separate `component_id` and plain per-component region_ids; on pan_graph_edges, the legacy "
        f"`directions` column has been replaced with `genomes_json` (if it was still present). The "
        f"per-edge component can be joined from pan_graph_nodes via either endpoint.",
        nl_before=1, nl_after=1, mc='green'
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=f"A simple script to upgrade an anvi'o pan-graph database from version {current_version} to version {next_version}"
    )
    parser.add_argument('pan_graph_db', metavar='PAN_GRAPH_DB',
                        help=f"An anvi'o pan-graph database of version {current_version}")
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.pan_graph_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
