#!/usr/bin/env python

import sys
import json
import argparse

import anvio.dbinfo as dbinfo
import anvio.terminal as terminal

from anvio.errors import ConfigError

current_version, next_version = [x[1:] for x in __name__.split('_to_')]

run = terminal.Run()
progress = terminal.Progress()

# The target (v7) pan-graph table schemas are declared inline on purpose. A
# migration must pin the schema exactly as it was at the time of the version
# bump and must NOT import anvio.tables -- that module always tracks the current
# head and would silently drift as the schema evolves in later versions,
# breaking this script retroactively.
#
# v7 turns the numeric, 0-based component_id and the plain numeric-string
# region_id into 1-based prefixed *names* that mirror the gene-cluster id
# convention: component_id becomes "CP_0001", "CP_0002", ... (CP_0001 is still
# the largest component) and region_id embeds its component plus a 1-based index
# ("CP_0001_1", "CP_0001_2", ...), so a region id is globally unique on its own.
# component_id is therefore a text column in v7.
pan_graph_nodes_table_name      = 'pan_graph_nodes'
pan_graph_nodes_table_structure = ['node_id', 'node_type', 'region_id', 'gene_cluster_id', 'gene_calls_json', 'synteny_position_json', 'node_x', 'node_y', 'alignment_summary', 'component_id']
pan_graph_nodes_table_types     = [  'str'  ,    'str'   ,    'str'   ,       'str'      ,       'str'      ,          'str'         , 'numeric', 'numeric',       'str'       ,    'str'      ]

pan_graph_regions_table_name      = 'pan_graph_regions'
pan_graph_regions_table_structure = ['component_id', 'region_id', 'region_type', 'x_min', 'x_max', 'num_synteny_gene_clusters', 'num_gene_clusters', 'num_gene_calls', 'max_expansion', 'min_expansion', 'complexity', 'complexity_normalized', 'diversity', 'diversity_normalized', 'weight', 'weight_normalized', 'composite_variability_score', 'complexity_mm_scaled', 'diversity_mm_scaled', 'expansion_mm_scaled', 'weight_mm_scaled']
pan_graph_regions_table_types     = [    'str'     ,    'str'   ,     'str'    , 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric']


def _to_component_name(value):
    """v6 stored component_id as a 0-based number (0 = largest). v7 stores a
    1-based, zero-padded prefixed name ("CP_0001", "CP_0002", ...).
    Already-migrated names pass through untouched so the migration is
    idempotent."""
    if value is None or value == '':
        return "CP_0001"
    if isinstance(value, str) and value.startswith('CP_'):
        return value
    return f"CP_{int(float(value)) + 1:04d}"


def _to_region_name(component_value, region_value):
    """v6 stored region_id as a plain 0-based numeric string ("0", "1", ...)
    unique only within its component. v7 embeds the (renamed) component and a
    1-based index ("CP_0001_1", "CP_0001_2", ...) so the id is globally unique.
    Empty/unset region ids (nodes outside any region) and already-migrated names
    (which contain the "CP_" prefix) pass through."""
    if region_value is None or region_value == '':
        return region_value
    if isinstance(region_value, str) and region_value.startswith('CP_'):
        return region_value
    return f"{_to_component_name(component_value)}_{int(float(region_value)) + 1}"


def _column_names(pan_graph_db, table_name):
    """Return the column names of an existing SQLite table without reading rows."""
    return [d[0] for d in pan_graph_db._exec(
        f"SELECT * FROM {table_name} LIMIT 0").description]


def _rewrite_table(pan_graph_db, table_name, structure, types, transform):
    """Drop and recreate ``table_name`` with the pinned v7 (structure, types),
    reinserting every row after passing its carryover dict through
    ``transform`` (which renames the component_id / region_id values in place)."""
    old_columns = _column_names(pan_graph_db, table_name)
    rows = pan_graph_db._exec(f"SELECT * FROM {table_name}").fetchall()

    new_rows = []
    for row in rows:
        carryover = {c: v for c, v in zip(old_columns, row)}
        transform(carryover)
        new_rows.append([carryover.get(col) for col in structure])

    pan_graph_db._exec(f"DROP TABLE {table_name}")
    pan_graph_db.create_table(table_name, structure, types)
    placeholders = ','.join(['?'] * len(structure))
    pan_graph_db._exec_many(
        f"INSERT INTO {table_name} VALUES ({placeholders})", new_rows)


def _transform_ids(carryover):
    # Build the region name from the ORIGINAL component + region values before
    # overwriting component_id (the region name embeds the component).
    region_new = _to_region_name(carryover.get('component_id'), carryover.get('region_id'))
    carryover['component_id'] = _to_component_name(carryover.get('component_id'))
    carryover['region_id'] = region_new


def _migrate_nodes(pan_graph_db):
    _rewrite_table(pan_graph_db, pan_graph_nodes_table_name,
                   pan_graph_nodes_table_structure, pan_graph_nodes_table_types, _transform_ids)


def _migrate_regions(pan_graph_db):
    _rewrite_table(pan_graph_db, pan_graph_regions_table_name,
                   pan_graph_regions_table_structure, pan_graph_regions_table_types, _transform_ids)


def _migrate_states(pan_graph_db):
    """Rename the active-component selector stored under
    ``graph_layout.component`` in every state so it matches the new component
    names. Returns the number of states rewritten."""
    states = pan_graph_db.get_table_as_dict('states')

    rewritten = 0
    for state_name, row in states.items():
        try:
            state = json.loads(row['content'])
        except (ValueError, TypeError, KeyError):
            run.warning(f"Could not parse pan-graph state '{state_name}' -- leaving it untouched.")
            continue

        graph_layout = state.get('graph_layout') if isinstance(state, dict) else None
        if not isinstance(graph_layout, dict) or 'component' not in graph_layout:
            continue

        graph_layout['component'] = _to_component_name(graph_layout['component'])
        pan_graph_db._exec(
            "UPDATE states SET content = ? WHERE name = ?",
            (json.dumps(state, indent=4), state_name),
        )
        rewritten += 1

    return rewritten


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

    progress.update("Renaming component/region ids in pan_graph_nodes (0 -> C_1 / R_1) ...")
    _migrate_nodes(pan_graph_db)

    progress.update("Renaming component/region ids in pan_graph_regions ...")
    _migrate_regions(pan_graph_db)

    progress.update("Renaming the active-component selector in states ...")
    _migrate_states(pan_graph_db)

    progress.update("Updating version")
    pan_graph_db.remove_meta_key_value_pair('version')
    pan_graph_db.set_version(next_version)

    progress.update("Committing changes")
    pan_graph_db.disconnect()

    progress.end()

    run.info_single(
        f"Your pan-graph database is now version {next_version}. Components and regions now use "
        f"1-based, named ids mirroring the gene-cluster convention: component_id \"0\" became "
        f"\"CP_0001\" (CP_0001 is still the largest component) and region \"0\" of component \"0\" "
        f"became \"CP_0001_1\" (the region id now embeds its component). Nothing else about the "
        f"graph changed.",
        nl_before=1, nl_after=1, mc='green')


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
