#!/usr/bin/env python

import sys
import json
import argparse

import networkx as nx

import anvio.dbinfo as dbinfo
import anvio.terminal as terminal

from anvio.errors import ConfigError

current_version, next_version = [x[1:] for x in __name__.split('_to_')]

run = terminal.Run()
progress = terminal.Progress()

# The target (v6) pan-graph table schemas are declared inline on purpose. A
# migration must pin the schema exactly as it was at the time of the version
# bump and must NOT import anvio.tables -- that module always tracks the current
# head and would silently drift as the schema evolves in later versions,
# breaking this script retroactively.
pan_graph_nodes_table_name      = 'pan_graph_nodes'
pan_graph_nodes_table_structure = ['node_id', 'node_type', 'region_id', 'gene_cluster_id', 'gene_calls_json', 'synteny_position_json', 'node_x', 'node_y', 'alignment_summary', 'component_id']
pan_graph_nodes_table_types     = [  'str'  ,    'str'   ,    'str'   ,       'str'      ,       'str'      ,          'str'         , 'numeric', 'numeric',       'str'       ,  'numeric'    ]

pan_graph_edges_table_name      = 'pan_graph_edges'
pan_graph_edges_table_structure = ['edge_id', 'source', 'target', 'weight' , 'genomes_json']
pan_graph_edges_table_types     = [  'str'  ,   'str' ,   'str' , 'numeric',     'str'      ]

pan_graph_regions_table_name      = 'pan_graph_regions'
pan_graph_regions_table_structure = ['component_id', 'region_id', 'region_type', 'x_min', 'x_max', 'num_synteny_gene_clusters', 'num_gene_clusters', 'num_gene_calls', 'max_expansion', 'min_expansion', 'complexity', 'complexity_normalized', 'diversity', 'diversity_normalized', 'weight', 'weight_normalized', 'composite_variability_score', 'complexity_mm_scaled', 'diversity_mm_scaled', 'expansion_mm_scaled', 'weight_mm_scaled']
pan_graph_regions_table_types     = [  'numeric'   ,    'str'   ,     'str'    , 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric']

# Every pan-graph produced by the master line was collapsed to its single
# largest connected component at generation time, so every v5 database is
# single-component -- component_id is 0 throughout.
SINGLE_COMPONENT_ID = 0

# The engine renamed several node-type values after v5 was cut, without a
# version bump, so a v5 database still carries the old names in two places that
# use *different* spellings for the same concept:
#   * the `node_type` column of pan_graph_nodes, and
#   * the `type_colors` keys inside each state's JSON.
# The v6 renderer only knows the new names, so normalise both here -- otherwise
# duplication/tRNA nodes render uncoloured and set_UI_settings crashes on the
# missing colour key. Both maps are applied idempotently (already-renamed
# values pass through untouched).
NODE_TYPE_RENAMES = {
    'paralog': 'duplication',
    'rearranged': 'rearrangement',
    'trna': 'rna',
}

STATE_TYPE_COLOR_KEY_RENAMES = {
    'multi_copy': 'duplication',
    'trna': 'rna',
}


def _column_names(pan_graph_db, table_name):
    """Return the column names of an existing SQLite table without reading rows."""
    return [d[0] for d in pan_graph_db._exec(
        f"SELECT * FROM {table_name} LIMIT 0").description]


def _norm_region_id(rid):
    """v5 stored region_id as a number; v6 stores a plain per-component string
    id. Normalise node.region_id and region.region_id through the same function
    so the (component_id, region_id) pairs still join after migration."""
    if rid is None or rid == '':
        return rid
    try:
        return str(int(float(rid)))
    except (ValueError, TypeError):
        return str(rid)


def _migrate_nodes(pan_graph_db):
    """v6 adds a `component_id` column to pan_graph_nodes, stores region_id as a
    plain string, and normalises legacy node_type values (paralog -> duplication,
    rearranged -> rearrangement, trna -> rna). Every v5 db is single-component,
    so component_id is 0."""
    old_columns = _column_names(pan_graph_db, pan_graph_nodes_table_name)
    if 'component_id' in old_columns:
        # Already migrated (e.g. a manual repair) -- skip silently.
        return

    rows = pan_graph_db._exec(f"SELECT * FROM {pan_graph_nodes_table_name}").fetchall()
    new_rows = []
    for row in rows:
        carryover = {c: v for c, v in zip(old_columns, row)}
        carryover['component_id'] = SINGLE_COMPONENT_ID
        carryover['region_id'] = _norm_region_id(carryover.get('region_id'))
        node_type = carryover.get('node_type')
        carryover['node_type'] = NODE_TYPE_RENAMES.get(node_type, node_type)
        new_rows.append([carryover.get(col) for col in pan_graph_nodes_table_structure])

    pan_graph_db._exec(f"DROP TABLE {pan_graph_nodes_table_name}")
    pan_graph_db.create_table(
        pan_graph_nodes_table_name,
        pan_graph_nodes_table_structure,
        pan_graph_nodes_table_types,
    )
    placeholders = ','.join(['?'] * len(pan_graph_nodes_table_structure))
    pan_graph_db._exec_many(
        f"INSERT INTO {pan_graph_nodes_table_name} VALUES ({placeholders})",
        new_rows,
    )


def _parse_directions(directions_raw):
    """The v5 `directions` column is a JSON object mapping each genome crossing
    the edge to its traversal orientation, 'L' (reverse) or 'R' (forward). The
    keys are therefore the edge's genome membership. Return a
    (genomes, orientations) tuple: `genomes` is the sorted list of genome names
    on the edge and `orientations` is the set of orientations present
    (upper-cased). Both are empty if the payload is missing or unparseable, and
    `genomes` is empty for the legacy list/scalar payloads that carried no
    genome names."""
    if not directions_raw:
        return [], set()
    try:
        parsed = json.loads(directions_raw)
    except (ValueError, TypeError):
        return [], set()
    if isinstance(parsed, dict):
        genomes = sorted(str(g) for g in parsed.keys())
        orientations = {str(v).upper() for v in parsed.values()}
    elif isinstance(parsed, (list, tuple)):
        genomes = []
        orientations = {str(v).upper() for v in parsed}
    else:
        genomes = []
        orientations = {str(parsed).upper()}
    return genomes, orientations


def _migrate_edges(pan_graph_db):
    """v6 has no per-edge orientation: it replaces the v5 `directions` column
    with `genomes_json` (the set of genomes crossing the edge) and treats every
    edge as forward source->target.

    A fully reversed edge -- one the v5 renderer drew reversed, i.e. with no
    forward ('R') orientation -- is traversed target->source in genome order.
    Since v6 cannot store orientation, such edges are RE-EXPRESSED as forward
    edges by swapping source and target rather than being dropped: dropping them
    would strand any node reachable only through a reversed edge, leaving it
    isolated in the migrated graph. Forward / mixed edges keep their source->
    target as-is.

    The v6 renderer requires a DAG, but flipping a reversed edge can close a
    directed cycle against the forward edges. So flipped edges are added on top
    of the forward graph one at a time, and a flipped edge is CUT only when
    adding it would actually create a cycle -- the minimum needed to stay
    acyclic while keeping the rest of each reversed run connected.

    Every surviving edge recovers its genome membership from the keys of the v5
    `directions` object, so genome-based edge colouring keeps working. (Legacy
    list/scalar `directions` payloads carried no genome names, so those edges get
    an empty genomes_json.)

    Returns (flipped_kept, cut) -- reversed edges re-oriented and kept, and
    reversed edges cut because they would have created a cycle."""
    old_columns = _column_names(pan_graph_db, pan_graph_edges_table_name)
    if 'genomes_json' in old_columns and 'directions' not in old_columns:
        # Already on the target shape.
        return 0, 0

    dir_idx = old_columns.index('directions') if 'directions' in old_columns else None
    rows = pan_graph_db._exec(f"SELECT * FROM {pan_graph_edges_table_name}").fetchall()

    # Split into forward/mixed edges (kept source->target) and fully reversed
    # edges (to be flipped). Each entry is the v6-shaped carryover dict.
    forward, reversed_edges = [], []
    for row in rows:
        carryover = {c: v for c, v in zip(old_columns, row)}
        genomes = []
        is_reversed = False
        if dir_idx is not None:
            genomes, orientations = _parse_directions(row[dir_idx])
            is_reversed = bool(orientations) and 'R' not in orientations
        carryover['genomes_json'] = json.dumps(genomes)
        (reversed_edges if is_reversed else forward).append(carryover)

    # Build the forward DAG, then graft flipped reversed edges onto it, cutting
    # only the ones that would introduce a cycle. `edge_id` order keeps the
    # choice of which edge to cut deterministic across runs.
    graph = nx.DiGraph()
    for e in forward:
        graph.add_edge(e['source'], e['target'])

    flipped_kept, cut = 0, 0
    kept_reversed = []
    for e in sorted(reversed_edges, key=lambda e: e['edge_id']):
        new_source, new_target = e['target'], e['source']  # flip endpoints
        # Adding new_source->new_target closes a cycle iff a path already runs
        # new_target -> ... -> new_source.
        if new_source in graph and new_target in graph and nx.has_path(graph, new_target, new_source):
            cut += 1
            continue
        graph.add_edge(new_source, new_target)
        e['source'], e['target'] = new_source, new_target
        kept_reversed.append(e)
        flipped_kept += 1

    new_rows = [[e.get(col) for col in pan_graph_edges_table_structure]
                for e in forward + kept_reversed]

    pan_graph_db._exec(f"DROP TABLE {pan_graph_edges_table_name}")
    pan_graph_db.create_table(
        pan_graph_edges_table_name,
        pan_graph_edges_table_structure,
        pan_graph_edges_table_types,
    )
    placeholders = ','.join(['?'] * len(pan_graph_edges_table_structure))
    pan_graph_db._exec_many(
        f"INSERT INTO {pan_graph_edges_table_name} VALUES ({placeholders})",
        new_rows,
    )
    return flipped_kept, cut


def _migrate_regions(pan_graph_db):
    """v6 prepends a `component_id` column to pan_graph_regions and stores
    region_id as a plain string. Single-component v5 dbs get component_id 0."""
    old_columns = _column_names(pan_graph_db, pan_graph_regions_table_name)
    if 'component_id' in old_columns:
        return

    rows = pan_graph_db._exec(f"SELECT * FROM {pan_graph_regions_table_name}").fetchall()
    new_rows = []
    for row in rows:
        carryover = {c: v for c, v in zip(old_columns, row)}
        carryover['component_id'] = SINGLE_COMPONENT_ID
        carryover['region_id'] = _norm_region_id(carryover.get('region_id'))
        new_rows.append([carryover.get(col) for col in pan_graph_regions_table_structure])

    pan_graph_db._exec(f"DROP TABLE {pan_graph_regions_table_name}")
    pan_graph_db.create_table(
        pan_graph_regions_table_name,
        pan_graph_regions_table_structure,
        pan_graph_regions_table_types,
    )
    placeholders = ','.join(['?'] * len(pan_graph_regions_table_structure))
    pan_graph_db._exec_many(
        f"INSERT INTO {pan_graph_regions_table_name} VALUES ({placeholders})",
        new_rows,
    )


def _migrate_states(pan_graph_db):
    """Rename the legacy node type_colors keys (multi_copy -> duplication,
    trna -> rna) inside every stored state so the v6 renderer finds a colour for
    each node type. States are left untouched if they cannot be parsed or carry
    no type_colors block, and keys that are already on the v6 names pass through.

    Returns the number of states whose type_colors were rewritten."""
    states = pan_graph_db.get_table_as_dict('states')

    rewritten = 0
    for state_name, row in states.items():
        try:
            state = json.loads(row['content'])
        except (ValueError, TypeError, KeyError):
            run.warning(f"Could not parse pan-graph state '{state_name}' -- leaving it untouched.")
            continue

        nodes_block = state.get('nodes') if isinstance(state, dict) else None
        type_colors = nodes_block.get('type_colors') if isinstance(nodes_block, dict) else None
        if not isinstance(type_colors, dict):
            continue

        changed = False
        for old_key, new_key in STATE_TYPE_COLOR_KEY_RENAMES.items():
            if old_key in type_colors:
                # Don't clobber a v6 key that already exists; just drop the legacy one.
                type_colors.setdefault(new_key, type_colors[old_key])
                del type_colors[old_key]
                changed = True

        if changed:
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
    else:
         raise ConfigError(
            "Due to major changes in the pan-graph engine, migrating your pangenome-graph-db is not supported."
            "We therefore kindly ask you to re-run the pan genome graph generation on your data instead.")

    pan_graph_db = pan_graph_db_info.load_db()

    progress.new("Migrating")

    progress.update("Rewriting pan_graph_nodes (adding component_id = 0; normalising node types) ...")
    _migrate_nodes(pan_graph_db)

    progress.update("Rewriting pan_graph_edges (flipping reversed edges; directions -> genomes_json) ...")
    reversed_edges_flipped, reversed_edges_cut = _migrate_edges(pan_graph_db)

    progress.update("Rewriting pan_graph_regions (adding component_id = 0) ...")
    _migrate_regions(pan_graph_db)

    progress.update("Renaming legacy node type_colors keys in states ...")
    _migrate_states(pan_graph_db)

    progress.update("Updating version")
    pan_graph_db.remove_meta_key_value_pair('version')
    pan_graph_db.set_version(next_version)

    progress.update("Committing changes")
    pan_graph_db.disconnect()

    progress.end()

    run.info('Reversed edges re-oriented', reversed_edges_flipped)
    run.info('Reversed edges cut (would have created a cycle)', reversed_edges_cut,
             mc='red' if reversed_edges_cut else 'green', nl_after=1)

    run.warning(
        f"Your pan-graph database is now version {next_version}, but please read this carefully. "
        f"The whole pan-graph has been assigned to a single component (component_id = 0): every "
        f"pan-graph produced by the older (master) engine was already reduced to its single largest "
        f"connected component, so this is faithful. However, {reversed_edges_flipped} reversed "
        f"edge(s) (edges the old engine traversed in the reverse orientation) were RE-ORIENTED during "
        f"migration by swapping their endpoints, because the new engine has no notion of edge "
        f"orientation -- this keeps every node connected instead of stranding nodes that were only "
        f"reachable in reverse. Of those, {reversed_edges_cut} edge(s) had to be CUT because "
        f"re-orienting them would have created a cycle (the new engine requires an acyclic graph); the "
        f"rest of each affected run stays connected. Edges recovered their per-edge genome membership "
        f"from the old `directions` payload, so genome-based edge colouring still works. Because the "
        f"two engine versions differ "
        f"substantially, this migrated database is a best-effort, structural upgrade only -- anvi'o "
        f"strongly recommends that you re-run the entire pan-graph generation on your data rather "
        f"than relying on this migrated database for analysis.",
        header=f"PAN-GRAPH DB MIGRATED TO v{next_version} -- PLEASE RE-RUN THE PAN-GRAPH",
        lc='yellow'
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
