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

KNOWN_FIXED_KEYS = frozenset({
    'rearranged_color', 'accessory_color', 'paralog_color', 'singleton_color', 'core_color',
    'trna_color', 'layer_color', 'non_back_color', 'back_color', 'flexsaturation',
    'arrow', 'flexarrow', 'globalbackbone', 'flexglobalbackbone', 'flextree',
    'tree_length', 'tree_offset', 'tree_thickness', 'distx', 'disty', 'size', 'circ',
    'edge', 'flexlinear', 'line', 'label', 'label_offset', 'search_hit',
    'inner_margin', 'outer_margin', 'inner', 'start_angle', 'end_angle', 'num_position',
    'flexcondtr', 'condtr', 'flexmaxlength', 'maxlength', 'flexgroupcompress',
    'groupcompress', 'track_line_width', 'region_label_size', 'region_label_min_width',
    'region_label_distance', 'flexbinlabels', 'bin_label_orientation', 'bin_label_size',
    'bin_ring_height', 'bin_ring_opacity', 'flexbinedges', 'bin_edge_thickness',
    'bin_edge_color', 'bin_edge_opacity'
})


def convert_state_v4_to_v5(old_state, genome_names):
    s = old_state

    drawing_type = 'linear' if s.get('flexlinear', False) else 'circular'

    genome_names = [g for g in genome_names if g]
    genome_related_keys = set()
    for genome in genome_names:
        genome_related_keys.update({genome, 'flex' + genome, genome + 'layer', 'flex' + genome + 'layer'})

    genomes = {}
    for genome in genome_names:
        genomes[genome] = {
            'color': s.get(genome, '#000000'),
            'show': s.get('flex' + genome, True),
            'track_height': s.get(genome + 'layer', 7),
            'show_track': s.get('flex' + genome + 'layer', True)
        }

    all_handled = KNOWN_FIXED_KEYS | genome_related_keys
    remaining_keys = set(s.keys()) - all_handled
    flex_layer_names = {k[4:] for k in remaining_keys if k.startswith('flex')}
    imported_layers = {}
    for layer in flex_layer_names:
        if layer in s:
            imported_layers[layer] = {
                'visible': s.get('flex' + layer, False),
                'height': s.get(layer, 0)
            }

    max_edge_length = s.get('maxlength', 1000)
    if max_edge_length == -1:
        max_edge_length = 1000

    return {
        'drawing': {
            'type': drawing_type,
            'inner_radius': s.get('inner', 41),
            'start_angle': s.get('start_angle', 0),
            'end_angle': s.get('end_angle', 270),
            'node_x_spacing': s.get('distx', 45),
            'node_y_spacing': s.get('disty', 120)
        },
        'nodes': {
            'radius': s.get('size', 15),
            'outline_width': s.get('circ', 5),
            'fade_by_prevalence': s.get('flexsaturation', True),
            'type_colors': {
                'core': s.get('core_color', '#BCBCBC'),
                'rearrangement': s.get('rearranged_color', '#8FF0A4'),
                'accessory': s.get('accessory_color', '#DC8ADD'),
                'multi_copy': s.get('paralog_color', '#FFA348'),
                'singleton': s.get('singleton_color', '#99C1F1'),
                'trna': s.get('trna_color', '#F66151')
            }
        },
        'edges': {
            'width': s.get('edge', 5)
        },
        'graph_layout': {
            'grouping_enabled': s.get('flexcondtr', False),
            'grouping_threshold': s.get('condtr', -1),
            'max_edge_length_enabled': s.get('flexmaxlength', True),
            'max_edge_length': max_edge_length,
            'group_compression_enabled': s.get('flexgroupcompress', False),
            'group_compression': s.get('groupcompress', 1.0)
        },
        'layers': {
            'backbone': {
                'visible': s.get('flexglobalbackbone', True),
                'height': s.get('globalbackbone', 3),
                'backbone_color': s.get('back_color', '#3D70A0'),
                'variable_region_color': s.get('non_back_color', '#F8E45C')
            },
            'orientation_arrow': {
                'visible': s.get('flexarrow', True),
                'height': s.get('arrow', 3)
            },
            'search': {
                'hit_height': s.get('search_hit', 3)
            }
        },
        'layers_tree': {
            'visible': s.get('flextree', True),
            'height': s.get('tree_length', 44),
            'offset': s.get('tree_offset', 1),
            'line_width': s.get('tree_thickness', 10)
        },
        'genome_tracks': {
            'line_width': s.get('track_line_width', 5),
            'background_color': s.get('layer_color', '#F5F5F5'),
            'genomes': genomes
        },
        'imported_layers': imported_layers,
        'labels': {
            'font_size': s.get('label', 0),
            'offset': s.get('label_offset', 1),
            'position_tick_count': s.get('num_position', 20)
        },
        'margins': {
            'inner': s.get('inner_margin', 3),
            'outer': s.get('outer_margin', 0)
        },
        'region_labels': {
            'font_size': s.get('region_label_size', 13),
            'min_width_px': s.get('region_label_min_width', 80),
            'distance': s.get('region_label_distance', 2)
        },
        'bins': {
            'show_labels': s.get('flexbinlabels', True),
            'label_orientation': s.get('bin_label_orientation', 'natural'),
            'label_font_size': s.get('bin_label_size', 19.5),
            'ring_height': s.get('bin_ring_height', 4),
            'ring_opacity': s.get('bin_ring_opacity', 0.8),
            'show_edges': s.get('flexbinedges', True),
            'edge_thickness': s.get('bin_edge_thickness', 4),
            'edge_color': s.get('bin_edge_color', '#FFFFFF'),
            'edge_opacity': s.get('bin_edge_opacity', 1.0)
        }
    }


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
    progress.update("Reading genome names ...")
    genome_names_str = pan_graph_db.get_meta_value('genome_names')
    genome_names = genome_names_str.split(',') if genome_names_str else []

    progress.update("Converting state JSON keys from flat (v4) to nested (v5) structure ...")
    states_table = pan_graph_db.get_table_as_dict('states')

    for state_name, row in states_table.items():
        try:
            old_state = json.loads(row['content'])
        except (json.JSONDecodeError, KeyError):
            run.warning(f"Could not parse state '{state_name}' — skipping.")
            continue

        new_state = convert_state_v4_to_v5(old_state, genome_names)
        pan_graph_db._exec(
            "UPDATE states SET content = ? WHERE name = ?",
            (json.dumps(new_state, indent=4), state_name)
        )

    progress.update("Updating version")
    pan_graph_db.remove_meta_key_value_pair('version')
    pan_graph_db.set_version(next_version)

    progress.update("Committing changes")
    pan_graph_db.disconnect()

    progress.end()

    run.info_single(
        f"Your pan-graph database is now version {next_version}. State variable names have been "
        f"restructured from a flat layout to a nested, human-readable format.",
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
