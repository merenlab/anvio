#!/usr/bin/env python
# -*- coding: utf-8

import sys
import os
import json
from urllib.parse import unquote

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
from anvio.errors import ConfigError, FilesNPathsError
import anvio.filesnpaths as filesnpaths

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['ahenoch']
__provides__ = ['config-file']
__requires__ = ['config-file']
__description__ = "Modify the config-json file or any other file in JSON format."

# Most of these functions wer created with the help of AI.

def main():
    args = get_args()

    try:  
        run_program()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def find_patch_leaf_path(patch):
    
    if not isinstance(patch, dict) or not patch:
        raise ConfigError("Patch must be a non-empty dict.")

    path = []
    current = patch

    while isinstance(current, dict):
        if len(current) != 1:
            raise ValueError("patch must have exactly one key at each level")
        key = next(iter(current))
        path.append(key)
        current = current[key]

    # current is the leaf value
    return path, current


def try_update_from_patch(node, patch_path, new_value, base_path):
    
    if not isinstance(node, dict):
        return None

    current = node
    last_i = len(patch_path) - 1

    for i, key in enumerate(patch_path):
        if not isinstance(current, dict) or key not in current:
            return None

        if i == last_i:
            current[key] = new_value
            return base_path + [key]

        current = current[key]

    return None


def global_update(root, patch, run):
    
    patch_path, new_value = find_patch_leaf_path(patch)

    changed_paths = []
    stack = [(root, [])]  # (node, path_to_node)

    while stack:
        node, node_path = stack.pop()

        changed = try_update_from_patch(node, patch_path, new_value, node_path)
        if changed is not None:
            changed_paths.append("/".join(str(p) for p in changed))

        if isinstance(node, dict):
            for k, v in node.items():
                if isinstance(v, (dict, list)):
                    stack.append((v, node_path + [k]))
        elif isinstance(node, list):
            for idx, item in enumerate(node):
                if isinstance(item, (dict, list)):
                    stack.append((item, node_path + [idx]))

    run.info_single(f"Successfully changed JSON at paths {changed_paths} to {new_value}")


def local_update(target, patch, run, only_existing=True, path=()):
    for key, value in patch.items():
        current_path = path + (key,)

        if isinstance(value, dict):
            if key not in target:
                if only_existing:
                    raise ConfigError(f"Missing key at path {'/'.join(current_path)}")
                target[key] = {}

            if not isinstance(target[key], dict):
                raise ConfigError(
                    f"Expected dict at path {'/'.join(current_path)}, "
                    f"found {type(target[key]).__name__}"
                )

            local_update(
                target[key],
                value,
                run,
                only_existing=only_existing,
                path=current_path,
            )

        else:
            if key not in target and only_existing:
                raise ConfigError(f"Missing key at path {'/'.join(current_path)}")

            old_value = target[key]
            target[key] = value
            run.info_single(f"Successfully changed JSON at path {'/'.join(current_path)} from {old_value} to {value}")


def run_program():
    args = get_args()
    run = terminal.Run()

    filesnpaths.is_file_json_formatted(args.config_file)
    filesnpaths.is_output_file_writable(args.config_file, ok_if_exists=True)

    config_file_dict = json.load(open(args.config_file, 'r'))

    if args.local_change:
        for string in args.local_change:
            json_string = unquote(string)
            json_dict = json.loads(json_string)
            
            local_update(config_file_dict, json_dict, run)

    if args.global_change:
        for string in args.global_change:
            json_string = unquote(string)
            json_dict = json.loads(json_string)
            
            global_update(config_file_dict, json_dict, run)

    json.dump(config_file_dict, open(args.config_file, 'w'), indent=4, ensure_ascii=False)
    run.info_single(f"Successfully saved JSON")


def get_args():
    from anvio.argparse import ArgumentParser
    parser = ArgumentParser(description=__description__)

    parser.add_argument(*anvio.A('config-file'), **anvio.K('config-file', {'required': True}))
    parser.add_argument('--local-change', nargs='+', type=str, help = "JSON configuration string necessary ()")
    parser.add_argument('--global-change', nargs='+', type=str, help = "JSON configuration string necessary ()")

    return parser.get_args(parser)


if __name__ == '__main__':
    main()