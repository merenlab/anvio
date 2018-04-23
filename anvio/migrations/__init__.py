#!/usr/bin/env python
# -*- coding: utf-8

import os
import importlib
from pathlib import Path

migration_scripts = {}

base_path = os.path.dirname(__file__)

for script_full_path in Path(base_path).glob('*/v*_to_v*.py'):
    script_full_path = str(script_full_path)
    script_path, script_filename = os.path.split(script_full_path)

    script_name = script_filename[:-3]
    db_type = os.path.basename(script_path)

    if not db_type in migration_scripts:
        migration_scripts[db_type] = {}

    spec = importlib.util.spec_from_file_location(script_name, script_full_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    
    migration_scripts[db_type][script_name] = module
