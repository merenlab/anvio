#!/usr/bin/env python
# -*- coding: utf-8

import os
import importlib
from pathlib import Path

migration_scripts = {}

base_path = os.path.dirname(__file__)

for script_file in Path(base_path).glob('*/v*_to_v*.py'):
    script_path, script_name = os.path.split(script_file)

    db_type = os.path.dirname(script_path)

    if not db_type in migration_scripts:
        migration_scripts[db_type] = {}

    spec = importlib.util.spec_from_file_location(script_name[:-3], script_file)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    
    migration_scripts[db_type][script_name] = module
