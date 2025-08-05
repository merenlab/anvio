#!/usr/bin/env python3

import os
import re
from collections import defaultdict
from typing import Dict, Set, List, Tuple

from utils_migration_map import FUNCTION_TO_MODULE

def find_all_releavnt_function_calls(content: str) -> Set[str]:
    """Find all utils.function_name calls in the content, excluding comments and docstrings."""
    lines = content.split('\n')
    pattern = r'\butils\.([a-zA-Z_][a-zA-Z0-9_]*)'
    matches = set()

    in_docstring = False
    docstring_delimiter = None

    for line in lines:
        stripped_line = line.strip()

        # Handle docstring detection
        if not in_docstring:
            # Check if line starts a docstring
            if (stripped_line.startswith('"""') or stripped_line.startswith("'''")):
                docstring_delimiter = stripped_line[:3]
                # Check if it's a single-line docstring
                if stripped_line.count(docstring_delimiter) >= 2 and len(stripped_line) > 3:
                    # Single line docstring, don't set in_docstring flag
                    continue
                else:
                    in_docstring = True
                    continue
        else:
            # We're in a docstring, check if this line ends it
            if docstring_delimiter in line:
                in_docstring = False
                docstring_delimiter = None
            continue

        # Skip comment lines
        if stripped_line.startswith('#'):
            continue

        # Process the line for utils calls, but remove inline comments first
        comment_pos = line.find('#')
        if comment_pos != -1:
            # Only process the part before the comment
            line_to_process = line[:comment_pos]
        else:
            line_to_process = line

        # Find utils function calls in this line
        line_matches = re.findall(pattern, line_to_process)
        matches.update(line_matches)

    return matches


def find_utils_imports(content: str) -> List[str]:
    """Find import statements that import anvio.utils as utils."""
    lines = content.split('\n')
    utils_import_patterns = [
        r'^\s*import\s+anvio\.utils\s+as\s+utils\s*$',
        r'^\s*from\s+anvio\s+import\s+utils\s*$',
        r'^\s*import\s+anvio\.utils\s*$'
    ]

    import_lines = []
    for i, line in enumerate(lines):
        for pattern in utils_import_patterns:
            if re.match(pattern, line):
                import_lines.append(i)
                break

    return import_lines

def convert_module_path_to_import(module_path: str) -> str:
    """Convert 'anvio/utils/sequences' to 'anvio.utils.sequences'."""
    return module_path.replace('/', '.')

def generate_new_imports(functions_used: Set[str]) -> Dict[str, List[str]]:
    """Generate new import statements grouped by module."""
    imports_by_module = defaultdict(list)

    for func in functions_used:
        if func in FUNCTION_TO_MODULE:
            module = FUNCTION_TO_MODULE[func]
            imports_by_module[module].append(func)
        else:
            print(f"Warning: Function '{func}' not found in FUNCTION_TO_MODULE mapping")

    return dict(imports_by_module)

def create_import_statements(imports_by_module: Dict[str, List[str]]) -> List[str]:
    """Create the actual import statement strings."""
    import_statements = []

    for module, functions in sorted(imports_by_module.items()):
        module_import = convert_module_path_to_import(module)
        if len(functions) == 1:
            import_statements.append(f"from {module_import} import {functions[0]}")
        else:
            # Sort functions for consistent output
            sorted_functions = sorted(functions)
            if len(sorted_functions) <= 3:  # Keep short imports on one line
                func_list = ", ".join(sorted_functions)
                import_statements.append(f"from {module_import} import {func_list}")
            else:  # Multi-line import for readability
                import_statements.append(f"from {module_import} import (")
                for i, func in enumerate(sorted_functions):
                    if i == len(sorted_functions) - 1:
                        import_statements.append(f"    {func}")
                    else:
                        import_statements.append(f"    {func},")
                import_statements.append(")")

    return import_statements

def find_import_insertion_point(lines: List[str]) -> int:
    """Find the best place to insert new imports (after existing imports)."""
    last_import_line = -1
    in_docstring = False
    docstring_quotes = None

    for i, line in enumerate(lines):
        stripped = line.strip()

        # Handle docstrings
        if not in_docstring:
            if stripped.startswith('"""') or stripped.startswith("'''"):
                docstring_quotes = stripped[:3]
                if not (stripped.endswith(docstring_quotes) and len(stripped) > 3):
                    in_docstring = True
                continue
        else:
            if docstring_quotes in line:
                in_docstring = False
                continue
            continue

        # Skip comments and empty lines at the beginning
        if not stripped or stripped.startswith('#'):
            continue

        # Check if this is an import statement
        if (stripped.startswith('import ') or
            stripped.startswith('from ') or
            (i > 0 and lines[i-1].strip().endswith('\\'))):  # Continuation of import
            last_import_line = i
        elif last_import_line >= 0:
            # We've found imports and now hit a non-import line
            break

    # Insert after the last import, or at the beginning if no imports found
    return last_import_line + 1 if last_import_line >= 0 else 0

def refactor_file(file_path: str, dry_run: bool = True) -> Tuple[bool, str]:
    """Refactor a single Python file."""
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()
    except Exception as e:
        return False, f"Error reading file: {e}"

    # Find utils function calls
    functions_used = find_all_releavnt_function_calls(content)

    if not functions_used:
        return False, "No utils function calls found"

    # Find utils import lines
    utils_import_lines = find_utils_imports(content)

    lines = content.split('\n')

    # Remove old utils import statements
    new_lines = []
    for i, line in enumerate(lines):
        if i not in utils_import_lines:
            new_lines.append(line)

    # Replace utils.function_name with function_name
    new_content = '\n'.join(new_lines)
    for func in functions_used:
        if func in FUNCTION_TO_MODULE:
            # Use word boundaries to avoid partial matches
            pattern = r'\butils\.' + re.escape(func) + r'\b'
            new_content = re.sub(pattern, func, new_content)

    # Generate new import statements
    imports_by_module = generate_new_imports(functions_used)
    new_import_statements = create_import_statements(imports_by_module)

    # Find where to insert new imports
    new_lines = new_content.split('\n')
    insertion_point = find_import_insertion_point(new_lines)

    # Insert new imports
    for i, import_stmt in enumerate(reversed(new_import_statements)):
        if import_stmt == ")" and i == 0:  # Handle multi-line imports
            new_lines.insert(insertion_point, import_stmt)
        elif import_stmt.startswith("    "):  # Handle multi-line import items
            new_lines.insert(insertion_point, import_stmt)
        elif import_stmt.startswith("from") and import_stmt.endswith("("):
            new_lines.insert(insertion_point, import_stmt)
        else:
            new_lines.insert(insertion_point, import_stmt)

    final_content = '\n'.join(new_lines)

    if not dry_run:
        with open(file_path, 'w', encoding='utf-8') as f:
            f.write(final_content)

    summary = f"Functions found: {', '.join(sorted(functions_used))}"
    if utils_import_lines:
        summary += f"\nRemoved {len(utils_import_lines)} utils import line(s)"
    summary += f"\nAdded {len(new_import_statements)} new import statement(s)"

    return True, summary

def find_python_files(directory: str, exclude_dirs: Set[str] = None) -> List[str]:
    """Find all Python files in directory, excluding specified directories."""
    if exclude_dirs is None:
        exclude_dirs = {'data', '__pycache__', '.git', 'build', 'dist'}

    exclude_dirs.add('utils')

    python_files = []
    for root, dirs, files in os.walk(directory):
        # Remove excluded directories from dirs in-place
        dirs[:] = [d for d in dirs if d not in exclude_dirs]

        for file in files:
            if file.endswith('.py'):
                python_files.append(os.path.join(root, file))

    return python_files

def main():
    import argparse

    parser = argparse.ArgumentParser(description='Refactor anvio utils imports')
    parser.add_argument('directory', help='Directory to process (default: current directory)',
                       nargs='?', default='.')
    parser.add_argument('--dry-run', action='store_true',
                       help='Show what would be changed without making changes')
    parser.add_argument('--exclude-dirs', nargs='*', default=['data'],
                       help='Directories to exclude (default: data)')

    args = parser.parse_args()

    exclude_dirs = set(args.exclude_dirs)
    python_files = find_python_files(args.directory, exclude_dirs)

    print(f"Found {len(python_files)} Python files to process")
    if args.dry_run:
        print("DRY RUN MODE - No files will be modified")
    print("-" * 80)

    processed = 0
    modified = 0

    for file_path in python_files:
        try:
            changed, summary = refactor_file(file_path, dry_run=args.dry_run)
            processed += 1

            if changed:
                modified += 1
                print(f"\n✓ {file_path}")
                print(f"  {summary}")
            else:
                if summary != "No utils function calls found":
                    print(f"\n✗ {file_path}")
                    print(f"  {summary}")

        except Exception as e:
            print(f"\n✗ Error processing {file_path}: {e}")

    print("-" * 80)
    print(f"Processed: {processed} files")
    print(f"Modified: {modified} files")
    if args.dry_run:
        print("\nRun without --dry-run to apply changes")

if __name__ == '__main__':
    main()
