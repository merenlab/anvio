# Prompt: make `anvio.TABULATE` decide on unicode the same way the taxonomy tree does

> This file is a note-to-self to kick off a fresh Claude Code session. It is not part of the
> codebase -- delete it once the work it describes is done.

---

## The task

Please read `ARCHITECTURE.md` first.

`anvio.TABULATE()` in `anvio/__init__.py:154` is supposed to fall back from a fancy unicode table to
an ASCII one when the output stream cannot encode box-drawing characters. It gets the check wrong,
so in practice it *always* falls back, and no anvi'o program has ever shown the fancy table. Fix it
by reusing the helper that was added to `anvio/terminal.py` for exactly this question, and make sure
the fancy tables actually look right everywhere `TABULATE` is used.

This is deliberately a separate PR from the `--tree-output` work that introduced the helper.

## The bug

```python
# anvio/__init__.py:154
def TABULATE(table, header, numalign="right", max_width=0):
    """Encoding-safe `tabulate`"""

    from tabulate import tabulate

    tablefmt = "fancy_grid" if sys.stdout.encoding == "UTF-8" else "grid"
```

`sys.stdout.encoding` is reported by Python as `'utf-8'` (lowercase, with the hyphen) on macOS and
Linux, so `== "UTF-8"` is False and `tablefmt` is always `"grid"`. Verify for yourself:

```bash
python -c "import sys; print(repr(sys.stdout.encoding))"                        # -> 'utf-8'
python -c "import sys; print('fancy_grid' if sys.stdout.encoding == 'UTF-8' else 'grid')"   # -> grid
```

So today every anvi'o table renders like the left, when the intent was the right:

```
+----------+-----+          ╒══════════╤═════╕
| taxon    |   n |          │ taxon    │   n │
+==========+=====+          ╞══════════╪═════╡
| Bacteria |  94 |          │ Bacteria │  94 │
+----------+-----+          ├──────────┼─────┤
| Archaea  |   2 |          │ Archaea  │   2 │
+----------+-----+          ╘══════════╧═════╛
```

## The fix

Use `terminal.stdout_supports_unicode()`, which normalizes the encoding spelling and tolerates
`sys.stdout.encoding` being `None`:

```python
tablefmt = "fancy_grid" if terminal.stdout_supports_unicode() else "grid"
```

**If that helper does not exist yet** (i.e. the `--tree-output` PR has not landed), add it to
`anvio/terminal.py` next to `get_terminal_size()`, and keep the docstring's warning that this is not
the same question as `sys.stdout.isatty()` -- see "Do not use isatty" below.

### Gotcha: the import has to be deferred

`anvio/terminal.py` does `import anvio` at its top, so a module-level `import anvio.terminal` in
`anvio/__init__.py` is circular and will not work. Put the import inside the function body, next to
the `from tabulate import tabulate` that is already deferred there for the same kind of reason, and
leave a short comment saying why. This has been confirmed to work:

```python
def TABULATE(table, header, numalign="right", max_width=0):
    """Encoding-safe `tabulate`"""

    from tabulate import tabulate

    # deferred because `anvio.terminal` imports `anvio`, so importing it at module level here
    # would be circular
    import anvio.terminal as terminal
```

### Do not use `isatty` here

It is tempting, and it is wrong. `sys.stdout.isatty()` answers "is a human watching a terminal?"
(which is what `ttycolors.color_text` and the `Progress` class correctly ask before emitting ANSI
escapes or live progress lines). Whether a stream can *encode* `╒` is a different question:

- `anvi-something > table.txt` is not a TTY, but the file holds box-drawing characters perfectly
  well, and the pretty table is what you want in it.
- An ASCII terminal *is* a TTY, and printing `╒` there raises `UnicodeEncodeError`.

## Blast radius

`TABULATE` has 9 call sites across 4 files, and this changes the terminal output of all of them:

```
anvio/taxonomyops/__init__.py:686, 989, 1631, 1677
anvio/cli/estimate_genome_completeness.py:104, 149, 201
anvio/cli/tabulate.py:26
```

Two things to check rather than assume:

1. **The `max_width` truncation path** (`anvio/__init__.py:162-167`) slices each line and re-appends
   `l[-2:]`. The border characters differ between `grid` and `fancy_grid`, so confirm truncated
   tables still close cleanly. Exercise it with `anvi-tabulate ... --max-width 60`.
2. **Column alignment with wide/CJK characters.** `terminal.tabulate()` (`anvio/terminal.py:1030`)
   exists specifically to patch `tabulate`'s width handling for ANSI codes and line breaks, and
   `TABULATE` does *not* go through it. Worth a look at whether it should, though that may belong in
   yet another PR -- don't scope-creep.

No test breakage expected: the `diff`/`cmp` calls in `anvio/tests/*.sh` all compare TAB-delimited
output files, not terminal tables.

## Optional companion fix

`anvio/programsearch.py:104` asks the same question a third way:

```python
if sys.stdout.encoding.lower() != "utf-8":
```

This handles the casing but raises `AttributeError` when `sys.stdout.encoding` is `None` (which
happens when stdout has been replaced with a stream that does not declare one). Switching it to the
helper folds all three call sites into one implementation. Small and worth doing here, but ask
before bundling it if you would rather keep the PR to one file.

## Verification

```bash
# the tables should now be fancy
anvi-estimate-scg-taxonomy -c Sample_01.db --metagenome-mode | head -20
anvi-estimate-genome-completeness -c Sample_01.db 2>&1 | tail -20

# truncation still closes cleanly
anvi-tabulate <some-tab-delimited-file> --max-width 60

# still degrades to ASCII when it must
python -c "
import io, sys, anvio
sys.stdout = io.TextIOWrapper(io.BytesIO(), encoding='ascii')
anvio.TABULATE([['a', 1]], ['x', 'y'])
sys.stdout.seek(0); sys.stderr.write(sys.stdout.buffer.getvalue().decode('ascii'))
"

ruff check .
```

## Already checked, so you don't have to

- `sys.stdout.encoding` is `'utf-8'` on this machine, in a terminal and through a pipe, and even
  under `LANG=C LC_ALL=C` -- so the current check is dead in every ordinary case.
- A deferred `import anvio.terminal as terminal` inside `TABULATE` works.
- A module-level one does not (`terminal.py` imports `anvio`).
- Nothing in `anvio/tests/` diffs terminal table output.
- `anvio.AS_MARKDOWN` does not flow through `TABULATE` at all -- markdown output is handled
  separately in `anvio/utils.py:1080` and `anvio/filesnpaths.py:719`, so it is out of scope here.
