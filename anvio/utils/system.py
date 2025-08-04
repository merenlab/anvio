import os
import shutil
import linecache
import tracemalloc

import anvio

from anvio.errors import ConfigError
from anvio.terminal import Run, SuppressAllOutput

from anvio.utils.files import human_readable_file_size

with SuppressAllOutput():
    try:
        import psutil
        PSUTIL_OK=True
    except:
        PSUTIL_OK=False


def get_total_memory_usage(keep_raw=False):
    """Get the total memory, including children

    Parameters
    ==========
    keep_raw : bool, False
        A human readable format is returned, e.g. "1.41 GB". If keep_raw, the raw number is
        returned, e.g. 1515601920
    """
    if not PSUTIL_OK:
        return None

    current_process = psutil.Process(os.getpid())
    mem = current_process.memory_info().rss
    for child in current_process.children(recursive=True):
        try:
            mem += child.memory_info().rss
        except:
            pass

    return mem if keep_raw else human_readable_file_size(mem)



def display_top_memory_usage(snapshot, key_type='lineno', limit=10):
    """A pretty-print for the tracemalloc memory usage module

    Modified from https://docs.python.org/3/library/tracemalloc.html

    Examples
    ========
    >>> import tracemalloc
    >>> import anvio.utils as utils
    >>> tracemalloc.start()
    >>> snap = tracemalloc.take_snapshot
    >>> utils.display_top_memory_usage(snap)
    Top 10 lines
    #1: anvio/bamops.py:160: 4671.3 KiB
        constants.cigar_consumption,
    #2: anvio/bamops.py:96: 2571.6 KiB
        self.cigartuples = np.array(read.cigartuples)
    #3: python3.6/linecache.py:137: 1100.0 KiB
        lines = fp.readlines()
    #4: <frozen importlib._bootstrap_external>:487: 961.4 KiB
    #5: typing/templates.py:627: 334.3 KiB
        return type(base)(name, (base,), dct)
    #6: typing/templates.py:923: 315.7 KiB
        class Template(cls):
    #7: python3.6/_weakrefset.py:84: 225.2 KiB
        self.data.add(ref(item, self._remove))
    #8: targets/npyimpl.py:411: 143.2 KiB
        class _KernelImpl(_Kernel):
    #9: _vendor/pyparsing.py:3349: 139.7 KiB
        self.errmsg = "Expected " + _ustr(self)
    #10: typing/context.py:456: 105.1 KiB
        def on_disposal(wr, pop=self._globals.pop):
    3212 other: 4611.9 KiB
    Total allocated size: 15179.4 KiB
    """

    snapshot = snapshot.filter_traces((
        tracemalloc.Filter(False, "<frozen importlib._bootstrap>"),
        tracemalloc.Filter(False, "<unknown>"),
    ))
    top_stats = snapshot.statistics(key_type)

    print("Top %s lines" % limit)
    for index, stat in enumerate(top_stats[:limit], 1):
        frame = stat.traceback[0]
        # replace "/path/to/module/file.py" with "module/file.py"
        filename = os.sep.join(frame.filename.split(os.sep)[-2:])
        print("#%s: %s:%s: %.1f KiB"
              % (index, filename, frame.lineno, stat.size / 1024))
        line = linecache.getline(frame.filename, frame.lineno).strip()
        if line:
            print('    %s' % line)

    other = top_stats[limit:]
    if other:
        size = sum(stat.size for stat in other)
        print("%s other: %.1f KiB" % (len(other), size / 1024))
    total = sum(stat.size for stat in top_stats)
    print("Total allocated size: %.1f KiB" % (total / 1024))



def get_available_program_names_in_active_environment(prefix=None, contains=None, postfix=None):
    """Find all executable programs in the current environment that match the given criteria.

    Parameters
    ==========
    prefix : str, optional
        The prefix to search for (e.g., 'anvi-')
    contains : str, optional
        String that must be contained in the program name (e.g., 'graph')
    postfix : str, optional
        The postfix/suffix to search for (e.g., '.py', '-dev')

    Returns
    =======
    program_names : set
        A set of program names that match all specified criteria
    """
    program_names = set()

    # Get all directories in PATH
    path_dirs = os.environ.get('PATH', '').split(os.pathsep)

    for path_dir in path_dirs:
        if not path_dir or not os.path.isdir(path_dir):
            continue

        try:
            # List all files in the directory
            for item in os.listdir(path_dir):
                item_path = os.path.join(path_dir, item)

                # Check if it's a file and executable
                if os.path.isfile(item_path) and os.access(item_path, os.X_OK):
                    # Check all specified criteria
                    matches = True

                    if prefix is not None:
                        if not item.lower().startswith(prefix.lower()):
                            matches = False

                    if contains is not None and matches:
                        if contains.lower() not in item.lower():
                            matches = False

                    if postfix is not None and matches:
                        if not item.lower().endswith(postfix.lower()):
                            matches = False

                    if matches:
                        program_names.add(item)

        except (PermissionError, OSError):
            # Skip directories we can't read
            continue

    return program_names



def is_program_exists(program, dont_raise=False):
    IsExe = lambda p: os.path.isfile(p) and os.access(p, os.X_OK)

    fpath, fname = os.path.split(program)

    if fpath:
        if IsExe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = os.path.expanduser(path).strip('"')
            exe_file = os.path.join(path, program)
            if IsExe(exe_file):
                return exe_file

    if dont_raise:
        return False

    raise ConfigError("An anvi'o function needs '%s' to be installed on your system, but it doesn't seem to appear "
                       "in your path :/ If you are certain you have it on your system (for instance you can run it "
                       "by typing '%s' in your terminal window), you may want to send a detailed bug report. Sorry!"\
                        % (program, program))



def check_h5py_module():
    """To make sure we do have the h5py module.

       The reason this function is here is becasue we removed h5py from anvi'o dependencies,
       but some migration scripts may still need it if the user has very old databases. In
       those cases the user must install it manually."""

    try:
        import h5py
        h5py.__version__
    except:
        raise ConfigError("There is an issue but it is easy to resolve and everything is fine! "
                          "To continue, please first install the newest version of the Python module `h5py` "
                          "by running `pip install h5py` in your anvi'o environment. "
                          "The reason why the standard anvi'o package does not include this module is both "
                          "complicated and really unimportant. Re-running the migration after `h5py` is installed "
                          "will make things go smoothly.")



def is_all_npm_packages_installed(run=Run()):
    """A function to test whether all npm packages are installed in the interactive directory.

    This check is for ensuring that necessary npm packages are installed in the
    anvio/data/interactive directory.
    """

    # find the root directory of anvi'o module
    anvio_module_path = os.path.dirname(os.path.abspath(anvio.__file__))
    interactive_dir_path = os.path.join(anvio_module_path, 'data', 'interactive')

    if not os.path.exists(interactive_dir_path):
        raise ConfigError("The interactive directory does not exist in the anvi'o module. "
                          "Please ensure the directory is present.")

    # Check if Node.js is installed
    if shutil.which("node") is None:
        run.warning("It seems your installation is missing Node.js, a recent requirement of anvi'o "
                    "environments. Please run the following command in your terminal, and you should "
                    "be good to go:", header="⚠️  YOUR ATTENTION PLEASE ⚠️", overwrite_verbose=True, lc='yellow')
        run.info_single("      conda install -c conda-forge nodejs", level=0, overwrite_verbose=True, nl_before=1)
        raise ConfigError("Node.js is not installed. Please install it using conda and try again.")

    # Check if node_modules exists and is not empty
    node_modules_path = os.path.join(interactive_dir_path, 'node_modules')

    if not os.path.exists(node_modules_path) or not os.listdir(node_modules_path):
        run.warning("Anvi'o recently changed its use of external libraries for interactive interfaces"
                    "from git submodules to npm packages. Your current setup does not seem to have the "
                    "necessary files in place, so the purpose of this warning is to help you match your "
                    "setup to most up-to-date anvi'o code. If you run the commands below in your terminal, "
                    "you will most likely be fine :) But if things don't work out, please reach out to us "
                    "on GitHub or Discord since this is a new feature and some hiccups may occur.",
                    header="⚠️  YOUR ATTENTION PLEASE ⚠️", overwrite_verbose=True,
                    lc='yellow')
        run.info_single(f"1) cd {interactive_dir_path}", level=0, overwrite_verbose=True)
        run.info_single("2) npm install", level=0, overwrite_verbose=True)
        run.info_single("3) cd -", level=0, overwrite_verbose=True)

        raise ConfigError("Some npm packages seem to be missing in your interactive directory. "
                          "Please run 'npm install' in the interactive directory and try again.")
    else:
        return True

