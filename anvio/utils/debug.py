import os
import shutil
import numpy as np
import pandas as pd

import anvio

from anvio.errors import ConfigError
from anvio.terminal import Run, TimeCode

# for full output
pd.options.display.max_columns=100
pd.options.display.max_rows=100

def compare_times(calls, as_matrix=False, iterations_per_call=1):
    """Compare times between function calls

    Parameters
    ==========
    calls : list of tuples
        Each element should be a (name, function, args, kwargs) tuples. If there are no args or
        kwargs, the element should look like (name, function, [], {})

    as_matrix : bool, False
        If True, results are output as a pandas matrix, where each element is a time difference between
        calls. Otherwise, a dictionary is returned

    iterations_per_call : int, 1
        How many times should each function call be ran? Time will be averaged

    Returns
    =======
    times : pd.DataFrame or dict
        If as_matrix, pd.DataFrame is returned, where times[i, j] is how much faster i is than j.
        Otherwise, dictionary of {name: time} is returned
    """

    call_times = np.zeros((len(calls), iterations_per_call))
    names, *_ = zip(*calls)
    for i, call in enumerate(calls):
        name, function, args, kwargs = call

        for j in range(iterations_per_call):
            try:
                with TimeCode(quiet=True) as t:
                    function(*args, **kwargs)
            except:
                raise ConfigError("compare_times :: function call with name '%s' failed." % name)

            call_times[i, j] = t.time.total_seconds()

    averaged_call_times = np.mean(call_times, axis=1)

    if not as_matrix:
        return dict(zip(names, averaged_call_times))

    matrix = []
    for i, _time in enumerate(call_times):
        row = []

        for j, _time in enumerate(call_times):
            row.append(averaged_call_times[j] - averaged_call_times[i] if i > j else 'NA')

        matrix.append(row)

    return pd.DataFrame(matrix, columns=names, index=names)


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

