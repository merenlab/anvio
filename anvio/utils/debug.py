import numpy as np
import pandas as pd

from anvio.errors import ConfigError
from anvio.terminal import TimeCode

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

