"""Validation helpers for user-supplied metabolism tables."""

import numpy as np
import pandas as pd

from anvio.errors import ConfigError


# These are the default missing-value markers accepted by pandas' text reader.
# Keep recognition scoped to numeric fields so identifiers such as the literal
# sample or gene ID "NA" are not changed.
PANDAS_DEFAULT_MISSING_MARKERS = frozenset({
    '', '#N/A', '#N/A N/A', '#NA', '-1.#IND', '-1.#QNAN', '-NaN', '-nan',
    '1.#IND', '1.#QNAN', '<NA>', 'N/A', 'NA', 'NaN', 'None', 'NULL',
    'n/a', 'nan', 'null',
})


def parse_numeric_column(dataframe, column, source_path, identifier_column=None):
    """Parse a numeric column, preserving standard missing markers and rejecting bad values."""
    raw_values = dataframe[column]
    stripped_values = raw_values.astype('string').str.strip()
    missing = raw_values.isna() | stripped_values.isin(PANDAS_DEFAULT_MISSING_MARKERS)
    numeric_values = pd.to_numeric(stripped_values.mask(missing, pd.NA), errors='coerce')

    missing_values = missing.to_numpy(dtype=bool)
    malformed = ((~missing) & numeric_values.isna()).to_numpy(dtype=bool)
    numeric_as_float = numeric_values.to_numpy(dtype='float64', na_value=np.nan)
    non_finite = (~missing_values) & numeric_values.notna().to_numpy(dtype=bool) & ~np.isfinite(numeric_as_float)
    invalid_positions = np.flatnonzero(malformed | non_finite)

    if len(invalid_positions):
        bad_rows = []
        for row_position in invalid_positions[:5]:
            row_number = row_position + 2
            row_details = [f"row {row_number}"]
            if identifier_column and identifier_column in dataframe.columns:
                row_details.append(f"{identifier_column}={dataframe.iloc[row_position][identifier_column]!r}")
            bad_value = raw_values.iloc[row_position]
            if hasattr(bad_value, 'item'):
                bad_value = bad_value.item()
            row_details.append(f"value={bad_value!r}")
            bad_rows.append(' (' + ', '.join(row_details) + ')')

        raise ConfigError(f"The '{column}' column in '{source_path}' must contain finite numbers or recognized "
                          f"missing-value markers. Invalid value(s): {', '.join(bad_rows)}.")

    return numeric_values.astype('float64')
