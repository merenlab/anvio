"""Shared method validation for anvi'o's supported ANI workflows."""

from anvio.errors import ConfigError


SUPPORTED_ANI_METHODS = ('ANIb', 'ANIm')


def validate_ani_method(method):
    """Reject ANI methods retired by the Python 3.13 port with guidance."""
    if method == 'ANIblastall':
        raise ConfigError("The ANI method 'ANIblastall' is retired in the Python 3.13 port because it depends on legacy BLAST. "
                          "Recalculate with ANIb if fragment-based BLASTN ANI fits your analysis; ANIb is a distinct method "
                          "and its values are not guaranteed to match ANIblastall.")
    if method == 'TETRA':
        raise ConfigError("The ANI method 'TETRA' is retired in the Python 3.13 port and has no equivalent in the supported "
                          "ANIb/ANIm methods. Choose one of those methods and recalculate if alignment-based ANI is suitable.")
    if method not in SUPPORTED_ANI_METHODS:
        raise ConfigError("Unsupported ANI method '%s'. Choose ANIb or ANIm." % method)
