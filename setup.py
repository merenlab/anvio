import glob
from setuptools import setup

# find all anvi'o programs
scripts = glob.glob("bin/anvi-*") + glob.glob("sandbox/anvi-script-*")

setup(
    scripts=scripts,
)
