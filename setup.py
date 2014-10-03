import os
import glob
from setuptools import setup, find_packages

os.chdir(os.path.normpath(os.path.join(os.path.abspath(__file__), os.pardir)))

with open(os.path.join(os.path.dirname(__file__), 'README.md')) as readme:
    README = readme.read()


setup(
    name = "PaPi",
    version = "0.1",
    packages = find_packages(),
    scripts = [script for script in glob.glob('bin/*') if not script.endswith('-OBSOLETE')],

    install_requires = ['bottle>=0.12.7', 'pysam>=0.7.5', 'hcluster>=0.2.0', 'ete2>=2.2', 'scipy>=0.14.0'],


    author = "A. Murat Eren",
    author_email = "a.murat.eren@gmail.com",
    description = "Post-assembly environmental genomics pipeline",
    longer_description=README,
    license = "GPL",
    keywords = "metagenomics microbiology shotgun genomics MBL pipeline sequencing bam",
    url = "http://meren.org/research/papi/",   # project home page, if any
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Environment :: Web Environment',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Natural Language :: English',
        'Operating System :: MacOS',
        'Operating System :: POSIX',
        'Programming Language :: JavaScript',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering',
    ],
)
