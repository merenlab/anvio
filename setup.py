import os
import glob
from setuptools import setup, find_packages

os.chdir(os.path.normpath(os.path.join(os.path.abspath(__file__), os.pardir)))

with open(os.path.join(os.path.dirname(__file__), 'README.md')) as readme:
    README = readme.read()


setup(
    name = "papi",
    version = open('VERSION').read().strip(),

    scripts = [script for script in glob.glob('bin/*') if not script.endswith('-OBSOLETE')],
    include_package_data = True,

    packages = find_packages(),

    install_requires = ['bottle>=0.12.7', 'pysam==0.7.7', 'hcluster>=0.2.0', 'ete2>=2.2', 'scipy>=0.14.0', 'scikit-learn>=0.15'],

    author = "PaPi Authors",
    author_email = "a.murat.eren@gmail.com",
    description = "Post-assembly Environmental (Meta)genomics Pipeline",
    longer_description=README,
    license = "GPLv3+",
    keywords = "metagenomics microbiology shotgun genomics MBL pipeline sequencing bam",
    url = "https://meren.github.io/projects/papi/",
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
