import os
import sys
import glob

init_py_path = os.path.normpath(os.path.dirname(os.path.abspath(__file__))) + '/anvio/__init__.py'
version_string = [l.strip() for l in open(init_py_path).readlines() if l.strip().startswith('anvio_version')][0]
anvio_version = version_string.split('=')[1].strip().strip("'").strip('"')

requirements = [req.strip() for req in open('requirements.txt', 'rU').readlines() if not req.startswith('#')]

try:
    if sys.version_info.major != 3:
        sys.stderr.write("Your active Python major version ('%d') is not compatible with what anvi'o expects :/ We recently switched to Python 3.\n" % sys.version_info.major)
        sys.exit(-1)
except Exception:
    sys.stderr.write("(anvi'o failed to learn about your Python version, but it will pretend as if nothing happened)\n\n")

try:
    import numpy
except ImportError:
    print("You need to have numpy installed on your system to run setup.py. Sorry!")
    sys.exit()

try:
    from Cython.Distutils import build_ext
except ImportError:
    print("You need to have Cython installed on your system to run setup.py. Sorry!")
    sys.exit()

from setuptools import setup, find_packages, Extension

if os.environ.get('USER','') == 'vagrant':
    del os.link

os.chdir(os.path.normpath(os.path.join(os.path.abspath(__file__), os.pardir)))

with open(os.path.join(os.path.dirname(__file__), 'README.md')) as readme:
    README = readme.read()

include_dirs_for_concoct = [numpy.get_include(), '/opt/local/include/', '/usr/local/include/']
library_dirs_for_concoct = ['/opt/local/lib/', '/usr/local/lib/']

setup(
    name = "anvio",
    version = anvio_version,

    scripts = [script for script in glob.glob('bin/*') + glob.glob('sandbox/*') if not script.endswith('-OBSOLETE')],
    include_package_data = True,

    packages = find_packages(),

    install_requires = requirements,

    cmdclass = {'build_ext': build_ext},
    ext_modules = [
                    Extension("anvio.vbgmm", sources=["./anvio/extensions/concoct/vbgmm.pyx", "./anvio/extensions/concoct/c_vbgmm_fit.c"],
                                libraries =['gsl', 'gslcblas'], include_dirs=include_dirs_for_concoct, library_dirs=library_dirs_for_concoct),
                  ],

    author = "anvi'o Authors",
    author_email = "a.murat.eren@gmail.com",
    description = "An interactive analysis and visualization platform for 'omics data. See https://merenlab.org/projects/anvio for more information",
    license = "GPLv3+",
    keywords = "metagenomics metatranscriptomics microbiology shotgun genomics MBL pipeline sequencing bam visualization SNP SNV",
    url = "https://merenlab.org/projects/anvio/",
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Environment :: Web Environment',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Natural Language :: English',
        'Operating System :: MacOS',
        'Operating System :: POSIX',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: JavaScript',
        'Programming Language :: C',
        'Topic :: Scientific/Engineering',
    ],
)
