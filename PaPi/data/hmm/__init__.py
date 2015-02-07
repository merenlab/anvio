import os
import glob
import textwrap

import PaPi.utils as u
import PaPi.terminal as terminal
from PaPi.constants import allowed_chars

run = terminal.Run()

dir_path = os.path.dirname(os.path.abspath(__file__))
sources = {}

allowed_chars = allowed_chars.replace('.', '').replace('-', '')
PROPER = lambda w: not len([c for c in w if c not in allowed_chars]) \
                   and len(w) >= 3 \
                   and w[0] not in '_0123456789'


for source in [s for s in glob.glob(os.path.join(dir_path, '*')) if s.find('.py') < 0 and s.find('.txt') < 0]:
    if not PROPER(os.path.basename(source)):
        raise u.ConfigError, "One of the search database directory ('%s') contains characters in its name\
                              PaPi does not like. Directory names should be at least three characters long\
                              and must not contain any characters but ASCII letters, digits and\
                              underscore" % os.path.basename(source)

    for f in ['reference.txt', 'kind.txt', 'genes.txt', 'genes.hmm.gz']:
        if not os.path.exists(os.path.join(source, f)):
            raise u.ConfigError, "Each search database directory must contain following files:\
                                  'kind.txt', 'reference.txt', 'genes.txt', and 'genes.hmm.gz'. %s does not seem\
                                  to be a proper source." % os.path.basename(source)

    kind = open(os.path.join(source, 'kind.txt')).readlines()[0].strip()
    if not PROPER(kind):
        raise u.ConfigError, "'kind.txt' defines the kind of search this database offers. This file must contain a single\
                              word that is at least three characters long, and must not contain any characters but\
                              ASCII letters, digits, and underscore. Here are some nice examples: 'singlecopy',\
                              or 'pathogenicity', or 'noras_selection'. But yours is '%s'." % (kind)

    genes = u.get_TAB_delimited_file_as_dictionary(os.path.join(source, 'genes.txt'), column_names = ['gene', 'accession'])

    sources[os.path.basename(source)] = {'ref': os.path.join(source, 'reference.txt'),
                                         'kind': kind,
                                         'genes': genes.keys(),
                                         'model': os.path.join(source, 'genes.hmm.gz')}


# lets make sure stuff we need is installed on this system
# these are necessary for single-copy gene analysis.
missing_programs = []
for p in ['prodigal', 'hmmscan']:
    try:
        u.is_program_exists(p)
    except u.ConfigError:
        missing_programs.append(p)

if len(sources):
    if len(missing_programs):
        run.info('WARNING', '', header = True)
        print textwrap.fill(u.remove_spaces('PaPi found one or more databases to perform single-copy gene analysis\
                                             that may be useful for downstream analyses. However, this process require\
                                             certain programs to be present in your system. Here is the command(s) PaPi\
                                             tried to access and failed: %s. Please see the documentation for system\
                                             requirements.' % ', '.join(missing_programs)), 80) + '\n'
        sources = {}
    else:
        run.info('HMM profiling data',
                 'Loaded from %d source%s; %s' % (len(sources),
                                                  's' if len(sources) > 1 else '',
                                                  ', '.join(['%s (%d genes)' % (s, len(sources[s]['genes'])) for s in sources])))
