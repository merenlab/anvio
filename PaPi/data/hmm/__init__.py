import os
import glob
import textwrap

import PaPi.utils as u
import PaPi.terminal as terminal
from PaPi.constants import allowed_chars

run = terminal.Run()

dir_path = os.path.dirname(os.path.abspath(__file__))
sources = {}

for source in [s for s in glob.glob(os.path.join(dir_path, '*')) if s.find('.py') < 0 and s.find('.txt') < 0]:
    if len([c for c in os.path.basename(source) if c not in allowed_chars]):
        raise u.ConfigError, "One of the directories for single-copy gene analysis ('%s') contains\
                              characters PaPi does not like. Becathese directory names will also be\
                              used as variable names, they must be composed of ASCII letters,\
                              digits, '_' and '.' alone." % os.path.basename(source)

    for f in ['reference.txt', 'genes.txt', 'genes.hmm.gz']:
        if not os.path.exists(os.path.join(source, f)):
            raise u.ConfigError, "Each directory with single-copy gene analysis must contain a\
                                  'reference.txt', 'genes.txt', and 'genes.hmm'. %s does not seem\
                                  to be a proper source." % os.path.basename(source)

    genes = u.get_TAB_delimited_file_as_dictionary(os.path.join(source, 'genes.txt'), column_names = ['pfam_id', 'gene'], indexing_field = 1)

    sources[os.path.basename(source)] = {'ref': os.path.join(source, 'reference.txt'),
                                         'genes': genes.keys(),
                                         'hmm': os.path.join(source, 'genes.hmm.gz')}


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
        run.info('Bacterial single-copy genes database',
                 'Loaded from %d source%s; %s' % (len(sources),
                                                  's' if len(sources) > 1 else '',
                                                  ', '.join(['%s (%d genes)' % (s, len(sources[s]['genes'])) for s in sources])))
