import os
import glob

import anvio.utils as u
import anvio.terminal as terminal

run = terminal.Run()

dir_path = os.path.dirname(os.path.abspath(__file__))
sources = u.get_HMM_sources_dictionary([s for s in glob.glob(os.path.join(dir_path, '*')) if s.find('.py') < 0 and s.find('.txt') < 0 and not os.path.basename(s).startswith('__')])

# lets make sure stuff we need is installed on this system
# to perform an HMM search
missing_programs = u.get_missing_programs_for_hmm_analysis()

if len(sources):
    if len(missing_programs):
        run.warning("Anvi'o found one or more databases to perform single-copy gene analysis\
                     that may be useful for downstream analyses. However, this process require\
                     certain programs to be present in your system. Here is the command(s) anvio\
                     tried to access and failed: %s. Please see the documentation for system\
                     requirements." % ', '.join(missing_programs))
        sources = {}
    else:
        run.info('HMM profiles',
                 '%d source%s been loaded: %s' % (len(sources),
                                                  's have' if len(sources) > 1 else ' has',
                                                  ', '.join(['%s (%d genes, domain: %s)' % (s, len(sources[s]['genes']), sources[s]['domain']) for s in sources])))
