import os
import glob

import anvio.utils as u
import anvio.terminal as terminal

run = terminal.Run()

dir_path = os.path.dirname(os.path.abspath(__file__))
sources = u.get_HMM_sources_dictionary([s for s in glob.glob(os.path.join(dir_path, '*')) if not s.endswith('.py') and s.find('.txt') < 0 and not os.path.basename(s).startswith('__')])
scg_domain_to_source = dict([(sources[s]['domain'], s) for s in sources if sources[s]['kind'] == 'singlecopy'])

if len(sources):
   run.info('HMM profiles',
            '%d source%s been loaded: %s' % (len(sources),
                                             's have' if len(sources) > 1 else ' has',
                                             ', '.join(['%s (%d genes, domain: %s)' % (s, len(sources[s]['genes']), sources[s]['domain']) for s in sources])))
