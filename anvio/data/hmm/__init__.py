import glob
import os

import anvio.terminal as terminal
import anvio.utils as u
from anvio.terminal import pluralize as P

run = terminal.Run()

dir_path = os.path.dirname(os.path.abspath(__file__))
sources = u.get_HMM_sources_dictionary([s for s in glob.glob(os.path.join(dir_path, '*')) if not s.endswith('.py') and s.find('.txt') < 0 and not os.path.basename(s).startswith('__')])
scg_domain_to_source = dict([(sources[s]['domain'], s) for s in sources if sources[s]['kind'] == 'singlecopy'])

if len(sources):
    run.info('HMM profiles', f"{len(sources)} {P('has', len(sources), alt='have')} been loaded:")
    for s in sources:
        run.info_single(f"{s} ({len(sources[s]['genes'])} genes, domain: {sources[s]['domain']})", level=2, mc='cyan')
