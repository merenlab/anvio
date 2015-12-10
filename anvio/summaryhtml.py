# coding: utf-8
"""Creates an HTML output to act as a front-end for the static summary directory."""

import os
import shutil

import anvio
import anvio.terminal as terminal

from anvio.errors import ConfigError

# get django imported
try:
    from django.conf import settings
    absolute = os.path.join(os.path.dirname(os.path.realpath(__file__)))
    template_dir = os.path.join(absolute, 'data/static/template')
    html_content_dir = os.path.join(absolute, 'data/static/content')
    settings.configure(DEBUG=True, TEMPLATE_DEBUG=True, DEFAULT_CHARSET='utf-8', TEMPLATE_DIRS = (template_dir,))

    try:
        import django
        django.setup()
    except:
        pass

    from django.template.loader import render_to_string
    from django.template.defaultfilters import register
except ImportError:
    raise ConfigError, 'You need to have Django module (http://djangoproject.com) installed on your system to generate HTML output.'


__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2015, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


pp = terminal.pretty_print
run = terminal.Run()
progress = terminal.Progress()


class SummaryHTMLOutput:
    def __init__(self, summary_dict = {}, r = run, p = progress):
        self.run = r
        self.progress = p
        self.summary_dict = summary_dict


    def generate(self, quick = False):
        self.progress.new('Copying static files')
        self.copy_files()

        self.progress.new('Rendering')
        index_html = self.render(quick)

        self.run.info('HTML Output', index_html)

        return index_html


    def copy_files(self):
        self.progress.update('...')
        shutil.copytree(html_content_dir, os.path.join(self.summary_dict['meta']['output_directory'], '.html'))
        self.progress.end()


    def render(self, quick = False):
        self.progress.update('Processing the template ...')

        if quick:
            rendered = render_to_string('index-mini.tmpl', self.summary_dict)
        else:
            rendered = render_to_string('index.tmpl', self.summary_dict)

        index_html = os.path.join(self.summary_dict['meta']['output_directory'], 'index.html')
        self.progress.update('Writing the index file ...')
        open(index_html, 'w').write(rendered.encode('utf-8'))

        self.progress.end()
        return index_html


@register.filter(name='lookup')
def lookup(d, index):
    if index in d:
        return d[index]
    return ''

@register.filter(name='humanize')
def humanize(s):
    return s.replace('_', ' ')

@register.filter(name='sumvals')
def sumvals(d):
    return sum(d.values())

@register.filter(name='humanize_n')
def humanize_n(n):
    if isinstance(n, str):
        return n

    for unit in ['', ' Kb', ' Mb']:
        if abs(n) < 1000.0:
            return "%3.2f%s" % (n, unit)
        n /= 1000.0
    return "%.2f%s" % (n, 'Gb')

@register.filter(name='pretty')
def pretty(n):
    if not n:
        return 'None'
    try:
        return pp(int(n))
    except ValueError:
        return n
