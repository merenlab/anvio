# coding: utf-8
# pylint: disable=line-too-long
"""Creates an HTML output to act as a front-end for the static summary directory."""

import os
import json
import zlib
import base64
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
    local_settings = {
        'DEBUG': True,
        'TEMPLATE_DEBUG': True,
        'DEFAULT_CHARSET': 'utf-8'
    }

    try:
        # Django 1.10+ will use TEMPLATES variable instead of TEMPLATE_DIRS.
        from django.template.backends.django import DjangoTemplates
        local_settings.update({
            'TEMPLATES': [
                {
                    'BACKEND': 'django.template.backends.django.DjangoTemplates',
                    'DIRS': (template_dir,),
                    'APP_DIRS': False,
                }
            ]
        })
    except ImportError:
        local_settings.update({'TEMPLATE_DIRS': (template_dir,)})

    try:
        import django
    except:
        pass

    from django.template.loader import render_to_string
    from django.template.defaultfilters import register
except ImportError:
    raise ConfigError('You need to have Django module (http://djangoproject.com) installed on your system to generate HTML output.')


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
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
    def __init__(self, summary_dict={}, r=run, p=progress):
        settings.configure(**local_settings)
        django.setup()

        self.run = r
        self.progress = p
        self.summary_dict = summary_dict

        self.summary_type = self.summary_dict['meta']['summary_type']

        if self.summary_type not in ['profile', 'pan', 'saav', 'vignette']:
            raise ConfigError("Unknown summary type '%s'" % self.summary_type)


    def generate(self, quick=False):
        self.progress.new('Generating the output')
        self.progress.update('Copying static files')
        self.copy_files()

        self.progress.update('Rendering')
        rendered = self.render(quick)

        self.progress.update('Writing')
        index_html = os.path.join(self.summary_dict['meta']['output_directory'], 'index.html')
        open(index_html, 'wb').write(rendered.encode('utf-8'))

        self.progress.end()

        self.run.info('HTML Output', index_html, nl_before=1, nl_after=1, mc='green')

        return index_html


    def copy_files(self):
        self.progress.update('...')
        shutil.copytree(html_content_dir, os.path.join(self.summary_dict['meta']['output_directory'], '.html'))


    def render(self, quick=False):
        if self.summary_type == 'saav':
                rendered = render_to_string('saavs-index.tmpl', self.summary_dict)
        elif self.summary_type == 'pan':
            rendered = render_to_string('pan-index.tmpl', self.summary_dict)
        elif self.summary_type == 'profile':
            if quick:
                rendered = render_to_string('profile-index-mini.tmpl', self.summary_dict)
            else:
                rendered = render_to_string('profile-index.tmpl', self.summary_dict)
        elif self.summary_type == 'vignette':
            rendered = render_to_string('vignette.tmpl', self.summary_dict)
        else:
            raise ConfigError("You cray...")

        return rendered


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

@register.filter(name='humanize_f')
def humanize_f(n):
    if isinstance(n, str):
        return n

    return "%.2f" % n

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

@register.filter(name='convert_to_json')
def convert_to_json(obj):
    return json.dumps(obj)

@register.filter(name='base64_encode')
def base64_encode(data):
    return base64.b64encode(data)

@register.filter(name='zlib_encode')
def zlib_encode(data):
    return zlib.compress(data.encode("utf-8"))

