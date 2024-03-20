# -*- coding: utf-8
# pylint: disable=line-too-long
"""A library to process anvi'o authors file"""

import os

import anvio
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.errors import ConfigError


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


class AnvioAuthors:
    def __init__(self, authors_yaml_file_path=os.path.join(os.path.dirname(anvio.__file__), 'data/misc/PEOPLE/DEVELOPERS.yaml'), skip_init=False, r=terminal.Run(), p=terminal.Progress()):
        self.run = r
        self.progress = p

        self.authors_yaml_file_path = authors_yaml_file_path
        self.author_avatars_directory = os.path.join(os.path.dirname(authors_yaml_file_path), 'AVATARS')

        self.essential_author_info_keys = ['github', 'name', 'email']

        self.authors = {}

        self.initialized = False

        if not skip_init:
            self.init_authors()


    def init_authors(self):
        if self.initialized:
            # so we are asked to re-initialize
            self.authors = {}

        """Initializes the `self.authors` dictionary."""

        if not os.path.exists(self.authors_yaml_file_path):
            raise ConfigError(f"The YAML file for authors is not found at {self.authors_yaml_file_path} :(")

        if not os.path.exists(self.author_avatars_directory):
            raise ConfigError(f"The directory that was supposed to contain all the author avatars is not "
                              f"seem to be at {self.author_avatars_directory}. Not good.")


        authors_yaml = utils.get_yaml_as_dict(self.authors_yaml_file_path)

        authors_missing_avatars_in_yaml = set([])
        authors_missing_avatars_on_disk = set([])

        for entry in authors_yaml:
            missing_keys = [info_key for info_key in self.essential_author_info_keys if info_key not in entry]
            if len(missing_keys):
                for key in entry:
                    self.run.info(key, entry[key])

                raise ConfigError(f"Author entries in the YAML file must at least describe the following "
                                  f"keys: {', '.join(self.essential_author_info_keys)}. The entry shown "
                                  f"above misses some of them: {', '.join(missing_keys)}.")

            # Usernames should be case insensitive
            entry['github'] = entry['github'].lowercase()   
            if entry['github'] in self.authors:
                raise ConfigError(f"The GitHub username '{entry['github']}' is used for multiple authors in "
                                  f"the YAML file :/")

            if 'avatar' not in entry:
                authors_missing_avatars_in_yaml.add(entry['github'])
                entry['avatar'] = 'no-avatar.png'

            entry['avatar'] = os.path.join(self.author_avatars_directory, entry['avatar'])

            if not os.path.exists(entry['avatar']):
                entry['avatar'] = os.path.join(self.author_avatars_directory, 'no-avatar.png')
                authors_missing_avatars_on_disk.add(entry['github'])

            self.authors[entry['github']] = entry

        if len(authors_missing_avatars_on_disk):
            self.run.warning(f"The following authors have avatars defined in the YAML file, but there are no "
                             f"corresponding image files in the avatars directory: {', '.join(authors_missing_avatars_on_disk)}. "
                             f"So anvi'o asssigned a generic profile image for these people. To make sure their avatars are "
                             f"associated with their user names in the YAML file, you need to put the matching file under this "
                             f"directory: `{self.author_avatars_directory}`).", header="MISSING AVATARS ON DISK")

        self.run.warning(f"Some authors in the YAML file didn't have any avatars defined for them, so anvi'o will "
                         f"assign a generic profile image for these people: {', '.join(authors_missing_avatars_in_yaml)}. "
                         f"You could make these people much prettier by adding an `avatar` key for them, and put the "
                         f"corresponding file under the '{self.author_avatars_directory}' directory. JUST SAYING.",
                         header="MISSING AVATARS IN THE YAML FILE")
