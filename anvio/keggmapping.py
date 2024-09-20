#!/usr/bin/env python
# -*- coding: utf-8
"""Make KEGG pathway maps incorporating data sourced from anvi'o databases."""

import os
import re
import fitz
import math
import shutil
import functools
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from argparse import Namespace
from itertools import combinations
from typing import Dict, Iterable, List, Literal, Tuple, Union

import anvio.kegg as kegg
import anvio.kgml as kgml
import anvio.utils as utils
import anvio.dbinfo as dbinfo
import anvio.terminal as terminal
import anvio.reactionnetwork as rn
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.genomestorage import GenomeStorage
from anvio.dbops import ContigsDatabase, PanSuperclass
from anvio import FORCE_OVERWRITE, QUIET, __version__ as VERSION


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2024, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = VERSION
__maintainer__ = "Samuel Miller"
__email__ = "samuelmiller10@gmail.com"
__status__ = "Development"


# The colors of qualitative and repeating colormaps are sampled in order, whereas other colormaps,
# including sequential colormaps, are sampled evenly.
qualitative_colormaps: List[str] = [
    'Pastel1',
    'Pastel2',
    'Paired',
    'Accent',
    'Dark2',
    'Set1',
    'Set2',
    'Set3',
    'tab10',
    'tab20',
    'tab20b',
    'tab20c'
]
repeating_colormaps: List[str] = [
    'flag',
    'prism'
]

class Mapper:
    """
    Make KEGG pathway maps incorporating data sourced from anvi'o databases.

    Attributes
    ==========
    kegg_context : anvio.kegg.KeggContext
        This contains anvi'o KEGG database attributes, such as filepaths.

    available_pathway_numbers : List[str]
        ID numbers of all pathways set up with PNG and KGML files in the KEGG data directory.

    pathway_names : Dict[str, str]
        The names of all KEGG pathways, including those without files in the KEGG data directory.
        Keys are pathway ID numbers and values are pathway names.

    rn_constructor : anvio.reactionnetwork.Constructor
        Used for loading reaction networks from anvi'o databases.

    xml_ops : anvio.kgml.XMLOps
        Used for loading KGML files as pathway objects.

    overwrite_output : bool
        If True, methods in this class overwrite existing output files.

    name_files : bool
        Include the pathway name along with the number in output map file names.

    run : anvio.terminal.Run
        This object prints run information to the terminal.

    progress : anvio.terminal.Progress
        This object prints transient progress information to the terminal.
    """
    def __init__(
        self,
        kegg_dir: str = None,
        overwrite_output: bool = FORCE_OVERWRITE,
        name_files: bool = False,
        run: terminal.Run = terminal.Run(),
        progress: terminal.Progress = terminal.Progress(),
        quiet: bool = QUIET
    ) -> None:
        """
        Parameters
        ==========
        kegg_dir : str, None
            Directory containing an anvi'o KEGG database. The default argument of None expects the
            KEGG database to be set up in the default directory used by the program
            anvi-setup-kegg-data.

        overwrite_output : bool, anvio.FORCE_OVERWRITE
            If True, methods in this class overwrite existing output files.

        name_files : bool, False
            Include the pathway name along with the number in output map file names.

        run : anvio.terminal.Run, anvio.terminal.Run()
            This object prints run information to the terminal.

        progress : anvio.terminal.Progress, anvio.terminal.Progress()
            This object prints transient progress information to the terminal.

        quiet : bool, anvio.QUIET
            If True, run and progress information is not printed to the terminal.
        """
        args = Namespace()
        args.kegg_data_dir = kegg_dir
        self.kegg_context = kegg.KeggContext(args)

        if not os.path.exists(self.kegg_context.kegg_map_image_kgml_file):
            raise ConfigError(
                "One of the key files required by KEGG pathway maps is missing in your active "
                "anvi'o installation. If your KEGG data are not stored at the default KEGG data "
                "location, include that path using the 'kegg_dir' argument. Otherwise, please "
                "consider using the program `anvi-setup-kegg-data` to set up the latest KEGG data "
                "that includes the necessary files for KEGG pathway maps."
            )

        available_pathway_numbers: List[str] = []
        for row in pd.read_csv(
            self.kegg_context.kegg_map_image_kgml_file, sep='\t', index_col=0
        ).itertuples():
            if row.KO + row.EC + row.RN == 0:
                continue
            available_pathway_numbers.append(row.Index[-5:])
        self.available_pathway_numbers = available_pathway_numbers

        pathway_names: Dict[str, str] = {}
        for pathway_number, pathway_name in pd.read_csv(
            self.kegg_context.kegg_pathway_list_file, sep='\t', header=None
        ).itertuples(index=False):
            pathway_names[pathway_number[3:]] = pathway_name
        self.pathway_names = pathway_names

        self.rn_constructor = rn.Constructor(kegg_dir=self.kegg_context.kegg_data_dir)

        self.xml_ops = kgml.XMLOps()
        self.drawer = kgml.Drawer(kegg_dir=self.kegg_context.kegg_data_dir)

        self.name_files = name_files
        self.overwrite_output = overwrite_output
        self.run = run
        self.progress = progress
        self.quiet = self._quiet = quiet

    def map_contigs_database_kos(
        self,
        contigs_db: str,
        output_dir: str,
        pathway_numbers: Iterable[str] = None,
        color_hexcode: str = '#2ca02c',
        draw_maps_lacking_kos: bool = False
    ) -> Dict[str, bool]:
        """
        Draw pathway maps, highlighting KOs present in the contigs database.

        Parameters
        ==========
        contigs_db : str
            File path to a contigs database containing KO annotations.

        output_dir : str
            Path to the output directory in which pathway map PDF files are drawn. The directory is
            created if it does not exist.

        pathway_numbers : Iterable[str], None
            Regex patterns to match the ID numbers of the drawn pathway maps. The default of None
            draws all available pathway maps in the KEGG data directory.

        color_hexcode : str, '#2ca02c'
            This is the color, by default green, for reactions containing contigs database KOs.
            Alternatively to a color hex code, the string, 'original', can be provided to use the
            original color scheme of the reference map. In global and overview maps, KOs are
            represented in reaction lines. The foreground color of lines is set. In standard maps,
            KOs are represented in boxes, the background color of which is set, or lines.

        draw_maps_lacking_kos : bool, False
            If False, by default, only draw maps containing any of the KOs in the contigs database.
            If True, draw maps regardless, meaning that nothing may be colored.

        Returns
        =======
        Dict[str, bool]
            Keys are pathway numbers. Values are True if the map was drawn, False if the map was not
            drawn because it did not contain any of the select KOs and 'draw_maps_lacking_kos' was
            False.
        """
        # Retrieve the IDs of all KO annotations in the contigs database.
        self.progress.new("Loading KO data from the contigs database")
        self.progress.update("...")

        self._check_contigs_db(contigs_db)
        self._check_contigs_db_ko_annotation(contigs_db)

        cdb = ContigsDatabase(contigs_db)
        ko_ids = cdb.db.get_single_column_from_table(
            'gene_functions',
            'accession',
            unique=True,
            where_clause='source = "KOfam"'
        )
        self.progress.end()

        drawn = self._map_kos_fixed_colors(
            ko_ids,
            output_dir,
            pathway_numbers=pathway_numbers,
            color_hexcode=color_hexcode,
            draw_maps_lacking_kos=draw_maps_lacking_kos
        )
        count = sum(drawn.values()) if drawn else 0
        self.run.info("Number of maps drawn", count)

        return drawn

    def map_genomes_storage_genome_kos(
        self,
        genomes_storage_db: str,
        genome_name: str,
        output_dir: str,
        pathway_numbers: Iterable[str] = None,
        color_hexcode: str = '#2ca02c',
        draw_maps_lacking_kos: bool = False
    ) -> Dict[str, bool]:
        """
        Draw pathway maps, highlighting KOs present in the genome.

        Parameters
        ==========
        genomes_storage_db : str
            File path to a genomes storage database containing KO annotations.

        genome_name : str
            Name of a genome in the genomes storage.

        output_dir : str
            Path to the output directory in which pathway map PDF files are drawn. The directory is
            created if it does not exist.

        pathway_numbers : Iterable[str], None
            Regex patterns to match the ID numbers of the drawn pathway maps. The default of None
            draws all available pathway maps in the KEGG data directory.

        color_hexcode : str, '#2ca02c'
            This is the color, by default green, for reactions containing KOs in the genome.
            Alternatively to a color hex code, the string, 'original', can be provided to use the
            original color scheme of the reference map. In global and overview maps, KOs are
            represented in reaction lines. The foreground color of lines is set. In standard maps,
            KOs are represented in boxes, the background color of which is set, or lines.

        draw_maps_lacking_kos : bool, False
            If False, by default, only draw maps containing any of the KOs in the genome. If True,
            draw maps regardless, meaning that nothing may be colored.

        Returns
        =======
        Dict[str, bool]
            Keys are pathway numbers. Values are True if the map was drawn, False if the map was not
            drawn because it did not contain any of the select KOs and 'draw_maps_lacking_kos' was
            False.
        """
        # Retrieve the IDs of all KO annotations for the genome.
        self.progress.new("Loading KO data from the genome")
        self.progress.update("...")

        self._check_genomes_storage_db(genomes_storage_db)
        self._check_genomes_storage_ko_annotation(genomes_storage_db)

        gsdb = GenomeStorage(
            genomes_storage_db,
            genome_names_to_focus=[genome_name],
            function_annotation_sources=['KOfam'],
            run=terminal.Run(verbose=False),
            progress=terminal.Progress(verbose=False)
        )
        ko_ids = gsdb.db.get_single_column_from_table(
            'gene_function_calls',
            'accession',
            unique=True,
            where_clause=f'genome_name = "{genome_name}" AND source = "KOfam"'
        )
        self.progress.end()

        drawn = self._map_kos_fixed_colors(
            ko_ids,
            output_dir,
            pathway_numbers=pathway_numbers,
            color_hexcode=color_hexcode,
            draw_maps_lacking_kos=draw_maps_lacking_kos
        )
        count = sum(drawn.values()) if drawn else 0
        self.run.info("Number of maps drawn", count)

        return drawn

    def map_contigs_databases_kos(
        self,
        contigs_dbs: Iterable[str],
        output_dir: str,
        groups_txt: str = None,
        group_threshold: float = None,
        pathway_numbers: Iterable[str] = None,
        draw_contigs_db_files: Union[Iterable[str], bool] = False,
        draw_grid: Union[Iterable[str], bool] = False,
        colormap: Union[bool, str, mcolors.Colormap] = True,
        colormap_limits: Tuple[float, float] = None,
        colormap_scheme: Literal['by_count', 'by_membership'] = None,
        reverse_overlay: bool = False,
        color_hexcode: str = '#2ca02c',
        colorbar: bool = True,
        draw_maps_lacking_kos: bool = False
    ) -> Dict[Literal['unified', 'individual', 'grid'], Dict]:
        """
        Draw pathway maps, highlighting KOs across contigs databases (representing, for example,
        genomes or metagenomes) or groups of databases (representing, for example, taxonomic
        groups of genomes or geographical groups of metagenomes).

        A reaction on a map can correspond to one or more KOs, and a KO can annotate one or more
        sequences in a contigs database. In global and overview maps, reaction lines are colored.
        In standard maps, reaction boxes or lines are colored.

        Parameters
        ==========
        contigs_dbs : Iterable[str]
            File paths to contigs databases containing KO annotations. Databases should have
            different project names, by which they are uniquely identified.

        output_dir : str
            Path to the output directory in which pathway map PDF files are drawn. The directory is
            created if it does not exist.

        groups_txt : str, None
            A tab-delimited text file specifying which group each contigs database belongs to. The
            first column must have the header, 'contigs_db', and the second column must have the
            header, 'group'. Items in the first column must be the file paths or project names of
            contigs databases. The second column contains group names, which are recommended to be
            single words without fancy characters, such as 'HIGH_TEMPERATURE' or 'LOW_FITNESS'
            rather than 'my group #1' or 'IS-THIS-OK?'. Each contigs database can only be associated
            with a single group. The 'group_threshold' argument must also be used for the groups to
            take effect, assigning colors based on group membership and drawing individual files
            ('draw_contigs_db_files') and map grids ('draw_grid') for groups rather than individual
            databases.

        group_threshold : float, None
            The proportion of contigs databases in a group containing data of interest for the group
            to be represented in terms of presence/absence in a map feature. Here is a concrete
            example. Say each contigs database represents a genome, and the 'groups_txt' argument,
            which must be used with this argument, groups these genomes by their species, 'A', 'B',
            and 'C'. You wish to understand the distribution of metabolic capabilities across the 3
            species from KO annotations of genes. Reaction colors are assigned based on the groups
            rather than individual genomic contigs databases containing the reaction. Thresholds
            between 0 and 1 can be set to define group membership: a threshold of 0.0 would mean
            that ANY genome in the group can contain the KO for the KO to be considered present in
            the group; a threshold of 0.75 means at least 75% of the genomes in the group must
            contain the KO for it to be present; a threshold of 1.0 means that ALL genomes in the
            group must contain the KO for it to be present. In our example, set the threshold to
            0.5. Reaction J on a map corresponds to KO X, and Reaction K on a map corresponds to KO
            Y. 90% of species A genomes, 50% of species B genomes, and 10% of species C genomes
            contain KO X, so Reaction J would be colored to indicate that it is represented in
            species A and B. 0% of species A genomes, 15% of species B genomes, and 40% of species C
            genomes contain KO Y, so Reaction K would not be colored.

        pathway_numbers : Iterable[str], None
            Regex patterns to match the ID numbers of the drawn pathway maps. The default of None
            draws all available pathway maps in the KEGG data directory.

        draw_contigs_db_files : Union[Iterable[str], bool], False
            Draw pathway maps for each contigs database if not False. If True, draw maps for all of
            the contigs databases. Alternatively, the project names of a subset of contigs databases
            can be provided.

            If groups are defined by the 'groups_txt' argument, then maps are instead drawn for
            individual groups. A subset of group names can be provided to draw maps for select
            groups.

        draw_grid : Union[Iterable[str], bool], False
            If not False, draw a grid for each pathway map showing both the unified map of input
            contigs databases and a map for each contigs database, facilitating identification of
            the contigs databases containing reactions highlighted in the unified map. If True,
            include all of the contigs databases in the grid. Alternatively, the project names of a
            subset of contigs databases can be provided.

            If groups are defined by the 'groups_txt' argument, then the grid instead shows the
            unified map of groups and maps for individual groups. A subset of group names can be
            provided to select maps in the grid.

        colormap : Union[bool, str, matplotlib.colors.Colormap], True
            Reactions are dynamically colored to reflect the contigs databases (or groups of
            databases) involving the reaction, unless the argument value is False. False overrides
            dynamic coloring via a colormap with the argument provided to 'color_hexcode', so that
            reactions represented by KOs in contigs databases are assigned predetermined colors.

            The default argument value of True automatically assigns a colormap given the colormap
            scheme (see the 'colormap_scheme' argument). The scheme, 'by_count', uses by default the
            sequential colormap, 'plasma_r', which spans yellow (fewer databases or groups) to
            blue-violet (more databases or groups). This accentuates reactions that are shared
            rather than unshared across databases. In contrast, a colormap spanning dark to light,
            such as 'plasma', is better for drawing attention to unshared reactions. The scheme,
            'by_membership', uses by default the qualitative colormap, 'tab10'; it contains distinct
            colors suitable for clearly differentiating the different databases or groups containing
            reactions.

            The name of a Matplotlib Colormap or a Colormap object itself can also be provided to be
            used in lieu of the default. See the following webpage for named colormaps:
            https://matplotlib.org/stable/users/explain/colors/colormaps.html#classes-of-colormaps

        colormap_limits : Tuple[float, float], None
            Limit the fraction of the colormap used in dynamically selecting colors. The first value
            is the lower cutoff and the second value is the upper cutoff, e.g., (0.2, 0.9) limits
            color selection to 70% of the colormap, trimming the bottom 20% and top 10%. By default,
            for the colormap scheme, 'by_count', the colormap is 'plasma_r', and the limits are set
            to (0.1, 0.9). By default, for the scheme, 'by_membership', the colormap is
            qualititative ('tab10'), and limits are set to (0.0, 1.0).

        colormap_scheme : Literal['by_count', 'by_membership'], None
            There are two ways of dynamically coloring reactions by inclusion in contigs databases
            or groups of databases: by count, or explicitly by database/group or combination of
            databases/groups. Given the default argument value of None, with 4 or more
            databases/groups, reactions are colored by count, and with 2 or 3, explicitly by
            membership. In coloring by count, the colormap should be sequential, such that the color
            of a reaction changes 'smoothly' with the count. In contrast, coloring by database/group
            means reaction color is determined by membership in a database/group or combination of
            databases/groups, so a qualitative colormap can be suitable for clearly differentiating
            each.

        reverse_overlay : bool, False
            By default, with False, reactions in more contigs databases or groups of databases are
            drawn on top of those in fewer databases/groups. With True, the opposite applies;
            especially in global maps with a non-default colormap spanning dark to light, this
            accentuates unshared rather than shared parts of a pathway.

        color_hexcode : str, '#2ca02c'
            This is the color, by default green, for reactions containing contigs database KOs.
            Alternatively to a color hex code, the string, 'original', can be provided to use the
            original color scheme of the reference map. The 'colormap' argument must be False for
            this argument to be used, overriding dynamic coloring based on database/group membership
            with static coloring based on presence/absence in any database.

        colorbar : bool, True
            If True and coloring by database/group membership, save a colorbar legend to the file,
            'colorbar.pdf', in the output directory.

        draw_maps_lacking_kos : bool, False
            If False, by default, only draw maps containing any of the select KOs. If True, draw
            maps regardless, meaning that nothing may be colored.

        Returns
        =======
        Dict[Literal['unified', 'individual', 'grid'], Dict]
            Keys in the outer dictionary are different types of files that can be drawn. 'unified'
            maps show data from all contigs databases or groups of databases. 'individual' maps show
            data from individual databases or groups. 'grid' images show both unified and individual
            maps. 'unified' and 'grid' values are Dict[str, bool], where keys are pathway numbers,
            and values are True if the map was drawn, False if the map was not drawn because it did
            not contain any of the select KOs and 'draw_maps_lacking_kos' was False. 'individual'
            values are Dict[str, Dict[str, bool]], where keys in the outer dictionary are contigs
            database project names or group names, keys in the inner dictionary are pathway numbers,
            and values in the inner dictionary are True if the map was drawn, False if the map was
            not drawn because it did not contain any of the select KOs and 'draw_maps_lacking_kos'
            was False.
        """
        # This method is similar to map_pan_database_kos, and almost identical after KOs are loaded.

        # Load groups.
        if (
            (groups_txt is None and group_threshold is not None) or
            (groups_txt is not None and group_threshold is None)
        ):
            raise ConfigError(
                "To group contigs databases, arguments to both 'groups_txt' and 'group_threshold' "
                "must be provided."
            )

        if groups_txt is None:
            source_group_dict = None
            group_source_dict = None
            membership = contigs_dbs
        else:
            source_group_dict, group_source_dict = utils.get_groups_txt_file_as_dict(
                groups_txt, run=self.run, progress=self.progress
            )
            if not 0 <= group_threshold <= 1:
                raise ConfigError(
                    f"'group_threshold' must be a number between 0 and 1, not {group_threshold}"
                )
            membership = list(group_source_dict)

        # Set the colormap scheme.
        if colormap is False:
            scheme = 'static'
            if groups_txt is not None:
                self.run.warning(
                    "Groups were provided by 'groups_txt', but these will be ignored, since "
                    "'colormap' was set to False, and dynamic coloring based on KO membership in "
                    "groups will be overridden by static coloring based on KO presence/absence in "
                    "any contigs database."
                )
        else:
            if colormap_scheme is None:
                if len(membership) < 4:
                    scheme = 'by_membership'
                else:
                    scheme = 'by_count'
            elif colormap_scheme == 'by_count':
                scheme = 'by_count'
            elif colormap_scheme == 'by_membership':
                scheme = 'by_membership'
            else:
                raise AssertionError

        # Set the colormap.
        if colormap is True:
            if scheme == 'by_count':
                cmap = plt.colormaps['plasma_r']
                if colormap_limits is None:
                    colormap_limits = (0.1, 0.9)
            elif scheme == 'by_membership':
                cmap = plt.colormaps['tab10']
                if colormap_limits is None:
                    colormap_limits = (0.0, 1.0)
            else:
                raise AssertionError
        elif colormap is False:
            cmap = None
        elif isinstance(colormap, str):
            cmap = plt.colormaps[colormap]
            if colormap_limits is None:
                colormap_limits = (0.0, 1.0)
        elif isinstance(colormap, mcolors.Colormap):
            cmap = colormap
            if colormap_limits is None:
                colormap_limits = (0.0, 1.0)
        else:
            raise AssertionError

        # Set how the colormap is sampled.
        if cmap is None:
            sampling = None
        else:
            if cmap.name in qualitative_colormaps + repeating_colormaps:
                sampling = 'in_order'
            else:
                sampling = 'even'

        # Trim the colormap.
        if cmap is not None and colormap_limits is not None and colormap_limits != (0.0, 1.0):
            assert 0.0 <= colormap_limits[0] <= colormap_limits[1] <= 1.0
            cmap = mcolors.LinearSegmentedColormap.from_list(
                f'trunc({cmap.name},{colormap_limits[0]:.2f},{colormap_limits[1]:.2f})',
                cmap(range(
                    int(colormap_limits[0] * cmap.N), math.ceil(colormap_limits[1] * cmap.N)
                ))
            )

        self.progress.new("Loading KO data from contigs databases")
        self.progress.update("...")

        self._check_contigs_dbs(contigs_dbs)
        self._check_contigs_dbs_ko_annotation(contigs_dbs)

        # Load contigs database metadata.
        project_names: Dict[str, str] = {}
        for contigs_db in contigs_dbs:
            contigs_db_info = dbinfo.ContigsDBInfo(contigs_db)
            self_table = contigs_db_info.get_self_table()

            annotation_sources = self_table['gene_function_sources']
            assert annotation_sources is not None and 'KOfam' in annotation_sources.split(',')

            project_name = self_table['project_name']
            assert project_name not in project_names
            project_names[project_name] = contigs_db

        # Check that groups include all contigs databases.
        if groups_txt is not None:
            source_paths = [os.path.abspath(source) for source in source_group_dict]
            ungrouped_contigs_dbs: List[str] = []
            for contigs_db in contigs_dbs:
                if os.path.abspath(contigs_db) not in source_paths:
                    ungrouped_contigs_dbs.append(contigs_db)
            message = ', '.join([f"'{contigs_db}'" for contigs_db in ungrouped_contigs_dbs])
            raise ConfigError(
                "The following 'contigs_dbs' were not found in the groups provided by "
                f"'groups_txt': {message}"
            )

        # Find which contigs databases contain each KO.
        ko_dbs: Dict[str, List[str]] = {}
        for project_name, contigs_db in project_names.items():
            cdb = ContigsDatabase(contigs_db)
            for ko_id in cdb.db.get_single_column_from_table(
                'gene_functions',
                'accession',
                unique=True,
                where_clause='source = "KOfam"'
            ):
                try:
                    ko_dbs[ko_id].append(project_name)
                except KeyError:
                    ko_dbs[ko_id] = [project_name]
        self.progress.end()

        # Find the numeric IDs of the maps to draw.
        pathway_numbers = self._find_maps(output_dir, 'kos', patterns=pathway_numbers)

        filesnpaths.gen_output_directory(output_dir, progress=self.progress, run=self.run)

        drawn: Dict[Literal['unified', 'individual', 'grid'], Dict] = {
            'unified': {},
            'individual': {},
            'grid': {}
        }

        self.progress.new("Drawing 'unified' map incorporating data from all contigs databases")
        exceeds_colors: Tuple[int, int] = None
        if scheme == 'static':
            # Draw unified maps of all contigs databases with a static reaction color.
            for pathway_number in pathway_numbers:
                self.progress.update(pathway_number)
                if color_hexcode == 'original':
                    drawn['unified'][pathway_number] = self._draw_map_kos_original_color(
                        pathway_number,
                        ko_dbs,
                        output_dir,
                        draw_map_lacking_kos=draw_maps_lacking_kos
                    )
                else:
                    drawn['unified'][pathway_number] = self._draw_map_kos_single_color(
                        pathway_number,
                        ko_dbs,
                        color_hexcode,
                        output_dir,
                        draw_map_lacking_kos=draw_maps_lacking_kos
                    )
        else:
            # Draw unified maps with dynamic coloring by number of contigs databases.
            assert cmap is not None
            color_priority: Dict[str, float] = {}
            if scheme == 'by_count':
                # Sample the colormap for colors representing each possible number of contigs
                # databases. Lower color values correspond to smaller numbers of databases.
                if sampling == 'in_order':
                    if len(contigs_dbs) == 1:
                        sample_points = range(1, 2)
                    else:
                        sample_points = range(len(contigs_dbs))
                elif sampling == 'even':
                    if len(contigs_dbs) == 1:
                        sample_points = np.linspace(1, 1, 1)
                    else:
                        sample_points = np.linspace(0, 1, len(contigs_dbs))
                else:
                    raise AssertionError

                if len(contigs_dbs) > cmap.N:
                    exceeds_colors = (cmap.N, len(contigs_dbs))

                for sample_point in sample_points:
                    if reverse_overlay:
                        color_priority[mcolors.rgb2hex(cmap(sample_point))] = 1 - sample_point
                    else:
                        color_priority[mcolors.rgb2hex(cmap(sample_point))] = sample_point
                db_combos = None
            elif scheme == 'by_membership':
                # Sample the colormap for colors representing the different contigs databases and
                # their combinations. Lower color values correspond to smaller numbers of databases.
                db_combos = []
                for db_count in range(1, len(contigs_dbs) + 1):
                    db_combos += list(combinations(project_names, db_count))

                if sampling == 'in_order':
                    sample_points = range(len(db_combos))
                elif sampling == 'even':
                    sample_points = np.linspace(0, 1, len(db_combos))
                else:
                    raise AssertionError

                if len(db_combos) > cmap.N:
                    exceeds_colors = (cmap.N, len(db_combos))

                for sample_point in sample_points:
                    if reverse_overlay:
                        color_priority[
                            mcolors.rgb2hex(cmap(sample_point))
                        ] = 1 - sample_point / cmap.N
                    else:
                        color_priority[
                            mcolors.rgb2hex(cmap(sample_point))
                        ] = (sample_point + 1) / cmap.N
            else:
                raise AssertionError

            if colorbar:
                # Draw a colorbar in a separate file.
                _draw_colorbar = self._draw_colorbar
                if scheme == 'by_count':
                    _draw_colorbar = functools.partial(
                        _draw_colorbar,
                        color_labels=range(1, len(contigs_dbs) + 1),
                        label='database count'
                    )
                elif scheme == 'by_membership':
                    _draw_colorbar = functools.partial(
                        _draw_colorbar,
                        color_labels=[', '.join(db_combo) for db_combo in db_combos],
                        label='databases'
                    )
                _draw_colorbar(
                    color_priority, os.path.join(output_dir, 'colorbar.pdf')
                )

            for pathway_number in pathway_numbers:
                self.progress.update(pathway_number)
                drawn['unified'][pathway_number] = self._draw_map_kos_membership(
                    pathway_number,
                    ko_dbs,
                    color_priority,
                    output_dir,
                    cmap,
                    source_combos=db_combos,
                    draw_map_lacking_kos=draw_maps_lacking_kos
                )
        self.progress.end()

        if exceeds_colors:
            self.run.warning(
                f"There were fewer distinct colors available in the colormap ({exceeds_colors[0]}) "
                f"than were needed ({exceeds_colors[1]}), so some colors were repeated in use."
            )

        if draw_contigs_db_files is False and draw_grid is False:
            count = sum(drawn['unified'].values()) if drawn['unified'] else 0
            self.run.info("Number of maps drawn", count)
            return

        # Determine the individual database maps to draw.
        if draw_contigs_db_files == True:
            draw_files_project_names = list(project_names)
        elif draw_contigs_db_files == False:
            draw_files_project_names = []
        else:
            for project_name in draw_contigs_db_files:
                assert project_name in project_names
            draw_files_project_names = draw_contigs_db_files
        seen = set()
        draw_files_project_names = [
            project_name for project_name in list(draw_files_project_names)
            if not (project_name in seen or seen.add(project_name))
        ]

        # Determine the map grids to draw.
        if draw_grid == True:
            draw_grid_project_names = list(project_names)
        elif draw_grid == False:
            draw_grid_project_names = []
        else:
            for project_name in draw_grid:
                assert project_name in project_names
            draw_grid_project_names = draw_grid
        seen = set()
        draw_grid_project_names = [
            project_name for project_name in list(draw_grid_project_names)
            if not (project_name in seen or seen.add(project_name))
        ]

        seen = set()
        draw_project_names = [
            project_name for project_name in draw_files_project_names + draw_grid_project_names
            if not (project_name in seen or seen.add(project_name))
        ]

        # Draw individual database maps needed as final outputs or for grids.
        for project_name in draw_project_names:
            self.progress.new(f"Drawing maps for contigs database '{project_name}'")
            self.progress.update("...")
            progress = self.progress
            self.progress = terminal.Progress(verbose=False)
            run = self.run
            self.run = terminal.Run(verbose=False)
            drawn['individual'][project_name] = self.map_contigs_database_kos(
                project_names[project_name],
                os.path.join(output_dir, project_name),
                pathway_numbers=pathway_numbers,
                color_hexcode=color_hexcode,
                draw_maps_lacking_kos=draw_maps_lacking_kos
            )
            self.progress = progress
            self.run = run
            self.progress.end()

        if draw_grid == False:
            count = sum(drawn['unified'].values()) if drawn['unified'] else 0
            self.run.info(
                "Number of 'unified' maps drawn incorporating data from all contigs databases",
                count
            )
            if not drawn['individual']:
                count = 0
            else:
                count = sum([sum(d.values()) if d else 0 for d in drawn['individual'].values()])
            self.run.info("Number of maps drawn for individual contigs databases", count)
            return

        self.progress.new("Drawing map grid")
        self.progress.update("...")

        # Draw empty maps needed to fill in grids.
        paths_to_remove: List[str] = []
        if not draw_maps_lacking_kos:
            # Make a new dictionary with outer keys being pathway numbers, inner dictionaries
            # indicating which maps were drawn per contigs database.
            drawn_pathway_number: Dict[str, Dict[str, bool]] = {}
            for project_name, drawn_project_name in drawn['individual'].items():
                for pathway_number, drawn_map in drawn_project_name.items():
                    try:
                        drawn_pathway_number[pathway_number][project_name] = drawn_map
                    except KeyError:
                        drawn_pathway_number[pathway_number] = {project_name: drawn_map}

            # Draw empty maps as needed, for pathways with some but not all maps drawn.
            progress = self.progress
            self.progress = terminal.Progress(verbose=False)
            run = self.run
            self.run = terminal.Run(verbose=False)
            for pathway_number, drawn_project_name in drawn_pathway_number.items():
                if set(drawn_project_name.values()) != set([True, False]):
                    continue
                for project_name, drawn_map in drawn_project_name.items():
                    if drawn_map:
                        continue
                    self.map_contigs_database_kos(
                        project_names[project_name],
                        os.path.join(output_dir, project_name),
                        pathway_numbers=[pathway_number],
                        color_hexcode=color_hexcode,
                        draw_maps_lacking_kos=True
                    )
                    if self.name_files:
                        pathway_name = '_' + self._get_filename_pathway_name(pathway_number)
                    else:
                        pathway_name = ''
                    paths_to_remove.append(os.path.join(
                        output_dir, project_name, f'kos_{pathway_number}{pathway_name}.pdf'
                    ))
            self.progress = progress
            self.run = run

        # Draw map grids.
        grid_dir = os.path.join(output_dir, 'grid')
        filesnpaths.gen_output_directory(grid_dir, progress=self.progress, run=self.run)
        for pathway_number in pathway_numbers:
            self.progress.update(pathway_number)
            if self.name_files:
                pathway_name = '_' + self._get_filename_pathway_name(pathway_number)
            else:
                pathway_name = ''
            unified_map_path = os.path.join(output_dir, f'kos_{pathway_number}{pathway_name}.pdf')
            if not os.path.exists(unified_map_path):
                continue
            in_paths = [unified_map_path]
            labels = ['all']

            pdf_doc = fitz.open(in_paths[0])
            page = pdf_doc.load_page(0)
            input_aspect_ratio = page.rect.width / page.rect.height
            landscape = True if input_aspect_ratio > 1 else False

            for project_name in draw_grid_project_names:
                if self.name_files:
                    pathway_name = '_' + self._get_filename_pathway_name(pathway_number)
                else:
                    pathway_name = ''
                individual_map_path = os.path.join(
                    output_dir, project_name, f'kos_{pathway_number}{pathway_name}.pdf'
                )
                if not os.path.exists(individual_map_path):
                    break
                in_paths.append(individual_map_path)
                labels.append(project_name)
            else:
                if self.name_files:
                    pathway_name = '_' + self._get_filename_pathway_name(pathway_number)
                else:
                    pathway_name = ''
                out_path = os.path.join(grid_dir, f'kos_{pathway_number}{pathway_name}.pdf')
                self.make_grid(in_paths, out_path, labels=labels, landscape=landscape)
                drawn['grid'][pathway_number] = True
        self.progress.end()

        # Remove individual database maps that were only needed for map grids.
        for path in paths_to_remove:
            os.remove(path)
        for project_name in set(draw_project_names).difference(set(draw_files_project_names)):
            shutil.rmtree(os.path.join(output_dir, project_name))
            drawn['individual'].pop(project_name)

        count = sum(drawn['unified'].values()) if drawn['unified'] else 0
        self.run.info(
            "Number of 'unified' maps drawn incorporating data from all contigs databases",
            count
        )
        if draw_contigs_db_files:
            if not drawn['individual']:
                count = 0
            else:
                count = sum([sum(d.values()) if d else 0 for d in drawn['individual'].values()])
            self.run.info("Number of maps drawn for individual contigs databases", count)
        count = sum(drawn['grid'].values()) if drawn['grid'] else 0
        self.run.info("Number of map grids drawn", count)

        return drawn

    def map_pan_database_kos(
        self,
        pan_db: str,
        genomes_storage_db: str,
        output_dir: str,
        pathway_numbers: Iterable[str] = None,
        draw_genome_files: Union[Iterable[str], bool] = False,
        draw_grid: Union[Iterable[str], bool] = False,
        colormap: Union[str, mcolors.Colormap, None] = 'plasma_r',
        colormap_limits: Tuple[float, float] = None,
        reverse_overlay: bool = False,
        color_hexcode: str = '#2ca02c',
        colorbar: bool = True,
        draw_maps_lacking_kos: bool = False,
        consensus_threshold: float = None,
        discard_ties: bool = None
    ) -> Dict[Literal['unified', 'individual', 'grid'], Dict]:
        """
        Draw pathway maps, highlighting consensus KOs from the pan database.

        A reaction on a map can correspond to one or more KOs, and a KO can annotate one or more
        gene clusters. In global and overview maps, reaction lines are colored. In standard maps,
        reaction boxes or lines are colored.

        Parameters
        ==========
        pan_db : str
            File path to a pangenomic database. If a reaction network was stored in the database,
            then consensus KOs are determined using the consensus_threshold and discard_ties
            parameters stored as database metadata unless explicitly given here as arguments. These
            parameters are only stored in the database when a reaction network is stored.

        genomes_storage_db : str
            Path to the genomes storage database associated with the pan database. This contains
            KO annotations.

        output_dir : str
            Path to the output directory in which pathway map PDF files are drawn. The directory is
            created if it does not exist.

        pathway_numbers : Iterable[str], None
            Regex patterns to match the ID numbers of the drawn pathway maps. The default of None
            draws all available pathway maps in the KEGG data directory.

        draw_genome_files : Union[Iterable[str], bool], False
            Draw pathway maps for genomes of the pangenome if not False. If True, draw maps for all
            of the genomes. Alternatively, the names of a subset of genomes can be provided.

        draw_grid : Union[Iterable[str], bool], False
            If not False, draw a grid for each pathway map showing both the pangenomic map and a map
            for each genome of the pangenome, facilitating identification of the genomes containing
            reactions highlighted in the pangenomic map. If True, include all of the genomes in the
            grid. Alternatively, the names of a subset of genomes can be provided.

        colormap : Union[str, matplotlib.colors.Colormap, None], 'plasma_r'
            Reactions are dynamically colored to reflect the number of genomes involving the
            reaction, unless the argument value is None. None overrides dynamic coloring via a
            colormap using the argument provided to 'color_hexcode', so that reactions in the
            pangenome are assigned predetermined colors.

            Here is how a reaction is assigned a genome count. A reaction element in a map can
            contain one or more KOs. Find corresponding consensus KOs from the anvi'o pangenomic
            database. Each consensus KO is assigned to one or more gene clusters. Counted genomes
            have one or more genes in gene clusters with these consensus KOs.

            This argument can take either be the name of a built-in matplotlib colormap or a
            Colormap object itself. The default sequential colormap, 'plasma_r', spans yellow (fewer
            genomes) to blue-violet (more genomes). This accentuates reactions that are shared
            rather than unshared across genomes. A colormap spanning dark (fewer genomes) to light
            (more genomes), such as 'plasma', is better for drawing attention to unshared reactions.

            See the following webpage for named colormaps:
            https://matplotlib.org/stable/users/explain/colors/colormaps.html#classes-of-colormaps

        colormap_limits : Tuple[float, float], (0.0, 1.0)
            Limit the fraction of the colormap used in dynamically selecting colors. The first value
            is the lower cutoff and the second value is the upper cutoff, e.g., (0.2, 0.9) limits
            color selection to 70% of the colormap, trimming the bottom 20% and top 10%. The default
            limits with the default colormap scheme, 'plasma_r', are set to (0.1, 0.9).

        reverse_overlay : bool, False
            By default, with False, reactions in more genomes are drawn on top of those in fewer
            genomes. With True, the opposite applies; especially in global maps with a non-default
            colormap spanning dark to light, this accentuates unshared rather than shared parts of a
            pathway.

        color_hexcode : str, '#2ca02c'
            This is the color, by default green, for reactions containing consensus KOs from the pan
            database. Alternatively to a color hex code, the string, 'original', can be provided to
            use the original color scheme of the reference map. The 'colormap' argument must be
            False for this argument to be used, overriding dynamic coloring based on quantitative
            data with static coloring based on presence/absence in the pangenome.

        colorbar : bool, True
            If True and coloring by number of genomes, save a colorbar legend to the file,
            'colorbar.pdf', in the output directory.

        draw_maps_lacking_kos : bool, False
            If False, by default, only draw maps containing any of the select KOs. If True, draw
            maps regardless, meaning that nothing may be colored.

        consensus_threshold : float, None
            With a value of None, if a reaction network was stored in the pan database, then the
            consensus_threshold metavalue that was also stored in the database is used to find
            consensus KOs. If a reaction network was not stored, then with a value of None, the KO
            annotation most frequent in a gene cluster is assigned to the cluster itself. If a
            numerical value is provided (must be on [0, 1]), at least this proportion of genes in
            the cluster must have the most frequent annotation for the cluster to be annotated.

        discard_ties : bool, None
            With a value of None, if a reaction network was stored in the pan database, then the
            discard_ties metavalue that was also stored in the database is used to find consensus
            KOs. If a reaction network was not stored, then with a value of None, discard_ties
            assumes a value of False. A value of True means that if multiple KO annotations are most
            frequent among genes in a cluster, then a consensus KO is not assigned to the cluster
            itself, whereas a value of False would cause one of the most frequent KOs to be
            arbitrarily chosen.

        Returns
        =======
        Dict[Literal['unified', 'individual', 'grid'], Dict]
            Keys in the outer dictionary are different types of files that can be drawn. 'unified'
            maps show data from all genomes. 'individual' maps show data from individual genomes.
            'grid' images show both unified and individual maps. 'unified' and 'grid' values are
            Dict[str, bool], where keys are pathway numbers, and values are True if the map was
            drawn, False if the map was not drawn because it did not contain any of the select KOs
            and 'draw_maps_lacking_kos' was False. 'individual' values are Dict[str, Dict[str,
            bool]], where keys in the outer dictionary are genome names, keys in the inner
            dictionary are pathway numbers, and values in the inner dictionary are True if the map
            was drawn, False if the map was not drawn because it did not contain any of the select
            KOs and 'draw_maps_lacking_kos' was False.
        """
        # This method is similar to map_contigs_databases_kos, and almost identical after KOs are
        # loaded.
        if isinstance(colormap, str):
            assert colormap in mpl.colormaps()

        self._check_pan_db(pan_db)
        self._check_genomes_storage_db(genomes_storage_db)
        self._check_genomes_storage_ko_annotation(genomes_storage_db)

        # Load pan database metadata.
        pan_db_info = dbinfo.PanDBInfo(pan_db)
        self_table = pan_db_info.get_self_table()

        # Parameterize how consensus KOs are found.
        if consensus_threshold is None:
            consensus_threshold = self_table['reaction_network_consensus_threshold']
            if consensus_threshold is not None:
                consensus_threshold = float(consensus_threshold)
                assert 0 <= consensus_threshold <= 1
                self.run.info_single(
                    "No consensus threshold was explicitly specified for consensus KO assignment "
                    f"to gene clusters, but there was a value of '{consensus_threshold}' stored in "
                    "the pan database from reaction network construction, so this was used. (The "
                    "default if this were not the case is 0, or no threshold.)"
                )

        if discard_ties is None:
            discard_ties = self_table['reaction_network_discard_ties']
            if discard_ties is None:
                discard_ties = False
            else:
                discard_ties = bool(int(discard_ties))
                self.run.info_single(
                    "It was not explicitly specified whether to discard ties in consensus KO "
                    f"assignment to gene clusters, but there was a value of '{discard_ties}' "
                    "stored in the pan database from reaction network construction, so this was "
                    "used. (The default if this were not the case is False, or do not discard "
                    "ties.)"
                )

        # Find consensus KOs from the loaded pan database.
        self.progress.new("Loading consensus KO data from pan database")
        self.progress.update("...")
        progress = self.progress
        self.progress = terminal.Progress(verbose=False)
        run = self.run
        self.run = terminal.Run(verbose=False)
        args = Namespace()
        args.pan_db = pan_db
        args.genomes_storage = genomes_storage_db
        args.consensus_threshold = consensus_threshold
        args.discard_ties = discard_ties
        pan_super = PanSuperclass(args, r=self.run, p=self.progress)
        pan_super.init_gene_clusters()
        pan_super.init_gene_clusters_functions()
        pan_super.init_gene_clusters_functions_summary_dict()
        gene_clusters: Dict[str, Dict[str, List[int]]] = pan_super.gene_clusters
        gene_clusters_functions_summary_dict: Dict = pan_super.gene_clusters_functions_summary_dict

        consensus_cluster_ids: List[str] = []
        consensus_ko_ids: List[str] = []
        for cluster_id, gene_cluster_functions_data in gene_clusters_functions_summary_dict.items():
            gene_cluster_ko_data = gene_cluster_functions_data['KOfam']
            if gene_cluster_ko_data == {'function': None, 'accession': None}:
                continue
            consensus_cluster_ids.append(cluster_id)
            consensus_ko_ids.append(gene_cluster_ko_data['accession'])
        self.progress = progress
        self.run = run
        self.progress.end()

        # Find the numeric IDs of the maps to draw.
        pathway_numbers = self._find_maps(output_dir, 'kos', patterns=pathway_numbers)

        filesnpaths.gen_output_directory(output_dir, progress=self.progress, run=self.run)

        genome_names = self_table['external_genome_names'].split(',')

        drawn: Dict[Literal['unified', 'individual', 'grid'], Dict] = {
            'unified': {},
            'individual': {},
            'grid': {}
        }

        self.progress.new("Drawing 'unified' map incorporating data from all genomes")
        exceeds_colors: Tuple[int, int] = None
        if colormap is None:
            # Draw pangenomic maps with a static reaction color.
            for pathway_number in pathway_numbers:
                if color_hexcode == 'original':
                    drawn['unified'][pathway_number] = self._draw_map_kos_original_color(
                        pathway_number,
                        set(consensus_ko_ids),
                        output_dir,
                        draw_map_lacking_kos=draw_maps_lacking_kos
                    )
                else:
                    drawn['unified'][pathway_number] = self._draw_map_kos_single_color(
                        pathway_number,
                        set(consensus_ko_ids),
                        color_hexcode,
                        output_dir,
                        draw_map_lacking_kos=draw_maps_lacking_kos
                    )
            cmap = None
            sampling = None
        else:
            # Draw pangenomic maps with dynamic coloring by number of genomes.
            if isinstance(colormap, str):
                cmap = plt.colormaps[colormap]
                if colormap_limits is None:
                    colormap_limits = (0.1, 0.9)
            else:
                cmap = colormap

            # Set how the colormap is sampled.
            if cmap.name in qualitative_colormaps + repeating_colormaps:
                sampling = 'in_order'
            else:
                sampling = 'even'

            # Trim the colormap.
            if cmap is not None and colormap_limits is not None and colormap_limits != (0.0, 1.0):
                assert 0.0 <= colormap_limits[0] <= colormap_limits[1] <= 1.0
                cmap = mcolors.LinearSegmentedColormap.from_list(
                    f'trunc({cmap.name},{colormap_limits[0]:.2f},{colormap_limits[1]:.2f})',
                    cmap(range(
                        int(colormap_limits[0] * cmap.N), math.ceil(colormap_limits[1] * cmap.N)
                    ))
                )

            # For each consensus KO -- which can annotate more than one gene cluster -- find which
            # genomes contribute genes to clusters represented by the KO.
            ko_genomes: Dict[str, List[str]] = {}
            for cluster_id, ko_id in zip(consensus_cluster_ids, consensus_ko_ids):
                for genome_name, gcids in gene_clusters[cluster_id].items():
                    if not gcids:
                        continue
                    try:
                        ko_genomes[ko_id].append(genome_name)
                    except KeyError:
                        ko_genomes[ko_id] = [genome_name]
            for ko_id, ko_genome_names in ko_genomes.items():
                ko_genomes[ko_id] = list(set(ko_genome_names))

            # Sample the colormap for colors representing each possible number of genomes. Lower
            # color values correspond to smaller numbers of databases.
            if sampling == 'in_order':
                if len(genome_names) == 1:
                    sample_points = range(1, 2)
                else:
                    sample_points = range(len(genome_names))
            elif sampling == 'even':
                if len(genome_names) == 1:
                    sample_points = np.linspace(1, 1, 1)
                else:
                    sample_points = np.linspace(0, 1, len(genome_names))
            else:
                raise AssertionError

            if len(genome_names) > cmap.N:
                exceeds_colors = (cmap.N, len(genome_names))

            color_priority: Dict[str, float] = {}
            for sample_point in sample_points:
                if reverse_overlay:
                    color_priority[mcolors.rgb2hex(cmap(sample_point))] = 1 - sample_point
                else:
                    color_priority[mcolors.rgb2hex(cmap(sample_point))] = sample_point

            if colorbar:
                self._draw_colorbar(
                    color_priority,
                    os.path.join(output_dir, 'colorbar.pdf'),
                    color_labels=range(1, len(genome_names) + 1),
                    label='genomes'
                )
            for pathway_number in pathway_numbers:
                self.progress.update(pathway_number)
                drawn['unified'][pathway_number] = self._draw_map_kos_membership(
                    pathway_number,
                    ko_genomes,
                    color_priority,
                    output_dir,
                    cmap,
                    draw_map_lacking_kos=draw_maps_lacking_kos
                )
        self.progress.end()

        if exceeds_colors:
            self.run.warning(
                f"There were fewer distinct colors available in the colormap ({exceeds_colors[0]}) "
                f"than were needed ({exceeds_colors[1]}), so some colors were repeated in use."
            )

        if draw_genome_files is False and draw_grid is False:
            count = sum(drawn['unified'].values()) if drawn['unified'] else 0
            self.run.info("Number of maps drawn", count)
            return

        # Determine the individual genome maps to draw.
        if draw_genome_files == True:
            draw_files_genome_names = genome_names
        elif draw_genome_files == False:
            draw_files_genome_names = []
        else:
            for genome_name in draw_genome_files:
                assert genome_name in genome_names
            draw_files_genome_names = draw_genome_files
        seen = set()
        draw_files_genome_names = [
            genome_name for genome_name in list(draw_files_genome_names)
            if not (genome_name in seen or seen.add(genome_name))
        ]

        # Determine the map grids to draw.
        if draw_grid == True:
            draw_grid_genome_names = genome_names
        elif draw_grid == False:
            draw_grid_genome_names = []
        else:
            for genome_name in draw_grid:
                assert genome_name in genome_names
            draw_grid_genome_names = draw_grid
        seen = set()
        draw_grid_genome_names = [
            genome_name for genome_name in list(draw_grid_genome_names)
            if not (genome_name in seen or seen.add(genome_name))
        ]

        seen = set()
        draw_genome_names = [
            genome_name for genome_name in draw_files_genome_names + draw_grid_genome_names
            if not (genome_name in seen or seen.add(genome_name))
        ]

        # Draw individual genome maps needed as final outputs or for grids.
        for genome_name in draw_genome_names:
            self.progress.new(f"Drawing maps for genome '{genome_name}'")
            self.progress.update("...")
            progress = self.progress
            self.progress = terminal.Progress(verbose=False)
            run = self.run
            self.run = terminal.Run(verbose=False)
            drawn['individual'][genome_name] = self.map_genomes_storage_genome_kos(
                genomes_storage_db,
                genome_name,
                os.path.join(output_dir, genome_name),
                pathway_numbers=pathway_numbers,
                color_hexcode=color_hexcode,
                draw_maps_lacking_kos=draw_maps_lacking_kos
            )
            self.progress = progress
            self.run = run
            self.progress.end()

        if draw_grid == False:
            count = sum(drawn['unified'].values()) if drawn['unified'] else 0
            self.run.info(
                "Number of 'unified' maps drawn incorporating data from all genomes",
                count
            )
            if not drawn['individual']:
                count = 0
            else:
                count = sum([sum(d.values()) if d else 0 for d in drawn['individual'].values()])
            self.run.info("Number of maps drawn for individual genomes", count)
            return

        self.progress.new("Drawing map grid")
        self.progress.update("...")

        # Draw empty maps needed to fill in grids.
        paths_to_remove: List[str] = []
        if not draw_maps_lacking_kos:
            # Make a new dictionary with outer keys being pathway numbers, inner dictionaries
            # indicating which maps were drawn per genome.
            drawn_pathway_number: Dict[str, Dict[str, bool]] = {}
            for genome_name, drawn_genome_name in drawn['individual'].items():
                for pathway_number, drawn_map in drawn_genome_name.items():
                    try:
                        drawn_pathway_number[pathway_number][genome_name] = drawn_map
                    except KeyError:
                        drawn_pathway_number[pathway_number] = {genome_name: drawn_map}

            # Draw empty maps as needed, for pathways with some but not all maps drawn.
            progress = self.progress
            self.progress = terminal.Progress(verbose=False)
            run = self.run
            self.run = terminal.Run(verbose=False)
            for pathway_number, drawn_genome_name in drawn_pathway_number.items():
                if set(drawn_genome_name.values()) != set([True, False]):
                    continue
                for genome_name, drawn_map in drawn_genome_name.items():
                    if drawn_map:
                        continue
                    self.map_genomes_storage_genome_kos(
                        genomes_storage_db,
                        genome_name,
                        os.path.join(output_dir, genome_name),
                        pathway_numbers=[pathway_number],
                        color_hexcode=color_hexcode,
                        draw_maps_lacking_kos=True
                    )
                    if self.name_files:
                        pathway_name = '_' + self._get_filename_pathway_name(pathway_number)
                    else:
                        pathway_name = ''
                    paths_to_remove.append(os.path.join(
                        output_dir, genome_name, f'kos_{pathway_number}{pathway_name}.pdf'
                    ))
            self.progress = progress
            self.run = run

        # Draw map grids.
        grid_dir = os.path.join(output_dir, 'grid')
        filesnpaths.gen_output_directory(grid_dir, progress=self.progress, run=self.run)
        for pathway_number in pathway_numbers:
            self.progress.update(pathway_number)
            if self.name_files:
                pathway_name = '_' + self._get_filename_pathway_name(pathway_number)
            else:
                pathway_name = ''
            unified_map_path = os.path.join(output_dir, f'kos_{pathway_number}{pathway_name}.pdf')
            if not os.path.exists(unified_map_path):
                continue
            in_paths = [unified_map_path]
            labels = ['pangenome']

            pdf_doc = fitz.open(in_paths[0])
            page = pdf_doc.load_page(0)
            input_aspect_ratio = page.rect.width / page.rect.height
            landscape = True if input_aspect_ratio > 1 else False

            for genome_name in draw_grid_genome_names:
                if self.name_files:
                    pathway_name = '_' + self._get_filename_pathway_name(pathway_number)
                else:
                    pathway_name = ''
                individual_map_path = os.path.join(
                    output_dir, genome_name, f'kos_{pathway_number}{pathway_name}.pdf'
                )
                if not os.path.exists(individual_map_path):
                    break
                in_paths.append(individual_map_path)
                labels.append(genome_name)
            else:
                if self.name_files:
                    pathway_name = '_' + self._get_filename_pathway_name(pathway_number)
                else:
                    pathway_name = ''
                out_path = os.path.join(grid_dir, f'kos_{pathway_number}{pathway_name}.pdf')
                self.make_grid(in_paths, out_path, labels=labels, landscape=landscape)
                drawn['grid'][pathway_number] = True
        self.progress.end()

        # Remove individual genome maps that were only needed for map grids.
        for path in paths_to_remove:
            os.remove(path)
        for genome_name in set(draw_genome_names).difference(set(draw_files_genome_names)):
            shutil.rmtree(os.path.join(output_dir, genome_name))
            drawn['individual'].pop(genome_name)

        count = sum(drawn['unified'].values()) if drawn['unified'] else 0
        self.run.info(
            "Number of 'unified' maps drawn incorporating data from all genomes",
            count
        )
        if draw_genome_files:
            if not drawn['individual']:
                count = 0
            else:
                count = sum([sum(d.values()) if d else 0 for d in drawn['individual'].values()])
            self.run.info("Number of maps drawn for individual genomes", count)
        count = sum(drawn['grid'].values()) if drawn['grid'] else 0
        self.run.info("Number of map grids drawn", count)

        return drawn

    def _map_kos_fixed_colors(
        self,
        ko_ids: Iterable[str],
        output_dir: str,
        pathway_numbers: List[str] = None,
        color_hexcode: str = '#2ca02c',
        draw_maps_lacking_kos: bool = False
    ) -> Dict[str, bool]:
        """
        Draw pathway maps, highlighting reactions containing select KOs in either a single color
        provided by a hex code or the colors originally used in the reference map.

        Parameters
        ==========
        ko_ids : Iterable[str]
            KO IDs to be highlighted in the maps.

        output_dir : str
            Path to the output directory in which pathway map PDF files are drawn. The directory is
            created if it does not exist.

        pathway_numbers : Iterable[str], None
            Regex patterns to match the ID numbers of the drawn pathway maps. The default of None
            draws all available pathway maps in the KEGG data directory.

        color_hexcode : str, '#2ca02c'
            This is the color, by default green, for reactions containing provided KOs.
            Alternatively to a color hex code, the string, 'original', can be provided to use the
            original color scheme of the reference map. In global maps, KOs are represented in
            reaction lines, and in overview maps, KOs are represented in reaction arrows. The
            foreground color of the lines and arrows is set. In standard maps, KOs are represented
            in boxes, the background color of which is set.

        draw_maps_lacking_kos : bool, False
            If False, by default, only draw maps containing any of the select KOs. If True, draw
            maps regardless, meaning that nothing may be colored.

        Returns
        =======
        Dict[str, bool]
            Keys are pathway numbers. Values are True if the map was drawn, False if the map was not
            drawn because it did not contain any of the select KOs and 'draw_maps_lacking_kos' was
            False.
        """
        # Find the numeric IDs of the maps to draw.
        pathway_numbers = self._find_maps(output_dir, 'kos', patterns=pathway_numbers)

        filesnpaths.gen_output_directory(output_dir, progress=self.progress, run=self.run)

        # Draw maps.
        self.progress.new("Drawing map")
        drawn: Dict[str, bool] = {}
        for pathway_number in pathway_numbers:
            self.progress.update(pathway_number)
            if color_hexcode == 'original':
                drawn[pathway_number] = self._draw_map_kos_original_color(
                    pathway_number,
                    ko_ids,
                    output_dir,
                    draw_map_lacking_kos=draw_maps_lacking_kos
                )
            else:
                drawn[pathway_number] = self._draw_map_kos_single_color(
                    pathway_number,
                    ko_ids,
                    color_hexcode,
                    output_dir,
                    draw_map_lacking_kos=draw_maps_lacking_kos
                )
        self.progress.end()

        return drawn

    @staticmethod
    def _check_contigs_db(contigs_db: str) -> None:
        """
        Check the validity of an expected contigs database.

        Parameters
        ==========
        contigs_db : str
            File path to an expected contigs database.
        """
        if not os.path.exists(contigs_db):
            raise ConfigError(
                f"There was no file at the following expected contigs database path: '{contigs_db}'"
            )

        contigs_db_info = dbinfo.ContigsDBInfo(contigs_db, dont_raise=True, expecting='contigs')
        if contigs_db_info is None:
            raise ConfigError(
                "The file at the following expected contigs database path is not a contigs "
                f"database: '{contigs_db}'"
            )

    @staticmethod
    def _check_contigs_db_ko_annotation(contigs_db: str) -> None:
        """
        Check that a contigs database was annotated with KOs.

        Parameters
        ==========
        contigs_db : str
            File path to a contigs database.
        """
        contigs_db_info = dbinfo.ContigsDBInfo(contigs_db, expecting='contigs')
        if 'KOfam' not in contigs_db_info.get_functional_annotation_sources():
            raise ConfigError(
                f"The contigs database, '{contigs_db}', was never annotated with KOs. This can be "
                "rectified by running `anvi-run-kegg-kofams` on the database."
            )

    @staticmethod
    def _check_genomes_storage_db(genomes_storage_db: str) -> None:
        """
        Check the validity of an expected genomes storage database.

        Parameters
        ==========
        genomes_storage_db : str
            File path to an expected genomes storage database.
        """
        if not os.path.exists(genomes_storage_db):
            raise ConfigError(
                "There was no file at the following expected genomes storage database path: "
                f"'{genomes_storage_db}'"
            )

        gsdb_info = dbinfo.GenomeStorageDBInfo(
            genomes_storage_db, dont_raise=True, expecting='genomestorage'
        )
        if gsdb_info is None:
            raise ConfigError(
                "The file at the following expected genomes storage database path is not a genomes "
                f"storage database: '{genomes_storage_db}'"
            )

    @staticmethod
    def _check_genomes_storage_ko_annotation(genomes_storage_db: str) -> None:
        """
        Check that a genomes storage database was annotated with KOs.

        Parameters
        ==========
        genomes_storage_db : str
            File path to a genomes storage database.
        """
        gsdb_info = dbinfo.GenomeStorageDBInfo(genomes_storage_db, expecting='genomestorage')
        if 'KOfam' not in gsdb_info.get_functional_annotation_sources():
            raise ConfigError(
                f"The genomes storage database, '{genomes_storage_db}', was never annotated with "
                "KOs. The genomes storage should be remade with annotated genomes, which can be "
                "rectified by running `anvi-run-kegg-kofams` on the genome databases."
            )

    @staticmethod
    def _check_contigs_dbs(contigs_dbs: Iterable[str]) -> None:
        """
        Check the validity of expected contigs databases.

        Parameters
        ==========
        contigs_dbs : Iterable[str]
            File paths to expected contigs databases.
        """
        invalid_paths: List[str] = []
        invalid_filetypes: List[str] = []
        for contigs_db in contigs_dbs:
            if not os.path.exists(contigs_db):
                invalid_paths.append(contigs_db)
            if invalid_paths:
                continue

            contigs_db_info = dbinfo.ContigsDBInfo(contigs_db, dont_raise=True, expecting='contigs')
            if contigs_db_info is None:
                invalid_filetypes.append(contigs_db)
            if invalid_filetypes:
                continue

        if invalid_paths:
            paths = ', '.join([f'{path}' for path in invalid_paths])
            raise ConfigError(
                f"There were no files at the following expected contigs database paths: {paths}"
            )

        if invalid_filetypes:
            paths = ', '.join([f'{path}' for path in invalid_filetypes])
            raise ConfigError(
                "The files at the following expected contigs database paths are not contigs "
                f"databases: {paths}"
            )

    @staticmethod
    def _check_contigs_dbs_ko_annotation(contigs_dbs: Iterable[str]) -> None:
        unannotated: List[str] = []
        for contigs_db in contigs_dbs:
            contigs_db_info = dbinfo.ContigsDBInfo(contigs_db, expecting='contigs')
            if 'KOfam' not in contigs_db_info.get_functional_annotation_sources():
                unannotated.append(contigs_db)
            if unannotated:
                continue

        if unannotated:
            paths = ', '.join([f'{path}' for path in unannotated])
            raise ConfigError(
                "The following contigs databases were never annotated with KOs, but this can be "
                f"rectified by running `anvi-run-kegg-kofams` on them: {paths}"
            )

    @staticmethod
    def _check_pan_db(pan_db: str) -> None:
        """
        Check the validity of an expected pan database.

        Parameters
        ==========
        pan_db : str
            File path to an expected pan database.
        """
        if not os.path.exists(pan_db):
            raise ConfigError(
                f"There was no file at the following expected pan database path: '{pan_db}'"
            )

        pan_db_info = dbinfo.PanDBInfo(pan_db, dont_raise=True, expecting='pan')
        if pan_db_info is None:
            raise ConfigError(
                "The file at the following expected pan database path is not a pan database: "
                f"'{pan_db}'"
            )

    def _find_maps(self, output_dir: str, prefix: str, patterns: List[str] = None) -> List[str]:
        """
        Find the numeric IDs of maps to draw given the file prefix, checking that the map can be
        drawn in the target output direcotry.

        Parameters
        ==========
        output_dir : str
            Path to the output directory in which pathway map PDF files are drawn. The directory is
            created if it does not exist.

        prefix : str
            Output filenames are formatted as <prefix>_<pathway_number>.pdf or
            <prefix>_<pathway_number>_<pathway_name>.pdf.

        patterns : List[str], None
            Regex patterns of pathway numbers, which are five digits.
        """
        if patterns is None:
            pathway_numbers = self.available_pathway_numbers
        else:
            pathway_numbers = self._get_pathway_numbers_from_patterns(patterns)

        if not self.overwrite_output:
            for pathway_number in pathway_numbers:
                if self.name_files:
                    pathway_name = '_' + self._get_filename_pathway_name(pathway_number)
                else:
                    pathway_name = ''
                out_path = os.path.join(output_dir, f'{prefix}_{pathway_number}{pathway_name}.pdf')
                if os.path.exists(out_path):
                    raise ConfigError(
                        f"Output files would be overwritten in the output directory, {output_dir}. "
                        "Either delete the contents of the directory, or use the option to "
                        "overwrite output destinations."
                    )

        return pathway_numbers

    def _get_pathway_numbers_from_patterns(self, patterns: Iterable[str]) -> List[str]:
        """
        Among pathways available in the KEGG data directory, get those with ID numbers matching the
        given regex patterns.

        Parameters
        ==========
        patterns : Iterable[str]
            Regex patterns of pathway numbers, which are five digits.

        Returns
        =======
        List[str]
            Pathway numbers matching the regex patterns.
        """
        pathway_numbers: List[str] = []
        for pattern in patterns:
            for available_pathway_number in self.available_pathway_numbers:
                if re.match(pattern, available_pathway_number):
                    pathway_numbers.append(available_pathway_number)

        # Maintain the order of pathway numbers recovered from patterns.
        seen = set()
        return [
            pathway_number for pathway_number in pathway_numbers
            if not (pathway_number in seen or seen.add(pathway_number))
        ]

    def _draw_map_kos_single_color(
        self,
        pathway_number: str,
        ko_ids: Iterable[str],
        color_hexcode: str,
        output_dir: str,
        draw_map_lacking_kos: bool = False
    ) -> bool:
        """
        Draw a pathway map, highlighting reactions containing select KOs in a single color.

        Parameters
        ==========
        pathway_number : str, None
            Numeric ID of the map to draw.

        ko_ids : Iterable[str]
            Select KOs, any of which in the map are colored.

        color_hexcode : str
            This is the color, by default green, for reactions containing provided KOs. A reaction
            on a map can correspond to one or more KOs, and a KO can annotate one or more sequences
            in a contigs database. In global and overview maps, reaction lines are colored. In
            standard maps, reaction boxes or lines are colored.

        output_dir : str
            Path to an existing output directory in which map PDF files are drawn.

        draw_map_lacking_kos : bool, False
            If False, by default, only draw the map if it contains any of the select KOs. If True,
            draw the map regardless, meaning that nothing may be highlighted.

        Returns
        =======
        bool
            True if the map was drawn, False if the map was not drawn because it did not contain any
            of the select KOs and 'draw_map_lacking_kos' was False.
        """
        pathway = self._get_pathway(pathway_number)

        select_entries = pathway.get_entries(kegg_ids=ko_ids)
        if not select_entries and not draw_map_lacking_kos:
            return False

        # Set the color of Graphics elements for reactions containing select KOs. For other Graphics
        # elements, change the 'fgcolor' attribute to a nonsense value of '0' to ensure that the
        # elements with the prioritized color can be distinguished from other elements. Also, in
        # overview and standard maps, widen lines from the base map default of 1.0.
        all_entries = pathway.get_entries(entry_type='ortholog')
        select_uuids = [entry.uuid for entry in select_entries]
        for entry in all_entries:
            if entry.uuid in select_uuids:
                for uuid in entry.children['graphics']:
                    graphics: kgml.Graphics = pathway.uuid_element_lookup[uuid]
                    if pathway.is_global_map:
                        assert graphics.type == 'line'
                        graphics.fgcolor = color_hexcode
                        graphics.bgcolor = '#FFFFFF'
                    elif pathway.is_overview_map:
                        assert graphics.type == 'line'
                        graphics.fgcolor = color_hexcode
                        graphics.bgcolor = '#FFFFFF'
                        graphics.width = 5.0
                    else:
                        if graphics.type == 'rectangle':
                            graphics.fgcolor = '#000000'
                            graphics.bgcolor = color_hexcode
                        elif graphics.type == 'line':
                            graphics.fgcolor = color_hexcode
                            graphics.bgcolor = '#FFFFFF'
                            graphics.width = 5.0
                        else:
                            raise AssertionError(
                                "Ortholog entries are assumed to have Graphics elements of type "
                                "'rectangle' or 'line', not the encountered type, "
                                f"'{graphics.type}'."
                            )
            else:
                for uuid in entry.children['graphics']:
                    graphics: kgml.Graphics = pathway.uuid_element_lookup[uuid]
                    graphics.fgcolor = '0'

        # Set the color priority so that the colored reactions are prioritized for display on top.
        # Recolor "unprioritized" reactions to a background color. In global and overview maps,
        # recolor circles to reflect the colors of prioritized reactions involving the compounds.
        color_priority: Dict[str, Dict[str, Dict[Tuple[str, str], float]]] = {}
        if pathway.is_global_map:
            color_priority['ortholog'] = {'line': {(color_hexcode, '#FFFFFF'): 1.0}}
            recolor_unprioritized_entries = 'g'
            color_associated_compounds = 'high'
        elif pathway.is_overview_map:
            color_priority['ortholog'] = {'line': {(color_hexcode, '#FFFFFF'): 1.0}}
            recolor_unprioritized_entries = 'w'
            color_associated_compounds = 'high'
        else:
            color_priority['ortholog'] = {
                'rectangle': {('#000000', color_hexcode): 1.0},
                'line': {(color_hexcode, '#FFFFFF'): 1.0}
            }
            recolor_unprioritized_entries = 'w'
            color_associated_compounds = None
        pathway.set_color_priority(
            color_priority,
            recolor_unprioritized_entries=recolor_unprioritized_entries,
            color_associated_compounds=color_associated_compounds
        )

        # Draw the map.
        if self.name_files:
            pathway_name = '_' + self._get_filename_pathway_name(pathway_number)
        else:
            pathway_name = ''
        out_path = os.path.join(output_dir, f'kos_{pathway_number}{pathway_name}.pdf')
        if os.path.exists(out_path) and self.overwrite_output:
            os.remove(out_path)
        else:
            filesnpaths.is_output_file_writable(out_path, ok_if_exists=False)
        self.drawer.draw_map(pathway, out_path)
        return True

    def _draw_map_kos_original_color(
        self,
        pathway_number: str,
        ko_ids: Iterable[str],
        output_dir: str,
        draw_map_lacking_kos: bool = False
    ) -> bool:
        """
        Draw a pathway map, highlighting reactions containing select KOs in the color or colors
        originally used in the reference map.

        Parameters
        ==========
        pathway_number : str, None
            Numeric ID of the map to draw.

        ko_ids : Iterable[str]
            Select KOs, any of which in the map are colored.

        output_dir : str
            Path to an existing output directory in which map PDF files are drawn.

        draw_map_lacking_kos : bool, False
            If False, by default, only draw the map if it contains any of the select KOs. If True,
            draw the map regardless, meaning that nothing may be highlighted.

        Returns
        =======
        bool
            True if the map was drawn, False if the map was not drawn because it did not contain any
            of the select KOs and 'draw_map_lacking_kos' was False.
        """
        pathway = self._get_pathway(pathway_number)

        select_entries = pathway.get_entries(kegg_ids=ko_ids)
        if not select_entries and not draw_map_lacking_kos:
            return False

        # Set "secondary" colors of ortholog Graphics elements for reactions containing select KOs:
        # white background color of lines or black foreground text of rectangles. For other Graphics
        # elements, change the 'fgcolor' attribute to a nonsense value to ensure that the elements
        # with prioritized colors can be distinguished from other elements. Also, in overview and
        # standard maps, widen lines from the base map default of 1.0.
        all_entries = pathway.get_entries(entry_type='ortholog')
        select_uuids = [entry.uuid for entry in select_entries]
        prioritized_colors: Dict[str, List[Tuple[str, str]]] = {}
        for entry in all_entries:
            if entry.uuid in select_uuids:
                for uuid in entry.children['graphics']:
                    graphics: kgml.Graphics = pathway.uuid_element_lookup[uuid]
                    if pathway.is_global_map:
                        assert graphics.type == 'line'
                        graphics.bgcolor = '#FFFFFF'
                    elif pathway.is_overview_map:
                        assert graphics.type == 'line'
                        graphics.bgcolor = '#FFFFFF'
                        graphics.width = 5.0
                    else:
                        if graphics.type == 'rectangle':
                            graphics.fgcolor = '#000000'
                        elif graphics.type == 'line':
                            graphics.bgcolor = '#FFFFFF'
                            graphics.width = 5.0
                        else:
                            raise AssertionError(
                                "Ortholog entries are assumed to have Graphics elements of type "
                                "'rectangle' or 'line', not the encountered type, "
                                f"'{graphics.type}'."
                            )
                    try:
                        graphics_type_prioritized_colors = prioritized_colors[graphics.type]
                    except:
                        prioritized_colors[graphics.type] = graphics_type_prioritized_colors = []
                    graphics_type_prioritized_colors.append((graphics.fgcolor, graphics.bgcolor))
            else:
                for uuid in entry.children['graphics']:
                    graphics: kgml.Graphics = pathway.uuid_element_lookup[uuid]
                    graphics.fgcolor = '0'

        # By default, global maps but not overview and standard maps display reaction graphics in
        # more than one color. Give higher priority to reaction entries that are encountered later
        # (occur further down in the KGML file), and would thus be rendered above earlier reactions.
        color_priority: Dict[str, Dict[str, Dict[Tuple[str, str], float]]] = {'ortholog': {}}
        for graphics_type, graphics_type_prioritized_colors in prioritized_colors.items():
            seen = set()
            unique_prioritized_colors = [
                colors for colors in graphics_type_prioritized_colors
                if not (colors in seen or seen.add(colors))
            ]
            priorities = np.linspace(0, 1, len(unique_prioritized_colors) + 1)[1: ]
            graphics_type_color_priority = {
                colors: priority for colors, priority in zip(unique_prioritized_colors, priorities)
            }
            color_priority['ortholog'][graphics_type] = graphics_type_color_priority

        # Recolor "unprioritized" reactions to a background color. In global and overview maps,
        # recolor circles to reflect the colors of prioritized reactions involving the compounds.
        if pathway.is_global_map:
            recolor_unprioritized_entries = 'g'
            color_associated_compounds = 'high'
        elif pathway.is_overview_map:
            recolor_unprioritized_entries = 'w'
            color_associated_compounds = 'high'
        else:
            recolor_unprioritized_entries = 'w'
            color_associated_compounds = None
        pathway.set_color_priority(
            color_priority,
            recolor_unprioritized_entries=recolor_unprioritized_entries,
            color_associated_compounds=color_associated_compounds
        )

        # Draw the map.
        if self.name_files:
            pathway_name = '_' + self._get_filename_pathway_name(pathway_number)
        else:
            pathway_name = ''
        out_path = os.path.join(output_dir, f'kos_{pathway_number}{pathway_name}.pdf')
        if os.path.exists(out_path) and self.overwrite_output:
            os.remove(out_path)
        else:
            filesnpaths.is_output_file_writable(out_path, ok_if_exists=False)
        self.drawer.draw_map(pathway, out_path)
        return True

    def _draw_map_kos_membership(
        self,
        pathway_number: str,
        ko_membership: Dict[str, List[str]],
        color_priority: Dict[str, float],
        output_dir: str,
        colormap: mcolors.Colormap,
        source_combos: List[Tuple[str]] = None,
        draw_map_lacking_kos: bool = False
    ) -> bool:
        """
        Draw a pathway map, coloring reactions by their membership in sources.

        For a pangenome, reactions are colored by genomes containing consensus KOs in the reaction.
        For contigs databases, reactions are colored by databases containing KOs in the reaction. By
        default, with 'source_combos' being None, coloring reflects the count of genomes or
        databases rather than actual genome or database membership.

        In global and overview maps, compounds involved in colored reactions are given the color of
        the reaction with the highest priority.

        Parameters
        ==========
        pathway_number : str
            Numeric ID of the map to draw.

        ko_membership : Dict[str, List[str]]
            Keys are KO IDs. Values are lists of "sources:" genome names or project names of contigs
            databases.

            A KO can annotate more than one gene cluster in a pangenome; a list contains the names
            of genomes contributing genes to clusters represented by the KO.

        color_priority : Dict[str, float]
            Keys are color hex codes. If 'by_count' is True, there should be one color for each
            possible number of genomes or databases. If 'by_count' is False, there should be one
            color for each individual genome or database and combination thereof. Values are
            priorities. KOs with higher priority colors are drawn over KOs with lower priority
            colors.

        output_dir : str
            Path to an existing output directory in which map PDF files are drawn.

        colormap : matplotlib.colors.Colormap
            This colormap is used to interpolate the colors of compounds involved in reactions with
            color-prioritized KOs. Colors in the color_priority arguments should be drawn from this
            colormap.

        source_combos : List[Tuple[str]], None
            With the default argument value of None, reactions are colored by number of pangenomic
            genomes or contigs databases containing the reaction. A list of "source combination"
            tuples can be provided instead to color explicitly by genome or database membership.
            Tuples should consist of source names (genome names or database project names) and their
            combinations, e.g., [('A', ), ('B', ), ('C', ), ('A', 'B'), ('A', 'C'), ('B', 'C'),
            ('A', 'B', 'C')].

        draw_map_lacking_kos : bool, False
            If False, by default, only draw the map if it contains any of the select KOs. If True,
            draw the map regardless, meaning that nothing may be highlighted.

        Returns
        =======
        bool
            True if the map was drawn, False if the map was not drawn because it did not contain any
            of the select KOs and 'draw_map_lacking_kos' was False.
        """
        pathway = self._get_pathway(pathway_number)

        combo_lookup: Dict[Tuple[str], Tuple[str]] = {}
        if source_combos is not None:
            for combo in source_combos:
                combo_lookup[tuple(sorted(combo))] = combo

        entries = pathway.get_entries(kegg_ids=ko_membership)
        if not entries and not draw_map_lacking_kos:
            return False

        # Change the colors of the KO graphics. A reaction Entry can represent multiple KOs. Also,
        # in overview and standard maps, widen lines from the base map default of 1.0.
        color_hexcodes = list(color_priority)
        for entry in entries:
            source_names = []
            for kegg_name in entry.name.split():
                split_kegg_name = kegg_name.split(':')
                kegg_id = split_kegg_name[1]
                try:
                    source_names += ko_membership[kegg_id]
                except KeyError:
                    continue
            assert len(source_names)

            if source_combos is None:
                color_hexcode = color_hexcodes[len(set(source_names)) - 1]
            else:
                source_combo = combo_lookup[tuple(sorted(set(source_names)))]
                color_hexcode = color_hexcodes[source_combos.index(source_combo)]
            for uuid in entry.children['graphics']:
                graphics: kgml.Graphics = pathway.uuid_element_lookup[uuid]
                if pathway.is_global_map:
                    assert graphics.type == 'line'
                    graphics.fgcolor = color_hexcode
                    graphics.bgcolor = '#FFFFFF'
                elif pathway.is_overview_map:
                    assert graphics.type == 'line'
                    graphics.fgcolor = color_hexcode
                    graphics.bgcolor = '#FFFFFF'
                    graphics.width = 5.0
                else:
                    if graphics.type == 'rectangle':
                        graphics.fgcolor = '#000000'
                        graphics.bgcolor = color_hexcode
                    elif graphics.type == 'line':
                        graphics.fgcolor = color_hexcode
                        graphics.bgcolor = '#FFFFFF'
                        graphics.width = 5.0
                    else:
                        raise AssertionError(
                            "Ortholog entries are assumed to have Graphics elements of type "
                            f"'rectangle' or 'line', not the encountered type, '{graphics.type}'."
                        )

        # Set the color priorities of entries for proper overlaying in the image. Recolor
        # "unprioritized" KO graphics to a background color. In global and overview maps, recolor
        # circles to reflect the colors of prioritized reactions involving the compounds.
        ortholog_color_priority: Dict[str, Dict[Tuple[str, str], float]] = {}
        if pathway.is_global_map:
            ortholog_color_priority['line'] = line_color_priority = {}
            for color_hexcode, priority in color_priority.items():
                line_color_priority[(color_hexcode, '#FFFFFF')] = priority
            pathway.set_color_priority(
                {'ortholog': ortholog_color_priority},
                recolor_unprioritized_entries='g',
                color_associated_compounds='high',
                colormap=colormap
            )
        elif pathway.is_overview_map:
            ortholog_color_priority['line'] = line_color_priority = {}
            for color_hexcode, priority in color_priority.items():
                line_color_priority[(color_hexcode, '#FFFFFF')] = priority
            pathway.set_color_priority(
                {'ortholog': ortholog_color_priority},
                recolor_unprioritized_entries='w',
                color_associated_compounds='high',
                colormap=colormap
            )
        else:
            ortholog_color_priority['rectangle'] = rectangle_color_priority = {}
            ortholog_color_priority['line'] = line_color_priority = {}
            for color_hexcode, priority in color_priority.items():
                rectangle_color_priority[('#000000', color_hexcode)] = priority
                line_color_priority[(color_hexcode, '#FFFFFF')] = priority
            pathway.set_color_priority(
                {'ortholog': ortholog_color_priority},
                recolor_unprioritized_entries='w'
            )

        # Draw the map.
        if self.name_files:
            pathway_name = '_' + self._get_filename_pathway_name(pathway_number)
        else:
            pathway_name = ''
        out_path = os.path.join(output_dir, f'kos_{pathway_number}{pathway_name}.pdf')
        if os.path.exists(out_path) and self.overwrite_output:
            os.remove(out_path)
        else:
            filesnpaths.is_output_file_writable(out_path, ok_if_exists=False)
        self.drawer.draw_map(pathway, out_path)
        return True

    def _get_pathway(self, pathway_number: str) -> kgml.Pathway:
        """
        Get a Pathway object for the KGML file used in drawing a pathway map.

        Parameters
        ==========
        pathway_number : str
            Numeric ID of the map to draw.

        Returns
        =======
        kgml.Pathway
            Representation of the KGML file as an object.
        """
        # KOs correspond to arrows rather than boxes in global and overview maps.
        is_global_map = False
        is_overview_map = False
        if re.match(kegg.GLOBAL_MAP_ID_PATTERN, pathway_number):
            is_global_map = True
        elif re.match(kegg.OVERVIEW_MAP_ID_PATTERN, pathway_number):
            is_overview_map = True

        # A 1x resolution global 'KO' image is used as the base of the drawing, whereas a 2x
        # overview or standard 'map' image is used as the base. The global 'KO' image grays out
        # all reaction arrows that are not annotated by KO ID. Select the KGML file accordingly.
        if is_global_map:
            kgml_path = os.path.join(
                self.kegg_context.kgml_1x_ko_dir, f'ko{pathway_number}.xml'
            )
        else:
            kgml_path = os.path.join(
                self.kegg_context.kgml_2x_ko_dir, f'ko{pathway_number}.xml'
            )
        pathway = self.xml_ops.load(kgml_path)

        return pathway

    def _get_filename_pathway_name(self, pathway_number: str) -> str:
        """
        Format the pathway name corresponding to the number to include in file paths.

        Replace all non-alphanumeric characters except parentheses, brackets, and curly braces with
        underscores. Replace multiple consecutive underscores with a single underscore. Strip
        leading and trailing underscores.

        Parameters
        ==========
        pathway_number : str
            Numeric ID of a pathway map.

        Returns
        =======
        str
            Altered version of the pathway name.
        """
        try:
            pathway_name = self.pathway_names[pathway_number]
        except KeyError:
            raise ConfigError(
                f"The pathway number, '{pathway_number}', is not recognized in the table of KEGG "
                "pathway names set up in the KEGG data directory, which can be found here: "
                f"'{self.kegg_context.kegg_pathway_list_file}'."
            )

        altered = re.sub(r'[^a-zA-Z0-9()\[\]\{\}]', '_', pathway_name)
        altered = re.sub(r'_+', '_', altered)
        altered = altered.strip('_')

        return altered

    def _draw_colorbar(
        self,
        colors: Iterable,
        out_path: str,
        color_labels: Iterable[str] = None,
        label: str = None
    ) -> None:
        """
        Save a standalone colorbar to a file.

        Parameters
        ==========
        colors : Iterable
            Sequence of Matplotlib color specifications for matplotlib.colors.ListedColormap color
            parameter.

        out_path : str
            Path to PDF output file.

        color_labels : Iterable[str], None
            Labels corresponding to each color.

        label : str, None
            Overall colorbar label.
        """
        if color_labels is not None:
            assert len(colors) == len(color_labels)

        fig, ax = plt.subplots(figsize=(1, 6))

        cmap = mcolors.ListedColormap(colors)
        norm = mcolors.BoundaryNorm(boundaries=range(len(colors) + 1), ncolors=len(colors))

        cb = plt.colorbar(
            plt.cm.ScalarMappable(norm=norm, cmap=cmap),
            cax=ax,
            orientation='vertical'
        )

        # Don't show tick marks.
        cb.ax.tick_params(size=0)

        if color_labels:
            # Calculate appropriate font size of tick labels based on color segment height.
            height_in_data_coords = 1 / len(colors)
            height_in_points = (
                ax.transData.transform((0, height_in_data_coords)) - ax.transData.transform((0, 0))
            )
            if height_in_points[1] < 10:
                tick_font_size = height_in_points[1] * 2
            else:
                tick_font_size = min(height_in_points[1], 24)

            cb.set_ticks(np.arange(len(colors)) + 0.5)
            cb.set_ticklabels(color_labels, fontsize=tick_font_size)

        if label:
            label_font_size = min(tick_font_size * 1.25, 30)
            cb.set_label(label, rotation=270, labelpad=label_font_size * 1.25, fontsize=label_font_size)

        if os.path.exists(out_path) and self.overwrite_output:
            os.remove(out_path)
        else:
            filesnpaths.is_output_file_writable(out_path, ok_if_exists=False)
        plt.savefig(out_path, format='pdf', bbox_inches='tight')
        plt.close()

    def make_grid(
        self,
        in_paths: Iterable[str],
        out_path: str,
        labels: Iterable[str] = None,
        landscape: bool = False,
        margin: float = 10.0
    ) -> None:
        """
        Write a PDF containing a grid of input PDF images.

        Parameters
        ==========
        in_paths : Iterable[str]
            Paths to input PDFs.

        out_path : str
            Path to output PDF.

        labels : Iterable[str], None
            Labels displayed over grid cells corresponding to input files.

        landscape : bool, False
            Page layout is portrait if False, landscape if True.

        margin : float, 10.0
            Minimum space between cells.
        """
        if labels:
            assert len(in_paths) == len(labels)

        # Find the number of rows and columns in the grid.
        cols = math.ceil(math.sqrt(len(in_paths)))
        rows = math.ceil(len(in_paths) / cols)

        # Find the width and height of each cell.
        width, height = fitz.paper_size(f'{"letter-l" if landscape else "letter"}')
        cell_width = (width - (cols + 1) * margin) / cols
        cell_height = (height - (rows + 1) * margin) / rows

        fontsize = margin * 0.8

        # Create a new PDF document.
        output_doc = fitz.open()
        output_page = output_doc.new_page(width=width, height=height)

        # Loop through input PDF files, placing them in the grid.
        for i, pdf_path in enumerate(in_paths):
            pdf_doc = fitz.open(pdf_path)
            page = pdf_doc.load_page(0)

            # Calculate position in the grid.
            row = i // cols
            col = i % cols
            x = margin + col * (cell_width + margin)
            y = margin + row * (cell_height + margin)

            # Resize the input PDF to the cell by the longest dimension, maintaining aspect ratio.
            input_aspect_ratio = page.rect.width / page.rect.height
            if input_aspect_ratio > 1:
                draw_width = cell_width
                draw_height = cell_width / input_aspect_ratio
            else:
                draw_height = cell_height
                draw_width = cell_height * input_aspect_ratio

            # If the resized shorter side still exceeds the cell size, resize by the shorter side.
            if draw_width > cell_width:
                draw_width = cell_width
                draw_height = cell_width / input_aspect_ratio
            if draw_height > cell_height:
                draw_height = cell_height
                draw_width = cell_height * input_aspect_ratio

            # Find upper left drawing coordinates.
            draw_x = x + (cell_width - draw_width) / 2
            draw_y = y + (cell_height - draw_height) / 2

            # Place the input PDF.
            rect = fitz.Rect(draw_x, draw_y, draw_x + draw_width, draw_y + draw_height)
            output_page.show_pdf_page(rect, pdf_doc, 0)

            if labels:
                # Draw labels above each image.
                label = labels[i]
                label_x = draw_x
                label_y = draw_y
                output_page.insert_text((label_x, label_y), label, fontsize=fontsize)

        output_doc.save(out_path)

    @property
    def quiet(self):
        return self._quiet

    @quiet.setter
    def quiet(self, new_value: bool):
        self._quiet = new_value
        self.run.verbose = not self.quiet
        self.progress.verbose = not self.quiet
