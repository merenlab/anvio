# -*- coding: utf-8 -*-
# pylint: disable=line-too-long

__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2020, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__maintainer__ = "Samuel Miller"
__email__ = "samuelmiller10@gmail.com"
__status__ = "Development"


THREEPRIME_VARIANTS = ['CCA', 'CC', 'C',
                       'CCAA', 'CCAC', 'CCAG', 'CCAT',
                       'CCAAA', 'CCAAC', 'CCAAG', 'CCAAT',
                       'CCACA', 'CCACC', 'CCACG', 'CCACT',
                       'CCAGA', 'CCAGC', 'CCAGG', 'CCAGT',
                       'CCATA', 'CCATC', 'CCATG', 'CCATT']

TRNA_FEATURE_NAMES = ['trna_his_position_0',
                      'acceptor_stem',
                      'fiveprime_acceptor_stem_sequence',
                      'position_8',
                      'position_9',
                      'd_arm',
                      'd_stem',
                      'fiveprime_d_stem_sequence',
                      'd_loop',
                      'threeprime_d_stem_sequence',
                      'position_26',
                      'anticodon_arm',
                      'anticodon_stem',
                      'fiveprime_anticodon_stem_sequence',
                      'anticodon_loop',
                      'threeprime_anticodon_stem_sequence',
                      'v_loop',
                      't_arm',
                      't_stem',
                      'fiveprime_t_stem_sequence',
                      't_loop',
                      'threeprime_t_stem_sequence',
                      'threeprime_acceptor_stem_sequence',
                      'discriminator',
                      'acceptor']

# the longest known tRNA (selenocysteine) is 101 bp
# ref: Santesmasses, Mariotti & Guigo, 2017, "Computational identification of the selenocysteine tRNA (tRNASec) in genomes," PLoS Comp. Biol., 13(2): e1005383
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5330540/
LONGEST_KNOWN_TRNA_LENGTH = 101
