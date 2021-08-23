/**
 * Javascript library for anvi'o genome view
 *
 *  Authors: Isaac Fink <iafink@uchicago.edu>
 *           Matthew Klein <mtt.l.kln@gmail.com>
 *           A. Murat Eren <a.murat.eren@gmail.com>
 *
 * Copyright 2021, The anvi'o project (http://anvio.org)
 *
 * Anvi'o is a free software. You can redistribute this program
 * and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 *
 * You should have received a copy of the GNU General Public License
 * along with anvi'o. If not, see <http://opensource.org/licenses/GPL-3.0>.
 *
 * @license GPL-3.0+ <http://opensource.org/licenses/GPL-3.0>
 */

/**
 * File Overview : This file contains temporary functions and processes used to test various functionality within 
 * genomeview. They are gathered here so as not to pollute the other genomeview files. 
 */

function drawTestShades() {
  shadeGeneClusters(["GC_00000034", "GC_00000097", "GC_00000002"], { "GC_00000034": "green", "GC_00000097": "red", "GC_00000002": "purple" });
}

function generateMockADL() {
  for (let i = 0; i < genomeData.genomes.length; i++) { // generate mock additional data layer content
    let gcContent = []
    let coverage = []
    for (let j = 0; j < genomeMax; j++) {
      gcContent.push(Math.floor(Math.random() * 45))
      coverage.push(Math.floor(Math.random() * 45))
    }
    let genomeLabel = Object.keys(genomeData.genomes[i][1]['contigs']['info'])[0];
    let additionalDataObject = {
      'genome': genomeLabel,
      'coverage': coverage,
      'coverage-color': 'blue',
      'gcContent': gcContent,
      'gcContent-color': 'purple',
      'ruler': true // TODO: store any genome-specific scale data here
    }
    stateData['additional-data-layers'].push(additionalDataObject)
  }
}
function generateMockGenomeOrder() {
  stateData['genome-order-method'] = [{
    'name': 'cats',
    'ordering': 'some order'
  }, {
    'name': 'dogs',
    'ordering': 'some other order'
  }, {
    'name': 'birds',
    'ordering': 'beaks to tails'
  }]
}