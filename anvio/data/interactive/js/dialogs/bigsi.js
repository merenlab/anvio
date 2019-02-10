/**
 *  Search random sequence on BIGSI (http://www.bigsi.io)
 *
 *  Author: Özcan Esen <ozcanesen@gmail.com>
 *  Credits: A. Murat Eren
 *  Copyright 2018, The anvio Project
 *
 * This file is part of anvi'o (<https://github.com/merenlab/anvio>).
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


function BIGSI(split_name, sequence) {
    this.dialog = document.createElement('div');
    this.dialog.setAttribute('class', 'modal fade in');

    this.split_name = split_name;
    this.sequence = sequence;

    this.dialog.innerHTML = `<div class="modal-dialog" style="width: 800px;">
            <div class="modal-content">
                <div class="modal-header">
                    <button class="close" data-dismiss="modal" type="button"><span>&times;</span></button>
                    <h4 class="modal-title">BIGSI Search</h4>
                </div>

                <div class="modal-body">
                    <div class='info-section'></div>

                    <hr>

                    <div class="results">
                        <center><img src="images/loading.gif" /> <b>Fetching results...</b><br /></center>
                    </div>
                </div>

                <div class="modal-footer">
                    <button class="btn btn-primary" type="button">Do it again</button>
                    <button class="btn btn-default" data-dismiss="modal" type="button">Close</button>
                </div>
            </div>
        </div>`;

    this.dialog.querySelector('.btn-primary').addEventListener('click', (event) => { this.Search(); });
    $(this.dialog).modal('show').on('hidden.bs.modal', () => this.dialog.remove());
};


BIGSI.prototype.Search = function() {
    this.dialog.querySelector('.results').innerHTML = '<center><img src="images/loading.gif" /> <b>Fetching results...</b><br /></center>';

    let sub_sequence = '';
    const search_length = 150;
    let start = Math.floor(Math.random() * (this.sequence.length - search_length));

    if (this.sequence.length < search_length) {
        sub_sequence = this.sequence;
    } else {
        sub_sequence = this.sequence.substring(start, start + search_length);
    }

    this.dialog.querySelector('.info-section').innerHTML = `<span style="word-wrap: break-word">Searching a random ${search_length} bp sequence from ${this.sequence.length} bp long split "<b style="word-wrap: break-word">${this.split_name}</b>" on <a target="_blank" href="http://www.bigsi.io">BIGSI</a></span>
                    <div style="font-family: monospace; word-break: break-all; background: #f1f1f1; margin: 10px 5px; padding: 6px;">${sub_sequence}</div>`;

    $.ajax({
        type: 'GET',
        cache: false,
        url: `http://api.bigsi.io/search?seq=${sub_sequence}&threshold=1&score=0`,
        success: (data) => {
            let table = '';
            let results = data[Object.keys(data)[0]]['results'];

            if (!Object.keys(results).length) {
                table += '<div style="background: #ffdfdf; margin: 10px 5px; padding: 6px;">Tataaaa :) No results for your search :(</div>`';
            } else {
                table += '<table class="table table-condensed table-striped" style="margin-top: 20px; margin-bottom: 20px;"><thead> <tr><th>% k-mers found</th> <th>Accession</th> <th>Hits</th> </tr> </thead> <tbody>';
                for (accession in results) {
                    let item = results[accession];
                    table += `<tr><td>${item['percent_kmers_found']}</td><td><a target="_blank" href="https://www.ebi.ac.uk/ena/data/view/${accession}">${accession}</a></td> <td>${item['species']}</td> </tr>`;
                }

                table += '</tbody></table>';
                table += '<div style="background: #f1e9de; margin: 10px 5px; padding: 6px;">If you find these results helpful for your research, please cite <i>"Ultrafast search of all deposited bacterial and viral genomic data"</i> by Bradley et al (<a target="_blank" href="https://doi.org/10.1038/s41587-018-0010-1">doi:10.1038/s41587-018-0010-1</a>).</div>`';
            }

            this.dialog.querySelector('.results').innerHTML = table;
        },
        error: () => {
            this.dialog.querySelector('.results').innerHTML = "<center><b>BIGSI returned nothing :/ (probably there was an error, but anvi'o is Jon Snow).</b></center>";
        }
    });
}


