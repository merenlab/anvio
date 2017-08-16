
function searchContigs() 
{
    var svalue = $('#searchValue').val();

    if (svalue == "")
    {
        alert("Search value shouldn't be empty.");
        return;
    }
    var column = $('#searchLayerList').val();
    search_column = layerdata[0][column];
    var operator = $('#searchOperator').val();
    
    if (operator < 6)
    {
        var operator_text = $('#searchOperator option:selected').text();

        // logical operator
        var _pre = "layerdata[";
        var _post = "][" + column + "] " + operator_text + " \"" + svalue.trim() + "\"";

    }
    else if (operator == 6)
    {
        // contains
        var _pre = "layerdata[";
        var _post = "][" + column + "].toString().indexOf(\"" + svalue + "\") != -1";
    }

    var _len = layerdata.length;
    var _counter = 0;
    search_results = [];

    $('#search_result_message').html("Searching...");

    for (var row=1; row < _len; row++)
    {
        if (layerdata[row][column]==null)
            continue;

        if (eval(_pre + row + _post)){
            search_results.push({'split': layerdata[row][0], 'value': layerdata[row][column]});
            _counter++;
        }
    }
    $('#search_result_message').html(_counter + " splits found.");
}

function searchFunctions() {
    $.ajax({
        type: 'POST',
        cache: false,
        url: '/data/search_functions?timestamp=' + new Date().getTime(),
        data: {terms: $('#searchFunctionsValue').val()},
        success: function(data) {
            if (data['status'] == 0) {
                $('.search-message').hide();
                $('.search-message').html('');

                search_results = [];
                search_column = 'Annotation';

                for (var i=0; i < data['results'].length; i++) {
                    var _gene_caller_id = data['results'][i][0];
                    var _source         = data['results'][i][1];
                    var _annotation     = data['results'][i][2];
                    var _search_term    = data['results'][i][3];
                    var _split_name     = data['results'][i][4];

                    var _beginning = _annotation.toLowerCase().indexOf(_search_term.toLowerCase());
                    _annotation = [_annotation.slice(0, _beginning), 
                                   '<mark>', 
                                   _annotation.slice(_beginning, _beginning + _search_term.length), 
                                   '</mark>', 
                                   _annotation.slice(_beginning + _search_term.length, _annotation.length)
                                   ].join("");

                    search_results.push({'split': data['results'][i][4], 'value': '<b>Gene caller id:</b> ' + _gene_caller_id +
                                                                                  '</br><b>Source:</b> ' + _source +
                                                                                  '</br><b>Annotation:</b> ' + _annotation});
                }
                $('#search_result_message_functions').html(data['results'].length + " results found.");
            } else {
                $('.search-message').show();
                $('.search-message').html(data['message']);
            };
        }
    });
}

function showSearchResult() {
    var clear_link = '<a href="#" onclick="$(\'.search-results-display, #search-results-table-search-item, #search-results-table-search-name, #search-results-table-header\').html(\'\');">(clear)</a>';
    $("#search-results-table-header").html('<h4>Search results ' + clear_link + '</h4>');
    $("#search-results-table-search-name").html('Split name');
    $("#search-results-table-search-item").html(search_column);

    var rows = "";
    var _len = search_results.length;
    for (var i=0; i < _len; i++)
    {
        var col1 = search_results[i]['split'];
        var col2 = search_results[i]['value'];

        rows = rows + `<tr><td data-value="${ col1 }""><a href='#' class='no-link' onclick='highlightSplit("${ col1 }");'>${col1}</a></td><td data-value="${ col2 }">${ col2 }</td></tr>`;
    }
    $("#tbody_search_body").html(rows);
}

function highlightResult() {
    // check if tree exists
    if ($.isEmptyObject(label_to_node_map)) {
        alert('Draw tree first.');
        return;
    }

    highlighted_splits = [];

    for (var i=0; i < search_results.length; i++) {
        var _contig_name = search_results[i]['split'];
        highlighted_splits.push(_contig_name);
    }

    redrawBins(); 
}

function highlightSplit(name) {
    // check if tree exists
    if ($.isEmptyObject(label_to_node_map)) {
        alert('Draw tree first.');
        return;
    }

    highlighted_splits = [name];
    redrawBins();
}

function appendResult() {
    // check if tree exists
    if ($.isEmptyObject(label_to_node_map)) {
        alert('Draw tree first.');
        return;
    }

    var bin_id = getBinId();

    if (bin_id === 'undefined')
        return;

    var bins_to_update = [];
    var _len = search_results.length;
    for (var i=0; i < _len; i++) {
        _contig_name = search_results[i]['split'];
        if (SELECTED[bin_id].indexOf(_contig_name) == -1) {
            SELECTED[bin_id].push(_contig_name);

            if (bins_to_update.indexOf(bin_id) == -1)
                bins_to_update.push(bin_id);
        }

        for (var bid = 1; bid <= bin_counter; bid++) {
            // don't remove nodes from current bin
            if (bid == bin_id)
                continue;

            var pos = SELECTED[bid].indexOf(_contig_name);
            if (pos > -1) {
                SELECTED[bid].splice(pos, 1);

                if (bins_to_update.indexOf(bid) == -1)
                    bins_to_update.push(bid);
            }
        }
    }

    updateBinsWindow(bins_to_update);
    redrawBins();
}

function removeResult() {
    // check if tree exists
    if ($.isEmptyObject(label_to_node_map)) {
        alert('Draw tree first.');
        return;
    }

    var bin_id = getBinId();

    if (bin_id === 'undefined')
        return;

    var bins_to_update = [];
    var _len = search_results.length;
    for (var i=0; i < _len; i++) {
        _contig_name = search_results[i]['split'];

        var pos = SELECTED[bin_id].indexOf(_contig_name);
        if (pos > -1) {
            SELECTED[bin_id].splice(pos, 1);
            
            if (bins_to_update.indexOf(bin_id) == -1)
                bins_to_update.push(bin_id);
        }
    }

    updateBinsWindow(bins_to_update);
    redrawBins();
}
