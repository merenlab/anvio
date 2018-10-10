
function searchContigs() 
{
    var svalue = $('#searchValue').val();

    if (svalue == "")
    {
        alert("Search value shouldn't be empty.");
        return;
    }
    var column = $('#searchLayerList').val();
    search_column = (column == 0) ? 'Item Name' : layerdata[0][column];
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
        var _post = "][" + column + "].toString().toLowerCase().indexOf(\"" + svalue + "\".toLowerCase()) != -1";
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
    $('#search_result_message').html(_counter + " result(s) found.");
}

function searchFunctions() {
    $('#search_functions_button').prop( "disabled", true );
    $('#search_result_message_functions').html('<img src="images/loading.gif" style="width: 16px;" />');

    let requested_sources = [];
    $('#functions_sources_list input:checkbox:checked').each((i, v) => {
        requested_sources.push($(v).val());
    });

    $.ajax({
        type: 'POST',
        cache: false,
        url: '/data/search_functions',
        data: {
            terms: $('#searchFunctionsValue').val(),
            sources: requested_sources
        },
        success: function(data) {
            if (data['status'] == 0) {
                $('.search-message').hide();
                $('.search-message').html('');

                search_results = [];
                search_column = 'Annotation';

                for (var i=0; i < data['results'].length; i++) {
                    console.log(mode);
                    if (mode == 'pan') {
                        var _gene_caller_id = data['results'][i][0];
                        var _genome_name    = data['results'][i][1];
                        var _source         = data['results'][i][2];
                        var _accession      = data['results'][i][3];
                        var _annotation     = data['results'][i][4];
                        var _search_term    = data['results'][i][5];
                        var _split_name     = data['results'][i][6];                        
                    }
                    else
                    {
                        var _gene_caller_id = data['results'][i][0];
                        var _source         = data['results'][i][1];
                        var _accession      = data['results'][i][2];
                        var _annotation     = data['results'][i][3];
                        var _search_term    = data['results'][i][4];
                        var _split_name     = data['results'][i][5];
                    }

                    var _beginning = _annotation.toLowerCase().indexOf(_search_term.toLowerCase());
                    _annotation = [_annotation.slice(0, _beginning), 
                                   '<mark>', 
                                   _annotation.slice(_beginning, _beginning + _search_term.length), 
                                   '</mark>', 
                                   _annotation.slice(_beginning + _search_term.length, _annotation.length)
                                   ].join("");

                    var reported_column = _split_name;
                    if (mode == 'gene') {
                        reported_column = _gene_caller_id;
                    }
                    search_results.push({'split': reported_column , 'value': '<b>Gene caller id:</b> ' + _gene_caller_id +
                                                                                  '</br><b>Source:</b> ' + _source +
                                                                                  '</br><b>Accession:</b> ' + _accession +
                                                                                  '</br><b>Annotation:</b> ' + _annotation});
                }
                $('#search_result_message_functions').html(data['results'].length + " result(s) found in " + data['item_count'] + ' unique item(s).');
            } else {
                $('.search-message').show();
                $('.search-message').html(data['message']);
                $('#search_result_message_functions').html('');
            }

            $('#search_functions_button').prop( "disabled", false );
        },
        done: function() {
        }
    });
}

function filterGeneClusters() {
    var parameters = {};

    $('.pan-filters input:text').each(function (index, input){
        if (!$(input).prop('disabled')) {
            parameters[$(input).attr('parameter')] = $(input).val();
        }
    });

    if (Object.keys(parameters).length == 0) {
        $('.pan-filter-error').show();
        $('.pan-filter-error').html("You need to select at least one filter.");
        $('#search_result_message_pan_filter').html('');
        return;
    }

    $.ajax({
        type: 'POST',
        cache: false,
        url: '/data/filter_gene_clusters',
        data: parameters,
        success: function(data) {
            if (data['status'] == 0) {
                $('.pan-filter-error').hide();
                $('.pan-filter-error').html('');

                search_results = [];
                search_column = '';

                for (var i=0; i < data['gene_clusters_list'].length; i++) {
                    search_results.push({'split': data['gene_clusters_list'][i], 'value': ''});
                }
                $('#search_result_message_pan_filter').html(data['gene_clusters_list'].length + " gene clusters passed the filter.");
            } else {
                $('.pan-filter-error').show();
                $('.pan-filter-error').html(data['message']);
                $('#search_result_message_pan_filter').html('');
            };
        }
    });
}

function showSearchResult() {
    var clear_link = '<a href="#" onclick="$(\'.search-results-display, #search-results-table-search-item, #search-results-table-search-name, #search-results-table-header\').html(\'\');">(clear)</a>';
    $("#search-results-table-header").html('<h4>Search results ' + clear_link + '</h4>');
    $("#search-results-table-search-name").html('Item Name');
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
    if (!drawer) {
        alert('Draw tree first.');
        return;
    }

    highlighted_splits = [];
    for (var i=0; i < search_results.length; i++) {
        highlighted_splits.push(search_results[i]['split']);
    }

    bins.HighlightItems(highlighted_splits); 
}

function highlightSplit(name) {
    bins.HighlightItems(name); 
}

function appendResult() {
    let node_list = [];
    for (const result of search_results) {
        let node = drawer.tree.GetLeafByName(result['split']);
        node_list.push(node);
    }

    bins.AppendNode(node_list);
}

function removeResult() {
    let node_list = [];
    for (const result of search_results) {
        let node = drawer.tree.GetLeafByName(result['split']);
        node_list.push(node);
    }

    bins.RemoveNode(node_list);
}
