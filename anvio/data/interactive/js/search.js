
function searchContigs() 
{
    var svalue = $('#searchValue').val();

    if (svalue == "")
    {
        alert("Search value shouldn't be empty.");
        return;
    }
    var column = $('#searchLayerList').val();
    search_column = column;
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
            search_results.push(row);
            _counter++;
        }
    }
    $('#search_result_message').html(_counter + " splits found.");
}

function showSearchResult() {
    var clear_link = '<a href="#" onclick="$(\'.search-results-display, #search-results-table-search-item, #search-results-table-search-name, #search-results-table-header\').html(\'\');">(clear)</a>';
    $("#search-results-table-header").html('<h4>Search results ' + clear_link + '</h4>');
    $("#search-results-table-search-name").html('Split name');
    $("#search-results-table-search-item").html(layerdata[0][search_column]);

    var rows = "";
    var _len = search_results.length;
    for (var i=0; i < _len; i++)
    {
        rows = rows + "<tr><td data-value=" + layerdata[search_results[i]][0] + "><a href='#' class='no-link' onclick='highlightSplit(\"" + layerdata[search_results[i]][0] + "\");'>" + layerdata[search_results[i]][0] + "</a></td><td data-value=" + layerdata[search_results[i]][search_column] + ">" + layerdata[search_results[i]][search_column] + "</td></tr>";
    }
    $(".search-results-display").html(rows);
}

function highlightResult() {
    // check if tree exists
    if ($.isEmptyObject(label_to_node_map)) {
        alert('Draw tree first.');
        return;
    }

    highlighted_splits = [];

    for (var i=0; i < search_results.length; i++) {
        var _contig_name = layerdata[search_results[i]][0];
        
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
        _contig_name = layerdata[search_results[i]][0];
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
        _contig_name = layerdata[search_results[i]][0];
        var _id = label_to_node_map[_contig_name].id;

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
