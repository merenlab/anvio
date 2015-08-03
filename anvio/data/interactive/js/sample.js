sampleOrganizationResponse = [JSON.parse('{"num_reads":{"order":"","newick":"((DAY_19:3.8927e-12,(DAY_15A:1.17876e-12,DAY_17A:1.17876e-12)Int16:3.8927e-12)Int18:1.58727e-11,((DAY_15B:1.14795e-12,(DAY_17B:4.17546e-14,DAY_23:4.17546e-14)Int12:1.14795e-12)Int15:5.55138e-12,(DAY_22B:1.20141e-12,((DAY_18:3.08426e-14,DAY_22A:3.08426e-14)Int11:7.76164e-13,(DAY_16:3.99589e-13,DAY_24:3.99589e-13)Int13:7.76164e-13)Int14:1.20141e-12)Int17:5.55138e-12)Int19:1.58727e-11);"},"even_odd":{"order":"","newick":"((DAY_22B:0.000332086,((DAY_18:8.54156e-06,DAY_22A:8.54156e-06)Int11:0.00021462,(DAY_16:0.000110463,DAY_24:0.000110463)Int13:0.00021462)Int14:0.000332086)Int15:0.0376259,((DAY_15B:0.000506205,(DAY_17B:1.8431e-05,DAY_23:1.8431e-05)Int12:0.000506205)Int16:0.00441781,(DAY_19:0.00170499,(DAY_15A:0.00051519,DAY_17A:0.00051519)Int17:0.00170499)Int18:0.00441781)Int19:0.0376259);"},"basic":{"order":"DAY_15A,DAY_15B,DAY_16,DAY_17A,DAY_17B,DAY_18,DAY_19,DAY_22A,DAY_22B,DAY_23,DAY_24","newick":""}, "mini_test": {"order":"s204_6M,s204_7M,s204_9M", "newick":""}}')];


function get_newick_leaf_order(newick)
{
    var _order_list = [];

    var t = new Tree();
    t.Parse(newick);

    var n = new NodeIterator(t.root);
    var q = n.Begin();

    order_counter = 0;
    while (q != null)
    {
        if (q.IsLeaf()) {
            _order_list.push(q.label);
        }
        q=n.Next();
    }

    return _order_list;
}

$(document).ready(function() {
    $('#sample_organization').change(function() {
        if (this.value == 'none') 
            return;

        var organization = sampleOrganizationResponse[0][this.value];
        var sample_order;

        // get new sample order
        if (organization['order'] != "")
        {
            sample_order = organization['order'].split(',');
        }
        else
        {
            sample_order = get_newick_leaf_order(organization['newick']);
        }
        $.map(sample_order, $.trim);

        for(var i=0; i < sample_order.length; i++)
        {
            var layer_id = getLayerId(sample_order[i]);
            var detached_row = $('#height' + layer_id).closest('tr').detach();
            $('#tbody_layers').append(detached_row);
        }
    });
});


function drawSampleLayers() {

}
