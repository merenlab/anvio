/**
 * Contigs db stats visualization
 *
 *  Authors: Ozcan Esen
 *
 * Copyright 2015-2021, The anvi'o project (http://anvio.org)
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

function draw_n_values_plot(container, stats) {
    var svg = d3.select(container)
        .append('svg')
        .attr('viewBox', '0 50 800 260')

    var g = svg.append("g")
        .attr("transform", "translate(0, 100)");

    var plot_width = 780;
    var plot_height = 180;

    var xscale = d3.scale.ordinal().rangeBands([120, plot_width]);
    var yscale = d3.scale.sqrt().rangeRound([plot_height, 0]);
    var color_scale = d3.scale.sqrt().range(['#fff7bc', '#d95f0e']);

    var data = stats['n_values'];
    xscale.domain(data.map(function(d, i) { return i; }));
    yscale.domain([0, d3.max(data, function(d) { return d.length; })]);
    color_scale.domain([0, d3.max(data, function(d) { return d.length; })]);

    var bar_group = g.selectAll(".bar")
                        .data(data)
                        .enter()
                        .append("g")
                        .on('mouseenter', function(d, i) {
                            d3.select(this)
                                .select('rect')
                                .attr('fill-opacity', '0.5');

                            var tooltip_pos = {
                                'x': Math.max(Math.min(xscale(i) - 100,  600), 120),
                                'y': yscale(d.length) - 42,
                                'width': 200,
                                'height': 40
                            };

                            var tooltip = d3.select(this)
                                .append('g')
                                    .attr('class', 'plot-tooltip')
                                    .attr('transform', 'translate(' + tooltip_pos.x + ',' + tooltip_pos.y + ')');

                            tooltip.append('rect')
                                    .attr('rx', '8')
                                    .attr('ry', '8')
                                    .attr('height', tooltip_pos.height)
                                    .attr('fill-opacity', '0.7')
                                    .attr('width', tooltip_pos.width);

                            var tooltip_text = tooltip.append('text')
                                                .attr('x', '10')
                                                .attr('y', '15')
                                                .attr('fill', '#FFFFFF')
                                                .attr('font-family', 'Montserrat', 'Helvetica')
                                                .attr('font-size', '12px');

                            tooltip_text.append('tspan')
                                            .text('L' + (i+1) + ': ' + d.num_contigs + ' contigs');

                            tooltip_text.append('tspan')
                                            .text('N' + (i+1) + ': ' + getCommafiedNumberString(d.length) + ' nts')
                                            .attr('x', '10')
                                            .attr('dy', '1.4em');

                        })
                        .on('mouseleave', function(d) {
                            d3.select(this)
                                .select('rect')
                                .attr('fill-opacity', '1');

                            d3.select(this)
                                .select('.plot-tooltip').remove();
                        });

    bar_group.append('rect')
            .attr("class", "bar")
            .attr("x", function(d,i) { return xscale(i); })
            .attr("y", function(d) { return yscale(d.length); })
            .attr("width", xscale.rangeBand())
            .attr("fill", function(d,i) { return (i == 49) ? '#DC143C' : color_scale(d.length); })
            .attr("height", function(d) { return plot_height - yscale(d.length); });

    g.append('text')
            .text('L50')
            .style("font-size", '9.5')
            .attr("x", xscale(49) - 2)
            .attr("y", yscale(0) + 10);

    g.append('text')
            .text('L0')
            .style("font-size", '9.5')
            .attr("x", xscale(0) - 2)
            .attr("y", yscale(0) + 10);

    g.append('text')
            .text('L100')
            .style("font-size", '9.5')
            .attr("x", xscale(99) - 2)
            .attr("y", yscale(0) + 10);

    var yAxis = d3.svg.axis()
        .scale(yscale)
        .orient("left")
        .tickFormat(function(d) { return getReadableSeqSizeString(d); })
        .tickSubdivide(0);

    g.append("g")
        .attr("class", "x_axis")
        .attr("transform", 'translate(118,0)')
        .style("font-size", '10')
        .call(yAxis);
};

function draw_gene_counts_chart(container, gene_counts) {
        var svg = d3.select(container)
            .append('svg')
            .attr('viewBox', '0 20 800 300')

        var g = svg.append("g")
            .attr("transform", "translate(0, 100)");

        var plot_width = 780;
        var plot_height = 180;

        var data = [];
        for (key in gene_counts) {
            data.push({'name': key, 'value': gene_counts[key]});
        }
        data.sort(function(a, b) { return (a['value'] < b['value']) - (a['value'] > b['value'])});

        var xscale = d3.scale.ordinal().rangeBands([120, plot_width]);
        var yscale = d3.scale.linear().rangeRound([plot_height, 0]);
        var color_scale = d3.scale.linear().range(['#fff7bc', '#d95f0e']);

        xscale.domain(data.map(function(d) { return d.name; }));
        yscale.domain([0, d3.max(data, function(d) { return d.value; })]);
        color_scale.domain([0, d3.max(data, function(d) { return d.value; })]);

        var bar_group = g.selectAll(".bar")
                            .data(data)
                            .enter()
                            .append("g")
                            .on('mouseenter', function(d) {
                                d3.select(this)
                                    .select('rect')
                                    .attr('fill-opacity', '0.5');

                                var tooltip_pos = {
                                    'x': Math.max(Math.min(xscale(d.name) - 80,  620), 120),
                                    'y': yscale(d.value) - 50,
                                    'width': 160,
                                    'height': 40
                                };

                                var tooltip = d3.select(this)
                                    .append('g')
                                        .attr('class', 'plot-tooltip')
                                        .attr('transform', 'translate(' + tooltip_pos.x + ',' + tooltip_pos.y + ')');

                                tooltip.append('rect')
                                        .attr('rx', '8')
                                        .attr('ry', '8')
                                        .attr('height', tooltip_pos.height)
                                        .attr('fill-opacity', '0.7')
                                        .attr('width', tooltip_pos.width);

                                var tooltip_text = tooltip.append('text')
                                                    .attr('x', '10')
                                                    .attr('y', '15')
                                                    .attr('fill', '#FFFFFF')
                                                    .attr('font-family', 'Montserrat', 'Helvetica')
                                                    .attr('font-size', '12px');

                                tooltip_text.append('tspan')
                                                .text('Gene: ' + d.name);

                                tooltip_text.append('tspan')
                                                .text('Count: ' + d.value)
                                                .attr('x', '10')
                                                .attr('dy', '1.4em');

                            })
                            .on('mouseleave', function(d) {
                                d3.select(this)
                                    .select('rect')
                                    .attr('fill-opacity', '1');

                                d3.select(this)
                                    .select('.plot-tooltip').remove();
                            });

    bar_group.append('rect')
            .attr("class", "bar")
            .attr("x", function(d) { return xscale(d.name); })
            .attr("y", function(d) { return yscale(d.value); })
            .attr("width", xscale.rangeBand())
            .attr("fill", function(d) { return color_scale(d.value); })
            .attr("height", function(d) { return plot_height - yscale(d.value); });

    var yAxis = d3.svg.axis()
        .scale(yscale)
        .orient("left")
        .tickFormat(d3.format("d"))
        .tickSubdivide(0);

    g.append("g")
        .attr("class", "x_axis")
        .attr("transform", 'translate(118,0)')
        .call(yAxis);

    var bins = d3.layout.histogram()
        .bins(data[0].value)
        (data.map(function(x) { return x.value }));

    var bins_y = d3.scale.linear()
        .domain(d3.extent(bins, function(d){return d.length }))
        .range([0, 70]);

    var bins_y_reversed = d3.scale.linear()
        .domain(d3.extent(bins, function(d){return d.length }).reverse())
        .range([0, 70]);

    var histogram = g.append('g').attr('class', 'histogram');

    histogram.selectAll('rect')
            .data(bins)
            .enter()
            .append('rect')
            .attr("class", "histogram-bar")
            .attr("y", function(d) { return yscale(d.x); })
            .attr("height", function(d) { return 1; })
            .attr("x", function(d) { return 80 - bins_y(d.y); })
            .attr("width", function(d) { return bins_y(d.y); } )
            .attr("fill", '#d95f0e');

    var axis_for_histogram = d3.svg.axis()
        .scale(bins_y_reversed)
        .orient("top")
        .ticks(5)
        .tickFormat(d3.format("d"))
        .tickSubdivide(0);

    g.append("g")
        .attr("class", "histogram_axis")
        .attr("transform", 'translate(10,0)')
        .style("font-size", '8')
        .call(axis_for_histogram);
};
