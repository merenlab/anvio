var VIEWER_WIDTH = window.innerWidth || document.documentElement.clientWidth || document.getElementsByTagName('body')[0].clientWidth;
var VIEWER_HEIGHT = window.innerHeight || document.documentElement.clientHeight || document.getElementsByTagName('body')[0].clientHeight;

var layers;
var coverage;

function getUrlVars() {
	var map = {};
	var parts = window.location.href.replace(/[?&]+([^=&]+)=([^&]*)/gi, function(m,key,value) {
		map[key] = value;
	});
	return map;
}

function loadAll() {
    contig_id = getUrlVars()["contig"];

    document.getElementById("header").innerHTML = "<strong>" + contig_id + "</strong> detailed";

    $.ajax({
        type: 'GET',
        cache: false,
        url: '/data/charts/' + contig_id + '?timestamp=' + new Date().getTime(),
 
        /* function prototype for interactive-binning::bottle

           @route('/data/charts/<contig_name>')
           def charts(contig_name):
               print contig_name
               data = {'layers': ['contig_1', 'contig_2', ..., 'contig_N'],
                       'coverage': 
                          [c_1, c_2, ..., c_n]}
               return json.dumps(data) */

        success: function(data) {
            d = JSON.parse(data);
            layers = d.layers;
            coverage = d.coverage;

            createCharts();
        }
    });
}


function createCharts(){
    /* Adapted from Tyler Craft's Multiple area charts with D3.js article:
    http://tympanus.net/codrops/2012/08/29/multiple-area-charts-with-d3-js/  */

    var margin = {top: 20, right: 10, bottom: 150, left: 50};
    var width = VIEWER_WIDTH * .85;
    var height = (150 * (layers.length + 1)) - margin.top - margin.bottom;
    var contextHeight = 50;
    var contextWidth = width;

    var svg = d3.select("#chart-container").append("svg")
            .attr("width", width + margin.left + margin.right)
            .attr("height", (height + margin.top + margin.bottom));

    $('#chart-container').css("width", (width + 100) + "px");
    
    var charts = [];
    var maxDataPoint = 0;
    
    var layersCount = layers.length;
    var chartHeight = height * (1 / layersCount);
    
    coverage.forEach(function(d) {
        for (var prop in d) {
            if (d.hasOwnProperty(prop)) {
                d[prop] = parseFloat(d[prop]);
                
                if (d[prop] > maxDataPoint) {
                    maxDataPoint = d[prop];
                }
            }
        }
        
    });

    for(var i = 0; i < layersCount; i++){
        charts.push(new Chart({
                        name: layers[i],
                        coverage: coverage[i],
                        id: i,
                        width: width,
                        height: height * (1 / layersCount),
                        maxDataPoint: maxDataPoint,
                        svg: svg,
                        margin: margin,
                        showBottomAxis: (i == layers.length - 1)
                }));
        
    }
    
    /* Let's create the context brush that will 
            let us zoom and pan the chart */
    var contextXScale = d3.scale.linear().range([0, contextWidth]).domain(charts[0].xScale.domain());    
    
    var contextAxis = d3.svg.axis()
                            .scale(contextXScale)
                            .tickSize(contextHeight)
                            .orient("bottom");

    var contextArea = d3.svg.area()
                            .interpolate("monotone")
                            .x(function(d) { return contextXScale(d); })
                            .y0(contextHeight)
                            .y1(0);

    var brush = d3.svg.brush()
                                        .x(contextXScale)
                                        .on("brush", onBrush);

    var context = svg.append("g")
                .attr("class","context")
                .attr("transform", "translate(" + (margin.left) + "," + (height + margin.top * 3) + ")");
    
    context.append("g")
                .attr("class", "x axis top")
                .attr("transform", "translate(0,0)")
                .call(contextAxis)
                                        
    context.append("g")
                .attr("class", "x brush")
                .call(brush)
                .selectAll("rect")
                .attr("y", 0)
                .attr("height", contextHeight);
    
    function onBrush(){
        /* this will return a date range to pass into the chart object */
        var b = brush.empty() ? contextXScale.domain() : brush.extent();
        b = [Math.floor(b[0]), Math.floor(b[1])];
        for(var i = 0; i < layersCount; i++){
            charts[i].showOnly(b);
        }
    }
}

function Chart(options){
    this.chartData = options.coverage;
    this.width = options.width;
    this.height = options.height;
    this.maxDataPoint = options.maxDataPoint;
    this.svg = options.svg;
    this.id = options.id;
    this.name = options.name;
    this.margin = options.margin;
    this.showBottomAxis = options.showBottomAxis;
    
    var localName = this.name;
    
    this.xScale = d3.scale.linear()
                            .range([0, this.width])
                            .domain([0, this.chartData.length]);
    
    this.yScale = d3.scale.linear()
                            .range([this.height,0])
                            .domain([0,this.maxDataPoint]);
    var xS = this.xScale;
    var yS = this.yScale;
    
    /* 
        This is what creates the chart.
        There are a number of interpolation options. 
        'basis' smooths it the most, however, when working with a lot of data, this will slow it down 
    */
    this.area = d3.svg.area()
                            .x(function(d, i) { return xS(i); })
                            .y0(this.height)
                            .y1(function(d) { return yS(d); });
    /*
        This isn't required - it simply creates a mask. If this wasn't here,
        when we zoom/panned, we'd see the chart go off to the left under the y-axis 
    */
    this.svg.append("defs").append("clipPath")
                                    .attr("id", "clip-" + this.id)
                                    .append("rect")
                                    .attr("width", this.width)
                                    .attr("height", this.height);
    /*
        Assign it a class so we can assign a fill color
        And position it on the page
    */
    this.chartContainer = this.svg.append("g")
                        .attr('class',this.name.toLowerCase())
                        .attr("transform", "translate(" + this.margin.left + "," + (this.margin.top + (this.height * this.id) + (10 * this.id)) + ")");

    /* We've created everything, let's actually add it to the page */
    this.chartContainer.append("path")
                                            .data([this.chartData])
                                            .attr("class", "chart")
                                            .attr("clip-path", "url(#clip-" + this.id + ")")
                                            .attr("d", this.area);
                                    
    this.xAxisTop = d3.svg.axis().scale(this.xScale).orient("top");

    if(this.id == 0){
        this.chartContainer.append("g")
                    .attr("class", "x axis top")
                    .attr("transform", "translate(0,0)")
                    .call(this.xAxisTop);
    }
    
        
    this.yAxis = d3.svg.axis().scale(this.yScale).orient("left").ticks(5);
        
    this.chartContainer.append("g")
                                            .attr("class", "y axis")
                                            .attr("transform", "translate(-15,0)")
                                            .call(this.yAxis);
                                            
    this.chartContainer.append("text")
                                            .attr("class","country-title")
                                            .attr("transform", "translate(0,20)")
                                            .text(this.name);
    
}
    
Chart.prototype.showOnly = function(b){
        this.xScale.domain(b);
        this.chartContainer.select("path").data([this.chartData]).attr("d", this.area);
        this.chartContainer.select(".x.axis.top").call(this.xAxisTop);
}

