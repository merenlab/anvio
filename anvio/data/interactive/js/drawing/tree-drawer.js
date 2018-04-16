var TreeDrawer = function(tree, container, settings) {
    this.tree = tree;
    this.settings = settings;
    this.container = container;

    createBin('svg', 'viewport');
    createBin('viewport', 'tree');

    if (this.settings['tree-type'] == 'circlephylogram') {
        this.drawer = new CirclePhylogramTreeDrawer(this.tree, 'tree');
    } else {
        this.drawer = new PhylogramTreeDrawer(this.tree, 'tree');
    }
};


TreeDrawer.prototype.Draw = function() {
    this.drawer.CalculateLeafSizes();

    for (const p of this.tree.root.IterateChildren()) {
        if (p.IsLeaf()) {
            this.drawer.CalculateLeafCoordinate(p);
        } else {
            this.drawer.CalculateInternalNodeCoordinate(p);
        }
    }

    if (this.tree.rooted) {
        this.drawer.DrawRoot();
    }
                
    for (const p of this.tree.root.IterateChildren()) {
        if (p.IsLeaf()) {
            this.drawer.DrawLeaf(p);
        } else {
            this.drawer.DrawInternalNode(p);
        }
    }
};

var PhylogramTreeDrawer = function(tree, container) {
    this.tree = tree;
    this.container = container;
    this.height = 500;
    this.width = 500;

    this.root_length = 0.1 * this.height;
};


PhylogramTreeDrawer.prototype.CalculateLeafSizes = function() {
    let sum_size = 0;

    for (const p of this.tree.leaves) {
        sum_size += p.size;
    }

    let size_unit = this.width / sum_size;

    for (const p of this.tree.leaves) {
        p.size = p.size * size_unit;
    }
}

PhylogramTreeDrawer.prototype.CalculateLeafCoordinate = function(p) {
    var pt = [];

    if (p.order == 0) {
        pt['x'] = p.size / 2;
    } else {
        var prev_leaf = this.tree.leaves[p.order - 1];
        pt['x'] = prev_leaf.xy['x'] + prev_leaf.size / 2 + p.size / 2;
    }

    if (this.tree.has_edge_lengths) {
        pt['y'] = this.root_length + (p.path_length / this.tree.max_path_length) * (this.height - this.root_length);
    } else {
        pt['y'] = this.root_length + ((this.tree.root.depth - p.depth) / this.tree.root.depth) * (this.height - this.root_length);
    }

    p.xy['y'] = pt['y'];
    p.xy['x'] = pt['x'];
};


PhylogramTreeDrawer.prototype.CalculateInternalNodeCoordinate = function(p) {
    var pt = [];

    if (this.tree.has_edge_lengths) {
        pt['y'] = this.root_length + (p.path_length / this.tree.max_path_length) * (this.height - this.root_length);
    } else {
        pt['y'] = this.root_length + ((this.tree.root.depth - p.depth) / this.tree.root.depth) * (this.height - this.root_length);
    }

    var pl = p.child.xy;
    var pr = p.child.GetRightMostSibling().xy;

    pt['x'] = pl['x'] + (pr['x'] - pl['x']) / 2;
    p.xy['y'] = pt['y'];
    p.xy['x'] = pt['x'];
};


PhylogramTreeDrawer.prototype.DrawRoot = function() {
    var p0 = this.tree.root.xy
    var p1 = {};

    p1['x'] = p0['x'];
    p1['y'] = 0;

    drawLine(this.container, this.tree.root, p0, p1);
};


PhylogramTreeDrawer.prototype.DrawLeaf = function(p) {
    var p0 = p.xy
    var p1 = [];
    var anc = p.ancestor;
    if (anc) {
        p1['y'] = anc.xy['y'];
        p1['x'] = p0['x'];

        drawLine(this.container, p, p0, p1);
    }
};


PhylogramTreeDrawer.prototype.DrawInternalNode = function(p) {
    var p0 = [];
    var p1 = [];

    p0['y'] = p.xy['y'];
    p0['x'] = p.xy['x'];

    var anc = p.ancestor;
    if (anc) {
        p1['y'] = anc.xy['y'];
        p1['x'] = p0['x'];

        drawLine(this.container, p, p0, p1);
    }

    var pl = p.child.xy;
    var pr = p.child.GetRightMostSibling().xy;

    p0['y'] = p0['y'];
    p0['x'] = pl['x'];
    p1['y'] = p0['y'];
    p1['x'] = pr['x'];

    drawLine(this.container, p, p0, p1, true);
};


var CirclePhylogramTreeDrawer = function(tree, container) {
    this.tree = tree;
    this.container = container;

    this.radius = 3000;
    this.root_length = 0.1 * (this.radius / 2);

    this.angle_min = Math.toRadians(0);
    this.angle_max = Math.toRadians(270);
};


CirclePhylogramTreeDrawer.prototype.CalculateLeafSizes = function() {
    let sum_size = 0;

    for (const p of this.tree.leaves) {
        sum_size += p.size;
    }

    let size_unit = Math.abs(this.angle_max - this.angle_min) / sum_size;

    for (const p of this.tree.leaves) {
        p.size = p.size * size_unit;
    }
};

CirclePhylogramTreeDrawer.prototype.CalculateLeafCoordinate = function(p) {
    if (p.order == 0) {
        p.angle = this.angle_min + p.size / 2;
    } else {
        let prev_leaf = this.tree.leaves[p.order - 1];
        p.angle = prev_leaf.angle + prev_leaf.size / 2 + p.size / 2;
    }

    if (this.tree.has_edge_lengths) {
        p.radius = this.root_length + (p.path_length / this.tree.max_path_length) * ((this.radius / 2) - this.root_length);
    }
    else
    {
        p.radius = this.root_length + ((this.tree.root.depth - p.depth) / this.tree.root.depth) * ((this.radius / 2) - this.root_length);
    }

    var pt = [];
    pt['x'] = p.radius * Math.cos(p.angle);
    pt['y'] = p.radius * Math.sin(p.angle);

    p.xy['x'] = pt['x'];
    p.xy['y'] = pt['y'];
};


CirclePhylogramTreeDrawer.prototype.CalculateInternalNodeCoordinate = function(p) {
    var left_angle = p.child.angle;
    var right_angle = p.child.GetRightMostSibling().angle;

    p.angle = left_angle + (right_angle - left_angle) / 2;
    
    if (this.tree.has_edge_lengths) {
        p.radius = this.root_length + (p.path_length / this.tree.max_path_length) * ((this.radius / 2) - this.root_length);
    }
    else
    {
        p.radius = this.root_length + ((this.tree.root.depth - p.depth) / this.tree.root.depth) * ((this.radius / 2) - this.root_length);
    }

    var pt = [];
    pt['x'] = p.radius * Math.cos(p.angle);
    pt['y'] = p.radius * Math.sin(p.angle);

    p.xy['x'] = pt['x'];
    p.xy['y'] = pt['y'];

    var q = p.child;
    while (q) {
        pt = [];

        pt['x'] = p.radius * Math.cos(q.angle);
        pt['y'] = p.radius * Math.sin(q.angle);

        q.backarc = [];
        q.backarc['x'] = pt['x'];
        q.backarc['y'] = pt['y'];

        q = q.sibling;
    }
};


CirclePhylogramTreeDrawer.prototype.DrawRoot = function() {
    var p0 = this.tree.root.xy
    var p1 = {};

    p1['x'] = 0;
    p1['y'] = 0;

    drawLine(this.container, this.tree.root, p0, p1);
}


CirclePhylogramTreeDrawer.prototype.DrawLeaf = function(p) {
    var p0 = p.xy
    var p1 = p.backarc;
    drawLine(this.container, p, p0, p1);
};


CirclePhylogramTreeDrawer.prototype.DrawInternalNode = function(p) {
    var p0 = [];
    var p1 = [];

    p0['x'] = p.xy['x'];
    p0['y'] = p.xy['y'];

    var anc = p.ancestor;
    if (anc) {
        p0 = p.xy;
        p1 = p.backarc;

        drawLine(this.container, p, p0, p1);
    }
    p0 = p.child.backarc;
    p1 = p.child.GetRightMostSibling().backarc;

    var large_arc_flag = (Math.abs(p.child.GetRightMostSibling().angle - p.child.angle) > Math.PI) ? true : false;
    drawCircleArc(this.container, p, p0, p1, p.radius, large_arc_flag);
};

