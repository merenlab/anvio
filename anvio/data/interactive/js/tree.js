/**
 * Javascript library to parse newick trees
 *
 *  Authors: Özcan Esen <ozcanesen@gmail.com>
 *           Matthew Klein <mtt.l.kln@gmail.com>
 *           A. Murat Eren <a.murat.eren@gmail.com>
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


function Node(label) {
    this.ancestor = null;
    this.child = null;
    this.sibling = null;
    this.collapsed = false;
    this.label = label;
    this.size = 1;
    this.id = 0;
    this.xy = [];
    this.original_edge_length = 0.0;
    this.edge_length = 0.0;
    this.path_length = 0.0;
    this.depth = 0;
    this.order = null;
    this.max_child_path = 0;
    this.branch_support = "";
}


Node.prototype.IsLeaf = function() {
    return (!this.child);
};


Node.prototype.GetRightMostSibling = function() {
    var p = this;
    while (p.sibling) {
        p = p.sibling;
    }
    return p;
};


Node.prototype.GetBorderNodes = function() {
    var p1 = this;
    while (p1.child) {
        p1 = p1.child;
    }

    var p2 = this;
    while (p2.child) {
        p2 = p2.child.GetRightMostSibling();
    }

    return [p1, p2];
};


Node.prototype.GetChildren = function() {
    var n = new NodeIterator(this);
    var q = n.Begin();
    var children = [];

    while (q != null)
    {
        children.push(q);
        q=n.Next();
    }

    return children;
}

Node.prototype.GetAncestors = function() {
    let q = this;
    let ancestors = [];

    while (q)
    {
        ancestors.push(q);
        q = q.ancestor;
    }

    return ancestors;
}


Node.prototype.IterateChildren = function*() {
    var n = new NodeIterator(this);
    var q = n.Begin();

    while (q != null)
    {
        yield q;
        q=n.Next();
    }
}

Node.prototype.SetColor = function(color, width=3) {
    let line = document.getElementById('line' + this.id);
    if (line) {
        line.style['stroke-width'] = width;
        line.style['stroke'] = color;
    }

    let arc = document.getElementById('arc' + this.id);
    if (arc) {
        arc.style['stroke-width'] = width;
        arc.style['stroke'] = color;
    }
}


Node.prototype.ResetColor = function() {
    this.SetColor(LINE_COLOR, width=1);
}


Node.prototype.Rotate = function() {
    if (this.child) {
        var siblings = [];
        var p = this.child;

        while (p) {
            p.Rotate();
            siblings.push(p);
            p = p.sibling;
        }

        p = siblings.pop();
        this.child = p;

        while (siblings.length > 0) {
            p.sibling = siblings.pop()
        }

        p.sibling.sibling = null;
    }
};


function Tree() {
    this.root = null;
    this.num_leaves = 0;
    this.num_nodes = 0;
    this.nodes = [];
    this.leaves = [];
    this.rooted = true;
    this.has_edge_lengths = false;
    this.has_branch_supports = false;
    this.error = 0;
    this.label_to_leaves = {};
}


Tree.prototype.GetLeafByName = function(name) {
    if (this.label_to_leaves.hasOwnProperty(name)) {
        return this.label_to_leaves[name];
    }
    return null;
}


Tree.prototype.NewNode = function() {
    var node = new Node();
    node.id = this.num_nodes++;
    this.nodes[node.id] = node;
    return node;
};


Tree.prototype.Parse = function(str, edge_length_norm) {
    str = str.replace(/\(/g, "|(|");
    str = str.replace(/\)/g, "|)|");
    str = str.replace(/,/g, "|,|");
    str = str.replace(/:/g, "|:|");
    str = str.replace(/;/g, "|;|");
    str = str.replace(/\|\|/g, "|");
    str = str.replace(/^\|/, "");
    str = str.replace(/\|$/, "");

    var token = str.split("|");
    var curnode = this.NewNode();
    this.root = curnode;

    var state = 0;
    var stack = [];
    var i = 0;
    var q = null;

    this.error = 0;

    while ((state != 99) && (this.error == 0)) {
        switch (state) {
            case 0:
                if (ctype_alnum(token[i].charAt(0)) || token[i].charAt(0) == "'" || token[i].charAt(0) == '"' || token[i].charAt(0) == '_') {
                    if (isNumber(token[i]) && mode != 'gene' || isNaN(token[i] && mode != 'gene')){
                        curnode.branch_support = token[i];
                        this.has_branch_supports = true;
                    } 
                    else {
                        this.label_to_leaves[token[i]] = curnode;
                        curnode.label = token[i];
                    }
                    i++;
                    state = 1;
                } else {
                    switch (token[i]) {
                        case '(':
                            state = 2;
                            break;

                        default:
                            state = 99;
                            this.error = 1; // syntax
                            break;
                    }
                }
                break;


            case 1: // getinternode
                switch (token[i]) {
                    case ':':
                    case ',':
                    case ')':
                        state = 2;
                        break;
                    default:
                        state = 99;
                        this.error = 1; // syntax
                        break;
                }
                break;

            case 2: // nextmove
                switch (token[i]) {
                    case ':':
                        i++;
                        if (isNumber(token[i])) {
                            curnode.original_edge_length = token[i];

                            // nnormalization of edge lengths
                            if (edge_length_norm) {
                                curnode.edge_length = Math.sqrt(parseFloat(token[i]) * 1000000) / 1000000;
                            } else {
                                curnode.edge_length = parseFloat(token[i]);
                            }
                            if (parseFloat(token[i]) > 0) {
                                this.has_edge_lengths = true;
                            }
                            i++;
                        }
                        break;

                    case ',':
                        q = this.NewNode();
                        curnode.sibling = q;
                        var c = stack.length;
                        if (c == 0) {
                            state = 99;
                            this.error = 2; // missing (
                        } else {
                            q.ancestor = stack[c - 1];
                            curnode = q;
                            state = 0;
                            i++;
                        }
                        break;

                    case '(':
                        stack.push(curnode);
                        q = this.NewNode();
                        curnode.child = q;
                        q.ancestor = curnode;
                        curnode = q;
                        state = 0;
                        i++;
                        break;

                    case ')':
                        if (stack.length == 0) {
                            state = 99;
                            this.error = 3; // unbalanced
                        } else {
                            curnode = stack.pop();
                            state = 3;
                            i++;
                        }
                        break;

                    case ';':
                        if (stack.length == 0) {
                            state = 99;
                        } else {
                            state = 99;
                            this.error = 4; // stack not empty
                        }
                        break;

                    default:
                        state = 99;
                        this.error = 1; // syntax
                        break;
                }
                break;

                case 3: // finishchildren
                if (ctype_alnum(token[i].charAt(0)) || token[i].charAt(0) == "'" || token[i].charAt(0) == '"' || token[i].charAt(0) == '_') {
                    if (isNaN(token[i]) && mode != 'gene') {
                        curnode.branch_support = token[i];
                        this.has_branch_supports = true;
                    } else {
                        curnode.branch_support = parseFloat(token[i]);
                        this.has_branch_supports = true;
                    }
                    i++;
                } else {
                    switch (token[i]) {
                        case ':':
                            i++;
                            if (isNumber(token[i])) {
                                curnode.original_edge_length = parseFloat(token[i]);

                                // normalization of edge lengths
                                if (edge_length_norm) {
                                    curnode.edge_length = Math.sqrt(parseFloat(token[i]) * 1000000) / 1000000;
                                } else {
                                    curnode.edge_length = parseFloat(token[i]);
                                }
                                if (parseFloat(token[i]) > 0) {
                                    this.has_edge_lengths = true;
                                }
                                i++;
                            }
                            break;

                        case ')':
                            if (stack.length == 0) {
                                state = 99;
                                this.error = 3; // unbalanced
                            } else {
                                curnode = stack.pop();
                                i++;
                            }
                            break;

                        case ',':
                            q = this.NewNode();
                            curnode.sibling = q;

                            if (stack.length == 0) {
                                state = 99;
                                this.error = 2; // missing (
                            } else {
                                q.ancestor = stack[stack.length - 1];
                                curnode = q;
                                state = 0;
                                i++;
                            }
                            break;

                        case ';':
                            state = 2;
                            break;

                        default:
                            state = 99;
                            this.error = 1; // syntax
                            break;
                    }
                }
                break;
        }
    }
};


Tree.prototype.Serialize = function() {
    return this.SerializeNode(this.root) + ";";
};


Tree.prototype.SerializeNode = function(node) {
    var text = "";

    if (node.child) {
        text += "(" + this.SerializeNode(node.child) + ")";
    }

    if (!node.IsLeaf() && this.has_branch_supports) {
        text += node.branch_support;
    } else {
        text += node.label;
    }

    if (node.collapsed) {
        text += '_collapsed';
    }

    if (this.has_edge_lengths) {
        text += ":" + node.original_edge_length;
    }

    if (node.sibling) {
        text += "," + this.SerializeNode(node.sibling);
    }

    return text;
};


Tree.prototype.ComputeDepths = function() {
    for (var i in this.nodes) {
        if (this.nodes[i].IsLeaf()) {
            p = this.nodes[i].ancestor;
            var count = 1;
            while (p) {
                p.depth = Math.max(p.depth, count);
                count++;
                p = p.ancestor;
            }
        }
    }
};


Tree.prototype.FindLowestCommonAncestor = function(p1, p2) {
    let ancestors = new Set([]);

    while (p1) {
        ancestors.add(p1);
        p1 = p1.ancestor;
    }

    while (p2) {
        if (ancestors.has(p2)) {
            return p2;
        }
        p2 = p2.ancestor;
    }

    return null;
}


function NodeIterator(root)
{
    this.root = root;
    this.cur = null;
    this.stack = [];
}


NodeIterator.prototype.Begin = function()
{
/*    if (this.root.constructor === Array)
    {
        this.cur = 0;
        return label_to_node_map[this.root[0]];
    }*/
    this.cur = this.root;
    while (this.cur.child)
    {
        this.stack.push(this.cur);
        this.cur = this.cur.child;
    }
    return this.cur;
};


NodeIterator.prototype.Next = function()
{
/*    if (this.root.constructor === Array)
    {
        this.cur = this.cur + 1;
        if (this.cur >= this.root.length)
        {
            return null;
        }
        return label_to_node_map[this.root[this.cur]];
    }*/
    if (this.stack.length == 0)
    {
        this.cur = null;
    }
    else
    {
        if (this.cur.sibling)
        {
            var p = this.cur.sibling;
            while (p.child)
            {
                this.stack.push(p);
                p = p.child;
            }
            this.cur = p;
        }
        else
        {
            this.cur = this.stack.pop();
        }
    }
    return this.cur;
};

PreorderIterator.prototype = new NodeIterator;


function PreorderIterator()
{
    NodeIterator.apply(this, arguments)
};


PreorderIterator.prototype.Begin = function()
{
    this.cur = this.root;
    return this.cur;
};


PreorderIterator.prototype.Next = function()
{
    if (this.cur.child)
    {
        this.stack.push(this.cur);
        this.cur = this.cur.child;
    }
    else
    {
        while (this.stack.length > 0 && this.cur.sibling == null)
        {
            this.cur = this.stack.pop();
        }
        if (this.stack.length == 0)
        {
            this.cur = null;
        }
        else
        {
            this.cur = this.cur.sibling;
        }
    }
    return this.cur;
};