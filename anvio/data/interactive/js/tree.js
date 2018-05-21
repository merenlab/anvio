/**
 * Javascript library to parse newick trees
 *
 *  Author: Özcan Esen <ozcanesen@gmail.com>
 *  Credits: A. Murat Eren
 *  Copyright 2015, The anvio Project
 *
 * This file is part of anvi'o (<https://github.com/meren/anvio>).
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
    this.label = null;
    this.id = 0;
    this.xy = [];
    this.original_edge_length = 0.0;
    this.edge_length = 0.0;
    this.path_length = 0.0;
    this.depth = 0;
    this.order = null;
    this.max_child_path = 0;
}


Node.prototype.SetLabel = function(label) {
    this.collapsed = (label && label.endsWith('_collapsed')) ? true : false;
    this.label = (label) ? label.replace(/_collapsed$/,'') : '';
};


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
    this.rooted = true;
    this.has_edge_lengths = false;
    this.error = 0;
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
                    this.num_leaves++;
                    label = token[i];
                    curnode.SetLabel(label);
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
                    curnode.SetLabel(token[i]);
                    i++;
                } else {
                    switch (token[i]) {
                        case ':':
                            i++;
                            if (isNumber(token[i])) {
                                curnode.original_edge_length = token[i];
                                
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

Tree.prototype.FindNode = function(label) {
    var n = new NodeIterator(this.root);
    var q = n.Begin();
    while (q != null)
    {
        if (q.label == label)
            return q;

        q=n.Next();
    }

    console.log("Couldn't find item with label '" + label + "'");
    return null;
};


Tree.prototype.Serialize = function() {
    return this.SerializeNode(this.root) + ";";
};


Tree.prototype.SerializeNode = function(node) {
    var text = "";

    if (node.child) {
        text += "(" + this.SerializeNode(node.child) + ")";
    }

    text += node.label; 

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


function NodeIterator(root)
{
    this.root = root;
    this.cur = null;
    this.stack = [];
}


NodeIterator.prototype.Begin = function() 
{
    if (this.root.constructor === Array)
    {
        this.cur = 0;
        return label_to_node_map[this.root[0]];
    }
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
    if (this.root.constructor === Array)
    {
        this.cur = this.cur + 1;
        if (this.cur >= this.root.length)
        {
            return null;
        }
        return label_to_node_map[this.root[this.cur]];
    }
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


