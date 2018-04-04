/**
 *  Handles right click menu
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


ContextMenu = function(options) {
    this.container = options.container;
    this.event = event;
    this.node = options.node;
    this.layer = options.layer;
}

ContextMenu.prototype.BuildMenu = function() {
    var menu = [];

    if (this.node.IsLeaf()) {
        menu.push({'type': 'item', 'title': 'Add item to bin'     , 'action': 'select'},
                  {'type': 'item', 'title': 'Remove item from bin', 'action': 'remove'},
                  {'type': 'divider'},
                  {'type': 'item', 'title': 'Select layer'  , 'action': 'select_layer'},
                  {'type': 'item', 'title': 'Unselect layer', 'action': 'unselect_layer'});

        if (mode == 'gene') {
            menu.push({'type': 'divider'},
                      {'type': 'item', 'title': 'Inspect Context', 'action': 'inspect', 'param': 'context'},
                      {'type': 'item', 'title': 'Inspect Gene'   , 'action': 'inspect', 'param': 'gene'});            
        }        
        else if (mode == 'pan') {
            menu.push({'type': 'divider'},
                      {'type': 'item', 'title': 'Inspect', 'action': 'inspect', 'param': 'geneclusters'});            
        } else {
            menu.push({'type': 'divider'},
                      {'type': 'item', 'title': 'Inspect', 'action': 'inspect', 'param': 'split'}); 
        }
    } else {
        if (this.node.collapsed) {
            menu.push({'type': 'item', 'title': 'Expand this branch', 'action': 'expand'});
        } else {
            menu.push({'type': 'item', 'title': 'Add branch to bin'     , 'action': 'select'},
                      {'type': 'item', 'title': 'Remove branch from bin', 'action': 'remove'},
                      {'type': 'divider'},
                      {'type': 'item', 'title': 'Collapse this branch', 'action': 'collapse'},
                      {'type': 'item', 'title': 'Rotate this branch'  , 'action': 'rotate'},
                      {'type': 'item', 'title': 'Reroot tree here'    , 'action': 'reroot'});   
        }
    }

    return menu;
};


ContextMenu.prototype.Show = function() {
    this.container.querySelectorAll('.context-menu').forEach((menu) => {
        menu.remove();
    });

    var list = document.createElement('ul');
    list.setAttribute('class', 'dropdown-menu context-menu');
    list.setAttribute('role', 'menu');
    list.style.display = 'block';
    list.style.left = this.event.clientX + 'px';
    list.style.top = this.event.clientY + 'px';

    for (const item of this.BuildMenu()) {
        if (item.type == 'divider') {
            list.innerHTML += `<li class="divider"></li>`;
        } else {
            list.innerHTML += `<li><a href="#" menu-action="${item.action}" menu-param="${item.param}">${item.title}</a></li>`;
        }
    }

    this.container.appendChild(list);

    list.addEventListener('click', (event) => {
        let action = event.target.getAttribute('menu-action');
        let param = event.target.getAttribute('menu-param');
        this.MenuClickHandler(action, param);
    });

    this.container.addEventListener('click', (event) => {
        list.remove();
    }, {once: true});
};


ContextMenu.prototype.MenuClickHandler = function(action, param) {
    if (action == 'select') {
        bins.AppendNode(this.node);
    } 
    else if (action == 'remove') {
        bins.RemoveNode(this.node);
    }
    else if (action == 'select_layer') {
        $('#tbody_layers tr:nth-child(' + this.layer + ') input:checkbox').prop('checked', true);
    }
    else if (action == 'unselect_layer') {
        $('#tbody_layers tr:nth-child(' + this.layer + ') input:checkbox').prop('checked', false);
    }
    else if (action == 'inspect') {
        localStorage.state = JSON.stringify(serializeSettings(true), null, 4);
        window.open(generate_inspect_link('inspect_' + param, item_name), '_blank');
    }
};
