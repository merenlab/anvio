/**
 * Javascript library for state migrations
 *
 *  Authors: Ozcan Esen
 *
 * Copyright 2019-2021, The anvi'o project (http://anvio.org)
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

function migrate_state(state) {
    let current_version = state['version'];

    if (current_version == '1') {
        ['samples-layers', 'samples-categorical-colors', 'samples-stack-bar-colors'].forEach((key) => {
            if (state.hasOwnProperty(key)) {
                let backup = $.extend(true, {}, state[key]); // deep copy
                delete state[key];

                state[key] = {'default': backup};
            }

        });

        if (state.hasOwnProperty('samples-layer-order')) {
            let new_order = [];

            state['samples-layer-order'].forEach((sample_layer_name) => {
                new_order.push({
                    'group': 'default',
                    'layer_name': sample_layer_name
                });
            });

            delete state['samples-layer-order'];
            state['samples-layer-order'] = new_order;
        }

        state['version'] = '2';
        current_version = '2';
    }

    if (current_version == '2') {
        let new_stack_bar_colors = {};

        for (let layer_name in state['stack_bar_colors']) {
            let bar_names = layer_name.split('!')[1].split(';');
            new_stack_bar_colors[layer_name] = {};

            for (let j=0; j < bar_names.length; j++) {
                new_stack_bar_colors[layer_name][bar_names[j]] = state['stack_bar_colors'][layer_name][j];
            }
        }

        delete state['stack_bar_colors'];
        state['stack_bar_colors'] = new_stack_bar_colors;

        let new_samples_stack_bar_colors = {};
        for (let group in state['samples-stack-bar-colors']) {
            new_samples_stack_bar_colors[group] = {};
            for (let layer_name in state['samples-stack-bar-colors'][group]) {
                let bar_names = layer_name.split('!')[1].split(';');
                new_samples_stack_bar_colors[group][layer_name] = {};

                for (let j=0; j < bar_names.length; j++) {
                    new_samples_stack_bar_colors[group][layer_name][bar_names[j]] = state['samples-stack-bar-colors'][group][layer_name][j];
                }
            }
        }

        delete state['samples-stack-bar-colors'];
        state['samples-stack-bar-colors'] = new_samples_stack_bar_colors;

        state['version'] = '3';
        current_version = '3';
    }

    if (parseInt(current_version) < parseInt(VERSION)) {
        toastr.error(`Anvi'o failed to upgrade the state. State will be ignored.`);
        return {}
    }

    return state;
}
