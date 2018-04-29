
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
            let new_order = convert_samples_order_to_array(state['samples-layer-order']);
            delete state['samples-layer-order'];
            state['samples-layer-order'] = new_order;
        }

        state['version'] = '2';
        current_version = '2';
    }

    if (parseInt(current_version) < parseInt(VERSION)) {
        toastr.info(`Anvi'o failed to upgrade the state.`);
        throw "";
    }

    return state;
}
