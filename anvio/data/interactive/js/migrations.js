
function migrate_state(state) {
    let current_version = state['version'];

    if (current_version == '1') {
        if (state.hasownProperty('samples-layers')) {
            let layer_info = $.extend(true, {}, state['samples-layers']);
            delete state['samples-layers']
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
