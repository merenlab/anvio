var ANIMATIONS_ENABLED = true;
var SLIDE_INTERVAL = 4;
var SLIDE_STEP_SIZE = 15;

var is_left_panel_sliding = false;
var is_right_panel_sliding = false;
var is_right_panel2_sliding = false;

function toggleLeftPanel() {
    if (is_left_panel_sliding)
        return;

    is_left_panel_sliding = true;

    if ($('#panel-left').is(':visible')) {
        var animation_frame = function(){ 
            if (ANIMATIONS_ENABLED && $('#panel-left')[0].getBoundingClientRect().right > 0) {
                $('#panel-left').css('left', parseInt($('#panel-left').css('left')) - SLIDE_STEP_SIZE);
                $('#toggle-panel-left').css('left', $('#sidebar')[0].getBoundingClientRect().right + 'px');
                setTimeout(animation_frame, SLIDE_INTERVAL);
            }
            else {
                $('#panel-left').hide();
                $('#toggle-panel-left').css('left', '');
                $('#toggle-panel-left').removeClass('toggle-panel-left-pos');
                $('#toggle-panel-left-inner').html('&#9658;');
                is_left_panel_sliding = false;
            }
        };
        setTimeout(animation_frame, 1);
    } else {
        $('#panel-left').show();
        var animation_frame = function(){ 
            if (ANIMATIONS_ENABLED && $('#panel-left')[0].getBoundingClientRect().left < 0) {
                $('#panel-left').css('left', parseInt($('#panel-left').css('left')) + SLIDE_STEP_SIZE);
                $('#toggle-panel-left').css('left', $('#sidebar')[0].getBoundingClientRect().right + 'px');
                setTimeout(animation_frame, SLIDE_INTERVAL);
            }
            else {
                $('#panel-left').css('left', '0px');
                $('#toggle-panel-left').css('left', '');
                $('#toggle-panel-left').addClass('toggle-panel-left-pos');
                $('#toggle-panel-left-inner').html("&#9664;");
                is_left_panel_sliding = false;
            }
        };
        setTimeout(animation_frame, 1);
    }
}

function toggleRightPanel() {
    if (is_right_panel_sliding)
        return;

    if ($('#description-panel').is(':visible')) {
        if (!is_right_panel2_sliding) {
            toggleRight2Panel();
        }
        setTimeout(toggleRightPanel, SLIDE_INTERVAL * 2);
        return;
    }

    is_right_panel_sliding = true;

    if ($('#mouse_hover_panel').is(':visible')) {
        var animation_frame = function(){ 
            if (ANIMATIONS_ENABLED && $('#mouse_hover_panel')[0].getBoundingClientRect().left < document.body.clientWidth) {
                $('#mouse_hover_panel').css('left', parseInt($('#mouse_hover_panel').css('left')) + SLIDE_STEP_SIZE);
                $('#toggle-panel-right').css('left', $('#mouse_hover_panel')[0].getBoundingClientRect().left - parseInt($('#toggle-panel-right').css('width')) + 'px');
                setTimeout(animation_frame, SLIDE_INTERVAL);
            }
            else {
                $('#mouse_hover_panel').hide();
                $('#toggle-panel-right').css('left', '');
                $('#toggle-panel-right').removeClass('toggle-panel-right-pos');
                $('#toggle-panel-right-inner').html('&#9664;');
                is_right_panel_sliding = false;
            }
        };
        setTimeout(animation_frame, 1);
    } else {
        $('#mouse_hover_panel').show();
        var animation_frame = function(){ 
            if (ANIMATIONS_ENABLED && $('#mouse_hover_panel')[0].getBoundingClientRect().right > document.body.clientWidth) {
                $('#mouse_hover_panel').css('left', parseInt($('#mouse_hover_panel').css('left')) - SLIDE_STEP_SIZE);
                $('#toggle-panel-right').css('left', $('#mouse_hover_panel')[0].getBoundingClientRect().left - parseInt($('#toggle-panel-right').css('width')) + 'px');
                setTimeout(animation_frame, SLIDE_INTERVAL);
            }
            else {
                $('#mouse_hover_panel').css('left', '');
                $('#toggle-panel-right').css('left', '');
                $('#toggle-panel-right').addClass('toggle-panel-right-pos');
                $('#toggle-panel-right-inner').html("&#9658;");
                is_right_panel_sliding = false;
            }
        };
        setTimeout(animation_frame, 1);
    }
}

function toggleRight2Panel() {
    if (is_right_panel2_sliding)
        return;

    if ($('#mouse_hover_panel').is(':visible')) {
        if (!is_right_panel_sliding) {
            toggleRightPanel();
        }
        setTimeout(toggleRight2Panel, SLIDE_INTERVAL * 2);
        return;
    }

    is_right_panel2_sliding = true;

    if ($('#description-panel').is(':visible')) {
        var animation_frame = function(){ 
            if (ANIMATIONS_ENABLED && $('#description-panel')[0].getBoundingClientRect().left < document.body.clientWidth) {
                $('#description-panel').css('left', parseInt($('#description-panel').css('left')) + SLIDE_STEP_SIZE);
                $('#toggle-panel-right-2').css('left', $('#description-panel')[0].getBoundingClientRect().left - parseInt($('#toggle-panel-right-2').css('width')) + 'px');
                setTimeout(animation_frame, SLIDE_INTERVAL);
            }
            else {
                $('#description-panel').hide();
                $('#toggle-panel-right-2').css('left', '');
                $('#toggle-panel-right-2').removeClass('toggle-panel-right-pos-2');
                $('#toggle-panel-right-2-inner').html('&#9664;');
                is_right_panel2_sliding = false;
            }
        };
        setTimeout(animation_frame, 1);
    } else {
        $('#description-panel').show();
        $('#description-panel').css('left', document.body.clientWidth); 
        var animation_frame = function(){ 
            if (ANIMATIONS_ENABLED && $('#description-panel')[0].getBoundingClientRect().right > document.body.clientWidth) {
                $('#description-panel').css('left', parseInt($('#description-panel').css('left')) - SLIDE_STEP_SIZE);
                $('#toggle-panel-right-2').css('left', $('#description-panel')[0].getBoundingClientRect().left - parseInt($('#toggle-panel-right-2').css('width')) + 'px');
                setTimeout(animation_frame, SLIDE_INTERVAL);
            }
            else {
                $('#description-panel').css('left', '');
                $('#toggle-panel-right-2').css('left', '');
                $('#toggle-panel-right-2').addClass('toggle-panel-right-pos-2');
                $('#toggle-panel-right-2-inner').html("&#9658;");
                is_right_panel2_sliding = false;
            }
        };
        setTimeout(animation_frame, 1);
    }
}