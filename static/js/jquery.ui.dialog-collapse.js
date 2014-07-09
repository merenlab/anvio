/*
 * Extends JQuery UI Dialog to add collapse button feature.
 *
 * Copyright 2013.  Marko Martinović
 * http://www.techytalk.info
 */
(function($) {
    // Add default options and event callbacks
    $.extend($.ui.dialog.prototype.options, {
        collapseEnabled: true,
        beforeCollapse: null,
        collapse: null,
        beforeCollapseRestore: null,
        collapseRestore: null
    });

    // Backup old _init
    var _init = $.ui.dialog.prototype._init;

    // New _init
    $.ui.dialog.prototype._init = function() {
        // Apply old _init
        _init.apply(this, arguments);

        // Holds original this.options.resizable
        var resizableOld = null;
        if(this.options.collapseEnabled) {
            this.addCollapseButton = function() {
                // Hide the restore button if it exists
                if(this.uiDialogTitlebarCollapseRestore)
                    this.uiDialogTitlebarCollapseRestore.hide();

                // Add the collapse button if it doesn't exists
                if(!this.uiDialogTitlebarCollapse) {
                    this.uiDialogTitlebarCollapse = $('<button></button>')
                        .button({
                            label: 'collapse',
                            icons: {
                                primary: 'ui-icon ui-icon-arrowthickstop-1-n'
                            },
                            text: false
                        })
                        .addClass('ui-dialog-titlebar-collapse')
                        .appendTo( this.uiDialogTitlebar );
                    this._on( this.uiDialogTitlebarCollapse, {
                        // Run this.collapse on click
                        click: function( event ) {
                            event.preventDefault();
                            this.collapse( event );
                        }
                    });
                } else {
                    this.uiDialogTitlebarCollapse.show();
                }
            }

            this.addCollapseRestoreButton = function() {
                // Hide the collapse button if it exists
                if(this.uiDialogTitlebarCollapse)
                    this.uiDialogTitlebarCollapse.hide();

                // Add the restore button if it doesn't exists
                if(!this.uiDialogTitlebarCollapseRestore){
                    this.uiDialogTitlebarCollapseRestore = $('<button></button>')
                        .button({
                            label: 'restore',
                            icons: {
                                primary: 'ui-icon ui-icon-arrowthickstop-1-s'
                            },
                            text: false
                        })
                        .addClass('ui-dialog-titlebar-collapse-restore')
                        .appendTo( this.uiDialogTitlebar );
                    this._on( this.uiDialogTitlebarCollapseRestore, {
                        // Run this.restore on click
                        click: function( event ) {
                            event.preventDefault();
                            this.restore( event );
                        }
                    });
                } else {
                    this.uiDialogTitlebarCollapseRestore.show();
                }
            }

            this.collapse = function(event) {
                var self = this;

                // Allow people to abort collapse event
                if (false === self._trigger('beforeCollapse')) {
                    return;
                }

                // slideUp the dialog element
                this.element.slideUp('fast', function() {
                    // Deal with the resizable option
                    if(self.options.resizable){
                        // Backup old resizable option
                        resizableOld = self.options.resizable;

                        // Destroy the resizable and set dialog height to auto
                        self.uiDialog.resizable('destroy').css('height', 'auto');

                        // Overwrite original resizable option to disable vertical resize
                        self.options.resizable = 'e, w';

                        // Make resizable with the new resizable option
                        self._makeResizable();
                    }

                    // Replace collapse button with restore button
                    self.addCollapseRestoreButton();

                    // Trigger collapse event
                    self._trigger('collapse');
                });

                return self;
            };

            this.restore = function(event) {
                var self = this;

                // Allow people to abort restore event
                if (false === self._trigger('beforeCollapseRestore')) {
                        return;
                }

                // slideDown the dialog element
                this.element.slideDown('fast', function() {
                    // Deal with the resizable option
                    if(self.options.resizable){
                        // Destroy our horizontal only resize
                        self.uiDialog.resizable('destroy');

                        // Restore original resizable option from backup
                        self.options.resizable = resizableOld;

                        // Make resizable with the original resizable option
                        self._makeResizable();
                    }

                    // Replace restore button with collapse button
                    self.addCollapseButton();

                    // Trigger collapse event
                    self._trigger('collapseRestore');
                });

                return self;
            };

            // By default add both buttons, collapse button will hide restore
            this.addCollapseRestoreButton();
            this.addCollapseButton();

            // Deal with collapse and restore buttons position if close button is visible
            if(this.uiDialogTitlebarClose && this.uiDialogTitlebarClose.is(':visible')) {
                var right = parseFloat(this.uiDialogTitlebarClose.css('right'));

                $('.ui-dialog-titlebar-collapse, .ui-dialog-titlebar-collapse-restore')
                    .css('right', 2*right+this.uiDialogTitlebarClose.outerWidth()+'px');
            }
        }
    };
}(jQuery));