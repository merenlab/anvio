/**
 * This is a "plugin" for the sortable package:
 * https://www.npmjs.com/package/sortable-tablesort
 * https://github.com/tofsjonas/sortable
 *
 * Enhances the accessibility of class="sortable" tables by adding ARIA attributes and keyboard event listeners.
 * @param tables - A list of HTML table elements to enhance.
 */
var enhanceSortableAccessibility = function (tables) {
    /**
     * Generates an aria-label attribute for a table header cell based on its content and current sort direction.
     * @param element - The table header cell to update.
     * @param default_direction - The default sort direction for the table.
     */
    function updateAriaLabel(element, default_direction) {
        var _a;
        if (default_direction === void 0) { default_direction = ''; }
        // Generate aria-label based on header content
        var header_text = element.textContent || 'element';
        var current_direction = (_a = element.getAttribute('aria-sort')) !== null && _a !== void 0 ? _a : '';
        var new_direction = 'descending';
        if (current_direction === 'descending' || (default_direction && current_direction !== 'ascending')) {
            new_direction = 'ascending';
        }
        var aria_label = "Click to sort table by ".concat(header_text, " in ").concat(new_direction, " order");
        element.setAttribute('aria-label', aria_label);
        // element.setAttribute('title', aria_label) REMEMBER TO COMMENT OUT WHEN NOT TESTING!!
    }
    /**
     * Handles keyboard events on table header cells and triggers a click event when the Enter key is pressed.
     * @param event - The keyboard event to handle.
     */
    function handleKeyDown(event) {
        if (event.key === 'Enter') {
            var element = event.target;
            element.click();
        }
    }
    // Iterate over each table in the input list
    tables.forEach(function (table) {
        var default_direction = table.classList.contains('asc') ? 'ascending' : '';
        var headers = table.querySelectorAll('th');
        // Iterate over each header cell in the table
        headers.forEach(function (header) {
            var element = header;
            // Skip if the header cell already has a tabindex attribute
            if (element.hasAttribute('tabindex'))
                return;
            var update = function () {
                updateAriaLabel(element, default_direction);
            };
            // Add tabindex attribute and generate initial aria-label attribute
            element.setAttribute('tabindex', '0');
            update();
            // Attach click event listener to update aria-label attribute
            element.addEventListener('click', function () {
                // Add a delay to allow the new sort order to be applied
                setTimeout(update, 50);
            });
            // Attach focus event listener to update aria-label attribute
            element.addEventListener('focus', update);
            // Attach keyboard event listener to trigger click event
            element.addEventListener('keydown', handleKeyDown);
        });
    });
};

// Attach function to DOMContentLoaded event to execute when page is loaded
document.addEventListener('DOMContentLoaded', function () {
    enhanceSortableAccessibility(document.querySelectorAll('.sortable'));
});
