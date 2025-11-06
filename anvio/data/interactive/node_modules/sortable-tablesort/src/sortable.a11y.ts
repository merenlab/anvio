import { enhanceSortableAccessibility } from './enhanceSortableAccessibility'

// Attach function to DOMContentLoaded event to execute when page is loaded
document.addEventListener('DOMContentLoaded', () => {
  enhanceSortableAccessibility(document.querySelectorAll<HTMLTableElement>('.sortable'))
})
