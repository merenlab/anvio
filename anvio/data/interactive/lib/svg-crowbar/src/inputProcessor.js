import {doctype, prefix} from './const'

function getEmptySvgDeclarationComputed() {
  let emptySvg = document.createElementNS(prefix.svg, 'svg')
  document.body.appendChild(emptySvg)
  emptySvg.style.all = 'initial'
  const emptySvgDeclarationComputed = getComputedStyle(emptySvg)
  document.body.removeChild(emptySvg)
  emptySvg = null
  return emptySvgDeclarationComputed
}

function getSource(svg, {css = 'inline'} = {}) {
  if (!(svg instanceof SVGElement)) {
    throw new Error('SVG element is required')
  }

  svg.setAttribute('version', '1.1')
  // removing attributes so they aren't doubled up
  svg.removeAttribute('xmlns')
  svg.removeAttribute('xlink')

  // These are needed for the svg
  if (!svg.hasAttributeNS(prefix.xmlns, 'xmlns')) {
    svg.setAttributeNS(prefix.xmlns, 'xmlns', prefix.svg)
  }

  if (!svg.hasAttributeNS(prefix.xmlns, 'xmlns:xlink')) {
    svg.setAttributeNS(prefix.xmlns, 'xmlns:xlink', prefix.xlink)
  }

  if (css === 'inline') {
    setInlineStyles(svg, getEmptySvgDeclarationComputed())
  } else if (css === 'internal') {
    setInternalStyles(svg)
  }

  const source = new XMLSerializer().serializeToString(svg)
  const rect = svg.getBoundingClientRect()

  const result = {
    top: rect.top,
    left: rect.left,
    width: rect.width,
    height: rect.height,
    class: svg.getAttribute('class'),
    id: svg.getAttribute('id'),
    name: svg.getAttribute('name'),
    childElementCount: svg.childElementCount,
    source: doctype + source,
  }

  return result
}

function setInlineStyles(svg, emptySvgDeclarationComputed) {
  function explicitlySetStyle(element) {
    const cSSStyleDeclarationComputed = getComputedStyle(element)
    let key
    let value
    let computedStyleStr = ''
    for (let i = 0, len = cSSStyleDeclarationComputed.length; i < len; i++) {
      key = cSSStyleDeclarationComputed[i]
      value = cSSStyleDeclarationComputed.getPropertyValue(key)
      if (value !== emptySvgDeclarationComputed.getPropertyValue(key)) {
        computedStyleStr += `${key}:${value};`
      }
    }
    element.setAttribute('style', computedStyleStr)
  }
  function traverse(obj) {
    const tree = []
    tree.push(obj)
    visit(obj)
    function visit(node) {
      if (node && node.hasChildNodes()) {
        let child = node.firstChild
        while (child) {
          if (child.nodeType === 1 && child.nodeName !== 'SCRIPT') {
            tree.push(child)
            visit(child)
          }
          child = child.nextSibling
        }
      }
    }
    return tree
  }
  // hardcode computed css styles inside svg
  const allElements = traverse(svg)
  let i = allElements.length
  while (i--) {
    explicitlySetStyle(allElements[i])
  }
}

function setInternalStyles(svg) {
  const style = document.createElement('style')
  style.innerHTML = Array.from(document.styleSheets)
    .filter(
      (styleSheet) =>
        // Prevent CORS errors
        !styleSheet.href || styleSheet.href.startsWith(document.location.origin),
    )
    .map((styleSheet) =>
      Array.from(styleSheet.cssRules)
        .map((rule) => rule.cssText)
        .join(' '),
    )
    .join(' ')
  svg.prepend(style)
}

export default getSource
