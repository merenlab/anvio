const doctype =
  '<?xml version="1.0" standalone="no"?><!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">'
const prefix = {
  xmlns: 'http://www.w3.org/2000/xmlns/',
  xlink: 'http://www.w3.org/1999/xlink',
  svg: 'http://www.w3.org/2000/svg',
}
const REMOVE_TIMEOUT = 10
const DEFAULT_FILENAME = 'untitled'

const downloadSvg = (svgElement, filename, options) => download(getSource(svgElement, options), filename || getFilename(svgElement))

function download(source, filename = DEFAULT_FILENAME) {
  const url = URL.createObjectURL(new Blob([source.source], {type: 'text/xml'}))

  commenceDownload(`${filename}.svg`, url, () => URL.revokeObjectURL(url))
}

function getFilename(source) {
  if (!(source instanceof SVGElement)) {
    throw new Error('SVG Element is required')
  }

  return (
    source.getAttribute('id') ||
    source.getAttribute('class') ||
    document.title.replace(/[^a-z0-9]/gi, '-').toLowerCase() ||
    DEFAULT_FILENAME
  )
}

function commenceDownload(filename, imgdata, callback) {
  const a = document.createElement('a')
  document.body.appendChild(a)
  a.setAttribute('class', 'svg-crowbar')
  a.setAttribute('download', filename)
  a.setAttribute('href', imgdata)
  a.style.display = 'none'
  a.click()

  setTimeout(() => {
    if (callback) {
      callback()
    }
    document.body.removeChild(a)
  }, REMOVE_TIMEOUT)
}

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

  // if (css === 'inline') {
  //   setInlineStyles(svg, getEmptySvgDeclarationComputed())
  // } else if (css === 'internal') {
  //   setInternalStyles(svg)
  // }

  const source = new XMLSerializer().serializeToString(svg)
  console.log(source);
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