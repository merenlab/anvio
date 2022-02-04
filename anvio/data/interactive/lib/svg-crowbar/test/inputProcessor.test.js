import inputProcessor from '../src/inputProcessor'
import {prefix, doctype} from '../src/const'

const createSVG = () => document.createElementNS(prefix.svg, 'svg')

describe('inputProcessor', () => {
  test('does not add styles to empty svg', () => {
    expect(inputProcessor).toThrow()
  })

  test('empty SVG', () => {
    const result = inputProcessor(createSVG())
    expect(result.id).toBe(null)
    expect(result.class).toBe(null)
    expect(result.childElementCount).toBe(0)
    expect(result.top).toBe(0)
    expect(result.left).toBe(0)
    expect(result.width).toBe(0)
    expect(result.height).toBe(0)
    const svgString = result.source.replace(doctype, '')
    const foo = document.createElement('div')
    foo.innerHTML = svgString
    expect(foo.firstElementChild.tagName).toBe('svg')
    expect(foo.firstElementChild.getAttribute('style')).toBe('')
    expect(foo.firstElementChild.children.length).toBe(0)
  })
})
