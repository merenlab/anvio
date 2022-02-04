import {prefix, DEFAULT_FILENAME} from '../src/const'
import {getFilename} from '../src/util'

const createSVG = () => document.createElementNS(prefix.svg, 'svg')
describe('getFilename', () => {
  beforeEach(() => (document.title = ''))
  test('throws when receves not  an SVG', () => {
    expect(getFilename).toThrow()
  })

  test('uses default ', () => {
    expect(getFilename(createSVG())).toBe(DEFAULT_FILENAME)
  })
  test('uses id', () => {
    const id = 'boom'
    const elem = createSVG()
    elem.id = id
    expect(getFilename(elem)).toEqual(id)
  })

  test('uses classname', () => {
    const classname = 'whatever'
    const elem = createSVG()
    elem.classList.add(classname)
    expect(getFilename(elem)).toEqual(classname)
  })

  test('uses doc title', () => {
    const title = 'boom'
    document.title = title
    expect(getFilename(createSVG())).toBe(title)
  })
})
