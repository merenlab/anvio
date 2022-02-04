import download, {_fixSource} from '../src/png'
import inputProcessor from '../src/inputProcessor'

jest.mock('../src/util')

test('png download requires source', () => {
  expect(download).toThrow()
})

test('download succeeds with empty SVG', () => {
  expect(() => download(inputProcessor(createSVG()))).not.toThrow()
})

test('download uses provided safeSource fn', () => {
  const safeFnMock = jest.fn()
  download(inputProcessor(createSVG()), undefined, {fixSource: safeFnMock})
  expect(safeFnMock).toHaveBeenCalled()
})

test('download succeeds with non-ascii chars ☸☹☺☻☼☾☿✓', () => {
  const source = `<svg><text>boom! ☸☹☺☻☼☾☿✓</text></svg>`
  expect(() => _fixSource(source)).not.toThrow()
})

test.todo('commenceDownload have been called')
