import download from '../src/svg'
import inputProcessor from '../src/inputProcessor'
import {commenceDownload} from '../src/util'
import {DEFAULT_FILENAME} from '../src/const'

jest.mock('../src/util')

test('download requires source', () => {
  expect(download).toThrow()
})

test('download succeeds with empty SVG', () => {
  expect(() => download(inputProcessor(createSVG()))).not.toThrow()
  expect(commenceDownload).toHaveBeenCalledTimes(1)
  expect(commenceDownload.mock.calls[0][0]).toEqual(`${DEFAULT_FILENAME}.svg`)
})
