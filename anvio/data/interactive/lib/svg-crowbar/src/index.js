import getSource from './inputProcessor'
import download from './svg'
import downloadPNG from './png'
import {getFilename} from './util'

const downloadSvg = (svgElement, filename, options) =>
  download(getSource(svgElement, options), filename || getFilename(svgElement))
export default downloadSvg
const downloadPng = (svgElement, filename, options) =>
  downloadPNG(
    getSource(svgElement, options),
    filename || getFilename(svgElement),
    options?.downloadPNGOptions,
  )
export {downloadSvg, downloadPng}
