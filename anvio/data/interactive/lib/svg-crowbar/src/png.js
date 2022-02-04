import {DEFAULT_FILENAME} from './const'
import {commenceDownload} from './util'

export function _fixSource(source) {
  return btoa(
    unescape(
      encodeURIComponent(source.replace(/[\u00A0-\u2666]/g, (c) => `&#${c.charCodeAt(0)};`)),
    ),
  )
}

const DEFAULT_OPTIONS = {
  debug: false,
  fixSource: _fixSource,
  scale: 1,
}

function downloadPng(
  source,
  filename = DEFAULT_FILENAME,
  {
    debug = DEFAULT_OPTIONS.debug,
    fixSource = DEFAULT_OPTIONS.fixSource,
    scale = DEFAULT_OPTIONS.scale,
  } = DEFAULT_OPTIONS,
) {
  const canvas = document.createElement('canvas')
  const dpr = window.devicePixelRatio || 1
  document.body.appendChild(canvas)
  canvas.setAttribute('id', 'svg-image')
  canvas.setAttribute('width', source.width * dpr * scale)
  canvas.setAttribute('height', source.height * dpr * scale)
  if (debug === false) {
    canvas.style.display = 'none'
  }

  const context = canvas.getContext('2d')
  const imgsrc = `data:image/svg+xml;base64,${fixSource(source.source)}`
  const image = new Image()

  function onLoad() {
    context.scale(dpr * scale, dpr * scale)
    context.drawImage(image, 0, 0)
    const canvasdata = canvas.toDataURL('image/png')

    if (debug === false) {
      commenceDownload(`${filename}.png`, canvasdata, () => document.body.removeChild(canvas))
    }
  }

  image.onload = onLoad
  image.src = imgsrc
  if (debug === true) {
    document.body.appendChild(image)
  }
}

export default downloadPng
