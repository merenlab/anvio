import {REMOVE_TIMEOUT, DEFAULT_FILENAME} from './const'

export function getFilename(source) {
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

export function commenceDownload(filename, imgdata, callback) {
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
