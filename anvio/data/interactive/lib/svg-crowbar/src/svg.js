import {commenceDownload} from './util'
import {DEFAULT_FILENAME} from './const'

function download(source, filename = DEFAULT_FILENAME) {
  const url = URL.createObjectURL(new Blob([source.source], {type: 'text/xml'}))

  commenceDownload(`${filename}.svg`, url, () => URL.revokeObjectURL(url))
}

export default download
