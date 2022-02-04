import {prefix} from '../src/const'

global.createSVG = () => document.createElementNS(prefix.svg, 'svg')
