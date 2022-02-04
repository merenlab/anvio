# SVG Crowbar Library
[![NPM version](https://img.shields.io/npm/v/svg-crowbar.svg)](https://www.npmjs.com/package/svg-crowbar)
[![code style: prettier](https://img.shields.io/badge/code_style-prettier-ff69b4.svg?style=flat-square)](https://github.com/prettier/prettier)
[![Build status](https://github.com/cy6erskunk/svg-crowbar/actions/workflows/node.js.yml/badge.svg)](https://github.com/cy6erskunk/svg-crowbar/actions/workflows/node.js.yml)
![Copy README to gh-pages branch](https://github.com/cy6erskunk/svg-crowbar/workflows/Copy%20README%20to%20gh-pages%20branch/badge.svg)

A standalone 3.5Kb JS client library based on Chrome [bookmarklet](https://nytimes.github.io/svg-crowbar/).

The library provides functionality to trigger a download of a given SVG file having all the styles inlined,
to make it look the same when opened independently from the original HTML page.

It is also possible to use this library to convert an SVG to a PNG before downloading.

## Usage
```javascript
import downloadSvg from 'svg-crowbar';

downloadSvg(document.querySelector('svg'));
```    
or
```javascript
import { downloadPng } from 'svg-crowbar';

downloadPng(document.querySelector('svg'), 'my_svg', { css: 'internal' });
```

The `downloadSVG`/`downloadPNG` functions each have three arguments:

```javascript
downloadSVG(svgElement, [filename], [options])
downloadPNG(svgElement, [filename], [options])
```

### Options

- **svgElement** *(required)*
  
  A DOM element selector for an SVG, e.g. `document.querySelector('svg')`. An error is thrown if no valid SVG element was provided.

- **filename** *(optional)*

  A string to set the filename. This is determined by element id, class or page title, when not provided explicitly.

- **options** *(optional)*

  An object literal. It presently has two configurable properties:

- **options.css** *(optional)*

  This setting determines how the SVG will be styled:

  - **`'inline'`**

    Default value. Inlines all computed styles on every element in the SVG. This setting best ensures that the exported SVG is accurate cross-browser.

  - **`'internal'`**

    Adds an internal block of styles containing only explicitly declared style rules (from `document.styleSheets`). This can drastically reduce file-sizes and build time in exported SVGs, but could be less accurate as it does not include styles from the browser's user agent stylesheet, or from cross-origin stylesheets (e.g. external webfonts).

  - **`'none'`**

    Doesn't add any CSS. This gives the smallest file-size, but you might need to manually add your own styles to exported SVGs to ensure an accurate output. You can do this by injecting a `<style>` block into the selected SVG before exporting.

  Example:
  ```javascript
  const svg = document.querySelector('svg');

  // Add inline styles on SVG elements:
  downloadSvg(svg, 'my_svg'); 
  downloadSvg(svg, 'my_svg', { css: 'inline' });

  // Add a <style> block in the SVG:
  downloadSvg(svg, 'my_svg', { css: 'internal' });

  // Do not add CSS:
  downloadSvg(svg, 'my_svg', { css: 'none' });
  ```

- **options.downloadPNGOptions.scale** *(optional)*

  This setting determines at which scale the final PNG image is created, for example when resolution is desired. The default scale is 1 (ie no scaling).

  Example:
  ```javascript
  const svg = document.querySelector('svg');

  // Download a normal-scaled PNG
  downloadPng(svg, 'my_svg');
  downloadPng(svg, 'my_svg', {downloadPNGOptions:{ scale: 1 }});

  // Download a double-scaled PNG
  downloadPng(svg, 'my_svg', {downloadPNGOptions:{ scale: 2 }});

  ```

  ### UMD bundle

Thanks to @richardwestenra there's UMD bundle available in the package: 
simply add 
```html
<script src="node_modules/svg-crowbar/dist/main.js"></script>
```
to get `downloadSvg` and `downloadPng` global function.