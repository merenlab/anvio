<!-- markdownlint-disable MD033 MD026 -->
<h1>sortable</h1>
<h5>- a tiny, vanilla/plain JavaScript table sorter</h5>

Makes any table with **class="sortable"**, er, sortable. The user can click on a table header and change the sorting of the table rows.

Just include the JavaScript and it will work. No function calls are needed, everything is handled by an **eventListener**.

<h2>Demo</h2>

You can find a simple demo on <https://tofsjonas.github.io/sortable/>

<h2>Table of Contents</h2>

- [Factoids](#factoids)
- ["Installation"](#installation)
  - [1. link to jsDelivr](#1-link-to-jsdelivr)
  - [2. copy file to assets folder](#2-copy-file-to-assets-folder)
  - [3. npm package](#3-npm-package)
    - [a) use links in the html](#a-use-links-in-the-html)
    - [b) import files in javascript](#b-import-files-in-javascript)
- [Non-sortable field](#non-sortable-field)
  - [...using `class="no-sort"`](#using-classno-sort)
  - [...using `CSS`](#using-css)
  - [...using `td` instead of `th`](#using-td-instead-of-th)
- [Indicators/arrows on the left side](#indicatorsarrows-on-the-left-side)
- [NOTE ABOUT CSS/SCSS](#note-about-cssscss)
- [Sorting sizes, dates and such](#sorting-sizes-dates-and-such)
- [Alternative sorting](#alternative-sorting)
- [Colspans/Sort on specific column](#colspanssort-on-specific-column)
- [Ascending sort](#ascending-sort)
- [Tiebreaker / secondary sort](#tiebreaker--secondary-sort)
- [Empty/null rows always last](#emptynull-rows-always-last)
- [Accessibility](#accessibility)
- [Sort on load](#sort-on-load)
- [Thank you...](#thank-you)

## Factoids

- **1180 bytes** minified. (617 bytes gzipped)

- Works with **JavaScript generated tables**. (since we are using an eventListener)

- **Lightning fast**. _Huge_ tables will make it slow and may freeze the browser, especially for mobiles, so you know...

- Requires **thead** and **tbody**.

- **cross browser**, ~~ie9+~~ No longer ie9 compatible. Then again, maybe it already wasn't 🤷

- NOT tested with React, Angular, Vue, etc.

- Works with [Svelte](https://svelte.dev/)!

## "Installation"

There are three ways to use sortable, all of which have their pros and cons. [S Anand](https://github.com/sanand0) and [dkhgh](https://github.com/dkhgh) had some [interesting thoughts](https://github.com/tofsjonas/sortable/issues/28) about it.

1. Include a link to [jsDelivr](https://www.jsdelivr.com/package/gh/tofsjonas/sortable). (easiest)

2. Copy the file from [jsDelivr](https://www.jsdelivr.com/package/gh/tofsjonas/sortable) or [Github](https://github.com/tofsjonas/sortable) and put it in your assets folder. (in between)

3. Install the [npm package](https://www.npmjs.com/package/sortable-tablesort). (most reliable)

### 1. link to jsDelivr

```html
<table class="sortable">
  <thead>
    <tr>
      <th><span>Role</span></th>
      <th>Name</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>Genius</td>
      <td>Rick</td>
    </tr>
    <tr>
      <td><a href="javascript:alert('Inline javascript works!');">Sidekick</a></td>
      <td>Morty</td>
    </tr>
  </tbody>
</table>
<link href="https://cdn.jsdelivr.net/gh/tofsjonas/sortable@latest/sortable.min.css" rel="stylesheet" />
<script src="https://cdn.jsdelivr.net/gh/tofsjonas/sortable@latest/sortable.min.js"></script>
```

The `span` on line four is just there to prove that you can have elements inside `th`!

⚠️ _If you are concerned about bugs, I recommend using version numbers instead of "latest"._

### 2. copy file to assets folder

Same as above, but link to your own files

```html
...
<link href="/assets/sortable.min.css" rel="stylesheet" />
<script src="/assets/sortable.min.js"></script>
...
```

### 3. npm package

First,

```bash
npm install sortable-tablesort
# yarn add sortable-tablesort
# pnpm install sortable-tablesort
```

Now you can

#### a) use links in the html

Same as above, with links to files

```html
...
<link href="./node_modules/sortable-tablesort/sortable.min.css" rel="stylesheet" />
<script src="./node_modules/sortable-tablesort/sortable.min.js"></script>
...
```

or

#### b) import files in javascript

```javascript
// main.js
import 'sortable-tablesort/sortable.min.css'
import 'sortable-tablesort/sortable.min.js'
```

## Non-sortable field

### ...using `class="no-sort"`

If you wish to disable sorting for a specific field, the easiest (and best) way is to add `class="no-sort"` to it, like so:

```html
<thead>
  <tr>
    <th class="no-sort">Role</th>
    <th>Name</th>
  </tr>
</thead>
```

Sorting will not be triggered if you click on "Role".

### ...using `CSS`

This is a bit trickier, but it doesn't require any changes to the html, so I guess it could be worth it in some cases.

```css
/* the first column in every sortable table should not be sortable*/
.sortable th:nth-child(1) {
  pointer-events: none;
}

/* the seventh column in the second .sortable table should not be sortable*/
.sortable:nth-of-type(2) th:nth-child(7) {
  pointer-events: none;
}
```

### ...using `td` instead of `th`

The eventListener only triggers on `th`, not `td`, so this would disable sorting for "Role":

```html
<thead>
  <tr>
    <td>Role</td>
    <th>Name</th>
  </tr>
</thead>
```

⚠️ _Since `th` and `td` are not the same thing, you would most likely still have to use CSS to make them look the way you want. (It might also mess with accessibility.) In **some** cases it could be worth it, but I recommend the `.no-sort` alternative_.

## Indicators/arrows on the left side

If you have text that is aligned on the right side, you may want to have the arrows on the left side.

This is solved by adding a class to the css and using `::before` instead of `::after`.

(You can of course use a pure css solution, without class names - just like with the [non-sortable field](#non-sortable-field) - but _that_ I will leave for you to figure out.)

```css
.sortable th.indicator-left::after {
  content: '';
}
.sortable th.indicator-left::before {
  margin-right: 3px;
  content: '▸';
}
/* etc. */
```

> _Full example: [CSS](https://github.com/tofsjonas/sortable/blob/main/sortable-base.css), [SCSS](https://github.com/tofsjonas/sortable/blob/main/src/sortable-base.scss)_

## NOTE ABOUT CSS/SCSS

The `css/scss` in this repo was only ever meant as an example. It was never intended to be actually _used_.

That said, if you're feeling lazy, here are two stylesheets you can use:

```html
<!-- This will add arrows, as well as support for .no-sort and .indicator-left -->
<link href="https://cdn.jsdelivr.net/gh/tofsjonas/sortable@latest/sortable-base.min.css" rel="stylesheet" />

<!-- This will make it look like the tables in the example, with arrows, striped rows etc. -->
<link href="https://cdn.jsdelivr.net/gh/tofsjonas/sortable@latest/sortable.min.css" rel="stylesheet" />
```

## Sorting sizes, dates and such

Using the `data-sort` attribute in `tbody` > `td` you can have one visible value and one sortable value. This is useful in case you have for instance sizes like kb, Mb, GB, or really weird date formats. 😉

```html
<table class="sortable">
  <thead>
    <tr>
      <th>Movie Name</th>
      <th>Size</th>
      <th>Release date</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>Zack Snyder's Justice League</td>
      <td data-sort="943718400">900MB</td>
      <td data-sort="20210318">03/18/2021</td>
    </tr>
    <tr>
      <td>The Sound of Music</td>
      <td data-sort="1610612736">1.5GB</td>
      <td data-sort="19651209">12/09/1965</td>
    </tr>
  </tbody>
</table>
```

## Alternative sorting

If you click on a table header while holding **shift** or **alt** an alternative
`data-sort-alt` attribute will override `data-sort`.

```html
<table class="sortable">
  <thead>
    <tr>
      <th>Movie Name</th>
      <th>Size</th>
      <th>Release date</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>Something</td>
      <td data-sort-alt="c" data-sort="a">A</td>
      <td data-sort-alt="b" data-sort="c">B</td>
      <td data-sort-alt="a" data-sort="b">C</td>
    </tr>
    <tr>
      <td>Something else</td>
      <td data-sort-alt="e" data-sort="f">D</td>
      <td data-sort-alt="f" data-sort="e">E</td>
      <td data-sort-alt="d" data-sort="d">F</td>
    </tr>
  </tbody>
</table>
```

## Colspans/Sort on specific column

Using the `data-sort-col` attribute in `thead` > `th`, you can sort on a different column than the one that was clicked. For instance if you want to have colspans. Like so:

```html
<thead>
  <tr>
    <th></th>
    <th>Category</th>
    <th class="show_name">Show</th>
    <th colspan="2">Overall</th>
    <th colspan="2" data-sort-col="5">On Our Dates</th>
    <th data-sort-col="7">First Sold Out</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td class="tags">&nbsp;</td>
    <td class="category">Comedy</td>
    <td class="show_name">Show 1</td>
    <td class="ratio all" data-sort="72">18/25</td>
    <td class="pct all">72%</td>
    <td class="ratio ours" data-sort="75">3/4</td>
    <td class="pct ours">75%</td>
    <td>2022-07-30</td>
  </tr>
  ...
</tbody>
```

## Ascending sort

By adding `asc` to `table`, the default sorting direction will be **ascending** instead of descending

```html
<table class="sortable asc">
  <thead>
    ...
  </thead>
  <tbody>
    ...
  </tbody>
</table>
```

## Tiebreaker / secondary sort

If you wish to sort by a different column when two values are equal, you can use the `data-sort-tbr` attribute, like so:

```html
<table class="sortable asc">
  <thead>
    <tr>
      <th data-sort-tbr="1">Year</th>
      <th>Month</th>
      <th>Day</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>2010</td>
      <td>07</td>
      <td>25</td>
    </tr>
    <tr>
      <td>2010</td>
      <td>11</td>
      <td>12</td>
    </tr>
    <tr>
      <td>2010</td>
      <td>04</td>
      <td>25</td>
    </tr>
  </tbody>
</table>
```

When clicking **Year**, if they are the same, we will sort on **Month**.

## Empty/null rows always last

Adding `class="n-last"` to `<table class="sortable">` will make empty/null values always be sorted last, similar to what SQL does with `ORDER BY foo NULLS LAST`.

```html
<table class="sortable n-last">
  <thead>
    <tr>
      <th>Text</th>
      <th class="indicator-left">Number</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>jkl</td>
      <td>0.4</td>
    </tr>
    <tr>
      <td>This will always be sorted after the others</td>
      <td></td>
    </tr>
    <tr>
      <td>abc</td>
      <td>0</td>
    </tr>
    <tr>
      <td>def</td>
      <td>0.2</td>
    </tr>
  </tbody>
</table>
```

⚠️ _Note that a string of blank spaces is **not** considered null/empty. `<td data-sort=" "></td>` will be sorted normally._

## Accessibility

Sortable is not very accessible in its raw form. It does not support screen readers, and it does not have any keyboard support. Including `sortable.a11y.min.js` in your project will add some basic accessibility features.

```html
<table class="sortable">
  ...
</table>
<link href="sortable.min.css" rel="stylesheet" />
<script src="sortable.min.js"></script>
<script src="sortable.a11y.min.js"></script>
```

By including the file the global function `enhanceSortableAccessibility` will automatically run through all existing `.sortable` tables, but you can also run it manually, like so:

```js
enhanceSortableAccessibility([table1, table2,...etc.])
```

The function adds an `aria-label` to each th, as well as `tabindex="0"` to each th in the thead of each table, making it possible to tab through the headers. It updates the `aria-label` depending on the direction.

if you want to import it instead this _should_ work: (I haven't tested it)

```ts
import { enhanceSortableAccessibility } from 'sortable-tablesort/enhanceSortableAccessibility'
enhanceSortableAccessibility([table1, table2,...etc.])
```

## Sort on load

If you wish to sort a table on load, I would recommend doing something like this:

```html
<table class="sortable">
  <thead>
    <tr>
      <th>Movie Name</th>
      <th id="movie-size">Size</th>
      <th>Release date</th>
    </tr>
  </thead>
  <tbody>
    ...
  </tbody>
</table>

<script>
  window.addEventListener('load', function () {
    const el = document.getElementById('movie-size')
    // without id:
    // const el = document.querySelector('.sortable th:first-child')
    // const el = document.querySelector('.sortable th:nth-child(2)')
    // const el = document.querySelectorAll('.sortable')[3].querySelector('th:nth-child(7)')
    // etc.
    if (el) {
      el.click()
    }
  })
</script>
```

Combine this with `<table class="sortable asc">` to reverse the sort order. Or do `el.click()` twice!

## Thank you...

- ...[Nikita Dunajevs](https://github.com/dunajevs) for the [ascending sort](#ascending-sort) idea!

- ...[wodny](https://github.com/wodny) for the [alternative sorting](#alternative-sorting) idea!

- ...[Nick Kocharhook](https://github.com/nk9) for the [colspan sorting](#colspanssort-on-specific-column) idea!

- ...[mxve](https://github.com/mxve) for the **nested elements** inside `th` fix!

- ...[Christian Petersson](https://github.com/Issen007) and [Abit Salihu](https://github.com/abitsalihu) for the [sort on load](#sort-on-load) example!

- ...[GazHay](https://github.com/gazhay) for the idea to sort multiple `<tbody />`!

- ...[Gordan Ratkovic](https://github.com/GordanRatkovic) for the [tiebreaker / secondary sort](#tiebreaker--secondary-sort) idea!

- ...[chatcoda](https://github.com/chatcoda) for the `<td></td>` / `<td>0</td>` sorting bug fix!

- ...[Christian Garbs](https://github.com/mmitch) for fixing the `dataset` bug!

- ...[Witold Baryluk](https://github.com/baryluk) for pointing out some inefficiencies!

- ...[Nick](https://github.com/data-handler) for raising a whole bunch of issues! 🤯️

- ...[James Pudson](https://github.com/nepsilon) for the [empty last](#emptynull-rows-always-last) suggestion, AND for finding the "`data-sort` ignored when empty" bug! 🥳️

- ...[Jojo-IO](https://github.com/Jojo-IO) for the finding the "`parseFloat` messes up time values" bug!

- ...[Steve Wirt](https://github.com/swirtSJW) for lifting the Accessibility issue! 🦾️

[![jsdelivr](https://data.jsdelivr.com/v1/package/gh/tofsjonas/sortable/badge)](https://www.jsdelivr.com/package/gh/tofsjonas/sortable)
