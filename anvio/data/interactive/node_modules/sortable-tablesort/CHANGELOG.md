<!-- markdownlint-disable MD024 -->

# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [3.2.2] - 2024-02-14

### Fixed

- Correct path for `sass` field in `package.json`.

## [3.2.1] - 2024-02-14

### Added

- Introduced `style` and `sass` fields in `package.json` to enhance compatibility with modern JavaScript tooling and environments.

## [3.2.0] - 2024-02-14

### Added

- Introduced `main` and `module` fields in `package.json` to enhance compatibility with modern JavaScript tooling and environments.

## [3.1.0] - 2023-11-16

### Changed

- Using `rollup` instead of `tsc` to generate javascript files.
- `src/` folder now included in npm package, so that `enhanceSortableAccessibility` can be used in TypeScript projects.

### Added

- Test minified files also.
- Added `focus` eventListener to `enhanceSortableAccessibility`, so that the `aria-label` is kept up to date.

## [3.0.0] - 2023-10-17

### Changed

- `aria-sort="ascending|descending"` used instead of `class="dir-d|dir-up"` to keep track of direction.

### Breaking Changes

- `class="dir-d|dir-up"` removed.

## [2.4.0] - 2023-10-17

### Added

- Simple accessibility introduced with **sortable.a11y.js**. `enhanceSortableAccessibility` adds `aria-label` to the headers of an array of tables.

## [2.3.2] - 2023-08-22

### Fixed

- `parseFloat()` turned time values like **12:11:11** and **12:23:56** into **12**, sorting them incorrectly. Enter: `Number()`! 🦸‍♂️️

## [2.3.1] - 2023-08-19

### Fixed

- `data-sort-alt` and `data-sort` were ignored if empty. No longer! 🦸‍♂️️

## [2.3.0] - 2023-08-13

### Added

- `class="n-last"` places empty cells always last, similar to what SQL does with ORDER BY foo NULLS LAST.

## [2.2.0] - 2023-06-30

### Changed

- `th` clicks only triggered in `thead`, not in `tbody` or `tfoot`
- `class="no-sort"` is now part of core JavaScript functionality, not CSS only like before

## [2.1.3] - 2023-03-24

### Fixed

- sortable-base.\* back in npm package

## [2.1.2] - 2023-03-24

### Changed

- Code quality bump

## [2.1.1] - 2023-03-24

### Fixed

- Bugfix tiebreaker column = 0

## [2.1.0] - 2023-03-23

### Added

- Tiebreaker/secondary sort

## [2.0.1] - 2023-03-21

### Fixed

- Bugfix dataset

## [2.0.0] - 2023-02-24

### Removed

- IE9 support dropped

## [1.80.1] - 2023-02-08

### Fixed

- Bugfix

## [1.80.0] - 2023-02-07

### Changed

- Typescript in src

## [1.70] - 2023-02-06

### Added

- First release

## Acknowledgments

This CHANGELOG.md was generated with the assistance of [ChatGPT by OpenAI](https://www.openai.com/research/chatgpt).

[3.2.2]: https://github.com/tofsjonas/sortable/releases/tag/3.2.2
[3.2.1]: https://github.com/tofsjonas/sortable/releases/tag/3.2.1
[3.2.0]: https://github.com/tofsjonas/sortable/releases/tag/3.2.0
[3.1.0]: https://github.com/tofsjonas/sortable/releases/tag/3.1.0
[3.0.0]: https://github.com/tofsjonas/sortable/releases/tag/3.0.0
[2.4.0]: https://github.com/tofsjonas/sortable/releases/tag/2.4.0
[2.3.2]: https://github.com/tofsjonas/sortable/releases/tag/2.3.2
[2.3.1]: https://github.com/tofsjonas/sortable/releases/tag/2.3.1
[2.3.0]: https://github.com/tofsjonas/sortable/releases/tag/2.3.0
[2.2.0]: https://github.com/tofsjonas/sortable/releases/tag/2.2.0
[2.1.3]: https://github.com/tofsjonas/sortable/releases/tag/2.1.3
[2.1.2]: https://github.com/tofsjonas/sortable/releases/tag/2.1.2
[2.1.1]: https://github.com/tofsjonas/sortable/releases/tag/2.1.1
[2.1.0]: https://github.com/tofsjonas/sortable/releases/tag/2.1.0
[2.0.1]: https://github.com/tofsjonas/sortable/releases/tag/2.0.1
[2.0.0]: https://github.com/tofsjonas/sortable/releases/tag/2.0.0
[1.80.1]: https://github.com/tofsjonas/sortable/releases/tag/1.80.1
[1.80.0]: https://github.com/tofsjonas/sortable/releases/tag/1.80.0
[1.70]: https://github.com/tofsjonas/sortable/releases/tag/1.70
