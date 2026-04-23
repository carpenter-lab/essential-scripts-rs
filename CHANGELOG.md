# Changelog

## [0.3.0](https://github.com/carpenter-lab/essential-scripts-rs/compare/v0.2.0...v0.3.0) (2026-04-23)


### Features

* add copy_cellranger_outs command and improve error handling ([#17](https://github.com/carpenter-lab/essential-scripts-rs/issues/17)) ([d0f2f25](https://github.com/carpenter-lab/essential-scripts-rs/commit/d0f2f25a215686cc9ec6b108727a675537c9283a))
* **cli:** add feature flag for terminal sizing in help ([30ec41e](https://github.com/carpenter-lab/essential-scripts-rs/commit/30ec41e68a9dd9272f1fe328d17a083361bcd266))
* **cli:** update output file argument to have a default value of "-" for better usability ([7639419](https://github.com/carpenter-lab/essential-scripts-rs/commit/7639419585a6f348a6c983e1569737f7fe288ebe))
* **split:** enhance CDR3 sequence splitting with gene schema resolution and optional grouping ([8899ce2](https://github.com/carpenter-lab/essential-scripts-rs/commit/8899ce2e87300eb92819307495a9199666880c03))


### Bug Fixes

* **dataframe:** enable keeping nulls in explode options for improved data handling ([dcb4fe6](https://github.com/carpenter-lab/essential-scripts-rs/commit/dcb4fe698afc1af2bc97112691f94f0d0843a191))
* update explode options to handle nulls and empty values with updated polars crate ([b63f9c2](https://github.com/carpenter-lab/essential-scripts-rs/commit/b63f9c2c42d67898ff306f685faa4cf096dcf61d))

## [0.2.0](https://github.com/carpenter-lab/essential-scripts-rs/compare/v0.1.0...v0.2.0) (2026-02-06)


### Features

* add command to copy cellranger outputs ([#7](https://github.com/carpenter-lab/essential-scripts-rs/issues/7)) ([8d7b91d](https://github.com/carpenter-lab/essential-scripts-rs/commit/8d7b91de3f7059a591cd1dbb28e569c315b18275))
* implement TCR alignment scoring command with Parasail integration ([#5](https://github.com/carpenter-lab/essential-scripts-rs/issues/5)) ([773919c](https://github.com/carpenter-lab/essential-scripts-rs/commit/773919c2a6444a51e565a03be0d7336979994ce3))
* **plate-reader:** enhance data processing with stride-based flattening and improved output handling ([5ca621b](https://github.com/carpenter-lab/essential-scripts-rs/commit/5ca621bcbff11b67885320241ffef73f3062dc6f))

## 0.1.0 (2025-12-10)


### Features

* add a command to copy results from a Cell Ranger pipestance ([#2](https://github.com/carpenter-lab/essential-scripts-rs/issues/2)) ([8a02ece](https://github.com/carpenter-lab/essential-scripts-rs/commit/8a02ece582e317559bec0b82068fe5d5b3175d74))
* **aggregate:** enhance TCR aggregation with alpha chain support and error handling ([c0f2754](https://github.com/carpenter-lab/essential-scripts-rs/commit/c0f27548bfc9634c47dc8492004b17525f735498))
* **aggregate:** filter out null and empty values in `CDR3b` column during processing ([5eef3fa](https://github.com/carpenter-lab/essential-scripts-rs/commit/5eef3fa1dd8542595397693341d95bff8fa93d72))
* **cli:** add customizable column name for splitting in `split_sample_id` ([2e90042](https://github.com/carpenter-lab/essential-scripts-rs/commit/2e900426ddeb095f1a8202654fb45d19a27b5a01))
* **cli:** setup initial CLI commands with Clap ([fd005ec](https://github.com/carpenter-lab/essential-scripts-rs/commit/fd005ecf25a901b5fa6e2489237420bcefcad7c7))
* **io:** add functionality to write DataFrame or LazyFrame to CSV or stdout ([fd005ec](https://github.com/carpenter-lab/essential-scripts-rs/commit/fd005ecf25a901b5fa6e2489237420bcefcad7c7))
* **io:** add support for multi-format output and simplify error handling ([e8df7b8](https://github.com/carpenter-lab/essential-scripts-rs/commit/e8df7b86dd49fe3e4ea3450c8d4d78c27e25518a))
* **plate_reader:** add placeholder handling for plate reader data reformatting ([fd005ec](https://github.com/carpenter-lab/essential-scripts-rs/commit/fd005ecf25a901b5fa6e2489237420bcefcad7c7))


### Bug Fixes

* ComputeError(ErrString("CSV format does not support nested data")) when running `split_sample_id` ([2e90042](https://github.com/carpenter-lab/essential-scripts-rs/commit/2e900426ddeb095f1a8202654fb45d19a27b5a01))


### Performance Improvements

* **aggregate:** improve performance with streaming and maintain input order ([5eef3fa](https://github.com/carpenter-lab/essential-scripts-rs/commit/5eef3fa1dd8542595397693341d95bff8fa93d72))
