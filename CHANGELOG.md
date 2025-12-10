# Changelog

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
