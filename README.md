# CLI tools to make data processing and preparation easier

[![Rust](https://github.com/carpenter-lab/essential-scripts-rs/actions/workflows/rust.yml/badge.svg)](https://github.com/carpenter-lab/essential-scripts-rs/actions/workflows/rust.yml)
[![codecov](https://codecov.io/gh/carpenter-lab/essential-scripts-rs/graph/badge.svg?token=QAHBPDKTW5)](https://codecov.io/gh/carpenter-lab/essential-scripts-rs)

## Installation

To install the latest version of the package, run:

```bash
cargo install --locked --git https://github.com/carpenter-lab/essential-scripts-rs
```

or

```bash
cargo install --locked --git ssh://git@github.com/carpenter-lab/essential-scripts-rs.git
```

If you need the cargo command (i.e. `cargo install`), you can install it by following the
instructions [to install Rust](https://rust-lang.org/tools/install/).
The package will be installed in `~/.cargo/bin` and will take between 5 and 15 minutes to
compile.

## The Tools

All tools have a `--help` option which explains how to use the tool.
This only explains how to access the help option.
This can be run by using any of the following commands:

```
essential-scripts-rs help
# or
essential-scripts-rs help << tool name >>
# or
essential-scripts-rs << tool name >> --help
```

To run any tool, use the following setup:

```
essential-scripts-rs << tool name >>  << Positinal Arguments >>  << Options >>
```

Replace `<< tool name >>` with the actual tool name (without angled brackets `<`/`>`).
