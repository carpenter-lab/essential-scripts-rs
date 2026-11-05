# CLI tools to make data processing and preparation easier

[![Rust](https://github.com/carpenter-lab/essential-scripts-rs/actions/workflows/rust.yml/badge.svg)](https://github.com/carpenter-lab/essential-scripts-rs/actions/workflows/rust.yml)
[![codecov](https://codecov.io/gh/carpenter-lab/essential-scripts-rs/graph/badge.svg?token=QAHBPDKTW5)](https://codecov.io/gh/carpenter-lab/essential-scripts-rs)

## Installation

To install the latest version of the package, run:

```bash
RUSTFLAGS="-C target-cpu=native" cargo install \
  --locked \
  --git https://github.com/carpenter-lab/essential-scripts-rs
```

or

```bash
RUSTFLAGS="-C target-cpu=native" cargo install \
  --locked \
  --git ssh://git@github.com/carpenter-lab/essential-scripts-rs.git
```

If you need the cargo command (i.e. `cargo install`), you can install it by following the
instructions [to install Rust](https://rust-lang.org/tools/install/). The package will be installed in `~/.cargo/bin`
and will take between 5 and 15 minutes to compile.

### Troubleshooting errors

The `RUSTFLAGS` enviornment variable is optional but recommended.
This allows the SIMD code for TCRdist to run at the fastest speed possible.
This could theoretically cause errors during installation depending on the system.
If you encounter any errors, try removing the `RUSTFLAGS` environment variable.

The `--locked` option ensures the dependencies used during installation are the same ones used during testing.
If a crate is yanked, this will prevent the installation from succeeding.
Try without the `--locked` option if you encounter this error.

## Documentation

The CLI documentation can be found in the docs directory [here](docs/cli.md).

## The Tools

All tools have a `--help` option which explains how to use the tool. This only explains how to access the help option.
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
