# AGENTS.md

Guidance for coding agents working in `essential-scripts-rs`: a single Rust (edition 2024) binary
with data-wrangling subcommands for immunology/bioinformatics work (CellRanger, GLIPH2, TCR
alignment and clustering, Enrichr, GEO and Dryad submission, plate-reader files).

## Layout

- `src/main.rs`: the clap `Cli`, the `Commands` enum (each variant is a `#[command(flatten)]` of a
  module's own `Commands`) and the dispatch to `handle_command`. The hidden `--markdown-help` flag
  writes `docs/cli.md` via `src/docs.rs`.
- `src/lib.rs`: declares each module.
- One module per tool family. Larger ones split into `mod.rs` (clap `Commands` enum plus
  `handle_command`) and `core.rs` (the logic); `tcr_cluster/` has deeper submodules.
- Shared modules:
  - `io.rs`: reading and writing CSV/TSV/XLSX tables as polars frames.
  - `cache.rs`: downloads a file once into a cache directory, writing to a temporary file and
    renaming it, with optional size and MD5 checks (`cache::cached_file`). Used for every
    reference database.
  - `cores.rs`: parses thread-count arguments (`thread_parser`, `thread_default`: `0` means all
    cores, negative numbers subtract from the available cores) and builds rayon pools
    (`build_pool`).
  - `progress.rs`: `indicatif` progress bars plus the shared `ProgressArg` flag that turns them
    on or off.
- `tests/cli.rs`: end-to-end CLI tests with `assert_cmd`. `tests/data/`: fixtures.
- `docs/` is generated and git-ignored; don't commit it.

## Feature flags

| Feature      | Enables                                  | Extra system libraries                         |
|--------------|------------------------------------------|------------------------------------------------|
| `tcr`        | `tcr_align` (`parasail-rs`, `rand`)      | `libparasail-dev`, `parasail`, `cmake`, `clang` |
| `enrichment` | `enrich` (plotters, cairo)               | `libglib2.0-dev`, `libcairo2-dev`, `libfontconfig1-dev`, `libfreetype6-dev` |

Feature-gated modules always compile and register their subcommand. Without the feature, the
subcommand's `about` text ends with `Requires re-installing with --features <name>` and
`handle_command` returns an error. Keep that pattern (`#[cfg(feature = ...)]` /
`#[cfg(not(feature = ...))]` pairs) when adding gated code. Code must build and pass tests with
**both** default features and `--all-features`.

Each gated module has a test for the disabled state (for example `disabled_tests` in
`src/enrich/mod.rs`) that checks the error kind is `MissingSubcommand`. Add the same test for new
gated code.

## Commands

Always pass `--locked`; CI fails if `Cargo.lock` doesn't match `Cargo.toml`.

```bash
cargo build --locked --all-features
cargo nextest run --locked --all-features --all-targets   # CI uses nextest
cargo nextest run --locked --all-targets                  # default features too
cargo clippy --locked --all-features -- -D warnings
cargo clippy --locked -- -D warnings
cargo fmt --all
```

`cargo test` works if nextest isn't installed. If `Cargo.lock` is stale (typically after merging
`main`), fix it with `cargo metadata --all-features --format-version=1 > /dev/null` and commit the
lockfile on its own as `chore: update Cargo.lock`.

Other useful commands:

```bash
cargo llvm-cov nextest --locked --all-features --all-targets   # coverage, as in CI
cargo build --locked --all-features --profile profiler         # release build with debug symbols
cargo run --locked -- --markdown-help                          # regenerate docs/cli.md
```

Building with `--all-features` needs the system libraries in the feature table. On Debian or
Ubuntu, CI installs them with `apt-get install libparasail-dev parasail cmake clang
libglib2.0-dev libcairo2-dev libfontconfig1-dev libfreetype6-dev`.

## Code conventions

- **Errors:** return `anyhow::Result` and add context with `.context(...)` / `bail!`.
  `thiserror` is only for small typed errors (see `cores.rs`). No `unwrap()` on user input or I/O
  outside tests. Validation clap can't express goes in a `validate` method returning
  `clap::Error::raw(ErrorKind::ValueValidation, ...)`, so it reads like clap's own errors (see
  `CommonArgs::validate` in `src/tcr_cluster/cli.rs`).
- **Tables:** use `polars` (lazy where practical) and the readers/writers in `src/io.rs`, not
  hand-rolled file handling. Output path `-` means stdout where supported.
- **Threads and progress:** take thread counts with `thread_parser()` / `thread_default()`,
  resolve them with `cores::process_cores`, and run parallel work inside the pool from
  `cores::build_pool`. Flatten `ProgressArg` into any long-running subcommand.
- **CLI text:** clap derive. Each subcommand's `///` doc comment is its full user documentation
  and ends up in `docs/cli.md`, so keep it accurate and user-facing; give every argument a
  `help = "..."`. Clap derives the command name from the enum variant (`AggregateCellRangerTCR`
  → `aggregate-cell-ranger-tcr`), so renaming a variant renames the user-facing command.
- **Dependencies:** the shared-dependency feature pins in `Cargo.toml` (e.g. `smallvec`,
  `simd-adler32`) keep default and all-features builds sharing compiled crates. Don't remove
  them.
- Run `cargo fmt` before committing; Clippy runs with warnings as errors.

## Adding a subcommand

1. Create `src/<name>/mod.rs` with a `#[derive(Subcommand)] pub enum Commands` and
   `pub fn handle_command(cmd: Commands) -> anyhow::Result<()>`; put logic in `core.rs`.
2. Add `pub mod <name>;` to `src/lib.rs`.
3. Add a `#[command(flatten)]` variant to `Commands` in `src/main.rs` and a match arm in
   `main_helper`.
4. Add unit tests in a `#[cfg(test)] mod tests` block and add the subcommand to the
   `command_cases` template in `tests/cli.rs`.

## TCR code invariants

- `run-gliph` and `run-tcr-dist`: every algorithm implements `ClusterFinder` and goes through
  the shared pipeline in `tcr_cluster/pipeline.rs` for scoring and output. Don't copy scoring or
  output code into an algorithm.
- Both commands write the same TSV columns, listed in the `RunGliph` doc comment. A change to
  `OutputRow` changes both commands; update that doc comment in the same change.
- `tcr-match`: the fast prefilter (`k3_sum_fast`) must stay exact-safe. `PREFILTER_MARGIN`
  guarantees it never drops a pair that the reference kernel `k3_sum` would keep; the
  equivalence tests must keep passing.
- `score-tcr-alignments`: each group's `SmallRng` is seeded from `base_seed + i`, so results
  don't depend on how rayon schedules groups.
- `tcr_cluster/tcrdist/distance.rs` has AVX2 and AVX-512 kernels behind `cfg(target_feature)`.
  Keep the scalar fallback and the test that they agree.

## Performance work

- Profile with the `profiler` build profile (release settings plus debug symbols). Users install
  with `RUSTFLAGS="-C target-cpu=native"`, so measure with it too, but never rely on it for
  correctness: SIMD code needs a portable fallback.
- A speed-up must not change output. Add a test that compares the fast path with the reference
  path, as the recent performance PRs did.
- Parallel results must not depend on scheduling: seed randomness per item rather than per
  thread, and keep output order stable.
- In hot loops, precompute per-query data once and reuse it, and avoid allocating per pair or
  per row.

## Network and credentials

- Downloads go through `cache::cached_file`. Pin URLs to a commit or record version, as the
  existing constants do, and check size and MD5 when they are known.
- Dryad credentials (`--client-id` / `--client-secret`, or `DRYAD_CLIENT_ID` / `DRYAD_SECRET`)
  must never be printed, logged, written to fixtures or committed.
- Tests never call live services (Enrichr, Dryad, or the reference-data hosts). See *Tests* for
  how to avoid them.

## Tests

- Unit tests sit at the bottom of each file in `#[cfg(test)] mod tests`, mostly `rstest` with
  `#[case]`; `rstest_reuse` templates for repeated cases.
- Compare dataframes with `polars-testing`, other values with `pretty_assertions`.
- Write temporary files with `tempfile`/`temp_testdir` (helpers in `src/test_helpers.rs`), never
  into the repo. `*.txt` and `*.tsv` are git-ignored, so new fixtures with those extensions must be
  force-added.
- `tests/cli.rs` runs every subcommand in the `command_cases` template and has extra end-to-end
  cases for `run-tcr-dist`. `dryad-upload` has no CLI test yet.
- The nextest `ci` profile (`.config/nextest.toml`) turns off fail-fast and writes JUnit XML. Use
  `--profile ci` to reproduce a CI run.
- Testing without the network:
  - GLIPH references: copy the fixtures in `tests/data/gliph_reference/` into a `tempdir` and pass
    it as the cache directory.
  - `tcr-match`: implement the `GetData` trait with a test database (see `TestDatabase` in
    `tcrmatch/helper.rs`).
  - `geo-fastq`: the traits in `geo_submission/traits.rs` let tests substitute MD5 results.
  - Enrichr and Dryad: unit-test request building, response parsing and error messages only.

## Commits, CI and releases

- Conventional Commits, enforced by commitizen in pre-commit: `feat:`, `fix:`, `perf:`,
  `refactor:`, `test:`, `ci:`, `build:`, `chore:`, with an optional scope such as
  `feat(tcr-align): ...`. The type drives the release-please version bump.
- **Do not edit `CHANGELOG.md` or the `version` in `Cargo.toml`.** release-please owns both: it
  keeps a `chore(main): release X.Y.Z` PR open, and merging it tags `vX.Y.Z` and moves the
  `stable` tag. `release-docs.yml` then attaches `docs/cli.md` to that release.
- PRs target `main`. On PRs and the merge queue:
  - `rust.yml` checks `Cargo.lock`, then runs nextest under `cargo llvm-cov` with all features
    and with default features (one shared build, thanks to the feature pins), uploads coverage
    to Codecov, and builds `docs/cli.md` from the default-features binary.
  - `clippy.yml` runs Clippy with `-Dwarnings` for both feature sets and uploads the results to
    code scanning. Its steps are `continue-on-error`, so the job can pass with warnings; fix
    them anyway.
- Dependabot opens weekly Cargo updates, with minor and patch updates grouped into one PR.
