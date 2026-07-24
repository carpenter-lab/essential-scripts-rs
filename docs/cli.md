# CLI Documentation

This document contains the help content for the `essential-scripts-rs` command-line program.

**Command Overview:**

* [`essential-scripts-rs`↴](#essential-scripts-rs)
* [`essential-scripts-rs aggregate-cell-ranger-tcr`↴](#essential-scripts-rs-aggregate-cell-ranger-tcr)
* [`essential-scripts-rs split-sample-id`↴](#essential-scripts-rs-split-sample-id)
* [`essential-scripts-rs split-cdr3-seq`↴](#essential-scripts-rs-split-cdr3-seq)
* [`essential-scripts-rs reformat-plate-reader-data`↴](#essential-scripts-rs-reformat-plate-reader-data)
* [`essential-scripts-rs copy-cell-ranger-outs`↴](#essential-scripts-rs-copy-cell-ranger-outs)
* [`essential-scripts-rs score-tcr-alignments`↴](#essential-scripts-rs-score-tcr-alignments)
* [`essential-scripts-rs geo-fastq`↴](#essential-scripts-rs-geo-fastq)
* [`essential-scripts-rs run-enrichr`↴](#essential-scripts-rs-run-enrichr)

## `essential-scripts-rs`

A set of useful tools for data wrangling

**Usage:** `essential-scripts-rs [COMMAND]`

###### **Subcommands:**

* `aggregate-cell-ranger-tcr` — A set of useful tools for data wrangling
* `split-sample-id` — A set of useful tools for data wrangling
* `split-cdr3-seq` — A set of useful tools for data wrangling
* `reformat-plate-reader-data` — 
* `copy-cell-ranger-outs` — A set of useful tools for data wrangling
* `score-tcr-alignments` — Score a GLIPH2 output file with a TCR alignment pipeline. Requires re-installing with --features tcr
* `geo-fastq` — Match FastQ files by lane and sample, and compute MD5 checksums
* `run-enrichr` — Run Enrichr via API interface. Requires re-installing with --features enrichment



## `essential-scripts-rs aggregate-cell-ranger-tcr`

Aggregate CellRanger TCR output from multiple samples

Parse a set of input files with the Cell Ranger TCR format (`filtered_contig_annotations.csv`) and aggregate them into a single output file. The output will contain one row per unique combination of sample, barcode, and TCR chain, with the corresponding CDR3 sequences and gene segments. Optionally, the alpha chain can be retained in the output. The internal `sample_id` column is used, so there is no need for unique filenames if running on direct outputs from Cell Ranger

**Usage:** `essential-scripts-rs aggregate-cell-ranger-tcr [OPTIONS] <INPUT_FILES>... <OUTPUT_FILE>`

###### **Arguments:**

* `<INPUT_FILES>` — Input CSV files to process
* `<OUTPUT_FILE>` — Output file path

###### **Options:**

* `-k`, `--keep-alpha` — Keep alpha chain in output [default: false]

  Default value: `false`



## `essential-scripts-rs split-sample-id`

Split sample ID into subject and condition from GLIPH2 output

Splits the `subject:condition` column (the default) into two columns. Can be used on any file with a separator to split by however, these columns will always be called subject and condition

**Usage:** `essential-scripts-rs split-sample-id [OPTIONS] <INPUT_FILE>... <OUTPUT_FILE>`

###### **Arguments:**

* `<INPUT_FILE>` — Input CSV file to process
* `<OUTPUT_FILE>` — Output file path

###### **Options:**

* `-c`, `--column-name <COLUMN_NAME>` — Column to split

  Default value: `subject:condition`



## `essential-scripts-rs split-cdr3-seq`

Split CDR3 sequences and genes if a semicolon is present

This potentially splits each row into up to 4 new rows. In the output from scRepertoire, the rows each belong to a single cell barcode. In the case when a cell has 2 detectable CDR3 beta or alpha sequences, the resulting CDR3 and V/J columns are concatinated with a ";". For downstream applications, this results in treating this chimeric sequence as a real, biological sequence. This tool will expand this into up to 4 different pairs of chains. "beta_1;beta2" and "alpha_1;alpha_2" will be split into 4 rows each containing a single alpha and beta chain.

The TCR columns must be named CDR3a and CDR3b. Requires either CTgeneA/CTgeneB columns or TRAV/TRAJ/TRBV/TRBJ columns for the TCR genes. If group columns are not provided, each input row is treated as its own group and alpha/beta splits are anchored to original rows.

**Usage:** `essential-scripts-rs split-cdr3-seq [OPTIONS] <INPUT_FILE> [OUTPUT_FILE]`

###### **Arguments:**

* `<INPUT_FILE>` — Input CSV file to process
* `<OUTPUT_FILE>` — Output file path

  Default value: `-`

###### **Options:**

* `-g`, `--group <GROUP>` — Optional columns to group by; if omitted, each input row is treated as its own group



## `essential-scripts-rs reformat-plate-reader-data`

**Usage:** `essential-scripts-rs reformat-plate-reader-data <INPUT_FILE> <OUTPUT_PATH>`

###### **Arguments:**

* `<INPUT_FILE>` — Input Excel file to process
* `<OUTPUT_PATH>` — Output directory path. Will create one CSV per sheet.



## `essential-scripts-rs copy-cell-ranger-outs`

Copy selected outputs from Cell Ranger pipestances

Scans through a directory searching for Cell Ranger Pipestances (subdirectories containing a `*.mri.tgz` marker file). For each pipestance, copies selected outputs (H5, MEX, VDJ) into a destination directory. Those files are renamed from their default (`sample_filtered_feature_bc_matrix.h5`) to `<sample>.h5` and stored as a flat directory.

**Usage:** `essential-scripts-rs copy-cell-ranger-outs [OPTIONS] --base-path <BASE_PATH>`

###### **Options:**

* `-b`, `--base-path <BASE_PATH>` — Base directory containing Cell Ranger pipestances
* `-d`, `--dest <DEST>` — Destination directory to copy to
* `--h5` — Copy filtered H5 matrix as <sample>.h5

  Default value: `false`
* `--mex` — Copy filtered MEX directory into <sample>/

  Default value: `false`
* `--vdj` — Copy VDJ annotations as <sample>.csv

  Default value: `false`
* `--check` — Check each pipestance for presence of a *.mri.tgz marker and print results. No copying is performed.

  Default value: `false`



## `essential-scripts-rs score-tcr-alignments`

Score a GLIPH2 output file with a TCR alignment pipeline. Requires re-installing with --features tcr

**Usage:** `essential-scripts-rs score-tcr-alignments [OPTIONS] [INPUT_FILE] [OUTPUT_FILE]`

###### **Arguments:**

* `<INPUT_FILE>` — Input CSV file to process
* `<OUTPUT_FILE>` — Output file path

###### **Options:**

* `-r`, `--replicates <REPLICATES>`

  Default value: `1000`



## `essential-scripts-rs geo-fastq`

Match FastQ files by lane and sample, and compute MD5 checksums

**Usage:** `essential-scripts-rs geo-fastq [OPTIONS] [INPUT_DIRECTORIES]...`

###### **Arguments:**

* `<INPUT_DIRECTORIES>` — Directories containing fastq.gz files

###### **Options:**

* `--paired-output <PAIRED_OUTPUT>` — Output file for paired files by lane and sample [default: stdout]
* `--sample-output <SAMPLE_OUTPUT>` — Output file for all files per sample [default: stdout]
* `--md5-output <MD5_OUTPUT>` — Output file for file paths with MD5 checksums [default: stdout]
* `--parallel-md5` — Compute MD5s in parallel

  Default value: `false`
* `--threads <THREADS>` — Number of threads to use for parallel MD5. 0 = all available cores

  Default value: `0`
* `--progress` — Show progress bar
* `--no-progress` — Hide progress bar



## `essential-scripts-rs run-enrichr`

Run Enrichr via API interface. Requires re-installing with --features enrichment

**Usage:** `essential-scripts-rs run-enrichr [OPTIONS] --library <LIBRARY> --gene-list <GENE_LIST> <OUTPUT_FILE> [OUTPUT_PLOT]...`

###### **Arguments:**

* `<OUTPUT_FILE>` — Output file path
* `<OUTPUT_PLOT>` — Output file paths

###### **Options:**

* `-l`, `--library <LIBRARY>` — Enrichr Library to use

  Possible values: `reactome-pathways2024`, `reactome`, `bio-carta2016`, `wiki-pathways2024-human`, `go-biological-process`

* `-g`, `--gene-list <GENE_LIST>` — Input gene list to process. One gene per line
* `-b`, `--background <BACKGROUND>` — Input gene list to process as background. One gene per line



<hr/>

<small><i>
    This document was generated automatically by
    <a href="https://crates.io/crates/clap-markdown"><code>clap-markdown</code></a>.
</i></small>
