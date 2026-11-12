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
* [`essential-scripts-rs dryad-upload`↴](#essential-scripts-rs-dryad-upload)
* [`essential-scripts-rs run-gliph`↴](#essential-scripts-rs-run-gliph)
* [`essential-scripts-rs run-tcr-dist`↴](#essential-scripts-rs-run-tcr-dist)
* [`essential-scripts-rs tcr-match`↴](#essential-scripts-rs-tcr-match)

## `essential-scripts-rs`

A set of useful tools for data wrangling

**Usage:** `essential-scripts-rs [COMMAND]`

###### **Subcommands:**

* `aggregate-cell-ranger-tcr` — Aggregate CellRanger TCR output from multiple samples
* `split-sample-id` — Split sample ID into subject and condition from GLIPH2 output
* `split-cdr3-seq` — Split CDR3 sequences and genes if a semicolon is present
* `reformat-plate-reader-data` — 
* `copy-cell-ranger-outs` — Copy selected outputs from Cell Ranger pipestances
* `score-tcr-alignments` — Score a GLIPH2 output file with a TCR alignment pipeline
* `geo-fastq` — Match FastQ files by lane and sample, and compute MD5 checksums
* `run-enrichr` — Run Enrichr via API interface. Requires re-installing with --features enrichment
* `dryad-upload` — Upload files in a directory to a Dryad dataset DOI
* `run-gliph` — A reimplementation of the GLIPH algorithms
* `run-tcr-dist` — Cluster chains using tcrdist distances and score the clusters
* `tcr-match` — A reimplementation of TCRMatch



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

Score a GLIPH2 output file with a TCR alignment pipeline

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

  Default value: `8`
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



## `essential-scripts-rs dryad-upload`

Upload files in a directory to a Dryad dataset DOI

**Usage:** `essential-scripts-rs dryad-upload [OPTIONS] --client-id <CLIENT_ID> --client-secret <CLIENT_SECRET> --doi <DOI> [DIRECTORY]`

###### **Arguments:**

* `<DIRECTORY>`

  Default value: `.`

###### **Options:**

* `--client-id <CLIENT_ID>`
* `--client-secret <CLIENT_SECRET>`
* `--doi <DOI>`
* `--progress` — Show progress bar
* `--no-progress` — Hide progress bar



## `essential-scripts-rs run-gliph`

A reimplementation of the GLIPH algorithms

`run-gliph` reads a chain table with required `subject:condition` and `TcRb` columns (and optional `TcRa`, `V`, `J`, and `counts` columns). Each `subject:condition` value is split at the first colon during processing. An optional HLA TSV can be supplied with `--hla-file`; it has a header row and each row contains the subject followed by left/right allele columns for `HLA-DPA1`, `HLA-DPB1`, `HLA-DQA1`, `HLA-DQB1`, `HLA-DRB1`, and `HLA-DRB3/4/5` (`NA` or empty means missing). Output has eight HLA columns (DRB3, DRB4 and DRB5 are split by allele name); distinct alleles are joined with `/`, e.g. `DPA1*01:03/DPA1*02:01`, and all of them are used for scoring. Alleles whose enrichment p-value within a cluster is at or below `--hla-threshold` (default 0.05) are suffixed with `!` in the output. When omitted, HLA output columns are left blank. It writes one row per chain that belongs to a GLIPH cluster. The output is a TSV with the fixed column order `index`, `pattern`, `number_subject`, `number_unique_cdr3`, `TcRb`, `Fisher_score`, `Freq`, `subject`, `condition`, `final_score`, `hla_score`, `vb_score`, `expansion_score`, `length_score`, `cluster_size_score`, `type`, `ulTcRb`, `V`, `J`, `TcRa`, and the eight HLA columns above. When the input has CDR3α sequences, a `pattern_alpha` column follows with the motif shared by each cluster's CDR3αs (found as for `run-tcr-dist`), then `VA` and `JA` when alpha genes are given.

**Usage:** `essential-scripts-rs run-gliph [OPTIONS] <CHAINS_FILE> <OUTPUT_FILE>`

###### **Arguments:**

* `<CHAINS_FILE>` — TSV containing subject:condition, TcRb, and optional TcRa/V/J/counts columns
* `<OUTPUT_FILE>` — Output TSV file path; use '-' for stdout

###### **Options:**

* `--hla-file <HLA_FILE>` — TSV with a header row, keyed by subject, followed by left/right allele columns for DPA1, DPB1, DQA1, DQB1, DRB1 and DRB3/4/5
* `--method <METHOD>` — GLIPH algorithm to use

  Default value: `gliph2`

  Possible values: `gliph1`, `gliph2`

* `--sim-depth <SIM_DEPTH>` — Number of simulations to run for enrichment scores. Minimum score is 1/sim-depth.

  Default value: `1000`
* `--hla-threshold <HLA_THRESHOLD>` — Mark (with '!') HLA alleles whose cluster enrichment p-value is at or below this threshold

  Default value: `0.05`
* `--reference <REFERENCE>` — Reference dataset to use for global similarity scoring

  Default value: `human-v2cd4cd8`

  Possible values: `human-v1cd4`, `human-v1cd8`, `human-v1cd4cd8`, `human-v2cd4`, `human-v2cd8`, `human-v2cd4cd8`, `mouse-v1cd4`, `mouse-v1cd8`, `mouse-v1cd4cd8`, `mouse-v2cd4`, `mouse-v2cd8`, `mouse-v2cd4cd8`

* `--reference-cache-dir <REFERENCE_CACHE_DIR>`
* `--refresh-reference-cache`
* `--cores <CORES>` — Number of cores to use. 0 = all available; negative values subtract from available cores

  Default value: `8`
* `--sub-matrix <SUB_MATRIX>` — Substitution matrix to use for creating groups

  Possible values: `blosum62`, `tcr`




## `essential-scripts-rs run-tcr-dist`

Cluster chains using tcrdist distances and score the clusters

Output matches `run-gliph`: `pattern` is the CDR3β motif of each cluster and, when the input has CDR3α sequences, `pattern_alpha` is the CDR3α motif.

**Usage:** `essential-scripts-rs run-tcr-dist [OPTIONS] <CHAINS_FILE> <OUTPUT_FILE>`

###### **Arguments:**

* `<CHAINS_FILE>` — TSV containing subject:condition, TcRb, and optional TcRa/V/J/counts columns
* `<OUTPUT_FILE>` — Output TSV file path; use '-' for stdout

###### **Options:**

* `--hla-file <HLA_FILE>` — TSV with a header row, keyed by subject, followed by left/right allele columns for DPA1, DPB1, DQA1, DQB1, DRB1 and DRB3/4/5
* `--method <METHOD>` — GLIPH algorithm to use

  Default value: `gliph2`

  Possible values: `gliph1`, `gliph2`

* `--sim-depth <SIM_DEPTH>` — Number of simulations to run for enrichment scores. Minimum score is 1/sim-depth.

  Default value: `1000`
* `--hla-threshold <HLA_THRESHOLD>` — Mark (with '!') HLA alleles whose cluster enrichment p-value is at or below this threshold

  Default value: `0.05`
* `--reference <REFERENCE>` — Reference dataset to use for global similarity scoring

  Default value: `human-v2cd4cd8`

  Possible values: `human-v1cd4`, `human-v1cd8`, `human-v1cd4cd8`, `human-v2cd4`, `human-v2cd8`, `human-v2cd4cd8`, `mouse-v1cd4`, `mouse-v1cd8`, `mouse-v1cd4cd8`, `mouse-v2cd4`, `mouse-v2cd8`, `mouse-v2cd4cd8`

* `--reference-cache-dir <REFERENCE_CACHE_DIR>`
* `--refresh-reference-cache`
* `--cores <CORES>` — Number of cores to use. 0 = all available; negative values subtract from available cores

  Default value: `8`
* `--sub-matrix <SUB_MATRIX>` — Substitution matrix to use for creating groups

  Possible values: `blosum62`, `tcr`

* `--distance-mode <DISTANCE_MODE>` — Chains and CDRs used for the distance. V genes without an allele are looked up as `*01` [default: alpha-beta if CDR3α is present, else beta]

  Possible values:
  - `beta`:
    CDR3β only
  - `alpha-beta`:
    CDR3α + CDR3β
  - `beta-full`:
    CDR3β + germline CDR1β, CDR2β, and CDR2.5β (from TRBV)
  - `alpha-beta-full`:
    CDR3α + CDR3β + germline CDR1, CDR2, and CDR2.5 of both chains (from TRAV and TRBV)

* `--organism <ORGANISM>` — Organism of the germline V genes for the *-full distance modes

  Default value: `human`

  Possible values: `human`, `mouse`

* `--cluster-method <CLUSTER_METHOD>` — How records within each CDR3-edit threshold are grouped into clusters

  Default value: `shared-neighbor`

  Possible values:
  - `shared-neighbor`:
    Records within the threshold are linked only when their neighborhoods overlap enough (Jaccard index of their closed neighborhoods at least `--min-neighbor-overlap`); clusters are the connected components of the kept links. A single link between two dense groups does not join them, so chaining is limited by the same rule everywhere instead of a cluster size cap
  - `single-linkage`:
    Connected components (single linkage) found with a union-find straight from the distance matrix. One pass covers every threshold and no graph is built
  - `connected-components`:
    Connected components of an igraph graph built for each threshold. Same clusters as `single-linkage` but slower and memory-hungry at large thresholds
  - `louvain`:
    Louvain communities (modularity optimization, see `--resolution`) of the graph linking records within the threshold. Records are first split into single-linkage components at the largest threshold and each component is clustered on its own, straight from the distance matrix, so only one component's edges are in memory at a time

* `--min-neighbor-overlap <MIN_NEIGHBOR_OVERLAP>` — For shared-neighbor clustering: link two records within the threshold only if the Jaccard overlap of their neighborhoods (each record plus its neighbors) is at least this. 0 is single linkage; higher values split chains more readily

  Default value: `0.4`
* `--resolution <RESOLUTION>` — For Louvain clustering: modularity resolution. 1 is standard modularity; higher values give smaller, tighter clusters and lower values larger ones

  Default value: `1`



## `essential-scripts-rs tcr-match`

A reimplementation of TCRMatch

**Usage:** `essential-scripts-rs tcr-match [OPTIONS] <INPUT_FILE>`

###### **Arguments:**

* `<INPUT_FILE>` — Input file to process

###### **Options:**

* `-d`, `--database <DATABASE>`

  Default value: `IEDB`

  Possible values: `iedb`, `cedar`

* `-t`, `--threads <THREADS>` — Number of cores to use. 0 = all available; negative values subtract from available cores

  Default value: `8`
* `-s`, `--threshold <THRESHOLD>`

  Default value: `0.97`
* `-a`, `--airr`
* `-k`, `--keep-flanking`
* `-r`, `--trust4`
* `-m`, `--max-ram-gb <MAX_RAM_GB>`



<hr/>

<small><i>
    This document was generated automatically by
    <a href="https://crates.io/crates/clap-markdown"><code>clap-markdown</code></a>.
</i></small>
