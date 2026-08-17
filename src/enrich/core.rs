use crate::enrich::api::{EnrichrAPI, EnrichrAPITrait};
use crate::io::WriteToCsvOrStdout;
use cairo;
use plotters::coord::Shift;
use plotters::prelude::*;
use plotters_cairo::CairoBackend;
use polars::polars_utils::itertools::Itertools;
use polars::prelude::*;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::error::Error;
use std::fs;
use std::path::{Path, PathBuf};
use tokio;

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct EnrichrResult {
    rank: Vec<i32>,
    term: Vec<String>,
    p_value: Vec<f64>,
    zscore: Vec<f64>,
    library: String,
    overlap_genes: Vec<String>,
    q_value: Vec<f64>,
    combined_score: Vec<f64>,
}

impl EnrichrResult {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        rank: Vec<i32>,
        term: Vec<String>,
        p_value: Vec<f64>,
        zscore: Vec<f64>,
        library: String,
        overlap_genes: Vec<String>,
        q_value: Vec<f64>,
        combined_score: Vec<f64>,
    ) -> Self {
        EnrichrResult {
            rank,
            term,
            p_value,
            zscore,
            library,
            overlap_genes,
            q_value,
            combined_score,
        }
    }

    pub fn empty(library: &str) -> Self {
        Self::new(
            vec![0],
            vec![String::new()],
            vec![1.0],
            vec![0.0],
            library.to_owned(),
            vec![String::new()],
            vec![1.0],
            vec![0.0],
        )
    }
    pub fn new_from_json(json: &Value, library: &str) -> Self {
        let mut ranks: Vec<i32> = Vec::new();
        let mut terms: Vec<String> = Vec::new();
        let mut p_values: Vec<f64> = Vec::new();
        let mut zscores: Vec<f64> = Vec::new();
        let mut combined_scores: Vec<f64> = Vec::new();
        let mut overlap_genes_str: Vec<String> = Vec::new();
        let mut q_values: Vec<f64> = Vec::new();

        if let Some(array) = json.as_array() {
            for item in array {
                if let Some(row) = item.as_array()
                    && row.len() >= 7
                {
                    ranks.push(row[0].as_i64().unwrap_or(0) as i32);
                    terms.push(row[1].as_str().unwrap_or_default().to_string());
                    p_values.push(row[2].as_f64().unwrap_or(1.0));
                    zscores.push(row[3].as_f64().unwrap_or(0.0));
                    combined_scores.push(row[4].as_f64().unwrap_or(0.0));

                    let genes_str = row[5]
                        .as_array()
                        .map(|a| a.iter().filter_map(|v| v.as_str()).collect_vec().join(", "))
                        .unwrap_or_default();
                    overlap_genes_str.push(genes_str);

                    q_values.push(row[6].as_f64().unwrap_or(1.0));
                }
            }
        }

        Self::new(
            ranks,
            terms,
            p_values,
            zscores,
            library.to_owned(),
            overlap_genes_str,
            q_values,
            combined_scores,
        )
    }

    pub fn to_dataframe(&self) -> PolarsResult<DataFrame> {
        df!(
        "rank" => self.rank.clone(),
        "term" => self.term.clone(),
        "p_value" => self.p_value.clone(),
        "zscore" => self.zscore.clone(),
        "combined_score" => self.combined_score.clone(),
        "overlap_genes" => self.overlap_genes.clone(),
        "q_value" => self.q_value.clone(),
        "library" => vec![self.library.clone(); self.rank.len()]
        )
    }

    pub fn get_all_rows_as_values(&self) -> Vec<Vec<Value>> {
        // Convert the internal fields back into row arrays for storage
        (0..self.rank.len())
            .map(|i| {
                vec![
                    Value::from(self.rank[i]),
                    Value::from(self.term[i].clone()),
                    Value::from(self.p_value[i]),
                    Value::from(self.zscore[i]),
                    Value::from(self.combined_score[i]),
                    Value::from(
                        self.overlap_genes[i]
                            .split(", ")
                            .map(|s| Value::from(s.to_string()))
                            .collect::<Vec<_>>(),
                    ),
                    Value::from(self.q_value[i]),
                ]
            })
            .collect()
    }
    pub fn get_top_n(&self, n: usize) -> Self {
        let EnrichrResult {
            rank,
            term,
            p_value,
            zscore,
            combined_score,
            overlap_genes,
            q_value,
            library,
        } = self;
        let rank = rank.clone().into_iter().take(n).collect();
        let term = term.clone().into_iter().take(n).collect();
        let p_value = p_value.clone().into_iter().take(n).collect();
        let zscore = zscore.clone().into_iter().take(n).collect();
        let combined_score = combined_score.clone().into_iter().take(n).collect();
        let overlap_genes = overlap_genes.clone().into_iter().take(n).collect();
        let q_value = q_value.clone().into_iter().take(n).collect();
        Self::new(
            rank,
            term,
            p_value,
            zscore,
            library.clone(),
            overlap_genes,
            q_value,
            combined_score,
        )
    }
}

#[derive(Debug, Clone)]
pub struct Enrichment {
    pub gene_list: Vec<String>,
    pub libraries: Vec<String>,
    pub background: Option<Vec<String>>,
    pub results: Vec<EnrichrResult>,
}

impl Enrichment {
    pub fn new(gene_list: Vec<String>, libraries: Vec<String>) -> Self {
        Self {
            gene_list,
            libraries,
            background: None,
            results: Vec::new(),
        }
    }

    pub fn with_background(&mut self, background: Vec<String>) -> &mut Self {
        self.background = Some(background);
        self
    }

    pub async fn run<A: EnrichrAPITrait>(
        &mut self,
        api: &mut A,
    ) -> Result<&mut Self, Box<dyn Error>> {
        api.send_genes(&self.gene_list, &self.libraries, false)
            .await?;
        if self.background.is_some() {
            api.send_genes(&self.gene_list, &self.libraries, true)
                .await?;
        }

        for lib in &self.libraries {
            self.results.push(api.enrich(lib).await?);
        }

        Ok(self)
    }
    pub fn get_short_id<A: EnrichrAPITrait>(&self, api: &A) -> Option<String> {
        api.get_short_id()
    }

    pub fn save_results(&self, path_buf: PathBuf) -> Result<(), Box<dyn Error>> {
        let df: Vec<LazyFrame> = self
            .results
            .iter()
            .map(|result| result.to_dataframe().map(|df| df.lazy()))
            .collect::<Result<Vec<_>, _>>()?;
        let combined_df = concat(df, UnionArgs::default())?;
        println!("{}", combined_df.clone().collect()?);
        tokio::task::block_in_place(|| combined_df.write_to_tsv_or_stdout(path_buf));
        Ok(())
    }

    fn parse_color(s: &str) -> RGBColor {
        // Named fallbacks
        match s {
            "lightskyblue" => RGBColor(135, 206, 250),
            "lightgrey" | "lightgray" => RGBColor(211, 211, 211),
            "black" => RGBColor(0, 0, 0),
            "blue" => RGBColor(0, 0, 255),
            "red" => RGBColor(255, 0, 0),
            _ => {
                // Accept hex like "#RRGGBB"
                if let Some(hex) = s.strip_prefix('#')
                    && hex.len() == 6
                    && let Ok(rgb) = u32::from_str_radix(hex, 16)
                {
                    let r = ((rgb >> 16) & 0xFF) as u8;
                    let g = ((rgb >> 8) & 0xFF) as u8;
                    let b = (rgb & 0xFF) as u8;
                    return RGBColor(r, g, b);
                }
                // default
                RGBColor(135, 206, 250)
            }
        }
    }
    fn draw_message_on_root<DB: DrawingBackend>(
        root: &DrawingArea<DB, Shift>,
        message: &str,
        width: u32,
        height: u32,
    ) -> Result<(), Box<dyn Error + Send + Sync>>
    where
        DB::ErrorType: 'static,
    {
        const CENTER_DIVISOR: u32 = 2;
        const TEXT_X_OFFSET: i32 = 10;
        root.fill(&WHITE)?;
        // Draw message near center (x offset a little to the left)

        let x = (width / CENTER_DIVISOR) as i32 - TEXT_X_OFFSET;
        let y = (height / CENTER_DIVISOR) as i32;
        root.draw(&Text::new(message, (x, y), ("sans-serif", 24).into_font()))?;
        root.present()?;
        Ok(())
    }
    fn render_empty_to_path(
        _library: &str,
        path: &Path,
        message: &str,
        width: u32,
        height: u32,
    ) -> Result<(), Box<dyn Error + Send + Sync>> {
        let ext = path
            .extension()
            .and_then(|s| s.to_str())
            .unwrap_or("")
            .to_lowercase();

        match ext.as_str() {
            "svg" => {
                let root = SVGBackend::new(path, (width, height)).into_drawing_area();
                Self::draw_message_on_root(&root, message, width, height)
            }
            "pdf" => {
                let surface = cairo::PdfSurface::new(f64::from(width), f64::from(height), path)?;
                let context = cairo::Context::new(&surface)?;
                let root = CairoBackend::new(&context, (width, height))?.into_drawing_area();
                Self::draw_message_on_root(&root, message, width, height)?;
                surface.flush();
                surface.finish();
                Ok(())
            }
            "png" | "jpg" | "jpeg" => {
                let path_str = path
                    .to_str()
                    .ok_or_else(|| format!("Invalid path: {}", path.display()))?;
                let root = BitMapBackend::new(path_str, (width, height)).into_drawing_area();
                Self::draw_message_on_root(&root, message, width, height)
            }
            other => Err(format!("Unsupported output extension: {other}").into()),
        }
    }
    #[allow(clippy::too_many_arguments)]
    fn render_bar_chart<DB: DrawingBackend>(
        root: &DrawingArea<DB, Shift>,
        library: &str,
        terms: &[String],
        p_vals: &[f64],
        neg_log_p: &[f64],
        max_x: f64,
        n: usize,
        color_primary: RGBColor,
        color_secondary: RGBColor,
    ) -> Result<(), Box<dyn Error + Send + Sync>>
    where
        DB::ErrorType: 'static,
    {
        root.fill(&WHITE)?;

        let mut chart = ChartBuilder::on(root)
            .caption(library.replace('_', " "), ("sans-serif", 24))
            .margin(10)
            .x_label_area_size(50)
            .y_label_area_size(220)
            .build_cartesian_2d(0f64..(max_x * 1.1), 0usize..n)?;

        chart
            .configure_mesh()
            .x_desc("-log10(p-value)")
            .disable_mesh()
            .y_labels(n)
            .y_label_formatter(&|idx| terms.get(*idx).cloned().unwrap_or_default())
            .draw()?;

        chart.draw_series((0..n).map(|i| {
            let x = neg_log_p.get(i).copied().unwrap_or(0.0);
            let pval = p_vals.get(i).copied().unwrap_or(1.0);
            let fill = if pval < 0.05 {
                color_primary.filled()
            } else {
                color_secondary.filled()
            };
            Rectangle::new([(0.0, i), (x, i + 1)], fill)
        }))?;

        root.present()?;
        Ok(())
    }
    fn draw_bar_plot_file(
        results: &EnrichrResult,
        library: &str,
        path: &Path,
        color_primary: RGBColor,
        color_secondary: RGBColor,
        top_n: usize,
    ) -> Result<(), Box<dyn Error + Send + Sync>> {
        const BAR_HEIGHT_PER_ITEM: u32 = 70;
        const MIN_BAR_ROWS: u32 = 5;

        // Build top-N
        let top = results.get_top_n(top_n);

        let terms: Vec<String> = top.term.clone();
        let p_vals: Vec<f64> = top.p_value.clone();
        let neg_log_p: Vec<f64> = p_vals
            .iter()
            .map(|p| if *p > 0.0 { -p.log10() } else { 0.0 })
            .collect();

        // use minimum length to avoid mismatches and panics
        let n = std::cmp::min(terms.len(), std::cmp::min(p_vals.len(), neg_log_p.len()));
        let width = 1000u32;
        let height = BAR_HEIGHT_PER_ITEM * u32::try_from(n)?.max(MIN_BAR_ROWS);

        if n == 0 {
            let message = format!("No rows to plot for {}", library);
            return Self::render_empty_to_path(library, path, &message, width, height);
        }

        let max_x = neg_log_p.iter().copied().fold(0.0_f64, f64::max).max(1.0);

        let ext = path
            .extension()
            .and_then(|s| s.to_str())
            .unwrap_or("")
            .to_lowercase();

        match ext.as_str() {
            "svg" => {
                let root = SVGBackend::new(path, (width, height)).into_drawing_area();
                Self::render_bar_chart(
                    &root,
                    library,
                    &terms,
                    &p_vals,
                    &neg_log_p,
                    max_x,
                    n,
                    color_primary,
                    color_secondary,
                )
            }
            "pdf" => {
                let surface = cairo::PdfSurface::new(f64::from(width), f64::from(height), path)?;
                let context = cairo::Context::new(&surface)?;
                let root = CairoBackend::new(&context, (width, height))?.into_drawing_area();

                Self::render_bar_chart(
                    &root,
                    library,
                    &terms,
                    &p_vals,
                    &neg_log_p,
                    max_x,
                    n,
                    color_primary,
                    color_secondary,
                )?;
                surface.flush();
                surface.finish();
                Ok(())
            }
            "png" | "jpg" | "jpeg" => {
                let path_str = path
                    .to_str()
                    .ok_or_else(|| format!("Invalid path: {}", path.display()))?;
                let root = BitMapBackend::new(path_str, (width, height)).into_drawing_area();
                Self::render_bar_chart(
                    &root,
                    library,
                    &terms,
                    &p_vals,
                    &neg_log_p,
                    max_x,
                    n,
                    color_primary,
                    color_secondary,
                )
            }
            other => Err(format!("Unsupported output extension: {other}").into()),
        }
    }

    pub async fn bar_plot(
        &self,
        paths: Vec<PathBuf>,
        library: Option<String>,
        color: Option<String>,
        color2: Option<String>,
    ) -> Result<(), Box<dyn Error>> {
        let lib = library.unwrap_or_else(|| self.libraries[0].clone());
        let color = color.unwrap_or_else(|| "lightskyblue".to_string());
        let color2 = color2.unwrap_or_else(|| "lightgrey".to_string());

        let result = self
            .results
            .iter()
            .find(|r| r.library == lib)
            .ok_or_else(|| {
                Box::new(std::io::Error::other(format!(
                    "No results found for library: {lib}"
                ))) as Box<dyn Error>
            })?;

        // Clone for blocking thread
        let result_clone = result.clone();
        let lib_clone = lib.clone();
        let color_primary = Enrichment::parse_color(&color);
        let color_secondary = Enrichment::parse_color(&color2);
        let top_n = 10usize;
        let paths_clone = paths.clone();

        // spawn_blocking returns JoinHandle<Result<(), Box<dyn Error + Send + Sync>>>
        let join_result = tokio::task::spawn_blocking(move || {
            // closure runs on blocking thread
            for path in &paths_clone {
                Enrichment::draw_bar_plot_file(
                    &result_clone,
                    &lib_clone,
                    path.as_path(),
                    color_primary,
                    color_secondary,
                    top_n,
                )?;
            }
            Ok::<(), Box<dyn Error + Send + Sync>>(())
        })
        .await;

        // Map JoinError -> plain boxed error, then handle inner boxed error explicitly
        let inner_res: Result<(), Box<dyn Error + Send + Sync>> = match join_result {
            Ok(res) => res,
            Err(join_err) => {
                return Err(Box::from(format!(
                    "Failed to execute bar plot rendering task: {join_err}"
                )));
            }
        };

        if let Err(inner_err) = inner_res {
            // convert inner boxed Send+Sync error into a plain boxed error for this API
            return Err(Box::from(inner_err.to_string()));
        }

        Ok(())
    }
}

pub async fn enrich_command(
    library: String,
    gene_list: PathBuf,
    background: Option<PathBuf>,
    output_file: PathBuf,
    output_plot: Vec<PathBuf>,
) -> Result<(), Box<dyn Error>> {
    let genes = get_genes(&gene_list)?;
    let libraries = vec![library.clone()];
    let mut enrich = Enrichment::new(genes, libraries);

    if let Some(path) = &background {
        let bg_genes: Vec<String> = get_genes(path)?;
        enrich.with_background(bg_genes);
    }
    let mut api = EnrichrAPI::new(enrich.background.clone());
    enrich.run(&mut api).await?;

    if let Some(short_id) = enrich.get_short_id(&api)
        && background.is_none()
    {
        println!(
            "Results can be found at: https://maayanlab.cloud/Enrichr/enrich?dataset={short_id}"
        );
    }
    enrich.save_results(output_file)?;
    enrich
        .bar_plot(output_plot, Some(library), None, None)
        .await?;
    Ok(())
}

fn get_genes(gene_list: &PathBuf) -> Result<Vec<String>, Box<dyn Error>> {
    let genes: Vec<String> = fs::read_to_string(gene_list)?
        .lines()
        .map(|line| line.trim().to_string())
        .filter(|line| !line.is_empty())
        .collect();
    Ok(genes)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::enrich::api::APIFailure;
    use crate::test_helpers::*;
    use pretty_assertions::{assert_eq, assert_ne};
    use rstest::{Context, fixture, rstest};
    use std::time::Duration;
    use temp_testdir::TempDir;

    #[fixture]
    fn create_test_results() -> EnrichrResult {
        EnrichrResult::new(
            vec![1, 2, 3, 4, 5],
            vec![
                "t1".to_string(),
                "t2".to_string(),
                "t3".to_string(),
                "t4".to_string(),
                "t5".to_string(),
            ],
            vec![0.01, 0.02, 0.03, 0.04, 0.05],
            vec![5.0, 4.0, 3.0, 2.0, 1.0],
            "lib1".to_string(),
            vec![
                "g1".to_string(),
                "g2".to_string(),
                "g3".to_string(),
                "g4".to_string(),
                "g5".to_string(),
            ],
            vec![0.01, 0.02, 0.03, 0.04, 0.05],
            vec![5.0, 4.0, 3.0, 2.0, 1.0],
        )
    }

    #[rstest]
    fn test_get_genes() {
        let gene_list = PathBuf::from("tests/data/example_gene_list.txt");
        let actual = get_genes(&gene_list).unwrap();
        assert_eq!(actual.len(), 614);
        assert_eq!(actual[0], "A1CF");
        assert_eq!(actual[613], "ZNF831");
    }

    #[rstest]
    #[timeout(Duration::from_secs(10))]
    #[tokio::test(flavor = "multi_thread", worker_threads = 1)]
    async fn test_enrich_command_to_api_backend(#[values(true, false)] use_background: bool) {
        let library = "Reactome_Pathways_2024".to_string();
        let gene_list = PathBuf::from("tests/data/example_gene_list.txt");

        let genes = get_genes(&gene_list).unwrap();
        let libraries = vec![library.clone()];
        let mut enrich = Enrichment::new(genes, libraries);

        if use_background {
            let background = PathBuf::from("tests/data/example_background.txt");
            let bg_genes: Vec<String> = get_genes(&background).unwrap();
            enrich.with_background(bg_genes);
        }

        let mut api = EnrichrAPI::new(enrich.background.clone());
        enrich.run(&mut api).await.unwrap();
        let temp = TempDir::default();
        let tsv_file = temp.with_file_name(get_timestamp()).with_extension("tsv");
        let svg_file = temp.with_file_name(get_timestamp()).with_extension("svg");
        let pdf_file = temp.with_file_name(get_timestamp()).with_extension("pdf");
        let png_file = temp.with_file_name(get_timestamp()).with_extension("png");
        assert!(enrich.save_results(tsv_file).is_ok());
        assert!(
            enrich
                .bar_plot(
                    vec![svg_file, pdf_file, png_file],
                    Some(library),
                    None,
                    None
                )
                .await
                .is_ok()
        );

        let df = enrich.results[0].to_dataframe().unwrap();
        assert_eq!(df.height(), 791);
        assert_eq!(df.width(), 8);
    }

    #[test]
    fn test_enrichr_result_from_json() {
        let json = serde_json::json!([
            [1, "term1", 0.01, 2.0, 3.0, ["gene1", "gene2"], 0.05],
            [2, "term2", 0.05, 1.0, 2.0, ["gene3"], 0.1]
        ]);
        let result = EnrichrResult::new_from_json(&json, "lib1");
        assert_eq!(result.rank, vec![1, 2]);
        assert_eq!(result.term, vec!["term1", "term2"]);
        assert_eq!(result.p_value, vec![0.01, 0.05]);
        assert_eq!(result.zscore, vec![2.0, 1.0]);
        assert_eq!(result.combined_score, vec![3.0, 2.0]);
        assert_eq!(result.overlap_genes, vec!["gene1, gene2", "gene3"]);
        assert_eq!(result.q_value, vec![0.05, 0.1]);
        assert_eq!(result.library, "lib1");
    }

    struct MockEnrichrAPI;
    impl EnrichrAPITrait for MockEnrichrAPI {
        async fn send_genes(
            &mut self,
            _gene_list: &[String],
            _libraries: &[String],
            _send_background: bool,
        ) -> Result<(), APIFailure> {
            Ok(())
        }
        async fn enrich(&mut self, library_name: &String) -> Result<EnrichrResult, APIFailure> {
            Ok(EnrichrResult::empty(library_name))
        }
        fn get_short_id(&self) -> Option<String> {
            Some("short_id".to_string())
        }
    }

    #[rstest]
    #[timeout(Duration::from_secs(10))]
    #[tokio::test(flavor = "multi_thread", worker_threads = 1)]
    async fn test_enrichment_run() {
        let mut enrich = Enrichment::new(vec!["gene1".to_string()], vec!["lib1".to_string()]);
        let mut api = MockEnrichrAPI;
        enrich.run(&mut api).await.unwrap();
        assert_eq!(enrich.results.len(), 1);
        assert_eq!(enrich.results[0].library, "lib1");

        let temp = TempDir::default();
        let tsv_file = temp.with_file_name(get_timestamp()).with_extension("tsv");
        let svg_file = temp.with_file_name(get_timestamp()).with_extension("svg");
        let pdf_file = temp.with_file_name(get_timestamp()).with_extension("pdf");
        let png_file = temp.with_file_name(get_timestamp()).with_extension("png");
        assert!(enrich.save_results(tsv_file).is_ok());
        assert!(
            enrich
                .bar_plot(
                    vec![svg_file, pdf_file, png_file],
                    Some("lib1".to_string()),
                    None,
                    None
                )
                .await
                .is_ok()
        );
    }

    #[rstest]
    fn test_enrichr_result_to_dataframe(#[from(create_test_results)] result: EnrichrResult) {
        let df = result.to_dataframe().unwrap();
        assert_eq!(df.height(), 5);
    }

    #[rstest]
    fn test_enrichr_result_get_all_rows_as_values(
        #[from(create_test_results)] result: EnrichrResult,
    ) {
        let rows = result.get_all_rows_as_values();
        assert_eq!(rows.len(), 5);
        assert_eq!(rows[0].len(), 7);
    }

    #[rstest]
    #[case(2)]
    #[case(3)]
    #[case(10)]
    fn test_enrichr_result_get_top_n(
        #[context] ctx: Context,
        #[from(create_test_results)] result: EnrichrResult,
        #[case] n: usize,
    ) {
        let top = result.get_top_n(n);
        let n = if ctx.case == Some(2) {
            // Test case 2: n is larger than the number of rows
            // Want to be sure get_top_n() does not add extra rows
            // Then set n to 5 to make sure the returned top-N the number of rows in result
            assert_ne!(top.rank.len(), n);
            5
        } else {
            n
        };
        assert_eq!(top.rank.len(), n);
        assert_eq!(top.rank, (1..=n as i32).collect::<Vec<i32>>());
    }

    #[test]
    fn test_enrichment_with_background() {
        let mut enrich = Enrichment::new(vec!["gene1".to_string()], vec!["lib1".to_string()]);
        assert!(enrich.background.is_none());
        enrich.with_background(vec!["bg1".to_string()]);
        assert!(enrich.background.is_some());
    }

    #[rstest]
    #[case::red("red", RGBColor(255, 0, 0))]
    #[case::blue("blue", RGBColor(0, 0, 255))]
    #[case::white_hex("#FFFFFF", RGBColor(255, 255, 255))]
    #[case::black_hex("#000000", RGBColor(0, 0, 0))]
    #[case::default_hex("#GGGGGG", RGBColor(135, 206, 250))]
    #[case::invalid("#FF", RGBColor(135, 206, 250))]
    fn test_parse_color(#[case] input: &str, #[case] expected: RGBColor) {
        assert_eq!(Enrichment::parse_color(input), expected);
    }

    #[rstest]
    #[timeout(Duration::from_secs(10))]
    #[tokio::test(flavor = "multi_thread", worker_threads = 1)]
    async fn test_enrichr_result_empty() {
        let libraries = vec!["test_lib".to_string()];
        let gene_list = vec!["gene1".to_string()];
        let result = EnrichrResult::empty("test_lib");
        assert_eq!(result.library, "test_lib");
        assert_eq!(result.rank.len(), 1);
        assert_eq!(result.p_value, vec![1.0]);

        let enrich = Enrichment {
            gene_list,
            libraries,
            background: None,
            results: vec![result],
        };

        let temp = TempDir::default();
        let tsv_file = temp.with_file_name(get_timestamp()).with_extension("tsv");
        let svg_file = temp.with_file_name(get_timestamp()).with_extension("svg");
        let pdf_file = temp.with_file_name(get_timestamp()).with_extension("pdf");
        let png_file = temp.with_file_name(get_timestamp()).with_extension("png");
        assert!(enrich.save_results(tsv_file).is_ok());
        assert!(
            enrich
                .bar_plot(
                    vec![svg_file, pdf_file, png_file],
                    Some("test_lib".to_string()),
                    None,
                    None
                )
                .await
                .is_ok()
        );
    }

    #[test]
    fn test_enrichr_result_new_from_json_empty() {
        let json = serde_json::json!([]);
        let result = EnrichrResult::new_from_json(&json, "lib1");
        assert_eq!(result.rank.len(), 0);
        assert_eq!(result.library, "lib1");
    }

    #[test]
    fn test_enrichr_result_new_from_json_incomplete_row() {
        let json = serde_json::json!([[1, "term1"], [2, "term2", 0.05, 1.0, 2.0, ["gene3"], 0.1]]);
        let result = EnrichrResult::new_from_json(&json, "lib1");
        assert_eq!(result.rank.len(), 1);
        assert_eq!(result.rank[0], 2);
    }

    #[test]
    fn test_enrichment_new() {
        let enrich = Enrichment::new(vec!["gene1".to_string()], vec!["lib1".to_string()]);
        assert_eq!(enrich.gene_list.len(), 1);
        assert_eq!(enrich.libraries.len(), 1);
        assert!(enrich.background.is_none());
        assert_eq!(enrich.results.len(), 0);
    }
}
