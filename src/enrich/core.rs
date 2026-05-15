use crate::enrich::api::EnrichrAPI;
use crate::io::WriteToCsvOrStdout;
use cairo;
use plotters::coord::Shift;
use plotters::prelude::*;
use plotters_cairo::CairoBackend;
use polars::polars_utils::itertools::Itertools;
use polars::prelude::*;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::path::{Path, PathBuf};

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
    api: Option<EnrichrAPI>,
    results: Vec<EnrichrResult>,
}

impl Enrichment {
    pub fn new(gene_list: Vec<String>, libraries: Vec<String>) -> Self {
        Self {
            gene_list,
            libraries,
            background: None,
            api: None,
            results: Vec::new(),
        }
    }

    pub fn with_background(&mut self, background: Vec<String>) -> &mut Self {
        self.background = Some(background);
        self
    }

    pub fn build(&mut self) -> &mut Self {
        self.api = Some(EnrichrAPI::new(self.to_owned()));
        self
    }

    pub async fn run(&mut self) -> Result<&mut Self, Box<dyn std::error::Error>> {
        match &mut self.api {
            Some(api) => {
                api.send_genes(&self.gene_list, &self.libraries, false)
                    .await?;
                if self.background.is_some() {
                    api.send_genes(&self.gene_list, &self.libraries, true)
                        .await?;
                }

                for lib in &self.libraries {
                    self.results.push(api.enrich(lib).await?);
                }
            }
            None => {
                return Err(Box::new(std::io::Error::other(
                    "API not initialized. Call build() before run().",
                )));
            }
        }

        Ok(self)
    }
    pub fn get_short_id(&self) -> Option<String> {
        match self.api {
            None => None,
            Some(ref api) => api.get_short_id(),
        }
    }

    pub fn save_results(&self, path_buf: PathBuf) -> Result<(), Box<dyn std::error::Error>> {
        let df: Vec<LazyFrame> = self
            .results
            .iter()
            .map(|result| result.to_dataframe().unwrap().lazy())
            .collect();
        let combined_df = concat(df, UnionArgs::default())?;
        println!("{}", combined_df.clone().collect()?);
        tokio::task::block_in_place(|| {
            combined_df.write_to_tsv_or_stdout(path_buf);
        });
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
    ) -> Result<(), Box<dyn std::error::Error + Send + Sync>>
    where
        DB::ErrorType: 'static,
    {
        root.fill(&WHITE)?;
        // Draw message near center (x offset a little to the left)
        let x = (width / 2) as i32 - 10;
        let y = (height / 2) as i32;
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
    ) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
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
    ) -> Result<(), Box<dyn std::error::Error + Send + Sync>>
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
    ) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
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
        let height = 70u32 * u32::try_from(n)?.max(5);

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
    ) -> Result<(), Box<dyn std::error::Error>> {
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
                ))) as Box<dyn std::error::Error>
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
            Ok::<(), Box<dyn std::error::Error + Send + Sync>>(())
        })
        .await;

        // Map JoinError -> plain boxed error, then handle inner boxed error explicitly
        let inner_res: Result<(), Box<dyn std::error::Error + Send + Sync>> = match join_result {
            Ok(res) => res,
            Err(join_err) => {
                return Err(Box::from(format!("JoinError: {join_err}")));
            }
        };

        if let Err(inner_err) = inner_res {
            // convert inner boxed Send+Sync error into a plain boxed error for this API
            return Err(Box::from(inner_err.to_string()));
        }

        Ok(())
    }
}
