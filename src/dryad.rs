use anyhow::{Context, anyhow};
use clap::Subcommand;
use reqwest::blocking::{Body, Client, Response};
use reqwest::header::{ACCEPT, AUTHORIZATION, CONTENT_TYPE};
use serde::Deserialize;
use std::fs::File;
use std::path::{Path, PathBuf};
use urlencoding::encode as url_encode;

pub type Doi = String;
pub type FilePath = PathBuf;

pub struct DryadApiConfig;

impl DryadApiConfig {
    pub const API_URL: &'static str = "https://datadryad.org/api/v2/datasets";
    pub const TOKEN_URL: &'static str = "https://datadryad.org/oauth/token";
}

#[derive(Subcommand)]
pub enum Commands {
    /// Upload files in a directory to a Dryad dataset DOI
    DryadUpload {
        #[arg(long, env = "DRYAD_CLIENT_ID")]
        client_id: String,
        #[arg(long = "client-secret", env = "DRYAD_SECRET")]
        client_secret: String,
        #[arg(long)]
        doi: Doi,
        #[arg(default_value = ".")]
        directory: FilePath,
    },
}

pub fn handle_command(cmd: Commands) -> Result<(), clap::Error> {
    match cmd {
        Commands::DryadUpload {
            client_id,
            client_secret,
            doi,
            directory,
        } => upload_to_dryad(&client_id, &client_secret, &directory, &doi)
            .map_err(|e| clap::Error::raw(clap::error::ErrorKind::Io, e.to_string())),
    }
}

#[derive(Debug)]
pub struct DryadClient {
    client_id: String,
    client_secret: String,
    http: Client,
    token: String,
}

#[derive(Debug, Deserialize)]
struct TokenResponse {
    access_token: String,
}

impl DryadClient {
    pub fn new(
        client_id: impl Into<String>,
        client_secret: impl Into<String>,
    ) -> anyhow::Result<Self> {
        let mut client = Self {
            client_id: client_id.into(),
            client_secret: client_secret.into(),
            http: Client::new(),
            token: String::new(),
        };
        client.refresh_token()?;
        Ok(client)
    }

    fn refresh_token(&mut self) -> anyhow::Result<()> {
        self.token = self.get_token()?;
        Ok(())
    }

    fn get_token(&self) -> anyhow::Result<String> {
        let response = self
            .http
            .post(DryadApiConfig::TOKEN_URL)
            .form(&[
                ("client_id", self.client_id.as_str()),
                ("client_secret", self.client_secret.as_str()),
                ("grant_type", "client_credentials"),
            ])
            .send()
            .context("failed to request Dryad token")?
            .error_for_status()
            .context("Dryad token endpoint returned error")?;

        let token: TokenResponse = response
            .json()
            .context("failed to parse Dryad token response")?;
        Ok(token.access_token)
    }

    fn authorized_put(
        &self,
        url: &str,
        content_type: &str,
        file: File,
    ) -> anyhow::Result<Response> {
        self.http
            .put(url)
            .header(AUTHORIZATION, format!("Bearer {}", self.token))
            .header(ACCEPT, "application/json")
            .header(CONTENT_TYPE, content_type)
            .body(Body::new(file))
            .send()
            .context("failed to upload file")
    }

    fn upload_single_file(&mut self, file: &Path, doi_encoded: &str) -> anyhow::Result<()> {
        let file_name = file
            .file_name()
            .and_then(|s| s.to_str())
            .ok_or_else(|| anyhow!("invalid file name: {}", file.display()))?;
        let file_name_encoded = Self::encode(file_name);
        let file_url = format!(
            "{}/{}/files/{}",
            DryadApiConfig::API_URL,
            doi_encoded,
            file_name_encoded
        );

        let mime = detect_mime_type(file)?;
        let handle =
            File::open(file).with_context(|| format!("failed to open {}", file.display()))?;
        let response = self.authorized_put(&file_url, mime, handle)?;

        if response.status().is_client_error() {
            self.refresh_token()?;
            let retry_handle =
                File::open(file).with_context(|| format!("failed to reopen {}", file.display()))?;
            self.authorized_put(&file_url, mime, retry_handle)?
                .error_for_status()
                .with_context(|| format!("failed upload for {}", file.display()))?;
        } else {
            response
                .error_for_status()
                .with_context(|| format!("failed upload for {}", file.display()))?;
        }
        Ok(())
    }

    pub fn upload_files(&mut self, files: &[FilePath], doi: &str) -> anyhow::Result<()> {
        let doi_encoded = Self::encode(doi);
        for file in files {
            self.upload_single_file(file, &doi_encoded)?;
        }
        Ok(())
    }

    pub fn encode(s: &str) -> String {
        url_encode(s).to_string()
    }
}

fn detect_mime_type(path: &Path) -> anyhow::Result<&'static str> {
    let file_name = path
        .file_name()
        .and_then(|s| s.to_str())
        .ok_or_else(|| anyhow!("invalid file name: {}", path.display()))?;
    if file_name.ends_with(".tar.gz") {
        return Ok("application/tar+gzip");
    }
    match path.extension().and_then(|s| s.to_str()) {
        Some("h5") => Ok("application/x-hdf"),
        Some("csv") => Ok("text/csv"),
        Some("jpg") | Some("jpeg") => Ok("image/jpeg"),
        Some("png") => Ok("image/png"),
        Some("pdf") => Ok("application/pdf"),
        Some("txt") => Ok("text/plain"),
        Some("tiff") | Some("tif") => Ok("image/tiff"),
        Some("svg") => Ok("image/svg+xml"),
        _ => Err(anyhow!("unsupported file type: {}", path.display())),
    }
}

pub fn upload_to_dryad(
    client_id: &str,
    client_secret: &str,
    directory: &Path,
    doi: &str,
) -> anyhow::Result<()> {
    let files: Vec<PathBuf> = std::fs::read_dir(directory)
        .with_context(|| format!("failed to read directory {}", directory.display()))?
        .filter_map(|entry| entry.ok().map(|e| e.path()))
        .filter(|path| path.is_file())
        .collect();

    let mut client = DryadClient::new(client_id, client_secret)?;
    client.upload_files(&files, doi)
}
