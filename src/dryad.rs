use crate::progress::{Progress, ProgressArg, byte_and_file_progress_bar};
use anyhow::{Context, Result, anyhow};
use clap::Subcommand;
use reqwest::Url;
use reqwest::blocking::{Body, Client, Response};
use reqwest::header::{ACCEPT, AUTHORIZATION, CONTENT_TYPE};
use serde::Deserialize;
use std::fs::File;
use std::path::{Path, PathBuf};

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
        #[command(flatten)]
        progress: ProgressArg,
    },
}

pub fn handle_command(cmd: Commands) -> Result<()> {
    match cmd {
        Commands::DryadUpload {
            client_id,
            client_secret,
            doi,
            directory,
            progress,
        } => upload_to_dryad(
            &client_id,
            &client_secret,
            &directory,
            &doi,
            progress.get().unwrap_or(Progress::Progress),
        ),
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
    pub fn new(client_id: impl Into<String>, client_secret: impl Into<String>) -> Result<Self> {
        let mut client = Self {
            client_id: client_id.into(),
            client_secret: client_secret.into(),
            http: Client::new(),
            token: String::new(),
        };
        client.refresh_token()?;
        Ok(client)
    }

    fn refresh_token(&mut self) -> Result<()> {
        self.token = self.get_token()?;
        Ok(())
    }

    fn get_token(&self) -> Result<String> {
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

    fn authorized_put(&self, url: &str, content_type: &str, file: File) -> Result<Response> {
        self.http
            .put(url)
            .header(AUTHORIZATION, format!("Bearer {}", self.token))
            .header(ACCEPT, "application/json")
            .header(CONTENT_TYPE, content_type)
            .body(Body::new(file))
            .send()
            .context("failed to upload file")
    }

    fn ensure_success(response: Response, file: &Path) -> Result<()> {
        let status = response.status();
        if status.is_success() {
            return Ok(());
        }

        let body = response
            .text()
            .context("failed to read Dryad error response")?;
        let body = body.trim();

        if body.is_empty() {
            Err(anyhow!(
                "failed upload for {} (status {})",
                file.display(),
                status
            ))
        } else {
            Err(anyhow!(
                "failed upload for {} (status {}): {}",
                file.display(),
                status,
                body
            ))
        }
    }

    fn upload_single_file(&mut self, file: &Path, doi_encoded: &str) -> Result<()> {
        let file_name = file
            .file_name()
            .map(|s| s.to_string_lossy().into_owned())
            .ok_or_else(|| anyhow!("invalid file name: {}", file.display()))?;
        let file_url = Self::file_upload_url(doi_encoded, &file_name)?;

        let mime = detect_mime_type(file)?;
        let handle =
            File::open(file).with_context(|| format!("failed to open {}", file.display()))?;
        let response = self.authorized_put(&file_url, mime, handle)?;

        if response.status().as_u16() == 401 {
            self.refresh_token()?;
            let retry_handle =
                File::open(file).with_context(|| format!("failed to reopen {}", file.display()))?;
            let retry_response = self.authorized_put(&file_url, mime, retry_handle)?;
            Self::ensure_success(retry_response, file)?;
        } else {
            Self::ensure_success(response, file)?;
        }
        Ok(())
    }

    pub fn upload_files(
        &mut self,
        files: &[FilePath],
        doi: &str,
        progress: Progress,
    ) -> Result<()> {
        reject_files_with_spaces(files)?;
        let doi_encoded = Self::encode(doi);
        let total_bytes = files
            .iter()
            .map(|file| {
                std::fs::metadata(file)
                    .with_context(|| format!("failed to read metadata for {}", file.display()))
                    .map(|metadata| metadata.len())
            })
            .collect::<Result<Vec<_>>>()?
            .into_iter()
            .sum();
        let pb = byte_and_file_progress_bar(total_bytes, files.len(), progress)?;
        for file in files {
            self.upload_single_file(file, &doi_encoded)?;
            let len = std::fs::metadata(file)
                .with_context(|| format!("failed to read metadata for {}", file.display()))?
                .len();
            pb.inc_bytes(len);
            pb.inc_file();
        }
        Ok(())
    }

    pub fn encode(s: &str) -> String {
        urlencoding::encode(s).to_string()
    }

    fn file_upload_url(doi_encoded: &str, file_name: &str) -> Result<String> {
        let mut url = Url::parse(&format!(
            "{}/{}/files",
            DryadApiConfig::API_URL,
            doi_encoded
        ))
        .context("failed to construct Dryad upload URL")?;
        url.path_segments_mut()
            .map_err(|_| anyhow!("failed to construct Dryad upload URL"))?
            .push(file_name);
        Ok(url.to_string())
    }
}

fn reject_files_with_spaces(files: &[FilePath]) -> Result<()> {
    let invalid_files = files
        .iter()
        .filter(|path| {
            path.file_name()
                .map(|name| name.to_string_lossy().contains(' '))
                .unwrap_or(false)
        })
        .map(|path| path.display().to_string())
        .collect::<Vec<_>>();

    if invalid_files.is_empty() {
        return Ok(());
    }

    Err(anyhow!(
        "Dryad upload does not support filenames with spaces. Rename these files before uploading:\n{}",
        invalid_files.join("\n")
    ))
}

fn detect_mime_type(path: &Path) -> Result<&'static str> {
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
    progress: Progress,
) -> Result<()> {
    let files: Vec<PathBuf> = std::fs::read_dir(directory)
        .with_context(|| format!("failed to read directory {}", directory.display()))?
        .filter_map(|entry| entry.ok().map(|e| e.path()))
        .filter(|path| path.is_file())
        .collect();

    let mut client = DryadClient::new(client_id, client_secret)?;
    client.upload_files(&files, doi, progress)
}

#[cfg(test)]
mod tests {
    use super::{DryadClient, reject_files_with_spaces};
    use std::path::PathBuf;

    #[test]
    fn encodes_spaces_in_file_names() {
        assert_eq!(DryadClient::encode("my file.csv"), "my%20file.csv");
        assert_eq!(
            DryadClient::encode("CL022 Infected.tar.gz"),
            "CL022%20Infected.tar.gz"
        );
    }

    #[test]
    fn builds_upload_url_with_path_encoding_for_spaces() {
        let url =
            DryadClient::file_upload_url("doi%3A10.5061%2Fdryad.12345", "CL022 Infected.tar.gz")
                .unwrap();
        assert_eq!(
            url,
            "https://datadryad.org/api/v2/datasets/doi%3A10.5061%2Fdryad.12345/files/CL022%20Infected.tar.gz"
        );
    }

    #[test]
    fn rejects_filenames_with_spaces_before_upload() {
        let err = reject_files_with_spaces(&[
            PathBuf::from("/tmp/no_spaces.csv"),
            PathBuf::from("/tmp/CL022 Infected.tar.gz"),
        ])
        .unwrap_err();

        assert_eq!(
            err.to_string(),
            "Dryad upload does not support filenames with spaces. Rename these files before uploading:\n/tmp/CL022 Infected.tar.gz"
        );
    }
}
