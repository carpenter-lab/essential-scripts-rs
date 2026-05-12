use crate::enrich::core::{Enrichment, EnrichrResult};
use reqwest;
use reqwest::multipart;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::HashMap;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ErichrResponse {
    #[serde(rename = "userListId")]
    pub user_list_id: i32,
    #[serde(rename = "shortId")]
    pub short_id: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BackgroundResponse {
    pub backgroundid: String,
}

#[derive(Debug)]
pub struct APIFailure {
    pub message: String,
}

impl std::fmt::Display for APIFailure {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "API Failure: {}", self.message)
    }
}

impl std::error::Error for APIFailure {}

impl APIFailure {
    pub fn new(message: String) -> Self {
        APIFailure { message }
    }
}

#[derive(Clone, Debug)]
pub struct EnrichrAPI {
    client: reqwest::Client,
    list_response: Option<ErichrResponse>,
    background_response: Option<BackgroundResponse>,
    enrich_response: HashMap<String, Option<Vec<Vec<Value>>>>,
    failure: HashMap<String, bool>,
    background: bool,
}

impl EnrichrAPI {
    const BASE_URL: &'static str = "https://maayanlab.cloud/Enrichr";
    const LIST_URL: &'static str = "addList";
    const ENRICH_URL: &'static str = "enrich";
    const SPEEDRICH_URL: &'static str = "speedrich/api";
    const BACKGROUND_URL: &'static str = "addbackground";
    const BACKGROUND_ENRICH_URL: &'static str = "backgroundenrich";

    pub fn new(enrichment: Enrichment) -> Self {
        match enrichment.background {
            Some(_background) => EnrichrAPI {
                client: reqwest::Client::new(),
                list_response: None,
                background_response: None,
                enrich_response: HashMap::new(),
                failure: HashMap::new(),
                background: true,
            },
            None => EnrichrAPI {
                client: reqwest::Client::new(),
                list_response: None,
                background_response: None,
                enrich_response: HashMap::new(),
                failure: HashMap::new(),
                background: false,
            },
        }
    }

    fn get_list_url(&self, background: bool) -> String {
        let url = match self.background {
            true => format!("{}/{}", Self::BASE_URL, Self::SPEEDRICH_URL),
            false => Self::BASE_URL.into(),
        };
        match background {
            true => format!("{}/{}", url, Self::BACKGROUND_URL),
            false => format!("{}/{}", url, Self::LIST_URL),
        }
    }

    pub async fn send_genes(
        &mut self,
        gene_list: &[String],
        libraries: &[String],
        send_background: bool,
    ) -> Result<(), APIFailure> {
        let genes_text = gene_list.join("\n");

        let form = match send_background {
            true => multipart::Form::new().text("background", genes_text),
            false => multipart::Form::new()
                .text("list", genes_text)
                .text("description", ""),
        };

        let response = self
            .client
            .post(self.get_list_url(send_background))
            .multipart(form)
            .send()
            .await
            .map_err(|e| APIFailure {
                message: format!("Failed to send genes: {}", e),
            })?;
        if !response.status().is_success() {
            for lib in libraries {
                self.failure.insert(lib.clone(), true);
            }
            return Err(APIFailure::new(format!(
                "Gene submission failed {}",
                response.text().await.unwrap_or_default()
            )));
        }
        let response_text = response.text().await.map_err(|e| APIFailure {
            message: format!("Failed to read response: {}", e),
        })?;
        match send_background {
            true => {
                self.background_response = Some(serde_json::from_str(&response_text).map_err(
                    |e| APIFailure {
                        message: format!("Failed to parse background response: {}", e),
                    },
                )?);
            }
            false => {
                self.list_response =
                    Some(
                        serde_json::from_str(&response_text).map_err(|e| APIFailure {
                            message: format!("Failed to parse list response: {}", e),
                        })?,
                    );
            }
        }
        Ok(())
    }

    pub async fn enrich(&mut self, library_name: &String) -> Result<EnrichrResult, APIFailure> {
        match self.failure.get(library_name) {
            Some(true) => {
                return Ok(EnrichrResult::empty(library_name));
            }
            Some(false) | None => {}
        }
        let user_list_id = &self
            .list_response
            .as_ref()
            .ok_or(APIFailure::new("No list response".to_string()))?
            .user_list_id;
        let url = match self.background_response.clone() {
            Some(background_response) => {
                format!(
                    "{}/{}/{}?userListId={}&backgroundId={}&backgroundType={}",
                    Self::BASE_URL,
                    Self::SPEEDRICH_URL,
                    Self::BACKGROUND_ENRICH_URL,
                    user_list_id,
                    background_response.backgroundid,
                    library_name
                )
            }
            None => {
                format!(
                    "{}/{}?userListId={}&backgroundType={}",
                    Self::BASE_URL,
                    Self::ENRICH_URL,
                    user_list_id,
                    library_name
                )
            }
        };

        let response = self.client.get(&url).send().await.map_err(|e| APIFailure {
            message: format!("Failed enrichment request: {}", e),
        })?;

        if !response.status().is_success() {
            self.failure.insert(library_name.clone(), true);
            return Err(APIFailure::new(format!(
                "Enrichment request failed for {}",
                library_name
            )));
        }

        let response_json: Value = response.json().await.map_err(|e| APIFailure {
            message: format!("Failed to parse enrichment response: {}", e),
        })?;

        if let Some(lib_data) = response_json.get(library_name).and_then(|v| v.as_array()) {
            if lib_data.is_empty() {
                self.failure.insert(library_name.clone(), true);
                return Ok(EnrichrResult::empty(library_name));
            } else {
                self.enrich_response.insert(
                    library_name.clone(),
                    Some(
                        lib_data
                            .iter()
                            .map(|v| v.as_array().map(|a| a.clone()).unwrap_or_default())
                            .collect(),
                    ),
                );
            }
        }

        if let Some(results) = response_json.get(library_name) {
            let r = EnrichrResult::new_from_json(results, library_name);
            return Ok(r);
        }

        Ok(EnrichrResult::empty(library_name))
    }
}
