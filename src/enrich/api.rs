use crate::enrich::core::EnrichrResult;
use reqwest;
use reqwest::multipart;
use serde::{Deserialize, Deserializer, Serialize};
use serde_json::Value;
use std::collections::HashMap;

pub trait EnrichrAPITrait: Send + Sync {
    async fn send_genes(
        &mut self,
        gene_list: &[String],
        libraries: &[String],
        send_background: bool,
    ) -> Result<(), APIFailure>;
    async fn enrich(&mut self, library_name: &str) -> Result<EnrichrResult, APIFailure>;
    fn get_short_id(&self) -> Option<String>;
}

impl EnrichrAPITrait for EnrichrAPI {
    async fn send_genes(
        &mut self,
        gene_list: &[String],
        libraries: &[String],
        send_background: bool,
    ) -> Result<(), APIFailure> {
        self.send_genes(gene_list, libraries, send_background).await
    }
    async fn enrich(&mut self, library_name: &str) -> Result<EnrichrResult, APIFailure> {
        let library_name = library_name.to_string();
        self.enrich(&library_name).await
    }
    fn get_short_id(&self) -> Option<String> {
        self.get_short_id()
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ListResponse {
    #[serde(rename = "userListId")]
    pub user_list_id: i32,
    #[serde(rename = "shortId")]
    pub short_id: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BackgroundResponse {
    pub backgroundid: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct EnrichrResponse {
    #[serde()]
    pub library: String,
    pub response: EnrichrResult,
}

impl<'de> Deserialize<'de> for EnrichrResponse {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let v = Value::deserialize(deserializer)?;
        let obj = v
            .as_object()
            .ok_or_else(|| serde::de::Error::custom("expected top-level object"))?;

        let (library, results) = obj
            .iter()
            .next()
            .ok_or_else(|| serde::de::Error::custom("expected one library key"))?;

        Ok(EnrichrResponse {
            library: library.clone(),
            response: EnrichrResult::new_from_json(results, library),
        })
    }
}

#[derive(Debug)]
pub struct APIFailure {
    pub message: String,
    short_id: Option<String>,
}

impl std::fmt::Display for APIFailure {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self.short_id {
            Some(ref id) => write!(
                f,
                "API Failure: {}\nResults can be found at: https://maayanlab.cloud/Enrichr/enrich?dataset={}",
                self.message, id
            ),
            None => write!(f, "API Failure: {}", self.message),
        }
    }
}

impl std::error::Error for APIFailure {}

impl APIFailure {
    pub fn new(message: String, short_id: Option<&String>) -> Self {
        let short_id = short_id.map(ToString::to_string);
        APIFailure { message, short_id }
    }
}

impl From<reqwest::Error> for APIFailure {
    fn from(e: reqwest::Error) -> Self {
        match e.url() {
            Some(url) => APIFailure::new(format!("{e} at {url}"), None),
            None => APIFailure::new(e.to_string(), None),
        }
    }
}

#[derive(Clone, Debug)]
pub struct EnrichrAPI {
    client: reqwest::Client,
    list_response: Option<ListResponse>,
    background_response: Option<BackgroundResponse>,
    enrich_response: HashMap<String, Option<Vec<Vec<Value>>>>,
    failure: HashMap<String, bool>,
    background: bool,
}

impl EnrichrAPI {
    const BASE_URL: &'static str = "https://maayanlab.cloud";
    const ENRICHR_URL: &'static str = "Enrichr";
    const LIST_URL: &'static str = "addList";
    const ENRICH_URL: &'static str = "enrich";
    const SPEEDRICH_URL: &'static str = "speedrichr/api";
    const BACKGROUND_URL: &'static str = "addbackground";
    const BACKGROUND_ENRICH_URL: &'static str = "backgroundenrich";

    pub fn new(background: Option<Vec<String>>) -> Self {
        match background {
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
        let url = if self.background {
            format!("{}/{}", Self::BASE_URL, Self::SPEEDRICH_URL)
        } else {
            format!("{}/{}", Self::BASE_URL, Self::ENRICHR_URL)
        };
        if background {
            format!("{}/{}", url, Self::BACKGROUND_URL)
        } else {
            format!("{}/{}", url, Self::LIST_URL)
        }
    }

    pub async fn send_genes(
        &mut self,
        gene_list: &[String],
        libraries: &[String],
        send_background: bool,
    ) -> Result<(), APIFailure> {
        let genes_text = gene_list.join("\n");
        let background_pre = if send_background { "Background " } else { "" };

        let form = if send_background {
            multipart::Form::new().text("background", genes_text)
        } else {
            multipart::Form::new()
                .text("list", genes_text)
                .text("description", "")
        };

        let response = self
            .client
            .post(self.get_list_url(send_background))
            .multipart(form)
            .send()
            .await?;
        if !response.status().is_success() {
            for lib in libraries {
                self.failure.insert(lib.clone(), true);
            }
            return Err(APIFailure::new(
                format!(
                    "{}Gene submission failed {}",
                    response.text().await.unwrap_or_default(),
                    background_pre
                ),
                None,
            ));
        }
        let response_text = response.text().await.map_err(|e| APIFailure {
            message: format!("Failed to read {background_pre}response: {e}"),
            short_id: None,
        })?;
        if send_background {
            self.background_response =
                Some(
                    serde_json::from_str(&response_text).map_err(|e| APIFailure {
                        message: format!("Failed to parse background response: {e}"),
                        short_id: None,
                    })?,
                );
        } else {
            self.list_response =
                Some(
                    serde_json::from_str(&response_text).map_err(|e| APIFailure {
                        message: format!("Failed to parse list response: {e}"),
                        short_id: None,
                    })?,
                );
        }
        Ok(())
    }

    pub fn get_short_id(&self) -> Option<String> {
        self.list_response.as_ref().map(|r| r.short_id.clone())
    }

    pub async fn enrich(&mut self, library_name: &String) -> Result<EnrichrResult, APIFailure> {
        if let Some(true) = self.failure.get(library_name) {
            return Ok(EnrichrResult::empty(library_name));
        }
        let short_id = &self.list_response.as_ref().unwrap().short_id;
        let user_list_id = &self
            .list_response
            .as_ref()
            .ok_or(APIFailure::new("No list response".to_string(), None))?
            .user_list_id;
        let request = if let Some(background_response) = self.background_response.clone() {
            let url = format!(
                "{}/{}/{}",
                Self::BASE_URL,
                Self::SPEEDRICH_URL,
                Self::BACKGROUND_ENRICH_URL
            );
            let form = multipart::Form::new()
                .text("backgroundid", background_response.backgroundid)
                .text("userListId", user_list_id.to_string())
                .text("backgroundType", library_name.clone());
            self.client.post(url).multipart(form)
        } else {
            let url = format!(
                "{}/{}/{}?userListId={}&backgroundType={}",
                Self::BASE_URL,
                Self::ENRICHR_URL,
                Self::ENRICH_URL,
                user_list_id,
                library_name
            );
            self.client.get(&url)
        };

        let response = request.send().await?;

        if !response.status().is_success() {
            self.failure.insert(library_name.clone(), true);
            return Err(APIFailure::new(
                format!("Enrichment request failed for {library_name}"),
                Some(short_id),
            ));
        }
        let response_text = match response.text().await {
            Ok(text) => text,
            Err(e) => {
                return Err(APIFailure::new(
                    format!("Failed to read enrichment response body for {library_name}: {e}"),
                    Some(short_id),
                ));
            }
        };

        let response_json: Value = match serde_json::from_str(&response_text.replace("NaN", "null"))
        {
            Ok(json) => json,
            Err(e) => {
                return Err(APIFailure::new(
                    format!(
                        "Failed to parse enrichment response for {library_name}: {e}\nResponse body:\n{response_text}"
                    ),
                    Some(short_id),
                ));
            }
        };
        match serde_json::from_value::<EnrichrResponse>(response_json) {
            Ok(enrichr_response) => {
                self.enrich_response.insert(
                    enrichr_response.library.clone(),
                    Some(
                        enrichr_response.response.get_all_rows_as_values(), // you'll need this getter
                    ),
                );
                Ok(enrichr_response.response)
            }
            Err(e) => {
                self.failure.insert(library_name.clone(), true);
                Err(APIFailure::new(
                    format!("Failed to deserialize enrichment response for {library_name}: {e}"),
                    Some(short_id),
                ))
            }
        }

        // if let Some(lib_data) = response_json.get(library_name).and_then(|v| v.as_array()) {
        //     if lib_data.is_empty() {
        //         self.failure.insert(library_name.clone(), true);
        //         return Ok(EnrichrResult::empty(library_name));
        //     } else {
        //         self.enrich_response.insert(
        //             library_name.clone(),
        //             Some(
        //                 lib_data
        //                     .iter()
        //                     .map(|v| v.as_array().map(|a| a.clone()).unwrap_or_default())
        //                     .collect(),
        //             ),
        //         );
        //     }
        // }
        //
        // if let Some(results) = response_json.get(library_name) {
        //     let r = EnrichrResult::new_from_json(results, library_name);
        //     return Ok(r);
        // }
        //
        // Ok(EnrichrResult::empty(library_name))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case::without_short_id("Something went wrong", None, "API Failure: Something went wrong")]
    #[case::with_short_id(
        "Enrichment request failed",
        Some("abc123"),
        "API Failure: Enrichment request failed\nResults can be found at: https://maayanlab.cloud/Enrichr/enrich?dataset=abc123"
    )]
    fn test_api_failure_display(
        #[case] message: &str,
        #[case] short_id: Option<&str>,
        #[case] expected: &str,
    ) {
        let short_id = short_id.map(String::from);
        let failure = APIFailure::new(message.to_string(), short_id.as_ref());

        assert_eq!(failure.to_string(), expected);
    }
}
