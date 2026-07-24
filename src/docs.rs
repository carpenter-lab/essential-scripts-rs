use clap::Parser;
use clap_markdown::{MarkdownOptions, help_markdown_custom};
use std::{fs, io};

fn md_opts() -> MarkdownOptions {
    MarkdownOptions::new().title("CLI Documentation".to_string())
}

pub fn write_docs_to_file<P: Parser>(path: &str) -> io::Result<()> {
    let md = help_markdown_custom::<P>(&md_opts());
    fs::write(path, md)
}
