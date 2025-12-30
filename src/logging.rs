//! Custom logging setup for AnnoRefine

use colored::*;
use log::{Level, LevelFilter, Metadata, Record};
use std::fs::OpenOptions;
use std::io::Write;
use std::path::Path;
use std::sync::Mutex;

/// Custom logger that handles both console and file output
pub struct AnnoRefineLogger {
    console_level: LevelFilter,
    file_writer: Option<Mutex<Box<dyn Write + Send>>>,
}

impl AnnoRefineLogger {
    pub fn new(verbose: bool, log_file: Option<&Path>) -> Result<Self, std::io::Error> {
        let console_level = if verbose {
            LevelFilter::Debug // Show all: ERROR, WARN, INFO, DEBUG
        } else {
            LevelFilter::Info // Show only: ERROR, INFO (no WARN, no DEBUG)
        };

        let file_writer = if let Some(log_path) = log_file {
            let file = OpenOptions::new()
                .create(true)
                .write(true)
                .truncate(true)
                .open(log_path)?;
            Some(Mutex::new(Box::new(file) as Box<dyn Write + Send>))
        } else {
            None
        };

        Ok(AnnoRefineLogger {
            console_level,
            file_writer,
        })
    }
}

impl log::Log for AnnoRefineLogger {
    fn enabled(&self, metadata: &Metadata) -> bool {
        // Always log to file if available, otherwise use console level
        self.file_writer.is_some() || metadata.level() <= self.console_level
    }

    fn log(&self, record: &Record) {
        if !self.enabled(record.metadata()) {
            return;
        }

        let timestamp = chrono::Utc::now().format("%H:%M:%S");
        let level = record.level();
        let target = record.target();
        let message = record.args();

        // Format message with colors for console
        let colored_level = match level {
            Level::Error => "ERROR".red().bold(),
            Level::Warn => "WARN".yellow().bold(),
            Level::Info => "INFO".green().bold(),
            Level::Debug => "DEBUG".blue().bold(),
            Level::Trace => "TRACE".purple().bold(),
        };

        let colored_message = format!(
            "[{} {} {}] {}",
            timestamp.to_string().dimmed(),
            colored_level,
            target.cyan(),
            message
        );

        // Plain message for file logging
        let plain_message = format!("[{} {} {}] {}", timestamp, level, target, message);

        // Log to console based on level (with colors)
        match level {
            Level::Error => eprintln!("{}", colored_message),
            Level::Warn => {
                if self.console_level >= LevelFilter::Debug {
                    eprintln!("{}", colored_message);
                }
            }
            Level::Info => {
                if self.console_level >= LevelFilter::Info {
                    println!("{}", colored_message);
                }
            }
            Level::Debug => {
                if self.console_level >= LevelFilter::Debug {
                    println!("{}", colored_message);
                }
            }
            Level::Trace => {
                if self.console_level >= LevelFilter::Trace {
                    println!("{}", colored_message);
                }
            }
        }

        // Log to file if available (all levels, no colors)
        if let Some(ref file_writer) = self.file_writer {
            if let Ok(mut writer) = file_writer.lock() {
                let _ = writeln!(writer, "{}", plain_message);
                let _ = writer.flush();
            }
        }
    }

    fn flush(&self) {
        if let Some(ref file_writer) = self.file_writer {
            if let Ok(mut writer) = file_writer.lock() {
                let _ = writer.flush();
            }
        }
    }
}

/// Initialize the custom logger
pub fn init_logger(verbose: bool, log_file: Option<&Path>) -> Result<(), anyhow::Error> {
    let logger = AnnoRefineLogger::new(verbose, log_file)
        .map_err(|e| anyhow::anyhow!("Failed to create logger: {}", e))?;

    log::set_boxed_logger(Box::new(logger))
        .map_err(|e| anyhow::anyhow!("Failed to set logger: {}", e))?;
    log::set_max_level(LevelFilter::Debug); // Allow all levels, filtering happens in the logger

    Ok(())
}

/// Log detailed gene model changes
pub fn log_gene_model_changes(
    gene_id: &str,
    transcript_id: &str,
    change_type: &str,
    details: &str,
) {
    log::warn!(
        target: "annorefine::changes",
        "CHANGE: Gene={}, Transcript={}, Type={}, Details={}",
        gene_id, transcript_id, change_type, details
    );
}

/// Log refinement statistics
pub fn log_refinement_stats(
    gene_id: &str,
    transcript_id: &str,
    original_exons: usize,
    refined_exons: usize,
    original_length: u64,
    refined_length: u64,
) {
    log::debug!(
        target: "annorefine::stats",
        "STATS: Gene={}, Transcript={}, Exons: {}→{}, Length: {}→{}",
        gene_id, transcript_id, original_exons, refined_exons, original_length, refined_length
    );
}

/// Log CDS validation results
pub fn log_cds_validation(gene_id: &str, transcript_id: &str, is_valid: bool, issues: &[String]) {
    if is_valid {
        log::debug!(
            target: "annorefine::cds",
            "CDS_VALID: Gene={}, Transcript={}",
            gene_id, transcript_id
        );
    } else {
        log::warn!(
            target: "annorefine::cds",
            "CDS_INVALID: Gene={}, Transcript={}, Issues={}",
            gene_id, transcript_id, issues.join("; ")
        );
    }
}

/// Log UTR extensions
pub fn log_utr_extension(
    gene_id: &str,
    transcript_id: &str,
    utr_type: &str,
    original_pos: u64,
    new_pos: u64,
    extension_length: u64,
) {
    log_gene_model_changes(
        gene_id,
        transcript_id,
        &format!("{}_UTR_EXTENSION", utr_type),
        &format!("{}→{} (+{}bp)", original_pos, new_pos, extension_length),
    );
}

/// Log splice junction refinements
pub fn log_splice_junction_refinement(
    gene_id: &str,
    transcript_id: &str,
    exon_index: usize,
    original_donor: u64,
    original_acceptor: u64,
    new_donor: u64,
    new_acceptor: u64,
    support_count: u32,
) {
    log_gene_model_changes(
        gene_id,
        transcript_id,
        "SPLICE_JUNCTION_REFINEMENT",
        &format!(
            "Exon{}: donor {}→{}, acceptor {}→{} (support={})",
            exon_index + 1,
            original_donor,
            new_donor,
            original_acceptor,
            new_acceptor,
            support_count
        ),
    );
}
