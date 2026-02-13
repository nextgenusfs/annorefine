use anyhow::Result;
use clap::Parser;
use log::info;
use std::path::PathBuf;

use annorefine::{
    bam::{get_alignment_stats_with_gene_models, get_basic_alignment_stats, AlignmentStats},
    fasta::{get_sequence_stats, parse_fasta_file},
    gff3::parse_gff3_file,
    logging::init_logger,
    output::{validate_gene_models, Gff3Writer},
    refinement::RefinementEngine,
    types::{AnnoRefineError, LibraryType, RefinementConfig, StrandBias},
};

/// Convert detected strand bias to library type (assumes paired-end for stranded data)
fn convert_strand_bias_to_library_type(stats: &AlignmentStats) -> Result<LibraryType> {
    match stats.strand_bias {
        StrandBias::Unstranded => Ok(LibraryType::PairedUnstranded), // Assume paired-end unstranded
        StrandBias::ForwardStranded => Ok(LibraryType::PairedFR),    // FR library
        StrandBias::ReverseStranded => Ok(LibraryType::PairedRF),    // RF library
    }
}

/// Extend and refine UTRs using RNA-seq data
#[derive(Parser)]
pub struct UtrsCommand {
    /// Input FASTA file containing genome sequences
    #[arg(short = 'f', long = "fasta", value_name = "FILE")]
    fasta_file: PathBuf,

    /// Input GFF3 file containing gene annotations
    #[arg(short = 'g', long = "gff3", value_name = "FILE")]
    gff3_file: PathBuf,

    /// Input BAM file containing RNA-seq alignments
    #[arg(short = 'b', long = "bam", value_name = "FILE")]
    bam_file: PathBuf,

    /// Output GFF3 file for refined annotations
    #[arg(short = 'o', long = "output", value_name = "FILE")]
    output_file: PathBuf,

    /// Minimum coverage threshold for UTR extension
    #[arg(long = "min-coverage", default_value = "5")]
    min_coverage: u32,

    /// Minimum number of supporting reads for splice junction
    #[arg(long = "min-splice-support", default_value = "3")]
    min_splice_support: u32,

    /// Maximum UTR extension length (bp)
    #[arg(long = "max-utr-extension", default_value = "1000")]
    max_utr_extension: u32,

    /// Verbose output (shows warnings and debug info)
    #[arg(short = 'v', long = "verbose")]
    verbose: bool,

    /// Enable detection of novel genes from RNA-seq evidence
    #[arg(long = "detect-novel-genes")]
    detect_novel_genes: bool,

    /// Minimum coverage for novel gene detection
    #[arg(long = "min-novel-coverage", default_value = "10")]
    min_novel_coverage: u32,

    /// Minimum length for novel genes (bp)
    #[arg(long = "min-novel-length", default_value = "300")]
    min_novel_length: u32,

    /// Minimum exon length (bp)
    #[arg(long = "min-exon-length", default_value = "50")]
    min_exon_length: u32,

    /// Number of threads to use for parallel processing
    #[arg(short = 't', long = "threads")]
    threads: Option<usize>,

    /// Skip splice site validation
    #[arg(long = "no-splice-validation")]
    no_splice_validation: bool,

    /// Strand bias threshold for detecting stranded RNA-seq (0.5-1.0)
    #[arg(
        long = "strand-bias-threshold",
        value_name = "THRESHOLD",
        default_value = "0.65"
    )]
    strand_bias_threshold: f64,

    /// Library strandedness specification
    /// Options: auto, F, R, U, RF, FR, UU
    /// auto = auto-detect, F = single-end forward, R = single-end reverse, U = single-end unstranded
    /// RF = paired-end reverse/forward, FR = paired-end forward/reverse, UU = paired-end unstranded
    #[arg(long = "stranded", value_name = "TYPE", default_value = "auto")]
    stranded: String,

    /// Maximum reads to sample for strand detection
    #[arg(
        long = "max-reads-strand-detection",
        value_name = "N",
        default_value = "10000"
    )]
    max_reads_strand_detection: u32,
}

impl UtrsCommand {
    pub fn run(self) -> Result<()> {
        // Initialize logging
        init_logger(self.verbose, None)?;

        info!("Starting annorefine utrs v{}", env!("CARGO_PKG_VERSION"));
        info!("Input FASTA: {}", self.fasta_file.display());
        info!("Input GFF3: {}", self.gff3_file.display());
        info!("Input BAM: {}", self.bam_file.display());
        info!("Output GFF3: {}", self.output_file.display());

        // Set thread count
        if let Some(threads) = self.threads {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build_global()
                .map_err(|e| {
                    AnnoRefineError::InvalidInput(format!("Failed to set thread count: {}", e))
                })?;
            info!("Using {} threads for parallel processing", threads);
        } else {
            let num_threads = rayon::current_num_threads();
            info!(
                "Using {} threads for parallel processing (auto-detected)",
                num_threads
            );
        }

        // Validate input files exist
        for file in [&self.fasta_file, &self.gff3_file, &self.bam_file] {
            if !file.exists() {
                return Err(AnnoRefineError::InvalidInput(format!(
                    "File not found: {}",
                    file.display()
                ))
                .into());
            }
        }
        info!("All input files validated successfully");

        info!("Starting annotation refinement pipeline");

        // 1. Load genome sequences
        info!("Step 1: Loading genome sequences");
        let genome = parse_fasta_file(&self.fasta_file)?;
        let genome_stats = get_sequence_stats(&genome);
        info!("Loaded genome: {}", genome_stats);

        // 2. Parse gene models
        info!("Step 2: Parsing gene models");
        let gene_models = parse_gff3_file(&self.gff3_file)?;
        info!("Loaded {} gene models", gene_models.len());

        // Parse library type specification
        let library_type: LibraryType = self.stranded.parse().map_err(|e| {
            AnnoRefineError::InvalidInput(format!("Invalid --stranded option: {}", e))
        })?;

        info!("Library type specification: {}", library_type);

        // 3. Analyze BAM file with strand detection (auto-detect if needed)
        let (_bam_stats, final_library_type) = if library_type == LibraryType::Auto {
            info!("Step 3: Auto-detecting RNA-seq library type with gene-model-based strand detection");
            let stats = get_alignment_stats_with_gene_models(
                &self.bam_file,
                &gene_models,
                self.strand_bias_threshold,
                self.max_reads_strand_detection,
            )?;
            info!("BAM statistics: {}", stats);
            info!("Detected strand bias: {}", stats.strand_bias);

            // Convert detected strand bias to library type
            let detected_type = convert_strand_bias_to_library_type(&stats)?;
            info!("Auto-detected library type: {}", detected_type);
            (stats, detected_type)
        } else {
            info!("Step 3: Using specified library type: {}", library_type);
            // Still get basic BAM stats but don't use strand detection
            let stats = get_basic_alignment_stats(&self.bam_file)?;
            info!("BAM statistics: {}", stats);
            (stats, library_type)
        };

        // 4. Configure refinement engine
        let config = RefinementConfig {
            min_coverage: self.min_coverage,
            min_splice_support: self.min_splice_support,
            max_utr_extension: self.max_utr_extension,
            enable_novel_gene_detection: self.detect_novel_genes,
            min_novel_gene_coverage: self.min_novel_coverage,
            min_novel_gene_length: self.min_novel_length as u64,
            min_exon_length: self.min_exon_length as u64,
            validate_splice_sites: !self.no_splice_validation,
            strand_bias_threshold: self.strand_bias_threshold,
            max_reads_for_strand_detection: self.max_reads_strand_detection,
            library_type: final_library_type,
        };

        // 5. Run refinement (now parallel)
        info!("Step 4: Refining gene models");
        let engine = RefinementEngine::new(config);
        let mut gene_models_mut = gene_models; // Make mutable for refinement
        let summary = engine.refine_gene_models(&mut gene_models_mut, &self.bam_file, &genome)?;

        info!("Refinement complete: {}", summary);

        // 6. Validate and write output
        info!("Step 5: Writing refined annotations");
        validate_gene_models(&gene_models_mut)?;

        let mut writer = Gff3Writer::new(&self.output_file)?;
        writer.write_gene_models(&gene_models_mut)?;

        info!("Output written to: {}", self.output_file.display());
        Ok(())
    }
}
