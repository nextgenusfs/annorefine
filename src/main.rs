use anyhow::Result;
use clap::Parser;
use log::info;
use std::path::PathBuf;

use annorefine::{
    bam::get_alignment_stats,
    fasta::{get_sequence_stats, parse_fasta_file},
    gff3::parse_gff3_file,
    logging::init_logger,
    output::{validate_gene_models, Gff3Writer},
    refinement::RefinementEngine,
    types::RefinementConfig,
};

/// AnnoRefine: Genome annotation refinement using RNA-seq data
#[derive(Parser)]
#[command(name = "annorefine")]
#[command(about = "Refine genome annotations using RNA-seq alignments")]
#[command(version)]

struct Args {
    /// Input genome FASTA file
    #[arg(short = 'f', long = "fasta", value_name = "FILE")]
    fasta_file: PathBuf,

    /// Input GFF3 annotation file
    #[arg(short = 'g', long = "gff3", value_name = "FILE")]
    gff3_file: PathBuf,

    /// Input BAM alignment file (RNA-seq)
    #[arg(short = 'b', long = "bam", value_name = "FILE")]
    bam_file: PathBuf,

    /// Output refined GFF3 file
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

    /// Log file path (optional, logs all messages including detailed changes)
    #[arg(long = "log-file", value_name = "FILE")]
    log_file: Option<PathBuf>,

    /// Number of threads to use for parallel processing (default: number of CPU cores)
    #[arg(short = 't', long = "threads", value_name = "N")]
    threads: Option<usize>,

    /// Enable novel gene detection from RNA-seq evidence
    #[arg(long = "detect-novel-genes")]
    detect_novel_genes: bool,

    /// Minimum coverage for novel gene detection
    #[arg(long = "min-novel-coverage", value_name = "N", default_value = "10")]
    min_novel_coverage: u32,

    /// Minimum length for novel genes (bp)
    #[arg(long = "min-novel-length", value_name = "N", default_value = "300")]
    min_novel_length: u64,

    /// Minimum exon length for novel genes (bp)
    #[arg(long = "min-exon-length", value_name = "N", default_value = "50")]
    min_exon_length: u64,

    /// Disable splice site validation (allow non-canonical splice sites)
    #[arg(long = "no-splice-validation")]
    no_splice_validation: bool,
}

fn main() -> Result<()> {
    let args = Args::parse();

    // Initialize custom logger
    init_logger(args.verbose, args.log_file.as_deref())?;

    // Set up thread pool
    if let Some(threads) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .map_err(|e| anyhow::anyhow!("Failed to set thread pool: {}", e))?;
        info!("Using {} threads for parallel processing", threads);
    } else {
        let num_cpus = rayon::current_num_threads();
        info!(
            "Using {} threads for parallel processing (auto-detected)",
            num_cpus
        );
    }

    info!("Starting AnnoRefine v2025.1.0");
    info!("Input FASTA: {}", args.fasta_file.display());
    info!("Input GFF3: {}", args.gff3_file.display());
    info!("Input BAM: {}", args.bam_file.display());
    info!("Output GFF3: {}", args.output_file.display());

    // Validate input files exist
    validate_input_files(&args)?;

    // Run the main processing pipeline
    run_refinement_pipeline(&args)?;

    Ok(())
}

fn validate_input_files(args: &Args) -> Result<()> {
    if !args.fasta_file.exists() {
        anyhow::bail!("FASTA file does not exist: {}", args.fasta_file.display());
    }
    if !args.gff3_file.exists() {
        anyhow::bail!("GFF3 file does not exist: {}", args.gff3_file.display());
    }
    if !args.bam_file.exists() {
        anyhow::bail!("BAM file does not exist: {}", args.bam_file.display());
    }

    info!("All input files validated successfully");
    Ok(())
}

fn run_refinement_pipeline(args: &Args) -> Result<()> {
    info!("Starting annotation refinement pipeline");

    // 1. Parse genome FASTA file
    info!("Step 1: Loading genome sequences");
    let genome = parse_fasta_file(&args.fasta_file)?;
    let genome_stats = get_sequence_stats(&genome);
    info!("Genome loaded: {}", genome_stats);

    // 2. Parse GFF3 annotation file
    info!("Step 2: Loading gene annotations");
    let mut gene_models = parse_gff3_file(&args.gff3_file)?;
    info!("Loaded {} gene models", gene_models.len());

    // 3. Analyze BAM file
    info!("Step 3: Analyzing RNA-seq alignments");
    let bam_stats = get_alignment_stats(&args.bam_file)?;
    info!("BAM statistics: {}", bam_stats);

    // 4. Configure refinement parameters
    let config = RefinementConfig {
        min_coverage: args.min_coverage,
        min_splice_support: args.min_splice_support,
        max_utr_extension: args.max_utr_extension,
        enable_novel_gene_detection: args.detect_novel_genes,
        min_novel_gene_coverage: args.min_novel_coverage,
        min_novel_gene_length: args.min_novel_length,
        min_exon_length: args.min_exon_length,
        validate_splice_sites: !args.no_splice_validation,
    };

    // 5. Run refinement (now parallel)
    info!("Step 4: Refining gene models");
    info!("Using BAM file: {}", args.bam_file.display());
    let refinement_engine = RefinementEngine::new(config);
    let refinement_summary =
        refinement_engine.refine_gene_models(&mut gene_models, &args.bam_file, &genome)?;
    info!("Refinement complete: {}", refinement_summary);

    // 7. Validate refined gene models
    info!("Step 5: Validating refined gene models");
    validate_gene_models(&gene_models)?;

    // 8. Write output
    info!("Step 6: Writing refined annotations");
    let mut gff3_writer = Gff3Writer::new(&args.output_file)?;
    gff3_writer.write_gene_models(&gene_models)?;

    info!("Annotation refinement completed successfully!");
    info!("Output written to: {}", args.output_file.display());

    Ok(())
}
