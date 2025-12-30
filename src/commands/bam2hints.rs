//! Command-line interface for bam2hints functionality

use annorefine::bam::{BamReader, parse_bam_record};
use annorefine::bam2hints::Bam2HintsConverter;
use annorefine::hints_output::{HintsWriter, HintsStatistics};
use annorefine::logging::init_logger;
use annorefine::types::{Bam2HintsConfig, GenomicInterval, Strand, StrandBias, LibraryType};
use anyhow::Context;
use clap::Parser;
use log::{info, debug, warn};
use rayon::prelude::*;
use rust_htslib::bam::Read;
use std::fs::File;
use std::io::{Write, BufWriter};
use std::path::PathBuf;


/// Convert BAM alignments to Augustus hints
#[derive(Parser, Debug)]
#[command(name = "bam2hints")]
#[command(about = "Convert mRNA-to-genome alignments in BAM format into hints file in GFF format")]
pub struct Bam2HintsCommand {
    /// Input BAM file (must be sorted by target sequence names and coordinates)
    #[arg(short = 'i', long = "bam", value_name = "FILE", required = true)]
    pub input: PathBuf,

    /// Output GFF hints file
    #[arg(short = 'o', long = "output", value_name = "FILE", required = true)]
    pub output: PathBuf,

    /// Priority of hint group
    #[arg(short = 'p', long = "priority", default_value = "4")]
    pub priority: u32,

    /// Gaps at most this length are simply closed
    #[arg(short = 'g', long = "maxgaplen", default_value = "14")]
    pub max_gap_len: u32,

    /// Alignments with gaps shorter than this and longer than maxgaplen are discarded
    #[arg(short = 'm', long = "minintronlen", default_value = "32")]
    pub min_intron_len: u32,

    /// Alignments with longer gaps are discarded
    #[arg(short = 'M', long = "maxintronlen", default_value = "350000")]
    pub max_intron_len: u32,

    /// Minimum length of a 'dangling' exon
    #[arg(short = 'b', long = "MinEndBlockLen", default_value = "8")]
    pub min_end_block_len: u32,

    /// Maximum length of gap in query (cDNA) sequence
    #[arg(short = 'q', long = "maxQgaplen", default_value = "5")]
    pub max_query_gap_len: u32,

    /// Compute exonpart, exon and splice site hints in addition to intron hints
    #[arg(short = 'x', long = "exonhints")]
    pub exon_hints: bool,

    /// This many bp are cut off of each exonpart hint at end of alignment
    #[arg(short = 'e', long = "ep_cutoff", default_value = "10")]
    pub exonpart_cutoff: u32,

    /// Source identifier
    #[arg(short = 's', long = "source", default_value = "E")]
    pub source: String,

    /// Only retrieve intron hints
    #[arg(short = 'I', long = "intronsonly")]
    pub introns_only: bool,

    /// Do not summarize multiple identical intron hints to a single one
    #[arg(short = 'n', long = "nomult")]
    pub no_multiplicity: bool,

    /// Only keep the strongest hint for a region
    #[arg(short = 'r', long = "remove_redundant")]
    pub remove_redundant: bool,

    /// Maximal number of hints at a given position (0: filtering deactivated)
    #[arg(short = 'C', long = "maxcoverage", default_value = "0")]
    pub max_coverage: u32,

    /// Force splice site (dss, ass) hints even when introns-only (default: auto-enabled when not introns-only)
    #[arg(short = 'S', long = "ssOn")]
    pub splice_sites_on: bool,

    /// Include splice sites hints from the ends of a truncated alignment
    #[arg(short = 'T', long = "trunkSS")]
    pub truncated_splice_sites: bool,

    /// Fill this number in the score column
    #[arg(short = 'v', long = "score", default_value = "0.0")]
    pub score: f64,

    /// Alignments that span more than this are ignored
    #[arg(short = 'G', long = "maxgenelen", default_value = "400000")]
    pub max_gene_len: u32,

    /// Number of threads to use for parallel processing
    #[arg(short = 't', long = "threads")]
    pub threads: Option<usize>,



    /// Library strandedness specification
    /// Options: FR, RF, UU
    /// FR = paired-end forward/reverse, RF = paired-end reverse/forward, UU = paired-end unstranded
    #[arg(long = "stranded", value_name = "TYPE", required = true)]
    pub stranded: String,

    /// Enable verbose output
    #[arg(long = "verbose")]
    pub verbose: bool,
}

impl Bam2HintsCommand {
    /// Execute the bam2hints command
    pub fn run(&self) -> anyhow::Result<()> {
        // Initialize logging using the same system as UTRs command
        init_logger(self.verbose, None)?;

        info!("Starting annorefine bam2hints v{}", env!("CARGO_PKG_VERSION"));
        info!("Input BAM: {}", self.input.display());
        info!("Output GFF: {}", self.output.display());

        // Set thread count - default to single-threaded
        let num_threads = self.threads.unwrap_or(1);
        if num_threads > 1 {
            rayon::ThreadPoolBuilder::new()
                .num_threads(num_threads)
                .build_global()
                .with_context(|| format!("Failed to set thread count: {}", num_threads))?;
            info!("Using {} threads for parallel processing", num_threads);
        } else {
            let num_threads = rayon::current_num_threads();
            info!("Using {} threads for parallel processing (auto-detected)", num_threads);
        }

        // Validate input parameters
        self.validate_parameters()?;
        info!("All input files validated successfully");

        // Input is now required, so we can use it directly
        let input_path = &self.input;

        // Step 1: Determine library type and strand bias (manual specification only)
        info!("Step 1: Using specified library strandedness");
        let (library_type, strand_bias) = match self.stranded.as_str() {
            "FR" => (LibraryType::PairedFR, StrandBias::ForwardStranded),
            "RF" => (LibraryType::PairedRF, StrandBias::ReverseStranded),
            "UU" => (LibraryType::PairedUnstranded, StrandBias::Unstranded),
            _ => {
                return Err(anyhow::anyhow!(
                    "Invalid stranded type '{}'. Valid options: FR, RF, UU",
                    self.stranded
                ));
            }
        };
        info!("Using specified library type: {}", library_type);
        info!("Library strandedness determined: {}", strand_bias);

        // Step 2: Create configuration and open files
        info!("Step 2: Setting up configuration and opening files");
        let config = self.create_config(library_type, strand_bias);
        debug!("Configuration: {:?}", config);

        let mut bam_reader = BamReader::new(input_path)
            .with_context(|| format!("Failed to open BAM file: {}", input_path.display()))?;
        info!("BAM file opened successfully");

        // Open output file
        let output: Box<dyn Write> = Box::new(BufWriter::new(File::create(&self.output)
            .with_context(|| format!("Failed to create output file: {}", self.output.display()))?));
        info!("Output file created: {}", self.output.display());

        // Step 3: Process BAM file and generate hints
        if num_threads > 1 {
            info!("Step 3: Processing BAM file with {} threads", num_threads);
            self.process_bam_file_parallel(&mut bam_reader, output, config)?;
        } else {
            info!("Step 3: Processing BAM file in single-threaded mode");
            self.process_bam_file_single_threaded(&mut bam_reader, output, config)?;
        }

        info!("BAM to hints conversion completed successfully");
        info!("Output written to: {}", self.output.display());
        Ok(())
    }

    /// Validate command-line parameters
    fn validate_parameters(&self) -> anyhow::Result<()> {
        if self.max_gap_len >= self.min_intron_len {
            return Err(anyhow::anyhow!(
                "maxgaplen ({}) must be less than minintronlen ({})",
                self.max_gap_len,
                self.min_intron_len
            ));
        }

        if self.min_intron_len > self.max_intron_len {
            return Err(anyhow::anyhow!(
                "minintronlen ({}) must be less than or equal to maxintronlen ({})",
                self.min_intron_len,
                self.max_intron_len
            ));
        }

        Ok(())
    }

    /// Create configuration from command-line arguments
    fn create_config(&self, library_type: LibraryType, strand_bias: StrandBias) -> Bam2HintsConfig {
        Bam2HintsConfig {
            priority: self.priority,
            max_gap_len: self.max_gap_len,
            min_intron_len: self.min_intron_len,
            max_intron_len: self.max_intron_len,
            min_end_block_len: self.min_end_block_len,
            max_query_gap_len: self.max_query_gap_len,
            exonpart_cutoff: self.exonpart_cutoff,
            source: self.source.clone(),
            introns_only: self.introns_only, // True if --intronsonly flag is set, false otherwise
            no_multiplicity: self.no_multiplicity,
            remove_redundant: self.remove_redundant,
            max_coverage: self.max_coverage,
            // Enable splice sites by default when not introns-only, unless explicitly disabled
            splice_sites_on: self.splice_sites_on || !self.introns_only,
            truncated_splice_sites: self.truncated_splice_sites,
            score: self.score,
            max_gene_len: self.max_gene_len,
            library_type,
            strand_bias,
        }
    }





    /// Process the BAM file in parallel and generate hints
    fn process_bam_file_parallel(
        &self,
        bam_reader: &mut BamReader,
        output: Box<dyn Write>,
        config: Bam2HintsConfig,
    ) -> anyhow::Result<()> {
        let mut hints_writer = HintsWriter::new(output);
        let _statistics = HintsStatistics::new();

        // Write header
        let input_name = self.input.display().to_string();
        let config_summary = format!(
            "priority={}, source={}, introns_only={}, max_coverage={}",
            config.priority, config.source, config.introns_only, config.max_coverage
        );
        hints_writer.write_header(&input_name, &config_summary)?;

        info!("Processing BAM alignments in parallel...");

        // Get all chromosome names from BAM header (preserving order)
        let chromosomes = self.get_chromosomes_from_bam(bam_reader)?;
        info!("Found {} chromosomes to process", chromosomes.len());

        // Process chromosomes in parallel but collect results in BAM header order
        let chromosome_results: anyhow::Result<Vec<_>> = chromosomes
            .par_iter()
            .map(|chromosome| {
                self.process_chromosome_parallel(chromosome, &config)
            })
            .collect();

        let chromosome_results = chromosome_results?;

        // Collect and merge all hints in BAM header order
        let mut all_hints = Vec::new();
        let mut total_processed = 0;
        let mut _total_with_hints = 0;

        // Process results in the same order as chromosomes (BAM header order)
        for (hints, processed, with_hints) in chromosome_results {
            all_hints.extend(hints);
            total_processed += processed;
            _total_with_hints += with_hints;
        }

        info!("Processed {} alignments total, {} generated hints", total_processed, all_hints.len());

        // Apply coverage filtering if enabled and write hints
        // Sort hints using BAM header chromosome order for deterministic output
        if config.max_coverage > 0 {
            hints_writer.write_hints_with_coverage_filter_ordered(all_hints, config.max_coverage, &chromosomes)?;
        } else {
            hints_writer.write_hints_sorted_by_chromosome_order(all_hints, &chromosomes)?;
        }

        Ok(())
    }

    /// Get chromosome names from BAM header
    fn get_chromosomes_from_bam(&self, bam_reader: &BamReader) -> anyhow::Result<Vec<String>> {
        let header = bam_reader.get_header();
        let mut chromosomes = Vec::new();

        for i in 0..header.target_count() {
            let chr_name = String::from_utf8_lossy(header.target_names()[i as usize]).to_string();
            chromosomes.push(chr_name);
        }

        Ok(chromosomes)
    }

    /// Process a single chromosome in parallel
    fn process_chromosome_parallel(
        &self,
        chromosome: &str,
        config: &Bam2HintsConfig,
    ) -> anyhow::Result<(Vec<annorefine::types::AugustusHint>, u64, u64)> {
        // Create a new BAM reader for this thread
        let input_path = &self.input;

        let mut thread_bam_reader = BamReader::new(input_path)
            .with_context(|| format!("Failed to open BAM file for chromosome {}", chromosome))?;

        let mut converter = Bam2HintsConverter::new(config.clone());
        let mut processed_count = 0u64;
        let mut hints_count = 0u64;

        debug!("Processing chromosome: {}", chromosome);

        // Create a region covering the entire chromosome
        // We'll get the chromosome length from the BAM header
        let chr_length = thread_bam_reader.get_chromosome_length(chromosome)
            .with_context(|| format!("Failed to get length for chromosome {}", chromosome))?;
        let region = GenomicInterval {
            chromosome: chromosome.to_string(),
            start: 1,
            end: chr_length,
            strand: Strand::Forward, // We'll process both strands
        };

        // Get all alignments for this chromosome
        match thread_bam_reader.get_alignments_in_region(&region) {
            Ok(alignments) => {
                debug!("Processing {} alignments from chromosome {}", alignments.len(), chromosome);

                for alignment in &alignments {
                    processed_count += 1;

                    if let Err(e) = converter.process_alignment(alignment) {
                        warn!("Failed to process alignment {}: {}", alignment.read_name, e);
                        continue;
                    }
                    hints_count += 1;
                }
            }
            Err(e) => {
                warn!("Failed to get alignments for chromosome {}: {}", chromosome, e);
                return Ok((Vec::new(), 0, 0));
            }
        }

        // Get generated hints
        let hints = converter.get_hints();
        debug!("Generated {} hints from chromosome {}", hints.len(), chromosome);

        Ok((hints, processed_count, hints_count))
    }

    /// Process the BAM file in single-threaded mode and generate hints
    fn process_bam_file_single_threaded(
        &self,
        _bam_reader: &mut BamReader,
        output: Box<dyn Write>,
        config: Bam2HintsConfig,
    ) -> anyhow::Result<()> {
        info!("Processing BAM file in single-threaded mode");

        // Create converter
        let mut converter = Bam2HintsConverter::new(config.clone());

        // Open BAM file using rust-htslib directly for record iteration
        let input_path = &self.input;
        let mut reader = rust_htslib::bam::Reader::from_path(input_path)
            .with_context(|| format!("Failed to open BAM file: {}", input_path.display()))?;

        // Get header before iterating
        let header = reader.header().clone();

        // Get chromosome order from BAM header for deterministic sorting
        let mut chromosomes = Vec::new();
        for i in 0..header.target_count() {
            let chr_name = String::from_utf8_lossy(header.target_names()[i as usize]).to_string();
            chromosomes.push(chr_name);
        }

        // Process all alignments
        let mut processed_count = 0;
        let mut hints_count = 0;
        info!("Starting to process BAM records...");

        for result in reader.records() {
            let record = result.with_context(|| "Failed to read BAM record")?;

            // Skip unmapped reads
            if record.is_unmapped() {
                continue;
            }

            // Skip secondary alignments and duplicates
            if record.is_secondary() || record.is_duplicate() {
                continue;
            }

            processed_count += 1;

            // Log progress every 100,000 records
            if processed_count % 100_000 == 0 {
                info!("Processed {} alignments so far...", processed_count);
            }

            // Parse BAM record into RnaSeqAlignment
            match parse_bam_record(&record, &header) {
                Ok(alignment) => {
                    match converter.process_alignment(&alignment) {
                        Ok(_) => hints_count += 1,
                        Err(_) => {
                            // Skip alignments that can't be processed
                            continue;
                        }
                    }
                }
                Err(_) => {
                    // Skip records that can't be parsed
                    continue;
                }
            }
        }

        info!("Processed {} alignments, {} generated hints", processed_count, hints_count);

        // Get generated hints
        let hints = converter.get_hints();
        info!("Total hints generated: {}", hints.len());

        // Write hints to output
        let mut hints_writer = HintsWriter::new(output);

        // Write header
        let config_summary = format!(
            "priority={}, source={}, introns_only={}, max_coverage={}",
            config.priority, config.source, config.introns_only, config.max_coverage
        );
        hints_writer.write_header(input_path.to_string_lossy().as_ref(), &config_summary)
            .with_context(|| "Failed to write hints header")?;

        // Write hints with coverage filter if specified, using BAM header chromosome order
        if config.max_coverage > 0 {
            hints_writer.write_hints_with_coverage_filter_ordered(hints, config.max_coverage, &chromosomes)
                .with_context(|| "Failed to write hints with coverage filter")?;
        } else {
            hints_writer.write_hints_sorted_by_chromosome_order(hints, &chromosomes)
                .with_context(|| "Failed to write hints")?;
        }

        Ok(())
    }
}
