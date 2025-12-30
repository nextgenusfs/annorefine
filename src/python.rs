//! Python bindings for AnnoRefine using PyO3

#[cfg(feature = "python")]
use pyo3::prelude::*;
#[cfg(feature = "python")]
use pyo3::types::PyDict;
#[cfg(feature = "python")]
use rust_htslib::bam::Read;

#[cfg(feature = "python")]
use crate::{
    bam::get_alignment_stats_with_config,
    bam2hints::Bam2HintsConverter,
    gff3::parse_gff3_file,
    hints_output::HintsWriter,
    output::Gff3Writer,
    refinement::RefinementEngine,
    types::{RefinementConfig, Bam2HintsConfig, LibraryType, StrandBias},
};

/// Python wrapper for RefinementConfig
#[cfg(feature = "python")]
#[pyclass]
#[derive(Clone)]
pub struct PyRefinementConfig {
    #[pyo3(get, set)]
    pub min_coverage: u32,
    #[pyo3(get, set)]
    pub min_splice_support: u32,
    #[pyo3(get, set)]
    pub max_utr_extension: u32,
    #[pyo3(get, set)]
    pub enable_novel_gene_detection: bool,
    #[pyo3(get, set)]
    pub min_novel_gene_coverage: u32,
    #[pyo3(get, set)]
    pub min_novel_gene_length: u64,
    #[pyo3(get, set)]
    pub min_exon_length: u64,
    #[pyo3(get, set)]
    pub validate_splice_sites: bool,
    #[pyo3(get, set)]
    pub strand_bias_threshold: f64,
    #[pyo3(get, set)]
    pub max_reads_for_strand_detection: u32,
    #[pyo3(get, set)]
    pub library_type: String,
}

#[cfg(feature = "python")]
#[pymethods]
impl PyRefinementConfig {
    #[new]
    #[pyo3(signature = (
        min_coverage = 5,
        min_splice_support = 3,
        max_utr_extension = 1000,
        enable_novel_gene_detection = false,
        min_novel_gene_coverage = 10,
        min_novel_gene_length = 300,
        min_exon_length = 50,
        validate_splice_sites = true,
        strand_bias_threshold = 0.65,
        max_reads_for_strand_detection = 10000,
        library_type = "auto".to_string()
    ))]
    pub fn new(
        min_coverage: u32,
        min_splice_support: u32,
        max_utr_extension: u32,
        enable_novel_gene_detection: bool,
        min_novel_gene_coverage: u32,
        min_novel_gene_length: u64,
        min_exon_length: u64,
        validate_splice_sites: bool,
        strand_bias_threshold: f64,
        max_reads_for_strand_detection: u32,
        library_type: String,
    ) -> Self {
        Self {
            min_coverage,
            min_splice_support,
            max_utr_extension,
            enable_novel_gene_detection,
            min_novel_gene_coverage,
            min_novel_gene_length,
            min_exon_length,
            validate_splice_sites,
            strand_bias_threshold,
            max_reads_for_strand_detection,
            library_type,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "RefinementConfig(min_coverage={}, min_splice_support={}, max_utr_extension={}, enable_novel_gene_detection={}, min_novel_gene_coverage={}, min_novel_gene_length={}, min_exon_length={}, validate_splice_sites={}, strand_bias_threshold={}, max_reads_for_strand_detection={})",
            self.min_coverage,
            self.min_splice_support,
            self.max_utr_extension,
            self.enable_novel_gene_detection,
            self.min_novel_gene_coverage,
            self.min_novel_gene_length,
            self.min_exon_length,
            self.validate_splice_sites,
            self.strand_bias_threshold,
            self.max_reads_for_strand_detection
        )
    }
}

#[cfg(feature = "python")]
impl From<PyRefinementConfig> for RefinementConfig {
    fn from(py_config: PyRefinementConfig) -> Self {
        RefinementConfig {
            min_coverage: py_config.min_coverage,
            min_splice_support: py_config.min_splice_support,
            max_utr_extension: py_config.max_utr_extension,
            enable_novel_gene_detection: py_config.enable_novel_gene_detection,
            min_novel_gene_coverage: py_config.min_novel_gene_coverage,
            min_novel_gene_length: py_config.min_novel_gene_length,
            min_exon_length: py_config.min_exon_length,
            validate_splice_sites: py_config.validate_splice_sites,
            strand_bias_threshold: py_config.strand_bias_threshold,
            max_reads_for_strand_detection: py_config.max_reads_for_strand_detection,
            library_type: py_config.library_type.parse().unwrap_or(crate::types::LibraryType::Auto),
        }
    }
}

/// Python wrapper for Bam2HintsConfig
#[cfg(feature = "python")]
#[pyclass]
#[derive(Clone)]
pub struct PyBam2HintsConfig {
    #[pyo3(get, set)]
    pub priority: u32,
    #[pyo3(get, set)]
    pub max_gap_len: u32,
    #[pyo3(get, set)]
    pub min_intron_len: u32,
    #[pyo3(get, set)]
    pub max_intron_len: u32,
    #[pyo3(get, set)]
    pub min_end_block_len: u32,
    #[pyo3(get, set)]
    pub max_query_gap_len: u32,
    #[pyo3(get, set)]
    pub exonpart_cutoff: u32,
    #[pyo3(get, set)]
    pub source: String,
    #[pyo3(get, set)]
    pub introns_only: bool,
    #[pyo3(get, set)]
    pub no_multiplicity: bool,
    #[pyo3(get, set)]
    pub remove_redundant: bool,
    #[pyo3(get, set)]
    pub max_coverage: u32,
    #[pyo3(get, set)]
    pub splice_sites_on: bool,
    #[pyo3(get, set)]
    pub truncated_splice_sites: bool,
    #[pyo3(get, set)]
    pub score: f64,
    #[pyo3(get, set)]
    pub max_gene_len: u32,
    #[pyo3(get, set)]
    pub library_type: String,
}

#[cfg(feature = "python")]
#[pymethods]
impl PyBam2HintsConfig {
    #[new]
    #[pyo3(signature = (
        library_type,
        priority = 4,
        max_gap_len = 14,
        min_intron_len = 32,
        max_intron_len = 350000,
        min_end_block_len = 8,
        max_query_gap_len = 5,
        exonpart_cutoff = 10,
        source = "E".to_string(),
        introns_only = false,
        no_multiplicity = false,
        remove_redundant = false,
        max_coverage = 0,
        splice_sites_on = false,
        truncated_splice_sites = false,
        score = 0.0,
        max_gene_len = 400000
    ))]
    pub fn new(
        library_type: String,
        priority: u32,
        max_gap_len: u32,
        min_intron_len: u32,
        max_intron_len: u32,
        min_end_block_len: u32,
        max_query_gap_len: u32,
        exonpart_cutoff: u32,
        source: String,
        introns_only: bool,
        no_multiplicity: bool,
        remove_redundant: bool,
        max_coverage: u32,
        splice_sites_on: bool,
        truncated_splice_sites: bool,
        score: f64,
        max_gene_len: u32,
    ) -> Self {
        Self {
            priority,
            max_gap_len,
            min_intron_len,
            max_intron_len,
            min_end_block_len,
            max_query_gap_len,
            exonpart_cutoff,
            source,
            introns_only,
            no_multiplicity,
            remove_redundant,
            max_coverage,
            splice_sites_on,
            truncated_splice_sites,
            score,
            max_gene_len,
            library_type,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "Bam2HintsConfig(priority={}, max_gap_len={}, min_intron_len={}, max_intron_len={}, source='{}', introns_only={}, splice_sites_on={}, max_coverage={})",
            self.priority,
            self.max_gap_len,
            self.min_intron_len,
            self.max_intron_len,
            self.source,
            self.introns_only,
            self.splice_sites_on,
            self.max_coverage
        )
    }
}

#[cfg(feature = "python")]
impl From<PyBam2HintsConfig> for Bam2HintsConfig {
    fn from(py_config: PyBam2HintsConfig) -> Self {
        // Parse library type string
        let library_type = match py_config.library_type.to_uppercase().as_str() {
            "FR" => LibraryType::PairedFR,
            "RF" => LibraryType::PairedRF,
            "UU" => LibraryType::PairedUnstranded,
            _ => LibraryType::PairedUnstranded, // Default to unstranded if invalid
        };

        // Determine strand bias from library type
        let strand_bias = match library_type {
            LibraryType::PairedFR => StrandBias::ForwardStranded,
            LibraryType::PairedRF => StrandBias::ReverseStranded,
            _ => StrandBias::Unstranded,
        };

        Bam2HintsConfig {
            priority: py_config.priority,
            max_gap_len: py_config.max_gap_len,
            min_intron_len: py_config.min_intron_len,
            max_intron_len: py_config.max_intron_len,
            min_end_block_len: py_config.min_end_block_len,
            max_query_gap_len: py_config.max_query_gap_len,
            exonpart_cutoff: py_config.exonpart_cutoff,
            source: py_config.source,
            introns_only: py_config.introns_only,
            no_multiplicity: py_config.no_multiplicity,
            remove_redundant: py_config.remove_redundant,
            max_coverage: py_config.max_coverage,
            splice_sites_on: py_config.splice_sites_on,
            truncated_splice_sites: py_config.truncated_splice_sites,
            score: py_config.score,
            max_gene_len: py_config.max_gene_len,
            library_type,
            strand_bias,
        }
    }
}

/// Python wrapper for gene model information
#[cfg(feature = "python")]
#[pyclass]
#[derive(Clone)]
pub struct PyGeneModel {
    #[pyo3(get)]
    pub id: String,
    #[pyo3(get)]
    pub chromosome: String,
    #[pyo3(get)]
    pub start: u64,
    #[pyo3(get)]
    pub end: u64,
    #[pyo3(get)]
    pub strand: String,
    #[pyo3(get)]
    pub has_structural_changes: bool,
    #[pyo3(get)]
    pub transcript_count: usize,
}

#[cfg(feature = "python")]
#[pymethods]
impl PyGeneModel {
    fn __repr__(&self) -> String {
        format!(
            "GeneModel(id='{}', chromosome='{}', start={}, end={}, strand='{}', has_structural_changes={}, transcript_count={})",
            self.id, self.chromosome, self.start, self.end, self.strand, self.has_structural_changes, self.transcript_count
        )
    }
}

/// Main Python interface for AnnoRefine
#[cfg(feature = "python")]
#[pyfunction]
#[pyo3(signature = (fasta_file, gff3_file, bam_file, output_file, config = None, threads = None))]
pub fn refine_annotations<'py>(
    py: Python<'py>,
    fasta_file: &str,
    gff3_file: &str,
    bam_file: &str,
    output_file: &str,
    config: Option<PyRefinementConfig>,
    threads: Option<usize>,
) -> PyResult<Bound<'py, PyDict>> {
    // Set up threading - create a custom thread pool if specified
    let thread_pool = if let Some(num_threads) = threads {
        Some(
            rayon::ThreadPoolBuilder::new()
                .num_threads(num_threads)
                .build()
                .map_err(|e| {
                    PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                        "Failed to create thread pool: {}",
                        e
                    ))
                })?,
        )
    } else {
        None
    };

    // Use provided config or default
    let rust_config: RefinementConfig = config
        .unwrap_or_else(|| PyRefinementConfig::new(5, 3, 1000, false, 10, 300, 50, true, 0.65, 10000, "auto".to_string()))
        .into();

    // Load genome
    let genome = crate::fasta::parse_fasta_file(fasta_file).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Failed to load genome: {:?}", e))
    })?;

    // Parse GFF3
    let mut gene_models = parse_gff3_file(gff3_file).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Failed to parse GFF3: {:?}", e))
    })?;

    // Detect strand bias from BAM file
    let bam_stats = get_alignment_stats_with_config(
        bam_file,
        rust_config.strand_bias_threshold,
        rust_config.max_reads_for_strand_detection,
    )
    .map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Failed to analyze BAM file: {:?}", e))
    })?;

    // Create refinement engine
    let engine = RefinementEngine::new(rust_config);

    // Refine annotations with interrupt handling
    let summary = py
        .allow_threads(|| {
            if let Some(pool) = thread_pool {
                pool.install(|| {
                    engine.refine_gene_models_with_strand_bias(
                        &mut gene_models,
                        std::path::Path::new(bam_file),
                        &genome,
                        bam_stats.strand_bias,
                    )
                })
            } else {
                engine.refine_gene_models_with_strand_bias(
                    &mut gene_models,
                    std::path::Path::new(bam_file),
                    &genome,
                    bam_stats.strand_bias,
                )
            }
        })
        .map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("Refinement failed: {:?}", e))
        })?;

    // Write output
    let mut writer = Gff3Writer::new(output_file).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyIOError, _>(format!(
            "Failed to create output file: {:?}",
            e
        ))
    })?;

    writer.write_gene_models(&gene_models).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Failed to write output: {:?}", e))
    })?;

    // Convert gene models to Python objects
    let py_gene_models: Vec<PyGeneModel> = gene_models
        .iter()
        .map(|gene| PyGeneModel {
            id: gene.id.clone(),
            chromosome: gene.chromosome.clone(),
            start: gene.start,
            end: gene.end,
            strand: format!("{:?}", gene.strand),
            has_structural_changes: gene.has_structural_changes,
            transcript_count: gene.transcripts.len(),
        })
        .collect();

    // Create result dictionary
    let result = PyDict::new_bound(py);
    result.set_item("genes_processed", summary.genes_processed)?;
    result.set_item("genes_failed", summary.genes_failed)?;
    result.set_item(
        "transcripts_with_structure_changes",
        summary.transcripts_with_structure_changes,
    )?;
    result.set_item(
        "transcripts_with_5utr_extension",
        summary.transcripts_with_5utr_extension,
    )?;
    result.set_item(
        "transcripts_with_3utr_extension",
        summary.transcripts_with_3utr_extension,
    )?;
    result.set_item("novel_genes_detected", summary.novel_genes_detected)?;
    result.set_item("gene_models", py_gene_models.into_py(py))?;
    result.set_item("output_file", output_file)?;

    Ok(result)
}

/// Python interface for bam2hints conversion
#[cfg(feature = "python")]
#[pyfunction]
#[pyo3(signature = (bam_file, output_file, library_type, config = None, threads = None, contig = None, region = None))]
pub fn bam2hints_convert<'py>(
    py: Python<'py>,
    bam_file: &str,
    output_file: &str,
    library_type: &str,
    config: Option<PyBam2HintsConfig>,
    threads: Option<usize>,
    contig: Option<&str>,
    region: Option<(String, u64, u64)>,
) -> PyResult<Bound<'py, PyDict>> {
    // Clone config for later use
    let config_clone = config.clone();

    // Use provided config or create default with required library_type
    let rust_config: Bam2HintsConfig = config
        .unwrap_or_else(|| PyBam2HintsConfig::new(library_type.to_string(), 4, 14, 32, 350000, 8, 5, 10, "E".to_string(), false, false, false, 0, false, false, 0.0, 400000))
        .into();

    // Set up threading - create a custom thread pool if specified
    let thread_pool = if let Some(num_threads) = threads {
        Some(
            rayon::ThreadPoolBuilder::new()
                .num_threads(num_threads)
                .build()
                .map_err(|e| {
                    PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                        "Failed to create thread pool: {}",
                        e
                    ))
                })?,
        )
    } else {
        None
    };

    // Convert contig/region parameters
    let filter_contig = if let Some(ref r) = region {
        Some(r.0.clone())
    } else {
        contig.map(|s| s.to_string())
    };

    let filter_region = region.clone();

    // Process BAM file with interrupt handling and optional threading
    let (processed_count, hints_count, total_hints) = py.allow_threads(|| {
        // Set up thread pool if specified
        if let Some(pool) = thread_pool {
            pool.install(|| process_bam_for_hints(bam_file, output_file, rust_config, filter_contig, filter_region))
        } else {
            process_bam_for_hints(bam_file, output_file, rust_config, filter_contig, filter_region)
        }
    })
    .map_err(|e: anyhow::Error| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("BAM processing failed: {:?}", e))
    })?;

    // Create result dictionary
    let result = PyDict::new_bound(py);
    result.set_item("alignments_processed", processed_count)?;
    result.set_item("alignments_with_hints", hints_count)?;
    result.set_item("total_hints_generated", total_hints)?;
    result.set_item("output_file", output_file)?;
    result.set_item("config", config_clone.unwrap_or_else(|| PyBam2HintsConfig::new(library_type.to_string(), 4, 14, 32, 350000, 8, 5, 10, "E".to_string(), false, false, false, 0, false, false, 0.0, 400000)).into_py(py))?;

    Ok(result)
}

/// Get version information
#[cfg(feature = "python")]
#[pyfunction]
pub fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

/// Get the current number of threads used by rayon
#[cfg(feature = "python")]
#[pyfunction]
pub fn current_num_threads() -> usize {
    rayon::current_num_threads()
}

/// Test function that can be interrupted (for testing Ctrl+C handling)
#[cfg(feature = "python")]
#[pyfunction]
pub fn test_interruptible_operation(py: Python, duration_seconds: u64) -> PyResult<String> {
    let start = std::time::Instant::now();
    let duration = std::time::Duration::from_secs(duration_seconds);

    while start.elapsed() < duration {
        // Check for interrupts every 100ms
        py.check_signals()?;
        std::thread::sleep(std::time::Duration::from_millis(100));
    }

    Ok(format!("Completed {} second operation", duration_seconds))
}

/// Join multiple hint files and merge identical hints
///
/// This function reads hints from multiple GFF files, sorts them, and merges
/// identical hints by summing their multiplicity values. This is equivalent to
/// the Augustus `join_mult_hints.pl` script.
///
/// Hints are considered identical if they have the same:
/// - Chromosome (column 1)
/// - Feature type (column 3)
/// - Start position (column 4)
/// - End position (column 5)
/// - Strand (column 7)
/// - Frame (column 8)
///
/// When identical hints are found, their multiplicity values are summed.
/// Group attributes (group=, grp=, gp=) are removed from the output.
///
/// Args:
///     input_files: List of input GFF hint files to join
///     output_file: Path for output joined hints file
///     introns_only: If True, only output intron hints (useful for GeneMark)
///
/// Returns:
///     Dictionary with joining statistics:
///     - input_files: Number of input files processed
///     - total_input_hints: Total hints read from all input files
///     - output_hints: Number of hints in output file (after merging)
///     - output_file: Path to output file
///
/// Example:
///     >>> # Join hints from BAM and protein alignments
///     >>> result = annorefine.join_hints(
///     ...     input_files=["bam_hints.gff", "protein_hints.gff"],
///     ...     output_file="joined_hints.gff"
///     ... )
///     >>> print(f"Merged {result['total_input_hints']} hints into {result['output_hints']}")
///
///     >>> # Join hints and filter for GeneMark (introns only)
///     >>> result = annorefine.join_hints(
///     ...     input_files=["bam_hints.gff", "protein_hints.gff"],
///     ...     output_file="genemark_hints.gff",
///     ...     introns_only=True
///     ... )
#[cfg(feature = "python")]
#[pyfunction]
#[pyo3(signature = (input_files, output_file, introns_only = false))]
pub fn join_hints<'py>(
    py: Python<'py>,
    input_files: Vec<String>,
    output_file: &str,
    introns_only: bool,
) -> PyResult<Bound<'py, PyDict>> {
    use std::io::{BufRead, BufReader, Write};
    use std::fs::File;

    // Read all hints from all input files
    let mut all_hints: Vec<String> = Vec::new();
    let mut total_input_hints = 0;

    for input_file in &input_files {
        let file = File::open(input_file)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(
                format!("Failed to open input file {}: {}", input_file, e)
            ))?;
        let reader = BufReader::new(file);

        for line in reader.lines() {
            let line = line.map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(
                format!("Failed to read line: {}", e)
            ))?;

            // Skip empty lines and comments
            if line.trim().is_empty() || line.starts_with('#') {
                continue;
            }

            total_input_hints += 1;
            all_hints.push(line);
        }
    }

    // Sort hints by: chromosome, feature type, start, end, strand, frame
    // This matches the sorting requirement of join_mult_hints.pl
    all_hints.sort_by(|a, b| {
        let fields_a: Vec<&str> = a.split('\t').collect();
        let fields_b: Vec<&str> = b.split('\t').collect();

        if fields_a.len() < 8 || fields_b.len() < 8 {
            return std::cmp::Ordering::Equal;
        }

        // Compare: chromosome (0), feature (2), start (3), end (4), strand (6), frame (7)
        fields_a[0].cmp(fields_b[0])
            .then_with(|| fields_a[2].cmp(fields_b[2]))
            .then_with(|| {
                let start_a = fields_a[3].parse::<u64>().unwrap_or(0);
                let start_b = fields_b[3].parse::<u64>().unwrap_or(0);
                start_a.cmp(&start_b)
            })
            .then_with(|| {
                let end_a = fields_a[4].parse::<u64>().unwrap_or(0);
                let end_b = fields_b[4].parse::<u64>().unwrap_or(0);
                end_a.cmp(&end_b)
            })
            .then_with(|| fields_a[6].cmp(fields_b[6]))
            .then_with(|| fields_a[7].cmp(fields_b[7]))
    });

    // Merge identical hints and sum multiplicities
    let mut merged_hints: Vec<(Vec<String>, u32)> = Vec::new();
    let mut last_fields: Option<Vec<String>> = None;
    let mut last_mult: u32 = 0;

    for hint_line in all_hints {
        let fields: Vec<String> = hint_line.split('\t').map(|s| s.to_string()).collect();

        if fields.len() < 9 {
            continue;
        }

        // Extract multiplicity from attributes (column 8)
        let mut mult = 1u32;
        if let Some(mult_match) = fields[8].split(';')
            .find(|attr| attr.starts_with("mult="))
        {
            if let Some(mult_str) = mult_match.strip_prefix("mult=") {
                mult = mult_str.parse::<u32>().unwrap_or(1);
            }
        }

        // Check if this hint is identical to the last one
        // Compare: chromosome (0), feature (2), start (3), end (4), strand (6), frame (7)
        let is_identical = if let Some(ref lf) = last_fields {
            fields[0] == lf[0] && fields[2] == lf[2] &&
            fields[3] == lf[3] && fields[4] == lf[4] &&
            fields[6] == lf[6] && fields[7] == lf[7]
        } else {
            false
        };

        if is_identical {
            // Merge with last hint by adding multiplicities
            last_mult += mult;
        } else {
            // Save the previous hint if it exists
            if let Some(lf) = last_fields.take() {
                merged_hints.push((lf, last_mult));
            }
            // Start new hint
            last_fields = Some(fields);
            last_mult = mult;
        }
    }

    // Don't forget the last hint
    if let Some(lf) = last_fields {
        merged_hints.push((lf, last_mult));
    }

    // Write merged hints to output file
    let mut writer = File::create(output_file)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(
            format!("Failed to create output file: {}", e)
        ))?;

    // Write header
    writeln!(writer, "# joined hints from {} input files", input_files.len())
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(
            format!("Failed to write header: {}", e)
        ))?;

    let mut output_hints = 0;
    for (mut fields, mult) in merged_hints {
        // Filter by hint type if introns_only is set
        if introns_only && fields.len() > 2 && fields[2] != "intron" {
            continue;
        }

        // Remove group attributes and old mult attribute
        fields[8] = fields[8]
            .split(';')
            .filter(|attr| !attr.starts_with("mult=") &&
                           !attr.starts_with("group=") &&
                           !attr.starts_with("grp=") &&
                           !attr.starts_with("gp="))
            .collect::<Vec<&str>>()
            .join(";");

        // Add new mult attribute at the beginning
        if mult > 1 {
            fields[8] = format!("mult={};{}", mult, fields[8]);
        } else {
            // Even for mult=1, we should include it for consistency
            fields[8] = format!("mult=1;{}", fields[8]);
        }

        // Ensure attributes end with semicolon if they don't already
        if !fields[8].ends_with(';') {
            fields[8].push(';');
        }

        writeln!(writer, "{}", fields.join("\t"))
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(
                format!("Failed to write hint: {}", e)
            ))?;
        output_hints += 1;
    }

    // Create result dictionary
    let result = PyDict::new_bound(py);
    result.set_item("input_files", input_files.len())?;
    result.set_item("total_input_hints", total_input_hints)?;
    result.set_item("output_hints", output_hints)?;
    result.set_item("output_file", output_file)?;

    Ok(result)
}

/// Filter hints by type and other criteria
#[cfg(feature = "python")]
#[pyfunction]
#[pyo3(signature = (input_file, output_file, hint_types = None, min_multiplicity = 1, contig = None))]
pub fn filter_hints<'py>(
    py: Python<'py>,
    input_file: &str,
    output_file: &str,
    hint_types: Option<Vec<String>>,
    min_multiplicity: u32,
    contig: Option<&str>,
) -> PyResult<Bound<'py, PyDict>> {
    use std::io::{BufRead, BufReader, Write};
    use std::fs::File;

    // Parse hint types if provided
    let filter_types: Option<Vec<String>> = hint_types.map(|types| {
        types.iter().map(|t| t.to_lowercase()).collect()
    });

    let mut total_hints = 0;
    let mut filtered_hints = 0;
    let mut kept_hints = 0;

    py.allow_threads(|| -> anyhow::Result<()> {
        let input = File::open(input_file)
            .map_err(|e| anyhow::anyhow!("Failed to open input file: {}", e))?;
        let reader = BufReader::new(input);

        let output = File::create(output_file)
            .map_err(|e| anyhow::anyhow!("Failed to create output file: {}", e))?;
        let mut writer = std::io::BufWriter::new(output);

        for line in reader.lines() {
            let line = line.map_err(|e| anyhow::anyhow!("Failed to read line: {}", e))?;

            // Skip comments and empty lines
            if line.starts_with('#') || line.trim().is_empty() {
                writeln!(writer, "{}", line)
                    .map_err(|e| anyhow::anyhow!("Failed to write line: {}", e))?;
                continue;
            }

            total_hints += 1;

            // Parse GFF line
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 9 {
                // Invalid GFF line, skip
                filtered_hints += 1;
                continue;
            }

            let hint_contig = fields[0];
            let hint_type = fields[2].to_lowercase();
            let attributes = fields[8];

            // Filter by contig if specified
            if let Some(filter_contig) = contig {
                if hint_contig != filter_contig {
                    filtered_hints += 1;
                    continue;
                }
            }

            // Filter by hint type if specified
            if let Some(ref types) = filter_types {
                if !types.contains(&hint_type) {
                    filtered_hints += 1;
                    continue;
                }
            }

            // Filter by multiplicity if specified
            if min_multiplicity > 1 {
                // Parse multiplicity from attributes
                let mult = if let Some(mult_str) = attributes.split(';')
                    .find(|attr| attr.starts_with("mult="))
                {
                    mult_str.trim_start_matches("mult=")
                        .parse::<u32>()
                        .unwrap_or(1)
                } else {
                    1
                };

                if mult < min_multiplicity {
                    filtered_hints += 1;
                    continue;
                }
            }

            // Keep this hint
            writeln!(writer, "{}", line)
                .map_err(|e| anyhow::anyhow!("Failed to write line: {}", e))?;
            kept_hints += 1;
        }

        Ok(())
    })
    .map_err(|e: anyhow::Error| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("Filtering failed: {:?}", e))
    })?;

    // Create result dictionary
    let result = PyDict::new_bound(py);
    result.set_item("total_hints", total_hints)?;
    result.set_item("filtered_hints", filtered_hints)?;
    result.set_item("kept_hints", kept_hints)?;
    result.set_item("input_file", input_file)?;
    result.set_item("output_file", output_file)?;

    Ok(result)
}

/// Python module definition
#[cfg(feature = "python")]
#[pymodule]
fn _annorefine(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(refine_annotations, m)?)?;
    m.add_function(wrap_pyfunction!(bam2hints_convert, m)?)?;
    m.add_function(wrap_pyfunction!(join_hints, m)?)?;
    m.add_function(wrap_pyfunction!(filter_hints, m)?)?;
    m.add_function(wrap_pyfunction!(version, m)?)?;
    m.add_function(wrap_pyfunction!(current_num_threads, m)?)?;
    m.add_function(wrap_pyfunction!(test_interruptible_operation, m)?)?;
    m.add_class::<PyRefinementConfig>()?;
    m.add_class::<PyBam2HintsConfig>()?;
    m.add_class::<PyGeneModel>()?;

    // Add module metadata
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add("__author__", "Jon Palmer")?;
    m.add(
        "__description__",
        "Genome annotation refinement using RNA-seq data",
    )?;

    Ok(())
}

/// Helper function to process BAM file for hints (can be called with or without thread pool)
#[cfg(feature = "python")]
fn process_bam_for_hints(
    bam_file: &str,
    output_file: &str,
    rust_config: Bam2HintsConfig,
    filter_contig: Option<String>,
    filter_region: Option<(String, u64, u64)>,
) -> anyhow::Result<(u64, u64, usize)> {
    // Create converter with provided config (no auto-detection)
    let mut converter = Bam2HintsConverter::new(rust_config.clone());

    // Open BAM file - use IndexedReader if filtering by region
    let mut processed_count = 0;
    let mut hints_count = 0;

    if let Some((chr, start, end)) = filter_region {
        // Use indexed reader for region-based filtering
        let mut reader = rust_htslib::bam::IndexedReader::from_path(bam_file)
            .map_err(|e| anyhow::anyhow!("Failed to open BAM file (indexed): {}", e))?;

        let header = reader.header().clone();

        // Get target ID for chromosome
        let tid = header.target_names()
            .iter()
            .position(|name| *name == chr.as_bytes())
            .ok_or_else(|| anyhow::anyhow!("Chromosome '{}' not found in BAM header", chr))?;

        // Fetch region (convert to 0-based for BAM)
        reader.fetch((tid as u32, start.saturating_sub(1), end))
            .map_err(|e| anyhow::anyhow!("Failed to fetch region {}:{}-{}: {}", chr, start, end, e))?;

        // Process alignments in region
        for result in reader.records() {
            let record = result.map_err(|e| anyhow::anyhow!("Failed to read BAM record: {}", e))?;

            // Skip unmapped reads
            if record.is_unmapped() {
                continue;
            }

            // Skip secondary alignments and duplicates
            if record.is_secondary() || record.is_duplicate() {
                continue;
            }

            processed_count += 1;

            // Parse BAM record into RnaSeqAlignment
            match crate::bam::parse_bam_record(&record, &header) {
                Ok(alignment) => {
                    match converter.process_alignment(&alignment) {
                        Ok(_) => hints_count += 1,
                        Err(_) => continue,
                    }
                }
                Err(_) => continue,
            }
        }
    } else {
        // Use regular reader for whole-file or contig-based filtering
        let mut reader = rust_htslib::bam::Reader::from_path(bam_file)
            .map_err(|e| anyhow::anyhow!("Failed to open BAM file: {}", e))?;

        let header = reader.header().clone();

        // If filtering by contig, use indexed reader with fetch
        if let Some(ref contig) = filter_contig {
            // Need to reopen as indexed reader
            drop(reader);
            let mut indexed_reader = rust_htslib::bam::IndexedReader::from_path(bam_file)
                .map_err(|e| anyhow::anyhow!("Failed to open BAM file (indexed): {}", e))?;

            let header = indexed_reader.header().clone();

            // Get target ID for chromosome
            let tid = header.target_names()
                .iter()
                .position(|name| *name == contig.as_bytes())
                .ok_or_else(|| anyhow::anyhow!("Contig '{}' not found in BAM header", contig))?;

            // Fetch entire contig
            indexed_reader.fetch(tid as u32)
                .map_err(|e| anyhow::anyhow!("Failed to fetch contig {}: {}", contig, e))?;

            // Process alignments in contig
            for result in indexed_reader.records() {
                let record = result.map_err(|e| anyhow::anyhow!("Failed to read BAM record: {}", e))?;

                if record.is_unmapped() || record.is_secondary() || record.is_duplicate() {
                    continue;
                }

                processed_count += 1;

                match crate::bam::parse_bam_record(&record, &header) {
                    Ok(alignment) => {
                        match converter.process_alignment(&alignment) {
                            Ok(_) => hints_count += 1,
                            Err(_) => continue,
                        }
                    }
                    Err(_) => continue,
                }
            }
        } else {
            // Process all alignments
            for result in reader.records() {
                let record = result.map_err(|e| anyhow::anyhow!("Failed to read BAM record: {}", e))?;

                if record.is_unmapped() || record.is_secondary() || record.is_duplicate() {
                    continue;
                }

                processed_count += 1;

                match crate::bam::parse_bam_record(&record, &header) {
                    Ok(alignment) => {
                        match converter.process_alignment(&alignment) {
                            Ok(_) => hints_count += 1,
                            Err(_) => continue,
                        }
                    }
                    Err(_) => continue,
                }
            }
        }
    }

    // Get all hints
    let hints = converter.get_hints();

    // Get chromosome order from BAM header for deterministic sorting
    let bam_reader = rust_htslib::bam::Reader::from_path(bam_file)
        .map_err(|e| anyhow::anyhow!("Failed to open BAM file for header: {}", e))?;
    let header = bam_reader.header();
    let chromosomes: Vec<String> = header
        .target_names()
        .iter()
        .map(|name| String::from_utf8_lossy(name).to_string())
        .collect();

    // Write hints to output file
    let output_file_handle = std::fs::File::create(output_file)
        .map_err(|e| anyhow::anyhow!("Failed to create output file: {}", e))?;
    let mut writer = HintsWriter::new(output_file_handle);

    // Write header (matching CLI behavior)
    let config_summary = format!(
        "priority={}, source={}, introns_only={}, max_coverage={}",
        rust_config.priority, rust_config.source, rust_config.introns_only, rust_config.max_coverage
    );
    writer.write_header(bam_file, &config_summary)
        .map_err(|e| anyhow::anyhow!("Failed to write header: {:?}", e))?;

    // Sort hints by chromosome order (matching CLI behavior)
    writer.write_hints_sorted_by_chromosome_order(hints, &chromosomes)
        .map_err(|e| anyhow::anyhow!("Failed to write hints: {:?}", e))?;

    Ok((processed_count, hints_count, writer.hints_written))
}
