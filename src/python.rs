//! Python bindings for AnnoRefine using PyO3

#[cfg(feature = "python")]
use pyo3::prelude::*;
#[cfg(feature = "python")]
use pyo3::types::PyDict;


#[cfg(feature = "python")]
use crate::{
    gff3::parse_gff3_file,
    refinement::RefinementEngine,
    output::Gff3Writer,
    types::RefinementConfig,
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
        validate_splice_sites = true
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
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "RefinementConfig(min_coverage={}, min_splice_support={}, max_utr_extension={}, enable_novel_gene_detection={}, min_novel_gene_coverage={}, min_novel_gene_length={}, min_exon_length={}, validate_splice_sites={})",
            self.min_coverage,
            self.min_splice_support,
            self.max_utr_extension,
            self.enable_novel_gene_detection,
            self.min_novel_gene_coverage,
            self.min_novel_gene_length,
            self.min_exon_length,
            self.validate_splice_sites
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
        Some(rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build()
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("Failed to create thread pool: {}", e)))?)
    } else {
        None
    };

    // Use provided config or default
    let rust_config: RefinementConfig = config.unwrap_or_else(|| PyRefinementConfig::new(5, 3, 1000, false, 10, 300, 50, true)).into();

    // Load genome
    let genome = crate::fasta::parse_fasta_file(fasta_file)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Failed to load genome: {:?}", e)))?;

    // Parse GFF3
    let mut gene_models = parse_gff3_file(gff3_file)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Failed to parse GFF3: {:?}", e)))?;

    // Create refinement engine
    let engine = RefinementEngine::new(rust_config);

    // Refine annotations with interrupt handling
    let summary = py.allow_threads(|| {
        if let Some(pool) = thread_pool {
            pool.install(|| {
                engine.refine_gene_models(&mut gene_models, std::path::Path::new(bam_file), &genome)
            })
        } else {
            engine.refine_gene_models(&mut gene_models, std::path::Path::new(bam_file), &genome)
        }
    })
    .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("Refinement failed: {:?}", e)))?;

    // Write output
    let mut writer = Gff3Writer::new(output_file)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Failed to create output file: {:?}", e)))?;

    writer.write_gene_models(&gene_models)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Failed to write output: {:?}", e)))?;

    // Convert gene models to Python objects
    let py_gene_models: Vec<PyGeneModel> = gene_models.iter().map(|gene| {
        PyGeneModel {
            id: gene.id.clone(),
            chromosome: gene.chromosome.clone(),
            start: gene.start,
            end: gene.end,
            strand: format!("{:?}", gene.strand),
            has_structural_changes: gene.has_structural_changes,
            transcript_count: gene.transcripts.len(),
        }
    }).collect();

    // Create result dictionary
    let result = PyDict::new_bound(py);
    result.set_item("genes_processed", summary.genes_processed)?;
    result.set_item("genes_failed", summary.genes_failed)?;
    result.set_item("transcripts_with_structure_changes", summary.transcripts_with_structure_changes)?;
    result.set_item("transcripts_with_5utr_extension", summary.transcripts_with_5utr_extension)?;
    result.set_item("transcripts_with_3utr_extension", summary.transcripts_with_3utr_extension)?;
    result.set_item("novel_genes_detected", summary.novel_genes_detected)?;
    result.set_item("gene_models", py_gene_models.into_py(py))?;
    result.set_item("output_file", output_file)?;

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

/// Python module definition
#[cfg(feature = "python")]
#[pymodule]
fn _annorefine(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(refine_annotations, m)?)?;
    m.add_function(wrap_pyfunction!(version, m)?)?;
    m.add_function(wrap_pyfunction!(current_num_threads, m)?)?;
    m.add_function(wrap_pyfunction!(test_interruptible_operation, m)?)?;
    m.add_class::<PyRefinementConfig>()?;
    m.add_class::<PyGeneModel>()?;

    // Add module metadata
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add("__author__", "Jon Palmer")?;
    m.add("__description__", "Genome annotation refinement using RNA-seq data")?;

    Ok(())
}
