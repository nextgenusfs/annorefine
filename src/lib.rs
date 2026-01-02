//! AnnoRefine: Genome annotation refinement using RNA-seq data
//!
//! This library provides functionality to refine genome annotations by incorporating
//! RNA-seq alignment evidence to improve gene model predictions.

pub mod bam;
pub mod bam2hints;
pub mod fasta;
pub mod gff3;
pub mod hints_output;
pub mod logging;
pub mod output;
pub mod refinement;
pub mod translation;
pub mod types;

#[cfg(feature = "python")]
pub mod python;

// Re-export Python module for PyO3
#[cfg(feature = "python")]
pub use python::*;

// Re-export main types for library usage
pub use types::*;
pub use gff3::*;
pub use fasta::*;
pub use refinement::*;
pub use output::*;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_types() {
        // Test that basic types can be created
        let config = types::RefinementConfig::default();
        assert_eq!(config.min_coverage, 5);
        assert_eq!(config.min_splice_support, 3);
    }

    #[test]
    fn test_strand_enum() {
        use types::Strand;

        let forward = Strand::Forward;
        let reverse = Strand::Reverse;
        let unknown = Strand::Unknown;

        assert_ne!(forward, reverse);
        assert_ne!(forward, unknown);
        assert_ne!(reverse, unknown);
    }

    #[test]
    fn test_version_info() {
        // Test that version info is available
        let version = env!("CARGO_PKG_VERSION");
        assert!(!version.is_empty());
        // Version should be in CalVer format (YYYY.MM.*)
        assert!(version.contains("202") && (version.contains("2025") || version.contains("2026")));
    }
}




