//! AnnoRefine: Genome annotation refinement using RNA-seq data
//!
//! This library provides functionality to refine genome annotations by incorporating
//! RNA-seq alignment evidence to improve gene model predictions.

pub mod types;
pub mod fasta;
pub mod gff3;
pub mod bam;
pub mod refinement;
pub mod output;
pub mod translation;
pub mod logging;

#[cfg(feature = "python")]
pub mod python;

// Re-export Python module for PyO3
#[cfg(feature = "python")]
pub use python::*;

pub use types::*;

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_integration_fasta_gff3() {
        // Create test FASTA file
        let mut fasta_file = NamedTempFile::new().unwrap();
        writeln!(fasta_file, ">chr1 Test chromosome").unwrap();
        writeln!(fasta_file, "ATGCATGCATGCATGCATGCATGCATGCATGC").unwrap();

        // Create test GFF3 file
        let mut gff3_file = NamedTempFile::new().unwrap();
        writeln!(gff3_file, "##gff-version 3").unwrap();
        writeln!(gff3_file, "chr1\ttest\tgene\t10\t30\t.\t+\t.\tID=gene1").unwrap();
        writeln!(gff3_file, "chr1\ttest\tmRNA\t10\t30\t.\t+\t.\tID=mRNA1;Parent=gene1").unwrap();
        writeln!(gff3_file, "chr1\ttest\texon\t10\t30\t.\t+\t.\tID=exon1;Parent=mRNA1").unwrap();

        // Test parsing
        let genome = fasta::parse_fasta_file(fasta_file.path()).unwrap();
        let gene_models = gff3::parse_gff3_file(gff3_file.path()).unwrap();

        assert_eq!(genome.sequences.len(), 1);
        assert_eq!(gene_models.len(), 1);
        assert_eq!(gene_models[0].transcripts.len(), 1);
        assert_eq!(gene_models[0].transcripts[0].exons.len(), 1);

        // Test validation
        output::validate_gene_models(&gene_models).unwrap();
    }

    #[test]
    fn test_strand_parsing() {
        use std::str::FromStr;

        assert_eq!(Strand::from_str("+").unwrap(), Strand::Forward);
        assert_eq!(Strand::from_str("-").unwrap(), Strand::Reverse);
        assert_eq!(Strand::from_str(".").unwrap(), Strand::Unknown);
        assert!(Strand::from_str("x").is_err());
    }

    #[test]
    fn test_feature_type_parsing() {
        use std::str::FromStr;

        assert_eq!(FeatureType::from_str("gene").unwrap(), FeatureType::Gene);
        assert_eq!(FeatureType::from_str("mRNA").unwrap(), FeatureType::Mrna);
        assert_eq!(FeatureType::from_str("exon").unwrap(), FeatureType::Exon);
        assert_eq!(FeatureType::from_str("CDS").unwrap(), FeatureType::Cds);

        if let FeatureType::Other(s) = FeatureType::from_str("custom").unwrap() {
            assert_eq!(s, "custom");
        } else {
            panic!("Expected Other variant");
        }
    }
}
