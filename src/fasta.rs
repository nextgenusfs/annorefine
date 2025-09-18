//! FASTA file parsing and genome sequence handling

use crate::types::{Genome, GenomeSequence, AnnoRefineError, Result};
use bio::io::fasta;
use std::path::Path;
use std::fs::File;
use std::io::BufReader;
use log::{info, debug};

/// Parse a FASTA file and return a Genome structure
pub fn parse_fasta_file<P: AsRef<Path>>(path: P) -> Result<Genome> {
    let path = path.as_ref();
    info!("Parsing FASTA file: {}", path.display());
    
    let file = File::open(path)
        .map_err(|e| AnnoRefineError::FastaParse(
            format!("Failed to open FASTA file {}: {}", path.display(), e)
        ))?;
    
    let reader = BufReader::new(file);
    let fasta_reader = fasta::Reader::new(reader);
    
    let mut genome = Genome::new();
    let mut sequence_count = 0;
    
    for result in fasta_reader.records() {
        let record = result.map_err(|e| AnnoRefineError::FastaParse(
            format!("Failed to parse FASTA record: {}", e)
        ))?;
        
        let id = record.id().to_string();
        let description = if record.desc().is_some() {
            Some(record.desc().unwrap().to_string())
        } else {
            None
        };
        
        let sequence = record.seq().to_vec();
        
        debug!("Loaded sequence: {} (length: {})", id, sequence.len());
        
        let genome_sequence = GenomeSequence {
            id: id.clone(),
            description,
            sequence,
        };
        
        genome.add_sequence(genome_sequence);
        sequence_count += 1;
    }
    
    info!("Successfully loaded {} sequences from FASTA file", sequence_count);
    
    if sequence_count == 0 {
        return Err(AnnoRefineError::FastaParse(
            "No sequences found in FASTA file".to_string()
        ));
    }
    
    Ok(genome)
}

/// Validate that a sequence contains only valid DNA bases
pub fn validate_dna_sequence(sequence: &[u8]) -> Result<()> {
    for (i, &base) in sequence.iter().enumerate() {
        match base.to_ascii_uppercase() {
            b'A' | b'T' | b'G' | b'C' | b'N' => continue,
            _ => {
                return Err(AnnoRefineError::FastaParse(
                    format!("Invalid DNA base '{}' at position {}", 
                        base as char, i + 1)
                ));
            }
        }
    }
    Ok(())
}

/// Get sequence statistics
pub fn get_sequence_stats(genome: &Genome) -> SequenceStats {
    let mut total_length = 0;
    let mut gc_count = 0;
    let mut n_count = 0;
    let mut sequence_count = 0;
    
    for sequence in genome.sequences.values() {
        sequence_count += 1;
        total_length += sequence.sequence.len();
        
        for &base in &sequence.sequence {
            match base.to_ascii_uppercase() {
                b'G' | b'C' => gc_count += 1,
                b'N' => n_count += 1,
                _ => {}
            }
        }
    }
    
    let gc_content = if total_length > 0 {
        (gc_count as f64 / total_length as f64) * 100.0
    } else {
        0.0
    };
    
    let n_content = if total_length > 0 {
        (n_count as f64 / total_length as f64) * 100.0
    } else {
        0.0
    };
    
    SequenceStats {
        sequence_count,
        total_length,
        gc_content,
        n_content,
    }
}

/// Statistics about genome sequences
#[derive(Debug)]
pub struct SequenceStats {
    pub sequence_count: usize,
    pub total_length: usize,
    pub gc_content: f64,
    pub n_content: f64,
}

impl std::fmt::Display for SequenceStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, 
            "Sequences: {}, Total length: {} bp, GC content: {:.2}%, N content: {:.2}%",
            self.sequence_count, self.total_length, self.gc_content, self.n_content
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;
    
    #[test]
    fn test_validate_dna_sequence() {
        assert!(validate_dna_sequence(b"ATGCN").is_ok());
        assert!(validate_dna_sequence(b"atgcn").is_ok());
        assert!(validate_dna_sequence(b"ATGCX").is_err());
    }
    
    #[test]
    fn test_parse_simple_fasta() {
        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, ">seq1 Test sequence 1").unwrap();
        writeln!(temp_file, "ATGCATGC").unwrap();
        writeln!(temp_file, ">seq2").unwrap();
        writeln!(temp_file, "GCTAGCTA").unwrap();
        
        let genome = parse_fasta_file(temp_file.path()).unwrap();
        
        assert_eq!(genome.sequences.len(), 2);
        assert!(genome.get_sequence("seq1").is_some());
        assert!(genome.get_sequence("seq2").is_some());
        
        let seq1 = genome.get_sequence("seq1").unwrap();
        assert_eq!(seq1.sequence, b"ATGCATGC");
        assert_eq!(seq1.description, Some("Test sequence 1".to_string()));
    }
}
