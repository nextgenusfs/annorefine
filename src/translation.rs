//! Translation and CDS validation using genetic code tables

use crate::types::{Genome, Result, Strand, Transcript};
use log::{debug, warn};
use std::collections::HashMap;

/// Standard genetic code table (NCBI table 1)
pub struct GeneticCode {
    codon_table: HashMap<String, char>,
    start_codons: Vec<String>,
    stop_codons: Vec<String>,
}

impl GeneticCode {
    /// Create the standard genetic code (NCBI table 1)
    pub fn standard() -> Self {
        let mut codon_table = HashMap::new();

        // Standard genetic code (table 1)
        let codons = [
            ("TTT", 'F'),
            ("TTC", 'F'),
            ("TTA", 'L'),
            ("TTG", 'L'),
            ("TCT", 'S'),
            ("TCC", 'S'),
            ("TCA", 'S'),
            ("TCG", 'S'),
            ("TAT", 'Y'),
            ("TAC", 'Y'),
            ("TAA", '*'),
            ("TAG", '*'),
            ("TGT", 'C'),
            ("TGC", 'C'),
            ("TGA", '*'),
            ("TGG", 'W'),
            ("CTT", 'L'),
            ("CTC", 'L'),
            ("CTA", 'L'),
            ("CTG", 'L'),
            ("CCT", 'P'),
            ("CCC", 'P'),
            ("CCA", 'P'),
            ("CCG", 'P'),
            ("CAT", 'H'),
            ("CAC", 'H'),
            ("CAA", 'Q'),
            ("CAG", 'Q'),
            ("CGT", 'R'),
            ("CGC", 'R'),
            ("CGA", 'R'),
            ("CGG", 'R'),
            ("ATT", 'I'),
            ("ATC", 'I'),
            ("ATA", 'I'),
            ("ATG", 'M'),
            ("ACT", 'T'),
            ("ACC", 'T'),
            ("ACA", 'T'),
            ("ACG", 'T'),
            ("AAT", 'N'),
            ("AAC", 'N'),
            ("AAA", 'K'),
            ("AAG", 'K'),
            ("AGT", 'S'),
            ("AGC", 'S'),
            ("AGA", 'R'),
            ("AGG", 'R'),
            ("GTT", 'V'),
            ("GTC", 'V'),
            ("GTA", 'V'),
            ("GTG", 'V'),
            ("GCT", 'A'),
            ("GCC", 'A'),
            ("GCA", 'A'),
            ("GCG", 'A'),
            ("GAT", 'D'),
            ("GAC", 'D'),
            ("GAA", 'E'),
            ("GAG", 'E'),
            ("GGT", 'G'),
            ("GGC", 'G'),
            ("GGA", 'G'),
            ("GGG", 'G'),
        ];

        for (codon, aa) in codons.iter() {
            codon_table.insert(codon.to_string(), *aa);
        }

        let start_codons = vec!["ATG".to_string()];
        let stop_codons = vec!["TAA".to_string(), "TAG".to_string(), "TGA".to_string()];

        Self {
            codon_table,
            start_codons,
            stop_codons,
        }
    }

    /// Translate a codon to amino acid
    pub fn translate_codon(&self, codon: &str) -> Option<char> {
        self.codon_table.get(&codon.to_uppercase()).copied()
    }

    /// Check if codon is a start codon
    pub fn is_start_codon(&self, codon: &str) -> bool {
        self.start_codons.contains(&codon.to_uppercase())
    }

    /// Check if codon is a stop codon
    pub fn is_stop_codon(&self, codon: &str) -> bool {
        self.stop_codons.contains(&codon.to_uppercase())
    }
}

/// CDS validation result
#[derive(Debug, Clone)]
pub struct CdsValidationResult {
    pub is_valid: bool,
    pub has_start_codon: bool,
    pub has_stop_codon: bool,
    pub has_premature_stop: bool,
    pub frame_preserved: bool,
    pub protein_sequence: Option<String>,
    pub issues: Vec<String>,
}

/// Validate CDS regions for a transcript
pub fn validate_cds(
    transcript: &Transcript,
    genome: &Genome,
    chromosome: &str,
    strand: Strand,
    genetic_code: &GeneticCode,
) -> Result<CdsValidationResult> {
    if transcript.cds_regions.is_empty() {
        return Ok(CdsValidationResult {
            is_valid: true, // Non-coding transcripts are valid
            has_start_codon: false,
            has_stop_codon: false,
            has_premature_stop: false,
            frame_preserved: true,
            protein_sequence: None,
            issues: vec!["No CDS regions (non-coding transcript)".to_string()],
        });
    }

    // Extract CDS sequence
    let cds_sequence = extract_cds_sequence(transcript, genome, chromosome, strand)?;

    // Validate the CDS
    validate_cds_sequence(&cds_sequence, genetic_code)
}

/// Extract the complete CDS sequence from transcript
fn extract_cds_sequence(
    transcript: &Transcript,
    genome: &Genome,
    chromosome: &str,
    strand: Strand,
) -> Result<Vec<u8>> {
    let mut cds_sequence = Vec::new();

    // Sort CDS regions by position
    let mut sorted_cds = transcript.cds_regions.clone();
    sorted_cds.sort_by_key(|cds| cds.start);

    for cds_region in &sorted_cds {
        // Get the genomic sequence for this CDS region
        let region_seq = genome.get_subsequence(&crate::types::GenomicInterval {
            chromosome: chromosome.to_string(),
            start: cds_region.start,
            end: cds_region.end,
            strand,
        })?;

        cds_sequence.extend_from_slice(&region_seq);
    }

    Ok(cds_sequence)
}

/// Validate a CDS sequence
fn validate_cds_sequence(
    cds_sequence: &[u8],
    genetic_code: &GeneticCode,
) -> Result<CdsValidationResult> {
    let mut result = CdsValidationResult {
        is_valid: true,
        has_start_codon: false,
        has_stop_codon: false,
        has_premature_stop: false,
        frame_preserved: true,
        protein_sequence: None,
        issues: Vec::new(),
    };

    // Check if sequence length is multiple of 3
    if cds_sequence.len() % 3 != 0 {
        result.is_valid = false;
        result.frame_preserved = false;
        result.issues.push(format!(
            "CDS length ({} bp) is not a multiple of 3",
            cds_sequence.len()
        ));
    }

    if cds_sequence.len() < 3 {
        result.is_valid = false;
        result.issues.push("CDS too short (< 3 bp)".to_string());
        return Ok(result);
    }

    // Convert to string for codon analysis
    let cds_string = String::from_utf8_lossy(cds_sequence).to_uppercase();

    // Check start codon
    let start_codon = &cds_string[0..3];
    result.has_start_codon = genetic_code.is_start_codon(start_codon);
    if !result.has_start_codon {
        result
            .issues
            .push(format!("Invalid start codon: {}", start_codon));
    }

    // Translate the sequence
    let mut protein = String::new();
    let mut _has_internal_stop = false;

    for i in (0..cds_string.len()).step_by(3) {
        if i + 3 > cds_string.len() {
            break; // Incomplete codon at the end
        }

        let codon = &cds_string[i..i + 3];

        if let Some(amino_acid) = genetic_code.translate_codon(codon) {
            if amino_acid == '*' {
                // Stop codon
                if i + 3 == cds_string.len() {
                    // Stop codon at the end - this is correct
                    result.has_stop_codon = true;
                } else {
                    // Premature stop codon
                    _has_internal_stop = true;
                    result.has_premature_stop = true;
                    result.issues.push(format!(
                        "Premature stop codon {} at position {}",
                        codon,
                        i + 1
                    ));
                }
            }
            protein.push(amino_acid);
        } else {
            result.is_valid = false;
            result.issues.push(format!("Invalid codon: {}", codon));
        }
    }

    // Check if we have a stop codon at the end
    if !result.has_stop_codon && cds_string.len() >= 3 {
        let last_codon = &cds_string[cds_string.len() - 3..];
        if !genetic_code.is_stop_codon(last_codon) {
            result
                .issues
                .push("Missing stop codon at end of CDS".to_string());
        }
    }

    result.protein_sequence = Some(protein);

    // Overall validity
    result.is_valid = result.has_start_codon
        && result.has_stop_codon
        && !result.has_premature_stop
        && result.frame_preserved
        && result
            .issues
            .iter()
            .all(|issue| issue.contains("non-coding") || issue.starts_with("Missing stop codon"));

    Ok(result)
}

/// Validate CDS consistency after exon boundary changes
pub fn validate_cds_after_refinement(
    original_transcript: &Transcript,
    refined_transcript: &Transcript,
    genome: &Genome,
    chromosome: &str,
    strand: Strand,
    genetic_code: &GeneticCode,
) -> Result<bool> {
    // If no CDS regions, no validation needed
    if original_transcript.cds_regions.is_empty() && refined_transcript.cds_regions.is_empty() {
        return Ok(true);
    }

    // Validate the refined CDS
    let validation_result =
        validate_cds(refined_transcript, genome, chromosome, strand, genetic_code)?;

    if !validation_result.is_valid {
        debug!(
            "CDS validation failed for transcript {}: {:?}",
            refined_transcript.id, validation_result.issues
        );
        return Ok(false);
    }

    // Additional check: ensure CDS regions are still within exons
    for cds_region in &refined_transcript.cds_regions {
        let mut cds_covered = false;
        for exon in &refined_transcript.exons {
            if cds_region.start >= exon.start && cds_region.end <= exon.end {
                cds_covered = true;
                break;
            }
        }
        if !cds_covered {
            debug!(
                "CDS region {}-{} not covered by any exon in transcript {}",
                cds_region.start, cds_region.end, refined_transcript.id
            );
            return Ok(false);
        }
    }

    Ok(true)
}

/// Update CDS regions after exon boundary changes
pub fn update_cds_after_exon_changes(
    transcript: &mut Transcript,
    _genome: &Genome,
    _chromosome: &str,
    _strand: Strand,
) -> Result<()> {
    // For now, we'll keep CDS regions unchanged if they're still within exon boundaries
    // A more sophisticated implementation would adjust CDS boundaries to maintain frame

    let mut valid_cds = Vec::new();

    for cds_region in &transcript.cds_regions {
        let mut cds_valid = false;

        // Check if CDS region is still within an exon
        for exon in &transcript.exons {
            if cds_region.start >= exon.start && cds_region.end <= exon.end {
                cds_valid = true;
                break;
            }
        }

        if cds_valid {
            valid_cds.push(cds_region.clone());
        } else {
            warn!(
                "Removing CDS region {}-{} that is no longer within exon boundaries",
                cds_region.start, cds_region.end
            );
        }
    }

    transcript.cds_regions = valid_cds;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_genetic_code() {
        let code = GeneticCode::standard();

        assert_eq!(code.translate_codon("ATG"), Some('M'));
        assert_eq!(code.translate_codon("TAA"), Some('*'));
        assert_eq!(code.translate_codon("TTT"), Some('F'));

        assert!(code.is_start_codon("ATG"));
        assert!(code.is_stop_codon("TAA"));
        assert!(code.is_stop_codon("TAG"));
        assert!(code.is_stop_codon("TGA"));
    }

    #[test]
    fn test_cds_validation() {
        let code = GeneticCode::standard();

        // Valid CDS: ATG (start) + TTT (F) + TAA (stop)
        let valid_cds = b"ATGTTTTAA";
        let result = validate_cds_sequence(valid_cds, &code).unwrap();
        assert!(result.is_valid);
        assert!(result.has_start_codon);
        assert!(result.has_stop_codon);
        assert!(!result.has_premature_stop);

        // Invalid CDS: wrong start codon
        let invalid_cds = b"TTGTTTTAA";
        let result = validate_cds_sequence(invalid_cds, &code).unwrap();
        assert!(!result.is_valid);
        assert!(!result.has_start_codon);
    }
}
