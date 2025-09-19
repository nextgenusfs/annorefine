//! Output generation for refined gene models

use crate::types::{
    AnnoRefineError, CdsRegion, Exon, FeatureType, GeneModel, Result, Strand, Transcript,
};
use log::{debug, info};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// GFF3 writer for outputting refined gene models
pub struct Gff3Writer {
    writer: BufWriter<File>,
}

impl Gff3Writer {
    /// Create a new GFF3 writer
    pub fn new<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        info!("Creating GFF3 output file: {}", path.display());

        let file = File::create(path).map_err(AnnoRefineError::Io)?;

        let mut writer = BufWriter::new(file);

        // Write GFF3 header
        writeln!(writer, "##gff-version 3")?;
        writeln!(writer, "##generated-by annorefine v{}", env!("CARGO_PKG_VERSION"))?;

        Ok(Gff3Writer { writer })
    }

    /// Write gene models to the GFF3 file
    pub fn write_gene_models(&mut self, gene_models: &[GeneModel]) -> Result<()> {
        info!("Writing {} gene models to GFF3", gene_models.len());

        // Sort gene models by chromosome and then by start position for proper GFF3 format
        let mut sorted_genes = gene_models.to_vec();
        sorted_genes.sort_by(|a, b| {
            // First sort by chromosome name (with smart numerical sorting)
            match compare_chromosome_names(&a.chromosome, &b.chromosome) {
                std::cmp::Ordering::Equal => {
                    // If chromosomes are the same, sort by start position
                    a.start.cmp(&b.start)
                }
                other => other,
            }
        });

        debug!(
            "Sorted {} gene models by chromosome and position",
            sorted_genes.len()
        );

        for gene_model in &sorted_genes {
            self.write_gene_model(gene_model)?;
        }

        self.writer.flush()?;
        info!("Successfully wrote gene models to GFF3");

        Ok(())
    }

    /// Write a single gene model
    fn write_gene_model(&mut self, gene_model: &GeneModel) -> Result<()> {
        debug!("Writing gene: {}", gene_model.id);

        // Write gene feature
        self.write_gene_feature(gene_model)?;

        // Write transcripts
        for transcript in &gene_model.transcripts {
            self.write_transcript(transcript, gene_model)?;
        }

        Ok(())
    }

    /// Write gene feature
    fn write_gene_feature(&mut self, gene_model: &GeneModel) -> Result<()> {
        // Start with original attributes and update key ones
        let mut attributes = gene_model.original_attributes.clone();

        // Ensure ID and Name are current (in case they were refined)
        attributes.insert("ID".to_string(), vec![gene_model.id.clone()]);
        if let Some(name) = &gene_model.name {
            attributes.insert("Name".to_string(), vec![name.clone()]);
        }

        // Add warning note if gene has structural changes
        if gene_model.has_structural_changes {
            // Calculate CDS and protein length changes for all transcripts
            let mut change_details = Vec::new();

            for transcript in &gene_model.transcripts {
                // Calculate current CDS length
                let current_cds_length: u64 = transcript
                    .cds_regions
                    .iter()
                    .map(|cds| cds.end - cds.start)
                    .sum();

                // Get original CDS length
                let original_cds_length = gene_model
                    .original_cds_lengths
                    .get(&transcript.id)
                    .copied()
                    .unwrap_or(0);

                if current_cds_length > 0 || original_cds_length > 0 {
                    let original_protein_length = original_cds_length / 3;
                    let current_protein_length = current_cds_length / 3;

                    if original_cds_length != current_cds_length {
                        change_details.push(format!(
                            "{}:{}bp/{}aa→{}bp/{}aa",
                            transcript.id,
                            original_cds_length,
                            original_protein_length,
                            current_cds_length,
                            current_protein_length
                        ));
                    } else {
                        // Structure changed but CDS length stayed the same
                        change_details.push(format!(
                            "{}:{}bp/{}aa(structure_modified)",
                            transcript.id, current_cds_length, current_protein_length
                        ));
                    }
                }
            }

            let warning_note = if !change_details.is_empty() {
                format!("CAUTION: Gene structure modified by AnnoRefine - CDS and protein sequence may have changed ({}), previous functional annotations may be inaccurate",
                    change_details.join(", "))
            } else {
                "CAUTION: Gene structure modified by AnnoRefine - CDS and protein sequence may have changed, previous functional annotations may be inaccurate".to_string()
            };

            // Add to existing Note attribute or create new one
            match attributes.get_mut("Note") {
                Some(notes) => {
                    notes.push(warning_note);
                }
                None => {
                    attributes.insert("Note".to_string(), vec![warning_note]);
                }
            }

            // Also add a specific refinement flag
            attributes.insert(
                "annorefine_structural_changes".to_string(),
                vec!["true".to_string()],
            );
        }

        self.write_feature(
            &gene_model.chromosome,
            "annorefine",
            &FeatureType::Gene,
            gene_model.start,
            gene_model.end,
            None,
            gene_model.strand,
            None,
            &attributes,
        )?;

        Ok(())
    }

    /// Write transcript and its features
    fn write_transcript(&mut self, transcript: &Transcript, gene_model: &GeneModel) -> Result<()> {
        debug!("Writing transcript: {}", transcript.id);

        // Write mRNA feature - start with original attributes
        let mut mrna_attributes = transcript.original_attributes.clone();

        // Ensure key attributes are current
        mrna_attributes.insert("ID".to_string(), vec![transcript.id.clone()]);
        mrna_attributes.insert("Parent".to_string(), vec![gene_model.id.clone()]);

        if let Some(name) = &transcript.name {
            mrna_attributes.insert("Name".to_string(), vec![name.clone()]);
        }

        // Add warning note if gene has structural changes
        if gene_model.has_structural_changes {
            // Calculate CDS change for this specific transcript
            let current_cds_length: u64 = transcript
                .cds_regions
                .iter()
                .map(|cds| cds.end - cds.start)
                .sum();

            let original_cds_length = gene_model
                .original_cds_lengths
                .get(&transcript.id)
                .copied()
                .unwrap_or(0);

            let warning_note = if original_cds_length != current_cds_length {
                let original_protein_length = original_cds_length / 3;
                let current_protein_length = current_cds_length / 3;
                format!("CAUTION: Transcript structure modified by AnnoRefine - CDS changed from {}bp/{}aa to {}bp/{}aa, protein sequence may have changed",
                    original_cds_length, original_protein_length,
                    current_cds_length, current_protein_length)
            } else if current_cds_length > 0 {
                let protein_length = current_cds_length / 3;
                format!("CAUTION: Transcript structure modified by AnnoRefine - CDS length unchanged ({}bp/{}aa) but structure modified, protein sequence may have changed",
                    current_cds_length, protein_length)
            } else {
                "CAUTION: Transcript structure modified by AnnoRefine - non-coding transcript structure changed".to_string()
            };

            // Add to existing Note attribute or create new one
            match mrna_attributes.get_mut("Note") {
                Some(notes) => {
                    notes.push(warning_note);
                }
                None => {
                    mrna_attributes.insert("Note".to_string(), vec![warning_note]);
                }
            }
        }

        self.write_feature(
            &gene_model.chromosome,
            "annorefine",
            &FeatureType::Mrna,
            transcript.start,
            transcript.end,
            None,
            gene_model.strand,
            None,
            &mrna_attributes,
        )?;

        // Write exons
        for (i, exon) in transcript.exons.iter().enumerate() {
            self.write_exon(exon, transcript, gene_model, i)?;
        }

        // Write CDS regions
        for (i, cds) in transcript.cds_regions.iter().enumerate() {
            self.write_cds(cds, transcript, gene_model, i)?;
        }

        // Write UTRs if present
        self.write_utrs(transcript, gene_model)?;

        Ok(())
    }

    /// Write exon feature
    fn write_exon(
        &mut self,
        exon: &Exon,
        transcript: &Transcript,
        gene_model: &GeneModel,
        index: usize,
    ) -> Result<()> {
        let mut attributes = HashMap::new();

        let exon_id = if let Some(id) = &exon.id {
            id.clone()
        } else {
            format!("{}.exon{}", transcript.id, index + 1)
        };

        attributes.insert("ID".to_string(), vec![exon_id]);
        attributes.insert("Parent".to_string(), vec![transcript.id.clone()]);

        self.write_feature(
            &gene_model.chromosome,
            "annorefine",
            &FeatureType::Exon,
            exon.start,
            exon.end,
            None,
            gene_model.strand,
            None,
            &attributes,
        )?;

        Ok(())
    }

    /// Write CDS feature
    fn write_cds(
        &mut self,
        cds: &CdsRegion,
        transcript: &Transcript,
        gene_model: &GeneModel,
        index: usize,
    ) -> Result<()> {
        let mut attributes = HashMap::new();

        let cds_id = if let Some(id) = &cds.id {
            id.clone()
        } else {
            format!("{}.cds{}", transcript.id, index + 1)
        };

        attributes.insert("ID".to_string(), vec![cds_id]);
        attributes.insert("Parent".to_string(), vec![transcript.id.clone()]);

        self.write_feature(
            &gene_model.chromosome,
            "annorefine",
            &FeatureType::Cds,
            cds.start,
            cds.end,
            None,
            gene_model.strand,
            cds.phase,
            &attributes,
        )?;

        Ok(())
    }

    /// Write UTR features
    fn write_utrs(&mut self, transcript: &Transcript, _gene_model: &GeneModel) -> Result<()> {
        // Write 5' UTR
        if let Some(five_prime_utrs) = &transcript.five_prime_utr {
            for (i, utr_region) in five_prime_utrs.iter().enumerate() {
                let mut attributes = HashMap::new();
                attributes.insert(
                    "ID".to_string(),
                    vec![format!("{}.5utr{}", transcript.id, i + 1)],
                );
                attributes.insert("Parent".to_string(), vec![transcript.id.clone()]);

                self.write_feature(
                    &utr_region.chromosome,
                    "annorefine",
                    &FeatureType::FivePrimeUtr,
                    utr_region.start,
                    utr_region.end,
                    None,
                    utr_region.strand,
                    None,
                    &attributes,
                )?;
            }
        }

        // Write 3' UTR
        if let Some(three_prime_utrs) = &transcript.three_prime_utr {
            for (i, utr_region) in three_prime_utrs.iter().enumerate() {
                let mut attributes = HashMap::new();
                attributes.insert(
                    "ID".to_string(),
                    vec![format!("{}.3utr{}", transcript.id, i + 1)],
                );
                attributes.insert("Parent".to_string(), vec![transcript.id.clone()]);

                self.write_feature(
                    &utr_region.chromosome,
                    "annorefine",
                    &FeatureType::ThreePrimeUtr,
                    utr_region.start,
                    utr_region.end,
                    None,
                    utr_region.strand,
                    None,
                    &attributes,
                )?;
            }
        }

        Ok(())
    }

    /// Write a generic GFF3 feature
    fn write_feature(
        &mut self,
        seqid: &str,
        source: &str,
        feature_type: &FeatureType,
        start: u64,
        end: u64,
        score: Option<f64>,
        strand: Strand,
        phase: Option<u8>,
        attributes: &HashMap<String, Vec<String>>,
    ) -> Result<()> {
        let score_str = match score {
            Some(s) => format!("{:.3}", s),
            None => ".".to_string(),
        };

        let phase_str = match phase {
            Some(p) => p.to_string(),
            None => ".".to_string(),
        };

        let attributes_str = format_attributes(attributes);

        writeln!(
            self.writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            seqid, source, feature_type, start, end, score_str, strand, phase_str, attributes_str
        )?;

        Ok(())
    }
}

/// Format GFF3 attributes
fn format_attributes(attributes: &HashMap<String, Vec<String>>) -> String {
    if attributes.is_empty() {
        return ".".to_string();
    }

    let mut attr_strings = Vec::new();

    // Ensure ID comes first if present
    if let Some(id_values) = attributes.get("ID") {
        attr_strings.push(format!("ID={}", id_values.join(",")));
    }

    // Add other attributes in alphabetical order
    let mut other_keys: Vec<_> = attributes.keys().filter(|k| *k != "ID").collect();
    other_keys.sort();

    for key in other_keys {
        let values = &attributes[key];
        attr_strings.push(format!("{}={}", key, values.join(",")));
    }

    attr_strings.join(";")
}

/// Validate gene models before output
pub fn validate_gene_models(gene_models: &[GeneModel]) -> Result<()> {
    info!("Validating {} gene models", gene_models.len());

    for gene_model in gene_models {
        validate_gene_model(gene_model)?;
    }

    info!("All gene models passed validation");
    Ok(())
}

/// Validate a single gene model
fn validate_gene_model(gene_model: &GeneModel) -> Result<()> {
    // Check that gene has at least one transcript
    if gene_model.transcripts.is_empty() {
        return Err(AnnoRefineError::InvalidGeneModel(format!(
            "Gene {} has no transcripts",
            gene_model.id
        )));
    }

    // Check gene boundaries
    if gene_model.start >= gene_model.end {
        return Err(AnnoRefineError::InvalidGeneModel(format!(
            "Gene {} has invalid coordinates: {}-{}",
            gene_model.id, gene_model.start, gene_model.end
        )));
    }

    // Validate each transcript
    for transcript in &gene_model.transcripts {
        validate_transcript(transcript, gene_model)?;
    }

    Ok(())
}

/// Validate a transcript
fn validate_transcript(transcript: &Transcript, gene_model: &GeneModel) -> Result<()> {
    // Check that transcript has at least one exon
    if transcript.exons.is_empty() {
        return Err(AnnoRefineError::InvalidGeneModel(format!(
            "Transcript {} has no exons",
            transcript.id
        )));
    }

    // Check transcript boundaries
    if transcript.start >= transcript.end {
        return Err(AnnoRefineError::InvalidGeneModel(format!(
            "Transcript {} has invalid coordinates: {}-{}",
            transcript.id, transcript.start, transcript.end
        )));
    }

    // Check that transcript is within gene boundaries
    if transcript.start < gene_model.start || transcript.end > gene_model.end {
        return Err(AnnoRefineError::InvalidGeneModel(format!(
            "Transcript {} extends beyond gene {} boundaries",
            transcript.id, gene_model.id
        )));
    }

    // Validate exons
    for (i, exon) in transcript.exons.iter().enumerate() {
        if exon.start >= exon.end {
            return Err(AnnoRefineError::InvalidGeneModel(format!(
                "Exon {} in transcript {} (gene {}) has invalid coordinates: {}-{}",
                i + 1,
                transcript.id,
                gene_model.id,
                exon.start,
                exon.end
            )));
        }

        if exon.start == 0 {
            return Err(AnnoRefineError::InvalidGeneModel(format!(
                "Exon {} in transcript {} (gene {}) has zero start coordinate",
                i + 1,
                transcript.id,
                gene_model.id
            )));
        }

        // Check that exons don't overlap (except for alternative splicing)
        for (j, other_exon) in transcript.exons.iter().enumerate() {
            if i != j && exons_overlap(exon, other_exon) {
                return Err(AnnoRefineError::InvalidGeneModel(
                    format!("Overlapping exons in transcript {} (gene {}): exon {} ({}-{}) and exon {} ({}-{})",
                        transcript.id, gene_model.id, i + 1, exon.start, exon.end,
                        j + 1, other_exon.start, other_exon.end)
                ));
            }
        }
    }

    Ok(())
}

/// Check if two exons overlap
fn exons_overlap(exon1: &Exon, exon2: &Exon) -> bool {
    exon1.start < exon2.end && exon2.start < exon1.end
}

/// Smart chromosome name comparison that handles numerical sorting
/// Examples: chr1 < chr2 < chr10, scaffold1 < scaffold2 < scaffold10
fn compare_chromosome_names(a: &str, b: &str) -> std::cmp::Ordering {
    // Try to extract numeric suffix for smart sorting
    let extract_parts = |s: &str| -> (String, Option<u64>) {
        // Find the last sequence of digits
        let mut last_digit_start = None;
        let chars: Vec<char> = s.chars().collect();

        for (i, c) in chars.iter().enumerate().rev() {
            if c.is_ascii_digit() {
                if last_digit_start.is_none() {
                    last_digit_start = Some(i);
                }
            } else if last_digit_start.is_some() {
                // Found the start of the prefix
                let prefix = chars[..=i].iter().collect::<String>();
                let number_str = chars[i + 1..].iter().collect::<String>();
                if let Ok(number) = number_str.parse::<u64>() {
                    return (prefix, Some(number));
                }
                break;
            }
        }

        // If we found digits at the start, handle that case
        if let Some(start) = last_digit_start {
            if start == 0 {
                // Entire string is digits
                if let Ok(number) = s.parse::<u64>() {
                    return (String::new(), Some(number));
                }
            } else {
                // Digits at the end
                let prefix = chars[..start].iter().collect::<String>();
                let number_str = chars[start..].iter().collect::<String>();
                if let Ok(number) = number_str.parse::<u64>() {
                    return (prefix, Some(number));
                }
            }
        }

        // No numeric suffix found, return the whole string
        (s.to_string(), None)
    };

    let (prefix_a, num_a) = extract_parts(a);
    let (prefix_b, num_b) = extract_parts(b);

    // First compare prefixes
    match prefix_a.cmp(&prefix_b) {
        std::cmp::Ordering::Equal => {
            // Prefixes are the same, compare numbers
            match (num_a, num_b) {
                (Some(n_a), Some(n_b)) => n_a.cmp(&n_b),
                (Some(_), None) => std::cmp::Ordering::Less, // Numbers come before non-numbers
                (None, Some(_)) => std::cmp::Ordering::Greater, // Numbers come before non-numbers
                (None, None) => std::cmp::Ordering::Equal,   // Both have no numbers
            }
        }
        other => other,
    }
}
