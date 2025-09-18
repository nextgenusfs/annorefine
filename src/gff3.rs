//! GFF3 file parsing and gene model construction

use crate::types::{
    Gff3Feature, GeneModel, Transcript, Exon, CdsRegion, FeatureType, Strand,
    AnnoRefineError, Result
};
use std::collections::HashMap;
use std::path::Path;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::str::FromStr;
use log::{info, debug};

/// Parse a GFF3 file and return gene models
pub fn parse_gff3_file<P: AsRef<Path>>(path: P) -> Result<Vec<GeneModel>> {
    let path = path.as_ref();
    info!("Parsing GFF3 file: {}", path.display());
    
    let file = File::open(path)
        .map_err(|e| AnnoRefineError::Gff3Parse(
            format!("Failed to open GFF3 file {}: {}", path.display(), e)
        ))?;
    
    let reader = BufReader::new(file);
    let mut features = Vec::new();
    let mut line_number = 0;
    
    for line in reader.lines() {
        line_number += 1;
        let line = line?;
        
        // Skip comments and empty lines
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        
        match parse_gff3_line(&line) {
            Ok(feature) => features.push(feature),
            Err(e) => {
                debug!("Error parsing line {}: {}", line_number, e);
                continue;
            }
        }
    }
    
    info!("Parsed {} features from GFF3 file", features.len());
    
    // Build gene models from features
    build_gene_models(features)
}

/// Parse a single GFF3 line into a feature
fn parse_gff3_line(line: &str) -> Result<Gff3Feature> {
    let fields: Vec<&str> = line.split('\t').collect();
    
    if fields.len() != 9 {
        return Err(AnnoRefineError::Gff3Parse(
            format!("Expected 9 fields, found {}", fields.len())
        ));
    }
    
    let seqid = fields[0].to_string();
    let source = fields[1].to_string();
    let feature_type = FeatureType::from_str(fields[2])?;
    
    let start = fields[3].parse::<u64>()
        .map_err(|_| AnnoRefineError::Gff3Parse(
            format!("Invalid start position: {}", fields[3])
        ))?;
    
    let end = fields[4].parse::<u64>()
        .map_err(|_| AnnoRefineError::Gff3Parse(
            format!("Invalid end position: {}", fields[4])
        ))?;
    
    if start > end {
        return Err(AnnoRefineError::Gff3Parse(
            format!("Start position {} is greater than end position {}", start, end)
        ));
    }
    
    let score = if fields[5] == "." {
        None
    } else {
        Some(fields[5].parse::<f64>()
            .map_err(|_| AnnoRefineError::Gff3Parse(
                format!("Invalid score: {}", fields[5])
            ))?)
    };
    
    let strand = Strand::from_str(fields[6])?;
    
    let phase = if fields[7] == "." {
        None
    } else {
        Some(fields[7].parse::<u8>()
            .map_err(|_| AnnoRefineError::Gff3Parse(
                format!("Invalid phase: {}", fields[7])
            ))?)
    };
    
    let attributes = parse_attributes(fields[8])?;
    
    Ok(Gff3Feature {
        seqid,
        source,
        feature_type,
        start,
        end,
        score,
        strand,
        phase,
        attributes,
    })
}

/// Parse GFF3 attributes field
fn parse_attributes(attr_str: &str) -> Result<HashMap<String, Vec<String>>> {
    let mut attributes = HashMap::new();
    
    if attr_str.trim().is_empty() || attr_str == "." {
        return Ok(attributes);
    }
    
    for pair in attr_str.split(';') {
        let pair = pair.trim();
        if pair.is_empty() {
            continue;
        }
        
        let parts: Vec<&str> = pair.splitn(2, '=').collect();
        if parts.len() != 2 {
            return Err(AnnoRefineError::Gff3Parse(
                format!("Invalid attribute format: {}", pair)
            ));
        }
        
        let key = url_decode(parts[0])?;
        let values: Vec<String> = parts[1]
            .split(',')
            .map(|v| url_decode(v))
            .collect::<Result<Vec<_>>>()?;
        
        attributes.insert(key, values);
    }
    
    Ok(attributes)
}

/// Simple URL decoding for GFF3 attributes
fn url_decode(s: &str) -> Result<String> {
    let mut result = String::new();
    let mut chars = s.chars();
    
    while let Some(c) = chars.next() {
        if c == '%' {
            let hex: String = chars.by_ref().take(2).collect();
            if hex.len() != 2 {
                return Err(AnnoRefineError::Gff3Parse(
                    format!("Invalid URL encoding: %{}", hex)
                ));
            }
            
            let byte = u8::from_str_radix(&hex, 16)
                .map_err(|_| AnnoRefineError::Gff3Parse(
                    format!("Invalid hex in URL encoding: {}", hex)
                ))?;
            
            result.push(byte as char);
        } else {
            result.push(c);
        }
    }
    
    Ok(result)
}

/// Build gene models from parsed features
fn build_gene_models(features: Vec<Gff3Feature>) -> Result<Vec<GeneModel>> {
    let mut genes = HashMap::new();
    let mut transcripts = HashMap::new();
    let mut exons = Vec::new();
    let mut cds_regions = Vec::new();
    
    // Separate features by type
    for feature in features {
        match feature.feature_type {
            FeatureType::Gene => {
                if let Some(id) = feature.get_id() {
                    genes.insert(id.to_string(), feature);
                }
            }
            FeatureType::Mrna => {
                if let Some(id) = feature.get_id() {
                    transcripts.insert(id.to_string(), feature);
                }
            }
            FeatureType::Exon => exons.push(feature),
            FeatureType::Cds => cds_regions.push(feature),
            _ => {} // Ignore other feature types for now
        }
    }
    
    debug!("Found {} genes, {} transcripts, {} exons, {} CDS regions", 
           genes.len(), transcripts.len(), exons.len(), cds_regions.len());
    
    // Build gene models
    let mut gene_models = Vec::new();
    
    for (gene_id, gene_feature) in genes {
        let mut gene_transcripts = Vec::new();
        
        // Find transcripts for this gene
        for (transcript_id, transcript_feature) in &transcripts {
            if let Some(parent) = transcript_feature.get_parent() {
                if parent == gene_id {
                    // Build transcript
                    let transcript = build_transcript(
                        transcript_id,
                        transcript_feature,
                        &exons,
                        &cds_regions,
                    )?;
                    gene_transcripts.push(transcript);
                }
            }
        }
        
        if gene_transcripts.is_empty() {
            debug!("Gene {} has no transcripts", gene_id);
            continue;
        }
        
        // Calculate original CDS lengths for each transcript
        let mut original_cds_lengths = HashMap::new();
        for transcript in &gene_transcripts {
            let cds_length: u64 = transcript.cds_regions.iter()
                .map(|cds| cds.end - cds.start)
                .sum();
            original_cds_lengths.insert(transcript.id.clone(), cds_length);
        }

        let gene_model = GeneModel {
            id: gene_id.clone(),
            name: gene_feature.get_name().map(|s| s.to_string()),
            chromosome: gene_feature.seqid.clone(),
            start: gene_feature.start,
            end: gene_feature.end,
            strand: gene_feature.strand,
            transcripts: gene_transcripts,
            original_attributes: gene_feature.attributes.clone(),
            has_structural_changes: false, // Will be set during refinement
            original_cds_lengths,
        };
        
        gene_models.push(gene_model);
    }
    
    info!("Built {} gene models", gene_models.len());
    Ok(gene_models)
}

/// Build a transcript from features
fn build_transcript(
    transcript_id: &str,
    transcript_feature: &Gff3Feature,
    exons: &[Gff3Feature],
    cds_regions: &[Gff3Feature],
) -> Result<Transcript> {
    // Find exons for this transcript
    let mut transcript_exons = Vec::new();
    for exon in exons {
        if let Some(parent) = exon.get_parent() {
            if parent == transcript_id {
                transcript_exons.push(Exon {
                    id: exon.get_id().map(|s| s.to_string()),
                    start: exon.start,
                    end: exon.end,
                });
            }
        }
    }
    
    // Sort exons by position
    transcript_exons.sort_by_key(|e| e.start);
    
    // Find CDS regions for this transcript
    let mut transcript_cds = Vec::new();
    for cds in cds_regions {
        if let Some(parent) = cds.get_parent() {
            if parent == transcript_id {
                transcript_cds.push(CdsRegion {
                    id: cds.get_id().map(|s| s.to_string()),
                    start: cds.start,
                    end: cds.end,
                    phase: cds.phase,
                });
            }
        }
    }
    
    // Sort CDS regions by position
    transcript_cds.sort_by_key(|c| c.start);
    
    Ok(Transcript {
        id: transcript_id.to_string(),
        name: transcript_feature.get_name().map(|s| s.to_string()),
        start: transcript_feature.start,
        end: transcript_feature.end,
        exons: transcript_exons,
        cds_regions: transcript_cds,
        five_prime_utr: None,  // Will be computed later
        three_prime_utr: None, // Will be computed later
        original_attributes: transcript_feature.attributes.clone(),
    })
}
