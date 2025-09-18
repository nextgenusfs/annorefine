//! BAM file reading and RNA-seq alignment processing

use crate::types::{
    AnnoRefineError, GenomicInterval, Result, RnaSeqAlignment, SpliceJunction, Strand,
};
use log::debug;
use rust_htslib::{bam, bam::Read};
use std::collections::HashMap;
use std::path::Path;

/// BAM reader for processing RNA-seq alignments
pub struct BamReader {
    reader: bam::IndexedReader,
    header: bam::HeaderView,
}

impl BamReader {
    /// Open a BAM file for reading
    pub fn new<P: AsRef<Path>>(bam_path: P) -> Result<Self> {
        let bam_path = bam_path.as_ref();
        debug!("Opening BAM file: {}", bam_path.display());

        let reader = bam::IndexedReader::from_path(bam_path).map_err(|e| {
            AnnoRefineError::BamParse(format!(
                "Failed to open BAM file {}: {}",
                bam_path.display(),
                e
            ))
        })?;

        let header = reader.header().clone();

        debug!(
            "Successfully opened BAM file with {} reference sequences",
            header.target_count()
        );

        Ok(BamReader { reader, header })
    }

    /// Get alignments in a genomic region
    pub fn get_alignments_in_region(
        &mut self,
        region: &GenomicInterval,
    ) -> Result<Vec<RnaSeqAlignment>> {
        let tid = self.get_target_id(&region.chromosome)?;

        // Convert to 0-based coordinates for BAM
        let start = region.start.saturating_sub(1);
        let end = region.end;

        self.reader.fetch((tid, start, end)).map_err(|e| {
            AnnoRefineError::BamParse(format!(
                "Failed to fetch region {}:{}-{}: {}",
                region.chromosome, region.start, region.end, e
            ))
        })?;

        let mut alignments = Vec::new();

        for result in self.reader.records() {
            let record = result.map_err(|e| {
                AnnoRefineError::BamParse(format!("Failed to read BAM record: {}", e))
            })?;

            // Skip unmapped reads
            if record.is_unmapped() {
                continue;
            }

            // Skip secondary alignments and duplicates
            if record.is_secondary() || record.is_duplicate() {
                continue;
            }

            let alignment = parse_bam_record(&record, &self.header)?;
            alignments.push(alignment);
        }

        debug!(
            "Found {} alignments in region {}:{}-{}",
            alignments.len(),
            region.chromosome,
            region.start,
            region.end
        );

        Ok(alignments)
    }

    /// Get splice junctions from alignments in a region
    pub fn get_splice_junctions_in_region(
        &mut self,
        region: &GenomicInterval,
    ) -> Result<Vec<SpliceJunction>> {
        let alignments = self.get_alignments_in_region(region)?;
        let mut junction_counts: HashMap<(u64, u64, Strand), u32> = HashMap::new();

        for alignment in alignments {
            for junction in alignment.splice_junctions {
                let key = (junction.donor_pos, junction.acceptor_pos, junction.strand);
                *junction_counts.entry(key).or_insert(0) += 1;
            }
        }

        let mut junctions = Vec::new();
        for ((donor_pos, acceptor_pos, strand), count) in junction_counts {
            junctions.push(SpliceJunction {
                chromosome: region.chromosome.clone(),
                donor_pos,
                acceptor_pos,
                strand, // Use the actual strand from the alignment
                support_count: count,
            });
        }

        // Sort by position
        junctions.sort_by_key(|j| j.donor_pos);

        debug!(
            "Found {} unique splice junctions in region",
            junctions.len()
        );
        Ok(junctions)
    }

    /// Get coverage depth at each position in a region
    pub fn get_coverage_in_region(&mut self, region: &GenomicInterval) -> Result<Vec<u32>> {
        let alignments = self.get_alignments_in_region(region)?;
        let region_length = (region.end - region.start + 1) as usize;
        let mut coverage = vec![0u32; region_length];

        for alignment in alignments {
            let align_start = alignment.start.max(region.start);
            let align_end = alignment.end.min(region.end);

            for pos in align_start..=align_end {
                let idx = (pos - region.start) as usize;
                if idx < coverage.len() {
                    coverage[idx] = coverage[idx].saturating_add(1);
                }
            }
        }

        Ok(coverage)
    }

    fn get_target_id(&self, chromosome: &str) -> Result<u32> {
        for (tid, target) in self.header.target_names().iter().enumerate() {
            if *target == chromosome.as_bytes() {
                return Ok(tid as u32);
            }
        }

        Err(AnnoRefineError::BamParse(format!(
            "Chromosome '{}' not found in BAM header",
            chromosome
        )))
    }
}

/// Parse a BAM record into an RnaSeqAlignment
fn parse_bam_record(record: &bam::Record, header: &bam::HeaderView) -> Result<RnaSeqAlignment> {
    let read_name = String::from_utf8_lossy(record.qname()).to_string();

    let tid = record.tid();
    let chromosome = if tid >= 0 {
        String::from_utf8_lossy(header.tid2name(tid as u32)).to_string()
    } else {
        return Err(AnnoRefineError::BamParse("Invalid target ID".to_string()));
    };

    let start = record.pos() as u64 + 1; // Convert to 1-based

    // Get CIGAR string - try cached first, then get directly
    let (end, cigar_string, splice_junctions) = if let Some(cigar_view) = record.cigar_cached() {
        let end = cigar_view.end_pos() as u64;
        let cigar_string = format!("{}", cigar_view);
        let splice_junctions = extract_splice_junctions_from_cigar_view(
            cigar_view,
            start,
            &chromosome,
            if record.is_reverse() {
                Strand::Reverse
            } else {
                Strand::Forward
            },
        );
        (end, cigar_string, splice_junctions)
    } else {
        // Get CIGAR directly if not cached
        let cigar_owned = record.cigar();
        let end = cigar_owned.end_pos() as u64;
        let cigar_string = format!("{}", cigar_owned);
        let splice_junctions = extract_splice_junctions_from_cigar_owned(
            &cigar_owned,
            start,
            &chromosome,
            if record.is_reverse() {
                Strand::Reverse
            } else {
                Strand::Forward
            },
        );
        (end, cigar_string, splice_junctions)
    };

    let strand = if record.is_reverse() {
        Strand::Reverse
    } else {
        Strand::Forward
    };

    let mapping_quality = record.mapq();

    Ok(RnaSeqAlignment {
        read_name,
        chromosome,
        start,
        end,
        strand,
        cigar: cigar_string,
        mapping_quality,
        splice_junctions,
    })
}

/// Extract splice junctions from cached CIGAR string
fn extract_splice_junctions_from_cigar_view(
    cigar: &bam::record::CigarStringView,
    alignment_start: u64,
    chromosome: &str,
    strand: Strand,
) -> Vec<SpliceJunction> {
    let mut junctions = Vec::new();
    let mut current_pos = alignment_start;

    for operation in cigar.iter() {
        match operation {
            bam::record::Cigar::Match(len)
            | bam::record::Cigar::Equal(len)
            | bam::record::Cigar::Diff(len) => {
                current_pos += *len as u64;
            }
            bam::record::Cigar::RefSkip(len) => {
                // This is an intron (splice junction)
                let donor_pos = current_pos - 1; // Last position of upstream exon
                current_pos += *len as u64;
                let acceptor_pos = current_pos; // First position of downstream exon

                junctions.push(SpliceJunction {
                    chromosome: chromosome.to_string(),
                    donor_pos,
                    acceptor_pos,
                    strand,
                    support_count: 1, // Will be aggregated later
                });
            }
            bam::record::Cigar::Del(len) => {
                current_pos += *len as u64;
            }
            bam::record::Cigar::Ins(_)
            | bam::record::Cigar::SoftClip(_)
            | bam::record::Cigar::HardClip(_)
            | bam::record::Cigar::Pad(_) => {
                // These don't advance the reference position
            }
        }
    }

    junctions
}

/// Extract splice junctions from owned CIGAR string
fn extract_splice_junctions_from_cigar_owned(
    cigar: &bam::record::CigarString,
    alignment_start: u64,
    chromosome: &str,
    strand: Strand,
) -> Vec<SpliceJunction> {
    let mut junctions = Vec::new();
    let mut current_pos = alignment_start;

    for operation in cigar.iter() {
        match operation {
            bam::record::Cigar::Match(len)
            | bam::record::Cigar::Equal(len)
            | bam::record::Cigar::Diff(len) => {
                current_pos += *len as u64;
            }
            bam::record::Cigar::RefSkip(len) => {
                // This is an intron (splice junction)
                let donor_pos = current_pos - 1; // Last position of upstream exon
                current_pos += *len as u64;
                let acceptor_pos = current_pos; // First position of downstream exon

                junctions.push(SpliceJunction {
                    chromosome: chromosome.to_string(),
                    donor_pos,
                    acceptor_pos,
                    strand,
                    support_count: 1, // Will be aggregated later
                });
            }
            bam::record::Cigar::Del(len) => {
                current_pos += *len as u64;
            }
            bam::record::Cigar::Ins(_)
            | bam::record::Cigar::SoftClip(_)
            | bam::record::Cigar::HardClip(_)
            | bam::record::Cigar::Pad(_) => {
                // These don't advance the reference position
            }
        }
    }

    junctions
}

/// Calculate alignment statistics for a BAM file
pub fn get_alignment_stats<P: AsRef<Path>>(bam_path: P) -> Result<AlignmentStats> {
    let bam_path = bam_path.as_ref();
    let mut reader = bam::Reader::from_path(bam_path).map_err(|e| {
        AnnoRefineError::BamParse(format!(
            "Failed to open BAM file {}: {}",
            bam_path.display(),
            e
        ))
    })?;

    let mut total_reads = 0;
    let mut mapped_reads = 0;
    let mut spliced_reads = 0;
    let mut total_splice_junctions = 0;

    for result in reader.records() {
        let record = result
            .map_err(|e| AnnoRefineError::BamParse(format!("Failed to read BAM record: {}", e)))?;

        total_reads += 1;

        if !record.is_unmapped() {
            mapped_reads += 1;

            // Check for splice junctions (RefSkip operations in CIGAR)
            let cigar = match record.cigar_cached() {
                Some(cigar) => cigar,
                None => {
                    // If CIGAR is not cached, get it directly
                    let cigar_owned = record.cigar();
                    // We need to work with the owned CIGAR directly
                    let mut has_splice = false;
                    for operation in cigar_owned.iter() {
                        if let bam::record::Cigar::RefSkip(_) = operation {
                            has_splice = true;
                            total_splice_junctions += 1;
                        }
                    }
                    if has_splice {
                        spliced_reads += 1;
                    }
                    continue; // Skip the rest of the loop since we handled it here
                }
            };

            let mut has_splice = false;
            for operation in cigar.iter() {
                if let bam::record::Cigar::RefSkip(_) = operation {
                    has_splice = true;
                    total_splice_junctions += 1;
                }
            }
            if has_splice {
                spliced_reads += 1;
            }
        }
    }

    Ok(AlignmentStats {
        total_reads,
        mapped_reads,
        spliced_reads,
        total_splice_junctions,
    })
}

/// Statistics about BAM alignments
#[derive(Debug)]
pub struct AlignmentStats {
    pub total_reads: u64,
    pub mapped_reads: u64,
    pub spliced_reads: u64,
    pub total_splice_junctions: u64,
}

impl std::fmt::Display for AlignmentStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mapping_rate = if self.total_reads > 0 {
            (self.mapped_reads as f64 / self.total_reads as f64) * 100.0
        } else {
            0.0
        };

        let splice_rate = if self.mapped_reads > 0 {
            (self.spliced_reads as f64 / self.mapped_reads as f64) * 100.0
        } else {
            0.0
        };

        write!(
            f,
            "Total reads: {}, Mapped: {} ({:.2}%), Spliced: {} ({:.2}%), Splice junctions: {}",
            self.total_reads,
            self.mapped_reads,
            mapping_rate,
            self.spliced_reads,
            splice_rate,
            self.total_splice_junctions
        )
    }
}
