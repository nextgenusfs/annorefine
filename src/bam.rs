//! BAM file reading and RNA-seq alignment processing

use crate::types::{
    AnnoRefineError, GenomicInterval, Result, RnaSeqAlignment, SpliceJunction, Strand, StrandBias,
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
        self.get_alignments_in_region_with_library_type(region, crate::types::LibraryType::Auto)
    }

    /// Get alignments in a genomic region with library type filtering
    pub fn get_alignments_in_region_with_library_type(
        &mut self,
        region: &GenomicInterval,
        library_type: crate::types::LibraryType,
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

        // Apply library type filtering based on gene strand
        let filtered_alignments =
            filter_alignments_by_library_type(alignments, region.strand, library_type);

        debug!(
            "After library type filtering: {} alignments (library type: {})",
            filtered_alignments.len(),
            library_type
        );

        Ok(filtered_alignments)
    }

    /// Get splice junctions from alignments in a region
    pub fn get_splice_junctions_in_region(
        &mut self,
        region: &GenomicInterval,
    ) -> Result<Vec<SpliceJunction>> {
        self.get_splice_junctions_in_region_with_library_type(
            region,
            crate::types::LibraryType::Auto,
        )
    }

    /// Get splice junctions from alignments in a region with library type filtering
    pub fn get_splice_junctions_in_region_with_library_type(
        &mut self,
        region: &GenomicInterval,
        library_type: crate::types::LibraryType,
    ) -> Result<Vec<SpliceJunction>> {
        let alignments = self.get_alignments_in_region_with_library_type(region, library_type)?;
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
        self.get_coverage_in_region_with_library_type(region, crate::types::LibraryType::Auto)
    }

    /// Get coverage array for a genomic region with library type filtering
    pub fn get_coverage_in_region_with_library_type(
        &mut self,
        region: &GenomicInterval,
        library_type: crate::types::LibraryType,
    ) -> Result<Vec<u32>> {
        let alignments = self.get_alignments_in_region_with_library_type(region, library_type)?;
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

    /// Get the BAM header
    pub fn get_header(&self) -> &bam::HeaderView {
        &self.header
    }

    /// Get chromosome length from BAM header
    pub fn get_chromosome_length(&self, chromosome: &str) -> Result<u64> {
        let tid = self.get_target_id(chromosome)?;
        Ok(self.header.target_len(tid).unwrap_or(0) as u64)
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
pub fn parse_bam_record(record: &bam::Record, header: &bam::HeaderView) -> Result<RnaSeqAlignment> {
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

    // Extract NH tag (number of reported alignments) if present
    let num_hits = record.aux(b"NH").ok().and_then(|aux| match aux {
        rust_htslib::bam::record::Aux::U8(v) => Some(v as u32),
        rust_htslib::bam::record::Aux::U16(v) => Some(v as u32),
        rust_htslib::bam::record::Aux::U32(v) => Some(v),
        rust_htslib::bam::record::Aux::I8(v) => Some(v as u32),
        rust_htslib::bam::record::Aux::I16(v) => Some(v as u32),
        rust_htslib::bam::record::Aux::I32(v) => Some(v as u32),
        _ => None,
    });

    // Extract paired-end information from BAM flags
    let is_paired = record.is_paired();
    let is_first_in_pair = record.is_first_in_template();

    Ok(RnaSeqAlignment {
        read_name,
        chromosome,
        start,
        end,
        strand,
        cigar: cigar_string,
        mapping_quality,
        num_hits,
        splice_junctions,
        is_paired,
        is_first_in_pair,
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

/// Calculate alignment statistics for a BAM file with strand detection
pub fn get_alignment_stats<P: AsRef<Path>>(bam_path: P) -> Result<AlignmentStats> {
    get_alignment_stats_with_config(bam_path, 0.65, 10000)
}

/// Get basic alignment statistics without strand detection
pub fn get_basic_alignment_stats<P: AsRef<Path>>(bam_path: P) -> Result<AlignmentStats> {
    let mut bam_reader = BamReader::new(bam_path)?;
    let mut total_reads = 0;
    let mut mapped_reads = 0;
    let mut spliced_reads = 0;
    let mut total_splice_junctions = 0;
    let mut forward_reads = 0;
    let mut reverse_reads = 0;

    // Sample a reasonable number of reads for basic stats
    const MAX_SAMPLE_READS: u32 = 10000;
    let mut sampled = 0;

    for result in bam_reader.reader.records() {
        if sampled >= MAX_SAMPLE_READS {
            break;
        }

        let record = result.map_err(|e| AnnoRefineError::BamParse(e.to_string()))?;

        if record.is_unmapped() {
            continue;
        }

        total_reads += 1;
        mapped_reads += 1;
        sampled += 1;

        // Parse alignment to get strand and splice info
        match parse_bam_record(&record, &bam_reader.header) {
            Ok(alignment) => {
                match alignment.strand {
                    Strand::Forward => forward_reads += 1,
                    Strand::Reverse => reverse_reads += 1,
                    Strand::Unknown => continue,
                }

                if !alignment.splice_junctions.is_empty() {
                    spliced_reads += 1;
                    total_splice_junctions += alignment.splice_junctions.len() as u64;
                }
            }
            Err(_) => continue, // Skip problematic records
        }
    }

    // For basic stats, assume unstranded
    Ok(AlignmentStats {
        total_reads,
        mapped_reads,
        spliced_reads,
        total_splice_junctions,
        forward_reads,
        reverse_reads,
        strand_bias: StrandBias::Unstranded,
        strand_bias_ratio: 0.5,
    })
}

/// Calculate alignment statistics for a BAM file with gene-model-based strand detection
pub fn get_alignment_stats_with_gene_models<P: AsRef<Path>>(
    bam_path: P,
    gene_models: &[crate::types::GeneModel],
    strand_bias_threshold: f64,
    max_reads_for_strand_detection: u32,
) -> Result<AlignmentStats> {
    let bam_path = bam_path.as_ref();

    // Use gene-model-based strand detection for accurate sampling
    let strand_stats = detect_strand_bias_with_gene_models(
        bam_path,
        gene_models,
        max_reads_for_strand_detection,
        strand_bias_threshold,
    )?;

    Ok(strand_stats)
}

/// Legacy function for backward compatibility
pub fn get_alignment_stats_with_config<P: AsRef<Path>>(
    bam_path: P,
    strand_bias_threshold: f64,
    max_reads_for_strand_detection: u32,
) -> Result<AlignmentStats> {
    // Fallback to simple sequential sampling when no gene models available
    detect_strand_bias_sequential(
        bam_path,
        max_reads_for_strand_detection,
        strand_bias_threshold,
    )
}

/// Gene-model-based strand detection that samples from known gene regions
fn detect_strand_bias_with_gene_models<P: AsRef<Path>>(
    bam_path: P,
    gene_models: &[crate::types::GeneModel],
    max_reads_for_strand_detection: u32,
    strand_bias_threshold: f64,
) -> Result<AlignmentStats> {
    let bam_path = bam_path.as_ref();

    if gene_models.is_empty() {
        return Err(AnnoRefineError::BamParse(
            "No gene models provided for strand detection".to_string(),
        ));
    }

    debug!(
        "Using gene-model-based strand detection with {} gene models",
        gene_models.len()
    );

    // Sample alignments from selected gene models
    let sampled_alignments =
        sample_alignments_from_gene_models(bam_path, gene_models, max_reads_for_strand_detection)?;

    // Calculate statistics from sampled alignments
    calculate_stats_from_alignments(sampled_alignments, strand_bias_threshold)
}

/// Sequential strand detection (fallback for backward compatibility)
fn detect_strand_bias_sequential<P: AsRef<Path>>(
    bam_path: P,
    max_reads_for_strand_detection: u32,
    strand_bias_threshold: f64,
) -> Result<AlignmentStats> {
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
    let mut forward_reads = 0;
    let mut reverse_reads = 0;

    for result in reader.records() {
        let record = result
            .map_err(|e| AnnoRefineError::BamParse(format!("Failed to read BAM record: {}", e)))?;

        total_reads += 1;

        if !record.is_unmapped() && !record.is_secondary() && !record.is_duplicate() {
            mapped_reads += 1;

            // Count strand orientation
            if record.is_reverse() {
                reverse_reads += 1;
            } else {
                forward_reads += 1;
            }

            // Check for splice junctions
            if let Some(cigar) = record.cigar_cached() {
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

            // Stop early if we've sampled enough reads
            if mapped_reads >= max_reads_for_strand_detection as u64 {
                debug!(
                    "Sampled {} mapped reads for strand detection, stopping early",
                    mapped_reads
                );
                break;
            }
        }
    }

    // Calculate strand bias
    let (strand_bias, strand_bias_ratio) =
        detect_strand_bias(forward_reads, reverse_reads, strand_bias_threshold);

    Ok(AlignmentStats {
        total_reads,
        mapped_reads,
        spliced_reads,
        total_splice_junctions,
        forward_reads,
        reverse_reads,
        strand_bias,
        strand_bias_ratio,
    })
}

/// Sample alignments from forward strand genes only for strand detection
fn sample_alignments_from_gene_models<P: AsRef<Path>>(
    bam_path: P,
    gene_models: &[crate::types::GeneModel],
    max_reads: u32,
) -> Result<Vec<(RnaSeqAlignment, crate::types::Strand)>> {
    let mut bam_reader = BamReader::new(bam_path)?;
    let mut sampled_alignments = Vec::new();

    // Only use forward strand genes for strand detection
    let mut forward_genes: Vec<&crate::types::GeneModel> = gene_models
        .iter()
        .filter(|gene| gene.strand == crate::types::Strand::Forward)
        .collect();

    // Shuffle to get random sampling
    use rand::seq::SliceRandom;
    use rand::thread_rng;
    let mut rng = thread_rng();
    forward_genes.shuffle(&mut rng);

    // Target: ~50 forward strand genes (or as many as available)
    let target_genes = 50;
    let forward_sample_size = target_genes.min(forward_genes.len());

    debug!(
        "Sampling from {} forward strand genes for strand detection",
        forward_sample_size
    );

    if forward_sample_size == 0 {
        return Err(AnnoRefineError::BamParse(
            "No forward strand genes found for strand detection".to_string(),
        ));
    }

    let reads_per_gene = (max_reads as usize / forward_sample_size).max(10); // At least 10 reads per gene

    // Sample from forward strand genes only
    for gene in forward_genes.iter().take(forward_sample_size) {
        if sampled_alignments.len() >= max_reads as usize {
            break;
        }

        let gene_region = GenomicInterval {
            chromosome: gene.chromosome.clone(),
            start: gene.start,
            end: gene.end,
            strand: gene.strand,
        };

        match bam_reader.get_alignments_in_region(&gene_region) {
            Ok(alignments) => {
                // Filter for first-in-pair reads only (for paired-end stranded libraries)
                let first_in_pair_alignments: Vec<_> = alignments
                    .into_iter()
                    .filter(|alignment| {
                        // For strand detection, we only want first-in-pair reads from paired-end data
                        // or all reads from single-end data
                        !alignment.is_paired || alignment.is_first_in_pair
                    })
                    .take(reads_per_gene)
                    .collect();

                let count = first_in_pair_alignments.len();

                // Add gene strand information to each alignment
                for alignment in first_in_pair_alignments {
                    sampled_alignments.push((alignment, gene.strand));
                }

                debug!(
                    "Sampled {} alignments from forward gene {} ({}:{}-{})",
                    count, gene.id, gene.chromosome, gene.start, gene.end
                );
            }
            Err(e) => {
                debug!("Failed to get alignments from gene {}: {}", gene.id, e);
            }
        }
    }

    debug!(
        "Total sampled alignments for strand detection: {} from {} forward strand genes (first-in-pair only)",
        sampled_alignments.len(),
        forward_sample_size
    );

    Ok(sampled_alignments)
}

/// Calculate alignment statistics from forward-gene-only sampling for RF/FR/unstranded detection
fn calculate_stats_from_alignments(
    alignments_with_gene_strand: Vec<(RnaSeqAlignment, crate::types::Strand)>,
    strand_bias_threshold: f64,
) -> Result<AlignmentStats> {
    let total_reads = alignments_with_gene_strand.len() as u64;
    let mapped_reads = alignments_with_gene_strand.len() as u64; // All sampled alignments are mapped
    let mut spliced_reads = 0;
    let mut total_splice_junctions = 0;

    // Read orientation counters (all from forward strand genes)
    let mut forward_reads = 0; // Reads aligned in forward direction
    let mut reverse_reads = 0; // Reads aligned in reverse direction

    for (alignment, _gene_strand) in alignments_with_gene_strand {
        // Count read orientations (all genes are forward strand)
        match alignment.strand {
            Strand::Forward => forward_reads += 1,
            Strand::Reverse => reverse_reads += 1,
            Strand::Unknown => continue, // Skip unknown strand alignments
        }

        // Count splice junctions
        if !alignment.splice_junctions.is_empty() {
            spliced_reads += 1;
            total_splice_junctions += alignment.splice_junctions.len() as u64;
        }
    }

    // Calculate strand bias based on read orientations from forward genes
    let total_strand_informative = forward_reads + reverse_reads;
    let forward_ratio = if total_strand_informative > 0 {
        forward_reads as f64 / total_strand_informative as f64
    } else {
        0.5 // Default to unstranded if no informative reads
    };

    // Determine library type based on read orientations from forward strand genes:
    // - FR libraries: reads align in same direction as gene (forward) -> high forward_ratio
    // - RF libraries: reads align opposite to gene (reverse) -> low forward_ratio
    // - Unstranded: reads are ~50/50 forward/reverse -> forward_ratio ~0.5
    let strand_bias = if forward_ratio >= strand_bias_threshold {
        StrandBias::ForwardStranded // FR library (reads same direction as gene)
    } else if forward_ratio <= (1.0 - strand_bias_threshold) {
        StrandBias::ReverseStranded // RF library (reads opposite to gene)
    } else {
        StrandBias::Unstranded // Unstranded library (~50/50)
    };

    // Determine library type string for logging
    let library_type = if forward_ratio >= strand_bias_threshold {
        "FR (forward/reverse)"
    } else if forward_ratio <= (1.0 - strand_bias_threshold) {
        "RF (reverse/forward)"
    } else {
        "unstranded"
    };

    debug!(
        "Forward-gene strand detection (first-in-pair only): Forward reads: {}, Reverse reads: {}, Forward ratio: {:.3}",
        forward_reads, reverse_reads, forward_ratio
    );

    debug!(
        "Library type detected: {} ({:.1}% forward reads from forward genes, first-in-pair only)",
        library_type,
        forward_ratio * 100.0
    );

    Ok(AlignmentStats {
        total_reads,
        mapped_reads,
        spliced_reads,
        total_splice_junctions,
        forward_reads,
        reverse_reads,
        strand_bias,
        strand_bias_ratio: forward_ratio, // Use forward ratio as the bias ratio
    })
}

/// Detect strand bias from forward and reverse read counts
fn detect_strand_bias(forward_reads: u64, reverse_reads: u64, threshold: f64) -> (StrandBias, f64) {
    let total_mapped = forward_reads + reverse_reads;

    if total_mapped == 0 {
        return (StrandBias::Unstranded, 0.5);
    }

    let forward_ratio = forward_reads as f64 / total_mapped as f64;
    let reverse_ratio = reverse_reads as f64 / total_mapped as f64;

    let bias_ratio = forward_ratio.max(reverse_ratio);

    if bias_ratio >= threshold {
        if forward_ratio > reverse_ratio {
            (StrandBias::ForwardStranded, forward_ratio)
        } else {
            (StrandBias::ReverseStranded, reverse_ratio)
        }
    } else {
        (StrandBias::Unstranded, bias_ratio)
    }
}

/// Statistics about BAM alignments
#[derive(Debug)]
pub struct AlignmentStats {
    pub total_reads: u64,
    pub mapped_reads: u64,
    pub spliced_reads: u64,
    pub total_splice_junctions: u64,
    pub forward_reads: u64,
    pub reverse_reads: u64,
    pub strand_bias: StrandBias,
    pub strand_bias_ratio: f64,
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
            "Total reads: {}, Mapped: {} ({:.2}%), Spliced: {} ({:.2}%), Splice junctions: {}, Forward: {}, Reverse: {}, Strand bias: {} ({:.2}%)",
            self.total_reads,
            self.mapped_reads,
            mapping_rate,
            self.spliced_reads,
            splice_rate,
            self.total_splice_junctions,
            self.forward_reads,
            self.reverse_reads,
            self.strand_bias,
            self.strand_bias_ratio * 100.0
        )
    }
}

/// Filter alignments based on library type and gene orientation
/// This ensures we only use the appropriate reads for refinement
pub fn filter_alignments_by_library_type(
    alignments: Vec<RnaSeqAlignment>,
    gene_strand: Strand,
    library_type: crate::types::LibraryType,
) -> Vec<RnaSeqAlignment> {
    alignments
        .into_iter()
        .filter(|alignment| should_use_alignment_for_gene(alignment, gene_strand, library_type))
        .collect()
}

/// Determine if an alignment should be used for a gene based on library type and gene orientation
fn should_use_alignment_for_gene(
    alignment: &RnaSeqAlignment,
    gene_strand: Strand,
    library_type: crate::types::LibraryType,
) -> bool {
    use crate::types::LibraryType;

    match library_type {
        // Auto-detection: use all reads (strand detection will handle this)
        LibraryType::Auto => true,

        // Single-end libraries
        LibraryType::SingleForward => {
            // For single-end forward: use forward reads for forward genes, reverse reads for reverse genes
            match (gene_strand, alignment.strand) {
                (Strand::Forward, Strand::Forward) => true,
                (Strand::Reverse, Strand::Reverse) => true,
                _ => false,
            }
        }
        LibraryType::SingleReverse => {
            // For single-end reverse: use reverse reads for forward genes, forward reads for reverse genes
            match (gene_strand, alignment.strand) {
                (Strand::Forward, Strand::Reverse) => true,
                (Strand::Reverse, Strand::Forward) => true,
                _ => false,
            }
        }
        LibraryType::SingleUnstranded => {
            // For single-end unstranded: use all reads
            true
        }

        // Paired-end libraries (only use first-in-pair reads for strand determination)
        LibraryType::PairedFR => {
            // FR library: first-in-pair reads align same direction as gene
            if alignment.is_paired && !alignment.is_first_in_pair {
                return false; // Skip second-in-pair reads
            }
            match (gene_strand, alignment.strand) {
                (Strand::Forward, Strand::Forward) => true,
                (Strand::Reverse, Strand::Reverse) => true,
                _ => false,
            }
        }
        LibraryType::PairedRF => {
            // RF library: first-in-pair reads align opposite direction to gene
            if alignment.is_paired && !alignment.is_first_in_pair {
                return false; // Skip second-in-pair reads
            }
            match (gene_strand, alignment.strand) {
                (Strand::Forward, Strand::Reverse) => true,
                (Strand::Reverse, Strand::Forward) => true,
                _ => false,
            }
        }
        LibraryType::PairedUnstranded => {
            // For paired-end unstranded: use all first-in-pair reads (or all reads if single-end)
            if alignment.is_paired {
                alignment.is_first_in_pair
            } else {
                true
            }
        }
    }
}
