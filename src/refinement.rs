//! Gene model refinement using RNA-seq evidence

use crate::bam::BamReader;
use crate::types::{
    AlignmentCluster, AnnoRefineError, Exon, GeneModel, Genome, GenomicInterval, RefinementConfig,
    Result, SpliceJunction, Strand, Transcript,
};

use crate::logging::{log_gene_model_changes, log_splice_junction_refinement, log_utr_extension};

use log::{debug, info, warn};
use rayon::prelude::*;
use std::collections::HashMap;

/// Main refinement engine
pub struct RefinementEngine {
    config: RefinementConfig,
}

impl RefinementEngine {
    pub fn new(config: RefinementConfig) -> Self {
        Self { config }
    }

    /// Refine a collection of gene models using RNA-seq evidence
    pub fn refine_gene_models(
        &self,
        gene_models: &mut Vec<GeneModel>,
        bam_file_path: &std::path::Path,
        genome: &Genome,
    ) -> Result<RefinementSummary> {
        info!(
            "Starting gene model refinement for {} genes",
            gene_models.len()
        );

        // Group genes by chromosome for parallel processing
        let mut genes_by_chromosome: HashMap<String, Vec<usize>> = HashMap::new();
        for (index, gene_model) in gene_models.iter().enumerate() {
            genes_by_chromosome
                .entry(gene_model.chromosome.clone())
                .or_insert_with(Vec::new)
                .push(index);
        }

        info!(
            "Processing {} chromosomes in parallel",
            genes_by_chromosome.len()
        );

        // Note: We collect summaries from each thread instead of using shared state

        // Process each chromosome in parallel
        let chromosome_results: Result<Vec<_>> = genes_by_chromosome
            .par_iter()
            .map(|(chromosome, gene_indices)| {
                self.refine_genes_on_chromosome(
                    gene_models,
                    gene_indices,
                    chromosome,
                    bam_file_path,
                    genome,
                )
            })
            .collect();

        // Collect results and update gene models
        let chromosome_results = chromosome_results?;
        let mut final_summary = RefinementSummary::new();

        for (refined_genes, chrom_summary) in chromosome_results {
            // Update the original gene models with refined versions
            for (original_index, refined_gene) in refined_genes {
                gene_models[original_index] = refined_gene;
            }
            final_summary.merge(chrom_summary);
        }

        // 6. Novel gene detection (if enabled)
        if self.config.enable_novel_gene_detection {
            info!("Step 5: Detecting novel genes from RNA-seq evidence");
            let novel_genes = self.detect_novel_genes(gene_models, bam_file_path, genome)?;

            if !novel_genes.is_empty() {
                info!("Detected {} novel gene candidates", novel_genes.len());

                // Add novel genes to the gene models list
                for novel_gene in novel_genes {
                    gene_models.push(novel_gene);
                    final_summary.novel_genes_detected += 1;
                }
            } else {
                info!("No novel genes detected");
            }
        }

        info!(
            "Refinement complete: {} genes processed, {} failed, {} novel genes detected",
            final_summary.genes_processed,
            final_summary.genes_failed,
            final_summary.novel_genes_detected
        );

        Ok(final_summary)
    }

    /// Refine genes on a single chromosome
    fn refine_genes_on_chromosome(
        &self,
        gene_models: &[GeneModel],
        gene_indices: &[usize],
        chromosome: &str,
        bam_file_path: &std::path::Path,
        genome: &Genome,
    ) -> Result<(Vec<(usize, GeneModel)>, RefinementSummary)> {
        debug!(
            "Processing chromosome {} with {} genes",
            chromosome,
            gene_indices.len()
        );

        // Create a separate BAM reader for this thread
        let mut bam_reader = BamReader::new(bam_file_path)?;
        let mut summary = RefinementSummary::new();
        let mut refined_genes = Vec::new();

        // Process each gene on this chromosome
        for &gene_index in gene_indices {
            // Work with a mutable copy
            let mut gene_model = gene_models[gene_index].clone();

            match self.refine_single_gene(&mut gene_model, &mut bam_reader, genome) {
                Ok(gene_summary) => {
                    summary.merge(gene_summary);
                    summary.genes_processed += 1;

                    // Store the refined gene model with its original index
                    refined_genes.push((gene_index, gene_model));
                }
                Err(e) => {
                    warn!("Failed to refine gene {}: {}", gene_model.id, e);
                    summary.genes_failed += 1;

                    // Keep the original gene model if refinement failed
                    refined_genes.push((gene_index, gene_model));
                }
            }
        }

        debug!("Completed chromosome {} processing", chromosome);
        Ok((refined_genes, summary))
    }

    /// Refine a single gene model
    fn refine_single_gene(
        &self,
        gene_model: &mut GeneModel,
        bam_reader: &mut BamReader,
        genome: &Genome,
    ) -> Result<RefinementSummary> {
        debug!("Refining gene: {}", gene_model.id);

        let mut summary = RefinementSummary::new();

        // Get RNA-seq evidence for the gene region
        let gene_interval = gene_model.interval();
        let extended_interval = self.extend_interval_for_analysis(&gene_interval);

        let splice_junctions = bam_reader.get_splice_junctions_in_region(&extended_interval)?;
        let coverage = bam_reader.get_coverage_in_region(&extended_interval)?;

        // Track if any transcript has structural changes
        let mut gene_has_structural_changes = false;

        // Refine each transcript
        for transcript in &mut gene_model.transcripts {
            let transcript_summary = self.refine_transcript(
                transcript,
                &splice_junctions,
                &coverage,
                &extended_interval,
                genome,
                &gene_model.chromosome,
                gene_model.strand,
                &gene_model.id, // Pass gene ID for logging
            )?;

            // Check if this transcript had structural changes
            if transcript_summary.transcripts_with_structure_changes > 0 {
                gene_has_structural_changes = true;
            }

            summary.merge(transcript_summary);
        }

        // Mark gene if any transcript had structural changes
        if gene_has_structural_changes {
            gene_model.has_structural_changes = true;
            debug!("Gene {} marked as having structural changes", gene_model.id);
        }

        // Update gene boundaries based on refined transcripts
        gene_model.update_boundaries();

        Ok(summary)
    }

    /// Refine a single transcript
    fn refine_transcript(
        &self,
        transcript: &mut Transcript,
        splice_junctions: &[SpliceJunction],
        coverage: &[u32],
        region: &GenomicInterval,
        genome: &Genome,
        chromosome: &str,
        strand: Strand,
        gene_id: &str,
    ) -> Result<RefinementSummary> {
        let mut summary = RefinementSummary::new();

        // Store original transcript for CDS validation
        let _original_transcript = transcript.clone();

        // 1. Refine intron/exon structure based on splice junctions
        let structure_changed =
            self.refine_exon_structure(transcript, splice_junctions, gene_id, genome)?;

        // 2. Extend UTRs based on coverage (extend gene model first)
        let utr_changes = self.extend_utrs(transcript, coverage, region, gene_id)?;

        // 3. Re-predict CDS in the extended gene model
        // This may result in the same CDS, longer CDS, or completely different CDS
        let cds_changed = if structure_changed
            || utr_changes.five_prime_extended
            || utr_changes.three_prime_extended
        {
            self.re_predict_cds(transcript, genome, chromosome, strand, gene_id)?
        } else {
            false
        };

        // Update summary
        if structure_changed {
            summary.transcripts_with_structure_changes += 1;
        }
        if utr_changes.five_prime_extended {
            summary.transcripts_with_5utr_extension += 1;
        }
        if utr_changes.three_prime_extended {
            summary.transcripts_with_3utr_extension += 1;
        }
        if cds_changed {
            debug!("CDS re-predicted for transcript {}", transcript.id);
        }

        // Update transcript boundaries
        transcript.update_boundaries();

        Ok(summary)
    }

    /// Refine exon structure based on splice junction evidence
    fn refine_exon_structure(
        &self,
        transcript: &mut Transcript,
        splice_junctions: &[SpliceJunction],
        gene_id: &str,
        genome: &Genome,
    ) -> Result<bool> {
        if transcript.exons.len() < 2 {
            return Ok(false); // Single exon transcripts don't have splice junctions
        }

        let mut structure_changed = false;

        // Get current splice junctions from transcript
        let current_junctions = transcript.get_splice_junctions(Strand::Forward); // Strand will be corrected

        // Find well-supported splice junctions that could improve the model
        let mut supported_junctions = Vec::new();
        for junction in splice_junctions {
            if junction.support_count >= self.config.min_splice_support {
                supported_junctions.push(junction.clone());
            }
        }

        // For each current junction, check if there's better evidence nearby
        for i in 0..current_junctions.len() {
            let current_junction = &current_junctions[i];

            // Look for nearby supported junctions
            for supported in &supported_junctions {
                if self.junctions_are_compatible(current_junction, supported) {
                    // Validate splice sites before updating (if enabled)
                    if self.config.validate_splice_sites
                        && !self.validate_splice_junction(supported, genome)
                    {
                        debug!(
                            "Skipping junction update for transcript {} - invalid splice sites",
                            transcript.id
                        );
                        continue;
                    }

                    // Update exon boundaries if the junction is different
                    if current_junction.donor_pos != supported.donor_pos
                        || current_junction.acceptor_pos != supported.acceptor_pos
                    {
                        self.update_exon_boundaries(transcript, i, supported, gene_id)?;
                        structure_changed = true;
                        debug!(
                            "Updated splice junction for transcript {} with validated splice sites",
                            transcript.id
                        );
                    }
                }
            }
        }

        Ok(structure_changed)
    }

    /// Check if two splice junctions are compatible (close enough to be the same)
    fn junctions_are_compatible(&self, j1: &SpliceJunction, j2: &SpliceJunction) -> bool {
        let max_distance = 10; // Allow up to 10bp difference

        (j1.donor_pos as i64 - j2.donor_pos as i64).abs() <= max_distance
            && (j1.acceptor_pos as i64 - j2.acceptor_pos as i64).abs() <= max_distance
    }

    /// Validate splice junction with splice site motif checking
    fn validate_splice_junction(&self, junction: &SpliceJunction, genome: &Genome) -> bool {
        validate_splice_sites(
            junction.donor_pos,
            junction.acceptor_pos,
            junction.strand,
            &junction.chromosome,
            genome,
        )
    }

    /// Update exon boundaries based on a new splice junction
    fn update_exon_boundaries(
        &self,
        transcript: &mut Transcript,
        junction_index: usize,
        new_junction: &SpliceJunction,
        gene_id: &str,
    ) -> Result<()> {
        if junction_index >= transcript.exons.len() - 1 {
            return Err(AnnoRefineError::Refinement(
                "Junction index out of bounds".to_string(),
            ));
        }

        // Validate that the new junction makes sense
        if new_junction.donor_pos >= new_junction.acceptor_pos {
            return Err(AnnoRefineError::Refinement(format!(
                "Invalid splice junction: donor_pos ({}) >= acceptor_pos ({})",
                new_junction.donor_pos, new_junction.acceptor_pos
            )));
        }

        // Validate that the new boundaries won't create invalid exons
        let upstream_exon = &transcript.exons[junction_index];
        let downstream_exon = &transcript.exons[junction_index + 1];

        // Check that the new end position for upstream exon is valid
        if upstream_exon.start > new_junction.donor_pos {
            debug!("Skipping junction update: would create invalid upstream exon (start {} > new end {})",
                   upstream_exon.start, new_junction.donor_pos);
            return Ok(()); // Skip this update rather than create invalid exon
        }

        // Check that the new start position for downstream exon is valid
        if new_junction.acceptor_pos > downstream_exon.end {
            debug!("Skipping junction update: would create invalid downstream exon (new start {} > end {})",
                   new_junction.acceptor_pos, downstream_exon.end);
            return Ok(()); // Skip this update rather than create invalid exon
        }

        // Store original positions for logging
        let original_donor = transcript.exons[junction_index].end;
        let original_acceptor = transcript.exons[junction_index + 1].start;

        // Update the end of the upstream exon
        transcript.exons[junction_index].end = new_junction.donor_pos;

        // Update the start of the downstream exon
        transcript.exons[junction_index + 1].start = new_junction.acceptor_pos;

        // Log the change
        log_splice_junction_refinement(
            gene_id,
            &transcript.id,
            junction_index,
            original_donor,
            original_acceptor,
            new_junction.donor_pos,
            new_junction.acceptor_pos,
            new_junction.support_count,
        );

        debug!(
            "Updated exon boundaries for transcript {}: exon {} end -> {}, exon {} start -> {}",
            transcript.id,
            junction_index + 1,
            new_junction.donor_pos,
            junction_index + 2,
            new_junction.acceptor_pos
        );

        Ok(())
    }

    /// Calculate dynamic coverage threshold based on existing exon coverage
    fn calculate_dynamic_coverage_threshold(
        &self,
        transcript: &Transcript,
        coverage: &[u32],
        region: &GenomicInterval,
    ) -> Result<u32> {
        if transcript.exons.is_empty() {
            return Ok(self.config.min_coverage);
        }

        let mut exon_coverages = Vec::new();

        // Calculate average coverage for each exon
        for exon in &transcript.exons {
            // Skip exons that are outside our coverage region
            if exon.end < region.start || exon.start > region.end {
                continue;
            }

            let exon_start_offset = if exon.start >= region.start {
                (exon.start - region.start) as usize
            } else {
                0
            };

            let exon_end_offset = if exon.end <= region.end {
                (exon.end - region.start) as usize
            } else {
                coverage.len().saturating_sub(1)
            };

            if exon_start_offset < coverage.len()
                && exon_end_offset < coverage.len()
                && exon_start_offset <= exon_end_offset
            {
                let exon_coverage: u32 = coverage[exon_start_offset..=exon_end_offset]
                    .iter()
                    .sum::<u32>()
                    / (exon_end_offset - exon_start_offset + 1) as u32;

                if exon_coverage > 0 {
                    exon_coverages.push(exon_coverage);
                }
            }
        }

        if exon_coverages.is_empty() {
            debug!(
                "No exon coverage data available for transcript {}, using default threshold",
                transcript.id
            );
            return Ok(self.config.min_coverage);
        }

        // Calculate average coverage across all exons
        let average_exon_coverage =
            exon_coverages.iter().sum::<u32>() / exon_coverages.len() as u32;

        // Use 10% lower than average, but not less than the configured minimum
        let dynamic_threshold =
            ((average_exon_coverage as f64 * 0.9) as u32).max(self.config.min_coverage);

        debug!(
            "Transcript {}: average exon coverage = {}, dynamic threshold = {}",
            transcript.id, average_exon_coverage, dynamic_threshold
        );

        Ok(dynamic_threshold)
    }

    /// Extend UTRs based on coverage evidence
    fn extend_utrs(
        &self,
        transcript: &mut Transcript,
        coverage: &[u32],
        region: &GenomicInterval,
        gene_id: &str,
    ) -> Result<UtrExtensionResult> {
        let mut result = UtrExtensionResult {
            five_prime_extended: false,
            three_prime_extended: false,
        };

        if transcript.exons.is_empty() {
            return Ok(result);
        }

        // Calculate dynamic coverage threshold based on existing exons
        let dynamic_threshold =
            self.calculate_dynamic_coverage_threshold(transcript, coverage, region)?;

        // Determine transcript orientation
        let is_forward = region.strand == Strand::Forward;

        // Extend 5' UTR
        let five_prime_extension = if is_forward {
            self.find_upstream_extension(transcript, coverage, region, dynamic_threshold)?
        } else {
            self.find_downstream_extension(transcript, coverage, region, dynamic_threshold)?
        };

        // Only apply extensions that are meaningful (>= 5bp to avoid noise)
        const MIN_EXTENSION_LENGTH: u64 = 5;

        if five_prime_extension >= MIN_EXTENSION_LENGTH {
            self.apply_five_prime_extension(transcript, five_prime_extension, is_forward, gene_id)?;
            result.five_prime_extended = true;
            debug!(
                "Extended 5' UTR for transcript {} by {} bp using threshold {}",
                transcript.id, five_prime_extension, dynamic_threshold
            );
        } else if five_prime_extension > 0 {
            debug!(
                "Skipped 5' UTR extension for transcript {} ({} bp < {} bp minimum)",
                transcript.id, five_prime_extension, MIN_EXTENSION_LENGTH
            );
        }

        // Extend 3' UTR
        let three_prime_extension = if is_forward {
            self.find_downstream_extension(transcript, coverage, region, dynamic_threshold)?
        } else {
            self.find_upstream_extension(transcript, coverage, region, dynamic_threshold)?
        };

        if three_prime_extension >= MIN_EXTENSION_LENGTH {
            self.apply_three_prime_extension(
                transcript,
                three_prime_extension,
                is_forward,
                gene_id,
            )?;
            result.three_prime_extended = true;
            debug!(
                "Extended 3' UTR for transcript {} by {} bp using threshold {}",
                transcript.id, three_prime_extension, dynamic_threshold
            );
        } else if three_prime_extension > 0 {
            debug!(
                "Skipped 3' UTR extension for transcript {} ({} bp < {} bp minimum)",
                transcript.id, three_prime_extension, MIN_EXTENSION_LENGTH
            );
        }

        Ok(result)
    }

    /// Find extension upstream of the transcript
    fn find_upstream_extension(
        &self,
        transcript: &Transcript,
        coverage: &[u32],
        region: &GenomicInterval,
        coverage_threshold: u32,
    ) -> Result<u64> {
        let transcript_start = transcript.start;
        let region_start = region.start;

        if transcript_start <= region_start {
            return Ok(0);
        }

        let start_offset = (transcript_start - region_start) as usize;
        if start_offset >= coverage.len() {
            return Ok(0);
        }

        // Look upstream for continuous coverage above dynamic threshold
        let mut extension = 0;
        let mut consecutive_low_coverage = 0;
        const MAX_LOW_COVERAGE_GAP: usize = 2; // Allow small gaps in coverage (reduced from 3)

        for i in (0..start_offset).rev() {
            if coverage[i] >= coverage_threshold {
                extension += 1;
                consecutive_low_coverage = 0; // Reset gap counter
                if extension >= self.config.max_utr_extension as usize {
                    break;
                }
            } else {
                consecutive_low_coverage += 1;
                if consecutive_low_coverage > MAX_LOW_COVERAGE_GAP {
                    break; // Stop if we hit too many consecutive low-coverage positions
                }
                // Don't extend for low coverage positions - only count them as gaps
            }
        }

        Ok(extension as u64)
    }

    /// Find extension downstream of the transcript
    fn find_downstream_extension(
        &self,
        transcript: &Transcript,
        coverage: &[u32],
        region: &GenomicInterval,
        coverage_threshold: u32,
    ) -> Result<u64> {
        let transcript_end = transcript.end;
        let region_start = region.start;

        let end_offset = (transcript_end - region_start) as usize;
        if end_offset >= coverage.len() {
            return Ok(0);
        }

        // Look downstream for continuous coverage above dynamic threshold
        let mut extension = 0;
        let mut consecutive_low_coverage = 0;
        const MAX_LOW_COVERAGE_GAP: usize = 2; // Allow small gaps in coverage (reduced from 3)

        for i in (end_offset + 1)..coverage.len() {
            if coverage[i] >= coverage_threshold {
                extension += 1;
                consecutive_low_coverage = 0; // Reset gap counter
                if extension >= self.config.max_utr_extension as usize {
                    break;
                }
            } else {
                consecutive_low_coverage += 1;
                if consecutive_low_coverage > MAX_LOW_COVERAGE_GAP {
                    break; // Stop if we hit too many consecutive low-coverage positions
                }
                // Don't extend for low coverage positions - only count them as gaps
            }
        }

        Ok(extension as u64)
    }

    /// Apply 5' UTR extension
    fn apply_five_prime_extension(
        &self,
        transcript: &mut Transcript,
        extension: u64,
        is_forward: bool,
        gene_id: &str,
    ) -> Result<()> {
        if transcript.exons.is_empty() {
            return Ok(());
        }

        if is_forward {
            // Extend first exon upstream (ensure we don't go below 1)
            let original_start = transcript.exons[0].start;
            transcript.exons[0].start = transcript.exons[0].start.saturating_sub(extension).max(1);
            let actual_extension = original_start - transcript.exons[0].start;

            if actual_extension > 0 {
                log_utr_extension(
                    gene_id,
                    &transcript.id,
                    "5PRIME",
                    original_start,
                    transcript.exons[0].start,
                    actual_extension,
                );
            }
        } else {
            // Extend last exon downstream (reverse strand 5' UTR)
            let last_idx = transcript.exons.len() - 1;
            let original_end = transcript.exons[last_idx].end;
            transcript.exons[last_idx].end += extension;

            log_utr_extension(
                gene_id,
                &transcript.id,
                "5PRIME",
                original_end,
                transcript.exons[last_idx].end,
                extension,
            );
        }

        Ok(())
    }

    /// Apply 3' UTR extension
    fn apply_three_prime_extension(
        &self,
        transcript: &mut Transcript,
        extension: u64,
        is_forward: bool,
        gene_id: &str,
    ) -> Result<()> {
        if transcript.exons.is_empty() {
            return Ok(());
        }

        if is_forward {
            // Extend last exon downstream (forward strand 3' UTR)
            let last_idx = transcript.exons.len() - 1;
            let original_end = transcript.exons[last_idx].end;
            transcript.exons[last_idx].end += extension;

            log_utr_extension(
                gene_id,
                &transcript.id,
                "3PRIME",
                original_end,
                transcript.exons[last_idx].end,
                extension,
            );
        } else {
            // Extend first exon upstream (reverse strand 3' UTR)
            let original_start = transcript.exons[0].start;
            transcript.exons[0].start = transcript.exons[0].start.saturating_sub(extension).max(1);
            let actual_extension = original_start - transcript.exons[0].start;

            if actual_extension > 0 {
                log_utr_extension(
                    gene_id,
                    &transcript.id,
                    "3PRIME",
                    original_start,
                    transcript.exons[0].start,
                    actual_extension,
                );
            }
        }

        Ok(())
    }

    /// Extend the analysis interval to capture potential UTR extensions
    fn extend_interval_for_analysis(&self, interval: &GenomicInterval) -> GenomicInterval {
        GenomicInterval {
            chromosome: interval.chromosome.clone(),
            start: interval
                .start
                .saturating_sub(self.config.max_utr_extension as u64)
                .max(1),
            end: interval.end + self.config.max_utr_extension as u64,
            strand: interval.strand,
        }
    }

    /// Detect novel genes from RNA-seq evidence in regions without existing annotations
    fn detect_novel_genes(
        &self,
        existing_genes: &[GeneModel],
        bam_file_path: &std::path::Path,
        genome: &Genome,
    ) -> Result<Vec<GeneModel>> {
        info!("Scanning for novel gene candidates...");

        // Create a BAM reader for novel gene detection
        let mut bam_reader = BamReader::new(bam_file_path)?;
        let mut novel_genes = Vec::new();

        // Process each chromosome
        for (chromosome, sequence) in &genome.sequences {
            debug!("Scanning chromosome {} for novel genes", chromosome);

            // Get existing gene intervals on this chromosome
            let existing_intervals = self.get_existing_gene_intervals(existing_genes, chromosome);

            // Find uncovered regions with spliced alignments
            let alignment_clusters = self.find_alignment_clusters(
                &mut bam_reader,
                chromosome,
                sequence.sequence.len() as u64,
                &existing_intervals,
                genome,
            )?;

            // Convert clusters to gene candidates
            for cluster in alignment_clusters {
                if let Some(gene_candidate) = self.cluster_to_gene_candidate(cluster, genome)? {
                    novel_genes.push(gene_candidate);
                }
            }
        }

        info!(
            "Novel gene detection complete: {} candidates found",
            novel_genes.len()
        );
        Ok(novel_genes)
    }

    /// Get intervals covered by existing genes on a chromosome
    fn get_existing_gene_intervals(
        &self,
        genes: &[GeneModel],
        chromosome: &str,
    ) -> Vec<GenomicInterval> {
        genes
            .iter()
            .filter(|gene| gene.chromosome == chromosome)
            .map(|gene| GenomicInterval {
                chromosome: gene.chromosome.clone(),
                start: gene.start,
                end: gene.end,
                strand: gene.strand,
            })
            .collect()
    }

    /// Find clusters of alignments in uncovered regions (simplified implementation)
    fn find_alignment_clusters(
        &self,
        bam_reader: &mut BamReader,
        chromosome: &str,
        chromosome_length: u64,
        existing_intervals: &[GenomicInterval],
        genome: &Genome,
    ) -> Result<Vec<AlignmentCluster>> {
        debug!("Finding alignment clusters on chromosome {}", chromosome);

        // Get all splice junctions on this chromosome (both strands)
        let region = GenomicInterval {
            chromosome: chromosome.to_string(),
            start: 1,
            end: chromosome_length,
            strand: Strand::Forward, // This parameter is not used for junction extraction
        };
        let splice_junctions = bam_reader.get_splice_junctions_in_region(&region)?;

        debug!(
            "Found {} total splice junctions on chromosome {}",
            splice_junctions.len(),
            chromosome
        );

        // Debug: Check strand distribution
        let forward_count = splice_junctions
            .iter()
            .filter(|j| j.strand == Strand::Forward)
            .count();
        let reverse_count = splice_junctions
            .iter()
            .filter(|j| j.strand == Strand::Reverse)
            .count();
        let unknown_count = splice_junctions
            .iter()
            .filter(|j| j.strand == Strand::Unknown)
            .count();
        debug!(
            "Strand distribution: Forward={}, Reverse={}, Unknown={}",
            forward_count, reverse_count, unknown_count
        );

        // Filter out junctions that overlap with existing genes or have invalid splice sites
        let novel_junctions: Vec<_> = splice_junctions
            .into_iter()
            .filter(|junction| {
                !self.junction_overlaps_existing_genes(junction, existing_intervals)
                    && (!self.config.validate_splice_sites
                        || self.validate_splice_junction(junction, genome))
            })
            .collect();

        debug!("Found {} novel splice junctions", novel_junctions.len());

        // Debug: Check strand distribution of novel junctions
        let novel_forward = novel_junctions
            .iter()
            .filter(|j| j.strand == Strand::Forward)
            .count();
        let novel_reverse = novel_junctions
            .iter()
            .filter(|j| j.strand == Strand::Reverse)
            .count();
        let novel_unknown = novel_junctions
            .iter()
            .filter(|j| j.strand == Strand::Unknown)
            .count();
        debug!(
            "Novel junction strands: Forward={}, Reverse={}, Unknown={}",
            novel_forward, novel_reverse, novel_unknown
        );

        // Create clusters and remove duplicates at same coordinates
        let mut clusters: Vec<AlignmentCluster> = novel_junctions
            .into_iter()
            .map(|junction| AlignmentCluster {
                chromosome: chromosome.to_string(),
                start: junction.donor_pos.min(junction.acceptor_pos),
                end: junction.donor_pos.max(junction.acceptor_pos),
                strand: junction.strand,
                splice_junctions: vec![junction],
                coverage_profile: Vec::new(),
                alignment_count: 1,
            })
            .collect();

        // Remove clusters that overlap at the same coordinates but different strands
        // Keep the one with higher support
        clusters.sort_by_key(|c| (c.start, c.end, c.alignment_count));
        let mut filtered_clusters = Vec::new();
        let mut i = 0;

        while i < clusters.len() {
            let current = &clusters[i];
            let mut best_cluster = current.clone();
            let mut j = i + 1;

            // Check for overlapping clusters at same coordinates
            while j < clusters.len()
                && clusters[j].start == current.start
                && clusters[j].end == current.end
            {
                if clusters[j].alignment_count > best_cluster.alignment_count {
                    best_cluster = clusters[j].clone();
                }
                j += 1;
            }

            filtered_clusters.push(best_cluster);
            i = j;
        }

        let clusters = filtered_clusters;

        debug!("Formed {} alignment clusters", clusters.len());
        Ok(clusters)
    }

    /// Check if a splice junction overlaps with existing gene intervals
    fn junction_overlaps_existing_genes(
        &self,
        junction: &SpliceJunction,
        existing_intervals: &[GenomicInterval],
    ) -> bool {
        let junction_start = junction.donor_pos.min(junction.acceptor_pos);
        let junction_end = junction.donor_pos.max(junction.acceptor_pos);

        existing_intervals
            .iter()
            .any(|interval| junction_start < interval.end && junction_end > interval.start)
    }

    /// Convert an alignment cluster to a gene candidate (simplified implementation)
    fn cluster_to_gene_candidate(
        &self,
        cluster: AlignmentCluster,
        _genome: &Genome,
    ) -> Result<Option<GeneModel>> {
        // Check if cluster meets minimum requirements
        if cluster.end - cluster.start < self.config.min_novel_gene_length {
            debug!(
                "Cluster too short: {} bp < {} bp minimum",
                cluster.end - cluster.start,
                self.config.min_novel_gene_length
            );
            return Ok(None);
        }

        // Create a simple single-exon gene for now
        let gene_id = format!(
            "novel_{}_{}_{}_{}",
            cluster.chromosome,
            cluster.start,
            cluster.end,
            match cluster.strand {
                Strand::Forward => "plus",
                Strand::Reverse => "minus",
                Strand::Unknown => "unknown",
            }
        );

        let exon = Exon {
            start: cluster.start,
            end: cluster.end,
            id: Some(format!("{}-T1.exon1", gene_id)),
        };

        // Predict CDS regions for the novel gene
        let cds_regions =
            self.predict_cds_for_novel_gene(&gene_id, &[exon.clone()], cluster.strand, _genome)?;

        // Only create gene if we found a valid CDS
        if cds_regions.is_empty() {
            debug!("No valid CDS found for novel gene candidate {}", gene_id);
            return Ok(None);
        }

        // Create transcript
        let transcript = Transcript {
            id: format!("{}-T1", gene_id),
            name: Some("Novel transcript".to_string()),
            start: cluster.start,
            end: cluster.end,
            exons: vec![exon],
            cds_regions,
            five_prime_utr: None,
            three_prime_utr: None,
            original_attributes: std::collections::HashMap::new(),
        };

        // Create gene model
        let mut original_cds_lengths = std::collections::HashMap::new();
        original_cds_lengths.insert(transcript.id.clone(), 0); // Novel genes start with no CDS

        let gene_model = GeneModel {
            id: gene_id,
            name: Some("Novel gene".to_string()),
            chromosome: cluster.chromosome,
            start: cluster.start,
            end: cluster.end,
            strand: cluster.strand,
            transcripts: vec![transcript],
            original_attributes: std::collections::HashMap::new(),
            has_structural_changes: false,
            original_cds_lengths,
        };

        info!(
            "Created novel gene candidate: {} ({}:{}-{})",
            gene_model.id, gene_model.chromosome, gene_model.start, gene_model.end
        );

        Ok(Some(gene_model))
    }

    /// Predict CDS regions for a novel gene by finding the longest ORF
    fn predict_cds_for_novel_gene(
        &self,
        gene_id: &str,
        exons: &[Exon],
        strand: Strand,
        genome: &Genome,
    ) -> Result<Vec<crate::types::CdsRegion>> {
        if exons.is_empty() {
            return Ok(Vec::new());
        }

        // Get the chromosome sequence
        let chr_name = gene_id.split('_').nth(1).unwrap_or("unknown");

        let genome_seq = match genome.sequences.get(chr_name) {
            Some(seq) => &seq.sequence,
            None => {
                debug!("Could not find chromosome sequence for {}", chr_name);
                return Ok(Vec::new());
            }
        };

        // Extract exon sequences and concatenate
        let mut transcript_seq = String::new();

        for exon in exons {
            let start_idx = (exon.start - 1) as usize; // Convert to 0-based
            let end_idx = exon.end as usize;

            if end_idx <= genome_seq.len() && start_idx < end_idx {
                let exon_seq = &genome_seq[start_idx..end_idx];
                let exon_str = std::str::from_utf8(exon_seq).unwrap_or("");
                transcript_seq.push_str(exon_str);
            }
        }

        if transcript_seq.is_empty() {
            return Ok(Vec::new());
        }

        // Reverse complement if on minus strand
        if strand == Strand::Reverse {
            transcript_seq = reverse_complement(&transcript_seq);
        }

        // Find the longest ORF
        let longest_orf = find_longest_orf(&transcript_seq);

        if let Some((orf_start, orf_end)) = longest_orf {
            // Check if ORF meets minimum length requirement
            let orf_length = orf_end - orf_start;
            if orf_length < 150 {
                // Minimum 50 amino acids
                debug!(
                    "ORF too short for novel gene {}: {} bp",
                    gene_id, orf_length
                );
                return Ok(Vec::new());
            }

            // Convert transcript coordinates back to genomic coordinates
            let cds_regions =
                self.transcript_to_genomic_coords(orf_start, orf_end, exons, strand)?;

            info!(
                "Predicted CDS for novel gene {}: {} bp ORF",
                gene_id, orf_length
            );
            Ok(cds_regions)
        } else {
            debug!("No valid ORF found for novel gene {}", gene_id);
            Ok(Vec::new())
        }
    }

    /// Convert transcript coordinates to genomic CDS regions
    fn transcript_to_genomic_coords(
        &self,
        transcript_start: usize,
        transcript_end: usize,
        exons: &[Exon],
        _strand: Strand,
    ) -> Result<Vec<crate::types::CdsRegion>> {
        use crate::types::CdsRegion;

        let mut cds_regions = Vec::new();
        let mut transcript_pos = 0;
        let mut phase = 0u8;

        for exon in exons {
            let exon_length = (exon.end - exon.start) as usize;
            let exon_end_pos = transcript_pos + exon_length;

            // Check if this exon overlaps with the ORF
            if transcript_pos < transcript_end && exon_end_pos > transcript_start {
                let cds_start_in_exon = if transcript_pos < transcript_start {
                    transcript_start - transcript_pos
                } else {
                    0
                };

                let cds_end_in_exon = if exon_end_pos > transcript_end {
                    transcript_end - transcript_pos
                } else {
                    exon_length
                };

                if cds_start_in_exon < cds_end_in_exon {
                    let genomic_start = exon.start + cds_start_in_exon as u64;
                    let genomic_end = exon.start + cds_end_in_exon as u64;

                    cds_regions.push(CdsRegion {
                        start: genomic_start,
                        end: genomic_end,
                        phase: Some(phase),
                        id: None,
                    });

                    // Update phase for next CDS region
                    let cds_length = (cds_end_in_exon - cds_start_in_exon) % 3;
                    phase = (3 - (cds_length as u8 % 3)) % 3;
                }
            }

            transcript_pos = exon_end_pos;
        }

        Ok(cds_regions)
    }

    /// Re-predict CDS in an extended transcript model
    fn re_predict_cds(
        &self,
        transcript: &mut Transcript,
        genome: &Genome,
        chromosome: &str,
        strand: Strand,
        gene_id: &str,
    ) -> Result<bool> {
        // Store original CDS for comparison
        let original_cds = transcript.cds_regions.clone();

        // Extract the full transcript sequence from the extended exons
        let transcript_sequence =
            self.extract_transcript_sequence(transcript, genome, chromosome, strand)?;

        if transcript_sequence.is_empty() {
            debug!(
                "Could not extract transcript sequence for {}",
                transcript.id
            );
            return Ok(false);
        }

        // Find the longest ORF in the extended transcript
        let longest_orf = find_longest_orf(&transcript_sequence);

        if let Some((orf_start, orf_end)) = longest_orf {
            // Check if ORF meets minimum length requirement
            let orf_length = orf_end - orf_start;
            if orf_length < 150 {
                // Minimum 50 amino acids
                debug!(
                    "ORF too short for transcript {}: {} bp",
                    transcript.id, orf_length
                );
                // Keep original CDS if new ORF is too short
                return Ok(false);
            }

            // Convert transcript coordinates back to genomic coordinates
            let new_cds_regions =
                self.transcript_to_genomic_coords(orf_start, orf_end, &transcript.exons, strand)?;

            // Compare with original CDS
            let cds_changed = !cds_regions_equal(&original_cds, &new_cds_regions);

            if cds_changed {
                // Calculate length difference
                let original_length: u64 = original_cds.iter().map(|cds| cds.end - cds.start).sum();
                let new_length: u64 = new_cds_regions.iter().map(|cds| cds.end - cds.start).sum();
                let length_diff = (new_length as i64 - original_length as i64).abs() as u64;

                // Only log meaningful changes (>= 10bp difference or >5% change)
                let percent_change = if original_length > 0 {
                    (length_diff * 100) / original_length
                } else {
                    100
                };

                if length_diff >= 10 || percent_change >= 5 {
                    debug!("Significant CDS change for transcript {} in gene {}: {} bp → {} bp ({}bp diff)",
                           transcript.id, gene_id, original_length, new_length, length_diff);

                    // Log detailed change for significant changes only
                    log_gene_model_changes(
                        gene_id,
                        &transcript.id,
                        "CDS_REPREDICTION",
                        &format!(
                            "Original: {} bp, New: {} bp, Diff: {}bp ({}%)",
                            original_length, new_length, length_diff, percent_change
                        ),
                    );
                } else {
                    debug!(
                        "Minor CDS adjustment for transcript {}: {} bp → {} bp ({}bp diff)",
                        transcript.id, original_length, new_length, length_diff
                    );
                }

                // Update transcript with new CDS regardless of size
                transcript.cds_regions = new_cds_regions;
            }

            Ok(cds_changed)
        } else {
            debug!(
                "No valid ORF found in extended transcript {}",
                transcript.id
            );
            // Keep original CDS if no valid ORF found
            Ok(false)
        }
    }

    /// Extract transcript sequence from exons
    fn extract_transcript_sequence(
        &self,
        transcript: &Transcript,
        genome: &Genome,
        chromosome: &str,
        strand: Strand,
    ) -> Result<String> {
        let genome_seq = match genome.sequences.get(chromosome) {
            Some(seq) => &seq.sequence,
            None => {
                debug!("Could not find chromosome sequence for {}", chromosome);
                return Ok(String::new());
            }
        };

        // Extract and concatenate exon sequences
        let mut transcript_seq = String::new();

        for exon in &transcript.exons {
            let start_idx = (exon.start - 1) as usize; // Convert to 0-based
            let end_idx = exon.end as usize;

            if end_idx <= genome_seq.len() && start_idx < end_idx {
                let exon_seq = &genome_seq[start_idx..end_idx];
                let exon_str = std::str::from_utf8(exon_seq).unwrap_or("");
                transcript_seq.push_str(exon_str);
            }
        }

        // Reverse complement if on minus strand
        if strand == Strand::Reverse {
            transcript_seq = reverse_complement(&transcript_seq);
        }

        Ok(transcript_seq)
    }
}

/// Result of UTR extension analysis
#[derive(Debug, Default)]
struct UtrExtensionResult {
    five_prime_extended: bool,
    three_prime_extended: bool,
}

/// Summary of refinement results
#[derive(Debug, Default)]
pub struct RefinementSummary {
    pub genes_processed: u32,
    pub genes_failed: u32,
    pub transcripts_with_structure_changes: u32,
    pub transcripts_with_5utr_extension: u32,
    pub transcripts_with_3utr_extension: u32,
    pub novel_genes_detected: u32,
}

impl RefinementSummary {
    fn new() -> Self {
        Self::default()
    }

    fn merge(&mut self, other: RefinementSummary) {
        self.transcripts_with_structure_changes += other.transcripts_with_structure_changes;
        self.transcripts_with_5utr_extension += other.transcripts_with_5utr_extension;
        self.transcripts_with_3utr_extension += other.transcripts_with_3utr_extension;
        self.novel_genes_detected += other.novel_genes_detected;
    }
}

impl std::fmt::Display for RefinementSummary {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f,
            "Genes processed: {}, Failed: {}, Structure changes: {}, 5'UTR extensions: {}, 3'UTR extensions: {}",
            self.genes_processed, self.genes_failed, self.transcripts_with_structure_changes,
            self.transcripts_with_5utr_extension, self.transcripts_with_3utr_extension
        )
    }
}

/// Find the longest open reading frame in a DNA sequence
fn find_longest_orf(sequence: &str) -> Option<(usize, usize)> {
    let seq_bytes = sequence.as_bytes();
    let mut longest_orf = None;
    let mut max_length = 0;

    // Check all three reading frames
    for frame in 0..3 {
        let mut start_pos = None;

        for i in (frame..seq_bytes.len().saturating_sub(2)).step_by(3) {
            if i + 2 >= seq_bytes.len() {
                break;
            }
            let codon = &seq_bytes[i..i + 3];

            // Start codon (ATG)
            if codon == b"ATG" && start_pos.is_none() {
                start_pos = Some(i);
            }

            // Stop codons (TAA, TAG, TGA)
            if (codon == b"TAA" || codon == b"TAG" || codon == b"TGA") && start_pos.is_some() {
                let start = start_pos.unwrap();
                let length = i - start;

                if length > max_length {
                    max_length = length;
                    longest_orf = Some((start, i + 3)); // Include stop codon
                }

                start_pos = None;
            }
        }

        // Handle ORF that extends to end of sequence
        if let Some(start) = start_pos {
            let length = seq_bytes.len() - start;
            if length > max_length {
                max_length = length;
                longest_orf = Some((start, seq_bytes.len()));
            }
        }
    }

    longest_orf
}

/// Reverse complement a DNA sequence
fn reverse_complement(sequence: &str) -> String {
    sequence
        .chars()
        .rev()
        .map(|c| match c.to_ascii_uppercase() {
            'A' => 'T',
            'T' => 'A',
            'G' => 'C',
            'C' => 'G',
            'N' => 'N',
            _ => c,
        })
        .collect()
}

/// Compare two sets of CDS regions for equality
fn cds_regions_equal(cds1: &[crate::types::CdsRegion], cds2: &[crate::types::CdsRegion]) -> bool {
    if cds1.len() != cds2.len() {
        return false;
    }

    // Sort both by start position for comparison
    let mut sorted_cds1 = cds1.to_vec();
    let mut sorted_cds2 = cds2.to_vec();
    sorted_cds1.sort_by_key(|cds| cds.start);
    sorted_cds2.sort_by_key(|cds| cds.start);

    // Compare each CDS region
    for (c1, c2) in sorted_cds1.iter().zip(sorted_cds2.iter()) {
        if c1.start != c2.start || c1.end != c2.end {
            return false;
        }
    }

    true
}

/// Validate splice site motifs for a splice junction
fn validate_splice_sites(
    donor_pos: u64,
    acceptor_pos: u64,
    strand: Strand,
    chromosome: &str,
    genome: &Genome,
) -> bool {
    let genome_seq = match genome.sequences.get(chromosome) {
        Some(seq) => &seq.sequence,
        None => return false,
    };

    // Convert to 0-based coordinates
    let donor_idx = (donor_pos - 1) as usize;
    let acceptor_idx = (acceptor_pos - 1) as usize;

    // Check bounds
    if donor_idx + 2 >= genome_seq.len() || acceptor_idx < 2 {
        return false;
    }

    match strand {
        Strand::Forward => {
            // For forward strand: donor site should be GT, acceptor site should be AG
            let donor_motif = get_sequence_at_position(genome_seq, donor_idx, 2);
            let acceptor_motif = get_sequence_at_position(genome_seq, acceptor_idx - 2, 2);

            // Check canonical splice sites
            is_valid_donor_site(&donor_motif) && is_valid_acceptor_site(&acceptor_motif)
        }
        Strand::Reverse => {
            // For reverse strand: coordinates are flipped
            // Donor site (3' end of intron) should be AC, acceptor site (5' end) should be CT
            let donor_motif = get_sequence_at_position(genome_seq, donor_idx - 2, 2);
            let acceptor_motif = get_sequence_at_position(genome_seq, acceptor_idx, 2);

            // Reverse complement to get the actual motifs
            let donor_rc = reverse_complement(&donor_motif);
            let acceptor_rc = reverse_complement(&acceptor_motif);

            is_valid_donor_site(&acceptor_rc) && is_valid_acceptor_site(&donor_rc)
        }
        Strand::Unknown => false, // Cannot validate unknown strand
    }
}

/// Get sequence at a specific position
fn get_sequence_at_position(genome_seq: &[u8], start: usize, length: usize) -> String {
    if start + length > genome_seq.len() {
        return String::new();
    }

    std::str::from_utf8(&genome_seq[start..start + length])
        .unwrap_or("")
        .to_uppercase()
}

/// Check if a sequence is a valid donor site (GT, GC, or AT)
fn is_valid_donor_site(motif: &str) -> bool {
    matches!(motif, "GT" | "GC" | "AT")
}

/// Check if a sequence is a valid acceptor site (AG or AC)
fn is_valid_acceptor_site(motif: &str) -> bool {
    matches!(motif, "AG" | "AC")
}
