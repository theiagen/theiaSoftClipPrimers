#!/usr/bin/env python3
# Thanh Le Viet - Theiagen Genomics - 2025
# This script soft clips primer regions in a BAM file based on the direction of the read and the primer.


import pysam
import argparse
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional
import logging
import os

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class PrimerRegion:
    """Dataclass to store primer information."""

    chrom: str
    start: int
    end: int
    name: str
    strand: str


def read_primer_bed(bed_file: str) -> Dict[str, List[PrimerRegion]]:
    """Read primer regions from BED file."""
    primers_by_chrom = defaultdict(list)
    with open(bed_file, "r") as f:
        for line in f:
            fields = line.strip().split("\t")
            chrom, start, end, name = fields[0:4]
            strand = fields[5] if len(fields) > 5 else "+"
            primer = PrimerRegion(chrom, int(start), int(end), name, strand)
            primers_by_chrom[chrom].append(primer)
    return primers_by_chrom


def get_matching_primers(
    read: pysam.AlignedSegment,
    primers: List[PrimerRegion],
    read_start: int,
    read_end: int,
) -> List[PrimerRegion]:
    """Find primers that match read direction and overlap the read."""
    read_is_forward = not read.is_reverse
    matching = []

    for primer in primers:
        # Check if read direction matches primer direction
        if read_is_forward == (primer.strand == "+"):
            if read_is_forward:
                # For forward reads: only clip if the start (5' end) overlaps with primer
                # Ignore if only the end overlaps
                if read_start <= primer.end and read_start >= primer.start:
                    matching.append(primer)
            else:
                # For reverse reads: only clip if the end (3' end) overlaps with primer
                # Ignore if only the start overlaps
                if read_end >= primer.start and read_end <= primer.end:
                    matching.append(primer)
    return matching


def get_soft_clip_info(read: pysam.AlignedSegment) -> Tuple[int, int, int]:
    """Get existing soft clips and aligned length from read."""
    left_clip = right_clip = 0
    if read.cigartuples:
        # Check for soft clips (CIGAR code 4) at start of read
        if read.cigartuples[0][0] == 4:
            left_clip = read.cigartuples[0][1]
        # Check for soft clips at end of read
        if read.cigartuples[-1][0] == 4:
            right_clip = read.cigartuples[-1][1]

    seq_len = len(read.query_sequence) if read.query_sequence else 0
    aligned_len = seq_len - (left_clip + right_clip)
    return left_clip, right_clip, aligned_len


def create_modified_read(
    read: pysam.AlignedSegment, header: Dict
) -> pysam.AlignedSegment:
    """Create a copy of the read with all fields."""
    modified = pysam.AlignedSegment(header)
    for attr in [
        "query_name",
        "flag",
        "reference_id",
        "reference_start",
        "mapping_quality",
        "cigartuples",
        "next_reference_id",
        "next_reference_start",
        "template_length",
        "query_sequence",
        "query_qualities",
    ]:
        setattr(modified, attr, getattr(read, attr))

    for tag, value in read.get_tags():
        try:
            modified.set_tag(tag, value)
        except ValueError as e:
            logger.warning(f"Could not set tag {tag} for read {read.query_name}: {e}")

    return modified


def calculate_new_cigar(
    read_is_forward: bool,
    aligned_len: int,
    additional_clip: int,
    left_clip: int,
    right_clip: int,
) -> List[Tuple[int, int]]:
    """Calculate new CIGAR string based on clipping parameters."""
    if read_is_forward:
        total_left = left_clip + additional_clip
        new_cigar = [(4, total_left)]
        remaining = aligned_len - additional_clip
        if remaining > 0:
            new_cigar.append((0, remaining))
        if right_clip > 0:
            new_cigar.append((4, right_clip))
    else:
        new_cigar = []
        if left_clip > 0:
            new_cigar.append((4, left_clip))
        remaining = aligned_len - additional_clip
        if remaining > 0:
            new_cigar.append((0, remaining))
        new_cigar.append((4, right_clip + additional_clip))

    return new_cigar


def process_read(
    read: pysam.AlignedSegment,
    bam_in: pysam.AlignmentFile,
    primers_by_chrom: Dict[str, List[PrimerRegion]],
    screened_reads: Optional[str] = None,
) -> Tuple[pysam.AlignedSegment, bool]:
    """Process a single read and return modified read if needed."""
    if read.is_unmapped:
        return read, False

    chrom = bam_in.get_reference_name(read.reference_id)
    read_start = read.reference_start
    read_end = read.reference_end

    matching_primers = get_matching_primers(
        read, primers_by_chrom.get(chrom, []), read_start, read_end
    )

    if not matching_primers:
        return read, False

    modified = create_modified_read(read, bam_in.header)
    left_clip, right_clip, aligned_len = get_soft_clip_info(read)
    read_is_forward = not read.is_reverse

    if read_is_forward:
        max_primer_end = max(p.end for p in matching_primers)
        additional_clip = min(max_primer_end - read_start, aligned_len)
        new_start = read_start + additional_clip
    else:
        min_primer_start = min(p.start for p in matching_primers)
        additional_clip = min(read_end - min_primer_start, aligned_len)
        new_start = read_start

    new_cigar = calculate_new_cigar(
        read_is_forward, aligned_len, additional_clip, left_clip, right_clip
    )

    # Verify CIGAR length
    cigar_length = sum(length for _, length in new_cigar)
    if cigar_length != len(read.query_sequence or ""):
        logger.warning(f"CIGAR length mismatch for read {read.query_name}")
        return read, False

    modified.cigartuples = new_cigar
    modified.reference_start = new_start

    if screened_reads:
        for primer in matching_primers:
            screened_reads.write(
                f"{read.query_name}\t{chrom}\t{read_start}\t{read_end}\t{primer.name}\t"
                f"{primer.strand}\t{'+' if read_is_forward else '-'}\tsoft_clipped\n"
            )

    return modified, True


def process_bam(
    input_bam: str,
    output_bam: str,
    primers_by_chrom: Dict[str, List[PrimerRegion]],
    screened_reads_file: str,
):
    """Process BAM file and soft clip primer regions."""
    # Initialize counters
    total_reads = 0
    soft_clipped_reads = 0

    with (
        pysam.AlignmentFile(input_bam, "rb") as bam_in,
        open(screened_reads_file, "w") as screened_reads,
        pysam.AlignmentFile(output_bam, "wb", header=bam_in.header) as bam_out,
    ):
        # Write header for screened reads
        screened_reads.write(
            "read_name\tchromosome\tstart\tend\tprimer_name\tprimer_direction\tread_direction\taction\n"
        )

        # Process reads in chunks to save memory
        chunk_size = 100000
        reads_chunk = []

        for read in bam_in.fetch(until_eof=True):
            total_reads += 1
            modified_read, was_modified = process_read(
                read, bam_in, primers_by_chrom, screened_reads
            )
            if was_modified:
                soft_clipped_reads += 1
            reads_chunk.append((modified_read, read.is_unmapped))

            if len(reads_chunk) >= chunk_size:
                write_sorted_chunk(reads_chunk, bam_out)
                reads_chunk = []

        # Write remaining reads
        if reads_chunk:
            write_sorted_chunk(reads_chunk, bam_out)

    # Calculate proportion
    proportion = (soft_clipped_reads / total_reads) * 100 if total_reads > 0 else 0

    # Print statistics
    logger.info(f"Total reads processed: {total_reads:,}")
    logger.info(f"Reads soft clipped: {soft_clipped_reads:,}")
    logger.info(f"Proportion of reads soft clipped: {proportion:.2f}%")

    # Sort and index the final BAM file
    try:
        pysam.sort("-o", output_bam + ".sorted.bam", output_bam)
        os.rename(output_bam + ".sorted.bam", output_bam)
        pysam.index(output_bam)
    except Exception as e:
        logger.error(f"Failed to sort/index BAM file: {e}")


def write_sorted_chunk(
    reads: List[Tuple[pysam.AlignedSegment, bool]], bam_out: pysam.AlignmentFile
):
    """Write a sorted chunk of reads to the output BAM file."""
    # Separate mapped and unmapped reads
    mapped = [(r, pos) for r, pos in reads if not pos]
    unmapped = [r for r, pos in reads if pos]

    # Sort mapped reads by position
    mapped.sort(key=lambda x: (x[0].reference_id, x[0].reference_start))

    # Write sorted mapped reads
    for read, _ in mapped:
        bam_out.write(read)

    # Write unmapped reads
    for read in unmapped:
        bam_out.write(read)


def main():
    parser = argparse.ArgumentParser(description="Soft clip primer regions in BAM file")
    parser.add_argument("-i", "--input_bam", required=True, help="Input BAM file")
    parser.add_argument("-o", "--output_bam", required=True, help="Output BAM file")
    parser.add_argument("-p", "--primers", required=True, help="Primers BED file")
    parser.add_argument(
        "-s", "--screened", required=True, help="Output file for screened reads"
    )

    args = parser.parse_args()

    logger.info("Reading primer regions...")
    primers_by_chrom = read_primer_bed(args.primers)

    logger.info("Processing BAM file...")
    process_bam(args.input_bam, args.output_bam, primers_by_chrom, args.screened)

    logger.info("Done!")


if __name__ == "__main__":
    main()
