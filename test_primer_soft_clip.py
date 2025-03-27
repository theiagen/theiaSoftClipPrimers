#!/usr/bin/env python3

import unittest
import os
import tempfile
import pysam
from primer_soft_clip import (
    PrimerRegion,
    read_primer_bed,
    get_matching_primers,
    process_bam
)

class TestPrimerSoftClip(unittest.TestCase):
    def setUp(self):
        # Create temporary files for testing
        self.temp_dir = tempfile.mkdtemp()
        
        # Define header for BAM files
        self.header = {'HD': {'VN': '1.0'},
                      'SQ': [{'LN': 4000000, 'SN': 'Chromosome'}]}  # Adjusted length for M. tuberculosis
        
        # Create a test BED file with real primer data
        self.bed_content = """Chromosome\t2725853\t2725873\tahpC_F\t20\t+
Chromosome\t2726437\t2726457\tahpC_R\t20\t-
Chromosome\t1460914\t1460934\tatpE_F\t20\t+
Chromosome\t1461381\t1461401\tatpE_R\t20\t-
Chromosome\t1046049\t1046071\teis_F\t21\t+
Chromosome\t2715455\t2715477\teis_F\t22\t-
Chromosome\t2714908\t2714928\tMDL_eis-R1\t20\t+"""
        self.bed_file = os.path.join(self.temp_dir, "test_primers.bed")
        with open(self.bed_file, "w") as f:
            f.write(self.bed_content)
            
        # Create a test BAM file
        self.unsorted_bam = os.path.join(self.temp_dir, "test.unsorted.bam")
        self.bam_file = os.path.join(self.temp_dir, "test.bam")
        
        # Write unsorted BAM
        with pysam.AlignmentFile(self.unsorted_bam, "wb", header=self.header) as bam:
            # Create test reads that overlap with primers
            def create_test_read(name, start, is_reverse=False):
                read = pysam.AlignedSegment()
                read.query_name = name
                read.query_sequence = "ATCG" * 10
                read.flag = 16 if is_reverse else 0
                read.reference_id = 0
                read.reference_start = start
                read.mapping_quality = 20
                read.cigartuples = [(0, 40)]  # 40M - this determines reference_end
                read.query_qualities = pysam.qualitystring_to_array("I" * 40)
                read.template_length = 0
                read.next_reference_id = -1
                read.next_reference_start = -1
                return read
            
            # Forward read overlapping with ahpC_F
            bam.write(create_test_read("read1", 2725853))
            
            # Reverse read overlapping with ahpC_R
            bam.write(create_test_read("read2", 2726437, is_reverse=True))
            
            # Forward read overlapping with atpE_F
            bam.write(create_test_read("read3", 1460914))
            
            # Read with existing soft clips
            soft_clipped_read = create_test_read("read4", 2725853)
            soft_clipped_read.query_sequence = "ATCG" * 25  # 100bp read
            soft_clipped_read.query_qualities = pysam.qualitystring_to_array("I" * 100)
            soft_clipped_read.cigartuples = [(4, 2), (0, 98)]  # 2S98M
            bam.write(soft_clipped_read)
        
        # Sort the BAM file
        pysam.sort("-o", self.bam_file, self.unsorted_bam)
        
        # Index the sorted BAM
        pysam.index(self.bam_file)
        
        # Remove unsorted BAM
        os.remove(self.unsorted_bam)

    def tearDown(self):
        # Clean up temporary files
        for f in os.listdir(self.temp_dir):
            os.remove(os.path.join(self.temp_dir, f))
        os.rmdir(self.temp_dir)

    def test_primer_region_creation(self):
        """Test PrimerRegion class creation with real primer data"""
        primer = PrimerRegion("Chromosome", 2725853, 2725873, "ahpC_F", "+")
        self.assertEqual(primer.chrom, "Chromosome")
        self.assertEqual(primer.start, 2725853)
        self.assertEqual(primer.end, 2725873)
        self.assertEqual(primer.name, "ahpC_F")
        self.assertEqual(primer.strand, "+")

    def test_read_primer_bed(self):
        """Test reading primer BED file with real primer data"""
        primers = read_primer_bed(self.bed_file)
        self.assertEqual(len(primers["Chromosome"]), 7)  # Total number of primers
        self.assertEqual(primers["Chromosome"][0].name, "ahpC_F")
        self.assertEqual(primers["Chromosome"][1].name, "ahpC_R")
        self.assertEqual(primers["Chromosome"][2].name, "atpE_F")

    def test_get_matching_primers(self):
        """Test finding matching primers with real primer coordinates"""
        primers = read_primer_bed(self.bed_file)
        
        def create_test_read(start, length, is_reverse=False):
            read = pysam.AlignedSegment()
            read.reference_id = 0
            read.reference_start = start
            read.is_reverse = is_reverse
            read.cigartuples = [(0, length)]  # length M (matches)
            read.query_sequence = "A" * length
            return read
        
        # Test cases with debug output
        def test_primer_match(read, primers, expected_match=None):
            read_end = read.reference_start + read.query_length
            matching = get_matching_primers(read, primers, read.reference_start, read_end)
            if expected_match:
                print(f"""
                Read: start={read.reference_start}, end={read_end}, is_reverse={read.is_reverse}
                Expected primer: {expected_match}
                Matching primers: {[p.name for p in matching]}
                """)
            return matching

        # Test forward read overlapping with ahpC_F (2725853-2725873)
        forward_read = create_test_read(2725853, 40, is_reverse=False)
        matching = test_primer_match(forward_read, primers["Chromosome"], "ahpC_F")
        self.assertEqual(len(matching), 1, "Should match ahpC_F primer")
        self.assertEqual(matching[0].name, "ahpC_F")
        
        # Test reverse read overlapping with ahpC_R (2726437-2726457)
        # For reverse reads, need to position so that read_end matches primer_end
        reverse_read = create_test_read(2726417, 40, is_reverse=True)  # Will end at 2726457
        matching = test_primer_match(reverse_read, primers["Chromosome"], "ahpC_R")
        self.assertEqual(len(matching), 1, "Should match ahpC_R primer")
        self.assertEqual(matching[0].name, "ahpC_R")
        
        # Test no matches for wrong direction
        wrong_dir_read = create_test_read(2725853, 40, is_reverse=True)
        matching = test_primer_match(wrong_dir_read, primers["Chromosome"])
        self.assertEqual(len(matching), 0, "Should not match when direction is wrong")

        # Print all primers for debugging
        print("\nAll primers:")
        for p in primers["Chromosome"]:
            print(f"Primer {p.name}: start={p.start}, end={p.end}, strand={p.strand}")

    def test_process_bam(self):
        """Test BAM processing with real primer data"""
        output_bam = os.path.join(self.temp_dir, "output.bam")
        screened_reads = os.path.join(self.temp_dir, "screened.txt")
        primers = read_primer_bed(self.bed_file)
        
        process_bam(self.bam_file, output_bam, primers, screened_reads)
        
        # Check if output files exist
        self.assertTrue(os.path.exists(output_bam))
        self.assertTrue(os.path.exists(screened_reads))
        
        # Read output BAM without using index
        with pysam.AlignmentFile(output_bam, "rb") as bam:
            reads = list(bam.fetch(until_eof=True))
            self.assertTrue(len(reads) > 0)

    def test_process_bam_with_existing_soft_clips(self):
        """Test BAM processing with reads that already have soft clips"""
        output_bam = os.path.join(self.temp_dir, "output_soft_clipped.bam")
        screened_reads = os.path.join(self.temp_dir, "screened_soft_clipped.txt")
        primers = read_primer_bed(self.bed_file)
        
        process_bam(self.bam_file, output_bam, primers, screened_reads)
        
        # Read the processed BAM and check the CIGAR string
        with pysam.AlignmentFile(output_bam, "rb") as bam:
            for read in bam.fetch():
                if read.query_name == "read4":
                    # Check if the soft clips were properly combined
                    self.assertEqual(read.cigartuples[0], (4, 22))  # 2 + 20 = 22 bases soft clipped
                    self.assertEqual(read.cigartuples[1], (0, 78))  # Remaining bases

if __name__ == '__main__':
    unittest.main() 