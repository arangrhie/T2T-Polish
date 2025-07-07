#!/bin/python3

"""
BAM file filter based on read length and alignment identity.

This script reads a BAM file, calculates alignment identity from CIGAR and NM tags,
and outputs only reads with length > 10kbp and identity > 85%.

Requirements:
    pysam

Usage:
    python filt_bam_len_identity.py input.bam output.bam
"""

import sys
import argparse
import pysam
from typing import Optional

def parse_cigar_for_aligned_length(cigar_tuples):
    """
    Parse CIGAR string to get aligned length and number of matches/mismatches.
    
    Args:
        cigar_tuples: List of (operation, length) tuples from pysam
        
    Returns:
        aligned_length
    """
    aligned_length = 0
    
    for operation, length in cigar_tuples:
        # CIGAR operations:
        # 0: M (match/mismatch)
        # 1: I (insertion)
        # 2: D (deletion) 
        # 3: N (skipped region)
        # 4: S (soft clipping)
        # 5: H (hard clipping)
        # 6: P (padding)
        # 7: = (sequence match)
        # 8: X (sequence mismatch)
        
        if operation in [0, 1, 2, 7, 8]:  # M, I, D, =, X : consume reference except I
            aligned_length += length
            
    return aligned_length

def calculate_identity(read) -> Optional[float]:
    """
    Calculate alignment identity using CIGAR and NM tag.
    
    Identity = (aligned_length - edit_distance) / aligned_length * 100
    
    Args:
        read: pysam AlignedSegment object
        
    Returns:
        Identity percentage or None if calculation not possible
    """
    if read.is_unmapped:
        return None
        
    # Get NM tag (edit distance)
    try:
        nm = read.get_tag('NM')
    except KeyError:
        print(f"Warning: Read {read.query_name} missing NM tag", file=sys.stderr)
        return None
    
    # Parse CIGAR to get aligned length
    if read.cigartuples is None:
        return None
        
    aligned_length = parse_cigar_for_aligned_length(read.cigartuples)
    
    if aligned_length == 0:
        return None
    
    # Calculate identity
    # Identity = (aligned bases - edit distance) / aligned bases * 100
    identity = ((aligned_length - nm) / aligned_length) * 100
    
    return max(0.0, min(100.0, identity))  # Clamp between 0-100%

def filter_bam(input_bam: str, output_bam: str, min_length: int = 10000, min_identity: float = 85.0):
    """
    Filter BAM file based on read length and alignment identity.
    
    Args:
        input_bam: Path to input BAM file
        output_bam: Path to output BAM file
        min_length: Minimum read length (default: 10000 bp)
        min_identity: Minimum identity percentage (default: 85.0%)
    """
    
    # Statistics counters
    total_reads = 0
    passed_length = 0
    passed_identity = 0
    passed_both = 0
    unmapped_reads = 0
    
    with pysam.AlignmentFile(input_bam, "rb") as infile:
        with pysam.AlignmentFile(output_bam, "wb", template=infile) as outfile:
            
            for read in infile:
                total_reads += 1
                
                # Skip unmapped reads
                if read.is_unmapped:
                    unmapped_reads += 1
                    continue
                
                # Check read length
                read_length = read.query_length
                if read_length is None or read_length < min_length:
                    continue
                passed_length += 1
                
                # Calculate identity
                identity = calculate_identity(read)
                if identity is None:
                    continue
                    
                if identity < min_identity:
                    continue
                passed_identity += 1
                
                # Read passes both filters
                passed_both += 1
                outfile.write(read)
                
                # Print progress for long files
                if total_reads % 100000 == 0:
                    print(f"Processed {total_reads:,} reads...", file=sys.stderr)
    
    # Print statistics
    print(f"\n=== Filtering Statistics ===")
    print(f"Total reads processed: {total_reads:,}")
    print(f"Unmapped reads: {unmapped_reads:,}")
    print(f"Reads ≥ {min_length:,} bp: {passed_length:,} ({passed_length/total_reads*100:.1f}%)")
    print(f"Reads ≥ {min_identity}% identity: {passed_identity:,}")
    print(f"Reads passing both filters: {passed_both:,} ({passed_both/total_reads*100:.1f}%)")
    print(f"Output written to: {output_bam}")

def main():
    parser = argparse.ArgumentParser(
        description="Filter BAM file by read length and alignment identity",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python filter_bam_by_identity.py input.bam output.bam
    python filter_bam_by_identity.py input.bam output.bam --min-length 15000 --min-identity 90
        """
    )
    
    parser.add_argument("input_bam", help="Input BAM file")
    parser.add_argument("output_bam", help="Output BAM file")
    parser.add_argument("--min-length", type=int, default=10000,
                       help="Minimum read length in bp (default: 10000)")
    parser.add_argument("--min-identity", type=float, default=85.0,
                       help="Minimum identity percentage (default: 85.0)")
    parser.add_argument("--verbose", "-v", action="store_true",
                       help="Print detailed statistics for each read")
    
    args = parser.parse_args()
    
    # Validate arguments
    if args.min_length <= 0:
        print("Error: Minimum length must be positive", file=sys.stderr)
        sys.exit(1)
        
    if not (0 <= args.min_identity <= 100):
        print("Error: Identity must be between 0 and 100", file=sys.stderr)
        sys.exit(1)
    
    print(f"Filtering BAM file with criteria:", file=sys.stderr)
    print(f"  Minimum length: {args.min_length:,} bp", file=sys.stderr)
    print(f"  Minimum identity: {args.min_identity}%", file=sys.stderr)
    print(f"  Input: {args.input_bam}", file=sys.stderr)
    print(f"  Output: {args.output_bam}", file=sys.stderr)
    print(f"", file=sys.stderr)
    
    try:
        filter_bam(args.input_bam, args.output_bam, args.min_length, args.min_identity)
        
        # Index the output BAM file
        print("Indexing output BAM file...", file=sys.stderr)
        pysam.index(args.output_bam)
        print(f"Index created: {args.output_bam}.bai", file=sys.stderr)
        
    except FileNotFoundError as e:
        print(f"Error: File not found - {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
