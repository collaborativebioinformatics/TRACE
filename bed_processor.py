#!/usr/bin/env python3
"""
BED File Processor for Structural Variants

Processes BED files containing structural variants (SVs) and modifies coordinates based on SVTYPE:
- DEL: Extends start by -1000bp and end by +1000bp
- DUP/INV: Creates two rows with flanking regions around breakpoints

Authors: Michal Izydorczyk, Nicola Wong, John Adedeji, Julian Chiu, and Thomas X. Garcia
"""

import argparse
import sys
import os
import re
import logging
from typing import List, Tuple, Optional

def setup_logging(verbose: bool = False, debug: bool = False) -> None:
    """Configure logging settings."""
    if debug:
        level = logging.DEBUG
    elif verbose:
        level = logging.INFO
    else:
        level = logging.WARNING
    
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="""
        Process BED files containing structural variants and modify coordinates based on SVTYPE.
        
        For SVTYPE=DEL: Extends the region by 1000bp on each side.
        For SVTYPE=DUP or SVTYPE=INV: Creates two separate 2000bp regions around breakpoints.
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
        Example usage:
            %(prog)s input.bed -o output.bed
            %(prog)s input.bed -o output.bed --verbose
            %(prog)s input.bed -o output.bed --flank-size 500
        """
    )
    
    parser.add_argument(
        'input_file',
        help='Input BED file with structural variants'
    )
    
    parser.add_argument(
        '-o', '--output',
        dest='output_file',
        default=None,
        help='Output BED file (default: input_processed.bed)'
    )
    
    parser.add_argument(
        '-f', '--flank-size', '--flank-length',
        type=int,
        default=1000,
        dest='flank_size',
        help='Size/length of flanking regions in bp (default: 1000)'
    )
    
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose output'
    )
    
    parser.add_argument(
        '-d', '--debug',
        action='store_true',
        help='Enable debug output (very verbose)'
    )
    
    parser.add_argument(
        '--min-coord',
        type=int,
        default=0,
        help='Minimum coordinate value (default: 0)'
    )
    
    parser.add_argument(
        '--skip-invalid',
        action='store_true',
        help='Skip invalid lines instead of exiting with error'
    )
    
    return parser.parse_args()

def extract_svtype(info_field: str) -> Optional[str]:
    """
    Extract SVTYPE from the INFO field of a BED line.
    
    Args:
        info_field: The INFO column containing SVTYPE information
    
    Returns:
        The SVTYPE value (e.g., 'DEL', 'DUP', 'INV') or None if not found
    """
    svtype_match = re.search(r'SVTYPE=([^;]+)', info_field)
    if svtype_match:
        return svtype_match.group(1)
    return None

def process_deletion(chrom: str, start: int, end: int, fields: List[str], 
                    flank_size: int, min_coord: int) -> List[str]:
    """
    Process deletion variant by extending coordinates.
    
    Args:
        chrom: Chromosome name
        start: Start coordinate
        end: End coordinate
        fields: All fields from the original BED line
        flank_size: Size of flanking region
        min_coord: Minimum allowed coordinate
    
    Returns:
        List containing the modified BED line
    """
    new_start = max(min_coord, start - flank_size)
    new_end = end + flank_size
    
    modified_fields = fields.copy()
    modified_fields[1] = str(new_start)
    modified_fields[2] = str(new_end)
    
    logging.debug(f"DEL: {chrom}:{start}-{end} -> {chrom}:{new_start}-{new_end}")
    
    return ['\t'.join(modified_fields)]

def process_duplication_inversion(chrom: str, start: int, end: int, fields: List[str],
                                 flank_size: int, min_coord: int, svtype: str) -> List[str]:
    """
    Process duplication or inversion variant by creating two flanking regions.
    
    Args:
        chrom: Chromosome name
        start: Start coordinate
        end: End coordinate
        fields: All fields from the original BED line
        flank_size: Size of flanking region
        min_coord: Minimum allowed coordinate
        svtype: Type of SV (DUP or INV)
    
    Returns:
        List containing two modified BED lines
    """
    results = []
    
    # First row: region around start breakpoint
    row1_start = max(min_coord, start - flank_size)
    row1_end = start + flank_size
    
    row1_fields = fields.copy()
    row1_fields[1] = str(row1_start)
    row1_fields[2] = str(row1_end)
    
    # Keep original variant ID for proper annotation lookup
    # (removed ID modification that was breaking annotation)
    
    results.append('\t'.join(row1_fields))
    logging.debug(f"{svtype} row1: {chrom}:{row1_start}-{row1_end}")
    
    # Second row: region around end breakpoint
    row2_start = max(min_coord, end - flank_size)
    row2_end = end + flank_size
    
    row2_fields = fields.copy()
    row2_fields[1] = str(row2_start)
    row2_fields[2] = str(row2_end)
    
    # Keep original variant ID for proper annotation lookup
    # (removed ID modification that was breaking annotation)
    
    results.append('\t'.join(row2_fields))
    logging.debug(f"{svtype} row2: {chrom}:{row2_start}-{row2_end}")
    
    return results

def process_bed_line(line: str, flank_size: int, min_coord: int, 
                     line_num: int, skip_invalid: bool) -> List[str]:
    """
    Process a single BED line based on its SVTYPE.
    
    Args:
        line: Input BED line
        flank_size: Size of flanking region
        min_coord: Minimum allowed coordinate
        line_num: Line number for error reporting
        skip_invalid: Whether to skip invalid lines
    
    Returns:
        List of processed BED lines
    
    Raises:
        ValueError: If the line format is invalid and skip_invalid is False
    """
    line = line.strip()
    if not line or line.startswith('#'):
        return [line] if line else []
    
    fields = line.split('\t')
    
    # Validate minimum required fields
    if len(fields) < 8:
        error_msg = f"Line {line_num}: Invalid BED format - expected at least 8 columns, got {len(fields)}"
        if skip_invalid:
            logging.warning(error_msg)
            return []
        else:
            raise ValueError(error_msg)
    
    try:
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        # Check if INFO field is in column 8 (index 7) or column 9 (index 8)
        # This handles both standard BED and VCF-like formats
        if len(fields) > 8 and 'SVTYPE=' in fields[8]:
            info_field = fields[8]
        elif len(fields) > 7 and 'SVTYPE=' in fields[7]:
            info_field = fields[7]
        else:
            # Try to find SVTYPE in any field
            info_field = None
            for i in range(7, len(fields)):
                if 'SVTYPE=' in fields[i]:
                    info_field = fields[i]
                    break
            if not info_field:
                error_msg = f"Line {line_num}: No field containing SVTYPE found"
                if skip_invalid:
                    logging.warning(error_msg)
                    return [line]
                else:
                    logging.warning(error_msg + " - returning unchanged")
                    return [line]
    except (IndexError, ValueError) as e:
        error_msg = f"Line {line_num}: Error parsing coordinates - {e}"
        if skip_invalid:
            logging.warning(error_msg)
            return []
        else:
            raise ValueError(error_msg)
    
    # Extract SVTYPE from INFO field
    svtype = extract_svtype(info_field)
    
    if not svtype:
        error_msg = f"Line {line_num}: No SVTYPE found in INFO field"
        if skip_invalid:
            logging.warning(error_msg)
            return [line]  # Return unchanged
        else:
            logging.warning(error_msg + " - returning unchanged")
            return [line]
    
    logging.info(f"Processing line {line_num}: {chrom}:{start}-{end} SVTYPE={svtype}")
    
    # Process based on SVTYPE
    if svtype == 'DEL':
        return process_deletion(chrom, start, end, fields, flank_size, min_coord)
    elif svtype in ['DUP', 'INV']:
        return process_duplication_inversion(chrom, start, end, fields, flank_size, min_coord, svtype)
    else:
        logging.info(f"Line {line_num}: SVTYPE={svtype} not in [DEL, DUP, INV] - returning unchanged")
        return [line]

def process_bed_file(input_file: str, output_file: str, flank_size: int, 
                    min_coord: int, skip_invalid: bool) -> Tuple[int, int, int]:
    """
    Process entire BED file.
    
    Args:
        input_file: Path to input BED file
        output_file: Path to output BED file
        flank_size: Size of flanking region
        min_coord: Minimum allowed coordinate
        skip_invalid: Whether to skip invalid lines
    
    Returns:
        Tuple of (lines_read, lines_written, lines_skipped)
    
    Raises:
        FileNotFoundError: If input file doesn't exist
        IOError: If there are problems reading/writing files
    """
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")
    
    lines_read = 0
    lines_written = 0
    lines_skipped = 0
    
    try:
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            for line_num, line in enumerate(infile, 1):
                lines_read += 1
                
                try:
                    processed_lines = process_bed_line(
                        line, flank_size, min_coord, line_num, skip_invalid
                    )
                    
                    for processed_line in processed_lines:
                        if processed_line:  # Skip empty lines
                            outfile.write(processed_line + '\n')
                            lines_written += 1
                    
                    if not processed_lines:
                        lines_skipped += 1
                        
                except ValueError as e:
                    logging.error(f"Error processing line {line_num}: {e}")
                    raise
                except Exception as e:
                    logging.error(f"Unexpected error at line {line_num}: {e}")
                    raise
                    
    except IOError as e:
        logging.error(f"I/O error: {e}")
        raise
    
    return lines_read, lines_written, lines_skipped

def main():
    """Main function to orchestrate BED file processing."""
    args = parse_arguments()
    
    # Setup logging
    setup_logging(verbose=args.verbose, debug=args.debug)
    
    # Determine output file name
    if args.output_file is None:
        base_name = os.path.splitext(args.input_file)[0]
        args.output_file = f"{base_name}_processed.bed"
    
    logging.info(f"Input file: {args.input_file}")
    logging.info(f"Output file: {args.output_file}")
    logging.info(f"Flank size: {args.flank_size} bp")
    logging.info(f"Minimum coordinate: {args.min_coord}")
    
    try:
        # Process the BED file
        lines_read, lines_written, lines_skipped = process_bed_file(
            args.input_file,
            args.output_file,
            args.flank_size,
            args.min_coord,
            args.skip_invalid
        )
        
        # Print summary
        print(f"\nProcessing complete!")
        print(f"Lines read: {lines_read}")
        print(f"Lines written: {lines_written}")
        print(f"Lines skipped: {lines_skipped}")
        print(f"Output written to: {args.output_file}")
        
        if args.verbose:
            if lines_written > lines_read:
                print(f"Note: Output has more lines than input due to DUP/INV variants creating two rows each")
        
        return 0
        
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        print("Use --skip-invalid to skip problematic lines", file=sys.stderr)
        return 2
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        logging.exception("Unexpected error occurred")
        return 3

if __name__ == "__main__":
    sys.exit(main())