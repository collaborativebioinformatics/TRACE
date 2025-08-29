#!/usr/bin/env python3
"""
squash_intersect.py - Squash multiple rows per variant into single rows

This script processes BED files from Sniffles intersections and combines
multiple rows with the same SnifflesID into a single row, preserving
data from columns 11-19.

Authors: Michal Izydorczyk, Nicola Wong, John Adedeji, Julian Chiu, and Thomas X. Garcia
"""

import argparse
import sys
import os
import logging
from collections import defaultdict
from typing import Dict, List, Tuple


def setup_logging(debug: bool) -> None:
    """
    Set up logging configuration.
    
    Args:
        debug: If True, set logging level to DEBUG, otherwise INFO
    """
    level = logging.DEBUG if debug else logging.INFO
    format_str = '%(asctime)s - %(levelname)s - %(message)s' if debug else '%(levelname)s: %(message)s'
    
    logging.basicConfig(
        level=level,
        format=format_str,
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def parse_arguments() -> argparse.Namespace:
    """
    Parse command line arguments.
    
    Returns:
        Parsed arguments namespace
    """
    parser = argparse.ArgumentParser(
        description='Squash multiple BED file rows with the same SnifflesID into single rows',
        epilog="""
Example usage:
  %(prog)s input.bed                # outputs to input_squashed.tsv
  %(prog)s input.bed -o output.tsv  # outputs to output.tsv
  %(prog)s input.bed -o -           # outputs to stdout
  %(prog)s input.bed --debug        # outputs to input_squashed.tsv with debug logging
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        'input_file',
        help='Input BED file with Sniffles data'
    )
    
    parser.add_argument(
        '-o', '--output',
        help='Output TSV file (default: input_squashed.tsv, use "-" for stdout)',
        default=None
    )
    
    parser.add_argument(
        '--debug',
        action='store_true',
        help='Enable debug logging'
    )
    
    parser.add_argument(
        '--skip-header',
        action='store_true',
        help='Skip header line in input file'
    )
    
    parser.add_argument(
        '--no-header',
        action='store_true',
        help='Do not write header line to output'
    )
    
    parser.add_argument(
        '--delimiter',
        default='\t',
        help='Input file delimiter (default: tab)'
    )
    
    return parser.parse_args()


def validate_input_file(filepath: str) -> None:
    """
    Validate that the input file exists and is readable.
    
    Args:
        filepath: Path to the input file
        
    Raises:
        FileNotFoundError: If file doesn't exist
        PermissionError: If file isn't readable
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Input file not found: {filepath}")
    
    if not os.path.isfile(filepath):
        raise ValueError(f"Input path is not a file: {filepath}")
    
    if not os.access(filepath, os.R_OK):
        raise PermissionError(f"Input file is not readable: {filepath}")
    
    logging.debug(f"Input file validated: {filepath}")


def process_bed_line(line: str, delimiter: str, line_num: int) -> Tuple[str, List[str]]:
    """
    Process a single line from the BED file.
    
    Args:
        line: Input line
        delimiter: Column delimiter
        line_num: Line number for error reporting
        
    Returns:
        Tuple of (SnifflesID, list of columns 11-19)
        
    Raises:
        ValueError: If line doesn't have enough columns
    """
    line = line.strip()
    
    if not line or line.startswith('#'):
        return None, None
    
    columns = line.split(delimiter)
    
    # BED files are 0-indexed, but we need at least 19 columns (0-18)
    if len(columns) < 19:
        logging.warning(f"Line {line_num}: Insufficient columns ({len(columns)} < 19), skipping")
        return None, None
    
    # SnifflesID is in column 4 (0-indexed column 3)
    sniffles_id = columns[3]
    
    # Extract columns 11-19 (0-indexed columns 10-18)
    elements = columns[10:19]
    
    logging.debug(f"Line {line_num}: SnifflesID={sniffles_id}, Elements={' '.join(elements)}")
    
    return sniffles_id, elements


def squash_variants(input_file: str, delimiter: str, skip_header: bool) -> Dict[str, List[List[str]]]:
    """
    Read the input file and group rows by SnifflesID.
    
    Args:
        input_file: Path to input BED file
        delimiter: Column delimiter
        skip_header: Whether to skip the first line
        
    Returns:
        Dictionary mapping SnifflesID to list of element lists
    """
    variants = defaultdict(list)
    line_count = 0
    skipped_count = 0
    
    try:
        with open(input_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                if line_num == 1 and skip_header:
                    logging.debug("Skipping header line")
                    continue
                
                sniffles_id, elements = process_bed_line(line, delimiter, line_num)
                
                if sniffles_id is None:
                    skipped_count += 1
                    continue
                
                variants[sniffles_id].append(elements)
                line_count += 1
                
                if line_count % 10000 == 0:
                    logging.info(f"Processed {line_count} valid lines...")
    
    except Exception as e:
        logging.error(f"Error reading input file: {e}")
        raise
    
    logging.info(f"Processed {line_count} valid lines, skipped {skipped_count} lines")
    logging.info(f"Found {len(variants)} unique SnifflesIDs")
    
    # Log statistics about duplicates
    dup_counts = [len(v) for v in variants.values() if len(v) > 1]
    if dup_counts:
        logging.info(f"SnifflesIDs with multiple rows: {len(dup_counts)}")
        logging.info(f"Max rows per SnifflesID: {max(dup_counts)}")
        logging.info(f"Average rows per duplicated SnifflesID: {sum(dup_counts)/len(dup_counts):.2f}")
    
    return variants


def format_output_line(sniffles_id: str, element_lists: List[List[str]]) -> str:
    """
    Format the output line for a single SnifflesID.
    
    Args:
        sniffles_id: The SnifflesID
        element_lists: List of element lists from different rows
        
    Returns:
        Formatted TSV line
    """
    # Join elements within each row with commas
    formatted_rows = [','.join(elements) for elements in element_lists]
    
    # Join multiple rows with pipes
    combined_elements = '|'.join(formatted_rows)
    
    return f"{sniffles_id}\t{combined_elements}"


def write_output(variants: Dict[str, List[List[str]]], output_file: str, no_header: bool) -> None:
    """
    Write the squashed variants to output.
    
    Args:
        variants: Dictionary of variants
        output_file: Output file path (None for stdout)
        no_header: Whether to skip writing header
    """
    output_handle = sys.stdout if output_file is None else open(output_file, 'w')
    
    try:
        # Write header if requested
        if not no_header:
            output_handle.write("SnifflesID\tElements\n")
        
        # Write data sorted by SnifflesID for consistency
        for sniffles_id in sorted(variants.keys()):
            line = format_output_line(sniffles_id, variants[sniffles_id])
            output_handle.write(line + '\n')
        
        logging.info(f"Wrote {len(variants)} variants to output")
        
    finally:
        if output_file is not None:
            output_handle.close()
            logging.info(f"Output written to: {output_file}")


def main():
    """Main function to coordinate the script execution."""
    try:
        # Parse arguments
        args = parse_arguments()
        
        # Set up logging
        setup_logging(args.debug)
        
        logging.debug(f"Arguments: {args}")
        
        # Validate input
        validate_input_file(args.input_file)
        
        # Determine output file
        output_file = args.output
        if output_file is None:
            # Default: replace .bed with _squashed.tsv
            if args.input_file.endswith('.bed'):
                output_file = args.input_file[:-4] + '_squashed.tsv'
            else:
                output_file = args.input_file + '_squashed.tsv'
            logging.info(f"No output file specified, using: {output_file}")
        elif output_file == '-':
            # User explicitly requested stdout
            output_file = None
        
        # Process the input file
        logging.info(f"Processing input file: {args.input_file}")
        variants = squash_variants(args.input_file, args.delimiter, args.skip_header)
        
        if not variants:
            logging.warning("No valid data found in input file")
            sys.exit(1)
        
        # Write output
        write_output(variants, output_file, args.no_header)
        
        logging.info("Processing complete")
        
    except KeyboardInterrupt:
        logging.error("Process interrupted by user")
        sys.exit(130)
    except FileNotFoundError as e:
        logging.error(f"File error: {e}")
        sys.exit(1)
    except PermissionError as e:
        logging.error(f"Permission error: {e}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Unexpected error: {e}")
        if args.debug:
            logging.exception("Full traceback:")
        sys.exit(1)


if __name__ == "__main__":
    main()