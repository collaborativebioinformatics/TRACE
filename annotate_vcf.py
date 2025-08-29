#!/usr/bin/env python3
"""
VCF Annotation Script - Adds transposable element annotations to VCF INFO field.

This script reads a VCF file and adds RM_TE annotations to the INFO column
from external TSV files by matching SnifflesID values. The annotations are
added as a new INFO tag rather than a separate column.

Authors: Michal Izydorczyk, Nicola Wong, John Adedeji, Julian Chiu, and Thomas X. Garcia
"""

import sys
import argparse
import os
from pathlib import Path
import logging
from typing import Dict, List, Tuple, Optional


def setup_logging(debug: bool = False) -> None:
    """
    Configure logging based on debug flag.
    
    Args:
        debug: If True, set logging to DEBUG level, otherwise INFO
    """
    level = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def generate_output_filename(input_file: str) -> str:
    """
    Generate output filename by replacing .vcf with _annotated.vcf.
    
    Args:
        input_file: Path to input VCF file
    
    Returns:
        Path to output file with _annotated.vcf suffix
    """
    input_path = Path(input_file)
    
    # Check if file ends with .vcf (case-insensitive)
    if input_path.suffix.lower() == '.vcf':
        # Replace .vcf with _annotated.vcf
        output_path = input_path.with_suffix('').with_name(
            input_path.stem + '_annotated.vcf'
        )
    else:
        # If it doesn't end with .vcf, just append _annotated.vcf
        output_path = input_path.with_name(
            input_path.name + '_annotated.vcf'
        )
    
    return str(output_path)


def parse_arguments() -> argparse.Namespace:
    """
    Parse command line arguments.
    
    Returns:
        Parsed arguments namespace
    """
    parser = argparse.ArgumentParser(
        description='Annotate VCF file INFO field with transposable element data from external TSV files.',
        epilog='''
Example usage:
    %(prog)s input.vcf --data data1.tsv data2.tsv
        # Creates input_annotated.vcf with RM_TE tag in INFO field
    
    %(prog)s input.vcf --data data1.tsv data2.tsv --output custom.vcf
        # Uses custom output filename
    
    %(prog)s input.vcf --data data1.tsv data2.tsv --column-name REPEAT_ANNOTATION
        # Uses custom INFO tag name
    
    %(prog)s input.vcf --data data1.tsv data2.tsv --stdout
        # Outputs to stdout instead of file
    
    %(prog)s input.vcf --data data1.tsv data2.tsv --debug
        # Enable debug logging
        ''',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        'input_file',
        type=str,
        help='Input tab-separated file (VCF format) with 10 columns'
    )
    
    parser.add_argument(
        '--data',
        type=str,
        nargs=2,
        required=True,
        metavar=('FILE1', 'FILE2'),
        help='Two data files containing SnifflesID and annotation data'
    )
    
    parser.add_argument(
        '--output',
        '-o',
        type=str,
        default=None,
        help='Override output file path (default: input_annotated.vcf)'
    )
    
    parser.add_argument(
        '--stdout',
        action='store_true',
        help='Output to stdout instead of file'
    )
    
    parser.add_argument(
        '--debug',
        action='store_true',
        help='Enable debug mode with verbose logging'
    )
    
    parser.add_argument(
        '--column-name',
        type=str,
        default='RM_TE',
        help='Name for the INFO tag (default: RM_TE)'
    )
    
    parser.add_argument(
        '--na-value',
        type=str,
        default='-',
        help='Value to use when no match is found (default: -)'
    )
    
    parser.add_argument(
        '--overwrite',
        action='store_true',
        help='Overwrite output file if it exists'
    )
    
    return parser.parse_args()


def validate_files(input_file: str, data_files: List[str]) -> None:
    """
    Validate that all input files exist and are readable.
    
    Args:
        input_file: Path to input VCF file
        data_files: List of paths to data files
    
    Raises:
        FileNotFoundError: If any file doesn't exist
        PermissionError: If any file is not readable
    """
    all_files = [input_file] + data_files
    
    for file_path in all_files:
        path = Path(file_path)
        if not path.exists():
            raise FileNotFoundError(f"File not found: {file_path}")
        if not path.is_file():
            raise ValueError(f"Not a regular file: {file_path}")
        if not os.access(file_path, os.R_OK):
            raise PermissionError(f"File is not readable: {file_path}")
    
    logging.debug(f"All input files validated successfully")


def validate_output_file(output_file: str, overwrite: bool) -> None:
    """
    Validate output file path.
    
    Args:
        output_file: Path to output file
        overwrite: Whether to overwrite existing file
    
    Raises:
        FileExistsError: If file exists and overwrite is False
        PermissionError: If directory is not writable
    """
    output_path = Path(output_file)
    
    # Check if file exists
    if output_path.exists() and not overwrite:
        raise FileExistsError(
            f"Output file already exists: {output_file}\n"
            f"Use --overwrite to replace it or specify a different output file"
        )
    
    # Check if parent directory is writable
    parent_dir = output_path.parent
    if not parent_dir.exists():
        raise FileNotFoundError(f"Output directory does not exist: {parent_dir}")
    if not os.access(parent_dir, os.W_OK):
        raise PermissionError(f"Output directory is not writable: {parent_dir}")
    
    logging.debug(f"Output file path validated: {output_file}")


def read_data_files(data_files: List[str]) -> Dict[str, str]:
    """
    Read data files and create a lookup dictionary.
    
    Args:
        data_files: List of paths to data files (TSV format)
    
    Returns:
        Dictionary mapping SnifflesID to annotation data
    
    Raises:
        ValueError: If file format is invalid
    """
    lookup_dict = {}
    
    for file_path in data_files:
        logging.info(f"Reading data file: {file_path}")
        
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                # Skip header line if present
                header = f.readline().strip()
                if not header:
                    raise ValueError(f"Empty data file: {file_path}")
                
                logging.debug(f"Data file header: {header}")
                
                line_num = 1
                for line in f:
                    line_num += 1
                    line = line.strip()
                    
                    if not line:  # Skip empty lines
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) < 2:
                        logging.warning(
                            f"Line {line_num} in {file_path} has fewer than 2 columns, skipping"
                        )
                        continue
                    
                    sniffles_id = parts[0]
                    annotation_data = '\t'.join(parts[1:])  # Join all columns after first
                    
                    if sniffles_id in lookup_dict:
                        logging.warning(
                            f"Duplicate SnifflesID '{sniffles_id}' found in {file_path}, "
                            f"overwriting previous value"
                        )
                    
                    lookup_dict[sniffles_id] = annotation_data
                    logging.debug(f"Added: {sniffles_id} -> {annotation_data[:50]}...")
        
        except Exception as e:
            raise ValueError(f"Error reading data file {file_path}: {str(e)}")
    
    logging.info(f"Loaded {len(lookup_dict)} annotations from data files")
    return lookup_dict


def process_vcf_file(
    input_file: str,
    output_file: Optional[str],
    lookup_dict: Dict[str, str],
    column_name: str,
    na_value: str,
    use_stdout: bool = False
) -> Tuple[int, int]:
    """
    Process the VCF file and add annotation to INFO column.
    
    Args:
        input_file: Path to input VCF file
        output_file: Path to output file (None for auto-generated)
        lookup_dict: Dictionary mapping SnifflesID to annotation data
        column_name: Name for the INFO tag (default: RM_TE)
        na_value: Value to use when no match is found
        use_stdout: If True, output to stdout instead of file
    
    Returns:
        Tuple of (total_records, annotated_records)
    
    Raises:
        ValueError: If file format is invalid
    """
    total_records = 0
    annotated_records = 0
    header_written = False
    
    # Determine output destination
    if use_stdout:
        output_handle = sys.stdout
        logging.info("Writing output to stdout")
        actual_output = "stdout"
    else:
        if output_file is None:
            output_file = generate_output_filename(input_file)
        output_handle = open(output_file, 'w', encoding='utf-8')
        logging.info(f"Writing output to: {output_file}")
        actual_output = output_file
    
    try:
        with open(input_file, 'r', encoding='utf-8') as f:
            line_num = 0
            last_info_line = False
            for line in f:
                line_num += 1
                line = line.rstrip('\n')
                
                # Handle header lines (starting with ##)
                if line.startswith('##'):
                    # Check if this is an INFO line
                    if line.startswith('##INFO'):
                        last_info_line = True
                        output_handle.write(line + '\n')
                    else:
                        # If we just finished INFO lines and haven't written header yet
                        if last_info_line and not header_written:
                            output_handle.write(f'##INFO=<ID={column_name},Number=1,Type=String,Description="Transposable element annotation from RepeatMasker and genomic intersections">\n')
                            header_written = True
                            last_info_line = False
                        output_handle.write(line + '\n')
                    continue
                
                # Handle column header line (starting with #CHROM)
                if line.startswith('#CHROM') or (line.startswith('#') and not line.startswith('##')):
                    # If we haven't written the header yet, write it now
                    if not header_written:
                        output_handle.write(f'##INFO=<ID={column_name},Number=1,Type=String,Description="Transposable element annotation from RepeatMasker and genomic intersections">\n')
                        header_written = True
                    # Don't add new column, keep original header
                    output_handle.write(line + '\n')
                    logging.debug(f"Processing column headers")
                    continue
                
                # Process data lines
                if not line:  # Skip empty lines
                    continue
                
                parts = line.split('\t')
                
                # Validate that line has expected number of columns (VCF should have at least 8)
                if len(parts) < 8:
                    logging.warning(
                        f"Line {line_num} has fewer than 8 columns, skipping: {line[:50]}..."
                    )
                    output_handle.write(line + '\n')
                    continue
                
                # Extract SnifflesID from 3rd column (index 2)
                sniffles_id = parts[2]
                total_records += 1
                
                # Look up annotation data
                if sniffles_id in lookup_dict:
                    annotation = lookup_dict[sniffles_id]
                    annotated_records += 1
                    logging.debug(f"Found annotation for {sniffles_id}")
                    
                    # Add annotation to INFO field (column 8, index 7)
                    info_field = parts[7]
                    if info_field == '.':
                        # Empty INFO field, replace with our annotation
                        parts[7] = f'{column_name}={annotation}'
                    else:
                        # Append to existing INFO field
                        parts[7] = f'{info_field};{column_name}={annotation}'
                    
                    # Write modified line
                    output_handle.write('\t'.join(parts) + '\n')
                else:
                    # No annotation found, write original line
                    logging.debug(f"No annotation found for {sniffles_id}")
                    output_handle.write(line + '\n')
        
        logging.info(f"Processed {total_records} records, annotated {annotated_records}")
        
    except Exception as e:
        # If we were writing to a file and an error occurred, try to clean up
        if not use_stdout and output_file:
            try:
                output_handle.close()
                os.remove(output_file)
                logging.debug(f"Removed incomplete output file: {output_file}")
            except:
                pass
        raise ValueError(f"Error processing VCF file: {str(e)}")
    
    finally:
        if not use_stdout:
            output_handle.close()
    
    return total_records, annotated_records, actual_output


def main():
    """Main function to orchestrate the annotation process."""
    try:
        # Parse arguments
        args = parse_arguments()
        
        # Setup logging
        setup_logging(args.debug)
        
        logging.info("Starting VCF annotation process")
        logging.debug(f"Arguments: {args}")
        
        # Validate input files
        validate_files(args.input_file, args.data)
        
        # Determine output file
        if args.stdout:
            output_file = None
            use_stdout = True
        else:
            if args.output:
                output_file = args.output
            else:
                output_file = generate_output_filename(args.input_file)
            use_stdout = False
            
            # Validate output file
            validate_output_file(output_file, args.overwrite)
        
        # Read data files and create lookup dictionary
        lookup_dict = read_data_files(args.data)
        
        # Process VCF file
        total_records, annotated_records, actual_output = process_vcf_file(
            args.input_file,
            output_file,
            lookup_dict,
            args.column_name,
            args.na_value,
            use_stdout
        )
        
        # Print summary statistics
        if use_stdout:
            # Print to stderr if outputting to stdout
            print(f"\nSummary:", file=sys.stderr)
            print(f"  Input file: {args.input_file}", file=sys.stderr)
            print(f"  Total records processed: {total_records}", file=sys.stderr)
            print(f"  Records annotated: {annotated_records}", file=sys.stderr)
            print(f"  Records without annotation: {total_records - annotated_records}", file=sys.stderr)
        else:
            # Print to stdout if outputting to file
            print(f"\nAnnotation completed successfully!")
            print(f"  Input file: {args.input_file}")
            print(f"  Output file: {actual_output}")
            print(f"  Total records processed: {total_records}")
            print(f"  Records annotated: {annotated_records}")
            print(f"  Records without annotation: {total_records - annotated_records}")
            
            logging.info(f"Output written to: {actual_output}")
        
        logging.info("VCF annotation completed successfully")
        
    except FileNotFoundError as e:
        logging.error(f"File error: {e}")
        sys.exit(1)
    except FileExistsError as e:
        logging.error(f"File exists: {e}")
        sys.exit(1)
    except PermissionError as e:
        logging.error(f"Permission error: {e}")
        sys.exit(1)
    except ValueError as e:
        logging.error(f"Value error: {e}")
        sys.exit(1)
    except KeyboardInterrupt:
        logging.info("Process interrupted by user")
        sys.exit(130)
    except Exception as e:
        logging.error(f"Unexpected error: {e}")
        if logging.getLogger().level == logging.DEBUG:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()