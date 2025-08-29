#!/usr/bin/env python3
"""
VCF to FASTA converter for insertion structural variants (INS only)

This script extracts only insertion sequences from Sniffles or CuteSV VCF files
and outputs them in FASTA format.

Authors: Michal Izydorczyk, Nicola Wong, John Adedeji, Julian Chiu, and Thomas X. Garcia
License: MIT
"""

import argparse
import sys
import os
import re
import gzip
from pathlib import Path
from typing import TextIO, Tuple, Optional
import logging


class VCFParseError(Exception):
    """Custom exception for VCF parsing errors"""
    pass


def setup_logging(verbose: bool = False) -> None:
    """
    Set up logging configuration
    
    Args:
        verbose: Enable verbose logging if True
    """
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s',
        level=level,
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def validate_vcf_file(filepath: Path) -> None:
    """
    Validate that the file exists and appears to be a VCF file
    
    Args:
        filepath: Path to the VCF file
        
    Raises:
        FileNotFoundError: If file doesn't exist
        VCFParseError: If file doesn't appear to be a valid VCF
    """
    if not filepath.exists():
        raise FileNotFoundError(f"Input file not found: {filepath}")
    
    if not filepath.is_file():
        raise VCFParseError(f"Input path is not a file: {filepath}")
    
    # Check if file is readable
    if not os.access(filepath, os.R_OK):
        raise PermissionError(f"Cannot read file: {filepath}")
    
    # Check file extension
    valid_extensions = ['.vcf', '.vcf.gz']
    if not any(str(filepath).endswith(ext) for ext in valid_extensions):
        logging.warning(f"File doesn't have typical VCF extension: {filepath}")


def open_vcf_file(filepath: Path) -> TextIO:
    """
    Open a VCF file, handling gzip compression if necessary
    
    Args:
        filepath: Path to the VCF file
        
    Returns:
        File handle for reading
    """
    if str(filepath).endswith('.gz'):
        return gzip.open(filepath, 'rt', encoding='utf-8')
    else:
        return open(filepath, 'r', encoding='utf-8')


def parse_info_field(info_str: str) -> dict:
    """
    Parse the INFO field of a VCF record
    
    Args:
        info_str: INFO field string from VCF
        
    Returns:
        Dictionary of INFO field key-value pairs
    """
    info_dict = {}
    for item in info_str.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info_dict[key] = value
        else:
            info_dict[item] = True
    return info_dict


def extract_sv_type(info_dict: dict) -> Optional[str]:
    """
    Extract the SV type from INFO field
    
    Args:
        info_dict: Parsed INFO field dictionary
        
    Returns:
        SV type (INS) or None if not found or not INS
    """
    if 'SVTYPE' in info_dict:
        return info_dict['SVTYPE']
    return None


def clean_sequence(seq: str) -> str:
    """
    Clean sequence by removing non-DNA characters
    
    Args:
        seq: Input sequence
        
    Returns:
        Cleaned sequence containing only A, T, G, C, N
    """
    # Remove any characters that aren't standard DNA bases
    cleaned = re.sub(r'[^ATGCNatgcn]', '', seq)
    return cleaned.upper()


def format_fasta_sequence(sequence: str, line_length: int = 80) -> str:
    """
    Format a sequence for FASTA output with proper line wrapping
    
    Args:
        sequence: DNA sequence
        line_length: Maximum characters per line
        
    Returns:
        Formatted sequence with line breaks
    """
    lines = []
    for i in range(0, len(sequence), line_length):
        lines.append(sequence[i:i + line_length])
    return '\n'.join(lines)


def process_vcf_record(
    chrom: str, 
    pos: str, 
    var_id: str, 
    ref: str, 
    alt: str, 
    info: str
) -> Tuple[Optional[str], Optional[str]]:
    """
    Process a single VCF record to extract insertion information
    
    Args:
        chrom: Chromosome
        pos: Position
        var_id: Variant ID
        ref: Reference allele
        alt: Alternative allele
        info: INFO field
        
    Returns:
        Tuple of (header, sequence) or (None, None) if not INS
    """
    try:
        # Parse INFO field
        info_dict = parse_info_field(info)
        sv_type = extract_sv_type(info_dict)
        
        # Only process INS variants
        if sv_type != 'INS':
            return None, None
        
        # Skip symbolic insertions
        if alt.startswith('<') and alt.endswith('>'):
            logging.debug(f"Skipping symbolic insertion {var_id}: {alt}")
            return None, None
        
        # Create header with variant ID and additional info
        header_parts = [var_id if var_id != '.' else f"{chrom}_{pos}_INS"]
        header_parts.append(f"CHR={chrom}")
        header_parts.append(f"POS={pos}")
        
        # Add SV length if available
        if 'SVLEN' in info_dict:
            header_parts.append(f"LEN={info_dict['SVLEN']}")
        
        header = ' '.join(header_parts)
        
        # Extract insertion sequence from ALT
        # Remove any bracket notations or special characters
        alt_clean = re.sub(r'[\[\]<>].*', '', alt)
        sequence = clean_sequence(alt_clean)
        
        # Skip if sequence is too short or empty
        if len(sequence) <= 1:
            logging.debug(f"Skipping variant {var_id}: sequence too short")
            return None, None
        
        return header, sequence
        
    except Exception as e:
        logging.warning(f"Error processing record at {chrom}:{pos}: {e}")
        return None, None


def process_vcf_file(
    input_file: Path, 
    output_file: Path, 
    min_length: int = 1,
    max_length: Optional[int] = None
) -> Tuple[int, int]:
    """
    Process VCF file and extract INS variants to FASTA
    
    Args:
        input_file: Path to input VCF file
        output_file: Path to output FASTA file
        min_length: Minimum sequence length to include
        max_length: Maximum sequence length to include (None for no limit)
        
    Returns:
        Tuple of (insertions_extracted, symbolic_insertions_skipped)
    """
    validate_vcf_file(input_file)
    
    insertions_extracted = 0
    symbolic_skipped = 0
    processed_lines = 0
    
    logging.info(f"Processing VCF file: {input_file}")
    logging.info(f"Extracting INS variants only")
    logging.info(f"Output will be written to: {output_file}")
    
    try:
        with open_vcf_file(input_file) as vcf_handle, \
             open(output_file, 'w') as fasta_handle:
            
            for line_num, line in enumerate(vcf_handle, 1):
                line = line.strip()
                
                # Skip empty lines
                if not line:
                    continue
                
                # Skip header lines
                if line.startswith('#'):
                    continue
                
                processed_lines += 1
                
                # Parse VCF record
                fields = line.split('\t')
                if len(fields) < 8:
                    logging.warning(f"Line {line_num}: Invalid VCF format (expected at least 8 fields)")
                    continue
                
                chrom = fields[0]
                pos = fields[1]
                var_id = fields[2]
                ref = fields[3]
                alt = fields[4]
                info = fields[7]
                
                # Count symbolic insertions
                if 'SVTYPE=INS' in info and alt.startswith('<') and alt.endswith('>'):
                    symbolic_skipped += 1
                
                # Process the record
                header, sequence = process_vcf_record(
                    chrom, pos, var_id, ref, alt, info
                )
                
                if header and sequence:
                    # Apply length filters
                    seq_len = len(sequence)
                    if seq_len < min_length:
                        continue
                    if max_length and seq_len > max_length:
                        continue
                    
                    # Write to FASTA
                    fasta_handle.write(f">{header}\n")
                    fasta_handle.write(f"{format_fasta_sequence(sequence)}\n")
                    
                    insertions_extracted += 1
                
                # Progress indicator for large files
                if processed_lines % 10000 == 0:
                    logging.info(f"Processed {processed_lines} variants...")
    
    except Exception as e:
        logging.error(f"Error processing VCF file: {e}")
        raise
    
    return insertions_extracted, symbolic_skipped


def generate_output_filename(input_file: Path, output_file: Optional[Path]) -> Path:
    """
    Generate output filename based on input file if not specified
    
    Args:
        input_file: Path to input VCF file
        output_file: User-specified output file path (optional)
        
    Returns:
        Path to output file
    """
    if output_file:
        return output_file
    
    # Remove .vcf or .vcf.gz extension and add _INS.fasta
    input_str = str(input_file)
    if input_str.endswith('.vcf.gz'):
        base = input_str[:-7]
    elif input_str.endswith('.vcf'):
        base = input_str[:-4]
    else:
        base = input_str
    
    return Path(f"{base}_INS.fasta")


def main():
    """Main function to handle command line arguments and run the conversion"""
    
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description='Extract INS (insertion) structural variants from VCF files to FASTA format',
        epilog='''
This script extracts only insertion variants (SVTYPE=INS) from VCF files.
Symbolic insertions (<INS>) are skipped as they don't contain sequence data.

Examples:
  %(prog)s input.vcf
  %(prog)s input.vcf -o output_insertions.fasta
  %(prog)s input.vcf.gz --min-length 50 --max-length 1000
  %(prog)s sniffles_output.vcf --verbose
        ''',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Positional arguments
    parser.add_argument(
        'input',
        type=Path,
        help='Input VCF file from Sniffles or CuteSV (can be gzipped)'
    )
    
    # Optional arguments
    parser.add_argument(
        '-o', '--output',
        type=Path,
        help='Output FASTA file (default: input_INS.fasta)'
    )
    
    parser.add_argument(
        '--min-length',
        type=int,
        default=1,
        help='Minimum sequence length to include (default: 1)'
    )
    
    parser.add_argument(
        '--max-length',
        type=int,
        default=None,
        help='Maximum sequence length to include (default: no limit)'
    )
    
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose output'
    )
    
    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s 1.0.0'
    )
    
    # Parse arguments
    args = parser.parse_args()
    
    # Set up logging
    setup_logging(args.verbose)
    
    # Generate output filename if not specified
    output_file = generate_output_filename(args.input, args.output)
    
    # Validate arguments
    if args.min_length < 1:
        parser.error("Minimum length must be at least 1")
    
    if args.max_length and args.max_length < args.min_length:
        parser.error("Maximum length must be greater than minimum length")
    
    try:
        # Process the VCF file
        ins_count, symbolic_count = process_vcf_file(
            args.input,
            output_file,
            args.min_length,
            args.max_length
        )
        
        # Print summary
        print(f"\nProcessing complete!")
        print(f"Insertions extracted: {ins_count}")
        if symbolic_count > 0:
            print(f"Symbolic insertions skipped: {symbolic_count}")
            print(f"  (Symbolic insertions like <INS> don't contain sequence data)")
        print(f"Output written to: {output_file}")
        
    except FileNotFoundError as e:
        logging.error(f"File not found: {e}")
        sys.exit(1)
    except PermissionError as e:
        logging.error(f"Permission denied: {e}")
        sys.exit(1)
    except VCFParseError as e:
        logging.error(f"VCF parsing error: {e}")
        sys.exit(1)
    except KeyboardInterrupt:
        logging.info("\nProcessing interrupted by user")
        sys.exit(130)
    except Exception as e:
        logging.error(f"Unexpected error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()