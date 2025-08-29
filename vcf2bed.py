#!/usr/bin/env python3
"""
Convert structural variant VCF file to BED format with enhanced BND parsing

Authors: Michal Izydorczyk, Nicola Wong, John Adedeji, Julian Chiu, and Thomas X. Garcia
"""

import argparse
import sys
import re
import os

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Convert structural variant VCF file to BED format with enhanced BND parsing',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Example usage:
  %(prog)s input.vcf                    # Creates input.bed in same directory
  %(prog)s input.vcf -o output.bed      # Creates output.bed
  %(prog)s /path/to/input.vcf           # Creates /path/to/input.bed
  %(prog)s input.vcf --no-BND-parse     # Skip BND ALT parsing
  %(prog)s input.vcf --include-INS      # Include INS variants in output

Output BED format:
  Column 1: Chromosome
  Column 2: Start position (0-based)
  Column 3: End position
  Column 4: Variant ID
  Column 5: REF allele (from VCF column 4)
  Column 6: ALT allele (from VCF column 5)
  Column 7: Quality score (from VCF column 6)
  Column 8: FILTER (from VCF column 7)
  Column 9: INFO (from VCF column 8)
  Column 10: FORMAT (from VCF column 9)
  
For BND variants with ALT parsing enabled (default):
  - Creates duplicate lines using breakpoint information from ALT column
  - ALT formats supported:
    * N]chr:pos] or N[chr:pos[
    * ]chr:pos]N or [chr:pos[N
    * REF]chr:pos] or REF[chr:pos[
    * ]chr:pos]REF or [chr:pos[REF
  - Duplicate line uses chr and pos from ALT field
        '''
    )
    
    parser.add_argument('input_vcf', 
                        help='Input VCF file containing structural variants')
    
    parser.add_argument('-o', '--output',
                        help='Output BED file (default: same path as input with .bed extension)',
                        default=None)
    
    parser.add_argument('--no-BND-parse',
                        action='store_true',
                        help='Do not parse BND ALT column for duplicate entries (default: parse BND)')
    
    parser.add_argument('--include-INS',
                        action='store_true',
                        help='Include INS variants in output (default: exclude INS)')
    
    return parser.parse_args()

def extract_info_field(info_string, field_name):
    """Extract a specific field value from VCF INFO column"""
    pattern = f'{field_name}=([^;]+)'
    match = re.search(pattern, info_string)
    return match.group(1) if match else None

def parse_bnd_alt(alt_field):
    """
    Parse BND ALT field to extract chromosome and position
    ALT formats supported:
    - N]chr:pos] or N[chr:pos[
    - ]chr:pos]N or [chr:pos[N  
    - REF]chr:pos] or REF[chr:pos[
    - ]chr:pos]REF or [chr:pos[REF
    Returns: (chromosome, position) or (None, None) if parsing fails
    """
    # Updated comprehensive patterns for all BND formats
    patterns = [
        # Pattern 1: Starts with nucleotide(s) followed by bracket notation
        # Matches: N]chr:pos], N[chr:pos[, T]chr:pos], G[chr:pos[, etc.
        r'^[ACGTN]+[\[\]]([^:]+):(\d+)[\[\]]$',
        
        # Pattern 2: Starts with bracket notation, ends with nucleotide(s)
        # Matches: ]chr:pos]N, [chr:pos[N, ]chr:pos]G, [chr:pos[T, etc.
        r'^[\[\]]([^:]+):(\d+)[\[\]][ACGTN]+$',
        
        # Pattern 3: Only bracket notation (no nucleotides on either side)
        # This might occur in some special cases
        r'^[\[\]]([^:]+):(\d+)[\[\]]$'
    ]
    
    for pattern in patterns:
        match = re.search(pattern, alt_field)
        if match:
            chrom = match.group(1)
            pos = int(match.group(2))
            return chrom, pos
    
    return None, None

def vcf_to_bed(input_file, output_file, parse_bnd=True, include_ins=False):
    """Convert VCF file to BED format with enhanced BND parsing"""
    
    try:
        with open(input_file, 'r') as vcf, open(output_file, 'w') as bed:
            line_count = 0
            variant_count = 0
            skipped_ins_count = 0
            parsed_bnd_count = 0
            unparsed_bnd_count = 0
            
            for line in vcf:
                line_count += 1
                
                # Skip header lines
                if line.startswith('#'):
                    continue
                
                # Parse VCF fields
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    print(f"Warning: Line {line_count} has insufficient fields, skipping", file=sys.stderr)
                    continue
                
                chrom = fields[0]
                pos = int(fields[1])
                var_id = fields[2]
                ref = fields[3]
                alt = fields[4]
                qual = fields[5]
                filter_field = fields[6]
                info = fields[7]
                format_field = fields[8] if len(fields) >= 9 else ''
                
                # Extract SVTYPE from INFO field
                svtype = extract_info_field(info, 'SVTYPE')
                
                # Skip INS variants if not included
                if svtype == 'INS' and not include_ins:
                    skipped_ins_count += 1
                    continue
                
                # Extract END and SVLEN from INFO field
                end_str = extract_info_field(info, 'END')
                svlen_str = extract_info_field(info, 'SVLEN')
                
                # Determine end position
                if end_str:
                    end_pos = int(end_str)
                elif svlen_str:
                    svlen = int(svlen_str)
                    if svlen < 0:
                        # For deletions, END = POS + abs(SVLEN)
                        end_pos = pos + abs(svlen)
                    else:
                        # For insertions, END = POS (single position)
                        end_pos = pos
                else:
                    # Default: single position variant
                    end_pos = pos
                
                # Convert to BED format (0-based start)
                bed_start = pos - 1
                bed_end = end_pos
                
                # Write primary BED line
                bed.write(f"{chrom}\t{bed_start}\t{bed_end}\t{var_id}\t{ref}\t{alt}\t{qual}\t{filter_field}\t{info}\t{format_field}\n")
                variant_count += 1
                
                # For BND variants, optionally create duplicate line with parsed ALT information
                if svtype == 'BND' and parse_bnd:
                    bnd_chrom, bnd_pos = parse_bnd_alt(alt)
                    if bnd_chrom and bnd_pos:
                        # Handle edge case where position is 0 (telomeric regions)
                        # In BED format, position 0 becomes 0 (not -1)
                        if bnd_pos == 0:
                            bnd_start = 0
                            bnd_end = 1  # Use 1 for end position to have valid interval
                        else:
                            # Normal case: convert to 0-based
                            bnd_start = bnd_pos - 1
                            bnd_end = bnd_pos
                        
                        bed.write(f"{bnd_chrom}\t{bnd_start}\t{bnd_end}\t{var_id}\t{ref}\t{alt}\t{qual}\t{filter_field}\t{info}\t{format_field}\n")
                        variant_count += 1
                        parsed_bnd_count += 1
                    else:
                        unparsed_bnd_count += 1
                        # Optionally log unparsed BND for debugging
                        # print(f"Warning: Could not parse BND ALT field: {var_id}\t{alt}", file=sys.stderr)
            
            print(f"Successfully converted {variant_count} variant lines from {input_file} to {output_file}")
            if skipped_ins_count > 0:
                print(f"Skipped {skipped_ins_count} INS variants (use --include-INS to include them)")
            if parse_bnd and (parsed_bnd_count > 0 or unparsed_bnd_count > 0):
                print(f"BND parsing: {parsed_bnd_count} parsed, {unparsed_bnd_count} unparsed")
            
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found", file=sys.stderr)
        sys.exit(1)
    except PermissionError:
        print(f"Error: Permission denied accessing '{input_file}' or '{output_file}'", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)

def main():
    args = parse_arguments()
    
    # Determine output file
    if args.output:
        output_file = args.output
    else:
        # Replace .vcf extension with .bed
        base_name = os.path.splitext(args.input_vcf)[0]
        output_file = base_name + '.bed'
    
    # Convert VCF to BED with enhanced options
    vcf_to_bed(args.input_vcf, output_file, 
               parse_bnd=(not args.no_BND_parse), 
               include_ins=args.include_INS)

if __name__ == '__main__':
    main()