#!/usr/bin/env python3
"""
Convert RepeatMasker output to BED format

Authors: Michal Izydorczyk, Nicola Wong, John Adedeji, Julian Chiu, and Thomas X. Garcia
"""

import sys

def parse_repeatmasker_to_bed(rm_file, bed_file):
    with open(rm_file) as infile, open(bed_file, 'w') as out:
        out.write("#SnifflesID\tref_start\tref_end\trepeat_name\trepeat_family\tstrand\trepeat_start\trepeat_end\tID_number\n")
        for line in infile:
            if line.strip() == "" or line.startswith("   SW"):
                continue  # skip headers/empty lines

            parts = line.strip().split()
            if len(parts) < 15:
                continue  # skip malformed lines

            chrom = parts[4]
            start = int(parts[5]) - 1  # BED is zero-based
            end = int(parts[6])
            strand = '+' if parts[8] == '+' else '-'
            repeat_name = parts[9]
            repeat_family = parts[10]
            id_ = parts[14]

            # Repeat positions in repeat consensus sequence
            # For + strand: col 11 = begin, col 12 = end, col 13 = (left)
            # For C strand: col 11 = (left), col 12 = end, col 13 = begin
            if strand == '+':
                # + strand: columns 11 and 12 are begin and end
                rep_start = parts[11].strip('()')
                rep_stop = parts[12].strip('()')
            else:
                # C strand: column 12 is end, column 13 is begin (reverse orientation)
                # Column 11 is (left) which we don't use
                rep_start = parts[12].strip('()')
                rep_stop = parts[13].strip('()')

            bed_line = f"{chrom}\t{start}\t{end}\t{repeat_name}\t{repeat_family}\t{strand}\t{rep_start}\t{rep_stop}\t{id_}\n"
            out.write(bed_line)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} input.repeatmasker.out output.bed")
        sys.exit(1)

    parse_repeatmasker_to_bed(sys.argv[1], sys.argv[2])