#!/usr/bin/env python3
"""
TRACE.py - Transposon Repeat Annotation and Characterization Engine

Orchestrates the entire pipeline for identifying and annotating transposon repeats
from a Sniffles structural variant VCF file.
"""

import argparse
import subprocess
import os
import sys
import logging
from pathlib import Path
from datetime import datetime
import shutil

__author__ = "Michal Izydorczyk, Nicola Wong, John Adedeji, Julian Chiu, and Thomas X. Garcia"

def setup_logging(debug: bool, log_file: str) -> None:
    """Set up logging configuration."""
    level = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )

def check_dependencies() -> None:
    """Check if required external tools are available."""
    required_tools = ['RepeatMasker', 'bedtools']
    missing_tools = []
    
    for tool in required_tools:
        if shutil.which(tool) is None:
            missing_tools.append(tool)
    
    if missing_tools:
        logging.error(f"Missing required tools: {', '.join(missing_tools)}")
        logging.error("Please install these tools and ensure they are in your PATH")
        sys.exit(1)

def validate_file_exists(file_path: str, description: str) -> None:
    """Validate that a file exists."""
    if not Path(file_path).exists():
        logging.error(f"{description} not found: {file_path}")
        sys.exit(1)

def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Run the TRACE pipeline for transposon repeat annotation.',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        'input_vcf',
        help='Input VCF file from Sniffles'
    )
    parser.add_argument(
        '--keep-intermediates',
        action='store_true',
        help='Keep all intermediate files generated during the pipeline'
    )
    parser.add_argument(
        '--remove-log',
        action='store_true',
        help='Remove the pipeline log file after completion (default: keep)'
    )
    parser.add_argument(
        '--intersect',
        required=True,
        help='Path to the file to intersect with (e.g., hg38.fa.bed)'
    )
    parser.add_argument(
        '--lib',
        required=True,
        help='Path to the RepeatMasker library (e.g., Dfam-RepeatMasker.lib)'
    )
    parser.add_argument(
        '--threads',
        type=int,
        default=1,
        help='Number of threads to use for RepeatMasker'
    )
    parser.add_argument(
        '--log-file',
        help='Path to log file (default: TRACE_pipeline_TIMESTAMP.log)'
    )
    parser.add_argument(
        '--debug',
        action='store_true',
        help='Enable debug logging'
    )
    return parser.parse_args()

def run_command(command: list, log_file: str):
    """Run a command and log its output."""
    logging.info(f"Running command: {' '.join(command)}")
    with open(log_file, 'a') as f:
        f.write(f"Running command: {' '.join(command)}\n")
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True
        )
        for line in process.stdout:
            f.write(line)
            sys.stdout.write(line)
        process.wait()
    if process.returncode != 0:
        logging.error(f"Command failed with exit code {process.returncode}")
        logging.error(f"Check log file for details: {log_file}")
        sys.exit(1)

def main():
    """Main function to orchestrate the TRACE pipeline."""
    args = parse_arguments()
    
    # Set up log file with timestamp if not specified
    if args.log_file:
        log_file = args.log_file
    else:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_file = f"TRACE_pipeline_{timestamp}.log"
    
    setup_logging(args.debug, log_file)
    
    # Check dependencies
    logging.info("Checking for required dependencies...")
    check_dependencies()
    
    # Validate input files
    logging.info("Validating input files...")
    validate_file_exists(args.input_vcf, "Input VCF file")
    validate_file_exists(args.intersect, "Intersect file")
    validate_file_exists(args.lib, "RepeatMasker library")

    # Get the directory of the TRACE.py script
    script_dir = Path(__file__).parent.resolve()
    
    # Create output directory for RepeatMasker
    repeatmasker_dir = Path("repeatmasker_output")
    repeatmasker_dir.mkdir(exist_ok=True)

    # --- INS Pipeline ---
    logging.info("--- Starting INS Pipeline ---")
    # 1. vcf2fasta.py
    ins_fasta = Path(args.input_vcf).stem + "_INS.fasta"
    cmd = [
        sys.executable,
        str(script_dir / "vcf2fasta.py"),
        args.input_vcf,
        "-o",
        ins_fasta
    ]
    run_command(cmd, log_file)
    
    # Validate output
    validate_file_exists(ins_fasta, "INS FASTA output")

    # 2. RepeatMasker
    cmd = [
        "RepeatMasker",
        "-no_is",
        "-lib",
        args.lib,
        "-pa",
        str(args.threads),
        "-dir",
        str(repeatmasker_dir),
        ins_fasta
    ]
    run_command(cmd, log_file)

    # 3. out2bed.py
    rm_out = repeatmasker_dir / f"{ins_fasta}.out"
    ins_bed = repeatmasker_dir / f"{ins_fasta}.bed"
    
    # Validate RepeatMasker output
    validate_file_exists(str(rm_out), "RepeatMasker output")
    
    cmd = [
        sys.executable,
        str(script_dir / "out2bed.py"),
        str(rm_out),
        str(ins_bed)
    ]
    run_command(cmd, log_file)
    
    # Validate output
    validate_file_exists(str(ins_bed), "out2bed output")

    # 4. squash_repeat_masker.py
    ins_squashed = str(ins_bed)[:-4] + "_squashed.tsv"
    cmd = [
        sys.executable,
        str(script_dir / "squash_repeat_masker.py"),
        str(ins_bed),
        "-o",
        ins_squashed
    ]
    run_command(cmd, log_file)
    
    # Validate output
    validate_file_exists(ins_squashed, "squash_repeat_masker output")

    # --- DEL, DUP, BND Pipeline ---
    logging.info("--- Starting DEL, DUP, BND Pipeline ---")
    # 5. vcf2bed.py
    sv_bed = Path(args.input_vcf).stem + ".bed"
    cmd = [
        sys.executable,
        str(script_dir / "vcf2bed.py"),
        args.input_vcf,
        "-o",
        sv_bed
    ]
    run_command(cmd, log_file)
    
    # Validate output
    validate_file_exists(sv_bed, "vcf2bed output")

    # 6. bed_processor.py
    processed_bed = Path(sv_bed).stem + "_processed.bed"
    cmd = [
        sys.executable,
        str(script_dir / "bed_processor.py"),
        sv_bed,
        "-o",
        processed_bed
    ]
    run_command(cmd, log_file)
    
    # Validate output
    validate_file_exists(processed_bed, "bed_processor output")

    # 7. bedtools intersect
    intersected_bed = Path(processed_bed).stem + "_intersected.bed"
    cmd = [
        "bedtools",
        "intersect",
        "-wa",
        "-wb",
        "-header",
        "-b",
        args.intersect,
        "-a",
        processed_bed
    ]
    # Redirect output to file
    logging.info(f"Running command: {' '.join(cmd)} > {intersected_bed}")
    with open(intersected_bed, 'w') as f_out, open(log_file, 'a') as f_log:
        f_log.write(f"Running command: {' '.join(cmd)} > {intersected_bed}\n")
        process = subprocess.Popen(cmd, stdout=f_out, stderr=subprocess.PIPE, text=True)
        _, stderr = process.communicate()
        if stderr:
            f_log.write(stderr)
            if process.returncode != 0:
                sys.stderr.write(stderr)
    if process.returncode != 0:
        logging.error(f"bedtools intersect failed with exit code {process.returncode}")
        sys.exit(1)
    
    # Validate output
    validate_file_exists(intersected_bed, "bedtools intersect output")

    # 8. squash_intersect.py
    intersect_squashed = Path(intersected_bed).stem + "_squashed.tsv"
    cmd = [
        sys.executable,
        str(script_dir / "squash_intersect.py"),
        intersected_bed,
        "-o",
        str(intersect_squashed)
    ]
    run_command(cmd, log_file)
    
    # Validate output
    validate_file_exists(str(intersect_squashed), "squash_intersect output")

    # --- Annotation ---
    logging.info("--- Annotating VCF ---")
    # 9. annotate_vcf.py
    annotated_vcf = Path(args.input_vcf).stem + "_annotated.vcf"
    cmd = [
        sys.executable,
        str(script_dir / "annotate_vcf.py"),
        args.input_vcf,
        "--data",
        ins_squashed,
        str(intersect_squashed),
        "-o",
        annotated_vcf
    ]
    run_command(cmd, log_file)
    
    # Validate final output
    validate_file_exists(annotated_vcf, "Final annotated VCF")

    # --- Cleanup ---
    if not args.keep_intermediates:
        logging.info("--- Cleaning up intermediate files ---")
        files_to_remove = [
            ins_fasta,
            str(rm_out),
            str(ins_bed),
            ins_squashed,
            sv_bed,
            processed_bed,
            intersected_bed,
            str(intersect_squashed)
        ]
        
        # Remove log file only if explicitly requested
        if args.remove_log:
            files_to_remove.append(log_file)
        
        for f in files_to_remove:
            try:
                Path(f).unlink()
                logging.info(f"Removed {f}")
            except OSError as e:
                logging.warning(f"Could not remove {f}: {e}")
        
        # Remove the repeatmasker directory
        try:
            shutil.rmtree(repeatmasker_dir)
            logging.info(f"Removed directory {repeatmasker_dir}")
        except OSError as e:
            logging.warning(f"Could not remove directory {repeatmasker_dir}: {e}")

    logging.info(f"--- Pipeline complete! Final output: {annotated_vcf} ---")
    if not args.remove_log:
        logging.info(f"Log file saved as: {log_file}")

if __name__ == "__main__":
    main()