import subprocess 
import os 

# Functions for detecting structural variations

def detect_deletions_from_bam(alignment_path: str, reference_genome_path: str, output_vcf_path: str, delly_executable: str = "delly", exclude_regions_bed: str = None) -> bool:
    """
    Detects deletions from an alignment file (BAM) using Delly.
    Assumes Delly is installed and in PATH or path provided.
    Reference genome must be FASTA format.
    """
    print(f"Detecting deletions from {alignment_path} using Delly -> {output_vcf_path}")
    
    output_dir = os.path.dirname(output_vcf_path)
    if output_dir: 
        os.makedirs(output_dir, exist_ok=True)

    command = [
        delly_executable, "call",
        "-t", "DEL",
        "-g", reference_genome_path,
        "-o", output_vcf_path,
        alignment_path
    ]
    if exclude_regions_bed and os.path.exists(exclude_regions_bed):
        command.extend(["-x", exclude_regions_bed])
        print(f"Excluding regions from: {exclude_regions_bed}")

    try:
        print(f"Running Delly (DEL): {' '.join(command)}")
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        if result.stdout: print(f"Delly (DEL) stdout:\n{result.stdout}")
        if result.stderr: print(f"Delly (DEL) stderr:\n{result.stderr}") 
        print(f"Deletion calling completed. Output at {output_vcf_path}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error during Delly (DEL) execution: {e}")
        print(f"Stdout: {e.stdout}")
        print(f"Stderr: {e.stderr}")
        return False
    except FileNotFoundError:
        print(f"Error: {delly_executable} not found.")
        return False

def detect_duplications_from_bam(alignment_path: str, reference_genome_path: str, output_vcf_path: str, delly_executable: str = "delly", exclude_regions_bed: str = None) -> bool:
    """
    Detects duplications from an alignment file (BAM) using Delly.
    """
    print(f"Detecting duplications from {alignment_path} using Delly -> {output_vcf_path}")
    output_dir = os.path.dirname(output_vcf_path)
    if output_dir: os.makedirs(output_dir, exist_ok=True)

    command = [
        delly_executable, "call",
        "-t", "DUP",
        "-g", reference_genome_path,
        "-o", output_vcf_path,
        alignment_path
    ]
    if exclude_regions_bed and os.path.exists(exclude_regions_bed):
        command.extend(["-x", exclude_regions_bed])
        print(f"Excluding regions from: {exclude_regions_bed}")

    try:
        print(f"Running Delly (DUP): {' '.join(command)}")
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        if result.stdout: print(f"Delly (DUP) stdout:\n{result.stdout}")
        if result.stderr: print(f"Delly (DUP) stderr:\n{result.stderr}")
        print(f"Duplication calling completed. Output at {output_vcf_path}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error during Delly (DUP) execution: {e}")
        print(f"Stdout: {e.stdout}")
        print(f"Stderr: {e.stderr}")
        return False
    except FileNotFoundError:
        print(f"Error: {delly_executable} not found.")
        return False

def detect_inversions_from_bam(alignment_path: str, reference_genome_path: str, output_vcf_path: str, delly_executable: str = "delly", exclude_regions_bed: str = None) -> bool:
    """
    Detects inversions from an alignment file (BAM) using Delly.
    """
    print(f"Detecting inversions from {alignment_path} using Delly -> {output_vcf_path}")
    output_dir = os.path.dirname(output_vcf_path)
    if output_dir: os.makedirs(output_dir, exist_ok=True)

    command = [
        delly_executable, "call",
        "-t", "INV",
        "-g", reference_genome_path,
        "-o", output_vcf_path,
        alignment_path
    ]
    if exclude_regions_bed and os.path.exists(exclude_regions_bed):
        command.extend(["-x", exclude_regions_bed])
        print(f"Excluding regions from: {exclude_regions_bed}")

    try:
        print(f"Running Delly (INV): {' '.join(command)}")
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        if result.stdout: print(f"Delly (INV) stdout:\n{result.stdout}")
        if result.stderr: print(f"Delly (INV) stderr:\n{result.stderr}")
        print(f"Inversion calling completed. Output at {output_vcf_path}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error during Delly (INV) execution: {e}")
        print(f"Stdout: {e.stdout}")
        print(f"Stderr: {e.stderr}")
        return False
    except FileNotFoundError:
        print(f"Error: {delly_executable} not found.")
        return False

def detect_translocations_from_bam(alignment_path: str, reference_genome_path: str, output_vcf_path: str, delly_executable: str = "delly", exclude_regions_bed: str = None) -> bool:
    """
    Detects translocations (inter-chromosomal) from an alignment file (BAM) using Delly.
    """
    print(f"Detecting translocations from {alignment_path} using Delly -> {output_vcf_path}")
    output_dir = os.path.dirname(output_vcf_path)
    if output_dir: os.makedirs(output_dir, exist_ok=True)

    command = [
        delly_executable, "call",
        "-t", "TRA", 
        "-g", reference_genome_path,
        "-o", output_vcf_path,
        alignment_path
    ]
    if exclude_regions_bed and os.path.exists(exclude_regions_bed):
        command.extend(["-x", exclude_regions_bed])
        print(f"Excluding regions from: {exclude_regions_bed}")

    try:
        print(f"Running Delly (TRA/BND): {' '.join(command)}")
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        if result.stdout: print(f"Delly (TRA/BND) stdout:\n{result.stdout}")
        if result.stderr: print(f"Delly (TRA/BND) stderr:\n{result.stderr}")
        print(f"Translocation calling completed. Output at {output_vcf_path}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error during Delly (TRA/BND) execution: {e}")
        print(f"Stdout: {e.stdout}")
        print(f"Stderr: {e.stderr}")
        return False
    except FileNotFoundError:
        print(f"Error: {delly_executable} not found. Please ensure it is installed and in your PATH or provide the correct path.")
        return False

def merge_sv_calls(input_vcf_list: list[str], output_merged_vcf_path: str, survivor_executable: str = "SURVIVOR", distance_bp: int = 1000, min_callers: int = 1, require_type_match: bool = True, require_strand_match: bool = True, min_sv_size: int = 30) -> bool:
    """
    Merges SV calls from multiple VCF files using SURVIVOR.
    Assumes SURVIVOR is installed and in PATH or path provided.
    """
    if not input_vcf_list:
        print("Error: No input VCF files provided for merging.")
        return False
    
    vcf_files_str = " ".join(input_vcf_list) # SURVIVOR expects VCF files as separate arguments
    print(f"Merging {len(input_vcf_list)} SV VCF files using SURVIVOR -> {output_merged_vcf_path}")
    output_dir = os.path.dirname(output_merged_vcf_path)
    if output_dir: os.makedirs(output_dir, exist_ok=True)
    command = [
        survivor_executable, "merge",
    ]
    
    temp_file_list_path = os.path.join(output_dir if output_dir else ".", "survivor_input_vcfs.txt")
    with open(temp_file_list_path, 'w') as f_list:
        for vcf_file in input_vcf_list:
            if not os.path.exists(vcf_file):
                print(f"Warning: Input VCF file for merging not found: {vcf_file}")
            f_list.write(vcf_file + '\n')
    
    command.append(temp_file_list_path) # Path to the file listing VCFs
    command.extend([
        str(distance_bp),
        str(min_callers),
        "1" if require_type_match else "0",
        "1" if require_strand_match else "0",
        "0",  # estimate_distance (0 = no, 1 = yes)
        str(min_sv_size),
        output_merged_vcf_path
    ])

    try:
        print(f"Running SURVIVOR merge: {' '.join(command)}")
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        if result.stdout: print(f"SURVIVOR merge stdout:\n{result.stdout}")
        if result.stderr: print(f"SURVIVOR merge stderr:\n{result.stderr}")
        print(f"SV merging completed. Output at {output_merged_vcf_path}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error during SURVIVOR merge execution: {e}")
        print(f"Stdout: {e.stdout}")
        print(f"Stderr: {e.stderr}")
        return False
    except FileNotFoundError:
        print(f"Error: {survivor_executable} not found.")
        return False
    finally:
        # Clean up the temporary file list
        if os.path.exists(temp_file_list_path):
            os.remove(temp_file_list_path)
