import subprocess 
import os 

# Functions for Single Nucleotide Polymorphism (SNP) and Insertion/Deletion (Indel) calling

def align_reads_to_reference(reads_file_path: str, reference_genome_path: str, output_alignment_path: str, threads: int = 4) -> bool:
    """
    Aligns sequencing reads to a reference genome using BWA-MEM.
    Assumes BWA is installed and in PATH.
    Assumes reference genome is indexed for BWA.
    """
    print(f"Aligning {reads_file_path} to {reference_genome_path} -> {output_alignment_path} using BWA-MEM")
    bwa_executable = "bwa"
    samtools_executable = "samtools" 

    # Ensure output directory exists
    output_dir = os.path.dirname(output_alignment_path)
    if output_dir: 
        os.makedirs(output_dir, exist_ok=True)

    # BWA MEM command to output SAM
    bwa_command = [
        bwa_executable, "mem",
        "-t", str(threads),
        reference_genome_path,
        reads_file_path
    ]
    # Samtools view command to convert SAM to BAM
    # Outputting directly to BAM is more efficient
    samtools_command = [
        samtools_executable, "view",
        "-@ ", str(threads - 1 if threads > 1 else 1), # samtools threads
        "-Sb", # Output BAM
        "-o", output_alignment_path,
        "-" # Read from stdin
    ]

    try:
        print(f"Running BWA: {' '.join(bwa_command)}")
        bwa_process = subprocess.Popen(bwa_command, stdout=subprocess.PIPE, text=False) # text=False for binary pipe to samtools
        
        print(f"Running Samtools view: {' '.join(samtools_command)}")
        samtools_process = subprocess.Popen(samtools_command, stdin=bwa_process.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=False)
        
        # Allow bwa_process to receive a SIGPIPE if samtools_process exits.
        if bwa_process.stdout:
            bwa_process.stdout.close()
        
        stdout, stderr = samtools_process.communicate()
        
        if bwa_process.wait() != 0:
            # bwa might have written to its stderr, which is not captured here directly
            # but samtools_process.stderr might contain some info if bwa failed early
            print(f"Error during BWA execution. BWA stderr might provide clues (not captured directly). Samtools stderr: {stderr.decode() if stderr else 'N/A'}")
            return False

        if samtools_process.returncode != 0:
            print(f"Error during Samtools view execution: {samtools_process.returncode}")
            if stdout: print(f"Samtools stdout: {stdout.decode()}")
            if stderr: print(f"Samtools stderr: {stderr.decode()}")
            return False

        print(f"Alignment and BAM conversion completed. Output at {output_alignment_path}")
        return True

    except subprocess.CalledProcessError as e:
        print(f"Error during alignment: {e}")
        if hasattr(e, 'stdout') and e.stdout: print(f"Stdout: {e.stdout.decode()}")
        if hasattr(e, 'stderr') and e.stderr: print(f"Stderr: {e.stderr.decode()}")
        return False
    except FileNotFoundError:
        print(f"Error: {bwa_executable} or {samtools_executable} not found.")
        return False
    except Exception as e:
        print(f"An unexpected error occurred during alignment: {e}")
        return False


def sort_and_index_alignment(input_bam_path: str, sorted_output_bam_path: str, threads: int = 4) -> bool:
    """
    Sorts and indexes a BAM file using Samtools.
    Assumes Samtools is installed and in PATH.
    """
    print(f"Sorting and indexing {input_bam_path} -> {sorted_output_bam_path}")
    samtools_executable = "samtools"

    # Ensure output directory exists
    output_dir = os.path.dirname(sorted_output_bam_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    sort_command = [
        samtools_executable, "sort",
        "-@ ", str(threads),
        "-o", sorted_output_bam_path,
        input_bam_path
    ]
    index_command = [
        samtools_executable, "index",
        "-@ ", str(threads),
        sorted_output_bam_path
    ]

    try:
        print(f"Running Samtools sort: {' '.join(sort_command)}")
        result_sort = subprocess.run(sort_command, check=True, capture_output=True, text=True)
        if result_sort.stdout: print(f"Samtools sort stdout:\n{result_sort.stdout}")
        if result_sort.stderr: print(f"Samtools sort stderr:\n{result_sort.stderr}")
        
        print(f"Running Samtools index: {' '.join(index_command)}")
        result_index = subprocess.run(index_command, check=True, capture_output=True, text=True)
        if result_index.stdout: print(f"Samtools index stdout:\n{result_index.stdout}")
        if result_index.stderr: print(f"Samtools index stderr:\n{result_index.stderr}")

        print(f"Sorting and indexing completed. Output at {sorted_output_bam_path} and its index file.")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error during Samtools sort/index: {e}")
        print(f"Stdout: {e.stdout}")
        print(f"Stderr: {e.stderr}")
        return False
    except FileNotFoundError:
        print(f"Error: {samtools_executable} not found. Please ensure it is installed and in your PATH.")
        return False

def call_snps_indels(sorted_bam_path: str, reference_genome_path: str, output_vcf_path: str, gatk_path: str = "gatk") -> bool:
    """
    Calls SNPs and Indels from a sorted BAM file using GATK HaplotypeCaller.
    Assumes GATK is installed and accessible via gatk_path.
    Reference genome must be indexed (e.g., .fai file and .dict file).
    """
    print(f"Calling variants from {sorted_bam_path} against {reference_genome_path} -> {output_vcf_path} using GATK HaplotypeCaller")
    
    # Ensure output directory exists
    output_dir = os.path.dirname(output_vcf_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    command = [
        gatk_path, "HaplotypeCaller",
        "-R", reference_genome_path,
        "-I", sorted_bam_path,
        "-O", output_vcf_path
    ]
    try:
        print(f"Running GATK HaplotypeCaller: {' '.join(command)}")
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        if result.stdout: print(f"GATK HaplotypeCaller stdout:\n{result.stdout}")
        if result.stderr: print(f"GATK HaplotypeCaller stderr:\n{result.stderr}") # Progress and potential info
        print(f"Variant calling completed. Output at {output_vcf_path}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error during GATK HaplotypeCaller execution: {e}")
        print(f"Stdout: {e.stdout}")
        print(f"Stderr: {e.stderr}") 
        return False
    except FileNotFoundError:
        print(f"Error: {gatk_path} not found.")
        return False

def filter_variants(input_vcf_path: str, output_filtered_vcf_path: str, reference_genome_path: str, gatk_path: str = "gatk",
                    filter_expression: str = "QD < 2.0 || FS > 60.0 || MQ < 40.0",
                    filter_name: str = "BasicFilters") -> bool:
    """
    Filters variants using GATK VariantFiltration.
    """
    print(f"Filtering {input_vcf_path} -> {output_filtered_vcf_path} using GATK VariantFiltration")
    
    output_dir = os.path.dirname(output_filtered_vcf_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    command = [
        gatk_path, "VariantFiltration",
        "-R", reference_genome_path,
        "-V", input_vcf_path,
        "-O", output_filtered_vcf_path,
        "--filter-expression", filter_expression,
        "--filter-name", filter_name
    ]
    try:
        print(f"Running GATK VariantFiltration: {' '.join(command)}")
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        if result.stdout: print(f"GATK VariantFiltration stdout:\n{result.stdout}")
        if result.stderr: print(f"GATK VariantFiltration stderr:\n{result.stderr}")
        print(f"Variant filtering completed. Output at {output_filtered_vcf_path}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error during GATK VariantFiltration execution: {e}")
        print(f"Stdout: {e.stdout}")
        print(f"Stderr: {e.stderr}")
        return False
    except FileNotFoundError:
        print(f"Error: {gatk_path} not found.")
        return False

