from typing import Dict, List, Any, Iterator, Tuple
import os 

VCFRecord = Dict[str, Any]
VCFHeader = List[str] 

def get_vcf_header(vcf_file_path: str) -> VCFHeader:
    """Gets only the header lines (starting with '#') from a VCF file."""
    if not os.path.exists(vcf_file_path):
        print(f"Error: VCF file not found at {vcf_file_path}")
        return []
        
    header_lines: VCFHeader = []
    try:
        with open(vcf_file_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    header_lines.append(line.strip())
                else:
                    break # Stop after header section
    except Exception as e:
        print(f"Error reading VCF header from {vcf_file_path}: {e}")
        return []
    return header_lines

def stream_vcf_records(vcf_file_path: str, include_header_obj_in_records: bool = False) -> Iterator[VCFRecord]:
    """
    Streams VCF records (data lines) after skipping the header.
    A proper implementation would use a robust VCF parsing library (e.g., PyVCF).
    This is a simplified parser.

    Args:
        vcf_file_path (str): Path to the VCF file.
        include_header_obj_in_records (bool): If True, attempts to include a parsed header object
                                             (list of column names) in each record under 'HEADER_COLS'.
                                             This is non-standard for VCF records but can be useful.
    Yields:
        Iterator[VCFRecord]: A dictionary representing each VCF record.
    """
    if not os.path.exists(vcf_file_path):
        print(f"Error: VCF file not found at {vcf_file_path}")
        return iter([]) 

    column_names: List[str] = []

    def parse_info_field(info_str: str) -> Dict[str, Any]:
        info_dict = {}
        if info_str == "." or not info_str: return info_dict
        for item in info_str.split(';'):
            parts = item.split('=', 1)
            if len(parts) == 2:
                key, value = parts
                if value.isdigit(): info_dict[key] = int(value)
                elif value.replace('.', '', 1).isdigit(): info_dict[key] = float(value)
                else: info_dict[key] = value
            elif parts[0]: # Flag present
                info_dict[parts[0]] = True
        return info_dict

    try:
        with open(vcf_file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line: continue

                if line.startswith('##'): continue 

                if line.startswith('#CHROM'): 
                    column_names = line.lstrip('#').split('\t')
                    continue
                
                if not column_names: 
                    continue

                fields = line.split('\t')
                record: VCFRecord = dict(zip(column_names, fields))
                
                if 'POS' in record:
                    try: record['POS'] = int(record['POS'])
                    except ValueError: print(f"Warning: Could not convert POS '{record['POS']}' to int. Line: {line}")
                if 'QUAL' in record and record['QUAL'] != '.':
                    try: record['QUAL'] = float(record['QUAL'])
                    except ValueError: print(f"Warning: Could not convert QUAL '{record['QUAL']}' to float. Line: {line}")
                
                # Parse INFO field into a dictionary
                if 'INFO' in record:
                    record['INFO'] = parse_info_field(record['INFO'])
                else: # Should not happen in a valid VCF if INFO column exists
                    record['INFO'] = {}
                
                if 'FORMAT' in record and record['FORMAT'] and len(column_names) > 8 and len(fields) > 8:
                    format_keys = record['FORMAT'].split(':')
                    num_fixed_cols = 8 # CHROM to INFO
                    
                    sample_data_list = []
                    for i in range(num_fixed_cols + 1, len(column_names)): 
                        sample_name = column_names[i]
                        if i < len(fields): 
                            sample_values_str = fields[i]
                            sample_values = sample_values_str.split(':')
                            
                            # Create dict for this sample, handling potential missing values
                            # if len(sample_values) < len(format_keys), some trailing format fields are missing for this sample
                            sample_dict = {}
                            for k_idx, key in enumerate(format_keys):
                                if k_idx < len(sample_values):
                                    sample_dict[key] = sample_values[k_idx]
                                else:
                                    sample_dict[key] = '.' # Placeholder for missing value
                            sample_data_list.append({sample_name: sample_dict})
                        else: # Field for this sample column is missing
                            sample_data_list.append({sample_name: {fk: '.' for fk in format_keys}}) # Empty data for sample
                    record['SAMPLES'] = sample_data_list
                else:
                    record['SAMPLES'] = [] # No samples or no FORMAT field

                if include_header_obj_in_records:
                    record['HEADER_COLS'] = column_names # Non-standard, but for convenience
                
                yield record
    except Exception as e:
        print(f"Error streaming VCF records from {vcf_file_path}: {e}")


def write_vcf_records(records: List[VCFRecord], header_lines: VCFHeader, output_vcf_path: str) -> bool:
    """
    Writes VCF records to a file, including the provided header lines.
    This is a simplified writer. Assumes records are dictionaries.
    """
    print(f"Writing {len(records)} VCF records to {output_vcf_path}")
    
    output_dir = os.path.dirname(output_vcf_path)
    if output_dir: os.makedirs(output_dir, exist_ok=True)

    try:
        with open(output_vcf_path, 'w') as f_out:
            # Write header lines
            for h_line in header_lines:
                f_out.write(h_line + '\n')
            
            chrom_header_line = None
            for h_line in reversed(header_lines):
                if h_line.startswith("#CHROM"):
                    chrom_header_line = h_line
                    break
            
            if not chrom_header_line:
                print("Error: #CHROM line not found in provided VCF header. Cannot determine column order.")
                return False
            
            column_order = chrom_header_line.lstrip("#").split('\t')

            for record in records:
                row_values = []
                for col_name in column_order:
                    if col_name == "INFO":
                        info_parts = []
                        for key, val in sorted(record.get('INFO', {}).items()): 
                            if isinstance(val, bool) and val: 
                                info_parts.append(key)
                            elif val is not None and val != '.':
                                info_parts.append(f"{key}={val}")
                        row_values.append(";".join(info_parts) if info_parts else ".")
                    elif col_name == "FORMAT" and 'SAMPLES' in record and record['SAMPLES']:
                        if record['SAMPLES'][0]: # Check if first sample entry is not empty
                            first_sample_data_dict = list(record['SAMPLES'][0].values())[0]
                            row_values.append(":".join(first_sample_data_dict.keys())) 
                        else:
                            row_values.append(".") # No format if no sample data
                    elif col_name in record: # Standard fixed columns
                        row_values.append(str(record[col_name]))
                    elif 'SAMPLES' in record and record['SAMPLES']: # Sample data columns
                        found_sample_in_record = False
                        for sample_entry_dict in record['SAMPLES']:
                            if col_name in sample_entry_dict: 
                                sample_data_dict = sample_entry_dict[col_name] 
                                format_keys_str = record.get("FORMAT", "")
                                if format_keys_str:
                                    format_keys_list = format_keys_str.split(':')
                                    ordered_sample_values = [str(sample_data_dict.get(fk, ".")) for fk in format_keys_list]
                                    row_values.append(":".join(ordered_sample_values))
                                else: # No FORMAT string, just join what's there (should not happen for valid VCF)
                                    row_values.append(":".join(str(v) for v in sample_data_dict.values()))
                                found_sample_in_record = True
                                break
                        if not found_sample_in_record:
                            row_values.append(".") 
                    else:
                        row_values.append(".") 
                f_out.write("\t".join(row_values) + '\n')
        return True
    except Exception as e:
        print(f"Error writing VCF records to {output_vcf_path}: {e}")
        return False


def load_fasta_sequences(fasta_file_path: str) -> Dict[str, str]:
    """
    Loads sequences from a FASTA file into a dictionary.
    Keys are sequence IDs (first word after '>', before any space).
    Values are the concatenated sequences.
    """
    if not os.path.exists(fasta_file_path):
        print(f"Error: FASTA file not found at {fasta_file_path}")
        return {}

    sequences: Dict[str, str] = {}
    current_sequence_id = None
    current_sequence_parts: List[str] = []

    print(f"Loading sequences from FASTA file: {fasta_file_path}")
    try:
        with open(fasta_file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if current_sequence_id: 
                        sequences[current_sequence_id] = "".join(current_sequence_parts)
                    
                    current_sequence_id = line[1:].split(None, 1)[0] 
                    current_sequence_parts = []
                    if not current_sequence_id: 
                        print(f"Warning: FASTA entry with no ID found. Line: {line}")
                        current_sequence_id = f"UnnamedSequence_{len(sequences)+1}"

                elif current_sequence_id:
                    current_sequence_parts.append(line.upper())
            
            if current_sequence_id:
                sequences[current_sequence_id] = "".join(current_sequence_parts)
                
        print(f"Loaded {len(sequences)} sequences from {fasta_file_path}.")
    except Exception as e:
        print(f"Error loading FASTA file {fasta_file_path}: {e}")
    
    return sequences
