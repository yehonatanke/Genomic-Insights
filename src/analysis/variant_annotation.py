from typing import Dict, Any, List
import subprocess 
import os 

Variant = Dict[str, Any] 
Annotation = Dict[str, Any] 
GeneFeature = Dict[str, Any] 

def load_gene_annotations_from_gff_gtf(annotation_file_path: str) -> List[GeneFeature]:
    """
    Loads gene and transcript annotations from a GFF3 or GTF file.
    This is a very simplified parser. For robust parsing, use libraries like 'gffutils'.
    It extracts gene_id, type, chr, start, end, strand.
    """
    print(f"Loading gene annotations from {annotation_file_path}")
    gene_annotations: List[GeneFeature] = []
    
    if not os.path.exists(annotation_file_path):
        print(f"Error: Annotation file not found: {annotation_file_path}")
        return gene_annotations

    try:
        with open(annotation_file_path, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue

                feature_chr = fields[0]
                feature_type = fields[2]
                feature_start = int(fields[3])
                feature_end = int(fields[4])
                feature_strand = fields[6]
                attributes_str = fields[8]
                
                gene_id = None
                attributes_parts = attributes_str.split(';')
                for part in attributes_parts:
                    part = part.strip()
                    if 'gene_id "' in part:
                        gene_id = part.split('gene_id "')[1].split('"')[0]
                        break
                    elif 'GeneID=' in part: 
                        gene_id = part.split('GeneID=')[1].split(';')[0]
                        break
                
                if gene_id: 
                    gene_annotations.append({
                        'gene_id': gene_id,
                        'type': feature_type,
                        'chr': feature_chr,
                        'start': feature_start,
                        'end': feature_end,
                        'strand': feature_strand
                    })

        print(f"Loaded {len(gene_annotations)} features with gene_id from {annotation_file_path}")
    except Exception as e:
        print(f"Error parsing annotation file {annotation_file_path}: {e}")
    
    return gene_annotations


def predict_variant_effect(variant: Variant, gene_annotations: List[GeneFeature]) -> List[Annotation]:
    """
    Predicts the effect of a single variant based on gene annotations.
    This is a highly simplified version. Real effect prediction is complex and
    considers coding sequences, translation, splice sites, etc.
    Tools like SnpEff or VEP are recommended for accurate annotation.
    """
    annotations: List[Annotation] = []
    chrom = variant.get('CHROM')
    pos = int(variant.get('POS', 0)) 
    ref_allele = variant.get('REF', '')
    alt_allele = variant.get('ALT', '')

    if not all([chrom, pos > 0, ref_allele, alt_allele]):
        return [{'error': 'Incomplete or invalid variant information', 'variant_data': variant}]

    variant_start = pos
    variant_end = pos + len(ref_allele) -1 
    if len(ref_allele) < len(alt_allele): 
        variant_end = pos + len(ref_allele) 

    found_overlap = False
    for feature in gene_annotations:
        if feature.get('chr') == chrom and \
           max(variant_start, feature.get('start', -1)) <= min(variant_end, feature.get('end', -1)):
            
            found_overlap = True
            effect = "unknown_genic"
            impact = "MODIFIER" 

            feature_type = feature.get('type', '').lower()

            if 'exon' in feature_type:
                effect = "exonic"
                # Simplified impact based on type of change
                if len(ref_allele) != len(alt_allele): # Indel
                    if abs(len(ref_allele) - len(alt_allele)) % 3 != 0:
                        effect = "frameshift_variant"
                        impact = "HIGH"
                    else:
                        effect = "inframe_indel"
                        impact = "MODERATE"
                else: # SNV
                    effect = "coding_sequence_variant" 
                    impact = "MODERATE" 
            elif 'intron' in feature_type:
                effect = "intronic"
                impact = "MODIFIER"
            elif 'utr' in feature_type or "UTR" in feature_type:
                effect = "UTR_variant"
                impact = "MODIFIER"
            elif 'gene' in feature_type and effect == "unknown_genic":
                 effect = "genic_region" 
                 impact = "MODIFIER"

            annotations.append({
                'gene_id': feature.get('gene_id', 'N/A'),
                'feature_type': feature.get('type'),
                'feature_coordinates': f"{feature.get('chr')}:{feature.get('start')}-{feature.get('end')}",
                'effect': effect,
                'impact': impact,
                'variant_ref': ref_allele,
                'variant_alt': alt_allele
            })

    if not found_overlap:
        annotations.append({'effect': 'intergenic_variant', 'impact': 'MODIFIER', 'variant_ref': ref_allele, 'variant_alt': alt_allele})

    return annotations


def annotate_variants_from_vcf_snpeff(input_vcf_path: str, output_annotated_vcf_path: str,
                                   snpeff_path: str = "snpEff", genome_version: str = "GRCh38.99",
                                   snpeff_config_path: str = None) -> bool:
    """
    Annotates variants in a VCF file using SnpEff.
    Assumes SnpEff is installed and the specified genome_version database is downloaded.
    """
    print(f"Annotating VCF {input_vcf_path} using SnpEff (genome: {genome_version}) -> {output_annotated_vcf_path}")
    
    output_dir = os.path.dirname(output_annotated_vcf_path)
    if output_dir: os.makedirs(output_dir, exist_ok=True)

    command = [
        "java", "-jar", snpeff_path, 
        "eff",
        "-v", # Verbose
        genome_version,
        input_vcf_path
    ]
    if snpeff_config_path and os.path.exists(snpeff_config_path):
        command.insert(3, "-c") 
        command.insert(4, snpeff_config_path)
    
    try:
        print(f"Running SnpEff: {' '.join(command)}")
        with open(output_annotated_vcf_path, 'w') as f_out:
            result = subprocess.run(command, stdout=f_out, stderr=subprocess.PIPE, check=True, text=True)
        
        if result.stderr: 
             print(f"SnpEff stderr (summary/logs):\n{result.stderr}")
        print(f"SnpEff annotation completed. Output at {output_annotated_vcf_path}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error during SnpEff execution: {e}")
        if hasattr(e, 'stdout') and e.stdout: print(f"SnpEff Stdout: {e.stdout}") 
        if hasattr(e, 'stderr') and e.stderr: print(f"SnpEff Stderr: {e.stderr}")
        return False
    except FileNotFoundError:
        print(f"Error: SnpEff JAR not found at '{snpeff_path}'.")
        return False
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return False
