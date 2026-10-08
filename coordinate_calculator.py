"""
Coordinate calculator for converting between gene-relative and genomic positions.
This utility bridges the gap between the catalogue master file and genomic coordinates.
"""
import os
import re
from typing import Dict, List, Optional

import pandas as pd

from data_utils import DataLoader, GeneInfo


class CoordinateCalculator:
    """
    Calculate and convert between gene-relative and genomic coordinates.
    
    The catalogue master file uses gene-relative positions (e.g., c.102G>A),
    while the genomic coordinates file uses absolute genomic positions.
    This class provides utilities to convert between them.
    """
    
    def __init__(self, data_loader: DataLoader):
        self.data_loader = data_loader
        self._chromosome_length: Optional[int] = None
        self._chromosome_length_read = False
    
    def chromosome_length(self) -> Optional[int]:
        """Length of the chromosome, read once from the FASTA index."""
        if not self._chromosome_length_read:
            self._chromosome_length_read = True
            fai_path = getattr(self.data_loader, "fai_path", None)
            if fai_path and os.path.exists(fai_path):
                with open(fai_path) as handle:
                    fields = handle.readline().split("\t")
                if len(fields) > 1 and fields[1].isdigit():
                    self._chromosome_length = int(fields[1])
        return self._chromosome_length
    
    def _shortest_offset(self, offset: int) -> int:
        """
        The chromosome is circular: dnaA starts at position 1, so its upstream
        region lies at the end of the sequence. Measure the offset the short
        way round the circle.
        """
        length = self.chromosome_length()
        if length:
            if offset > length // 2:
                offset -= length
            elif offset < -(length // 2):
                offset += length
        return offset
    
    def calculate_relative_position(self, genomic_position: int, gene_info: GeneInfo) -> int:
        """
        Calculate the HGVS c. position of a genomic position.
        For + strand: relative_pos = genomic_pos - gene_start + 1
        For - strand: relative_pos = gene_end - genomic_position + 1
        Upstream of the gene the + 1 is dropped, because HGVS has no c.0:
        the base just before c.1 is c.-1.
        """
        if gene_info.strand == '+':
            offset = genomic_position - gene_info.start
        else:
            offset = gene_info.end - genomic_position
        offset = self._shortest_offset(offset)
        return offset + 1 if offset >= 0 else offset
    
    def parse_c_dot_notation(self, c_dot: str) -> Optional[int]:
        """
        Parse HGVS c. notation to extract position.
        Examples:
            c.102G>A -> 102
            c.-52A>C -> -52 (upstream)
            c.1044G>A -> 1044
            c.315C>T -> 315
        """
        if not c_dot or pd.isna(c_dot):
            return None
        
        # Match patterns like c.102, c.-52, c.1044
        match = re.search(r'c\.(-?\d+)', str(c_dot))
        if match:
            return int(match.group(1))
        
        return None
    
    def parse_p_notation(self, p_dot: str) -> Optional[int]:
        """
        Parse HGVS p. notation to extract amino acid position.
        Examples:
            p.Asp3Ala -> 3
            p.Glu107_Asp109del -> 107 (start of deletion)
            p.Thr112del -> 112
        """
        if not p_dot or pd.isna(p_dot):
            return None
        
        # Match patterns like p.Asp3, p.Glu107
        match = re.search(r'p\.[A-Z][a-z]+(\d+)', str(p_dot))
        if match:
            return int(match.group(1))
        
        return None
    
    def nucleotide_to_amino_acid_position(self, nucleotide_position: int) -> int:
        """Convert nucleotide position to amino acid position (1-based)."""
        return (nucleotide_position - 1) // 3 + 1
    
    def amino_acid_to_nucleotide_position(self, amino_acid_position: int) -> int:
        """Convert amino acid position to first nucleotide position (1-based)."""
        return (amino_acid_position - 1) * 3 + 1
    
    def calculate_genomic_from_c_dot(
        self, 
        c_dot: str, 
        gene_info: GeneInfo
    ) -> Optional[int]:
        """
        Calculate genomic position from c. notation.
        
        Args:
            c_dot: HGVS c. notation (e.g., "c.102G>A")
            gene_info: Gene information with start, end, strand
        
        Returns:
            Genomic position (1-based)
        """
        rel_pos = self.parse_c_dot_notation(c_dot)
        if rel_pos is None or rel_pos == 0:
            return None
        
        # Upstream variants (negative positions) have no c.0 to skip over,
        # so their offset from c.1 is the position itself.
        offset = rel_pos - 1 if rel_pos > 0 else rel_pos
        if gene_info.strand == '+':
            position = gene_info.start + offset
        else:
            position = gene_info.end - offset
        length = self.chromosome_length()
        return (position - 1) % length + 1 if length else position
    
    def calculate_c_dot_from_genomic(
        self,
        genomic_position: int,
        gene_info: GeneInfo
    ) -> Optional[int]:
        """
        Calculate c. position from genomic position.
        
        Args:
            genomic_position: Genomic position (1-based)
            gene_info: Gene information with start, end, strand
        
        Returns:
            Relative nucleotide position (c. notation number)
        """
        return self.calculate_relative_position(genomic_position, gene_info)
    
    def calculate_full_coordinates(
        self,
        variant: str,
        gene_name: str
    ) -> Dict:
        """
        Calculate all coordinate information for a variant.
        
        Args:
            variant: Variant name (e.g., "dnaA_c.102G>A" or "dnaA_p.Asp3Ala")
            gene_name: Gene name
        
        Returns:
            Dictionary with all coordinate information
        """
        gene_info = self.data_loader.get_gene_info(gene_name)
        if not gene_info:
            return {"error": f"Gene {gene_name} not found"}
        
        # Parse variant to extract mutation type and position
        parts = variant.split('_', 1)
        if len(parts) != 2:
            return {"error": f"Invalid variant format: {variant}"}
        
        mutation = parts[1]
        
        result = {
            "variant": variant,
            "gene": gene_name,
            "gene_info": {
                "name": gene_info.gene_name,
                "locus_tag": gene_info.locus_tag,
                "start": gene_info.start,
                "end": gene_info.end,
                "strand": gene_info.strand,
                "product": gene_info.product
            },
            "mutation": mutation
        }
        
        # Handle c. notation
        if mutation.startswith('c.'):
            rel_pos = self.parse_c_dot_notation(mutation)
            result["relative_nucleotide_position"] = rel_pos
            
            if rel_pos is not None:
                genomic_pos = self.calculate_genomic_from_c_dot(mutation, gene_info)
                result["genomic_position"] = genomic_pos
                # Upstream positions lie outside the coding sequence: no codon.
                result["amino_acid_position"] = (
                    self.nucleotide_to_amino_acid_position(rel_pos) if rel_pos > 0 else None
                )
        
        # Handle p. notation
        elif mutation.startswith('p.'):
            aa_pos = self.parse_p_notation(mutation)
            result["amino_acid_position"] = aa_pos
            
            if aa_pos is not None:
                nuc_pos = self.amino_acid_to_nucleotide_position(aa_pos)
                result["relative_nucleotide_position"] = nuc_pos
                
                genomic_pos = self.calculate_genomic_from_c_dot(f"c.{nuc_pos}", gene_info)
                result["genomic_position"] = genomic_pos
        
        return result
    
    def get_mutation_with_coordinates(
        self,
        variant: str,
        gene_name: str,
        drug: Optional[str] = None
    ) -> Dict:
        """
        Get complete mutation information with all coordinates and drug resistance data.
        
        Args:
            variant: Variant name
            gene_name: Gene name
            drug: Optional drug name to filter resistance data
        
        Returns:
            Complete mutation information dictionary
        """
        # Get coordinate information
        coord_info = self.calculate_full_coordinates(variant, gene_name)
        
        # Get drug resistance information from catalogue
        if drug:
            resistance_info = self.data_loader.get_drug_resistance_info(gene_name, variant)
        else:
            resistance_info = self.data_loader.get_drug_resistance_info(gene_name)
        
        if len(resistance_info) > 0:
            row = resistance_info.iloc[0]
            coord_info["drug"] = row.get('drug')
            coord_info["tier"] = row.get('tier')
            coord_info["effect"] = row.get('effect')
            coord_info["confidence"] = row.get('FINAL CONFIDENCE GRADING')
            coord_info["comment"] = row.get('Comment')
            
            # Add sensitivity/specificity if available
            if 'Sens_DATASET ALL' in row:
                coord_info["sensitivity"] = row.get('Sens_DATASET ALL')
            if 'Spec_DATASET ALL' in row:
                coord_info["specificity"] = row.get('Spec_DATASET ALL')
            if 'PPV_DATASET ALL' in row:
                coord_info["ppv"] = row.get('PPV_DATASET ALL')
        
        # Get all genomic coordinate entries for this variant
        genomic_coords = self.data_loader.search_mutations_by_gene(gene_name)
        variant_mask = genomic_coords['variant'] == variant
        variant_coords = genomic_coords[variant_mask]
        
        coord_info["all_nucleotide_changes"] = variant_coords.to_dict('records') if len(variant_coords) > 0 else []
        
        return coord_info
    
    def batch_calculate_coordinates(
        self,
        variants: List[str],
        gene_name: str
    ) -> List[Dict]:
        """
        Calculate coordinates for multiple variants.
        
        Args:
            variants: List of variant names
            gene_name: Gene name
        
        Returns:
            List of coordinate information dictionaries
        """
        results = []
        for variant in variants:
            result = self.get_mutation_with_coordinates(variant, gene_name)
            results.append(result)
        return results
