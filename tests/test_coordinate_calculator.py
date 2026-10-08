import os
import re
import sys

# Ensure parent directory is in path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import pandas as pd
import pytest

from data_utils import DataLoader, GeneInfo
from coordinate_calculator import CoordinateCalculator

DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'data'))


@pytest.fixture(scope="module")
def loader():
    instance = DataLoader(DATA_DIR)
    instance.load_gff3()
    return instance


@pytest.fixture(scope="module")
def calculator(loader):
    return CoordinateCalculator(loader)


def _gene(start, end, strand):
    return GeneInfo(gene_id="g", gene_name="g", locus_tag="g", start=start, end=end,
                    strand=strand, product="", chromosome="NC_000962.3")


def test_calculate_genomic_position():
    # Initialize DataLoader with the data directory
    loader = DataLoader("data")

    # Initialize CoordinateCalculator with the data loader
    calculator = CoordinateCalculator(loader)

    # Calculate the full coordinates for variant dnaA_c.7G>A
    result = calculator.calculate_full_coordinates("dnaA_c.7G>A", "dnaA")

    # Assert that genomic_position equals 7
    assert result.get("genomic_position") == 7


# ----------------------------------------------------------------------
# Genomic -> c. on both strands, either side of c.1 (HGVS has no c.0)
# ----------------------------------------------------------------------
@pytest.mark.parametrize("genomic, expected", [
    (1000, 1),     # first coding base
    (1002, 3),     # last base of codon 1
    (1003, 4),     # first base of codon 2
    (1999, 1000),  # last coding base
    (999, -1),     # base just before c.1
    (985, -15),
])
def test_relative_position_on_the_plus_strand(calculator, genomic, expected):
    assert calculator.calculate_relative_position(genomic, _gene(1000, 1999, "+")) == expected


@pytest.mark.parametrize("genomic, expected", [
    (1999, 1),     # first coding base is the gene end
    (1997, 3),
    (1000, 1000),  # last coding base is the gene start
    (2000, -1),    # base just before c.1
    (2014, -15),
])
def test_relative_position_on_the_minus_strand(calculator, genomic, expected):
    assert calculator.calculate_relative_position(genomic, _gene(1000, 1999, "-")) == expected


@pytest.mark.parametrize("strand", ["+", "-"])
def test_relative_position_never_returns_zero(calculator, strand):
    gene = _gene(1000, 1999, strand)
    positions = [calculator.calculate_relative_position(g, gene) for g in range(950, 2050)]
    assert 0 not in positions


@pytest.mark.parametrize("strand", ["+", "-"])
def test_genomic_and_relative_positions_round_trip(calculator, strand):
    gene = _gene(1000, 1999, strand)
    for genomic in range(950, 2050):
        relative = calculator.calculate_relative_position(genomic, gene)
        assert calculator.calculate_genomic_from_c_dot(f"c.{relative}A>G", gene) == genomic


def test_c_dot_zero_is_not_a_position(calculator):
    assert calculator.calculate_genomic_from_c_dot("c.0A>G", _gene(1000, 1999, "+")) is None


def test_c_dot_from_genomic_matches_relative_position(calculator):
    gene = _gene(1000, 1999, "-")
    assert calculator.calculate_c_dot_from_genomic(2000, gene) == -1
    assert calculator.calculate_c_dot_from_genomic(1999, gene) == 1


@pytest.mark.parametrize("nucleotide, codon", [(1, 1), (3, 1), (4, 2), (361, 121), (945, 315)])
def test_nucleotide_to_codon(calculator, nucleotide, codon):
    assert calculator.nucleotide_to_amino_acid_position(nucleotide) == codon


# ----------------------------------------------------------------------
# Real genes, against the positions in the WHO genomic coordinates file
# ----------------------------------------------------------------------
@pytest.mark.parametrize("gene, genomic, expected", [
    ("Rv0678", 779350, 361),    # + strand, p.Gly121Arg
    ("katG", 2155167, 945),     # - strand, p.Ser315Thr
    ("katG", 2156112, -1),      # - strand, katG_c.-1T>C
    ("fabG1", 1673425, -15),    # + strand, fabG1_c.-15C>T
    ("dnaA", 4411221, -312),    # dnaA starts at 1: its upstream wraps the circular chromosome
])
def test_relative_position_of_catalogue_variants(loader, calculator, gene, genomic, expected):
    gene_info = loader.get_gene_info(gene)
    assert calculator.calculate_relative_position(genomic, gene_info) == expected
    assert calculator.calculate_genomic_from_c_dot(f"c.{expected}", gene_info) == genomic


def test_upstream_variant_full_coordinates(calculator):
    result = calculator.calculate_full_coordinates("fabG1_c.-15C>T", "fabG1")
    assert result["genomic_position"] == 1673425
    assert result["relative_nucleotide_position"] == -15
    # Upstream of the start codon there is no amino acid
    assert result["amino_acid_position"] is None


def test_chromosome_length_comes_from_the_fasta_index(calculator):
    assert calculator.chromosome_length() == 4411532


def test_every_single_base_c_dot_variant_matches_the_who_coordinates(loader, calculator):
    coordinates = pd.read_csv(os.path.join(DATA_DIR, "genomic_coordinates.txt"), sep="\t")
    pattern = re.compile(r"^(.+?)_c\.(-?\d+)[ACGT]>[ACGT]$")
    mismatches = []
    checked = 0
    columns = ["variant", "position", "reference_nucleotide", "alternative_nucleotide"]
    for variant, position, reference, alternative in coordinates[columns].itertuples(index=False):
        match = pattern.match(str(variant))
        # The file also lists multi-base encodings of a variant; only a
        # single-base change sits exactly at the variant's c. position.
        if not match or len(str(reference)) != 1 or len(str(alternative)) != 1:
            continue
        gene_info = loader.get_gene_info(match.group(1))
        expected = int(match.group(2))
        checked += 1
        if (calculator.calculate_relative_position(int(position), gene_info) != expected
                or calculator.calculate_genomic_from_c_dot(f"c.{expected}", gene_info) != int(position)):
            mismatches.append((variant, position))
    assert checked > 10000
    assert mismatches == []
