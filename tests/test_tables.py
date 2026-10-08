import json
import os
import sys

# Ensure parent directory is in path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import pandas as pd
import pytest

import tables
from coordinate_calculator import CoordinateCalculator
from data_utils import DataLoader, RELAXED_THRESHOLDS_COLUMN

DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'data'))


@pytest.fixture(scope="module")
def loader():
    instance = DataLoader(DATA_DIR)
    instance.load_gff3()
    return instance


@pytest.fixture(scope="module")
def calculator(loader):
    return CoordinateCalculator(loader)


# ----------------------------------------------------------------------
# TASK-12: full WHO catalogue schema
# ----------------------------------------------------------------------
def test_all_who_metadata_columns_are_offered():
    names = [column["name"] for column in tables.build_resistance_columns()]
    for expected in ["Mutation", "Tier", "Confidence", "Effect", "Comment",
                     "CHANGES vs ver1", "Relaxed thresholds", "Silent mutation"]:
        assert expected in names


def test_required_columns_survive_hiding_everything():
    """TASK-12.4: optional columns hide, mutation identity and grading stay."""
    names = [column["name"] for column in tables.build_resistance_columns([])]
    assert names == ["Mutation", "Tier", "Confidence"]


def test_resistance_data_is_flattened_onto_stable_column_ids(loader):
    rows = loader.get_drug_resistance_info("katG")
    prepared = tables.prepare_resistance_data(rows)

    for key in [column["key"] for column in tables.RESISTANCE_COLUMNS]:
        assert key in prepared.columns
    assert "variant" in prepared.columns

    # The long WHO header is reachable through a short id.
    source = rows.set_index("mutation")[RELAXED_THRESHOLDS_COLUMN].dropna()
    if len(source) > 0:
        mutation = source.index[0]
        value = prepared.loc[prepared["mutation"] == mutation, "relaxed_thresholds"].iloc[0]
        assert value == source.iloc[0]


def test_silent_mutation_flag_is_rendered_as_yes_or_blank(loader):
    prepared = tables.prepare_resistance_data(loader.get_drug_resistance_info("katG"))
    assert set(prepared["silent_mutation"].unique()) <= {"Yes", ""}
    assert (prepared["silent_mutation"] == "Yes").any()


def test_tier_one_rows_sort_first(loader):
    prepared = tables.prepare_resistance_data(loader.get_drug_resistance_info("Rv0678"))
    tiers = pd.to_numeric(prepared["tier"], errors="coerce").dropna()
    assert list(tiers) == sorted(tiers)


def test_comment_and_version_changes_are_flagged_in_red():
    """TASK-12.3: remarks and version changes get their own styling."""
    styles = json.dumps(tables.resistance_style_conditional())
    assert "comment" in styles and "#B42318" in styles
    assert "changes_vs_ver1" in styles


def test_selected_mutation_is_highlighted_in_either_notation():
    styles = json.dumps(tables.resistance_style_conditional("p.Ser315Thr"))
    assert '{mutation} eq \\"p.Ser315Thr\\"' in styles
    assert '{mutation} eq \\"Ser315Thr\\"' in styles


@pytest.mark.parametrize("mutation,page", [
    ("p.Ala0Gly", 0),
    ("p.Ala8Gly", 1),
    ("p.Ala16Gly", 2),
    ("p.Missing", 0),
])
def test_page_of_mutation(mutation, page):
    """TASK-04.4: the table opens on the page holding the match."""
    prepared = pd.DataFrame({"mutation": [f"p.Ala{i}Gly" for i in range(24)]})
    assert tables.page_of_mutation(prepared, mutation, page_size=8) == page


# ----------------------------------------------------------------------
# TASK-07: long alleles and copyable sequence cells
# ----------------------------------------------------------------------
def test_long_alleles_are_truncated_with_their_size(loader, calculator):
    gene = loader.get_gene_info("katG")
    variants = loader.search_mutations_by_gene("katG")
    display = tables.prepare_coordinates_data(variants, gene, calculator)

    longest = display.loc[display["reference_nucleotide"].str.len().idxmax()]
    assert len(longest["reference_nucleotide"]) > 1000
    # The cell stays short, and still says how long the real allele is.
    assert len(longest["ref_display"]) < 60
    assert "bp)" in longest["ref_display"]


def test_short_alleles_are_shown_verbatim(loader, calculator):
    gene = loader.get_gene_info("katG")
    display = tables.prepare_coordinates_data(
        loader.search_mutations_by_gene("katG"), gene, calculator
    )
    snvs = display[display["reference_nucleotide"].str.len() == 1]
    assert (snvs["ref_display"] == snvs["reference_nucleotide"]).all()


@pytest.mark.parametrize("reference,alternative,expected", [
    ("A", "G", "SNV"),
    ("AT", "A", "Deletion (1 bp)"),
    ("A", "ATG", "Insertion (2 bp)"),
    ("AT", "GC", "MNV (2 bp)"),
])
def test_change_description(reference, alternative, expected):
    assert tables._describe_change(reference, alternative) == expected


def test_coordinates_table_enables_text_selection(loader, calculator):
    """TASK-07.3: sequence text must remain selectable so it can be copied."""
    gene = loader.get_gene_info("mmpS5")
    container = tables.create_genomic_coords_table(
        loader.search_mutations_by_gene("mmpS5"), gene, calculator
    )
    table = container.children
    css = json.dumps(table.css)
    assert "user-select: text" in css
    assert table.style_table["overflowX"] == "auto"
    assert "table-responsive-container" in container.className


def test_truncated_cells_get_a_full_sequence_tooltip(loader, calculator):
    gene = loader.get_gene_info("katG")
    container = tables.create_genomic_coords_table(
        loader.search_mutations_by_gene("katG"), gene, calculator
    )
    tooltips = container.children.tooltip_data
    assert any("ref_display" in entry for entry in tooltips)


def test_relative_positions_are_computed_for_the_reverse_strand(loader, calculator):
    gene = loader.get_gene_info("katG")
    display = tables.prepare_coordinates_data(
        loader.search_mutations_by_gene("katG"), gene, calculator
    )
    row = display[display["position_value"] == 2155167].iloc[0]
    # katG is on the minus strand: c. = gene_end - genomic + 1
    assert row["Gene Relative"] == f"c.{gene.end - 2155167 + 1}"


def test_upstream_positions_skip_c_zero(loader, calculator):
    gene = loader.get_gene_info("katG")
    display = tables.prepare_coordinates_data(
        loader.search_mutations_by_gene("katG"), gene, calculator
    )
    # katG_c.-1T>C sits one base past the gene end on the minus strand
    row = display[display["variant"] == "katG_c.-1T>C"].iloc[0]
    assert row["position_value"] == gene.end + 1
    assert row["Gene Relative"] == "c.-1"


# ----------------------------------------------------------------------
# TASK-14: summary tables
# ----------------------------------------------------------------------
def test_summary_tables_carry_a_totals_row(loader):
    summary = tables._with_totals(loader.get_catalogue_summary(), "Gene", "TOTAL")
    assert summary.iloc[-1]["Gene"] == "TOTAL"
    assert summary.iloc[-1]["Associations"] == summary.iloc[:-1]["Associations"].sum()


def test_summary_view_builds(loader):
    view = tables.catalogue_summary_view(loader)
    assert view is not None
