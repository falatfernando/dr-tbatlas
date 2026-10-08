import os
import sys

# Ensure parent directory is in path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import pandas as pd
import pytest

import mutation_filters as mf
import tables
from coordinate_calculator import CoordinateCalculator
from data_utils import DataLoader

DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'data'))


@pytest.fixture(scope="module")
def loader():
    instance = DataLoader(DATA_DIR)
    instance.load_gff3()
    return instance


@pytest.fixture(scope="module")
def rpob(loader):
    return tables.prepare_resistance_data(loader.get_drug_resistance_info("rpoB"))


@pytest.fixture(scope="module")
def rv0678_by_drug(loader):
    rows = loader.get_drug_resistance_info("Rv0678")
    return {
        drug: tables.prepare_resistance_data(rows[rows["drug"] == drug])
        for drug in sorted(rows["drug"].unique())
    }


def filters(**overrides):
    return dict(mf.RESISTANCE_DEFAULTS, **overrides)


# ----------------------------------------------------------------------
# Mutation parsing
# ----------------------------------------------------------------------
@pytest.mark.parametrize("mutation,level", [
    ("p.Ser450Leu", "protein"),
    ("c.-15C>T", "coding"),
    ("n.1401A>G", "noncoding"),
    ("LoF", "gene"),
    ("deletion", "gene"),
])
def test_mutation_level(mutation, level):
    assert mf.mutation_level(mutation) == level


@pytest.mark.parametrize("mutation,position", [
    ("p.Ser450Leu", 450),
    ("c.-15C>T", -15),
    ("n.1401A>G", 1401),
    ("p.Met97_Ser98dup", 97),
    ("LoF", None),
])
def test_mutation_position(mutation, position):
    assert mf.mutation_position(mutation) == position


# ----------------------------------------------------------------------
# Free-text search
# ----------------------------------------------------------------------
MUTATIONS = pd.Series([
    "p.Ser450Leu", "p.Ser1450Leu", "p.Ser45Leu", "p.His445Tyr", "c.-15C>T", "p.Thr400fs", "LoF",
])


def _matches(query):
    return MUTATIONS[mf.query_mask(MUTATIONS, query)].tolist()


@pytest.mark.parametrize("query", ["S450L", "Ser450Leu", "p.Ser450Leu", "ser450leu", "rpoB_S450L"])
def test_search_accepts_one_and_three_letter_residues(query):
    assert _matches(query) == ["p.Ser450Leu"]


def test_search_matches_positions_whole():
    assert _matches("450") == ["p.Ser450Leu"]
    assert _matches("45") == ["p.Ser45Leu"]


def test_search_accepts_several_terms():
    assert _matches("S450L, H445Y") == ["p.Ser450Leu", "p.His445Tyr"]


def test_search_finds_nucleotide_and_consequence_terms():
    assert _matches("c.-15C>T") == ["c.-15C>T"]
    assert _matches("fs") == ["p.Thr400fs"]
    assert _matches("lof") == ["LoF"]


def test_blank_search_keeps_everything():
    assert _matches("  ") == MUTATIONS.tolist()


# ----------------------------------------------------------------------
# Drug Resistance Profile
# ----------------------------------------------------------------------
def test_no_filters_keep_every_row(rpob):
    assert len(mf.filter_resistance(rpob, filters())) == len(rpob)


def test_grading_filter(rpob):
    kept = mf.filter_resistance(rpob, filters(grades=mf.ASSOCIATED_GRADES))
    assert len(kept) > 0
    assert set(kept["confidence"]) <= set(mf.ASSOCIATED_GRADES)


def test_rrdr_codon_range_of_protein_changes(rpob):
    """The rifampicin resistance-determining region, codons 426-452."""
    kept = mf.filter_resistance(rpob, filters(levels=["protein"], pos_min=426, pos_max=452))
    positions = kept["mutation"].apply(mf.mutation_position)
    assert len(kept) > 0
    assert positions.between(426, 452).all()
    assert "p.Ser450Leu" in kept["mutation"].tolist()


def test_position_range_drops_gene_level_entries(rv0678_by_drug):
    rows = next(iter(rv0678_by_drug.values()))
    kept = mf.filter_resistance(rows, filters(pos_min=-1000))
    assert (kept["mutation"].apply(mf.mutation_level) != "gene").all()


def test_flags(rv0678_by_drug):
    rows = next(iter(rv0678_by_drug.values()))
    with_comment = mf.filter_resistance(rows, filters(flags=["comment"]))
    assert len(with_comment) > 0
    assert (with_comment["comment"] != "").all()

    no_silent = mf.filter_resistance(rows, filters(flags=["hide_silent"]))
    assert (no_silent["silent_mutation"] != "Yes").all()


def test_drug_filter_empties_the_other_drugs(rv0678_by_drug):
    first, second = list(rv0678_by_drug)[:2]
    filtered = mf.filter_by_drug(rv0678_by_drug, filters(drugs=[first]))
    assert len(filtered[first]) == len(rv0678_by_drug[first])
    assert filtered[second].empty


def test_options_carry_counts_and_follow_catalogue_order(rv0678_by_drug):
    options = mf.resistance_filter_options(rv0678_by_drug)
    assert [o["value"] for o in options["drugs"]] == list(rv0678_by_drug)

    grades = [o["value"] for o in options["grades"]]
    assert grades == sorted(grades)
    assert all(o["label"].endswith(")") for o in options["effects"])


def test_flags_absent_from_a_gene_are_disabled(rpob):
    options = mf.resistance_filter_options({"Rifampicin": rpob})
    by_value = {o["value"]: o for o in options["flags"]}
    assert by_value["comment"]["disabled"] == (rpob["comment"] == "").all()
    assert not by_value["hide_silent"]["disabled"]


def test_filter_state_is_read_from_pattern_matching_inputs():
    inputs = [
        {"id": {"type": mf.RESISTANCE_FILTER, "field": "query"}, "value": "S450L"},
        {"id": {"type": mf.RESISTANCE_FILTER, "field": "drugs"}, "value": None},
    ]
    state = mf.collect(inputs, mf.RESISTANCE_DEFAULTS)
    assert state["query"] == "S450L"
    assert state["drugs"] == []
    assert mf.is_filtered(state)
    assert not mf.is_filtered(mf.collect([], mf.RESISTANCE_DEFAULTS))


def test_advanced_count_includes_hidden_filters_only():
    assert mf.advanced_count(filters(query="S450L", grades=mf.ASSOCIATED_GRADES)) == 0
    assert mf.advanced_count(filters(effects=["missense_variant"], pos_min=1, pos_max=9)) == 2
    assert mf.advanced_count(filters(flags=["comment", "hide_silent"])) == 2


# ----------------------------------------------------------------------
# Genomic Coordinates
# ----------------------------------------------------------------------
@pytest.mark.parametrize("ref,alt,kind", [
    ("C", "T", "SNV"),
    ("CA", "C", "Deletion"),
    ("C", "CA", "Insertion"),
    ("CAG", "TAC", "MNV"),
])
def test_change_kind(ref, alt, kind):
    assert mf.change_kind(ref, alt) == kind


@pytest.fixture(scope="module")
def katg_coordinates(loader):
    gene = loader.get_gene_info("katG")
    return tables.prepare_coordinates_data(
        loader.search_mutations_by_gene("katG"), gene, CoordinateCalculator(loader)
    )


def test_coordinates_filter_by_change_type(katg_coordinates):
    kept = mf.filter_coordinates(katg_coordinates, dict(mf.COORDS_DEFAULTS, changes=["Deletion"]))
    assert len(kept) > 0
    assert set(kept["change_kind"]) == {"Deletion"}


def test_coordinates_search_by_variant_or_position(katg_coordinates):
    by_variant = mf.filter_coordinates(katg_coordinates, dict(mf.COORDS_DEFAULTS, query="S315T"))
    assert set(by_variant["variant"]) == {"katG_p.Ser315Thr"}

    position = str(by_variant["position"].iloc[0])
    by_position = mf.filter_coordinates(katg_coordinates, dict(mf.COORDS_DEFAULTS, query=position))
    assert position in by_position["position"].astype(str).tolist()


def test_coordinates_linked_to_resistance_variants(katg_coordinates):
    kept = mf.filter_coordinates(katg_coordinates, mf.COORDS_DEFAULTS, variants={"katG_p.Ser315Thr"})
    assert set(kept["variant"]) == {"katG_p.Ser315Thr"}


def test_coordinate_tooltips_follow_filtered_rows(katg_coordinates):
    kept = mf.filter_coordinates(katg_coordinates, dict(mf.COORDS_DEFAULTS, changes=["Deletion"]))
    assert len(tables.coordinates_tooltips(kept)) == len(kept)
