import gzip
import json
import os
import sys
import urllib.parse

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import pytest

import browser_tracks
import genome_view
from data_utils import DataLoader

DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'data'))


@pytest.fixture(scope="module")
def loader():
    instance = DataLoader(DATA_DIR)
    instance.load_gff3()
    return instance


@pytest.fixture(scope="module")
def tracks(loader, tmp_path_factory):
    out_dir = tmp_path_factory.mktemp("tracks")
    return browser_tracks.BrowserTracks(loader, tracks_dir=str(out_dir)).build()


def _features(tracks, name):
    with gzip.open(os.path.join(tracks.tracks_dir, name), "rt", encoding="utf-8") as handle:
        return [line.rstrip("\n").split("\t") for line in handle if not line.startswith("#")]


def _attributes(fields):
    return dict(chunk.split("=", 1) for chunk in fields[8].split(";"))


# ----------------------------------------------------------------------
# TASK-09: flattened prokaryotic annotation
# ----------------------------------------------------------------------
def test_gene_track_has_one_feature_per_locus_and_no_subfeatures(tracks, loader):
    """
    The browser must see one feature per locus.

    With no CDS or exon children, JBrowse offers no intron-based sequence
    options at all, which is the point of the flattening.
    """
    features = _features(tracks, "genes.gff3.gz")
    assert {fields[2] for fields in features} <= {"gene", "pseudogene"}
    assert not any("Parent=" in fields[8] for fields in features)
    assert len(features) == len(loader._genes_sorted)


def test_gene_track_escapes_attribute_separators(tracks):
    """A product containing a comma must not be split into two values."""
    with_comma = [
        _attributes(fields)["product"]
        for fields in _features(tracks, "genes.gff3.gz")
        if "%2C" in _attributes(fields).get("product", "")
    ]
    assert with_comma
    assert all("," not in product for product in with_comma)


def test_gene_track_carries_annotation_for_feature_details(tracks):
    """TASK-11: product, note, catalogue drugs and a Mycobrowser link."""
    katg = next(
        _attributes(fields) for fields in _features(tracks, "genes.gff3.gz")
        if _attributes(fields)["ID"] == "Rv1908c"
    )
    assert katg["Name"] == "katG"
    assert katg["product"] == "catalase-peroxidase"
    assert "Note" in katg
    assert "Isoniazid" in katg["who_catalogue_drugs"].split(",")
    assert "mycobrowser.epfl.ch/genes/Rv1908c" in urllib.parse.unquote(katg["mycobrowser"])


# ----------------------------------------------------------------------
# WHO catalogue tracks
# ----------------------------------------------------------------------
def test_variant_tracks_are_split_by_grading(tracks):
    for key, grades in {"assoc": {"1", "2"}, "uncertain": {"3"}, "not_assoc": {"4", "5"}}.items():
        features = _features(tracks, f"{browser_tracks.VARIANT_TRACKS[key]}.gff3.gz")
        assert features
        assert {_attributes(fields)["grade_group"] for fields in features} <= grades


def test_katg_s315t_is_resistance_associated_at_its_codon(tracks):
    features = _features(tracks, f"{browser_tracks.VARIANT_TRACKS['assoc']}.gff3.gz")
    fields = next(f for f in features if _attributes(f)["ID"] == "katG_p.Ser315Thr")
    attributes = _attributes(fields)
    assert attributes["grade_group"] == "1"
    assert "Isoniazid: 1) Assoc w R" in urllib.parse.unquote(attributes["drug_gradings"]).split(",")
    # Codon 315 of the minus-strand katG is 2,155,167-2,155,169; one of the
    # recorded alternative changes reaches into the neighbouring codon.
    assert int(fields[3]) == 2155167 and int(fields[4]) == 2155170
    assert fields[6] == "-"


def test_drug_track_holds_only_that_drugs_associated_variants(tracks):
    features = _features(tracks, f"{browser_tracks.drug_track_id('Isoniazid')}.gff3.gz")
    ids = {_attributes(fields)["variant"] for fields in features}
    assert "katG_p.Ser315Thr" in ids
    assert {_attributes(fields)["grade_group"] for fields in features} <= {"1", "2"}
    gradings = {
        grading
        for fields in features
        for grading in urllib.parse.unquote(_attributes(fields)["drug_gradings"]).split(",")
    }
    assert all(grading.startswith("Isoniazid:") for grading in gradings)


def test_density_counts_every_variant(tracks):
    with open(os.path.join(tracks.tracks_dir, "who_variant_density.bed"), encoding="utf-8") as handle:
        total = sum(int(line.split("\t")[4]) for line in handle)
    assert total == len(tracks.variant_features())


# ----------------------------------------------------------------------
# Text search index
# ----------------------------------------------------------------------
def _search(tracks, term):
    """Resolve a term the way JBrowse's Trix adapter does."""
    ix_path = os.path.join(tracks.tracks_dir, "search.ix")
    with open(os.path.join(tracks.tracks_dir, "search.ixx"), encoding="utf-8") as handle:
        index = [(line[:-11], int(line[-11:-1], 16)) for line in handle if line.strip()]
    offset = 0
    for prefix, block_offset in index:
        if prefix[:len(term)] < term:
            offset = block_offset
    with open(ix_path, "rb") as handle:
        handle.seek(offset)
        for raw in handle:
            word, *docs = raw.decode("utf-8").rstrip("\n").split(" ")
            if word == term:
                return [
                    [urllib.parse.unquote(v) for v in json.loads(doc.replace("|", ","))]
                    for doc in docs
                ]
            if word > term:
                return []
    return []


def test_text_index_blocks_start_at_prefix_changes(tracks):
    """Every term sharing a prefix must be reachable from that prefix's block."""
    with open(os.path.join(tracks.tracks_dir, "search.ix"), "rb") as handle:
        raw = handle.read()
    with open(os.path.join(tracks.tracks_dir, "search.ixx"), encoding="utf-8") as handle:
        for line in handle:
            offset = int(line[-11:-1], 16)
            if offset:
                before = raw[:offset - 1].rsplit(b"\n", 1)[-1].split(b" ", 1)[0]
                after = raw[offset:].split(b" ", 1)[0]
                assert before[:5] != after[:5]


def test_text_index_is_sorted_with_valid_offsets(tracks):
    with open(os.path.join(tracks.tracks_dir, "search.ix"), "rb") as handle:
        raw = handle.read()
    words = [line.split(b" ", 1)[0] for line in raw.splitlines()]
    assert words == sorted(words)
    with open(os.path.join(tracks.tracks_dir, "search.ixx"), encoding="utf-8") as handle:
        for line in handle:
            offset = int(line[-11:-1], 16)
            assert offset == 0 or raw[offset - 1:offset] == b"\n"


def test_aliased_variant_carries_its_alias(tracks):
    """#5: a fabG1 promoter variant filed under inhA shows its fabG1 name."""
    features = _features(tracks, f"{browser_tracks.VARIANT_TRACKS['assoc']}.gff3.gz")
    fields = next(f for f in features if _attributes(f)["ID"] == "inhA_c.-777C>T")
    attributes = _attributes(fields)
    assert urllib.parse.unquote(attributes["alias"]) == "fabG1_c.-15C>T"
    assert attributes["gene"] == "inhA"
    # Recorded with equivalent multi-base changes, so the span starts upstream.
    assert int(fields[3]) <= 1673425 <= int(fields[4])


def test_unaliased_variant_has_no_alias_attribute(tracks):
    features = _features(tracks, f"{browser_tracks.VARIANT_TRACKS['assoc']}.gff3.gz")
    fields = next(f for f in features if _attributes(f)["ID"] == "katG_p.Ser315Thr")
    assert "alias" not in _attributes(fields)


@pytest.mark.parametrize("term", ["katg", "rv1908c", "catalase-peroxidase"])
def test_text_index_finds_genes(tracks, term):
    hits = _search(tracks, term)
    assert any(hit[1] == browser_tracks.GENES_TRACK and hit[2] == "katG" for hit in hits)
    assert hits[0][0].startswith("NC_000962.3:")


def test_text_index_finds_variants_on_their_track(tracks):
    hits = _search(tracks, "katg_p.ser315thr")
    assert hits and hits[0][1] == browser_tracks.VARIANT_TRACKS["assoc"]
    assert hits[0][2] == "katG_p.Ser315Thr"


def test_text_index_finds_variants_by_alias(tracks):
    hits = _search(tracks, "fabg1_c.-15c>t")
    assert hits and hits[0][2] == "inhA_c.-777C>T"


# ----------------------------------------------------------------------
# Browser configuration
# ----------------------------------------------------------------------
def test_session_only_references_configured_tracks(tracks):
    config = tracks.config(
        visible_tracks=[choice["value"] for choice in browser_tracks.TRACK_CHOICES]
        + [browser_tracks.drug_track_id(drug) for drug in tracks.drugs]
    )
    track_ids = [track["trackId"] for track in config["tracks"]]
    assert len(track_ids) == len(set(track_ids))

    known = set(track_ids) | {browser_tracks.SEQUENCE_TRACK}
    session = config["defaultSession"]["view"]["tracks"]
    assert {track["configuration"] for track in session} <= known
    # Every toggle, the reference and each drug track are shown once.
    assert len(session) == len(browser_tracks.TRACK_CHOICES) + 1 + len(tracks.drugs)


def test_hidden_tracks_stay_in_the_track_selector(tracks):
    config = tracks.config(visible_tracks=[browser_tracks.GENES_TRACK])
    session = {track["configuration"] for track in config["defaultSession"]["view"]["tracks"]}
    assert session == {browser_tracks.GENES_TRACK, browser_tracks.SEQUENCE_TRACK}
    assert len(config["tracks"]) > len(session)


def test_selection_track_marks_gene_and_mutation(loader):
    gene = loader.get_gene_info("katG")
    selection = {"variant": "katG_p.Ser315Thr", "position": 2155168,
                 "reference": "C", "alternative": "G"}
    gene_feature, variant_feature = browser_tracks.selection_features(gene, selection)
    assert (gene_feature["start"], gene_feature["end"]) == (gene.start - 1, gene.end)
    assert (variant_feature["start"], variant_feature["end"]) == (2155167, 2155168)
    assert variant_feature["name"] == "katG_p.Ser315Thr"


def test_focus_location(loader):
    gene = loader.get_gene_info("katG")
    selection = {"position": 2155168}
    assert genome_view.focus_location(gene, "2000") == (
        f"NC_000962.3:{gene.start - 2000}..{gene.end + 2000}"
    )
    assert genome_view.focus_location(gene, genome_view.FOCUS_MUTATION, selection) == (
        "NC_000962.3:2155128..2155208"
    )
    # Without coordinates the mutation focus falls back to the gene.
    assert genome_view.focus_location(gene, genome_view.FOCUS_MUTATION, {}) == (
        genome_view.focus_location(gene, genome_view.FOCUS_GENE)
    )
