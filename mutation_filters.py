"""
Mutation filters for the Drug Resistance Profile and Genomic Coordinates.

Replaces the DataTable native filter row, which expected query syntax typed
into an unlabelled cell under each header, with labelled controls: a free-text
mutation search, drop-down lists built from the values actually present for
the active gene (with their counts), a position range and a few on/off flags.

Filtering runs server-side, so every per-drug table, the per-drug counts and
the coordinates table stay in step with one set of filters.
"""
import re
from typing import Dict, Iterable, List, Optional

import dash_bootstrap_components as dbc
import pandas as pd
from dash import dcc, html

import search_utils
from data_utils import CONFIDENCE_GRADES

# Every control carries a pattern-matching id so that one callback can read
# them all with ``ALL``, and so that callbacks fire safely for genes whose
# panels are absent.
RESISTANCE_FILTER = "resistance-filter"
COORDS_FILTER = "coords-filter"
COORDS_FILTER_SUMMARY_ID = {"type": "coords-filter-summary", "index": "main"}


def resistance_filter_id(field: str) -> Dict[str, str]:
    return {"type": RESISTANCE_FILTER, "field": field}


def coords_filter_id(field: str) -> Dict[str, str]:
    return {"type": COORDS_FILTER, "field": field}


# Value of each resistance filter when nothing is filtered.
RESISTANCE_DEFAULTS: Dict = {
    "query": "",
    "drugs": [],
    "grades": [],
    "effects": [],
    "levels": [],
    "tiers": [],
    "changes": [],
    "pos_min": None,
    "pos_max": None,
    "flags": [],
}

# Filters folded away under "More filters"; counted on its button so a hidden
# filter is never forgotten.
ADVANCED_FIELDS = ("effects", "levels", "tiers", "changes", "flags")

COORDS_DEFAULTS: Dict = {"query": "", "changes": [], "linked": []}

ASSOCIATED_GRADES = list(CONFIDENCE_GRADES[:2])
UNCERTAIN_GRADES = [CONFIDENCE_GRADES[2]]
NOT_ASSOCIATED_GRADES = list(CONFIDENCE_GRADES[3:])

GRADE_LABELS = {
    "1) Assoc w R": "1) Associated with resistance",
    "2) Assoc w R - Interim": "2) Associated with resistance – interim",
    "3) Uncertain significance": "3) Uncertain significance",
    "4) Not assoc w R - Interim": "4) Not associated – interim",
    "5) Not assoc w R": "5) Not associated with resistance",
}

# Notation level of a catalogue mutation, from its HGVS prefix.
MUTATION_LEVELS = {
    "protein": "Protein change (p.)",
    "coding": "Coding or promoter nucleotide (c.)",
    "noncoding": "rRNA / non-coding nucleotide (n.)",
    "gene": "Gene-level (LoF, deletion)",
}

EFFECT_LABELS = {
    "upstream_gene_variant": "Upstream gene variant / promoter",
    "non_coding_transcript_exon_variant": "Non-coding transcript exon variant",
    "LoF": "Loss of function, LoF",
    "feature_ablation": "Feature ablation / gene deletion",
}

FLAG_LABELS = {
    "comment": "Has a WHO comment",
    "relaxed": "Flagged in the relaxed-thresholds simulation",
    "hide_silent": "Hide silent mutations",
}

CHANGE_KINDS = ("SNV", "Insertion", "Deletion", "MNV")

QUERY_PLACEHOLDER = "e.g. S450L, Ser450Leu, 450, c.-15C>T, fs"

_POSITION_REGEX = re.compile(r"^[pcn]\.\D*?(-?\d+)", re.IGNORECASE)
_QUERY_SPLIT = re.compile(r"[,;\s]+")


# ----------------------------------------------------------------------
# Mutation parsing
# ----------------------------------------------------------------------
def mutation_level(mutation: str) -> str:
    """Notation level of a catalogue mutation: protein, coding, noncoding or gene."""
    text = str(mutation or "").lower()
    if text.startswith("p."):
        return "protein"
    if text.startswith("c."):
        return "coding"
    if text.startswith("n."):
        return "noncoding"
    return "gene"


def mutation_position(mutation: str) -> Optional[int]:
    """
    First position in a mutation: the codon of a ``p.`` change, the nucleotide
    of a ``c.`` or ``n.`` change (negative upstream of the start codon).

    ``p.Ser450Leu`` -> 450, ``c.-15C>T`` -> -15, ``p.Met97_Ser98dup`` -> 97.
    Gene-level entries such as ``LoF`` have no position.
    """
    match = _POSITION_REGEX.match(str(mutation or ""))
    return int(match.group(1)) if match else None


def _term_pattern(term: str) -> Optional[re.Pattern]:
    """
    Pattern for one search term, accepting one- or three-letter residues.

    Numbers are matched whole, so ``450`` finds ``p.Ser450Leu`` but not
    ``p.Ser1450Leu`` or ``p.Ser45Leu``.
    """
    term = term.strip()
    if "_" in term:
        # A full variant such as "rpoB_S450L": the gene is already chosen.
        term = term.split("_", 1)[1]
    if not term:
        return None

    candidates = {term.lower(), search_utils.normalize_mutation(term).lower()}
    alternatives = []
    for candidate in sorted(candidates):
        if not candidate:
            continue
        pattern = re.escape(candidate)
        if candidate[0].isdigit():
            pattern = r"(?<!\d)" + pattern
        if candidate[-1].isdigit():
            pattern += r"(?!\d)"
        alternatives.append(pattern)
    return re.compile("|".join(alternatives)) if alternatives else None


def query_mask(values: pd.Series, query: Optional[str]) -> pd.Series:
    """Rows whose text matches any of the comma- or space-separated terms."""
    patterns = [
        pattern for pattern in (_term_pattern(term) for term in _QUERY_SPLIT.split(query or ""))
        if pattern is not None
    ]
    if not patterns:
        return pd.Series(True, index=values.index)

    lowered = values.fillna("").astype(str).str.lower()
    return lowered.apply(lambda text: any(p.search(text) for p in patterns))


# ----------------------------------------------------------------------
# Filter state
# ----------------------------------------------------------------------
def collect(inputs_list: List[Dict], defaults: Dict) -> Dict:
    """Filter state from the ``ALL`` inputs of a pattern-matching callback."""
    state = dict(defaults)
    for item in inputs_list:
        field = item["id"]["field"]
        if field in state:
            value = item.get("value")
            state[field] = defaults[field] if value is None else value
    return state


def _is_set(value) -> bool:
    return value not in (None, "", [])


def active_count(filters: Dict, fields: Iterable[str]) -> int:
    """Number of filters in use among ``fields``; each flag counts once."""
    count = 0
    for field in fields:
        value = filters.get(field)
        if field == "flags":
            count += len(value or [])
        elif _is_set(value):
            count += 1
    return count


def advanced_count(filters: Dict) -> int:
    """Active filters folded away under "More filters", the range counting once."""
    has_range = _is_set(filters.get("pos_min")) or _is_set(filters.get("pos_max"))
    return active_count(filters, ADVANCED_FIELDS) + int(has_range)


def is_filtered(filters: Dict, defaults: Dict = RESISTANCE_DEFAULTS) -> bool:
    return any(
        (str(filters.get(field)).strip() if field == "query" else filters.get(field))
        not in (None, "", [])
        for field in defaults
    )


# ----------------------------------------------------------------------
# Drug Resistance Profile
# ----------------------------------------------------------------------
def filter_resistance(prepared: pd.DataFrame, filters: Dict) -> pd.DataFrame:
    """
    Rows of a prepared resistance table that pass the filters.

    The drug filter is applied per table by the caller, since each drug has its
    own table.
    """
    if prepared.empty:
        return prepared

    mask = query_mask(prepared["mutation"], filters.get("query"))

    for field, column in (("grades", "confidence"), ("effects", "effect"),
                          ("changes", "changes_vs_ver1")):
        if filters.get(field):
            mask &= prepared[column].isin(filters[field])

    if filters.get("tiers"):
        mask &= prepared["tier"].astype(str).isin([str(t) for t in filters["tiers"]])

    if filters.get("levels"):
        mask &= prepared["mutation"].apply(mutation_level).isin(filters["levels"])

    low, high = filters.get("pos_min"), filters.get("pos_max")
    if _is_set(low) or _is_set(high):
        positions = prepared["mutation"].apply(mutation_position)
        in_range = positions.notna()
        if _is_set(low):
            in_range &= positions >= float(low)
        if _is_set(high):
            in_range &= positions <= float(high)
        mask &= in_range

    flags = set(filters.get("flags") or [])
    if "comment" in flags:
        mask &= prepared["comment"].astype(str).str.strip() != ""
    if "relaxed" in flags:
        mask &= prepared["relaxed_thresholds"].astype(str).str.strip() != ""
    if "hide_silent" in flags:
        mask &= prepared["silent_mutation"] != "Yes"

    return prepared[mask].reset_index(drop=True)


def filter_by_drug(prepared_by_drug: Dict[str, pd.DataFrame], filters: Dict) -> Dict[str, pd.DataFrame]:
    """Filtered rows for every drug; deselected drugs come back empty."""
    drugs = set(filters.get("drugs") or prepared_by_drug)
    return {
        drug: filter_resistance(prepared, filters) if drug in drugs else prepared.iloc[0:0]
        for drug, prepared in prepared_by_drug.items()
    }


def _counted_options(values: pd.Series, labels: Optional[Dict] = None,
                     order: Optional[Iterable] = None) -> List[Dict]:
    """Drop-down options for the values present, each with its row count."""
    counts = values[values.astype(str).str.strip() != ""].value_counts()
    keys = [key for key in order if key in counts.index] if order is not None else list(counts.index)
    labels = labels or {}
    return [
        {"label": f"{labels.get(key, key)} ({counts[key]:,})", "value": key}
        for key in keys
    ]


def _effect_label(effect: str) -> str:
    return EFFECT_LABELS.get(effect, effect.replace("_", " ").capitalize())


def resistance_filter_options(prepared_by_drug: Dict[str, pd.DataFrame]) -> Dict[str, List[Dict]]:
    """Options for each drop-down, drawn from the active gene's catalogue rows."""
    frames = [frame.assign(drug=drug) for drug, frame in prepared_by_drug.items()]
    rows = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(
        columns=["drug", "mutation", "confidence", "effect", "tier", "changes_vs_ver1"]
    )

    effects = rows["effect"]
    return {
        "drugs": _counted_options(rows["drug"], order=list(prepared_by_drug)),
        "grades": _counted_options(rows["confidence"], GRADE_LABELS, order=CONFIDENCE_GRADES),
        "effects": _counted_options(effects, {e: _effect_label(e) for e in effects.unique()}),
        "levels": _counted_options(
            rows["mutation"].apply(mutation_level), MUTATION_LEVELS, order=MUTATION_LEVELS
        ),
        "tiers": _counted_options(
            rows["tier"].astype(str),
            {t: f"Tier {t}" for t in rows["tier"].astype(str).unique()},
            order=sorted(rows["tier"].astype(str).unique()),
        ),
        "changes": _counted_options(rows["changes_vs_ver1"]),
        "flags": _flag_options(rows),
    }


def _flag_options(rows: pd.DataFrame) -> List[Dict]:
    """
    On/off flags with the number of rows carrying each one; a flag no row of
    this gene carries is disabled rather than hidden.
    """
    def present(column: str) -> int:
        if column not in rows:
            return 0
        values = rows[column].astype(str).str.strip()
        return int(((values != "") & (values != "nan")).sum())

    counts = {
        "comment": present("comment"),
        "relaxed": present("relaxed_thresholds"),
        "hide_silent": present("silent_mutation"),
    }
    return [
        {"label": f"{label} ({counts[value]:,})", "value": value, "disabled": counts[value] == 0}
        for value, label in FLAG_LABELS.items()
    ]


def _labelled(label: str, control, help_text: Optional[str] = None) -> html.Div:
    return html.Div([
        html.Label(label, className="filter-label"),
        control,
        html.Div(help_text, className="filter-help") if help_text else None,
    ], className="filter-field")


def _dropdown(field: str, options: List[Dict], placeholder: str) -> dcc.Dropdown:
    return dcc.Dropdown(
        id=resistance_filter_id(field),
        options=options,
        value=[],
        multi=True,
        placeholder=placeholder,
        className="filter-dropdown",
    )


def resistance_filter_panel(prepared_by_drug: Dict[str, pd.DataFrame]) -> html.Div:
    """Filter controls shown above the per-drug resistance tables."""
    options = resistance_filter_options(prepared_by_drug)
    total = sum(len(frame) for frame in prepared_by_drug.values())

    search = dbc.InputGroup([
        dbc.InputGroupText(html.I(className="bi bi-search")),
        dbc.Input(
            id=resistance_filter_id("query"),
            type="search",
            value="",
            placeholder=QUERY_PLACEHOLDER,
            debounce=350,
            autoComplete="off",
            className="filter-search",
        ),
    ], size="sm")

    grade_presets = html.Div([
        html.Span("Quick pick:", className="filter-help me-1"),
        dbc.Button("Associated (1–2)", id="filter-preset-associated",
                   className="filter-preset", size="sm"),
        dbc.Button("Uncertain (3)", id="filter-preset-uncertain",
                   className="filter-preset", size="sm"),
        dbc.Button("Not associated (4–5)", id="filter-preset-not-associated",
                   className="filter-preset", size="sm"),
    ], className="d-flex flex-wrap align-items-center gap-1 mt-1")

    primary = dbc.Row([
        dbc.Col(_labelled(
            "Find mutation", search,
            "One- or three-letter residues, a codon or nucleotide position, or "
            "several terms separated by commas.",
        ), lg=5),
        dbc.Col(_labelled(
            "Drug", _dropdown("drugs", options["drugs"], "All drugs")
        ), lg=3, md=6),
        dbc.Col([
            _labelled(
                "Confidence grading",
                _dropdown("grades", options["grades"], "All gradings")
            ),
            grade_presets,
        ], lg=4, md=6),
    ], className="g-3")

    position_range = dbc.InputGroup([
        dbc.Input(id=resistance_filter_id("pos_min"), type="number", placeholder="From",
                  debounce=True, step=1),
        dbc.InputGroupText("–"),
        dbc.Input(id=resistance_filter_id("pos_max"), type="number", placeholder="To",
                  debounce=True, step=1),
    ], size="sm")

    advanced = dbc.Collapse(
        html.Div([
            dbc.Row([
                dbc.Col(_labelled(
                    "Effect", _dropdown("effects", options["effects"], "Any effect")
                ), lg=4, md=6),
                dbc.Col(_labelled(
                    "Mutation type", _dropdown("levels", options["levels"], "Any type")
                ), lg=4, md=6),
                dbc.Col(_labelled(
                    "Change vs catalogue v1",
                    _dropdown("changes", options["changes"], "Any change")
                ), lg=4, md=6),
                dbc.Col(_labelled(
                    "Position range", position_range,
                    "Codon for p. changes; nucleotide for c./n. changes "
                    "(negative = upstream). Excludes gene-level entries.",
                ), lg=4, md=6),
                dbc.Col(_labelled(
                    "Gene tier", _dropdown("tiers", options["tiers"], "Any tier")
                ), lg=2, md=6),
                dbc.Col(_labelled(
                    "Flags",
                    dbc.Checklist(
                        id=resistance_filter_id("flags"),
                        options=options["flags"],
                        value=[],
                        switch=True,
                        className="filter-flags small",
                    )
                ), lg=6),
            ], className="g-3"),
        ], className="filter-advanced"),
        id="filter-advanced-collapse",
        is_open=False,
    )

    toolbar = html.Div([
        dbc.Button([
            html.I(className="bi bi-sliders me-2"),
            "More filters",
            dbc.Badge("", id="filter-advanced-count", className="filter-count ms-2",
                      style={"display": "none"}),
            html.I(className="bi bi-chevron-down ms-2", id="filter-advanced-chevron"),
        ], id="filter-advanced-toggle", className="btn-soft", size="sm"),
        html.Div(
            resistance_summary(total, total, len(prepared_by_drug), len(prepared_by_drug), False),
            id="resistance-filter-summary",
            className="filter-summary",
            role="status",
            **{"aria-live": "polite"},
        ),
        dbc.Button(
            [html.I(className="bi bi-x-circle me-1"), "Clear filters"],
            id="filter-clear",
            className="btn-soft",
            size="sm",
            disabled=True,
        ),
    ], className="filter-toolbar")

    return html.Div([
        html.Div([
            html.I(className="bi bi-funnel me-2"),
            "Filter mutations",
        ], className="filter-title"),
        primary,
        toolbar,
        advanced,
    ], className="filter-panel")


def resistance_summary(shown: int, total: int, drugs_shown: int, drugs_total: int,
                       filtered: bool) -> List:
    """Plain-language result line under the filters."""
    drug_word = "drug" if drugs_total == 1 else "drugs"
    if not filtered:
        return [
            "Showing all ", html.Strong(f"{total:,}"),
            f" catalogue entries across {drugs_total} {drug_word}.",
        ]
    return [
        "Showing ", html.Strong(f"{shown:,}"), f" of {total:,} catalogue entries in ",
        html.Strong(f"{drugs_shown}"), f" of {drugs_total} {drug_word}.",
    ]


def no_matches_alert() -> html.Div:
    return html.Div([
        html.I(className="bi bi-search me-2"),
        "No mutations match these filters. Loosen a filter or use ",
        html.Strong("Clear filters"),
        ".",
    ], className="filter-empty")


# ----------------------------------------------------------------------
# Genomic Coordinates
# ----------------------------------------------------------------------
def change_kind(reference: str, alternative: str) -> str:
    """SNV, insertion, deletion or MNV, from the allele lengths."""
    ref_len, alt_len = len(reference or ""), len(alternative or "")
    if ref_len == alt_len == 1:
        return "SNV"
    if ref_len > alt_len:
        return "Deletion"
    if alt_len > ref_len:
        return "Insertion"
    return "MNV"


def filter_coordinates(display: pd.DataFrame, filters: Dict,
                       variants: Optional[Iterable[str]] = None) -> pd.DataFrame:
    """
    Rows of the prepared coordinates table that pass the filters.

    ``variants``, when given, keeps only those catalogue variants: it links the
    table to the Drug Resistance Profile filters.
    """
    if display.empty:
        return display

    searchable = display["variant"].fillna("").astype(str) + " " + \
        display["position"].fillna("").astype(str)
    mask = query_mask(searchable, filters.get("query"))

    if filters.get("changes"):
        mask &= display["change_kind"].isin(filters["changes"])
    if variants is not None:
        mask &= display["variant"].isin(set(variants))

    return display[mask].reset_index(drop=True)


def coordinates_filter_bar(display: pd.DataFrame, linkable: bool) -> html.Div:
    """Compact filters above the Genomic Coordinates table."""
    counts = display["change_kind"].value_counts() if not display.empty else pd.Series(dtype=int)
    change_options = [
        {"label": f"{kind} ({counts[kind]:,})", "value": kind}
        for kind in CHANGE_KINDS if kind in counts.index
    ]

    columns = [
        dbc.Col(_labelled(
            "Find variant or position",
            dbc.InputGroup([
                dbc.InputGroupText(html.I(className="bi bi-search")),
                dbc.Input(
                    id=coords_filter_id("query"),
                    type="search",
                    value="",
                    placeholder="e.g. S450L, 761155, c.-15C>T",
                    debounce=350,
                    autoComplete="off",
                ),
            ], size="sm")
        ), lg=5),
        dbc.Col(_labelled(
            "Change type",
            dcc.Dropdown(
                id=coords_filter_id("changes"),
                options=change_options,
                value=[],
                multi=True,
                placeholder="Any change type",
                className="filter-dropdown",
            )
        ), lg=3, md=6),
    ]
    if linkable:
        columns.append(dbc.Col(_labelled(
            "Linked filters",
            dbc.Checklist(
                id=coords_filter_id("linked"),
                options=[{"label": "Match the Drug Resistance Profile filters",
                          "value": "linked"}],
                value=[],
                switch=True,
                className="filter-flags small",
            )
        ), lg=4, md=6))

    return html.Div([
        dbc.Row(columns, className="g-3"),
        html.Div(
            coordinates_summary(len(display), len(display), False),
            id=dict(COORDS_FILTER_SUMMARY_ID),
            className="filter-summary mt-2",
            role="status",
            **{"aria-live": "polite"},
        ),
    ], className="filter-panel filter-panel-compact")


def coordinates_summary(shown: int, total: int, filtered: bool) -> List:
    if not filtered:
        return ["Showing all ", html.Strong(f"{total:,}"), " nucleotide changes."]
    return ["Showing ", html.Strong(f"{shown:,}"), f" of {total:,} nucleotide changes."]
