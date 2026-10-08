"""
Table builders for the DR-TBAtlas.

Covers the Drug Resistance Profile with the full WHO catalogue schema
(TASK-12), the Genomic Coordinates table with horizontal overflow and
copyable sequence cells (TASK-07), and the catalogue summary tables
(TASK-14).
"""
from typing import Dict, List, Optional

import pandas as pd
from dash import dash_table, dcc, html
import dash_bootstrap_components as dbc

from data_utils import CONFIDENCE_GRADES, RELAXED_THRESHOLDS_COLUMN, GeneInfo
from mutation_filters import change_kind

# Longest sequence rendered inline in a table cell before truncation.
SEQ_DISPLAY_CHARS = 28
# Longest sequence rendered inside a hover tooltip.
SEQ_TOOLTIP_CHARS = 500

RESISTANCE_PAGE_SIZE = 8

# The coordinates table is created only once a gene is loaded. Dash resolves a
# plain string id eagerly whenever a callback that references it fires, so the
# dynamically created components use pattern-matching ids, which resolve to an
# empty list while they are absent.
COORDS_TABLE_ID = {"type": "coords-table", "index": "main"}

HIGHLIGHT_COLOR = "rgba(67, 142, 142, 0.18)"
HIGHLIGHT_BORDER = "#438E8E"

# Columns of the Drug Resistance Profile. Optional ones can be toggled off by
# the user; the mutation identity and grading are always shown.
RESISTANCE_COLUMNS: List[Dict] = [
    {"key": "mutation", "name": "Mutation", "source": "mutation", "required": True},
    {"key": "tier", "name": "Tier", "source": "tier", "required": True},
    {"key": "confidence", "name": "Confidence", "source": "FINAL CONFIDENCE GRADING", "required": True},
    {"key": "effect", "name": "Effect", "source": "effect"},
    {"key": "comment", "name": "Comment", "source": "Comment"},
    {"key": "changes_vs_ver1", "name": "CHANGES vs ver1", "source": "CHANGES vs ver1"},
    {"key": "relaxed_thresholds", "name": "Relaxed thresholds", "source": RELAXED_THRESHOLDS_COLUMN},
    {"key": "silent_mutation", "name": "Silent mutation", "source": "Silent mutation"},
]

OPTIONAL_RESISTANCE_COLUMNS = [c for c in RESISTANCE_COLUMNS if not c.get("required")]
DEFAULT_VISIBLE_COLUMNS = [c["key"] for c in OPTIONAL_RESISTANCE_COLUMNS]

_COLUMN_TOOLTIPS = {
    "tier": "WHO gene tier: 1 = established resistance gene, 2 = candidate gene.",
    "effect": "Predicted molecular consequence of the variant.",
    "confidence": "WHO final confidence grading for the drug/mutation pair.",
    "comment": "WHO remarks qualifying how the grading should be interpreted.",
    "changes_vs_ver1": "How this entry changed relative to catalogue version 1.",
    "relaxed_thresholds": RELAXED_THRESHOLDS_COLUMN,
    "silent_mutation": "Flagged by WHO as a silent (synonymous) mutation.",
}

_BASE_CELL_STYLE = {
    "textAlign": "left",
    "padding": "10px 12px",
    "fontFamily": "inherit",
    "fontSize": "0.875rem",
    "whiteSpace": "normal",
    "height": "auto",
    "minWidth": "90px",
    "maxWidth": "320px",
}

_BASE_HEADER_STYLE = {
    "backgroundColor": "#F8FAFC",
    "fontWeight": "700",
    "borderBottom": "2px solid #E2E8F0",
    "textTransform": "uppercase",
    "fontSize": "0.72rem",
    "letterSpacing": "0.05em",
    "whiteSpace": "normal",
    "height": "auto",
}

# Allows text selection inside DataTable cells so sequences can be copied.
_SELECTABLE_CSS = [
    {"selector": ".dash-cell-value", "rule": "user-select: text !important; -webkit-user-select: text !important;"},
    {"selector": "td, th", "rule": "user-select: text !important; -webkit-user-select: text !important;"},
    {"selector": ".dash-spreadsheet-container", "rule": "user-select: text !important;"},
]


def resistance_column_toggle(visible: Optional[List[str]] = None) -> html.Div:
    """Column visibility control for the resistance table (TASK-12.4)."""
    return html.Div([
        html.Label(
            [html.I(className="bi bi-layout-three-columns me-2"), "Columns"],
            className="small text-muted mb-1 d-block",
        ),
        dcc.Dropdown(
            id="resistance-column-toggle",
            options=[
                {"label": column["name"], "value": column["key"]}
                for column in OPTIONAL_RESISTANCE_COLUMNS
            ],
            value=list(visible if visible is not None else DEFAULT_VISIBLE_COLUMNS),
            multi=True,
            clearable=False,
            placeholder="Show optional WHO metadata columns…",
            className="column-toggle",
        ),
    ], className="mb-3")


def build_resistance_columns(visible_keys: Optional[List[str]] = None) -> List[Dict]:
    """DataTable column definitions for the selected optional columns."""
    visible = set(visible_keys if visible_keys is not None else DEFAULT_VISIBLE_COLUMNS)
    return [
        {"name": column["name"], "id": column["key"]}
        for column in RESISTANCE_COLUMNS
        if column.get("required") or column["key"] in visible
    ]


def prepare_resistance_data(drug_data: pd.DataFrame) -> pd.DataFrame:
    """
    Flatten catalogue rows into the table's own column ids.

    Working on stable ids rather than the catalogue's long header strings keeps
    the column-visibility callback and the styling rules readable.
    """
    prepared = pd.DataFrame(index=drug_data.index)
    for column in RESISTANCE_COLUMNS:
        source = column["source"]
        if source in drug_data.columns:
            prepared[column["key"]] = drug_data[source].fillna("")
        else:
            prepared[column["key"]] = ""

    # "Silent mutation" is a presence flag in the master file.
    prepared["silent_mutation"] = prepared["silent_mutation"].apply(
        lambda value: "Yes" if str(value).strip() else ""
    )
    prepared["variant"] = drug_data.get("variant", "")

    # Tier 1 first, then unnumbered tiers last.
    prepared["_tier_sort"] = pd.to_numeric(prepared["tier"], errors="coerce").fillna(99)
    prepared["_grade_sort"] = prepared["confidence"].apply(
        lambda value: next(
            (i for i, grade in enumerate(CONFIDENCE_GRADES) if grade == value), 99
        )
    )
    return (
        prepared.sort_values(["_tier_sort", "_grade_sort", "mutation"])
        .drop(columns=["_tier_sort", "_grade_sort"])
        .reset_index(drop=True)
    )


def resistance_style_conditional(selected_mutation: Optional[str] = None) -> List[Dict]:
    """Conditional styling for the resistance table."""
    styles: List[Dict] = [
        # Remarks and version changes are flagged in red (TASK-12.3).
        {
            "if": {"column_id": "comment", "filter_query": '{comment} ne ""'},
            "color": "#B42318",
            "fontStyle": "italic",
        },
        {
            "if": {
                "column_id": "changes_vs_ver1",
                "filter_query": '{changes_vs_ver1} ne "" && {changes_vs_ver1} ne "No change"',
            },
            "color": "#B42318",
            "fontWeight": "600",
        },
        {
            "if": {"column_id": "relaxed_thresholds", "filter_query": '{relaxed_thresholds} ne ""'},
            "color": "#B54708",
        },
        {
            "if": {"column_id": "confidence", "filter_query": '{confidence} contains "1) Assoc w R"'},
            "fontWeight": "700",
            "color": "#0B5D5D",
        },
        {"if": {"column_id": "mutation"}, "fontFamily": "ui-monospace, SFMono-Regular, Menlo, monospace"},
    ]

    if selected_mutation:
        # Highlight every row matching the active mutation (TASK-04, TASK-08).
        for candidate in _mutation_variants(selected_mutation):
            styles.append({
                "if": {"filter_query": f'{{mutation}} eq "{candidate}"'},
                "backgroundColor": HIGHLIGHT_COLOR,
                "borderTop": f"1px solid {HIGHLIGHT_BORDER}",
                "borderBottom": f"1px solid {HIGHLIGHT_BORDER}",
                "fontWeight": "600",
            })

    styles.append({
        "if": {"state": "active"},
        "backgroundColor": HIGHLIGHT_COLOR,
        "border": f"1px solid {HIGHLIGHT_BORDER}",
    })
    return styles


def _mutation_variants(mutation: str) -> List[str]:
    """Mutation spellings to match against, with and without a ``p.`` prefix."""
    bare = mutation
    for prefix in ("p.", "c.", "n."):
        if bare.lower().startswith(prefix):
            bare = bare[2:]
            break
    return list(dict.fromkeys([mutation, bare, f"p.{bare}", f"c.{bare}"]))


def create_drug_detail_table(
    drug_data: pd.DataFrame,
    drug: str,
    visible_columns: Optional[List[str]] = None,
    selected_mutation: Optional[str] = None,
    page_size: int = RESISTANCE_PAGE_SIZE,
) -> html.Div:
    """Create the detailed drug resistance table for one drug."""
    prepared = prepare_resistance_data(drug_data)

    table = dash_table.DataTable(
        id={"type": "drug-table", "index": drug},
        data=prepared.to_dict("records"),
        columns=build_resistance_columns(visible_columns),
        cell_selectable=True,
        sort_action="native",
        page_size=page_size,
        page_current=page_of_mutation(prepared, selected_mutation, page_size),
        style_table={"overflowX": "auto", "minWidth": "100%"},
        style_cell=_BASE_CELL_STYLE,
        style_header=_BASE_HEADER_STYLE,
        style_data={"borderBottom": "1px solid #F1F5F9"},
        style_data_conditional=resistance_style_conditional(selected_mutation),
        style_cell_conditional=[
            {"if": {"column_id": "tier"}, "width": "60px", "textAlign": "center"},
            {"if": {"column_id": "silent_mutation"}, "width": "90px", "textAlign": "center"},
            {"if": {"column_id": "comment"}, "minWidth": "200px"},
            {"if": {"column_id": "relaxed_thresholds"}, "minWidth": "180px"},
        ],
        tooltip_header={
            column["key"]: _COLUMN_TOOLTIPS[column["key"]]
            for column in RESISTANCE_COLUMNS
            if column["key"] in _COLUMN_TOOLTIPS
        },
        tooltip_delay=300,
        tooltip_duration=None,
        css=_SELECTABLE_CSS,
    )

    return html.Div(table, className="dash-table-container table-responsive-container")


def page_of_mutation(
    prepared: pd.DataFrame, mutation: Optional[str], page_size: int = RESISTANCE_PAGE_SIZE
) -> int:
    """Page index holding the first row matching a mutation (TASK-04.4)."""
    if not mutation or prepared.empty:
        return 0
    candidates = {value.lower() for value in _mutation_variants(mutation)}
    matches = prepared.index[prepared["mutation"].str.lower().isin(candidates)]
    if len(matches) == 0:
        return 0
    return int(matches[0]) // page_size


# ----------------------------------------------------------------------
# Genomic coordinates table
# ----------------------------------------------------------------------
def _truncate_sequence(value: str) -> str:
    """Shorten a long Ref/Alt allele for display, keeping its size visible."""
    text = "" if value is None else str(value)
    if len(text) <= SEQ_DISPLAY_CHARS:
        return text
    return f"{text[:SEQ_DISPLAY_CHARS]}… ({len(text):,} bp)"


def _sequence_tooltip(value: str) -> Optional[Dict]:
    """Markdown tooltip carrying the full (or leading) allele sequence."""
    text = "" if value is None else str(value)
    if len(text) <= SEQ_DISPLAY_CHARS:
        return None
    shown = text[:SEQ_TOOLTIP_CHARS]
    suffix = (
        f"\n\n_First {SEQ_TOOLTIP_CHARS:,} of {len(text):,} bp. "
        "Select the row to copy the full allele from Coordinate Analysis._"
        if len(text) > SEQ_TOOLTIP_CHARS
        else f"\n\n_{len(text):,} bp_"
    )
    return {"value": f"```\n{shown}\n```{suffix}", "type": "markdown"}


def prepare_coordinates_data(
    mutations_df: pd.DataFrame, gene_info: GeneInfo, coord_calculator
) -> pd.DataFrame:
    """Derive the display columns for the Genomic Coordinates table."""
    display = mutations_df.copy()

    numeric_position = pd.to_numeric(display["position"], errors="coerce")
    display["position_value"] = numeric_position

    display["relative_value"] = numeric_position.apply(
        lambda position: coord_calculator.calculate_relative_position(int(position), gene_info)
        if pd.notna(position) else None
    )
    display["Gene Relative"] = display["relative_value"].apply(
        lambda value: f"c.{int(value)}" if pd.notna(value) else "N/A"
    )
    display["Position Display"] = numeric_position.apply(
        lambda value: f"{int(value):,}" if pd.notna(value) else "N/A"
    )

    reference = display["reference_nucleotide"].fillna("")
    alternative = display["alternative_nucleotide"].fillna("")
    display["ref_display"] = reference.apply(_truncate_sequence)
    display["alt_display"] = alternative.apply(_truncate_sequence)
    display["Change"] = [
        _describe_change(ref, alt) for ref, alt in zip(reference, alternative)
    ]
    display["change_kind"] = [
        change_kind(ref, alt) for ref, alt in zip(reference, alternative)
    ]
    return display


def _describe_change(reference: str, alternative: str) -> str:
    """Compact human-readable description of an allele change."""
    ref_len, alt_len = len(reference), len(alternative)
    if ref_len == alt_len == 1:
        return "SNV"
    if ref_len > alt_len:
        return f"Deletion ({ref_len - alt_len:,} bp)"
    if alt_len > ref_len:
        return f"Insertion ({alt_len - ref_len:,} bp)"
    return f"MNV ({ref_len:,} bp)"


def coordinates_style_conditional(selected_variant: Optional[str] = None) -> List[Dict]:
    """Conditional styling for the coordinates table."""
    styles: List[Dict] = [
        {
            "if": {"column_id": c},
            "fontFamily": "ui-monospace, SFMono-Regular, Menlo, monospace",
        }
        for c in ("variant", "ref_display", "alt_display", "Gene Relative")
    ]
    if selected_variant:
        styles.append({
            "if": {"filter_query": f'{{variant}} eq "{selected_variant}"'},
            "backgroundColor": HIGHLIGHT_COLOR,
            "borderTop": f"1px solid {HIGHLIGHT_BORDER}",
            "borderBottom": f"1px solid {HIGHLIGHT_BORDER}",
            "fontWeight": "600",
        })
    styles.append({
        "if": {"state": "active"},
        "backgroundColor": HIGHLIGHT_COLOR,
        "border": f"1px solid {HIGHLIGHT_BORDER}",
    })
    return styles


COORDINATES_PAGE_SIZE = 10


def create_genomic_coords_table(
    mutations_df: pd.DataFrame,
    gene_info: GeneInfo,
    coord_calculator,
    selected_variant: Optional[str] = None,
) -> html.Div:
    """
    Create the genomic coordinates drill-down table.

    Long indel alleles are truncated to a fixed width with the full sequence
    available on hover, and the table scrolls horizontally rather than pushing
    later columns off-screen (TASK-07).
    """
    display = prepare_coordinates_data(mutations_df, gene_info, coord_calculator)

    columns = [
        {"name": "Variant", "id": "variant"},
        {"name": "Genomic Position", "id": "Position Display"},
        {"name": "Ref", "id": "ref_display"},
        {"name": "Alt", "id": "alt_display"},
        {"name": "Change", "id": "Change"},
        {"name": "Gene Relative", "id": "Gene Relative"},
    ]

    table = dash_table.DataTable(
        id=dict(COORDS_TABLE_ID),
        data=display.to_dict("records"),
        columns=columns,
        cell_selectable=True,
        sort_action="native",
        page_size=COORDINATES_PAGE_SIZE,
        page_current=page_of_variant(display, selected_variant),
        style_table={"overflowX": "auto", "minWidth": "100%"},
        style_cell=_BASE_CELL_STYLE,
        style_header=_BASE_HEADER_STYLE,
        style_data={"borderBottom": "1px solid #F1F5F9"},
        style_data_conditional=coordinates_style_conditional(selected_variant),
        style_cell_conditional=[
            {"if": {"column_id": "variant"}, "minWidth": "170px"},
            {"if": {"column_id": "ref_display"}, "minWidth": "150px", "maxWidth": "240px"},
            {"if": {"column_id": "alt_display"}, "minWidth": "150px", "maxWidth": "240px"},
            {"if": {"column_id": "Change"}, "minWidth": "130px"},
        ],
        tooltip_data=coordinates_tooltips(display),
        tooltip_delay=200,
        tooltip_duration=None,
        css=_SELECTABLE_CSS,
    )

    return html.Div(table, className="dash-table-container table-responsive-container")


def coordinates_tooltips(display: pd.DataFrame) -> List[Dict]:
    """Full-allele hover tooltips, one per row of ``display``."""
    tooltip_data = []
    for _, row in display.iterrows():
        tooltip = {}
        for column, source in (("ref_display", "reference_nucleotide"), ("alt_display", "alternative_nucleotide")):
            hint = _sequence_tooltip(row.get(source, ""))
            if hint:
                tooltip[column] = hint
        tooltip_data.append(tooltip)
    return tooltip_data


def page_of_variant(display: pd.DataFrame, variant: Optional[str]) -> int:
    """Page index holding the first row for a variant."""
    if not variant or display.empty:
        return 0
    matches = display.reset_index(drop=True).index[display["variant"].values == variant]
    if len(matches) == 0:
        return 0
    return int(matches[0]) // COORDINATES_PAGE_SIZE


# ----------------------------------------------------------------------
# Catalogue summary tables
# ----------------------------------------------------------------------
def _with_totals(summary: pd.DataFrame, label_column: str, label: str) -> pd.DataFrame:
    """Append an aggregate totals row to a summary table (TASK-14.2)."""
    if summary.empty:
        return summary

    totals = {}
    for column in summary.columns:
        if column == label_column:
            totals[column] = label
        elif pd.api.types.is_numeric_dtype(summary[column]):
            totals[column] = int(summary[column].sum())
        else:
            totals[column] = ""
    return pd.concat([summary, pd.DataFrame([totals])], ignore_index=True)


# WHO grading names are too long for a table header; the full text is offered
# as a header tooltip instead.
_SUMMARY_HEADER_LABELS = {
    "Associations": "Assoc.",
    "Drug count": "Drugs (n)",
    "1) Assoc w R": "1) AwR",
    "2) Assoc w R - Interim": "2) AwR-i",
    "3) Uncertain significance": "3) Uncertain",
    "4) Not assoc w R - Interim": "4) NotAwR-i",
    "5) Not assoc w R": "5) NotAwR",
}

_SUMMARY_HEADER_TOOLTIPS = {
    "Variants": "Distinct mutations catalogued for this row.",
    "Associations": "Catalogue rows, i.e. drug-variant associations.",
    "Drug count": "Number of drugs this gene is catalogued for.",
    "Tier 1": "Associations in a WHO tier 1 (established resistance) gene.",
    "Tier 2": "Associations in a WHO tier 2 (candidate) gene.",
    "LoF": "Distinct loss-of-function mutations.",
    "1) Assoc w R": "Final WHO grading: 1) Assoc w R",
    "2) Assoc w R - Interim": "Final WHO grading: 2) Assoc w R - Interim",
    "3) Uncertain significance": "Final WHO grading: 3) Uncertain significance",
    "4) Not assoc w R - Interim": "Final WHO grading: 4) Not assoc w R - Interim",
    "5) Not assoc w R": "Final WHO grading: 5) Not assoc w R",
}


def _summary_table(table_id: str, summary: pd.DataFrame, label_column: str) -> dash_table.DataTable:
    text_columns = (label_column, "Drugs", "Locus tag")
    columns = [
        {
            "name": _SUMMARY_HEADER_LABELS.get(column, column),
            "id": column,
            "type": "text" if column in text_columns else "numeric",
        }
        for column in summary.columns
    ]
    return dash_table.DataTable(
        id=table_id,
        data=summary.to_dict("records"),
        columns=columns,
        sort_action="native",
        filter_action="native",
        filter_options={"case": "insensitive", "placeholder_text": "Filter…"},
        page_action="none",
        fixed_rows={"headers": True},
        style_table={"overflowX": "auto", "maxHeight": "58vh", "overflowY": "auto"},
        style_cell=dict(_BASE_CELL_STYLE, fontSize="0.82rem"),
        style_header=_BASE_HEADER_STYLE,
        style_data={"borderBottom": "1px solid #F1F5F9"},
        style_cell_conditional=[
            {"if": {"column_id": "Drugs"}, "minWidth": "220px", "maxWidth": "320px"},
        ],
        style_data_conditional=[
            {
                "if": {"filter_query": f'{{{label_column}}} eq "TOTAL"'},
                "fontWeight": "700",
                "backgroundColor": "#F1F7F7",
                "borderTop": f"2px solid {HIGHLIGHT_BORDER}",
            },
            {"if": {"column_id": label_column}, "fontWeight": "600"},
        ],
        tooltip_header={
            column: tooltip
            for column, tooltip in _SUMMARY_HEADER_TOOLTIPS.items()
            if column in summary.columns
        },
        tooltip_delay=300,
        tooltip_duration=None,
        css=_SELECTABLE_CSS,
    )


def catalogue_summary_view(data_loader) -> html.Div:
    """Summary statistics view over the whole WHO catalogue (TASK-14)."""
    totals = data_loader.get_catalogue_totals()
    gene_summary = _with_totals(data_loader.get_catalogue_summary(), "Gene", "TOTAL")
    drug_summary = _with_totals(data_loader.get_drug_summary(), "Drug", "TOTAL")

    stat_cards = dbc.Row([
        _stat_card("Genes", totals["genes"], "bi-dna"),
        _stat_card("Drugs", totals["drugs"], "bi-capsule"),
        _stat_card("Distinct variants", totals["variants"], "bi-list-columns"),
        _stat_card("Drug–variant associations", totals["associations"], "bi-diagram-2"),
    ], className="g-3 mb-4")

    return html.Div([
        stat_cards,
        html.P([
            "A catalogue row is one drug–variant association, so a variant graded for "
            "several drugs is counted once per drug. ",
            html.Strong("Variants"),
            " therefore counts distinct mutations while ",
            html.Strong("Associations"),
            " counts catalogue rows. Tier counts follow the WHO gene tier "
            "(1 = established resistance gene, 2 = candidate), and the numbered "
            "columns are final WHO confidence gradings.",
        ], className="text-muted small"),
        dbc.Tabs([
            dbc.Tab(
                html.Div(
                    _summary_table("summary-gene-table", gene_summary, "Gene"),
                    className="dash-table-container table-responsive-container pt-3",
                ),
                label="By gene",
                tab_id="summary-by-gene",
            ),
            dbc.Tab(
                html.Div(
                    _summary_table("summary-drug-table", drug_summary, "Drug"),
                    className="dash-table-container table-responsive-container pt-3",
                ),
                label="By drug",
                tab_id="summary-by-drug",
            ),
        ], active_tab="summary-by-gene"),
    ])


def _stat_card(label: str, value: int, icon: str) -> dbc.Col:
    return dbc.Col(
        html.Div([
            html.I(className=f"bi {icon} stat-icon"),
            html.Div([
                html.Div(f"{value:,}", className="stat-value"),
                html.Div(label, className="stat-label"),
            ]),
        ], className="stat-card"),
        xs=6, lg=3,
    )
