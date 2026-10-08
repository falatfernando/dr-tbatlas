"""
DR-TBAtlas - A Dash application for visualizing TB genomic data and drug resistance.
"""
import os
from functools import lru_cache
from typing import Dict, List, Optional

import dash
import dash_bootstrap_components as dbc
import pandas as pd
from dash import ALL, Input, Output, State, callback_context, dcc, html
from dash.exceptions import PreventUpdate
from flask import send_from_directory

import browser_tracks
import genome_view
import layout as ui
import mutation_filters as mf
import search_utils
import tables
from coordinate_calculator import CoordinateCalculator
from data_utils import DataLoader, GeneInfo

# Get absolute path for data directory
DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')

# Initialize data loader
data_loader = DataLoader(data_dir=DATA_DIR)
data_loader.load_gff3()
data_loader.load_catalogue()
data_loader.load_genomic_coordinates()
# Flattened annotation and the derived catalogue, GC and search tracks used
# by the genome browser (TASK-09). The desktop launcher moves them to a
# per-user cache, since a packaged app's own folder may be read-only.
TRACKS_DIR = os.environ.get("DRTBATLAS_TRACKS_DIR") or os.path.join(
    DATA_DIR, browser_tracks.TRACKS_DIRNAME
)
browser = browser_tracks.BrowserTracks(data_loader, tracks_dir=TRACKS_DIR).build()

coord_calculator = CoordinateCalculator(data_loader)

DRUG_GENE_MAP = data_loader.get_drug_gene_map()
CATALOGUE_TOTALS = data_loader.get_catalogue_totals()

# Like the coordinates table, the analysis panel is created only once a gene
# is loaded, so it carries a pattern-matching id.
ANALYSIS_DISPLAY_ID = {"type": "analysis-display", "index": "main"}

# Initialize Dash app
app = dash.Dash(
    __name__,
    external_stylesheets=[dbc.themes.BOOTSTRAP, dbc.icons.BOOTSTRAP],
    suppress_callback_exceptions=True,
    meta_tags=[{"name": "viewport", "content": "width=device-width, initial-scale=1"}]
)
app.title = f"{ui.APP_NAME} | LaPAM"

# Expose the Flask server for gunicorn
server = app.server


# Flask route to serve data files for JBrowse
@server.route('/data/<path:path>')
def serve_data(path):
    tracks_prefix = browser_tracks.TRACKS_DIRNAME + "/"
    if path.startswith(tracks_prefix):
        return send_from_directory(TRACKS_DIR, path[len(tracks_prefix):])
    return send_from_directory(DATA_DIR, path)


# App layout
app.layout = html.Div([
    # Active gene and mutation requested by the user, from any entry point.
    dcc.Store(id="nav-store"),
    dcc.Store(id="current-gene-store"),
    # Single source of truth for the selected variant, shared by the
    # resistance table, the coordinates table and the analysis card.
    dcc.Store(id="selection-store"),

    ui.navbar(app),
    ui.hero_search(app),

    dbc.Container([
        dcc.Loading(
            id="loading-results",
            type="default",
            children=html.Div(id="search-results-container"),
        ),
        ui.footer(app),
    ], fluid=False),

    ui.browse_by_drug_modal(DRUG_GENE_MAP),
    ui.summary_modal(),
    ui.cite_modal(),
])


# ----------------------------------------------------------------------
# Chrome: navbar collapse and modals
# ----------------------------------------------------------------------
@app.callback(
    Output("navbar-collapse", "is_open"),
    Input("navbar-toggler", "n_clicks"),
    State("navbar-collapse", "is_open"),
    prevent_initial_call=True,
)
def toggle_navbar(n_clicks, is_open):
    return not is_open


@app.callback(
    Output("drugs-modal", "is_open"),
    [Input("open-drugs-modal", "n_clicks"),
     Input("close-drugs-modal", "n_clicks"),
     Input({"type": "gene-link", "drug": ALL, "index": ALL}, "n_clicks")],
    State("drugs-modal", "is_open"),
    prevent_initial_call=True,
)
def toggle_drugs_modal(open_clicks, close_clicks, gene_clicks, is_open):
    """Open from the header; close on explicit close or after picking a gene."""
    triggered = callback_context.triggered_id
    if triggered == "open-drugs-modal":
        return True
    return False


@app.callback(
    [Output("summary-modal", "is_open"),
     Output("summary-modal-body", "children")],
    [Input("open-summary-modal", "n_clicks"),
     Input("close-summary-modal", "n_clicks")],
    State("summary-modal", "is_open"),
    prevent_initial_call=True,
)
def toggle_summary_modal(open_clicks, close_clicks, is_open):
    """Build the summary tables on first open rather than at start-up."""
    if callback_context.triggered_id == "open-summary-modal":
        return True, tables.catalogue_summary_view(data_loader)
    return False, dash.no_update


@app.callback(
    Output("cite-modal", "is_open"),
    [Input("open-cite-modal", "n_clicks"),
     Input("open-cite-modal-footer", "n_clicks"),
     Input("close-cite-modal", "n_clicks")],
    prevent_initial_call=True,
)
def toggle_cite_modal(header_clicks, footer_clicks, close_clicks):
    return callback_context.triggered_id != "close-cite-modal"


# ----------------------------------------------------------------------
# Routing: every entry point resolves to one active gene
# ----------------------------------------------------------------------
@app.callback(
    [Output("nav-store", "data"),
     Output("search-feedback", "children"),
     Output("search-input", "value")],
    [Input("search-form", "n_submit"),
     Input({"type": "quick-search", "index": ALL}, "n_clicks"),
     Input({"type": "gene-link", "drug": ALL, "index": ALL}, "n_clicks"),
     Input({"type": "neighbor-link", "index": ALL}, "n_clicks")],
    State("search-input", "value"),
    prevent_initial_call=True,
)
def route_query(
    form_submit,
    quick_clicks,
    gene_link_clicks,
    neighbor_clicks,
    search_value,
):
    """
    Resolve a query from any entry point into the active gene.

    Handles the search form (submitted with the button or the Enter key,
    TASK-01), the quick-access pills, the Browse by Drug chips and the
    neighbour shortcuts (TASK-10).
    """
    triggered = callback_context.triggered_id
    if triggered is None:
        raise PreventUpdate

    # Clicks elsewhere in the interface carry their target in the component id:
    # the Browse by Drug chips and neighbour shortcuts carry a locus tag, while
    # a quick-access pill carries a query, which may be a full variant.
    if isinstance(triggered, dict):
        if not _triggered_click_value():
            raise PreventUpdate

        token = triggered.get("index")
        if triggered.get("type") == "quick-search":
            parsed = search_utils.parse_query(
                token, gene_resolver=data_loader.is_known_gene
            )
            return _navigate(parsed.gene_token, parsed.mutation, query=token)
        return _navigate(token)

    # Search submission.
    parsed = search_utils.parse_query(search_value, gene_resolver=data_loader.is_known_gene)

    if parsed.kind == "empty":
        return None, None, dash.no_update

    if parsed.kind == "standalone_mutation":
        # A bare mutation cannot be resolved uniquely (TASK-03).
        return dash.no_update, ui.mutation_advisory(parsed.message), dash.no_update

    return _navigate(parsed.gene_token, parsed.mutation, query=parsed.raw)


def _triggered_click_value():
    """The n_clicks value that triggered a pattern-matching callback."""
    triggered = callback_context.triggered
    return triggered[0]["value"] if triggered else None


def _navigate(token: Optional[str], mutation: Optional[str] = None, query: Optional[str] = None):
    """Build the routing outputs for a gene token, resolving it first."""
    gene_info = data_loader.get_gene_info(token)

    if gene_info is None:
        # Fall back to a partial match before giving up.
        matches = data_loader.search_genes(token) if token else []
        if matches:
            gene_info = matches[0]
        else:
            return dash.no_update, ui.gene_not_found_alert(query or token or ""), dash.no_update

    nav = {"locus_tag": gene_info.locus_tag, "mutation": mutation}
    display_query = (
        f"{gene_info.display_name}_{mutation}" if mutation else gene_info.display_name
    )
    return nav, None, display_query


# ----------------------------------------------------------------------
# Rendering the active gene
# ----------------------------------------------------------------------
@app.callback(
    [Output("search-results-container", "children"),
     Output("current-gene-store", "data"),
     Output("selection-store", "data")],
    Input("nav-store", "data"),
)
def render_gene(nav: Optional[Dict]):
    """Render every panel for the active gene."""
    if not nav or not nav.get("locus_tag"):
        return ui.welcome_panel(CATALOGUE_TOTALS), None, None

    gene_info = data_loader.get_gene_info(nav["locus_tag"])
    if gene_info is None:
        return ui.gene_not_found_alert(nav["locus_tag"]), None, None

    requested_mutation = nav.get("mutation")
    in_catalogue = data_loader.is_gene_in_catalogue(gene_info)
    catalogue_gene = data_loader.resolve_catalogue_gene(gene_info.locus_tag) or \
        data_loader.resolve_catalogue_gene(gene_info.gene_name)

    selected_variant = (
        f"{catalogue_gene}_{requested_mutation}"
        if catalogue_gene and requested_mutation else None
    )

    gene_data = {
        "gene_name": gene_info.gene_name,
        "locus_tag": gene_info.locus_tag,
        "start": gene_info.start,
        "end": gene_info.end,
        "strand": gene_info.strand,
        "chromosome": gene_info.chromosome,
        "product": gene_info.product,
        "protein_length": gene_info.protein_length,
        "catalogue_gene": catalogue_gene,
        "in_catalogue": in_catalogue,
    }

    results: List = [_gene_overview_card(gene_info, in_catalogue)]

    # Embedded JBrowse view: the genome explorer.
    results.append(genome_view.jbrowse_card(
        genome_view.jbrowse_view(browser, gene_info), browser.drugs
    ))
    results.append(genome_view.sequence_card(gene_info))

    if not in_catalogue:
        # Annotated but not catalogued: no resistance or coordinate tables
        # would have anything to show (TASK-05).
        results.append(ui.non_catalogue_card(gene_info))
        return results, gene_data, None

    drug_resistance = data_loader.get_drug_resistance_info(gene_info.locus_tag)
    if len(drug_resistance) > 0:
        results.append(_resistance_card(
            drug_resistance, gene_info.locus_tag, requested_mutation
        ))

    mutations_df = data_loader.search_mutations_by_gene(gene_info.locus_tag)
    if len(mutations_df) > 0:
        results.append(_coordinates_card(
            mutations_df, gene_info, selected_variant,
            linkable=len(drug_resistance) > 0,
        ))

    results.append(_analysis_card())

    initial_selection = None
    if requested_mutation:
        # Deep link straight to the requested variant (TASK-04).
        initial_selection = _build_selection(
            variant=selected_variant,
            mutation=requested_mutation,
            gene_data=gene_data,
        )

    return results, gene_data, initial_selection


def _gene_overview_card(gene_info: GeneInfo, in_catalogue: bool) -> dbc.Card:
    """Gene overview header card, including protein length (TASK-06)."""
    previous_gene, next_gene = data_loader.get_neighbor_genes(gene_info)

    if gene_info.protein_length:
        length_display = html.Span([
            f"{gene_info.length:,} bp",
            html.Span(f" ({gene_info.protein_length:,} aa)", className="text-muted"),
        ], className="fw-500")
    else:
        length_display = html.Span([
            f"{gene_info.length:,} bp",
            html.Span(
                f" ({gene_info.biotype.replace('_', ' ') or 'non-coding'})",
                className="text-muted",
            ),
        ], className="fw-500")

    return dbc.Card([
        dbc.CardHeader([
            dbc.Row([
                dbc.Col([
                    html.H4(gene_info.display_name, className="mb-0 fw-bold text-primary"),
                    html.Small(f"Locus: {gene_info.locus_tag}", className="text-muted"),
                    dbc.Badge(
                        "In WHO catalogue" if in_catalogue else "Not in WHO catalogue",
                        className=("badge-tb ms-2" if in_catalogue else "badge-muted ms-2"),
                    ),
                ], md=7),
                dbc.Col([
                    html.Div([
                        html.Span("Product: ", className="fw-bold"),
                        html.Span(gene_info.product or "Unknown"),
                    ], className="text-md-end text-muted small"),
                ], md=5),
            ], align="center"),
        ], className="card-header-custom"),
        dbc.CardBody([
            dbc.Row([
                _overview_field(
                    "Genomic Region",
                    f"{gene_info.chromosome}:{gene_info.start:,}-{gene_info.end:,}",
                ),
                dbc.Col([
                    html.Div([
                        html.Label("Strand", className="text-muted small d-block"),
                        dbc.Badge(gene_info.strand, className="badge-custom"),
                    ])
                ], md=3, xs=6),
                dbc.Col([
                    html.Div([
                        html.Label("Length", className="text-muted small d-block"),
                        length_display,
                    ])
                ], md=3, xs=6),
                _overview_field(
                    "Biotype",
                    (gene_info.biotype or "unknown").replace("_", " "),
                    width=3,
                ),
            ], className="gy-3"),

            html.Hr(className="my-3"),
            genome_view.neighbor_navigation(previous_gene, next_gene),

            html.Div([
                html.Label("Functional note", className="text-muted small d-block mt-3"),
                html.P(gene_info.note, className="small mb-0 gene-note"),
            ]) if gene_info.note else None,
        ]),
    ], className="result-card")


def _overview_field(label: str, value: str, width: int = 3) -> dbc.Col:
    return dbc.Col([
        html.Div([
            html.Label(label, className="text-muted small d-block"),
            html.Span(value, className="fw-500"),
        ])
    ], md=width, xs=6)


@lru_cache(maxsize=8)
def _prepared_resistance(locus_tag: str) -> Dict[str, pd.DataFrame]:
    """
    Prepared resistance rows of a gene, per drug, in table order.

    Cached because every change to the mutation filters re-filters them.
    Callers must not modify the frames.
    """
    drug_resistance = data_loader.get_drug_resistance_info(locus_tag)
    return {
        drug: tables.prepare_resistance_data(drug_resistance[drug_resistance["drug"] == drug])
        for drug in sorted(drug_resistance["drug"].dropna().unique())
    }


@lru_cache(maxsize=8)
def _prepared_coordinates(locus_tag: str) -> pd.DataFrame:
    """Prepared Genomic Coordinates rows of a gene; callers must not modify it."""
    gene_info = data_loader.get_gene_info(locus_tag)
    mutations_df = data_loader.search_mutations_by_gene(locus_tag)
    return tables.prepare_coordinates_data(mutations_df, gene_info, coord_calculator)


def _drug_count_label(shown: int, total: int) -> str:
    if shown == total:
        return f"{total:,} mutations"
    return f"{shown:,} of {total:,} mutations"


def _resistance_card(
    drug_resistance: pd.DataFrame, locus_tag: str, selected_mutation: Optional[str]
) -> dbc.Card:
    """Drug Resistance Profile card, one table per drug (TASK-12)."""
    drug_cards = []
    for drug in sorted(drug_resistance["drug"].dropna().unique()):
        drug_data = drug_resistance[drug_resistance["drug"] == drug]
        tiers = sorted({str(t) for t in drug_data["tier"].dropna().unique()})
        graded = drug_data["FINAL CONFIDENCE GRADING"].fillna("")
        associated = int(graded.str.startswith(("1)", "2)")).sum())

        drug_cards.append(
            dbc.Card([
                dbc.CardHeader([
                    dbc.Row([
                        dbc.Col(html.H6(drug, className="mb-0 fw-bold"), width="auto"),
                        dbc.Col(dbc.Badge(
                            _drug_count_label(len(drug_data), len(drug_data)),
                            id={"type": "drug-count", "index": drug},
                            className="badge-custom ms-2",
                        ), width="auto"),
                        dbc.Col(dbc.Badge(
                            f"{associated:,} assoc. with resistance",
                            className="badge-assoc ms-1",
                        ), width="auto"),
                        dbc.Col(html.Small(
                            f"Tier {', '.join(tiers) if tiers else 'n/a'}",
                            className="text-muted",
                        ), width="auto"),
                    ], align="center", className="g-2"),
                ], className="bg-light border-0 py-2"),
                dbc.CardBody([
                    tables.create_drug_detail_table(
                        drug_data, drug, selected_mutation=selected_mutation
                    ),
                ], className="pt-2"),
            ], id={"type": "drug-card", "index": drug}, className="mb-3 border-0 shadow-sm")
        )

    return dbc.Card([
        dbc.CardHeader([
            html.H5(
                [html.I(className="bi bi-capsule me-2"), "Drug Resistance Profile"],
                className="card-header-title",
            )
        ], className="card-header-custom"),
        dbc.CardBody([
            html.P(
                "Click any mutation to link it to its genomic coordinates and "
                "coordinate analysis below.",
                className="text-muted small mb-2",
            ),
            mf.resistance_filter_panel(_prepared_resistance(locus_tag)),
            tables.resistance_column_toggle(),
            html.Div(
                mf.no_matches_alert(),
                id="resistance-filter-empty",
                style={"display": "none"},
            ),
            html.Div(drug_cards),
        ]),
    ], className="result-card")


def _coordinates_card(
    mutations_df: pd.DataFrame,
    gene_info: GeneInfo,
    selected_variant: Optional[str],
    linkable: bool = False,
) -> dbc.Card:
    """
    Genomic Coordinates card (TASK-07).

    ``linkable`` offers to restrict the table to the variants left by the Drug
    Resistance Profile filters, when that card is shown.
    """
    return dbc.Card([
        dbc.CardHeader([
            dbc.Row([
                dbc.Col(html.H5(
                    [html.I(className="bi bi-map me-2"), "Genomic Coordinates"],
                    className="card-header-title",
                )),
                dbc.Col(html.Small(
                    f"{len(mutations_df):,} nucleotide changes",
                    className="text-muted text-end",
                ), width="auto"),
            ], align="center"),
        ], className="card-header-custom"),
        dbc.CardBody([
            html.P([
                "Select a row to view detailed coordinate calculations. Long "
                "insertion and deletion alleles are shortened — hover a cell to "
                "read the sequence, or select the row to copy it in full.",
            ], className="text-muted small mb-3"),
            mf.coordinates_filter_bar(_prepared_coordinates(gene_info.locus_tag), linkable),
            tables.create_genomic_coords_table(
                mutations_df, gene_info, coord_calculator, selected_variant
            ),
        ]),
    ], className="result-card")


def _analysis_card() -> dbc.Card:
    return dbc.Card([
        dbc.CardHeader([
            html.H5(
                [html.I(className="bi bi-calculator me-2"), "Coordinate Analysis"],
                className="card-header-title",
            )
        ], className="card-header-custom"),
        dbc.CardBody([
            html.Div(id=dict(ANALYSIS_DISPLAY_ID), children=_analysis_placeholder()),
        ]),
    ], id="calculator-card", className="result-card mb-5")


def _analysis_placeholder() -> html.Div:
    return html.Div([
        html.I(className="bi bi-info-circle me-2"),
        "Select a mutation in the Drug Resistance Profile or a row in the "
        "Genomic Coordinates table to see its coordinate analysis.",
    ], className="text-muted text-center py-4")


# ----------------------------------------------------------------------
# Shared selection state
# ----------------------------------------------------------------------
def _build_selection(
    variant: Optional[str],
    mutation: Optional[str],
    gene_data: Dict,
    position: Optional[int] = None,
    reference: Optional[str] = None,
    alternative: Optional[str] = None,
) -> Dict:
    """
    Normalise a selection, filling in coordinates when they are not known.

    A click in the resistance table identifies a mutation but no specific
    nucleotide change, so the first recorded coordinate for that variant is
    resolved here to drive the analysis card (TASK-08).
    """
    selection = {
        "variant": variant,
        "mutation": mutation,
        "position": position,
        "reference": reference,
        "alternative": alternative,
        "coordinate_count": None,
    }

    if variant:
        coordinates = data_loader.get_variant_coordinates(variant)
        selection["coordinate_count"] = int(len(coordinates))
        if position is None and len(coordinates) > 0:
            first = coordinates.iloc[0]
            selection["position"] = (
                int(first["position"]) if str(first["position"]).isdigit() else None
            )
            selection["reference"] = first.get("reference_nucleotide")
            selection["alternative"] = first.get("alternative_nucleotide")
            selection["resolved_from_catalogue"] = True

    return selection


@app.callback(
    Output("selection-store", "data", allow_duplicate=True),
    Input(dict(tables.COORDS_TABLE_ID, index=ALL), "active_cell"),
    [State(dict(tables.COORDS_TABLE_ID, index=ALL), "derived_viewport_data"),
     State("current-gene-store", "data")],
    prevent_initial_call=True,
)
def select_from_coordinates(active_cells, viewport_datas, gene_data):
    """A click in the Genomic Coordinates table selects that exact change."""
    active_cell = active_cells[0] if active_cells else None
    viewport_data = viewport_datas[0] if viewport_datas else None
    if not active_cell or not viewport_data or not gene_data:
        raise PreventUpdate

    row_index = active_cell["row"]
    if row_index >= len(viewport_data):
        raise PreventUpdate

    row = viewport_data[row_index]
    variant = row.get("variant")
    position = row.get("position")
    mutation = variant.split("_", 1)[1] if variant and "_" in variant else None

    return _build_selection(
        variant=variant,
        mutation=mutation,
        gene_data=gene_data,
        position=int(position) if str(position).isdigit() else None,
        reference=row.get("reference_nucleotide"),
        alternative=row.get("alternative_nucleotide"),
    )


@app.callback(
    Output("selection-store", "data", allow_duplicate=True),
    Input({"type": "drug-table", "index": ALL}, "active_cell"),
    [State({"type": "drug-table", "index": ALL}, "derived_viewport_data"),
     State("current-gene-store", "data")],
    prevent_initial_call=True,
)
def select_from_resistance(active_cells, viewports, gene_data):
    """A click in the Drug Resistance Profile selects that mutation (TASK-08)."""
    triggered = callback_context.triggered_id
    if not triggered or not gene_data:
        raise PreventUpdate

    ids = [item["id"] for item in callback_context.inputs_list[0]]
    try:
        position = next(
            i for i, item in enumerate(ids) if dict(item) == dict(triggered)
        )
    except StopIteration:
        raise PreventUpdate

    active_cell = active_cells[position]
    viewport = viewports[position] if position < len(viewports) else None
    if not active_cell or not viewport or active_cell["row"] >= len(viewport):
        raise PreventUpdate

    row = viewport[active_cell["row"]]
    mutation = row.get("mutation")
    variant = row.get("variant") or (
        f"{gene_data.get('catalogue_gene')}_{mutation}"
        if gene_data.get("catalogue_gene") and mutation else None
    )

    return _build_selection(variant=variant, mutation=mutation, gene_data=gene_data)


# ----------------------------------------------------------------------
# Selection consumers
# ----------------------------------------------------------------------
@app.callback(
    [Output(dict(tables.COORDS_TABLE_ID, index=ALL), "style_data_conditional"),
     Output(dict(tables.COORDS_TABLE_ID, index=ALL), "page_current")],
    Input("selection-store", "data"),
    State(dict(tables.COORDS_TABLE_ID, index=ALL), "data"),
    prevent_initial_call=True,
)
def sync_coordinates_table(selection, table_datas):
    """Highlight and page to the selected variant in the coordinates table."""
    variant = (selection or {}).get("variant")
    styles = tables.coordinates_style_conditional(variant)

    pages = []
    for rows in table_datas:
        page = 0
        if variant and rows:
            page = next(
                (i // tables.COORDINATES_PAGE_SIZE for i, row in enumerate(rows)
                 if row.get("variant") == variant),
                0,
            )
        pages.append(page)

    return [styles] * len(table_datas), pages


@app.callback(
    [Output({"type": "drug-table", "index": ALL}, "style_data_conditional"),
     Output({"type": "drug-table", "index": ALL}, "page_current")],
    Input("selection-store", "data"),
    State({"type": "drug-table", "index": ALL}, "data"),
    prevent_initial_call=True,
)
def sync_resistance_tables(selection, table_data):
    """Highlight the selected mutation in every drug table (TASK-04, TASK-08)."""
    mutation = (selection or {}).get("mutation")
    styles = tables.resistance_style_conditional(mutation)
    count = len(table_data)

    pages = [
        tables.page_of_mutation(pd.DataFrame(rows or []), mutation)
        if mutation and rows else 0
        for rows in table_data
    ]
    return [styles] * count, pages


@app.callback(
    Output({"type": "drug-table", "index": ALL}, "columns"),
    Input("resistance-column-toggle", "value"),
    prevent_initial_call=True,
)
def toggle_resistance_columns(visible_keys):
    """Show or hide the optional WHO metadata columns (TASK-12.4)."""
    columns = tables.build_resistance_columns(visible_keys)
    return [columns] * len(callback_context.outputs_list)


# ----------------------------------------------------------------------
# Mutation filters
# ----------------------------------------------------------------------
@app.callback(
    [Output("filter-advanced-collapse", "is_open"),
     Output("filter-advanced-chevron", "className")],
    Input("filter-advanced-toggle", "n_clicks"),
    State("filter-advanced-collapse", "is_open"),
    prevent_initial_call=True,
)
def toggle_advanced_filters(n_clicks, is_open):
    is_open = not is_open
    return is_open, f"bi bi-chevron-{'up' if is_open else 'down'} ms-2"


_GRADE_PRESETS = {
    "filter-preset-associated": mf.ASSOCIATED_GRADES,
    "filter-preset-uncertain": mf.UNCERTAIN_GRADES,
    "filter-preset-not-associated": mf.NOT_ASSOCIATED_GRADES,
}


@app.callback(
    Output({"type": mf.RESISTANCE_FILTER, "field": ALL}, "value"),
    [Input("filter-clear", "n_clicks"),
     *[Input(button, "n_clicks") for button in _GRADE_PRESETS]],
    State({"type": mf.RESISTANCE_FILTER, "field": ALL}, "value"),
    prevent_initial_call=True,
)
def reset_or_preset_filters(clear_clicks, *args):
    """Clear every filter, or set the confidence grading to a quick pick."""
    current = args[-1]
    triggered = callback_context.triggered_id
    if not _triggered_click_value():
        raise PreventUpdate

    fields = [item["id"]["field"] for item in callback_context.outputs_list]
    if triggered == "filter-clear":
        return [mf.RESISTANCE_DEFAULTS[field] for field in fields]

    grades = _GRADE_PRESETS[triggered]
    return [
        list(grades) if field == "grades" else value
        for field, value in zip(fields, current)
    ]


@app.callback(
    [Output({"type": "drug-table", "index": ALL}, "data"),
     Output({"type": "drug-table", "index": ALL}, "page_current", allow_duplicate=True),
     Output({"type": "drug-table", "index": ALL}, "active_cell"),
     Output({"type": "drug-table", "index": ALL}, "selected_cells"),
     Output({"type": "drug-card", "index": ALL}, "style"),
     Output({"type": "drug-count", "index": ALL}, "children"),
     Output("resistance-filter-summary", "children"),
     Output("resistance-filter-empty", "style"),
     Output("filter-clear", "disabled"),
     Output("filter-advanced-count", "children"),
     Output("filter-advanced-count", "style")],
    Input({"type": mf.RESISTANCE_FILTER, "field": ALL}, "value"),
    [State("current-gene-store", "data"),
     State("selection-store", "data")],
    prevent_initial_call=True,
)
def apply_resistance_filters(values, gene_data, selection):
    """Filter every drug table of the Drug Resistance Profile at once."""
    if not gene_data:
        raise PreventUpdate

    filters = mf.collect(callback_context.inputs_list[0], mf.RESISTANCE_DEFAULTS)
    prepared = _prepared_resistance(gene_data["locus_tag"])
    filtered = mf.filter_by_drug(prepared, filters)
    mutation = (selection or {}).get("mutation")

    table_ids = [item["id"]["index"] for item in callback_context.outputs_list[0]]
    card_ids = [item["id"]["index"] for item in callback_context.outputs_list[4]]
    count_ids = [item["id"]["index"] for item in callback_context.outputs_list[5]]
    empty = pd.DataFrame(columns=["mutation"])

    datas, pages = [], []
    for drug in table_ids:
        rows = filtered.get(drug, empty)
        datas.append(rows.to_dict("records"))
        pages.append(tables.page_of_mutation(rows, mutation))

    shown = sum(len(rows) for rows in filtered.values())
    total = sum(len(rows) for rows in prepared.values())
    drugs_shown = sum(1 for rows in filtered.values() if len(rows))
    active = mf.is_filtered(filters)
    advanced = mf.advanced_count(filters)

    return (
        datas,
        pages,
        [None] * len(table_ids),
        [[]] * len(table_ids),
        [{"display": "block" if len(filtered.get(drug, empty)) else "none"} for drug in card_ids],
        [_drug_count_label(len(filtered.get(drug, empty)), len(prepared.get(drug, empty)))
         for drug in count_ids],
        mf.resistance_summary(shown, total, drugs_shown, len(prepared), active),
        {"display": "none" if shown else "block"},
        not active,
        str(advanced),
        {"display": "inline-block" if advanced else "none"},
    )


@app.callback(
    [Output(dict(tables.COORDS_TABLE_ID, index=ALL), "data"),
     Output(dict(tables.COORDS_TABLE_ID, index=ALL), "tooltip_data"),
     Output(dict(tables.COORDS_TABLE_ID, index=ALL), "page_current", allow_duplicate=True),
     Output(dict(tables.COORDS_TABLE_ID, index=ALL), "active_cell"),
     Output(dict(tables.COORDS_TABLE_ID, index=ALL), "selected_cells"),
     Output(dict(mf.COORDS_FILTER_SUMMARY_ID, index=ALL), "children")],
    [Input({"type": mf.COORDS_FILTER, "field": ALL}, "value"),
     Input({"type": mf.RESISTANCE_FILTER, "field": ALL}, "value")],
    [State("current-gene-store", "data"),
     State("selection-store", "data")],
    prevent_initial_call=True,
)
def apply_coordinate_filters(coords_values, resistance_values, gene_data, selection):
    """
    Filter the Genomic Coordinates table, optionally to the variants left by
    the Drug Resistance Profile filters.
    """
    tables_count = len(callback_context.outputs_list[0])
    if not gene_data or not tables_count:
        raise PreventUpdate

    filters = mf.collect(callback_context.inputs_list[0], mf.COORDS_DEFAULTS)
    linked = "linked" in (filters.get("linked") or [])
    triggered = callback_context.triggered_id
    if isinstance(triggered, dict) and triggered.get("type") == mf.RESISTANCE_FILTER and not linked:
        raise PreventUpdate

    variants = None
    if linked:
        resistance_filters = mf.collect(callback_context.inputs_list[1], mf.RESISTANCE_DEFAULTS)
        filtered = mf.filter_by_drug(_prepared_resistance(gene_data["locus_tag"]), resistance_filters)
        variants = {
            variant for rows in filtered.values() for variant in rows.get("variant", [])
        }

    display = _prepared_coordinates(gene_data["locus_tag"])
    rows = mf.filter_coordinates(display, filters, variants)
    page = tables.page_of_variant(rows, (selection or {}).get("variant"))
    summary = mf.coordinates_summary(
        len(rows), len(display), linked or mf.is_filtered(filters, mf.COORDS_DEFAULTS)
    )

    return (
        [rows.to_dict("records")] * tables_count,
        [tables.coordinates_tooltips(rows)] * tables_count,
        [page] * tables_count,
        [None] * tables_count,
        [[]] * tables_count,
        [summary] * len(callback_context.outputs_list[5]),
    )


@app.callback(
    Output("genome-focus", "value"),
    Input("selection-store", "data"),
    prevent_initial_call=True,
)
def focus_selected_mutation(selection):
    """Selecting a mutation with known coordinates brings it into view."""
    if not (selection or {}).get("position"):
        raise PreventUpdate
    return genome_view.FOCUS_MUTATION


@app.callback(
    Output("jbrowse-view", "children"),
    [Input("genome-focus", "value"),
     Input("genome-tracks", "value"),
     Input("genome-drug-tracks", "value")],
    [State("selection-store", "data"),
     State("current-gene-store", "data")],
    prevent_initial_call=True,
)
def update_genome_browser(focus, visible_tracks, drug_tracks, selection, gene_data):
    """Rebuild the browser for a new focus or track selection."""
    if not gene_data:
        raise PreventUpdate

    gene_info = data_loader.get_gene_info(gene_data["locus_tag"])
    if gene_info is None:
        raise PreventUpdate

    return genome_view.jbrowse_view(
        browser,
        gene_info,
        focus=focus,
        visible_tracks=visible_tracks or [],
        drug_tracks=drug_tracks,
        selection=selection,
    )


@app.callback(
    [Output("sequence-output", "children"),
     Output("sequence-meta", "children"),
     Output("sequence-clipboard", "content")],
    [Input("sequence-kind", "value"),
     Input("sequence-flank", "value")],
    State("current-gene-store", "data"),
)
def render_sequence(kind, flank, gene_data):
    """Retrieve the coding sequence, protein or flanking regions (TASK-09.3)."""
    if not gene_data:
        raise PreventUpdate

    gene_info = data_loader.get_gene_info(gene_data["locus_tag"])
    if gene_info is None:
        raise PreventUpdate

    flank = int(flank or 0)

    if kind == "protein":
        protein = (data_loader.translate_gene(gene_info) or "").rstrip("*")
        meta = f"Protein translation · {len(protein):,} aa · bacterial codon table 11"
        return genome_view.format_sequence(protein, block=10, per_line=60), meta, protein

    if kind == "cds":
        sequence = data_loader.get_gene_sequence(gene_info)
        meta = (
            f"Coding sequence · {len(sequence):,} bp · "
            f"{gene_info.chromosome}:{gene_info.start:,}-{gene_info.end:,} "
            f"({gene_info.strand} strand, shown 5'→3')"
        )
    elif kind == "cds_flank":
        sequence = data_loader.get_gene_sequence(gene_info, flank=flank)
        meta = (
            f"Coding sequence ± {flank:,} bp · {len(sequence):,} bp · "
            f"{gene_info.chromosome}:{max(1, gene_info.start - flank):,}-"
            f"{gene_info.end + flank:,} ({gene_info.strand} strand, shown 5'→3')"
        )
    elif kind == "upstream":
        if gene_info.strand == "+":
            start, end = gene_info.start - flank, gene_info.start - 1
        else:
            start, end = gene_info.end + 1, gene_info.end + flank
        sequence = data_loader.get_sequence(start, end, gene_info.chromosome)
        if gene_info.strand == "-":
            sequence = _reverse_complement(sequence)
        meta = f"Upstream flank · {len(sequence):,} bp · {flank:,} bp requested"
    elif kind == "downstream":
        if gene_info.strand == "+":
            start, end = gene_info.end + 1, gene_info.end + flank
        else:
            start, end = gene_info.start - flank, gene_info.start - 1
        sequence = data_loader.get_sequence(start, end, gene_info.chromosome)
        if gene_info.strand == "-":
            sequence = _reverse_complement(sequence)
        meta = f"Downstream flank · {len(sequence):,} bp · {flank:,} bp requested"
    else:
        raise PreventUpdate

    return genome_view.format_sequence(sequence), meta, sequence


def _reverse_complement(sequence: str) -> str:
    from data_utils import _reverse_complement as rc

    return rc(sequence)


# ----------------------------------------------------------------------
# Coordinate analysis
# ----------------------------------------------------------------------
@app.callback(
    Output(dict(ANALYSIS_DISPLAY_ID, index=ALL), "children"),
    Input("selection-store", "data"),
    State("current-gene-store", "data"),
)
def display_coordinate_calculation(selection, gene_data):
    """Coordinate analysis for the shared selection (TASK-04, TASK-08)."""
    panels = len(callback_context.outputs_list)
    if not panels:
        raise PreventUpdate
    if not selection or not gene_data:
        return [_analysis_placeholder()] * panels

    variant = selection.get("variant")
    mutation = selection.get("mutation")
    if not variant and not mutation:
        return [_analysis_placeholder()] * panels

    gene_name = gene_data["gene_name"]
    gene_info = GeneInfo(
        gene_id="",
        gene_name=gene_name,
        locus_tag=gene_data["locus_tag"],
        start=gene_data["start"],
        end=gene_data["end"],
        strand=gene_data["strand"],
        product=gene_data.get("product", ""),
        chromosome=gene_data["chromosome"],
    )

    position = selection.get("position")
    reference = selection.get("reference") or "N/A"
    alternative = selection.get("alternative") or "N/A"
    relative = (
        coord_calculator.calculate_relative_position(int(position), gene_info)
        if position else None
    )
    amino_acid = coord_calculator.parse_p_notation(mutation) if mutation else None

    resistance_df = data_loader.get_drug_resistance_info(
        gene_data["locus_tag"], variant or mutation
    )
    if len(resistance_df) == 0 and mutation:
        if any(term in mutation.lower() for term in ("fs", "stop", "ins", "del", "?")):
            resistance_df = data_loader.get_drug_resistance_info(
                gene_data["locus_tag"], "LoF"
            )

    panel = html.Div([
        dbc.Row([
            dbc.Col([
                html.H5(variant or mutation, className="text-primary fw-bold mb-1"),
                html.Small(
                    _selection_provenance(selection),
                    className="text-muted",
                ),
            ]),
            dbc.Col(
                # The browser has already moved to the mutation; this only
                # scrolls back up to it.
                html.A(
                    [html.I(className="bi bi-eye me-1"), "View in genome browser"],
                    href="#genome-browser-card",
                    className="btn btn-sm btn-soft",
                ),
                width="auto",
            ) if position else None,
        ], className="mb-3 align-items-start"),

        dbc.Row([
            dbc.Col([
                _analysis_field("Genomic Position", f"{position:,}" if position else "N/A"),
                _analysis_field(
                    "Gene Context", "{} ({} strand)".format(gene_name, gene_data["strand"])
                ),
                _allele_field(reference, alternative),
            ], md=6),
            dbc.Col([
                _analysis_field(
                    "Gene Relative (c.)",
                    f"c.{relative}" if relative else "N/A",
                    value_class="fw-bold text-success",
                ),
                _analysis_field(
                    "Amino Acid Position",
                    f"{amino_acid}" if amino_acid else "—",
                ),
            ], md=6),
        ], className="mb-4 bg-light p-3 rounded-3 gy-3"),

        html.H6(
            "Resistance Association",
            className="fw-bold mb-2 small text-uppercase letter-spacing-1",
        ),
        _resistance_summary(resistance_df),

        dbc.Accordion([
            dbc.AccordionItem(
                html.Code(
                    _derivation_text(gene_data, position, relative),
                    style={
                        "whiteSpace": "pre-wrap",
                        "display": "block",
                        "padding": "15px",
                        "fontSize": "0.85rem",
                    },
                ),
                title="View Mathematical Derivation",
                className="mt-4 border-0",
            )
        ], start_collapsed=True, className="mt-4"),
    ])

    return [panel] * panels


def _selection_provenance(selection: Dict) -> str:
    """Explain which nucleotide change is being analysed."""
    count = selection.get("coordinate_count")
    if selection.get("resolved_from_catalogue") and count:
        if count > 1:
            return (
                f"Showing the first of {count:,} nucleotide changes recorded for "
                "this variant — select a row in Genomic Coordinates for a specific one."
            )
        return "Single nucleotide change recorded for this variant."
    if count and count > 1:
        return f"{count:,} nucleotide changes recorded for this variant."
    return ""


def _analysis_field(label: str, value: str, value_class: str = "fw-bold") -> html.Div:
    return html.Div([
        html.Label(label, className="text-muted small d-block"),
        html.Span(value, className=value_class),
    ], className="mb-3")


def _allele_field(reference: str, alternative: str) -> html.Div:
    """
    Nucleotide change, with a copy button for long indel alleles.

    Large deletions carry alleles thousands of bases long, so the full
    sequence is offered through the clipboard rather than being rendered
    inline (TASK-07).
    """
    long_allele = max(len(str(reference)), len(str(alternative))) > 40
    summary = f"{_shorten(reference)} > {_shorten(alternative)}"

    children = [
        html.Label("Nucleotide Change", className="text-muted small d-block"),
        html.Span(summary, className="fw-bold allele-summary"),
    ]

    if long_allele:
        children.extend([
            html.Div([
                html.Details([
                    html.Summary("Full alleles", className="small text-primary"),
                    html.Div([
                        html.Div("REF", className="small text-muted mt-2"),
                        html.Pre(str(reference), id="allele-ref", className="sequence-output allele-block"),
                        dcc.Clipboard(
                            target_id="allele-ref",
                            title="Copy reference allele",
                            className="copy-btn",
                        ),
                        html.Div("ALT", className="small text-muted mt-2"),
                        html.Pre(str(alternative), id="allele-alt", className="sequence-output allele-block"),
                        dcc.Clipboard(
                            target_id="allele-alt",
                            title="Copy alternative allele",
                            className="copy-btn",
                        ),
                    ]),
                ]),
            ], className="mt-2"),
        ])

    return html.Div(children, className="mb-3")


def _shorten(allele: str, limit: int = 20) -> str:
    text = str(allele)
    if len(text) <= limit:
        return text
    return f"{text[:limit]}… [{len(text):,} bp]"


def _resistance_summary(resistance_df: pd.DataFrame) -> html.Div:
    """Drug-by-drug grading summary for the selected variant."""
    if len(resistance_df) == 0:
        return html.P(
            "No drug resistance data found in the catalogue for this variant.",
            className="text-muted small",
        )

    rows = []
    for _, entry in resistance_df.iterrows():
        grading = str(entry.get("FINAL CONFIDENCE GRADING", "") or "N/A")
        comment = entry.get("Comment")
        changes = entry.get("CHANGES vs ver1")

        rows.append(html.Div([
            dbc.Row([
                dbc.Col([
                    html.Strong(entry.get("drug", "—"), className="text-dark"),
                    html.Br(),
                    html.Small([
                        f"Tier {entry.get('tier', 'N/A')}",
                        html.Span(
                            f" · {entry.get('effect')}" if pd.notna(entry.get("effect")) else "",
                            className="text-muted",
                        ),
                    ], className="text-muted"),
                ], md=5),
                dbc.Col([
                    dbc.Badge(
                        grading,
                        className="badge-tb" if grading.startswith(("1)", "2)")) else "badge-custom",
                    ),
                    html.Div(
                        html.Small(changes, className="text-danger"),
                    ) if pd.notna(changes) and str(changes) not in ("", "No change") else None,
                ], md=7, className="text-md-end"),
            ], className="align-items-center g-2"),
            html.Div(
                html.Small([html.I(className="bi bi-exclamation-circle me-1"), comment],
                           className="text-danger fst-italic"),
                className="mt-1",
            ) if pd.notna(comment) and str(comment).strip() else None,
        ], className="py-2 border-bottom"))

    return html.Div(rows, className="mt-2")


def _derivation_text(gene_data: Dict, position: Optional[int], relative: Optional[int]) -> str:
    """Human-readable derivation of the gene-relative coordinate."""
    strand = gene_data["strand"]
    start = gene_data["start"]
    end = gene_data["end"]

    # HGVS has no c.0, so upstream positions drop the + 1.
    upstream = relative is not None and relative < 0
    plus_one = "" if upstream else " + 1"
    if strand == "+":
        formula = "c. position = genomic - gene_start" + plus_one
        calculation = (
            "{:,} - {:,}{} = {}".format(position, start, plus_one, relative)
            if position and relative is not None else "n/a"
        )
    else:
        formula = "c. position = gene_end - genomic" + plus_one
        calculation = (
            "{:,} - {:,}{} = {}".format(end, position, plus_one, relative)
            if position and relative is not None else "n/a"
        )
    if upstream:
        formula += "\n(upstream of the gene: HGVS has no c.0, so c.-1 is the base before c.1)"

    return (
        "Gene: {} ({})\n"
        "Span: {:,}-{:,} ({} strand)\n\n"
        "Formula ({} strand):\n{}\n\n"
        "Calculation:\n{}\n\n"
        "Note: M. tuberculosis is a prokaryote, so the coding sequence is\n"
        "contiguous and no intron offsets apply."
    ).format(
        gene_data.get("gene_name", ""), gene_data.get("locus_tag", ""),
        start, end, strand, strand, formula, calculation,
    )


if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=8050)
