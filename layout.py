"""
Static layout, branding and modal content for the DR-TBAtlas.

Institutional identity, attribution and citation material live here so that
they are defined once and reused by the header and the footer (TASK-15).
"""
import os
from typing import Dict, List, Optional

import dash_bootstrap_components as dbc
from dash import dcc, html

from data_utils import FIRST_LINE_DRUGS, GeneInfo

APP_NAME = "DR-TBAtlas"
APP_VERSION = "1.5.1"

# --- Institutional identity -------------------------------------------------
LAPAM_NAME = "Laboratory of Applied Research in Mycobacteria (LaPAM)"
LAPAM_URL = "https://lapam-usp.github.io/"
USP_NAME = "University of São Paulo"
USP_URL = "https://www.usp.br/"
ATTRIBUTION = (
    "Developed by Fernando Falat, Laboratory of Applied Research in "
    "Mycobacteria (LaPAM), University of São Paulo, Brazil."
)

# --- Repositories and licence ----------------------------------------------
PROJECT_REPO_URL = "https://github.com/falatfernando/dr-tbatlas"
PERSONAL_GITHUB_URL = "https://github.com/falatfernando"
# Set to the laboratory's GitHub organisation to surface a second repository
# icon in the header and footer.
LAB_GITHUB_URL: Optional[str] = None
LICENSE_NAME = "GPL-3.0"
LICENSE_URL = "https://www.gnu.org/licenses/gpl-3.0.html"

# --- Citation ---------------------------------------------------------------
DOI = "10.5281/zenodo.21035590"
DOI_URL = f"https://doi.org/{DOI}"
CITATION_YEAR = "2026"
CITATION_TEXT = (
    f"Falat Rangel, F., & Guimarães, A. M. de S. ({CITATION_YEAR}). "
    f"DR-TBAtlas: an interactive web-based interface for visual exploration of the "
    f"WHO Mycobacterium tuberculosis drug resistance catalogue (Version {APP_VERSION}) "
    f"[Computer software]. Zenodo. {DOI_URL}"
)
CITATION_BIBTEX = f"""@software{{falat_drtbatlas,
  author    = {{Falat Rangel, Fernando and Guimar\\~aes, Ana Marcia de S\\'a}},
  title     = {{DR-TBAtlas: an interactive web-based interface for visual exploration
               of the WHO \\emph{{Mycobacterium tuberculosis}} drug resistance catalogue}},
  version   = {{{APP_VERSION}}},
  year      = {{{CITATION_YEAR}}},
  publisher = {{Zenodo}},
  doi       = {{{DOI}}},
  url       = {{{DOI_URL}}}
}}"""

REFERENCES = [
    {
        "label": "WHO mutation catalogue (2nd edition, 2023)",
        "text": (
            "World Health Organization. Catalogue of mutations in Mycobacterium "
            "tuberculosis complex and their association with drug resistance, "
            "second edition. Geneva: WHO; 2023."
        ),
        "url": "https://www.who.int/publications/i/item/9789240082410",
    },
    {
        "label": "JBrowse 2",
        "text": (
            "Diesh C, et al. JBrowse 2: a modular genome browser with views of "
            "synteny and structural variation. Genome Biology. 2023;24:74."
        ),
        "url": "https://doi.org/10.1186/s13059-023-02914-z",
    },
    {
        "label": "H37Rv reference genome",
        "text": "Mycobacterium tuberculosis H37Rv, RefSeq assembly NC_000962.3.",
        "url": "https://www.ncbi.nlm.nih.gov/nuccore/NC_000962.3",
    },
]

SEARCH_PLACEHOLDER = (
    "Search gene ID (e.g., Rv0677c), gene name (e.g., mmpS5), "
    "or variant (e.g., katG_Ser315Thr)..."
)

QUICK_SEARCH_GENES = ("rpoB", "gyrA", "katG", "Rv0678")


def _asset_exists(app, filename: str) -> bool:
    """Whether an optional asset has been provided in the assets folder."""
    assets_dir = getattr(app, "config", None) and app.config.get("assets_folder")
    if not assets_dir:
        return False
    return os.path.exists(os.path.join(assets_dir, filename))


def usp_mark(app, height: int = 38) -> html.A:
    """
    USP institutional mark.

    Renders ``assets/usp.png`` when the official file has been added; falls
    back to a typographic wordmark so the attribution is never broken.
    """
    if _asset_exists(app, "usp.png"):
        content = html.Img(
            src=app.get_asset_url("usp.png"),
            alt=f"{USP_NAME} logo",
            style={"height": f"{height}px"},
            className="usp-logo",
        )
    else:
        content = html.Span([
            html.Span("USP", className="usp-wordmark-initials"),
            html.Span("Universidade de São Paulo", className="usp-wordmark-name"),
        ], className="usp-wordmark")

    return html.A(
        content,
        href=USP_URL,
        target="_blank",
        rel="noopener noreferrer",
        title=USP_NAME,
        className="text-decoration-none",
    )


def _github_links(app, icon_height: int = 22) -> List[dbc.NavItem]:
    """Repository links for the header."""
    links = [
        (PROJECT_REPO_URL, "code.png", "DR-TBAtlas source code"),
        (PERSONAL_GITHUB_URL, "github.png", "Developer profile on GitHub"),
    ]
    if LAB_GITHUB_URL:
        links.append((LAB_GITHUB_URL, "github.png", f"{LAPAM_NAME} on GitHub"))

    return [
        dbc.NavItem(
            dbc.NavLink(
                html.Img(
                    src=app.get_asset_url(icon),
                    alt=title,
                    title=title,
                    style={"height": f"{icon_height}px"},
                ),
                href=url,
                target="_blank",
                external_link=True,
                className="d-flex align-items-center",
            )
        )
        for url, icon, title in links
    ]


def navbar(app) -> dbc.Navbar:
    """Application navigation bar with branding and the new entry points."""
    return dbc.Navbar(
        dbc.Container([
            html.Div([
                html.A(
                    html.Img(
                        src=app.get_asset_url("lapam.png"),
                        alt=f"{LAPAM_NAME} logo",
                        className="lab-logo",
                    ),
                    href=LAPAM_URL,
                    target="_blank",
                    rel="noopener noreferrer",
                    title=f"{LAPAM_NAME} — official website",
                    className="d-inline-flex align-items-center",
                ),
                html.A(
                    dbc.NavbarBrand(
                        [html.I("Mycobacterium tuberculosis"),
                         " Genomic Resistance Explorer"],
                        className="navbar-brand-text",
                    ),
                    href="/",
                    style={"textDecoration": "none"},
                ),
            ], className="d-flex align-items-center gap-3"),

            dbc.NavbarToggler(id="navbar-toggler", n_clicks=0),
            dbc.Collapse(
                dbc.Nav([
                    dbc.NavItem(dbc.NavLink("Home", href="/")),
                    dbc.NavItem(
                        dbc.NavLink(
                            [html.I(className="bi bi-capsule me-2"), "Browse by Drug"],
                            id="open-drugs-modal",
                            href="#",
                            className="nav-action",
                        )
                    ),
                    dbc.NavItem(
                        dbc.NavLink(
                            [html.I(className="bi bi-bar-chart-line me-2"), "Catalogue Summary"],
                            id="open-summary-modal",
                            href="#",
                            className="nav-action",
                        )
                    ),
                    dbc.NavItem(
                        dbc.Button(
                            [html.I(className="bi bi-quote me-2"), "How to Cite"],
                            id="open-cite-modal",
                            className="cite-btn ms-lg-2",
                            size="sm",
                        ),
                        className="d-flex align-items-center",
                    ),
                    *_github_links(app),
                ], className="ms-auto align-items-lg-center", navbar=True),
                id="navbar-collapse",
                is_open=False,
                navbar=True,
            ),
        ], fluid=True),
        className="navbar-custom sticky-top mb-0",
    )


def hero_search(app) -> html.Div:
    """
    Hero section with the search form.

    The input sits inside a form so that pressing ``Enter`` submits the query
    without reloading the page (TASK-01).
    """
    return html.Div([
        dbc.Container([
            dbc.Row([
                dbc.Col([
                    # The logo carries the wordmark, so its alt text is the
                    # heading's accessible name.
                    html.H1(
                        html.Img(
                            src=app.get_asset_url("dr-tbatlas-logo.png"),
                            alt=APP_NAME,
                            className="hero-logo",
                        ),
                        className="mb-3",
                    ),
                    html.P([
                        "Explore ",
                        html.I("M. tuberculosis"),
                        " WHO catalogue mutations and visualize the genomic "
                        "region with precision.",
                    ], className="lead text-muted mb-5"),

                    dbc.Card([
                        dbc.CardBody([
                            dbc.Form([
                                dbc.InputGroup([
                                    dbc.Input(
                                        id="search-input",
                                        type="search",
                                        placeholder=SEARCH_PLACEHOLDER,
                                        className="search-input",
                                        autoFocus=True,
                                    ),
                                    dbc.Button(
                                        html.I(className="bi bi-search"),
                                        id="search-button",
                                        type="submit",
                                        className="search-btn",
                                        title="Search",
                                    ),
                                ]),
                            ], id="search-form", prevent_default_on_submit=True),

                            html.Div([
                                html.Span("Accepts: ", className="text-muted small"),
                                html.Span("gene ID ", className="small fw-semibold"),
                                html.Code("Rv0677c", className="search-hint-code"),
                                html.Span(" · gene name ", className="small fw-semibold"),
                                html.Code("mmpS5", className="search-hint-code"),
                                html.Span(" · variant ", className="small fw-semibold"),
                                html.Code("katG_Ser315Thr", className="search-hint-code"),
                            ], className="mt-3 small text-center search-hint"),

                            html.Div([
                                html.Span("Quick access: ", className="text-muted small me-2"),
                                *[
                                    dbc.Button(
                                        gene,
                                        id={"type": "quick-search", "index": gene},
                                        size="sm",
                                        className="quick-search-btn",
                                    )
                                    for gene in QUICK_SEARCH_GENES
                                ],
                            ], className="mt-3 d-flex align-items-center justify-content-center flex-wrap"),

                            html.Div(id="search-feedback", className="mt-3"),
                        ])
                    ], className="search-card"),
                ], lg=9, className="mx-auto"),
            ])
        ], className="py-5"),
    ], className="hero-section")


def welcome_panel(totals: Dict[str, int]) -> html.Div:
    """
    Landing state shown before a gene is selected.

    Gives the three accepted query formats a worked example each and points at
    the two catalogue-wide entry points, rather than leaving the page empty.
    """
    examples = [
        ("bi-upc-scan", "Gene ID / locus tag", "Rv0677c",
         "The H37Rv identifier, as used in the reference annotation."),
        ("bi-tag", "Gene name", "mmpS5",
         "The gene symbol, as used in the WHO catalogue."),
        ("bi-crosshair", "Variant", "katG_Ser315Thr",
         "Gene, underscore, then the mutation in three-letter notation."),
    ]

    return html.Div([
        dbc.Row([
            dbc.Col(
                html.Div([
                    html.I(className=f"bi {icon} welcome-icon"),
                    html.H6(title, className="fw-bold mb-1"),
                    dbc.Button(
                        html.Code(example, className="welcome-example"),
                        id={"type": "quick-search", "index": example},
                        className="welcome-chip mb-2",
                        size="sm",
                    ),
                    html.P(description, className="text-muted small mb-0"),
                ], className="welcome-card"),
                md=4,
            )
            for icon, title, example, description in examples
        ], className="g-3"),

        dbc.Row([
            dbc.Col(html.Div([
                html.P([
                    html.Strong("Not sure where to start? "),
                    "Open ",
                    html.Strong("Browse by Drug"),
                    " to pick a gene from a drug's catalogue entries, or ",
                    html.Strong("Catalogue Summary"),
                    " for counts across the whole dataset.",
                ], className="small mb-2"),
                html.P([
                    f"The loaded catalogue covers {totals['variants']:,} distinct "
                    f"variants across {totals['genes']:,} genes and "
                    f"{totals['drugs']:,} drugs, of which "
                    f"{totals['1) Assoc w R'] + totals['2) Assoc w R - Interim']:,} "
                    "are graded as associated with resistance.",
                ], className="text-muted small mb-0"),
            ], className="welcome-note"))
        ], className="mt-3"),
    ], className="welcome-panel mb-5")


def mutation_advisory(message: str) -> dbc.Alert:
    """Advisory shown when a mutation is searched without its gene (TASK-03)."""
    return dbc.Alert([
        html.I(className="bi bi-info-circle-fill me-2"),
        html.Span(message),
    ], color="warning", className="rounded-3 shadow-sm mb-0 text-start")


def gene_not_found_alert(query: str) -> dbc.Alert:
    """Alert shown when a query matches no annotated gene."""
    return dbc.Alert([
        html.I(className="bi bi-exclamation-triangle-fill me-2"),
        f"'{query}' did not match any gene in the H37Rv annotation. ",
        "Try a gene name (katG), a locus tag (Rv1908c) or a variant "
        "(katG_Ser315Thr).",
    ], color="warning", className="rounded-3 shadow-sm mb-0 text-start")


def non_catalogue_card(gene_info: GeneInfo) -> dbc.Card:
    """
    Status card for a gene that is annotated but not catalogued (TASK-05).

    Replaces the variant and coordinate tables, which would otherwise prompt
    for a row selection that cannot exist.
    """
    return dbc.Card([
        dbc.CardHeader([
            html.H5(
                [html.I(className="bi bi-info-circle me-2"), "Gene Not in WHO Catalogue"],
                className="card-header-title",
            )
        ], className="card-header-custom"),
        dbc.CardBody([
            html.P([
                html.Strong(f"{gene_info.display_name} ({gene_info.locus_tag})"),
                " is part of the ",
                html.I("M. tuberculosis"),
                " reference genome (H37Rv) but is not currently indexed in the "
                "WHO catalogue of drug resistance-associated mutations.",
            ], className="mb-2"),
            html.P([
                "The genome browser and gene annotation above remain fully "
                "available. Drug resistance and coordinate analysis are only "
                "offered for catalogued genes — use ",
                html.Strong("Browse by Drug"),
                " in the header to see which genes the catalogue covers.",
            ], className="text-muted small mb-0"),
        ]),
    ], className="result-card info-card")


def browse_by_drug_modal(drug_map: Dict[str, List[Dict]]) -> dbc.Modal:
    """Browse the catalogue starting from a drug name (TASK-13)."""
    items = []
    for drug, genes in drug_map.items():
        is_first_line = drug in FIRST_LINE_DRUGS
        chips = [
            dbc.Button([
                html.Span(entry["gene"], className="fw-semibold"),
                html.Span(
                    f" · {entry['locus_tag']}" if entry["locus_tag"] != entry["gene"] else "",
                    className="text-muted small",
                ),
                dbc.Badge(
                    f"{entry['variant_count']:,}",
                    className="ms-2 badge-count",
                    title=f"{entry['variant_count']:,} catalogued mutations for {drug}",
                ),
            ],
                id={
                    "type": "gene-link",
                    "drug": drug,
                    "index": entry["locus_tag"] or entry["gene"],
                },
                className="gene-chip",
                size="sm",
            )
            for entry in genes
        ]

        items.append(
            dbc.AccordionItem(
                html.Div(chips, className="d-flex flex-wrap gap-2"),
                title=html.Span([
                    drug,
                    dbc.Badge(
                        "first-line" if is_first_line else "second-line / new",
                        className="ms-2 badge-line",
                    ),
                    html.Span(f"{len(genes)} genes", className="ms-2 text-muted small"),
                ]),
                item_id=f"drug-{drug}",
            )
        )

    return dbc.Modal([
        dbc.ModalHeader(dbc.ModalTitle([
            html.I(className="bi bi-capsule me-2"), "Browse by Drug",
        ])),
        dbc.ModalBody([
            html.P(
                "Select a drug to see the genes the WHO catalogue associates with "
                "it, then click a gene to load it in the explorer.",
                className="text-muted small",
            ),
            dbc.Accordion(items, start_collapsed=True, always_open=False, flush=True),
        ]),
        dbc.ModalFooter(
            dbc.Button("Close", id="close-drugs-modal", className="btn-soft", size="sm")
        ),
    ], id="drugs-modal", size="lg", scrollable=True, is_open=False)


def summary_modal() -> dbc.Modal:
    """Catalogue overview and summary statistics (TASK-14)."""
    return dbc.Modal([
        dbc.ModalHeader(dbc.ModalTitle([
            html.I(className="bi bi-bar-chart-line me-2"), "WHO Catalogue Summary",
        ])),
        dbc.ModalBody(
            dcc.Loading(html.Div(id="summary-modal-body"), type="default"),
        ),
        dbc.ModalFooter(
            dbc.Button("Close", id="close-summary-modal", className="btn-soft", size="sm")
        ),
    ], id="summary-modal", size="xl", scrollable=True, is_open=False)


def cite_modal() -> dbc.Modal:
    """Citation instructions with copy-to-clipboard support (TASK-15.5)."""
    return dbc.Modal([
        dbc.ModalHeader(dbc.ModalTitle([
            html.I(className="bi bi-quote me-2"), "How to Cite",
        ])),
        dbc.ModalBody([
            html.P(
                "If DR-TBAtlas supported your work, please cite the software and "
                "the underlying data sources.",
                className="text-muted small",
            ),

            html.H6("Software citation", className="fw-bold small text-uppercase mt-3"),
            html.Div([
                html.Div(CITATION_TEXT, id="citation-text", className="citation-block"),
                html.Div([
                    dcc.Clipboard(
                        target_id="citation-text",
                        title="Copy citation",
                        className="copy-btn",
                    ),
                    html.Span("Copy citation", className="small text-muted ms-2"),
                ], className="d-flex align-items-center mt-2"),
            ]),

            html.H6("BibTeX", className="fw-bold small text-uppercase mt-4"),
            html.Div([
                html.Pre(CITATION_BIBTEX, id="citation-bibtex", className="citation-block mb-0"),
                html.Div([
                    dcc.Clipboard(
                        target_id="citation-bibtex",
                        title="Copy BibTeX",
                        className="copy-btn",
                    ),
                    html.Span("Copy BibTeX", className="small text-muted ms-2"),
                ], className="d-flex align-items-center mt-2"),
            ]),

            html.H6("Data and software sources", className="fw-bold small text-uppercase mt-4"),
            html.Ul([
                html.Li([
                    html.Strong(reference["label"]),
                    ": ",
                    reference["text"],
                    " ",
                    html.A(
                        "link",
                        href=reference["url"],
                        target="_blank",
                        rel="noopener noreferrer",
                    ),
                ], className="small mb-2")
                for reference in REFERENCES
            ], className="mb-0 ps-3"),
        ]),
        dbc.ModalFooter([
            html.A(
                dbc.Badge(f"DOI {DOI}", className="badge-custom me-2"),
                href=DOI_URL,
                target="_blank",
                rel="noopener noreferrer",
                className="text-decoration-none",
            ),
            dbc.Button("Close", id="close-cite-modal", className="btn-soft", size="sm"),
        ]),
    ], id="cite-modal", size="lg", scrollable=True, is_open=False)


def footer(app) -> html.Footer:
    """Institutional attribution, repositories and licence notice (TASK-15)."""
    repo_links = [
        html.A([
            html.Img(src=app.get_asset_url("code.png"), alt="", className="footer-icon"),
            "Project repository",
        ], href=PROJECT_REPO_URL, target="_blank", rel="noopener noreferrer",
            className="footer-link"),
        html.A([
            html.Img(src=app.get_asset_url("github.png"), alt="", className="footer-icon"),
            "Developer profile",
        ], href=PERSONAL_GITHUB_URL, target="_blank", rel="noopener noreferrer",
            className="footer-link"),
    ]
    if LAB_GITHUB_URL:
        repo_links.append(
            html.A([
                html.Img(src=app.get_asset_url("github.png"), alt="", className="footer-icon"),
                "Laboratory repository",
            ], href=LAB_GITHUB_URL, target="_blank", rel="noopener noreferrer",
                className="footer-link")
        )

    return html.Footer([
        html.Hr(className="mb-4"),
        dbc.Row([
            dbc.Col([
                html.Div([
                    html.A(
                        html.Img(
                            src=app.get_asset_url("lapam.png"),
                            alt=f"{LAPAM_NAME} logo",
                            className="footer-lab-logo",
                        ),
                        href=LAPAM_URL,
                        target="_blank",
                        rel="noopener noreferrer",
                        title=f"{LAPAM_NAME} — official website",
                    ),
                    usp_mark(app),
                ], className="d-flex align-items-center gap-4 mb-3 flex-wrap"),
                html.P(ATTRIBUTION, className="fw-semibold mb-2 footer-attribution"),
                html.P([
                    f"DR-TBAtlas v{APP_VERSION} · ",
                    html.A(
                        dbc.Badge(LICENSE_NAME, className="badge-license"),
                        href=LICENSE_URL,
                        target="_blank",
                        rel="noopener noreferrer",
                        className="text-decoration-none",
                    ),
                    " · ",
                    html.A(
                        dbc.Badge(f"DOI {DOI}", className="badge-doi"),
                        href=DOI_URL,
                        target="_blank",
                        rel="noopener noreferrer",
                        className="text-decoration-none",
                    ),
                ], className="small mb-2 d-flex align-items-center gap-2 flex-wrap"),
            ], lg=7),

            dbc.Col([
                html.Div(repo_links, className="d-flex flex-column align-items-lg-end gap-2"),
                dbc.Button(
                    [html.I(className="bi bi-quote me-2"), "How to cite"],
                    id="open-cite-modal-footer",
                    size="sm",
                    className="btn-soft mt-2",
                ),
            ], lg=5, className="text-lg-end"),
        ], className="gy-3"),

        html.Hr(className="my-4"),
        html.P([
            html.Strong("Data sources: "),
            "H37Rv reference genome (NC_000962.3); WHO Catalogue of mutations in ",
            html.I("Mycobacterium tuberculosis"),
            " complex and their association with drug resistance, 2nd edition. "
            "Visualization powered by JBrowse 2.",
        ], className="text-center small mb-0 text-muted"),
    ], className="footer")
