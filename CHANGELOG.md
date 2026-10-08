# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.5.1] - 2026-10-08

### Fixed
- **Upstream c. positions were off by one.** HGVS has no c.0, so the base
  just before c.1 is c.-1, but the coordinate calculator counted it as c.0:
  every upstream variant in the Genomic Coordinates table and the Coordinate
  Analysis panel read one position too high (`katG_c.-1T>C` showed as
  `c.0`, `fabG1_c.-15C>T` as `c.-14`), and c. → genomic conversion was off by
  one in the other direction. The derivation text now shows the upstream
  formula.
- Upstream positions of dnaA, which starts at position 1, wrap round the
  circular chromosome instead of being reported as c.4411221 and so on.
- Upstream c. variants no longer report an amino acid position.

### Added
- Coordinate calculator tests on both strands, around c.1, across the
  chromosome origin, and against every single-base c. variant in the WHO
  genomic coordinates file (11,587 changes).

## [1.5.0] - 2026-10-08

### Added
- **Mutation filters for the Drug Resistance Profile.** A filter panel above
  the drug tables replaces the per-column filter row, which expected query
  syntax typed into an unlabelled cell. It offers:
  - a *Find mutation* box that accepts one- or three-letter residues
    (`S450L`, `Ser450Leu`), a whole codon or nucleotide position (`450`),
    nucleotide changes (`c.-15C>T`), consequences (`fs`, `LoF`) and several
    terms separated by commas;
  - drop-down lists for drug and confidence grading, with quick picks for
    associated (1–2), uncertain (3) and not associated (4–5);
  - under *More filters*: effect, mutation type (p., c., n., gene-level),
    change vs catalogue v1, gene tier, a codon/nucleotide position range and
    switches for WHO comments, the relaxed-thresholds simulation and silent
    mutations.
  Every option lists only the values present for the active gene, with their
  counts. Filters apply to all drug tables at once; drugs left with no match
  are hidden, each drug's badge reads "n of N mutations", a summary line
  states what is shown, and *Clear filters* resets everything. The *More
  filters* button shows how many hidden filters are active.
- **Genomic Coordinates filters**: search by variant or genomic position, a
  change-type list (SNV, insertion, deletion, MNV) and a switch that limits
  the table to the variants left by the Drug Resistance Profile filters.

### Changed
- The Catalogue Summary tables keep their per-column filter row, now
  case-insensitive and labelled "Filter…".

## [1.4.0] - 2026-10-06

### Added
- **Desktop app for Windows, macOS and Linux.** `launcher.py` starts the
  server, opens the explorer in the default browser and shows a small window
  to reopen or quit it; `packaging/drtbatlas.spec` bundles it with PyInstaller
  into a double-click executable that needs no Python installation.
- `.github/workflows/release.yml` builds and smoke-tests every platform on
  GitHub Actions and attaches the builds to a draft release for version tags.

### Changed
- The genome browser track folder can be moved with the
  `DRTBATLAS_TRACKS_DIR` environment variable; the desktop app keeps it in a
  per-user cache, since its own folder may be read-only.

## [1.3.0] - 2026-09-28

### Changed
- **Renamed the project from MtbRx to DR-TBAtlas.** The application title,
  hero banner, footer, citation text and BibTeX key (`falat_drtbatlas`),
  `CITATION.cff`, manuscript, preprint abstract, README, contributing guide
  and package name (`dr-tbatlas`) all use the new name, and the repository
  URL is now `github.com/falatfernando/dr-tbatlas`.
- The hero banner shows the DR-TBAtlas logo, and the site has a favicon.
- `assets/mtbrx.js` is now `assets/dr-tbatlas.js`, and its marker attribute
  is `data-drtbatlas-hidden`.
- **The embedded JBrowse 2 view is now the genome explorer**, opened on the
  whole chromosome rather than locked to the gene window, with:
  - genes coloured by strand and biotype, whose details carry the product,
    functional note, protein length, WHO catalogue drugs and tier, and a link
    to the Mycobrowser gene page;
  - WHO catalogue variant tracks by confidence grading (groups 1–2, 3 and
    4–5), coloured by grading, and a resistance-associated variant track for
    each drug;
  - WHO catalogue genes by tier, catalogue variant density (1 kb bins) and GC
    content (100 bp windows);
  - a "Current selection" track marking the active gene and mutation;
  - a text index, so the browser's location box finds gene names, locus
    tags, product words and catalogue variants;
  - focus presets (gene, ±500 bp to ±50 kb, selected mutation), track
    toggles and a per-drug track picker above the browser, and a colour
    legend and usage guide below it;
  - the application's colour theme.
- Selecting a mutation in the resistance or coordinate tables zooms the
  browser to its codon; the Coordinate Analysis panel links back up to it.
- Track files are generated into `data/tracks/` at start-up and rebuilt only
  when the annotation, catalogue or track format changes. They are gzipped,
  which cuts the largest download from 8.6 MB to under 1 MB.

### Removed
- The Plotly gene neighbourhood track, which the genome browser supersedes.
  The previous/next gene shortcuts in the gene overview remain.
- `DataLoader.get_prokaryotic_gff3_path` and `DataLoader.get_jbrowse_config`,
  replaced by `browser_tracks.BrowserTracks`. The direct `plotly` dependency
  is dropped (Dash still installs it).

### Fixed
- Start-up with gffutils 0.14, which no longer provides `gffutils.Database`.

## [1.2.0] - 2026-09-23

### Changed
- **Renamed the project from TBDashboard to MtbRx**, to set it apart from the
  many epidemiological "TB dashboards" in the literature. The application
  title, header, hero banner, footer, citation text and BibTeX key
  (`falat_mtbrx`), manuscript, preprint abstract, README, contributing guide
  and package name (`mtbrx`) all use the new name.
- The repository moved to `github.com/falatfernando/mtbrx`; the old
  `tbdashboard` URL redirects.
- The header now reads "*Mycobacterium tuberculosis* Genomic Resistance
  Explorer" and wraps on narrow screens instead of overflowing.
- `assets/tbdashboard.js` is now `assets/mtbrx.js`, and its marker attribute
  is `data-mtbrx-hidden`.

### Added
- `CITATION.cff`, so GitHub shows a "Cite this repository" button and Zenodo
  archives releases with the correct title and authors.

## [1.1.0] - 2026-09-21

### Added
- **Multimodal search**: the search bar accepts a locus tag (`Rv0677c`), a gene
  symbol (`mmpS5`) or a full variant (`katG_Ser315Thr`), normalises input, and
  up-converts legacy one-letter notation (`katG_S315T`) to the three-letter
  form used by the WHO catalogue.
- **Enter to search**: the input is wrapped in a form, so pressing `Enter`
  submits the query without reloading the page.
- **Mutation-only guard**: a mutation typed without its gene (`Ser315Thr`) is
  blocked with an advisory explaining the required `gene_mutation` format.
- **Variant deep linking**: searching a full variant loads the gene and
  automatically selects, highlights and pages to that mutation in both the
  Drug Resistance Profile and the Genomic Coordinates tables.
- **Reactive selection**: the resistance table, coordinates table and
  Coordinate Analysis card share one selection, linked in both directions.
- **Gene Neighbourhood track**: an interactive, prokaryote-native track where
  clicking any gene loads it into the whole dashboard, and hovering shows gene
  symbol, locus tag, product, functional note and coordinates.
- **Sequence Retrieval panel**: coding sequence, protein translation
  (bacterial codon table 11) and adjustable upstream/downstream flanks
  (default ±500 bp), with copy-to-clipboard.
- **Non-catalogue notice**: genes annotated in H37Rv but absent from the WHO
  catalogue now show an explicit status card instead of an empty table prompt.
- **Protein length**: the gene overview reports length as `2,223 bp (740 aa)`.
- **Full WHO schema**: the resistance table gained `Effect`, `Comment`,
  `CHANGES vs ver1`, `Relaxed thresholds simulation` and `Silent mutation`,
  with a column-visibility control and emphasis on remarks and version changes.
- **Browse by Drug**: a catalogue browser grouping genes by drug, first-line
  drugs first, with per-gene mutation counts; clicking a gene loads it.
- **Catalogue Summary**: per-gene and per-drug summary tables with tier,
  confidence-grading and loss-of-function counts, plus aggregate totals.
- **Branding and citation**: enlarged, linked LaPAM logo, USP institutional
  mark, formal attribution, licence and DOI badges, and a "How to Cite" modal
  with copyable citation text and BibTeX.
- Tests for query parsing, the data layer and the table builders.

### Changed
- **Prokaryotic genome browser**: the annotation served to JBrowse is
  flattened to one gene-level feature per locus. *M. tuberculosis* has no
  splicing, so the redundant gene/CDS hierarchy and every intron control are
  removed at the source rather than hidden.
- **Genomic Coordinates table**: scrolls horizontally, shortens long indel
  alleles to a fixed width with the full sequence on hover, keeps cell text
  selectable, and offers clipboard copy for multi-kilobase alleles.
- Gene lookup is served from an in-memory annotation index instead of a
  per-query SQLite scan, and catalogue gene aliases are derived from the
  annotation rather than a hand-maintained list, so all 65 catalogue genes
  resolve from either identifier.

### Fixed
- `coordinate_calculator` referenced `DataLoader`, `GeneInfo` and `typing`
  names it never imported, which raised `NameError` on Python below 3.14.
- The coordinate derivation text used nested same-quote f-strings, a syntax
  error before Python 3.12.
- Commas in GFF3 attribute values are escaped, so a product such as
  "KatG,catalase-peroxidase" is no longer parsed as two values.

## [1.0.0] - 2026-06-29

### Added
- **Gene Search**: Interactive lookup for *Mycobacterium tuberculosis* genes with full annotation display (locus tag, coordinates, product).
- **JBrowse Integration**: Built-in genomic browser leveraging `dash-jbrowse` to visualize gene models, genomic tracks, and custom sequence alignments.
- **Coordinate Calculator**: Utility for translating between gene-relative (HGVS c. and p. notation) and absolute genomic coordinates.
- **WHO Catalogue Integration**: Preloaded WHO tuberculosis mutation catalogue for drug resistance annotations.
- **Packaging and CI Setup**: Included standard `setup.py` packaging, unit tests, and GitHub Actions CI workflow.
