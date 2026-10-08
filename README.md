# DR-TBAtlas - *Mycobacterium tuberculosis* Genomic Resistance Explorer

DR-TBAtlas is a web-based genomic explorer for *Mycobacterium tuberculosis*, designed to bridge the gap between complex genomic data and clinical drug resistance interpretation. Built with Python Dash, it provides researchers and clinicians with an integrated platform to visualize genomic regions, calculate variant coordinates, and cross-reference mutations with the official WHO drug resistance catalogue. This app has it's on [![DOI](https://zenodo.org/badge/1199917151.svg)](https://doi.org/10.5281/zenodo.21035590) and a free hosted version is available at [Render](https://tbdashboard.onrender.com/).
***Note: the free Render instance may take 30-60 seconds to wake up on first load.***

## References

- **JBrowse 2**: Diesh et al, 2023. JBrowse 2: a modular genome browser with views of synteny and structural variation. *Genome Biology* 24:74. [https://doi.org/10.1186/s13059-023-02914-z](https://doi.org/10.1186/s13059-023-02914-z)
- **WHO Mutation Catalogue (2023)**: WHO catalogue of mutations in *Mycobacterium tuberculosis* complex and their association with drug resistance, second edition. [Publication](https://www.who.int/publications/i/item/9789240082410) | [GitHub Repository](https://github.com/GTB-tbsequencing/mutation-catalogue-2023/tree/main)

## Features

- **Multimodal Search**: One box accepts three query formats, and `Enter` submits:
  - Locus tag / gene ID — `Rv0677c`
  - Gene name / symbol — `mmpS5`
  - Variant — `katG_Ser315Thr` (three-letter notation; `katG_S315T` is accepted and normalised)
- **Variant Deep Linking**: Searching a variant loads the gene *and* selects, highlights and pages to that mutation in the resistance and coordinate tables
- **Genome Browser**: An embedded JBrowse 2 view of the whole H37Rv chromosome, opened on the active gene, with:
  - Genes coloured by strand and biotype, each with product, functional note, WHO catalogue drugs and tier, and a link to its Mycobrowser page in the feature details
  - WHO catalogue variant tracks split by confidence grading (associated, uncertain, not associated), plus one resistance-associated variant track per drug
  - WHO catalogue genes by tier, catalogue variant density and GC content tracks
  - A "Current selection" track that marks the active gene and the selected mutation
  - Reference sequence with six-frame translation at base level
  - Focus presets (gene, ±500 bp to ±50 kb, selected mutation) and track toggles above the browser; selecting a mutation in any table zooms the browser to its codon
  - Search from the browser's own location box by gene name, locus tag, product word or catalogue variant (`katG_p.Ser315Thr`)
  - JBrowse's built-in tools: track selector, motif search, SVG export, horizontal flip and region sequence retrieval
- **Drug Resistance Profiles**: The full WHO catalogue schema — mutation, tier, final confidence grading, effect, comment, `CHANGES vs ver1`, relaxed-thresholds simulation and silent-mutation flag, with a column-visibility control
- **Mutation Filters**: A labelled filter panel above the resistance tables:
  - Free-text mutation search in one- or three-letter notation (`S450L` or `Ser450Leu`), by position (`450`) or with several terms at once
  - Drop-down lists for drug, confidence grading (with quick picks), effect, mutation type, change vs catalogue v1 and gene tier, each listing only the values present for the gene, with counts
  - A codon/nucleotide position range, e.g. the rpoB RRDR (codons 426–452)
  - Switches for WHO comments, the relaxed-thresholds simulation and silent mutations
  - Coordinates table filters by variant, genomic position and change type, optionally linked to the resistance filters
- **Coordinate Calculator**: Automatically convert between:
  - Genomic coordinates (absolute position on chromosome)
  - Gene-relative coordinates (c. notation, e.g., c.102G>A)
  - Amino acid positions (p. notation, e.g., p.Asp3Ala)
- **Sequence Retrieval**: Coding sequence, protein translation, and adjustable upstream/downstream flanks (default ±500 bp), with copy-to-clipboard
- **Browse by Drug**: Start from a drug name and jump to any of its catalogue genes
- **Catalogue Summary**: Per-gene and per-drug counts by tier, confidence grading and loss of function, with aggregate totals
- **Drill-down Details**: Explore all nucleotide changes for each mutation, including multi-kilobase indel alleles

## Desktop App (no command line)

Ready-to-run versions for Windows, macOS and Linux are attached to each
[GitHub Release](https://github.com/falatfernando/dr-tbatlas/releases).
Download the file for your system, unzip it and double-click **DR-TBAtlas**:
a small window opens and the explorer appears in your web browser. Step-by-step
instructions, including the one-time security prompt on Windows and macOS,
are on the release page.

### Building the desktop app

The executables are made with [PyInstaller](https://pyinstaller.org) from
`launcher.py`, which starts the server, opens the browser and shows a small
control window. PyInstaller builds only for the system it runs on, so
`.github/workflows/release.yml` builds every platform on GitHub Actions:
pushing a version tag (`git tag 1.4.0 && git push origin 1.4.0`) creates a draft
release with all builds attached. To build locally:

```bash
pip install -r requirements.txt -r packaging/requirements-build.txt
pyinstaller packaging/drtbatlas.spec
```

The result is in `dist/`. The packaged app writes its genome browser tracks
and `launcher.log` to a per-user cache folder (`%LOCALAPPDATA%\DR-TBAtlas`,
`~/Library/Caches/DR-TBAtlas` or `~/.cache/DR-TBAtlas`) instead of `data/`.

## Installation

### Prerequisites

- Python 3.8 or higher
- pip package manager

### Setup

1. **Create virtual environment** (recommended):
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

2. **Install dependencies**:
```bash
pip install -r requirements.txt
```

3. **Verify data files**:
Ensure the following files exist in the `data/` directory:
- `h37rv.gff3` - Genome annotation file
- `h37rv.fasta` - Reference genome sequence
- `h37rv.fasta.fai` - FASTA index file
- `catalogue_master_file.txt` - WHO drug resistance catalogue
- `genomic_coordinates.txt` - Mutation coordinates

## Running Tests

To run the automated tests, ensure you have installed the dependencies (including `pytest`), and execute the following command from the root directory:

```bash
pytest tests/
```

## Usage

### Running the Application

```bash
python app.py
```

The application will start on `http://localhost:8050`

### Quick Start

1. Open your browser to `http://localhost:8050`
2. Search by gene name (`katG`), locus tag (`Rv1908c`) or variant (`katG_Ser315Thr`)
3. Click the search button or press `Enter`
4. Explore the results:
   - **Gene Info**: Genomic coordinates, strand, length in bp and amino acids, functional note
   - **Genome Browser**: JBrowse view of the gene with catalogue, GC and selection tracks
   - **Sequence Retrieval**: CDS, protein and flanking sequence
   - **Drug Resistance**: Associated drugs, tiers and WHO confidence gradings
   - **Genomic Coordinates**: All nucleotide changes
   - **Coordinate Analysis**: Detailed position calculations for the selected mutation

Alternatively, use **Browse by Drug** in the header to reach a gene from a drug
name, or **Catalogue Summary** for counts across the whole dataset.

### A note on prokaryotic annotation

*M. tuberculosis* is a bacterium: it has no splicing, and a coding sequence is
equivalent to its gene. The annotation handed to the genome browser is
therefore flattened to a single gene-level feature per locus (see
`browser_tracks.py`), which removes
both the redundant gene/CDS pair and every intron-based control from the
browser's sequence tools. Protein translation and flanking-sequence retrieval
are provided by the **Sequence Retrieval** panel instead.

### Coordinate Calculation Feature

The application automatically calculates the relationship between:

- **Genomic Position**: Absolute position on the chromosome (e.g., position 8 on NC_000962.3)
- **Gene Relative Position**: Position relative to gene start (e.g., c.7 for dnaA)
- **Amino Acid Position**: Codon position in the protein (e.g., p.Asp3)

#### Example Calculation

For `dnaA_p.Asp3Ala`:
```
Gene: dnaA
  - Genomic coordinates: NC_000962.3:1-1524 (+ strand)
  - Gene start: 1
  
Mutation: p.Asp3Ala (Aspartic acid → Alanine at amino acid position 3)
  - Amino acid position: 3
  - Nucleotide position (relative): (3-1) × 3 + 1 = 7
  - Genomic position: 1 + 7 - 1 = 7
  
Result: The mutation is at genomic position 7, which is c.7 in HGVS notation
```

For genes on the **minus strand**, the calculation is reversed:
```
Gene on - strand:
  - Genomic position = Gene end - Relative position + 1
```

## Project Structure

```
dr-tbatlas/
├── app.py                      # Dash application and callbacks
├── launcher.py                 # Desktop app entry point: server, browser and control window
├── layout.py                   # Static layout, branding, modals and citation
├── genome_view.py              # Genome browser card, controls and sequence panel
├── browser_tracks.py           # JBrowse track files, text index and configuration
├── tables.py                   # Resistance, coordinate and summary tables
├── search_utils.py             # Query parsing and mutation normalisation
├── data_utils.py               # Data loading and parsing utilities
├── coordinate_calculator.py    # Coordinate conversion utilities
├── requirements.txt            # Python dependencies
├── README.md                   # This file
├── assets/
│   ├── style.css              # Application styling
│   ├── dr-tbatlas.js          # Removes residual intron controls from JBrowse dialogs
│   ├── dr-tbatlas-logo.png    # DR-TBAtlas logo (hero), from icon.png
│   ├── favicon.ico            # Browser tab icon
│   ├── lapam.png              # Laboratory logo
│   └── usp.png                # Optional: USP crest, shown in the footer when present
├── packaging/                  # PyInstaller spec, build requirements and release notes
├── tests/                      # pytest suite
└── data/
    ├── h37rv.gff3             # Genome annotation
    ├── tracks/                # Generated at start-up: genome browser tracks and search index
    ├── h37rv.fasta            # Reference sequence
    ├── h37rv.fasta.fai        # FASTA index
    ├── catalogue_master_file.txt  # Drug resistance catalogue
    └── genomic_coordinates.txt    # Mutation coordinates
```

## Data Sources

- **H37Rv Reference Genome**: NC_000962.3
- **GFF3 Annotation**: RefSeq annotations via NCBI
- **Drug Resistance Catalogue**: [WHO Catalogue (2023) second edition](https://www.who.int/publications/i/item/9789240082410) for *M. tuberculosis* drug resistance.

### DataLoader Class

```python
from data_utils import DataLoader

loader = DataLoader(data_dir="data")
loader.load_gff3()
loader.load_catalogue()
loader.load_genomic_coordinates()

# Get gene information
gene_info = loader.get_gene_info("dnaA")

# Search for mutations
mutations = loader.search_mutations_by_gene("dnaA")

# Get drug resistance info
resistance = loader.get_drug_resistance_info("dnaA")
```

### CoordinateCalculator Class

```python
from coordinate_calculator import CoordinateCalculator

calc = CoordinateCalculator(data_loader)

# Calculate all coordinates for a variant
result = calc.calculate_full_coordinates("dnaA_c.102G>A", "dnaA")
# Returns: {
#   'variant': 'dnaA_c.102G>A',
#   'gene': 'dnaA',
#   'gene_info': {...},
#   'mutation': 'c.102G>A',
#   'relative_nucleotide_position': 102,
#   'genomic_position': ...,
#   'amino_acid_position': ...
# }

# Get complete mutation information
mutation_info = calc.get_mutation_with_coordinates(
    variant="dnaA_c.102G>A",
    gene_name="dnaA",
    drug="Amikacin"
)
```

## Troubleshooting

### Common Issues

1. **"Gene not found" error**:
   - Check that the GFF3 file is properly loaded
   - Try using the locus tag instead of gene name (e.g., `Rv0001` for `dnaA`)

2. **JBrowse not displaying**:
   - Ensure the GFF3 file path is correct
   - Check browser console for JavaScript errors
   - Verify that dash-jbrowse is properly installed

3. **Slow performance**:
   - The catalogue file is large (~35MB), initial load may take time
   - Consider filtering data for specific drugs/genes

### Debug Mode

Run with debug output:
```bash
python app.py
# Check console for error messages
```

## Future Enhancements
- [ ] Export results to PDF/Excel
- [ ] Batch search for multiple genes/mutations

## License

This project is licensed under the GNU General Public License v3.0.

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## Support

For questions or issues, please open an issue on the repository.
