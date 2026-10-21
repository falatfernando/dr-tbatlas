"""
Data loading and parsing utilities for DR-TBAtlas.
Handles GFF3, catalogue master file, and genomic coordinates.

*Mycobacterium tuberculosis* is a prokaryote: it has no splicing and no
introns, so a coding sequence (CDS) is equivalent to its gene. The annotation
is therefore flattened to a single gene-level feature per locus before it
reaches the genome browser (see :mod:`browser_tracks`).
"""
import os
import re
import threading
import urllib.parse
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import gffutils
import pandas as pd

# Gene-level feature types in the H37Rv RefSeq annotation.
GENE_FEATURE_TYPES = ("gene", "pseudogene")

# WHO confidence gradings, in catalogue order.
CONFIDENCE_GRADES = (
    "1) Assoc w R",
    "2) Assoc w R - Interim",
    "3) Uncertain significance",
    "4) Not assoc w R - Interim",
    "5) Not assoc w R",
)

# Optional WHO catalogue metadata columns surfaced in the resistance table
# (TASK-12). Keys are the column ids in the master file.
RELAXED_THRESHOLDS_COLUMN = (
    "Relaxed thresholds simulation "
    "(BDQ_Rv0678, CFZ_Rv0678, INH_katG, DLM_ddn/fbiA/fbiB/fbiC/fgd1/Rv2983)"
)

# First-line drugs are listed before the rest in the Browse by Drug module.
FIRST_LINE_DRUGS = ("Isoniazid", "Rifampicin", "Pyrazinamide", "Ethambutol")

# The catalogue files some variants under a neighbouring gene and keeps their
# other name in the Comment column: fabG1 and inhA share one promoter, so
# fabG1_c.-15C>T is catalogued as inhA_c.-777C>T with "Alias fabG1_c.-15C>T".
# The trailing lookahead drops the full stop that ends the sentence without
# cutting mutations such as "p.Leu203Leu".
ALIAS_COMMENT_REGEX = re.compile(r"\bAlias\s+([A-Za-z0-9]+_\S+?)(?=[.,;]?(?:\s|$))")


@dataclass
class GeneInfo:
    """Gene information from GFF3."""
    gene_id: str
    gene_name: str
    locus_tag: str
    start: int
    end: int
    strand: str
    product: str
    chromosome: str
    note: str = ""
    biotype: str = ""
    cds_length: Optional[int] = None
    protein_length: Optional[int] = None

    @property
    def length(self) -> int:
        """Genomic span of the gene in base pairs."""
        return self.end - self.start + 1

    @property
    def display_name(self) -> str:
        """Gene symbol when available, otherwise the locus tag."""
        return self.gene_name or self.locus_tag

    @property
    def is_protein_coding(self) -> bool:
        return self.biotype == "protein_coding"


@dataclass
class CatalogueAlias:
    """A catalogue variant that the WHO also names after another gene."""
    alias: str            # e.g. "fabG1_c.-15C>T"
    alias_gene: str       # e.g. "fabG1"
    alias_mutation: str   # e.g. "c.-15C>T"
    variant: str          # catalogue identifier, e.g. "inhA_c.-777C>T"
    gene: str             # catalogue gene, e.g. "inhA"
    mutation: str         # catalogue mutation, e.g. "c.-777C>T"
    drugs: List[Tuple[str, str]]  # (drug, final confidence grading)


@dataclass
class MutationInfo:
    """Mutation information with calculated coordinates."""
    variant: str
    gene: str
    mutation: str
    genomic_position: int
    gene_relative_position: Optional[int]  # Position relative to gene start
    gene_start: int
    gene_end: int
    strand: str
    drug: Optional[str]
    resistance_profile: Optional[str]
    tier: Optional[str]
    effect: Optional[str]
    confidence: Optional[str]


def _parse_gff_attributes(raw: str) -> Dict[str, str]:
    """Parse a GFF3 attribute column into a dict with values unescaped."""
    attributes: Dict[str, str] = {}
    for chunk in raw.strip().rstrip(";").split(";"):
        if not chunk or "=" not in chunk:
            continue
        key, _, value = chunk.partition("=")
        attributes[key.strip()] = urllib.parse.unquote(value.strip())
    return attributes


def protein_length_from_cds(cds_length: Optional[int]) -> Optional[int]:
    """
    Amino acid count for a CDS, excluding the termination codon (TASK-06).

    ``L_aa = (L_bp - 3) // 3``
    """
    if not cds_length or cds_length < 6:
        return None
    return (cds_length - 3) // 3


class DataLoader:
    """Load and parse TB genomic data files."""

    def __init__(self, data_dir: str = "data"):
        self.data_dir = data_dir
        self.gff3_path = os.path.join(data_dir, "h37rv.gff3")
        self.catalogue_path = os.path.join(data_dir, "catalogue_master_file.txt")
        self.genomic_coords_path = os.path.join(data_dir, "genomic_coordinates.txt")
        self.fasta_path = os.path.join(data_dir, "h37rv.fasta")
        self.fai_path = os.path.join(data_dir, "h37rv.fasta.fai")

        self.gene_db: Optional[gffutils.FeatureDB] = None
        self.genes_cache: Dict[str, GeneInfo] = {}
        self.catalogue_df: Optional[pd.DataFrame] = None
        self.genomic_coords_df: Optional[pd.DataFrame] = None

        # Gene index, built once from the GFF3 text (see _load_gene_index).
        self._genes_by_locus: Dict[str, GeneInfo] = {}
        self._genes_by_name: Dict[str, GeneInfo] = {}
        self._genes_sorted: List[GeneInfo] = []
        self._index_lock = threading.Lock()

        # Catalogue lookups, built lazily from the master file.
        self._catalogue_alias_map: Optional[Dict[str, str]] = None
        self._variant_aliases: Optional[List[CatalogueAlias]] = None
        self._fasta = None

    # ------------------------------------------------------------------
    # GFF3 / gene annotation
    # ------------------------------------------------------------------
    def get_gene_db(self) -> gffutils.FeatureDB:
        """Get gene database, ensuring thread safety for SQLite."""
        db_path = self.gff3_path + ".db"

        # In Dash/Flask, multiple threads access this data.
        # SQLite connections can't be shared across threads.
        # We'll create a new connection if needed.
        if not hasattr(self, '_thread_local_db'):
            self._thread_local_db = threading.local()

        if not hasattr(self._thread_local_db, 'db'):
            if not os.path.exists(db_path):
                self.load_gff3()
            self._thread_local_db.db = gffutils.FeatureDB(db_path, keep_order=True)

        return self._thread_local_db.db

    def load_gff3(self) -> gffutils.FeatureDB:
        """Load GFF3 file and create gene database."""
        db_path = self.gff3_path + ".db"

        # Create database if it doesn't exist
        if not os.path.exists(db_path):
            self.gene_db = gffutils.create_db(
                self.gff3_path,
                dbfn=db_path,
                force=True,
                keep_order=True,
                merge_strategy='merge',
                sort_attribute_values=True
            )

        self._load_gene_index()

        # For the current thread, initialize the DB
        return self.get_gene_db()

    def _load_gene_index(self) -> None:
        """
        Build an in-memory index of every gene-level feature in the GFF3.

        Reading the annotation once up front replaces a full SQLite scan per
        lookup and gives every gene its product, functional note and CDS
        length, which the overview card and the genome browser tracks need.
        """
        if self._genes_sorted:
            return

        with self._index_lock:
            if self._genes_sorted:
                return

            genes: Dict[str, GeneInfo] = {}
            cds_spans: Dict[str, int] = {}
            cds_products: Dict[str, str] = {}
            cds_notes: Dict[str, str] = {}

            with open(self.gff3_path, "r", encoding="utf-8") as handle:
                for line in handle:
                    if line.startswith("#"):
                        continue
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 9:
                        continue

                    seqid, _, feature_type, start, end = fields[0], fields[1], fields[2], fields[3], fields[4]
                    strand, attributes_raw = fields[6], fields[8]

                    if feature_type not in GENE_FEATURE_TYPES and feature_type != "CDS":
                        continue

                    attributes = _parse_gff_attributes(attributes_raw)
                    locus_tag = attributes.get("locus_tag", "")
                    if not locus_tag:
                        continue

                    if feature_type == "CDS":
                        # A bacterial gene has a single CDS, but summing the
                        # spans keeps the arithmetic correct either way.
                        cds_spans[locus_tag] = cds_spans.get(locus_tag, 0) + (
                            int(end) - int(start) + 1
                        )
                        if attributes.get("product"):
                            cds_products.setdefault(locus_tag, attributes["product"])
                        if attributes.get("Note"):
                            cds_notes.setdefault(locus_tag, attributes["Note"])
                        continue

                    genes[locus_tag] = GeneInfo(
                        gene_id=attributes.get("ID", ""),
                        gene_name=attributes.get("Name") or attributes.get("gene") or locus_tag,
                        locus_tag=locus_tag,
                        start=int(start),
                        end=int(end),
                        strand=strand,
                        product="",
                        chromosome=seqid,
                        note=attributes.get("Note", ""),
                        biotype=attributes.get("gene_biotype", ""),
                    )

            for locus_tag, gene in genes.items():
                gene.product = cds_products.get(locus_tag, "")
                if not gene.note:
                    gene.note = cds_notes.get(locus_tag, "")
                gene.cds_length = cds_spans.get(locus_tag)
                gene.protein_length = protein_length_from_cds(gene.cds_length)

                self._genes_by_locus[locus_tag.lower()] = gene
                # Locus tags win over symbols on collision, so a symbol never
                # shadows a real gene ID.
                self._genes_by_name.setdefault(gene.gene_name.lower(), gene)

            self._genes_sorted = sorted(genes.values(), key=lambda g: (g.chromosome, g.start))

    def get_gene_info(self, gene_name: str) -> Optional[GeneInfo]:
        """Get gene information by name or locus tag."""
        if not gene_name:
            return None

        token = str(gene_name).strip()
        if token in self.genes_cache:
            return self.genes_cache[token]

        self._load_gene_index()
        key = token.lower()
        gene_info = self._genes_by_locus.get(key) or self._genes_by_name.get(key)

        if gene_info:
            self.genes_cache[token] = gene_info
        return gene_info

    def is_known_gene(self, token: str) -> bool:
        """True when the token is an annotated gene symbol or locus tag."""
        return self.get_gene_info(token) is not None

    def search_genes(self, query: str) -> List[GeneInfo]:
        """Search for genes by name or locus tag (partial match)."""
        if not query:
            return []

        self._load_gene_index()
        needle = str(query).strip().lower()

        exact = [
            gene for gene in self._genes_sorted
            if gene.gene_name.lower() == needle or gene.locus_tag.lower() == needle
        ]
        partial = [
            gene for gene in self._genes_sorted
            if gene not in exact
            and (needle in gene.gene_name.lower() or needle in gene.locus_tag.lower())
        ]
        return exact + partial

    def get_neighbor_genes(
        self, gene_info: GeneInfo
    ) -> Tuple[Optional[GeneInfo], Optional[GeneInfo]]:
        """The gene immediately upstream and downstream in genomic order."""
        self._load_gene_index()
        same_contig = [g for g in self._genes_sorted if g.chromosome == gene_info.chromosome]
        try:
            position = next(
                i for i, g in enumerate(same_contig) if g.locus_tag == gene_info.locus_tag
            )
        except StopIteration:
            return None, None

        previous = same_contig[position - 1] if position > 0 else None
        following = same_contig[position + 1] if position + 1 < len(same_contig) else None
        return previous, following

    # ------------------------------------------------------------------
    # Reference sequence
    # ------------------------------------------------------------------
    def get_sequence(self, start: int, end: int, chromosome: Optional[str] = None) -> str:
        """Reference sequence for a 1-based inclusive genomic interval."""
        if self._fasta is None:
            import pyfaidx

            self._fasta = pyfaidx.Fasta(self.fasta_path)

        contig = chromosome or next(iter(self._fasta.keys()))
        start = max(1, int(start))
        end = min(len(self._fasta[contig]), int(end))
        if end < start:
            return ""
        return str(self._fasta[contig][start - 1:end]).upper()

    def get_gene_sequence(
        self, gene_info: GeneInfo, flank: int = 0, reverse_complement: bool = True
    ) -> str:
        """
        Coding sequence of a gene, optionally with symmetric flanks.

        With ``reverse_complement`` the sequence is returned in the gene's own
        orientation, so a minus-strand gene reads 5'->3' from its start codon.
        """
        sequence = self.get_sequence(
            gene_info.start - flank, gene_info.end + flank, gene_info.chromosome
        )
        if reverse_complement and gene_info.strand == "-":
            sequence = _reverse_complement(sequence)
        return sequence

    def translate_gene(self, gene_info: GeneInfo) -> Optional[str]:
        """Protein translation of a protein-coding gene (bacterial code 11)."""
        if not gene_info.is_protein_coding:
            return None

        from Bio.Seq import Seq

        cds = self.get_gene_sequence(gene_info)
        usable = len(cds) - (len(cds) % 3)
        if usable < 3:
            return None
        return str(Seq(cds[:usable]).translate(table=11, to_stop=False))

    # ------------------------------------------------------------------
    # WHO catalogue
    # ------------------------------------------------------------------
    def load_catalogue(self) -> pd.DataFrame:
        """Load catalogue master file."""
        if self.catalogue_df is None:
            catalogue = pd.read_csv(self.catalogue_path, sep='\t', dtype=str)
            # The master file repeats the "CHANGES vs ver1" header; pandas
            # disambiguates the second one and it carries unrelated numeric
            # codes, so it is dropped to keep the schema unambiguous.
            catalogue = catalogue.drop(
                columns=[c for c in catalogue.columns if re.fullmatch(r"CHANGES vs ver1\.\d+", c)]
            )
            self.catalogue_df = catalogue
        return self.catalogue_df

    def load_genomic_coordinates(self) -> pd.DataFrame:
        """Load genomic coordinates file."""
        if self.genomic_coords_df is None:
            coords = pd.read_csv(self.genomic_coords_path, sep='\t', dtype=str)
            # Pre-split the gene prefix so per-gene lookups are a cheap
            # equality test rather than a string scan over 145k rows.
            coords["gene"] = coords["variant"].str.split("_", n=1).str[0]
            self.genomic_coords_df = coords
        return self.genomic_coords_df

    def _build_catalogue_aliases(self) -> Dict[str, str]:
        """
        Map every identifier of a catalogue gene to its catalogue spelling.

        Built from the annotation rather than a hand-maintained list, so all
        catalogue genes resolve from either their symbol or their locus tag.
        """
        if self._catalogue_alias_map is not None:
            return self._catalogue_alias_map

        catalogue = self.load_catalogue()
        aliases: Dict[str, str] = {}

        for catalogue_gene in catalogue["gene"].dropna().unique():
            aliases[catalogue_gene.lower()] = catalogue_gene

            gene_info = self.get_gene_info(catalogue_gene)
            if gene_info:
                aliases[gene_info.gene_name.lower()] = catalogue_gene
                aliases[gene_info.locus_tag.lower()] = catalogue_gene

        self._catalogue_alias_map = aliases
        return aliases

    def resolve_catalogue_gene(self, gene_token: str) -> Optional[str]:
        """Catalogue spelling for a gene symbol or locus tag, if catalogued."""
        if not gene_token:
            return None
        aliases = self._build_catalogue_aliases()
        return aliases.get(str(gene_token).strip().lower())

    def is_gene_in_catalogue(self, gene_info: GeneInfo) -> bool:
        """Whether a gene has any entry in the WHO catalogue (TASK-05)."""
        if gene_info is None:
            return False
        return bool(
            self.resolve_catalogue_gene(gene_info.locus_tag)
            or self.resolve_catalogue_gene(gene_info.gene_name)
        )

    def get_variant_aliases(self) -> List[CatalogueAlias]:
        """
        Catalogue variants that carry an alias in another gene, in genomic order.

        Read from the "Alias ..." notes of the Comment column, so the list
        follows the catalogue rather than a hand-maintained table.
        """
        if self._variant_aliases is not None:
            return self._variant_aliases

        catalogue = self.load_catalogue()
        if "Comment" not in catalogue.columns:
            self._variant_aliases = []
            return self._variant_aliases

        rows = catalogue.assign(
            alias=catalogue["Comment"].str.extract(ALIAS_COMMENT_REGEX, expand=False)
        ).dropna(subset=["alias", "variant"])

        found = []
        for (alias, variant), group in rows.groupby(["alias", "variant"], sort=False):
            alias_gene, alias_mutation = alias.split("_", 1)
            first = group.iloc[0]
            position = pd.to_numeric(first.get("genomic position"), errors="coerce")
            found.append((position, CatalogueAlias(
                alias=alias,
                alias_gene=alias_gene,
                alias_mutation=alias_mutation,
                variant=variant,
                gene=first["gene"],
                mutation=first["mutation"],
                drugs=list(zip(
                    group["drug"].fillna(""),
                    group["FINAL CONFIDENCE GRADING"].fillna(""),
                )),
            )))

        found.sort(key=lambda item: (pd.isna(item[0]), item[0] if pd.notna(item[0]) else 0))
        self._variant_aliases = [alias for _, alias in found]
        return self._variant_aliases

    def get_aliases_for_gene(self, gene_info: GeneInfo) -> List[CatalogueAlias]:
        """Catalogue variants named after this gene but filed under another."""
        if gene_info is None:
            return []
        aliases = []
        for alias in self.get_variant_aliases():
            alias_gene = self.get_gene_info(alias.alias_gene)
            if alias_gene is not None and alias_gene.locus_tag == gene_info.locus_tag:
                aliases.append(alias)
        return aliases

    def resolve_variant_alias(
        self, gene_token: str, mutation: Optional[str]
    ) -> Optional[CatalogueAlias]:
        """The catalogue variant a gene/mutation pair is an alias of, if any."""
        gene_info = self.get_gene_info(gene_token) if gene_token else None
        if gene_info is None or not mutation:
            return None

        def bare(value: str) -> str:
            return re.sub(r"^[pcn]\.", "", str(value).strip().lower())

        wanted = bare(mutation)
        for alias in self.get_aliases_for_gene(gene_info):
            if bare(alias.alias_mutation) == wanted:
                return alias
        return None

    def get_drug_resistance_info(
        self, gene_name: str, mutation: Optional[str] = None
    ) -> pd.DataFrame:
        """Get drug resistance information for a gene/mutation."""
        catalogue = self.load_catalogue()

        catalogue_gene = self.resolve_catalogue_gene(gene_name)
        if catalogue_gene is None:
            return catalogue.iloc[0:0].copy()

        mask = catalogue['gene'] == catalogue_gene

        if mutation:
            # mutation might be a full variant like "Rv0678_p.Asp47fs"
            # or just the mutation part like "p.Asp47fs"
            mutation_lower = str(mutation).lower()

            # 1. Try matching full variant in the 'variant' column
            variant_mask = mask & (catalogue['variant'].str.lower() == mutation_lower)

            # 2. Try matching mutation part in the 'mutation' column
            mut_part = mutation_lower
            if '_' in mutation_lower:
                mut_part = mutation_lower.split('_', 1)[1]

            mutation_column_mask = mask & (catalogue['mutation'].str.lower() == mut_part)

            # 3. Accept the mutation with or without its "p."/"c." prefix.
            bare = re.sub(r"^[pcn]\.", "", mut_part)
            prefixed_mask = mask & catalogue['mutation'].str.lower().isin(
                {bare, f"p.{bare}", f"c.{bare}", f"n.{bare}"}
            )

            final_mask = variant_mask | mutation_column_mask | prefixed_mask

            # 4. Special case for LoF
            if "lof" in mutation_lower:
                final_mask = final_mask | (mask & (catalogue['mutation'].str.lower() == "lof"))

            return catalogue[final_mask].copy()

        return catalogue[mask].copy()

    def get_drug_gene_map(self) -> Dict[str, List[Dict]]:
        """
        Drugs mapped to their catalogue genes, for the Browse by Drug module.

        Each entry carries the gene's catalogue spelling, its locus tag and the
        number of distinct catalogued mutations for that drug/gene pair.
        """
        catalogue = self.load_catalogue()
        grouped = (
            catalogue.dropna(subset=["drug", "gene"])
            .groupby(["drug", "gene"])["mutation"]
            .nunique()
            .reset_index(name="variant_count")
        )

        drug_map: Dict[str, List[Dict]] = {}
        for drug, rows in grouped.groupby("drug"):
            entries = []
            for _, row in rows.sort_values("gene", key=lambda s: s.str.lower()).iterrows():
                gene_info = self.get_gene_info(row["gene"])
                entries.append({
                    "gene": row["gene"],
                    "locus_tag": gene_info.locus_tag if gene_info else "",
                    "variant_count": int(row["variant_count"]),
                })
            drug_map[drug] = entries

        # First-line drugs first, then the remainder alphabetically.
        ordered = [d for d in FIRST_LINE_DRUGS if d in drug_map]
        ordered += sorted(d for d in drug_map if d not in FIRST_LINE_DRUGS)
        return {drug: drug_map[drug] for drug in ordered}

    def get_catalogue_summary(self) -> pd.DataFrame:
        """
        Per-gene summary of the WHO catalogue (TASK-14).

        A catalogue row is a drug/gene/mutation association, so mutations are
        counted as distinct values rather than as rows; both figures are
        reported because they answer different questions.
        """
        catalogue = self.load_catalogue()
        records = []

        for catalogue_gene, rows in catalogue.dropna(subset=["gene"]).groupby("gene"):
            gene_info = self.get_gene_info(catalogue_gene)
            mutations = rows["mutation"].dropna()
            grades = rows["FINAL CONFIDENCE GRADING"].fillna("")

            record = {
                "Gene": catalogue_gene,
                "Locus tag": gene_info.locus_tag if gene_info else "",
                "Drugs": ", ".join(sorted(rows["drug"].dropna().unique())),
                "Drug count": int(rows["drug"].nunique()),
                "Variants": int(mutations.nunique()),
                "Associations": int(len(rows)),
                "Tier 1": int((rows["tier"] == "1").sum()),
                "Tier 2": int((rows["tier"] == "2").sum()),
                "LoF": int(
                    mutations[
                        mutations.str.contains("lof", case=False, na=False)
                    ].nunique()
                ),
            }
            for grade in CONFIDENCE_GRADES:
                record[grade] = int((grades == grade).sum())
            records.append(record)

        summary = pd.DataFrame(records)
        if summary.empty:
            return summary
        return summary.sort_values("Variants", ascending=False).reset_index(drop=True)

    def get_drug_summary(self) -> pd.DataFrame:
        """Per-drug summary of the WHO catalogue (TASK-14)."""
        catalogue = self.load_catalogue()
        records = []

        for drug, rows in catalogue.dropna(subset=["drug"]).groupby("drug"):
            grades = rows["FINAL CONFIDENCE GRADING"].fillna("")
            record = {
                "Drug": drug,
                "Genes": int(rows["gene"].nunique()),
                "Variants": int(rows["variant"].nunique()),
                "Associations": int(len(rows)),
                "Tier 1": int((rows["tier"] == "1").sum()),
                "Tier 2": int((rows["tier"] == "2").sum()),
            }
            for grade in CONFIDENCE_GRADES:
                record[grade] = int((grades == grade).sum())
            records.append(record)

        summary = pd.DataFrame(records)
        if summary.empty:
            return summary

        order = {drug: i for i, drug in enumerate(FIRST_LINE_DRUGS)}
        summary["_order"] = summary["Drug"].map(lambda d: order.get(d, len(order)))
        return (
            summary.sort_values(["_order", "Drug"])
            .drop(columns="_order")
            .reset_index(drop=True)
        )

    def get_catalogue_totals(self) -> Dict[str, int]:
        """Headline catalogue counts for the summary view."""
        catalogue = self.load_catalogue()
        grades = catalogue["FINAL CONFIDENCE GRADING"].fillna("")
        totals = {
            "genes": int(catalogue["gene"].nunique()),
            "drugs": int(catalogue["drug"].nunique()),
            "variants": int(catalogue["variant"].nunique()),
            "associations": int(len(catalogue)),
            "tier_1": int((catalogue["tier"] == "1").sum()),
            "tier_2": int((catalogue["tier"] == "2").sum()),
        }
        for grade in CONFIDENCE_GRADES:
            totals[grade] = int((grades == grade).sum())
        return totals

    # ------------------------------------------------------------------
    # Coordinates
    # ------------------------------------------------------------------
    def calculate_relative_position(self, genomic_position: int, gene_info: GeneInfo) -> int:
        """
        Calculate position relative to gene start.
        For + strand: relative_pos = genomic_pos - gene_start
        For - strand: relative_pos = gene_end - genomic_pos
        """
        if gene_info.strand == '+':
            return genomic_position - gene_info.start
        else:
            return gene_info.end - genomic_position

    def calculate_genomic_position(self, relative_position: int, gene_info: GeneInfo) -> int:
        """
        Calculate genomic position from gene-relative position.
        For + strand: genomic_pos = gene_start + relative_pos
        For - strand: genomic_pos = gene_end - relative_pos
        """
        if gene_info.strand == '+':
            return gene_info.start + relative_position
        else:
            return gene_info.end - relative_position

    def search_mutations_by_gene(self, gene_name: str) -> pd.DataFrame:
        """Search for all mutations related to a gene."""
        genomic_coords = self.load_genomic_coordinates()

        # Variants are prefixed with the catalogue spelling of the gene, so a
        # locus-tag search has to be translated first.
        catalogue_gene = self.resolve_catalogue_gene(gene_name) or gene_name
        return genomic_coords[genomic_coords["gene"] == catalogue_gene].copy()

    def get_variant_coordinates(self, variant: str) -> pd.DataFrame:
        """Every nucleotide change recorded for one variant identifier."""
        genomic_coords = self.load_genomic_coordinates()
        return genomic_coords[genomic_coords["variant"] == variant].copy()

    def get_mutation_details(self, variant: str) -> Dict:
        """Get detailed information about a specific mutation."""
        genomic_coords = self.load_genomic_coordinates()
        catalogue = self.load_catalogue()

        # Find in genomic coordinates
        gc_result = genomic_coords[genomic_coords['variant'] == variant]

        # Find in catalogue
        cat_result = catalogue[catalogue['variant'] == variant]

        return {
            'genomic_coordinates': gc_result.to_dict('records') if len(gc_result) > 0 else [],
            'catalogue': cat_result.to_dict('records') if len(cat_result) > 0 else []
        }


_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]
