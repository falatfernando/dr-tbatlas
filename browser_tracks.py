"""
Track data and configuration for the embedded JBrowse 2 genome browser.

The browser is given, besides the reference sequence and the flattened gene
annotation, a set of tracks derived from the WHO catalogue and the reference
itself, in the spirit of the Mycobrowser JBrowse instance:

* the reference sequence, with six-frame translation when zoomed in;
* genes coloured by strand and biotype, with product, functional note,
  catalogue drugs and a Mycobrowser link in the feature details;
* the catalogue genes, coloured by tier;
* catalogue variants split by confidence grading, and one
  resistance-associated variant track per drug;
* catalogue variant density and GC content, for the whole-genome overview;
* a "current selection" track marking the active gene and mutation;
* a text-search index, so the browser's location box accepts gene names,
  locus tags, products and variants.

Every file is derived from the tracked data at start-up and rebuilt only when
its sources change, like the flattened annotation.
"""
import gzip
import json
import os
import re
import urllib.parse
from collections import defaultdict
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd

from data_utils import FIRST_LINE_DRUGS, DataLoader, GeneInfo

ASSEMBLY_NAME = "H37Rv NC_000962.3"
TRACKS_DIRNAME = "tracks"
# Bump when the content of the generated files changes, so existing
# installations rebuild them on the next start.
TRACKS_FORMAT_VERSION = 2
DATA_URL = "/data"

MYCOBROWSER_GENE_URL = "https://mycobrowser.epfl.ch/genes/{locus_tag}"

# Bin width of the catalogue variant density track, and GC content window.
DENSITY_BIN_BP = 1000
GC_WINDOW_BP = 100

# Trix index layout: prefix length and target block size of the .ixx.
_TRIX_PREFIX = 5
_TRIX_BLOCK_BYTES = 16 * 1024

# Longest allele written out in full in a feature's details.
_MAX_ALLELE_CHARS = 30

# Grade group (the leading digit of the WHO grading) -> colour. Resistance
# reads red, uncertain grey and not associated green.
GRADE_COLORS = {
    "1": "#B42318",
    "2": "#DC6803",
    "3": "#98A2B3",
    "4": "#6FAE8F",
    "5": "#2E7D5B",
}
GRADE_LABELS = {
    "1": "Associated with resistance",
    "2": "Associated with resistance (interim)",
    "3": "Uncertain significance",
    "4": "Not associated (interim)",
    "5": "Not associated",
}

GENE_COLORS = {
    "forward": "#438E8E",
    "reverse": "#5B7DB1",
    "noncoding": "#C98B3F",
    "pseudogene": "#A0AEC0",
}
TIER_COLORS = {"1": "#7A1F5C", "2": "#B56FA0"}
SELECTION_COLOR = "#F2B705"
GC_COLOR = "#3A7D7D"
DENSITY_COLOR = "#B42318"

# Track ids.
SEQUENCE_TRACK = "h37rv-seq"
GENES_TRACK = "genes"
CATALOGUE_GENES_TRACK = "who-genes"
SELECTION_TRACK = "current-selection"
GC_TRACK = "gc-content"
DENSITY_TRACK = "who-variant-density"
VARIANT_TRACKS = {
    "assoc": "who-variants-assoc",
    "uncertain": "who-variants-uncertain",
    "not_assoc": "who-variants-not-assoc",
}
_VARIANT_GROUPS = {
    "assoc": ("1", "2"),
    "uncertain": ("3",),
    "not_assoc": ("4", "5"),
}
_GROUP_TRACK = {
    grade: VARIANT_TRACKS[key] for key, grades in _VARIANT_GROUPS.items() for grade in grades
}

# Tracks the user can toggle from the dashboard, in display order.
TRACK_CHOICES = [
    {"value": SELECTION_TRACK, "label": "Current selection"},
    {"value": GENES_TRACK, "label": "Genes"},
    {"value": CATALOGUE_GENES_TRACK, "label": "WHO catalogue genes"},
    {"value": VARIANT_TRACKS["assoc"], "label": "Variants: associated with R"},
    {"value": VARIANT_TRACKS["uncertain"], "label": "Variants: uncertain"},
    {"value": VARIANT_TRACKS["not_assoc"], "label": "Variants: not associated"},
    {"value": DENSITY_TRACK, "label": "Variant density"},
    {"value": GC_TRACK, "label": "GC content"},
]
DEFAULT_TRACKS = [
    SELECTION_TRACK,
    GENES_TRACK,
    VARIANT_TRACKS["assoc"],
    GC_TRACK,
]

THEME = {
    "palette": {
        "primary": {"main": "#438E8E"},
        "secondary": {"main": "#3A7D7D"},
        "tertiary": {"main": "#9CCBCB"},
        "quaternary": {"main": "#2D3748"},
    },
}


def drug_track_id(drug: str) -> str:
    return "who-drug-" + re.sub(r"[^a-z0-9]+", "-", drug.lower()).strip("-")


def grade_group(grading: Optional[str]) -> Optional[str]:
    """The leading WHO group number of a grading, e.g. "1) Assoc w R" -> "1"."""
    match = re.match(r"\s*([1-5])\)", str(grading or ""))
    return match.group(1) if match else None


def _uri_component(value: str) -> str:
    """Python equivalent of JavaScript's ``encodeURIComponent``."""
    return urllib.parse.quote(str(value), safe="-_.!~*'()")


def _is_fresh(target: str, sources: Iterable[str]) -> bool:
    if not os.path.exists(target):
        return False
    built = os.path.getmtime(target)
    return all(os.path.getmtime(source) <= built for source in sources)


def _open_track(path: str):
    """Text handle for a gzipped track file; JBrowse inflates it on load."""
    return gzip.open(path, "wt", encoding="utf-8", newline="\n")


def _gff_escape(value: str) -> str:
    """
    Escape a value for use inside a GFF3 attribute field.

    ``,``, ``;``, ``=``, ``&`` and ``%`` carry structural meaning in the GFF3
    attribute column, so they must stay percent-encoded — a product name such
    as "KatG,catalase-peroxidase" would otherwise be read as two values.
    """
    return urllib.parse.quote(str(value), safe=" ()[]{}<>/'+*:.-_")


def _gff_line(seqid, source, feature_type, start, end, strand, attributes) -> str:
    fields = []
    for key, value in attributes:
        if value is None or value == "" or value == []:
            continue
        if isinstance(value, (list, tuple)):
            value = ",".join(_gff_escape(v) for v in value)
        else:
            value = _gff_escape(value)
        fields.append(f"{key}={value}")
    return "\t".join([
        seqid, source, feature_type, str(start), str(end), ".", strand or ".", ".",
        ";".join(fields),
    ]) + "\n"


def _shorten_allele(allele: str) -> str:
    allele = str(allele)
    if len(allele) <= _MAX_ALLELE_CHARS:
        return allele
    return f"{allele[:_MAX_ALLELE_CHARS]}… [{len(allele):,} bp]"


def _mycobrowser_link(locus_tag: str) -> str:
    url = MYCOBROWSER_GENE_URL.format(locus_tag=urllib.parse.quote(locus_tag))
    return f'<a href="{url}" target="_blank" rel="noopener noreferrer">{locus_tag}</a>'


class BrowserTracks:
    """Builds the derived track files and the JBrowse 2 configuration."""

    def __init__(self, data_loader: DataLoader, tracks_dir: Optional[str] = None):
        self.data_loader = data_loader
        # Served under /data/tracks, so it only moves for tests.
        self.tracks_dir = tracks_dir or os.path.join(data_loader.data_dir, TRACKS_DIRNAME)
        self._variants: Optional[pd.DataFrame] = None
        self._drugs: List[str] = []

    # ------------------------------------------------------------------
    # Paths
    # ------------------------------------------------------------------
    def _path(self, name: str) -> str:
        return os.path.join(self.tracks_dir, name)

    @staticmethod
    def _url(name: str) -> str:
        return f"{DATA_URL}/{TRACKS_DIRNAME}/{name}"

    @property
    def _catalogue_sources(self) -> List[str]:
        loader = self.data_loader
        return [loader.catalogue_path, loader.genomic_coords_path, loader.gff3_path]

    # ------------------------------------------------------------------
    # Build
    # ------------------------------------------------------------------
    def build(self) -> "BrowserTracks":
        """Write every derived track file whose sources or format changed."""
        os.makedirs(self.tracks_dir, exist_ok=True)

        self._drugs = self._ordered_drugs()
        sources = self._catalogue_sources
        stamp = self._path("build.json")

        if (
            not _is_fresh(stamp, sources)
            or self._stamp_version(stamp) != TRACKS_FORMAT_VERSION
            or not self._all_outputs_exist()
        ):
            variants = self.variant_features()
            self._write_genes(self._path("genes.gff3.gz"))
            self._write_catalogue_genes(self._path("who_genes.gff3.gz"))
            for group, track_id in VARIANT_TRACKS.items():
                subset = variants[variants["grade_group"].isin(_VARIANT_GROUPS[group])]
                self._write_variants(self._path(f"{track_id}.gff3.gz"), subset)
            for drug in self._drugs:
                self._write_drug_variants(self._path(f"{drug_track_id(drug)}.gff3.gz"), variants, drug)
            self._write_density(self._path("who_variant_density.bed"), variants)
            self._write_text_index(variants)
            with open(stamp, "w", encoding="utf-8") as handle:
                json.dump({"version": TRACKS_FORMAT_VERSION, "drugs": self._drugs}, handle)
        return self

    @staticmethod
    def _stamp_version(stamp: str) -> Optional[int]:
        try:
            with open(stamp, encoding="utf-8") as handle:
                return json.load(handle).get("version")
        except (OSError, ValueError):
            return None

    def _output_names(self) -> List[str]:
        names = ["genes.gff3.gz", "who_genes.gff3.gz", "who_variant_density.bed",
                 "search.ix", "search.ixx", "search_meta.json"]
        names += [f"{track_id}.gff3.gz" for track_id in VARIANT_TRACKS.values()]
        names += [f"{drug_track_id(drug)}.gff3.gz" for drug in self._drugs]
        return names

    def _all_outputs_exist(self) -> bool:
        return all(os.path.exists(self._path(name)) for name in self._output_names())

    def _ordered_drugs(self) -> List[str]:
        drugs = set(self.data_loader.load_catalogue()["drug"].dropna().unique())
        ordered = [d for d in FIRST_LINE_DRUGS if d in drugs]
        return ordered + sorted(drugs - set(ordered))

    @property
    def drugs(self) -> List[str]:
        if not self._drugs:
            self._drugs = self._ordered_drugs()
        return self._drugs

    # ------------------------------------------------------------------
    # Catalogue variants
    # ------------------------------------------------------------------
    def variant_features(self) -> pd.DataFrame:
        """
        One row per graded catalogue variant that has genomic coordinates.

        A variant can be recorded as several nucleotide changes (alternative
        codons, or equivalent indel representations), so its span covers all
        of them. Its colour follows the strongest grading across drugs.
        """
        if self._variants is not None:
            return self._variants

        loader = self.data_loader
        catalogue = loader.load_catalogue()
        coordinates = loader.load_genomic_coordinates()

        graded = catalogue.dropna(subset=["variant", "drug"]).copy()
        graded["grade_group"] = graded["FINAL CONFIDENCE GRADING"].map(grade_group)
        graded = graded.dropna(subset=["grade_group"])

        coords = coordinates[coordinates["variant"].isin(graded["variant"])].copy()
        coords["position"] = coords["position"].astype(int)
        ref = coords["reference_nucleotide"].fillna("")
        coords["end"] = coords["position"] + ref.str.len().clip(lower=1) - 1
        coords["change"] = (
            coords["position"].astype(str) + " "
            + ref.map(_shorten_allele) + ">"
            + coords["alternative_nucleotide"].fillna("").map(_shorten_allele)
        )

        spans = coords.groupby("variant").agg(
            chromosome=("chromosome", "first"),
            start=("position", "min"),
            end=("end", "max"),
            changes=("change", lambda c: list(dict.fromkeys(c))),
        )

        per_variant = graded.groupby("variant").agg(
            gene=("gene", "first"),
            mutation=("mutation", "first"),
            tier=("tier", lambda t: sorted({str(v) for v in t.dropna()})),
            effect=("effect", "first"),
            grade_group=("grade_group", "min"),
        )
        per_drug: Dict[str, List[Tuple[str, str, str]]] = defaultdict(list)
        for variant, drug, grading, group in graded[
            ["variant", "drug", "FINAL CONFIDENCE GRADING", "grade_group"]
        ].itertuples(index=False):
            per_drug[variant].append((drug, grading, group))

        variants = spans.join(per_variant, how="inner").reset_index()
        variants["drugs"] = variants["variant"].map(per_drug)
        # Names the catalogue also gives a variant, e.g. fabG1_c.-15C>T (#5).
        aliases: Dict[str, List[str]] = defaultdict(list)
        for alias in loader.get_variant_aliases():
            aliases[alias.variant].append(alias.alias)
        variants["aliases"] = variants["variant"].map(lambda v: aliases.get(v, []))
        variants["strand"] = variants["gene"].map(self._gene_strand)
        self._variants = variants.sort_values(["chromosome", "start"]).reset_index(drop=True)
        return self._variants

    def _gene_strand(self, gene: str) -> str:
        gene_info = self.data_loader.get_gene_info(gene)
        return gene_info.strand if gene_info else "."

    def _variant_attributes(self, row, drug_filter: Optional[str] = None) -> List[Tuple]:
        drugs = row.drugs
        if drug_filter:
            drugs = [entry for entry in drugs if entry[0] == drug_filter]
        group = min(entry[2] for entry in drugs)
        gene_info = self.data_loader.get_gene_info(row.gene)
        return [
            ("ID", row.variant if not drug_filter else f"{row.variant}:{drug_filter}"),
            # A bare "LoF" or "deletion" label says nothing out of context.
            ("Name", row.mutation if row.mutation.startswith(("p.", "c.", "n.")) else row.variant),
            ("variant", row.variant),
            ("alias", row.aliases or None),
            ("gene", row.gene),
            ("locus_tag", gene_info.locus_tag if gene_info else None),
            ("grade_group", group),
            ("who_grading", GRADE_LABELS[group]),
            ("drug_gradings", [f"{drug}: {grading}" for drug, grading, _ in drugs]),
            ("tier", row.tier),
            ("effect", row.effect.replace("_", " ") if isinstance(row.effect, str) else None),
            ("nucleotide_changes", row.changes),
        ]

    def _write_variants(self, path: str, variants: pd.DataFrame) -> None:
        with _open_track(path) as out:
            out.write("##gff-version 3\n")
            for row in variants.itertuples(index=False):
                out.write(_gff_line(
                    row.chromosome, "WHO_catalogue_v2", "sequence_alteration",
                    row.start, row.end, row.strand, self._variant_attributes(row),
                ))

    def _write_drug_variants(self, path: str, variants: pd.DataFrame, drug: str) -> None:
        """Variants graded as resistance-associated (groups 1-2) for one drug."""
        with _open_track(path) as out:
            out.write("##gff-version 3\n")
            for row in variants.itertuples(index=False):
                if not any(d == drug and g in ("1", "2") for d, _, g in row.drugs):
                    continue
                out.write(_gff_line(
                    row.chromosome, "WHO_catalogue_v2", "sequence_alteration",
                    row.start, row.end, row.strand, self._variant_attributes(row, drug),
                ))

    def _write_density(self, path: str, variants: pd.DataFrame) -> None:
        """Graded catalogue variants per fixed-width bin, as scored BED."""
        counts: Dict[Tuple[str, int], int] = defaultdict(int)
        for chromosome, start in variants[["chromosome", "start"]].itertuples(index=False):
            counts[(chromosome, (start - 1) // DENSITY_BIN_BP)] += 1
        with open(path, "w", encoding="utf-8") as out:
            for (chromosome, bin_index), count in sorted(counts.items()):
                start = bin_index * DENSITY_BIN_BP
                out.write(
                    f"{chromosome}\t{start}\t{start + DENSITY_BIN_BP}\t"
                    f"{count}_variants\t{count}\n"
                )

    # ------------------------------------------------------------------
    # Genes
    # ------------------------------------------------------------------
    def _catalogue_gene_drugs(self) -> Dict[str, Dict]:
        """Drugs and tiers of every catalogue gene, keyed by locus tag."""
        catalogue = self.data_loader.load_catalogue().dropna(subset=["gene", "drug"])
        summary: Dict[str, Dict] = {}
        for gene, rows in catalogue.groupby("gene"):
            gene_info = self.data_loader.get_gene_info(gene)
            if gene_info is None:
                continue
            tiers = sorted({str(t) for t in rows["tier"].dropna()})
            summary[gene_info.locus_tag] = {
                "drugs": sorted(rows["drug"].unique()),
                "tier": tiers[0] if tiers else None,
                "tiers": tiers,
                "variants": int(rows["variant"].nunique()),
            }
        return summary

    def _gene_attributes(self, gene: GeneInfo, catalogue: Optional[Dict]) -> List[Tuple]:
        attributes = [
            ("ID", gene.locus_tag),
            ("Name", gene.display_name),
            ("locus_tag", gene.locus_tag),
            ("product", gene.product),
            ("Note", gene.note),
            ("gene_biotype", gene.biotype),
            ("length_bp", f"{gene.length:,}"),
            ("protein_length_aa", f"{gene.protein_length:,}" if gene.protein_length else None),
            ("mycobrowser", _mycobrowser_link(gene.locus_tag)),
        ]
        if catalogue:
            attributes += [
                ("who_catalogue_drugs", catalogue["drugs"]),
                ("who_tier", catalogue["tier"]),
                ("who_catalogue_variants", str(catalogue["variants"])),
            ]
        return attributes

    def _write_genes(self, path: str) -> None:
        catalogue = self._catalogue_gene_drugs()
        with _open_track(path) as out:
            out.write("##gff-version 3\n")
            for gene in self.data_loader._genes_sorted:
                feature_type = "pseudogene" if gene.biotype == "pseudogene" else "gene"
                out.write(_gff_line(
                    gene.chromosome, "RefSeq", feature_type, gene.start, gene.end,
                    gene.strand, self._gene_attributes(gene, catalogue.get(gene.locus_tag)),
                ))

    def _write_catalogue_genes(self, path: str) -> None:
        catalogue = self._catalogue_gene_drugs()
        with _open_track(path) as out:
            out.write("##gff-version 3\n")
            for gene in self.data_loader._genes_sorted:
                entry = catalogue.get(gene.locus_tag)
                if entry is None:
                    continue
                out.write(_gff_line(
                    gene.chromosome, "WHO_catalogue_v2", "gene", gene.start, gene.end,
                    gene.strand, self._gene_attributes(gene, entry),
                ))

    # ------------------------------------------------------------------
    # Text search (Trix index)
    # ------------------------------------------------------------------
    def _write_text_index(self, variants: pd.DataFrame) -> None:
        """
        Write a Trix index (``.ix`` + ``.ixx``) for the browser's search box.

        Each ``.ix`` line is a lower-case search term followed by the records
        it matches. A record is the JSON array ``[locString, trackId, label,
        ...]`` with every element URI-encoded and commas replaced by ``|``,
        which is the layout JBrowse's ``TrixTextSearchAdapter`` decodes.
        """
        terms: Dict[str, List[str]] = defaultdict(list)

        def record(loc: str, track_id: str, *labels: str) -> str:
            values = [_uri_component(v) for v in (loc, track_id, *labels) if v]
            return json.dumps(values, separators=(",", ":")).replace(",", "|")

        def add(term: str, doc: str) -> None:
            term = term.lower().strip()
            if term and doc not in terms[term]:
                terms[term].append(doc)

        for gene in self.data_loader._genes_sorted:
            flank = max(100, gene.length // 10)
            loc = f"{gene.chromosome}:{max(1, gene.start - flank)}..{gene.end + flank}"
            doc = record(loc, GENES_TRACK, gene.display_name, gene.locus_tag, gene.product)
            add(gene.display_name, doc)
            add(gene.locus_tag, doc)
            for word in re.findall(r"[A-Za-z0-9][A-Za-z0-9\-]{3,}", gene.product or ""):
                add(word, doc)

        for row in variants.itertuples(index=False):
            loc = f"{row.chromosome}:{max(1, row.start - 30)}..{row.end + 30}"
            doc = record(loc, _GROUP_TRACK[row.grade_group], row.variant, GRADE_LABELS[row.grade_group])
            add(row.variant, doc)
            for alias in row.aliases:
                add(alias, doc)

        index_lines = [f"{term} {' '.join(docs)}\n" for term, docs in sorted(terms.items())]

        # The .ixx maps a term prefix to the byte offset of its block. As in
        # UCSC ixIxx, a block only starts where the prefix changes: a lookup
        # starts at the last block whose prefix sorts before the query, so a
        # block starting midway through a prefix would skip earlier matches.
        ixx_lines = []
        offset, block_start, previous_prefix = 0, None, None
        for line in index_lines:
            prefix = line.split(" ", 1)[0][:_TRIX_PREFIX]
            if prefix != previous_prefix and (
                block_start is None or offset - block_start >= _TRIX_BLOCK_BYTES
            ):
                ixx_lines.append(f"{prefix}{offset:010X}\n")
                block_start = offset
            previous_prefix = prefix
            offset += len(line.encode("utf-8"))

        with open(self._path("search.ix"), "w", encoding="utf-8", newline="\n") as out:
            out.writelines(index_lines)
        with open(self._path("search.ixx"), "w", encoding="utf-8", newline="\n") as out:
            out.writelines(ixx_lines)
        with open(self._path("search_meta.json"), "w", encoding="utf-8") as out:
            json.dump({
                "tracks": [
                    {"trackId": GENES_TRACK, "attributesIndexed": ["Name", "locus_tag", "product"]},
                    {"trackId": VARIANT_TRACKS["assoc"], "attributesIndexed": ["variant"]},
                ],
            }, out)

    # ------------------------------------------------------------------
    # JBrowse configuration
    # ------------------------------------------------------------------
    def assembly(self) -> Dict:
        return {
            "name": ASSEMBLY_NAME,
            "aliases": ["H37Rv", "NC_000962.3"],
            "sequence": {
                "type": "ReferenceSequenceTrack",
                "trackId": SEQUENCE_TRACK,
                "adapter": self._fasta_adapter(),
            },
        }

    @staticmethod
    def _fasta_adapter() -> Dict:
        return {
            "type": "IndexedFastaAdapter",
            "fastaLocation": {"uri": f"{DATA_URL}/h37rv.fasta", "locationType": "UriLocation"},
            "faiLocation": {"uri": f"{DATA_URL}/h37rv.fasta.fai", "locationType": "UriLocation"},
        }

    def _gff_track(
        self, track_id: str, name: str, file_name: str, category: List[str],
        color: str, description: Optional[str] = None, display: Optional[Dict] = None,
    ) -> Dict:
        renderer = {"type": "SvgFeatureRenderer", "color1": color}
        if description:
            renderer["labels"] = {"description": description}
        return {
            "type": "FeatureTrack",
            "trackId": track_id,
            "name": name,
            "category": category,
            "assemblyNames": [ASSEMBLY_NAME],
            "adapter": {
                "type": "Gff3Adapter",
                "gffLocation": {"uri": self._url(file_name), "locationType": "UriLocation"},
            },
            "displays": [{
                "type": "LinearBasicDisplay",
                "displayId": f"{track_id}-LinearBasicDisplay",
                "renderer": renderer,
                **(display or {}),
            }],
        }

    @staticmethod
    def _grade_color_expression(attribute: str = "grade_group") -> str:
        expression = f"'{GRADE_COLORS['3']}'"
        for group in ("5", "4", "2", "1"):
            expression = (
                f"get(feature,'{attribute}')=='{group}'?'{GRADE_COLORS[group]}':({expression})"
            )
        return "jexl:" + expression

    def tracks(self, selection_features: Optional[List[Dict]] = None) -> List[Dict]:
        """Every track offered by the browser, including the track selector."""
        gene_color = (
            "jexl:get(feature,'type')=='pseudogene'?'{pseudo}':"
            "(get(feature,'gene_biotype')!='protein_coding'?'{nc}':"
            "(get(feature,'strand')>0?'{fwd}':'{rev}'))"
        ).format(
            pseudo=GENE_COLORS["pseudogene"], nc=GENE_COLORS["noncoding"],
            fwd=GENE_COLORS["forward"], rev=GENE_COLORS["reverse"],
        )
        tier_color = (
            f"jexl:get(feature,'who_tier')=='1'?'{TIER_COLORS['1']}':'{TIER_COLORS['2']}'"
        )
        grade_color = self._grade_color_expression()
        # Hotspots such as rpoB's RRDR stack dozens of variants deep.
        variant_display = {"maxFeatureScreenDensity": 6, "height": 180}

        tracks = [
            {
                "type": "FeatureTrack",
                "trackId": SELECTION_TRACK,
                "name": "Current selection",
                "category": ["DR-TBAtlas"],
                "assemblyNames": [ASSEMBLY_NAME],
                "adapter": {"type": "FromConfigAdapter", "features": selection_features or []},
                "displays": [{
                    "type": "LinearBasicDisplay",
                    "displayId": f"{SELECTION_TRACK}-LinearBasicDisplay",
                    "height": 90,
                    "renderer": {
                        "type": "SvgFeatureRenderer",
                        "color1": "jexl:get(feature,'color')",
                        "labels": {"description": "jexl:get(feature,'description')"},
                    },
                }],
            },
            self._gff_track(
                GENES_TRACK, "Genes (H37Rv RefSeq)", "genes.gff3.gz", ["Annotation"],
                gene_color, "jexl:get(feature,'product')", {"height": 130},
            ),
            self._gff_track(
                CATALOGUE_GENES_TRACK, "WHO catalogue genes (by tier)", "who_genes.gff3.gz",
                ["WHO catalogue"], tier_color,
                "jexl:'Tier '+get(feature,'who_tier')+' · '+get(feature,'who_catalogue_drugs')",
                {"height": 70},
            ),
            self._gff_track(
                VARIANT_TRACKS["assoc"], "Variants: associated with resistance (groups 1-2)",
                f"{VARIANT_TRACKS['assoc']}.gff3.gz", ["WHO catalogue", "Variants by grading"],
                grade_color, display=variant_display,
            ),
            self._gff_track(
                VARIANT_TRACKS["uncertain"], "Variants: uncertain significance (group 3)",
                f"{VARIANT_TRACKS['uncertain']}.gff3.gz", ["WHO catalogue", "Variants by grading"],
                grade_color, display=variant_display,
            ),
            self._gff_track(
                VARIANT_TRACKS["not_assoc"], "Variants: not associated with resistance (groups 4-5)",
                f"{VARIANT_TRACKS['not_assoc']}.gff3.gz", ["WHO catalogue", "Variants by grading"],
                grade_color, display=variant_display,
            ),
        ]

        for drug in self.drugs:
            track_id = drug_track_id(drug)
            tracks.append(self._gff_track(
                track_id, f"{drug}: resistance-associated variants",
                f"{track_id}.gff3.gz", ["WHO catalogue", "Resistance variants by drug"],
                grade_color, display=variant_display,
            ))

        tracks += [
            {
                "type": "QuantitativeTrack",
                "trackId": DENSITY_TRACK,
                "name": f"Catalogue variant density (per {DENSITY_BIN_BP // 1000} kb)",
                "category": ["WHO catalogue"],
                "assemblyNames": [ASSEMBLY_NAME],
                "adapter": {
                    "type": "BedAdapter",
                    "bedLocation": {
                        "uri": self._url("who_variant_density.bed"),
                        "locationType": "UriLocation",
                    },
                },
                "displays": [{
                    "type": "LinearWiggleDisplay",
                    "displayId": f"{DENSITY_TRACK}-LinearWiggleDisplay",
                    "height": 40,
                    "defaultRendering": "density",
                    "renderers": {"DensityRenderer": {"type": "DensityRenderer", "posColor": DENSITY_COLOR}},
                }],
            },
            {
                # A quantitative track over the GC content adapter: a second
                # ReferenceSequenceTrack cannot carry its own assembly names.
                "type": "QuantitativeTrack",
                "trackId": GC_TRACK,
                "name": f"GC content ({GC_WINDOW_BP} bp windows)",
                "category": ["Reference"],
                "assemblyNames": [ASSEMBLY_NAME],
                "adapter": {
                    "type": "GCContentAdapter",
                    "sequenceAdapter": self._fasta_adapter(),
                    "windowSize": GC_WINDOW_BP,
                    "windowDelta": GC_WINDOW_BP,
                },
                "displays": [{
                    "type": "LinearWiggleDisplay",
                    "displayId": f"{GC_TRACK}-LinearWiggleDisplay",
                    "height": 70,
                    "autoscale": "local",
                    "renderers": {"XYPlotRenderer": {"type": "XYPlotRenderer", "color": GC_COLOR}},
                }],
            },
        ]
        return tracks

    def text_search_adapters(self) -> List[Dict]:
        return [{
            "type": "TrixTextSearchAdapter",
            "textSearchAdapterId": "dr-tbatlas-search",
            "ixFilePath": {"uri": self._url("search.ix"), "locationType": "UriLocation"},
            "ixxFilePath": {"uri": self._url("search.ixx"), "locationType": "UriLocation"},
            "metaFilePath": {"uri": self._url("search_meta.json"), "locationType": "UriLocation"},
            "assemblyNames": [ASSEMBLY_NAME],
        }]

    def default_session(self, visible_tracks: List[str], tracks: List[Dict]) -> Dict:
        """A session that opens with the requested tracks, in track order."""
        by_id = {track["trackId"]: track for track in tracks}
        session_tracks = []
        # The reference sequence always sits under the selection and genes.
        ordered = [SELECTION_TRACK, GENES_TRACK, CATALOGUE_GENES_TRACK]
        ordered += [t for t in visible_tracks if t not in ordered]
        for track_id in ordered:
            if track_id not in visible_tracks or track_id not in by_id:
                continue
            track = by_id[track_id]
            session_tracks.append({
                "type": track["type"],
                "configuration": track_id,
                "displays": [{
                    "type": track["displays"][0]["type"],
                    "configuration": track["displays"][0]["displayId"],
                }],
            })
            if track_id == CATALOGUE_GENES_TRACK or (
                track_id == GENES_TRACK and CATALOGUE_GENES_TRACK not in visible_tracks
            ):
                session_tracks.append(self._sequence_session_track())
        if not any(t["configuration"] == SEQUENCE_TRACK for t in session_tracks):
            session_tracks.insert(0, self._sequence_session_track())

        return {
            "name": "DR-TBAtlas session",
            "view": {
                "id": "linear-genome-view",
                "type": "LinearGenomeView",
                "trackLabels": "offset",
                "tracks": session_tracks,
            },
        }

    @staticmethod
    def _sequence_session_track() -> Dict:
        return {
            "type": "ReferenceSequenceTrack",
            "configuration": SEQUENCE_TRACK,
            "displays": [{
                "type": "LinearReferenceSequenceDisplay",
                "configuration": f"{SEQUENCE_TRACK}-LinearReferenceSequenceDisplay",
            }],
        }

    def config(
        self,
        visible_tracks: Optional[List[str]] = None,
        selection_features: Optional[List[Dict]] = None,
    ) -> Dict:
        """Props for ``dash_jbrowse.LinearGenomeView``."""
        visible = list(visible_tracks if visible_tracks is not None else DEFAULT_TRACKS)
        tracks = self.tracks(selection_features)
        return {
            "assembly": self.assembly(),
            "tracks": tracks,
            "defaultSession": self.default_session(visible, tracks),
            "aggregateTextSearchAdapters": self.text_search_adapters(),
            "configuration": {"theme": THEME},
        }


def selection_features(
    gene: Optional[GeneInfo], selection: Optional[Dict] = None
) -> List[Dict]:
    """Features for the "Current selection" track: active gene and variant."""
    features: List[Dict] = []
    if gene is not None:
        features.append({
            "uniqueId": f"selected-gene-{gene.locus_tag}",
            "refName": gene.chromosome,
            "start": gene.start - 1,
            "end": gene.end,
            "strand": 1 if gene.strand == "+" else -1,
            "name": gene.display_name,
            "type": "gene",
            "description": f"Active gene · {gene.locus_tag}",
            "color": GENE_COLORS["forward"],
        })

    position = (selection or {}).get("position")
    if gene is not None and position:
        reference = str(selection.get("reference") or "N")
        label = selection.get("variant") or selection.get("mutation") or "Selected variant"
        features.append({
            "uniqueId": f"selected-variant-{label}-{position}",
            "refName": gene.chromosome,
            "start": int(position) - 1,
            "end": int(position) - 1 + max(1, len(reference)),
            "name": label,
            "type": "sequence_alteration",
            "description": (
                f"Selected · {int(position):,} {_shorten_allele(reference)}>"
                f"{_shorten_allele(selection.get('alternative') or 'N')}"
            ),
            "color": SELECTION_COLOR,
        })
    return features


def locus_string(chromosome: str, start: int, end: int) -> str:
    return f"{chromosome}:{max(1, int(start))}..{int(end)}"

