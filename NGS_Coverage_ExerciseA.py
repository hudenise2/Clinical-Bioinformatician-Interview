#!/usr/bin/env python3

"""
NGS Coverage Report exercise (part A)

Usage:
    python coverage_report.py \
        --input  <path to NGS148_34_139558_CB_CMCMD_S33_R1_001.sambamba_output.txt>
        --genelist <path to panel_genes.txt>
        --threshold <30 by default>
        --output  <namer of file for output>

Justification of Aggregation strategy:
    The test description indicates that every coding base of every gene has to be
    covered at least 30x. We have access to the % of bases in each exon region
    covered at ≥ 30× (percentage30). Therefore, a gene should be labelled FAIL
    if any of its exons gas value less than 100% in this data column.
    For the exercise, I chose to aggregate using the MINIMUM percentage30 across the
    gene's exons. Using the mean value would let one short 'bad' exon be hidden by 'good'
    exons, which could be clinically unsafe as decisions would be based on incomplete
    information. Using the minimum is the safe choice.
"""

import argparse
import logging
import re
import sys
from pathlib import Path

class Exon:
    """Handle data from one row from the sambamba output"""

    def __init__(self, chrom, start, end, gene_symbol, accession, percentage30):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.gene_symbol = gene_symbol
        self.accession = accession
        self.percentage30 = percentage30

    @property
    def region(self):
        return f"{self.chrom}:{self.start}-{self.end}"

class Gene:
    """Aggregate the exon data and the coverage for a single gene.
    A gene is built from a list of its Exons and the closs' methods report:
      - how many unique exons it has
      - the minimum coverage %
      - the exons that are failing (< 100% coverage)
      -  whether the gene is classified as a PASS or a FAIL from the data
    """

    def __init__(self, symbol):
        self.symbol = symbol
        self.exons = []       # unique exons (deduplicated by region)

    def add_exon(self, exon):
        """Add an exon to the list for one gene
        The lowest percentage30 is kept if a region is already present in the list
        """
        for i, present in enumerate(self.exons):
            if present.region == exon.region:
                if exon.percentage30 < present.percentage30:
                    self.exons[i] = exon
                return
        self.exons.append(exon)

    @property
    def failing_exons(self):
        """return exons below 100% coverage."""
        return [e for e in self.exons if e.percentage30 < 100.0]

    @property
    def min_coverage(self):
        """return the minimum percentage30 across all exons for a gene (or None)."""
        if not self.exons:
            return None
        return min(e.percentage30 for e in self.exons)

    @property
    def worst_exon(self):
        """return the failing exon and its lowest coverage (or None)."""
        if not self.failing_exons:
            return None
        return min(self.failing_exons, key=lambda e: e.percentage30)

    @property
    def status(self):
        """return the status for a gene."""
        if self.failing_exons:
            return "FAIL"
        return "PASS"

    def report_row(self):
        """return a ne tab-separated row per gene for the summary table."""
        min_cov = "NA" if self.min_coverage is None else f"{self.min_coverage:.2f}"
        worst = self.worst_exon
        worst_region = worst.region if worst else "NA"
        return "\t".join([
            self.symbol,
            str(len(self.exons)),
            str(len(self.failing_exons)),
            min_cov,
            worst_region,
            self.status,
        ])

class CoverageReport:
    """Generate the final report: a collection of Gene objects with status."""

    def __init__(self, sample_name, threshold):
        self.sample_name = sample_name
        self.threshold = threshold
        self.genes = {}          # Gene objects named by gene symbol
        self.panel_genes = []    # ordered list obtained from the panel file

    def add_exon(self, exon):
        """Add an exon to Gene objects or creating the Gene object first if not already existing."""
        if exon.gene_symbol in self.genes:
            gene = self.genes[exon.gene_symbol]
        else :
            gene = Gene(exon.gene_symbol)
            self.genes[exon.gene_symbol] = gene
        gene.add_exon(exon)

    def check_panel(self, panel_genes):
        """"return the gene(s) from the panel list not present in the Gene objects created from the sambamba file."""
        self.panel_genes = list(panel_genes)
        self.missing_genes = [symbol for symbol in self.panel_genes if symbol not in self.genes]
        return self.missing_genes

    @property
    def failing(self):
        """return as sorted list of Gene objects according to their min coverage"""
        failing = [g for g in self.genes.values() if g.status == "FAIL"]
        return sorted( failing,
            key=lambda g: g.min_coverage if g.min_coverage is not None else 999.0)

    @property
    def total(self):
        """return the total number of Gene objects"""
        return len(self.panel_genes) if self.panel_genes else len(self.genes)

    def summary_line(self):
        """return the run-level summary line"""
        line = (
            f"Sample {self.sample_name}: {len(self.failing)}/{self.total} "
            f"genes below {self.threshold:g}x threshold")
        return line

    def create_report(self, missing_genes):
        """return the header line and the data lines constituting the report"""
        lines = [
            f"# NGS Coverage Report — Sample {self.sample_name}",
            f"# Threshold: every coding base >= {self.threshold:g}x",
            "",
            "GeneSymbol\tTotalExons\tFailingExons\tMinCoverage%\t"
            "WorstExon(chr:start-end)\tStatus"]

        for gene in self.failing:
            lines.append(gene.report_row())
        lines.append("")
        lines.append(self.summary_line())

        # Missing genes
        if missing_genes:
            lines.append("")
            lines.append("# Genes in panel with no coverage data (MISSING)")
            lines.append("GeneSymbol\tTotalExons\tFailingExons\tStatus")
            for symbol in missing_genes:
                lines.append(symbol+"\t-\t-\ttMISSING")
        return "\n".join(lines) + "\n"

# Helper functions
def parse_sambamba(path):
    """Read the sambamba file and return a list of Exon objects.
    Malformed rows are logged and excluded.
    - input path (str): path to the sambamba file
    - return exons (list): list of Exon objects
    """
    exons = []
    with open(path, "r", encoding="utf-8") as f:
        line_count=0
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#chr"):
                # Using whitespace here as header line contains a mix of space and tab.
                # separating the field names.
                # Note that this is not the case for the data lines where only tab were used
                nber_fields=len(line.split())
            else:
                parts = line.split("\t")
                line_count += 1
                if len(parts) < nber_fields:
                    logging.warning(
                        f"Line {line_count}: too few columns to process: {len(parts)}/{nber_fields}; skipping")
                    continue
                gene_field = parts[6]
                if not gene_field or gene_field in (".", "NA", "N/A", "-"):
                    logging.warning(
                        f"Line {line_count}: missing gene symbol; skipping")
                    continue
                elif ";" in gene_field:
                    gene, accession = gene_field.split(";", 1)
                else:
                    gene = gene_field
                    accession = ""
                gene=gene.strip()
                if len(gene) == 0:
                    logging.warning(f"Line {line_count}: missing gene symbol; skipping")
                    continue
                try:
                    pct_field = float(parts[10].strip())
                except:
                    logging.warning(f"Line {line_count}: non numeric percentage30 value; skipping")
                    continue
                chr = parts[0].strip()
                start = int(parts[1])
                end = int(parts[2])
                exons.append(
                    Exon(
                        chrom=chr,
                        start=start,
                        end=end,
                        gene_symbol=gene,
                        accession=accession,
                        percentage30=pct_field))
    return exons

def load_panel_genes(path):
    """Read the panel gene list (contains one gene symbol per line).
    - input path (str): path to the panel gene file
    - return genes (list): list of gene symbol
    """
    genes = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            symbol = line.strip()
            if symbol and not symbol.startswith("#"):
                genes.append(symbol)
    return genes


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Generate an NGS coverage report from sambamba output provided." )
    parser.add_argument("--input", required=True,
                        help="path to sambamba coverage output file")
    parser.add_argument("--genelist", required=True,
                        help="path to panel gene list")
    parser.add_argument("--threshold", type=float, default=30.0,
                        help="coverage threshold in x (default: 30)")
    parser.add_argument("--output", required=True,
                        help="output report file name")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.WARNING, format="%(levelname)s: %(message)s")

    sample_name = "_".join(args.input.split("/")[-1].split("_")[:2])
    # Parse the sambamba file into Exon objects.
    exons = parse_sambamba(args.input)

    # Build the report, feeding each exon into the appropriate Gene object.
    report = CoverageReport(sample_name=sample_name, threshold=args.threshold)
    for exon in exons:
        report.add_exon(exon)

    # Cross-reference the panel list (identify the missing genes).
    panel_genes = load_panel_genes(args.genelist)
    missing_genes= report.check_panel(panel_genes)

    # Generate the report lines and write them to output file.
    text = report.create_report(missing_genes)
    Path(args.output).write_text(text, encoding="utf-8")
    print(f"Report written to {args.output}")

    # Print the run-level summary line to stdout as required.
    print(report.summary_line())

    # Optional: give the possibility of changing the script return code if particular
    # condition (too many failed or missing genes for example) is encountered.
    # Could be used to check the status of the run before continuing to downstrram analysis.
    return 0


if __name__ == "__main__":
    sys.exit(main())
