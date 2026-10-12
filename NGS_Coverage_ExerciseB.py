#!/usr/bin/env python3
"""
NGS Coverage Report exercise (part B)

I- For any failing gene, report which **specific exon(s)** are below threshold, including their
   genomic coordinates (`chr:start-end`).
   In NGS Coverage Report exercise (part A), I already report the coordinates of the worst failing exon(s).
   The task is an extension of that so sumary to the changes to the script for exercice A:
   - Addition of a method to the Gene class:
   
        def detail_lines(self):
            # generating lines for the failed exons.
            lines = []
            for exon in sorted(self.failing_exons, key=lambda e: e.percentage30):
                lines.append(f"    {exon.region}\t{exon.percentage30:.2f}%\t"
                    f"accession={exon.accession or 'NA'}")
        return lines

    - addition of a section to the create_report method from the CoverageReport class:

        def create_report(self):
            # return the header line and the data lines consistuting the report
            lines = [
                f"# NGS Coverage Report — Sample {self.sample_name}",
                f"# Threshold: every coding base >= {self.threshold:g}x",
                "",
                "GeneSymbol\tTotalExons\tFailingExons\tMinCoverage%\t"
                "WorstExon(chr:start-end)\tStatus",
            ]
    
            for gene in self.failing:
                lines.append(gene.report_row())
    
    
            # ADDED SECTION: add exon detail for each failing gene to report lines
            # Output of these changes had been added to the GitHub repository as outputB.txt
            detail = []
            for gene in self.failing:
                detail.append(f"- {gene.symbol}:")
                detail.extend(gene.detail_lines())
            if detail:
                lines.append("")
                lines.append("# Failing exons (Part B)")
                lines.extend(detail)
    
    
            lines.append("")
            lines.append(self.summary_line())
    
            # Missing genes (separate section)
            if self.missing_genes:
                lines.append("")
                lines.append("# Genes in panel with no coverage data (MISSING)")
                for symbol in self.missing_genes:
                    lines.append(symbol+"\t-\t-\t-\t-\tMISSING")
            return "\n".join(lines) + "\n"


    II- Briefly **discuss** (in a comment block or a short accompanying note) the difference between
   using a minimum vs. a weighted-mean aggregation strategy, and why one is more appropriate
   in a diagnostic context
   
   Using a weighted-mean aggregation strategy means computing an average value in which exons, in our case,
   would contribute according to their lengths or read count (the 'weight' criteria present in the data).
   It would have the effect of masking 'bad' exon(s) as the gene coverage average could still be close enough
   to 100% to be included despite having a region where we could not report on.
   In contrast, using the minimum ensure that only genes with full coverage above the threshold are included.
   
   III- Write at least **two unit tests** covering your core logic (e.g. aggregation, gene-symbol
    parsing, missing-gene detection)

    Two very simple unit tests can be run with the script below.
    It used classes and methods from the part A script and dummy files to test their behaviour.
    Due to time constraint, I did not script an aggregation test but the principle will be the
    same: create an example as a dummy file, process it using the part A appropriate components
    and compare the output with expectation using the check method.

"""

import os
# import classes and methods from the part A script
from NGS_Coverage_ExerciseA import (
    Exon,
    Gene,
    CoverageReport,
    parse_sambamba,
    load_panel_genes,
)

# Helper functions

def check(label, actual, expected):
    """Compare actual and expected value and print PASS/FAIL status.
    - input label (str): name of the test performed
    - input actual: output from the mock file processing (could be different datatypes)
    - input expected: expected output (could be different datatypes)
    - return None
    """
    if actual == expected:
        print(f"PASS: {label}")
    else:
        print(f"FAIL: {label}")
        print(f"      expected: {expected!r} and  got: {actual!r}")


def write_tmp(text, file):
    """Write text to a temp file and return its path.
    - input text (str): lines of text to be written in dummy file
    - input file (TextIOWrapper): dummy file open for writing the text.
     - return file.name (str): name of the dummy file
     """
    file.write(text)
    file.close()
    return file.name

# Test 1: parsing the compound 'GeneSymbol;Accession' column appropriately
def test_gene_annotation_splitting():
    """The parser function must split 'SYMBOL;ACCESSION' and store both parts."""
    # create mock data file
    content = (
        "#chromosome\tStart\tEnd\tSize\tReadCount\t\t"
        "GeneSymbol;Accession\t\t\t\tpercentage30\n"
        "chr1\t100\t200\t100\t50\t\tMYH7;NM_000257.3\t\t\t\t100.0\n"
        "chr1\t300\t400\t100\t50\t\tMYH7.;NM_000257.4\t\t\t\t95.0\n"
        "chr1\t500\t600\t100\t50\t\tTTN;NM_001267550.2\t\t\t\t80.0\n"
    )
    tempfile = open("temp_test", "w")
    path = write_tmp(content, tempfile)
    try:
        exons = parse_sambamba(path)
    finally:
        os.unlink(path)

    # Three rows should yield three Exon objects.
    check("parser: number of exons", len(exons), 3)
    # Extract the Gene symbol (part before ';').
    check("parser: first gene symbol", exons[0].gene_symbol, "MYH7")
    check("parser: second gene symbol", exons[1].gene_symbol, "MYH7.1")
    check("parser: third gene symbol", exons[2].gene_symbol, "TTN")
    # Extract the accession (part after ';').
    check("parser: first accession", exons[0].accession, "NM_000257.3")
    check("parser: second accession", exons[1].accession, "NM_000257.4")
    check("parser: third accession", exons[2].accession, "NM_001267550.2")

# Test 2: identifying the missing genes (in panel but absent from sambamba output)
def test_missing_gene_is_detected():
    """Full pipeline on a tiny synthetic input."""
    content = (
        "#chromosome\tStart\tEnd\tSize\tReadCount\t\t"
        "GeneSymbol;Accession\t\t\t\tpercentage30\n"
        "chr1\t100\t200\t100\t50\t\tMYH7;NM_000257.3\t\t\t\t100.0\n"
        "chr1\t300\t400\t100\t50\t\tMYH7.;NM_000257.4\t\t\t\t95.0\n"
        "chr1\t500\t600\t100\t50\t\tTTN;NM_001267550.2\t\t\t\t80.0\n"
    )
    panel_text = "MYH7\nTTN\nGHOST\n"
    temp_data = open("temp_data_test", "w")
    input_path = write_tmp(content, temp_data)
    temp_panel = open("temp_panel_test", "w")
    panel_path = write_tmp(panel_text, temp_panel)
    try:
        report = CoverageReport("NGS148_34", threshold=30.0)
        for exon in parse_sambamba(temp_data.name):
            report.add_exon(exon)
        report.check_panel(load_panel_genes(panel_path))
    finally:
        os.unlink(input_path)
        os.unlink(panel_path)
    check("check_panel: missing genes", report.missing_genes, ['GHOST'])

if __name__ == "__main__":
    tests = [
        test_gene_annotation_splitting,
        test_missing_gene_is_detected ]
    for t in tests:
        print(f"\n--- {t.__name__} ---")
        t()
    print("\nAll tests completed.")
