# Pre-Interview Exercise: NGS Coverage Report

## Background

The laboratory offers a diagnostic test for patients with suspected congenital myopathy or
congenital muscular dystrophy. The test sequences **84 genes** associated with these conditions
using a capture-based NGS approach (Agilent SureSelect pull-down, Illumina NextSeq).

Performance metrics for this test require that **every coding base of each gene is covered by
at least 30 reads** (≥ 30×). The analytical pipeline uses **sambamba** to generate per-exon
coverage statistics. The `percentage30` column in the sambamba output records the percentage
of bases in that region covered at ≥ 30×.

Your task is to write a Python script that processes this output and produces a coverage report
highlighting any genes that do not meet the 100% at 30× threshold.

---

## Files provided

| File | Description |
|---|---|
| `NGS148_34_139558_CB_CMCMD_S33_R1_001.sambamba_output.txt` | Example sambamba coverage output for one patient sample |
| `panel_genes.txt` | The 84 genes expected to be covered by this panel (one symbol per line) |

---

## Time guidance

The total time allowed for this exercise is **2.5 hours**.

- **Part A** should be completable in approximately **90 minutes**. Most of that time will be
  spent understanding the input format, deciding on an aggregation strategy, and handling
  edge cases — not on the volume of code.
- **Part B** is an extension. You are not expected to finish it within the time limit. It is
  there to give you room to demonstrate deeper thinking if you have time remaining.

If you find yourself running over on Part A, prioritise a working solution over a polished
one — correct output with rough edges is more informative than incomplete code that doesn't run.

---

## Part A — Core report (required)

Write a Python script that:

1. **Parses** the sambamba output file.

2. **Handles the gene name column.** The `GeneSymbol;Accession` column (column 7) is a compound
   field containing both a gene symbol and a transcript accession separated by a semicolon.
   You must split these and use only the gene symbol. Note that some genes are represented by
   more than one transcript accession in the panel — these should be treated as a single gene
   in your report.

3. **Aggregates** per-exon coverage to a per-gene result. You should justify your aggregation
   strategy in a comment in your code. Consider carefully: which approach correctly identifies
   a gene as failing when only a single short exon is below threshold?

4. **Reports** any gene with less than 100% coverage at 30×, using this output format:

   ```
   GeneSymbol  TotalExons  FailingExons  MinCoverage%  WorstExon(chr:start-end)  Status
   ```

   Additionally, print a **run-level summary line**, for example:
   ```
   Sample NGS148_34: 3/84 genes below 30x threshold
   ```

5. **Cross-references** the results against `panel_genes.txt`. Any gene in the panel list that
   has no rows in the sambamba output should be reported as `MISSING` — this can indicate a
   capture failure and must not be silently omitted.

6. **Handles malformed input gracefully.** The input may contain rows with missing gene symbols
   or non-numeric values in the `percentage30` column. These should be logged as warnings and
   skipped rather than causing a crash.

7. **Uses `argparse`** so the script can be run as follows:

   ```bash
   python coverage_report.py \
       --input  NGS148_34_139558_CB_CMCMD_S33_R1_001.sambamba_output.txt \
       --genelist panel_genes.txt \
       --threshold 30 \
       --output  report.txt
   ```

   The coverage threshold (`--threshold`) must be configurable — do not hardcode `30`.

---

## Part B — Extension (if time permits)

8. For any failing gene, report which **specific exon(s)** are below threshold, including their
   genomic coordinates (`chr:start-end`). A clinical scientist reviewing the report needs to
   know exactly which regions require follow-up (manual review, orthogonal assay, or Sanger
   sequencing).

9. Briefly **discuss** (in a comment block or a short accompanying note) the difference between
   using a minimum vs. a weighted-mean aggregation strategy, and why one is more appropriate
   in a diagnostic context.

10. Write at least **two unit tests** covering your core logic (e.g. aggregation, gene-symbol
    parsing, missing-gene detection).

---

## Notes on the input data

- The `Size` column is provided but its meaning is not guaranteed — inspect it carefully before
  using it in any calculation.
- The input contains a small number of deliberately anomalous rows; your script should handle
  these without crashing.
- Some genes are represented by more than one transcript accession; your report should still
  produce one row per gene.
