# TODO

## Rules
- Keep each item under a category (`BUG FIX` or `FEATURE ADDITION`).
- Use a unique ID for each item.
- Track `Report date`, `Status`, and `Status update date`.
- Allowed `Status` values: `planned`, `in-progress`, `blocked`, `done`.
- Use ISO 8601 date format: `YYYY-MM-DD` (example: `2026-04-29`).

## BUG FIX

### [BUG-006] SUMMARY_REPORT heredoc parsing failure
- Report date: 2026-04-30
- Status: done
- Status update date: 2026-04-30
- Affected files: `module/summary_report.nf`
- Root cause: `SUMMARY_REPORT` process script generated a malformed here-document (`python3 << 'EOF' ... EOF`) where shell parsing failed (`wanted \`EOF\``), leading to Python code being parsed with invalid leading indentation.
- Symptom:
  - `.command.sh: line 79: warning: here-document at line 2 delimited by end-of-file (wanted \`EOF\`)`
  - `IndentationError: unexpected indent`
- Error log:
```text
Error executing process > 'SUMMARY_REPORT'

Caused by:
  Process `SUMMARY_REPORT` terminated with an error exit status (1)

Command error:
  .command.sh: line 79: warning: here-document at line 2 delimited by end-of-file (wanted `EOF')
    File "<stdin>", line 1
      import os
  IndentationError: unexpected indent
```
- Progress:
  - 2026-04-30: adjusted `module/summary_report.nf` heredoc Python block to start at column 0 and avoid shell heredoc termination mismatch.
  - 2026-04-30: fixed embedded Python newline escaping (`\n` -> `\\n`) so generated `.command.sh` does not break string literals at runtime.
  - 2026-04-30: resume rerun completed `SUMMARY_REPORT` successfully; related embedded Python newline issue in `BIOM_GENERATE` also fixed.

### [BUG-005] OTU_MERGE input cardinality warning
- Report date: 2026-04-30
- Status: done
- Status update date: 2026-04-30
- Affected files: `main.nf`, `module/otu_merge.nf`
- Root cause: `main.nf` passed collected file lists to `OTU_MERGE`, but `OTU_MERGE` declared tuple inputs (`tuple val(...), path(...)`), causing input-set cardinality mismatch warnings.
- Fix:
  - Updated `OTU_MERGE` inputs to accept collected file lists directly as `path(consensus_fastas)` and `path(count_files)`.
  - Removed dependency on external `sample_ids` value and derived sample IDs from consensus FASTA filenames inside the Python script.
- Warning log:
```text
WARN: Input tuple does not match input set cardinality declared by process `OTU_MERGE` -- offending value: [/.../barcode31_F02_R01_clustered_consensus.fasta, ...]
```

### [BUG-004] OTU_COUNT_TABLE publish filename collision
- Report date: 2026-04-30
- Status: done
- Status update date: 2026-04-30
- Affected files: `module/otu_count.nf`
- Root cause: `OTU_COUNT_TABLE` published fixed output names (`otu_table_final.tsv`, `otu_table_summary.txt`, `versions_otu_count.yml`) for every sample into one directory (`04_otu_merge`), causing copy collisions and publish warnings.
- Fix:
  - Changed per-sample output filenames to include `sample_id`:
    - `${sample_id}_otu_table_final.tsv`
    - `${sample_id}_otu_table_summary.txt`
    - `${sample_id}_versions_otu_count.yml`
  - Updated process `output:` declarations and script write targets accordingly.
- Warning log:
```text
WARN: Failed to publish file: .../otu_table_final.tsv; to: .../results/04_otu_merge/otu_table_final.tsv [copy] -- See log file for details
```

### [BUG-003] GENERATE_PROVENANCE input filename collision
- Report date: 2026-04-29
- Status: done
- Status update date: 2026-04-29
- Affected files: `module/generate_provenance.nf`
- Root cause: `GENERATE_PROVENANCE` staged a collected list of `versions_*.yml` files with duplicate basenames into one work directory.
- Fix: set `stageAs: 'versions??/*'` on `v_files` input so each file is staged in a separate numbered subdirectory.
- Error log:
```text
Error executing process > 'GENERATE_PROVENANCE'

Caused by:
  Process `GENERATE_PROVENANCE` input file name collision -- There are multiple input files for each of the following file names: versions_cutadapt.yml, versions_otu_count.yml, versions_porechop.yml, versions_seqkit.yml, versions_biopython.yml, versions_krakentools.yml, versions_kraken2.yml, versions_as_main.yml, versions_as_post.yml, versions_as_pre.yml
```

### [BUG-002] OTU_COUNT_TABLE invalid path value
- Report date: 2026-04-28
- Status: done
- Status update date: 2026-04-28
- Affected files: `main.nf`, `module/otu_count.nf`
- Root cause: `OTU_COUNT_TABLE` received `[sample_id, path]` tuple channels but declared bare `path` inputs.
- Fix: join channels by `sample_id` in `main.nf`, and update process input to `tuple val(sample_id), path(consensus_fasta), path(read_count_file)`.
- Error log:
```text
Error executing process > 'OTU_COUNT_TABLE (1)'

Caused by:
  Not a valid path value: 'barcode31_F04_R01'
```

### [BUG-001] AS_DEDUPLICATE unbound raw_fastas
- Report date: 2026-04-28
- Status: done
- Status update date: 2026-04-28
- Affected files: `module/amplicon_sorter.nf`
- Root cause: `\${raw_fastas}` was escaped in the script block, so Bash with `-u` treated `raw_fastas` as unbound.
- Fix: use `${raw_fastas}` so Nextflow interpolates staged input FASTA paths.
- Error log:
```text
Error executing process > 'AS_DEDUPLICATE (barcode31_F02_R01)'

Caused by:
  Process `AS_DEDUPLICATE (barcode31_F02_R01)` terminated with an error exit status (1)

Command error:
  INFO:    Converting SIF file to temporary sandbox...
  .command.sh: line 5: raw_fastas: unbound variable
  INFO:    Cleaning up image...
```

- Executed `.command.sh`:
```bash
#!/bin/bash -ue
# Merge all cluster consensus files and remove exact duplicates.
# If all inputs are empty, create an empty consensus file for downstream zero-count handling.
has_content=0
for f in ${raw_fastas}; do
    if [[ -s "$f" ]]; then
        has_content=1
        break
    fi
done

if [[ "${has_content}" -eq 1 ]]; then
    cat ${raw_fastas} > merged_tmp.fasta
    seqkit rmdup -s merged_tmp.fasta -o barcode31_F02_R01_clustered_consensus.fasta
else
    touch barcode31_F02_R01_clustered_consensus.fasta
fi

echo "SeqKit (AS_POST): $(seqkit version | awk '{print $2}')" > versions_as_post.yml
```

## FEATURE ADDITION

### [FEAT-002] Container image maintenance workflow
  - Centralize container reference URLs in one config source.
  - Add pre-run connectivity checks for configured image URLs.
  - Add candidate tag suggestions when container pull fails.
- Notes:
  - Keep `https://depot.galaxyproject.org/singularity/...` style for Singularity pulls.


## INFRASTRUCTURE / CLEANUP

### [INFRA-001] Prepare trash/ and .gitignore for public release
- Report date: 2026-04-29
- Status: planned
- Status update date: N/A
- Summary: maintain trash/ folder for development work without committing it to public repository.
- Scope:
  - `trash/` directory created to store obsolete files (e.g., `nextflowPipelineGenerationLog.md`).
  - Files moved to `trash/` are tracked in git during development.
  - Before converting private → public repo, add `trash/` to `.gitignore`.
  - Add other development artifacts (e.g., `testContainer.sh` if only for dev testing, `.nextflow/` working dirs).
  - Add `helperScript/runPipeline.ipynb` to `.gitignore` before public release (Linux server specific; requires Jupyter Bash kernel).
- Files currently in trash/:
  - `nextflowPipelineGenerationLog.md` (superseded by README + TODO)
- Next step: Add entry to `.gitignore` when approaching public release.
### [FEAT-001] Zero/partial BLAST-hit handling and NA taxonomy fallback
- Report date: N/A
- Status: done
- Status update date: N/A
- Summary: added error handling for zero or partial BLAST hits and NA taxonomy assignment for unmatched OTUs.
- Changes:
  - `module/blast_annotate.nf`: skip `blastn` and emit empty `all_otus_blast_results.tsv` when input FASTA is empty.
  - `module/create_biom.nf` (`TAXONKIT_LINEAGE`): skip TaxonKit when BLAST is empty and filter invalid `staxids` (``, `0`, `N/A`).
  - `module/create_biom.nf` (`BIOM_GENERATE`): left-join OTU counts with BLAST/TaxonKit; assign `NA;NA;NA;NA;NA;NA;NA` for OTUs with no BLAST hit.
  - `module/summary_report.nf`: fallback taxonomy for no BLAST hit changed from `Unclassified` to `NA`.

### [FEAT-003] BIOM format output generation
- Report date: 2026-04-29
- Status: planned
- Status update date: N/A
- Summary: integrate BIOM_GENERATE and TAXONKIT_LINEAGE processes into main.nf workflow.
- Scope:
  - Include BIOM_GENERATE call in main.nf after BLAST_ANNOTATE.
  - Collect merged OTU count matrix and BLAST annotation results.
  - Generate BIOM-formatted output with integrated taxonomic and abundance data.
  - Output to `06_biom/` directory.
- Notes:
  - TAXONKIT_LINEAGE processes BLAST results to extract taxonomy lineage.
  - BIOM_GENERATE combines counts + lineage into HDF5-based BIOM format.
  - Requires pandas + biom-format containers.

### [FEAT-004] Output folder hierarchy restructuring
- Report date: 2026-04-29
- Status: planned
- Status update date: N/A
- Summary: standardize output folder naming and numbering to match pipeline flow.
- Current structure (irregular):
  ```
  results/
  ├── 00_preprocess
  ├── 01_demux
  ├── 03_amplicon_sorter
  ├── 04_kraken_filter      ← should be 02
  ├── 04_otu_merge          ← duplicate number
  ├── 05_annotation         ← should be 06, rename to blast_annotation
  ├── 00_provenance
  ├── 01_Final_Analysis, 02_Sequences, 03_Annotation (from FINALIZE_RESULTS)
  ```
- Target structure (sequential):
  ```
  results/
  ├── 00_preprocess
  ├── 01_demux
  ├── 02_kraken_filter
  ├── 03_amplicon_sorter
  ├── 04_otu_merge
  ├── 05_blast_annotation
  ├── 06_final_results       ← consolidate from FINALIZE_RESULTS
  ├── 07_biom               ← from BIOM_GENERATE
  ├── 08_summary_report     ← from SUMMARY_REPORT
  └── 99_provenance         ← move to end
  ```
- Affected files: all `module/*.nf` publishDir directives.

### [FEAT-005] SUMMARY_REPORT integration
- Report date: 2026-04-29
- Status: planned
- Status update date: N/A
- Summary: integrate SUMMARY_REPORT process into main.nf workflow.
- Scope:
  - Include SUMMARY_REPORT call after OTU_MERGE and BLAST_ANNOTATE.
  - Collect Kraken2 reports + OTU count matrix + BLAST results.
  - Generate phylum-level summary tables.
  - Output to `08_summary_report/` directory.
- Notes:
  - Produces two outputs: `summary_phylum_kraken2.tsv` and `summary_phylum_otu.tsv`.
  - Used for contamination and taxonomy verification.