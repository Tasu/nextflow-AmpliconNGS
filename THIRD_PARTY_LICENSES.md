# Third-Party Tools and Data Licensing Notes

This repository (workflow code and helper scripts) is distributed under the MIT License.

The workflow executes third-party tools and uses third-party databases. Those components are licensed separately by their original authors, and their terms remain in effect.

## Important points

- Using this workflow does not relicense third-party tools or databases.
- Running tools via Singularity/Apptainer containers does not remove or weaken third-party license obligations.
- If you redistribute container images, you are responsible for complying with each included component license (for example, notices, source code availability requirements, and attribution where required).
- If users pull containers directly from upstream registries and you do not redistribute them, each user still needs to comply with the relevant upstream licenses.

## Components referenced by this pipeline

The default configuration references tools such as:

- Porechop
- SeqKit
- Cutadapt
- Biopython
- Kraken2
- KrakenTools
- BLAST+
- TaxonKit
- BIOM format tooling

Database content is also separately licensed, including but not limited to NCBI BLAST and taxonomy datasets.

## Recommended compliance practice

- Keep this file and the upstream project links in your documentation.
- Record exact tool versions (the pipeline already writes versions files).
- Before publication or redistribution, verify current upstream license terms for each tool and dataset version in use.

This document is informational and not legal advice.
