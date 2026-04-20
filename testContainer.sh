Porechop
singularity run https://depot.galaxyproject.org/singularity/porechop:0.2.4--py39h2de1943_9 porechop -h

SeqKit
singularity run https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0 seqkit version

Cutadapt
singularity run https://depot.galaxyproject.org/singularity/cutadapt:5.0--py39hbcbf7aa_0 cutadapt --version

biopython
singularity run https://depot.galaxyproject.org/singularity/biopython:1.79

Kraken2
singularity run https://depot.galaxyproject.org/singularity/kraken2:2.1.3--pl5321h077b44d_4 kraken2 -version

extract_kraken_reads.py
singularity run docker://nanozoo/krakentools:1.2--13d5ba5 extract_kraken_reads.py -h

Pandas
singularity run docker://quay.io/biocontainers/pandas:2.2.1 

BLAST
singularity run https://depot.galaxyproject.org/singularity/blast:2.16.0--h66d330f_4 blastn -help

TaxonKit
singularity run https://depot.galaxyproject.org/singularity/taxonkit:0.18.0--h9ee0642_0 taxonkit -h

biom-format
singularity run https://depot.galaxyproject.org/singularity/biom-format:2.1.15
