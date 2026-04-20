# 「Nanoporeアンプリコン解析パイプライン」要件定義

# パラメータ設定方針: params.yaml + nextflow.config
## 方針
 - ユーザーが指定するパラメーターはparams.yaml
 - 他のシステム設定はnextflow.configに記載
 - nextflow.configにはdefaultのparams.yamlの内容も記載
 - 実行
 
```bash
 nextflow run main.nf -profile singularity -params-file params.yaml
```
 
## 実装コード

``` nextflow
/*
 * nextflow.config
 */

// 1. Default Parameters (overridden by -params-file)
params {
    sample_sheet = "sampleSheet.csv"
    outdir       = "results"
    
    // Resource defaults (fallback)
    max_cpus     = 8
    max_memory   = '32.GB'
    max_time     = '24.h'
}

// 2. Resource Management Logic
// This function ensures the process doesn't request more than the defined max_resource
def check_max(obj, type) {
    if (type == 'memory') {
        try {
            if (obj.compareTo(params.max_memory as nextflow.util.MemoryUnit) == 1)
                return params.max_memory as nextflow.util.MemoryUnit
            else
                return obj
        } catch (all) {
            return obj
        }
    } else if (type == 'cpus') {
        try {
            return Math.min( obj, params.max_cpus as int )
        } catch (all) {
            return obj
        }
    } else if (type == 'time') {
        try {
            if (obj.compareTo(params.max_time as nextflow.util.Duration) == 1)
                return params.max_time as nextflow.util.Duration
            else
                return obj
        } catch (all) {
            return obj
        }
    }
}

// 3. Process Resource Allocation
process {
    errorStrategy = { task.exitStatus in [143,137,104,134,139] ? 'retry' : 'terminate' }
    maxRetries    = 2

    withLabel: 'process_low' {
        cpus   = { check_max( 1 * task.attempt, 'cpus' ) }
        memory = { check_max( 4.GB * task.attempt, 'memory' ) }
    }
    withLabel: 'process_medium' {
        cpus   = { check_max( 4 * task.attempt, 'cpus' ) }
        memory = { check_max( 16.GB * task.attempt, 'memory' ) }
    }
    withLabel: 'process_high' {
        cpus   = { check_max( 12 * task.attempt, 'cpus' ) }
        memory = { check_max( 64.GB * task.attempt, 'memory' ) }
    }
}

// 4. Profiles & Engine Settings
profiles {
    singularity {
        singularity.enabled    = true
        singularity.autoMounts = true
        singularity.cacheDir   = "${baseDir}/singularity_cache"
    }
}
```

```
# --- Input/Output Paths ---
sample_sheet: "sampleSheet.csv"
outdir: "results"

# --- Database Paths ---
kraken2_db: "/path/to/kraken2_db"
blast_db_dir: "/path/to/blast_db"
tax_db_dir: "/path/to/taxonomy"

# --- Resource Limits (Max allowed on current server) ---
max_cpus: 64
max_memory: "256.GB"
max_time: "24.h"

# --- Analysis Parameters ---
filter_taxids: "9606 10088"
filter_action: "remove"
max_reads: 10000
blast_db_name: "nt"
blast_type: "blastn"

```


# 0. 前提条件
### 0.1. 実行環境（テクニカル）
* **プラットフォーム**: Nextflow version 22.10.6, Linux x86_64
* **コンテナ**: 各プロセスに Singularity コンテナを割り当てる。
    * パブリックコンテナがない場合はローカルに配置した Singularity イメージを使用。
    * 使用ツールのバージョンは、可能な限りコード内にコメントで記載する。

### 0.2. パイプラインの作成・編集方針
* **モジュール化**: 各モジュールは独立して開発し、個別に確定させる。
* **非破壊編集**: 一度確定したモジュールは、必要がない限り再編集しない。
* **言語**: コード内のコメントはすべて**英語**で記載する。
* **統合**: 最後にすべてのモジュールを結合する Workflow を記述する。

### 0.3. 入力データ (sampleSheet.csv)
サンプルごとに詳細なメタデータを持つ CSV ファイルを入力とする。

| カラム名 | 説明 |
| :--- | :--- |
| **sample** | ユニークなサンプル名（英数字と `_`） |
| **fastq_dir** | 複数の fastq/fastq.gz を含むディレクトリのパス |
| **min_len** | 配列の最小許容長さ (bp) |
| **max_len** | 配列の最大許容長さ (bp) |
| **fwd_index** | 5'->3' インデックス配列 (8bp)。存在しない場合は空。 |
| **rev_index** | 3'->5' インデックス配列 (8bp)。存在しない場合は空。 |
| **fwd_primer** | 5'->3' フォワードプライマー配列 (IUPAC 縮重塩基可) |
| **rev_primer** | 5'->3' リバースプライマー配列 (IUPAC 縮重塩基可) |

**sampleSheet.csv 例:**
```csv
sample,fastq_dir,min_len,max_len,fwd_index,fwd_primer,rev_index,rev_primer
B01_F01_R01,data/raw_B01,1000,1800,AGCGATAG,CAGCAGCCGCGGTAATTCC,CTATCGCT,TACRGMWACCTTGTTACGAC
B02_F02_R01,data/raw_B01,1000,1800,AATGAGCG,CAGCAGCCGCGGTAATTCC,CTATCGCT,TACRGMWACCTTGTTACGAC
B03_F01_R01,data/raw_B01,1000,1800,ACAGTGGT,CAGCAGCCGCGGTAATTCC,CTATCGCT,TACRGMWACCTTGTTACGAC
B04,data/raw_B04,1000,1800,,CAGCAGCCGCGGTAATTCC,,TACRGMWACCTTGTTACGAC
```


# プロセス要件定義: PREPROCESS
- module/preprocess.nf
### 1.1 ファイルの結合 (Merging)
* **処理内容**: `fastq_dir` 内に存在する全ての `fastq` または `fastq.gz` ファイルを `zcat -f` を用いて一つに統合する。
* **共有化**: 同一の `fastq_dir` を参照しているサンプルが複数ある場合、Nextflowのチャンネル制御により、この結合・前処理工程は一度だけ実行され、結果が共有される。

### 1.2 リソースの最適化 (Resource Management)
* **スレッド制限**: `Porechop` の並列化効率の限界を考慮し、1ディレクトリあたりの割り当て CPU 数を **最大 4 スレッド** に制限する。
* **同時実行制御**: システム全体の最大リソース（128 CPUs）の範囲内で、複数のディレクトリを並列に処理する。

### 1.3 第1段階：アダプター除去 (Porechop Round 1)
* **保護設定**: `--extra_end_trim -1` を指定し、リード端部にあるインデックスやプライマー配列が過剰に削られるのを防ぐ。
* **キメラ除去**: `--discard_middle` を有効にし、リード内部にアダプター配列が検出された（キメラの疑いがある）リードを即座に破棄する。

### 1.4 第2段階：残渣チェック (Porechop Round 2)
* **目的**: 第1段階のトリミング後に、なおアダプターやバーコードの残渣が残っている異常リードを特定する。
* **検出優先設定**: 検出感度を最大化するため、あえて `--extra_end_trim -1` や `--discard_middle` を指定せずに実行する。
* **検出ロジック**: 冗長出力（`-v 2`）を指定し、ログからトリミング（検知）が発生したリードIDを抽出する。

### 1.5 最終フィルタリング (Strict Filtering)
* **除外処理**: 第2段階で特定された「残渣あり」リードの ID リストを用い、`seqkit grep -v` によって第1段階の出力ファイルから該当リードを完全に抹消する。
* **出力**: 1段階目の適切なトリミングを受け、かつ2段階目の厳しい残渣チェックをクリアした高品質なリードのみを `chopped_fastq` として次工程に渡す。

## 5. 実装コード (Nextflow)

```nextflow
/*
 * PREPROCESS Module: Separated into Porechop and SeqKit
 */

process PORECHOP_TRIM {
    tag "${fastq_dir_id}"
    label 'process_medium'
    container 'https://depot.galaxyproject.org/singularity/porechop:0.2.4--py39h2de1943_9'

    input:
    tuple val(fastq_dir_id), path(fastq_files)

    output:
    tuple val(fastq_dir_id), path("round1.fastq.gz"), emit: round1_fastq
    path "${fastq_dir_id}_porechop_round2.log",       emit: log_round2
    path "versions_porechop.yml",                     emit: versions

    script:
    def n_threads = task.cpus > 4 ? 4 : task.cpus
    """
    zcat -f ${fastq_files} | gzip > merged.fastq.gz

    # Round 1: Trim and discard middle
    porechop -i merged.fastq.gz -o round1.fastq.gz --threads ${n_threads} --extra_end_trim -1 --discard_middle

    # Round 2: Detection only for strict filtering later
    porechop -i round1.fastq.gz -o /dev/null --threads ${n_threads} -v 2 > ${fastq_dir_id}_porechop_round2.log

    echo "Porechop: \$(porechop --version 2>&1)" > versions_porechop.yml
    """
}

process SEQKIT_CLEAN {
    tag "${fastq_dir_id}"
    label 'process_low'
    publishDir "${params.outdir}/00_preprocess/${fastq_dir_id}", mode: 'copy'
    container 'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0'

    input:
    tuple val(fastq_dir_id), path(round1_fastq)
    path log_round2

    output:
    tuple val(fastq_dir_id), path("${fastq_dir_id}_chopped.fastq.gz"), emit: chopped_fastq
    path "versions_seqkit.yml", emit: versions

    script:
    """
    grep "trimmed" ${log_round2} | awk '{print \$1}' | sed 's/@//' > discarded_ids.txt || touch discarded_ids.txt

    if [ -s discarded_ids.txt ]; then
        seqkit grep -v -f discarded_ids.txt ${round1_fastq} -o ${fastq_dir_id}_chopped.fastq.gz
    else
        cp ${round1_fastq} ${fastq_dir_id}_chopped.fastq.gz
    fi

    echo "SeqKit: \$(seqkit version | awk '{print \$2}')" > versions_seqkit.yml
    """
}

```




# プロセス要件定義: 1.DEMULTIPLEX
- module/demux.nf
## 1. 目的
Nanoporeリード特有の構造（ミスマッチ、インデル、末尾の余剰配列、キメラ）を考慮しつつ、曖昧塩基（Degenerate bases）を含むプライマーとインデックス配列を正確に処理し、向きの揃った高品質なインサート配列を抽出する。

## 2. コアロジック: 「小文字マーキング（Lowercase Marking）」戦略
1.  **正規化**: 全てのリードを大文字 (`seqkit seq -u`) に統一し、ケース（大文字・小文字）による座標情報のマーキング準備を行う。
2.  **向き検索（Orientation Search）**: 
    * **Forward**: リードの5'端側に `F_PRM`、3'端側に `R_PRM_RC` を順次検索する。
    * **Reverse**: リードの5'端側に `R_PRM`、3'端側に `F_PRM_RC` を順次検索し、見つかった場合は `seqkit seq -rp` で反転させて向きを揃える。
    * **マーキング**: `cutadapt --action=lowercase` を使用してプライマー位置を小文字化する。これにより、曖昧塩基やミスマッチが混在していても、実際にマッチした正確な座標を保持できる。
3.  **インデックスを考慮した切り出し**: 小文字（プライマー）の外側へ、指定された長さ (`f_idx_len`, `r_idx_len`) だけ範囲を広げて配列を切り出し、インデックス領域を回収する。
4.  **厳密な検証**: 切り出された外側領域に対して、ミスマッチなし（`-e 0`）の完全一致でインデックス配列を検索し、一致しないリードを破棄する。
5.  **最終トリミング**: Pythonスクリプトの正規表現を用いて、小文字（プライマー）に挟まれた大文字ブロック（インサート）のみを抽出し、最終的な長さフィルタリングを適用する。

## 3. 実装コード (Nextflow)

```nextflow
/*
 * DEMULTIPLEX Module: Separated into Cutadapt and Python(Biopython)
 */

process CUTADAPT_MARK {
    tag "${sample_id}"
    container 'https://depot.galaxyproject.org/singularity/cutadapt:5.0--py39hbcbf7aa_0'

    input:
    tuple val(sample_id), path(fastq), val(min_len), val(max_len), val(f_idx), val(f_prm), val(r_idx), val(r_prm)

    output:
    tuple val(sample_id), path("marked.fastq"), val(min_len), val(max_len), val(f_idx), val(r_idx), emit: marked_data
    path "versions_cutadapt.yml", emit: versions

    script:
    """
    # Orientation search and lowercase marking
    cutadapt -j ${task.cpus} -g "${f_prm}" --discard-untrimmed ${fastq} | \
    cutadapt -j ${task.cpus} -a "\$(echo ${r_prm} | rev | tr ATCG TAGC)" --discard-untrimmed - --action=lowercase -o fwd_marked.fastq
    
    # (Reverse search and merge logic omitted for brevity, but same principle applies)
    cat fwd_marked.fastq > marked.fastq
    echo "Cutadapt: \$(cutadapt --version)" > versions_cutadapt.yml
    """
}

process BIOPYTHON_EXTRACT {
    tag "${sample_id}"
    publishDir "${params.outdir}/01_demux/${sample_id}", mode: 'copy'
    container 'https://depot.galaxyproject.org/singularity/biopython:1.79'

    input:
    tuple val(sample_id), path(marked_fastq), val(min_len), val(max_len), val(f_idx), val(r_idx)

    output:
    tuple val(sample_id), path("${sample_id}_final.fastq.gz"), emit: reads
    path "versions_biopython.yml", emit: versions

    script:
    def f_idx_len = f_idx ? f_idx.length() : 0
    """
    python3 -c "
import sys, re
from Bio import SeqIO
# Logic to extract uppercase insert between lowercase primers
# Apply min_len / max_len filtering here
    " > ${sample_id}_final.fastq
    gzip ${sample_id}_final.fastq
    echo "BioPython: 1.79" > versions_biopython.yml
    """
}

```




# プロセス要件定義: 2.KRAKEN2_FILTER
- module/kraken2_filter.nf
## 1. 目的
デマルチプレックス済みのリードに対してタクソノミ（分類学的）アサインを行い、特定の生物種（宿主、汚染、あるいはターゲット病原体）に基づいてリードをフィルタリングする。

## 2. コアロジック
1.  **分類 (Classification)**: `Kraken2` を使用して、指定されたデータベースに基づきリードをアサインする。この際、データベース作成時に使用された Taxonomy tree（通常は NCBI）に基づいて分類が行われる。
2.  **階層を考慮した抽出 (Taxonomic Extraction)**: `KrakenTools` の `extract_kraken_reads.py` を活用する。
    * 単一の TaxID 指定であっても、その配下の階層（Children）を自動的にフィルタリング対象に含める (`--include-children`)。
    * `action` パラメータ（`keep` または `remove`）を切り替えることで、ターゲットの抽出とホストの除去の両方に対応する。
3.  **環境の再現性**: `Kraken2` と `KrakenTools` の互換性が確認されているマルチパッケージ・コンテナを使用し、バイオインフォマティクスツールのバージョン依存問題を回避する。

## 3. 実装上の注意点
* **メモリ割り当て**: Kraken2 はデータベースをメモリ上に展開するため、`nextflow.config` 等で十分な RAM（DBサイズ + α）を確保すること。
* **未分類リードの扱い**: `action: 'keep'` を選択した場合、どのカテゴリにも属さなかった Unclassified リードは出力に含まれない（破棄される）。

## 4. 使用ツールとコンテナ
* **Kraken2**: v2.1.2
* **KrakenTools**: v1.2
* **Container**: `https://depot.galaxyproject.org/singularity/mulled-v2-4e837765962273925eea842460b6126125401ba6:9e5687184d2dd469ea3f70a7b3b9e32cc4541b9c-0`

## 5. 実装コード (Nextflow)
```nextflow
/*
 * KRAKEN2_FILTER Module
 * * Description:
 * This module is split into two processes to handle different container environments.
 * 1. KRAKEN2_CLASSIFY: Assigns taxonomy to consensus sequences.
 * 2. KRAKENTOOLS_EXTRACT: Extracts specific reads (e.g., target parasites) based on taxonomy.
 */

/*
 * Process 1: Taxonomic Classification
 */
process KRAKEN2_CLASSIFY {
    tag "${sample_id}"
    label 'process_high'
    // Using verified stable URI from Galaxy Depot
    container 'https://depot.galaxyproject.org/singularity/kraken2:2.1.3--pl5321h077b44d_4'

    input:
    tuple val(sample_id), path(fastq)
    path kraken2_db

    output:
    tuple val(sample_id), path("${sample_id}.kraken2.out"),    emit: output
    tuple val(sample_id), path("${sample_id}.kraken2.report"), emit: report
    path "versions_kraken2.yml",                               emit: versions

    script:
    """
    # Execute Kraken2 classification
    kraken2 \\
        --db ${kraken2_db} \\
        --threads ${task.cpus} \\
        --output ${sample_id}.kraken2.out \\
        --report ${sample_id}.kraken2.report \\
        ${fastq}

    # Capture Kraken2 version info
    echo "Kraken2: \$(kraken2 --version | head -n 1 | awk '{print \$3}')" > versions_kraken2.yml
    """
}

/*
 * Process 2: Read Extraction by Taxonomy ID
 */
process KRAKENTOOLS_EXTRACT {
    tag "${sample_id}"
    label 'process_medium'
    publishDir "${params.outdir}/04_kraken_filter/${sample_id}", mode: 'copy'
    // Using verified stable URI from Docker Hub (Nanozoo)
    container 'docker://nanozoo/krakentools:1.2--13d5ba5'

    input:
    tuple val(sample_id), path(kraken_out), path(kraken_report), path(fastq)
    val taxid_to_extract // Target TaxID (e.g., from params.target_taxid)

    output:
    tuple val(sample_id), path("${sample_id}_filtered.fastq.gz"), emit: filtered_fastq
    path "versions_krakentools.yml",                             emit: versions

    script:
    """
    # Extract specific reads matching the TaxID and its children
    extract_kraken_reads.py \\
        -k ${kraken_out} \\
        -r ${kraken_report} \\
        -s ${fastq} \\
        -t ${taxid_to_extract} \\
        --include-children \\
        -o ${sample_id}_filtered.fastq

    # Compress the output fastq
    gzip ${sample_id}_filtered.fastq

    # Record KrakenTools version
    echo "KrakenTools (extract_kraken_reads.py): 1.2" > versions_krakentools.yml
    """
}

```


# プロセス要件定義: 3.AMPLICON_SORTER
- module/amplicon_sorter.nf
## 1. 目的
Nanoporeアンプリコンリードを類似度に基づきクラスタリングし、非冗長なコンセンサス配列セットを構築する。

## 2. コアロジック
1.  **サンプリングと長さ選別**:
    * 感度 0.1% を担保するため、解析リード数を `10,000` (`max_reads`) に制限する。
    * サンプルごとに指定された `min_length`, `max_length` を用いてノイズリードを除去する。
2.  **非冗長なコンセンサス回収**:
    * `Amplicon_Sorter` が出力する個別のクラスターファイル（`*_n_consensussequences.fasta`）を収集する。
    * ツールが自動生成する統合ファイル（`*_consensussequences.fasta`）をあえて使用せず、手動でマージした後に `seqkit rmdup -s` を実行することで、配列ベースでの厳密な重複除去（De-duplication）を行う。
3.  **結果の整理**:
    * 最終的に重複を除去した `${sample_id}_clustered_consensus.fasta` を後続の解析（アノテーション等）の入力とする。
    * 全ての出力ファイルを `all_outputs/` に保持し、詳細なクラスタリング結果（`results.txt`）をサマリーとして抽出する。

## 3. 使用ツールとコンテナ
* **Amplicon_Sorter**: v2024.05.28
* **Seqkit**: (Container内蔵)
* **Container**: `docker://quay.io/biocontainers/amplicon_sorter:2024.05.28--pyhdfd78af_1`

## 4. 実装コード (Nextflow)

```nextflow
/*
 * AMPLICON_SORTER Workflow Components
 * 1. AS_PRE_STATS: Count reads before clustering (using SeqKit)
 * 2. AMPLICON_SORTER: Core clustering (using Local SIF + Python Script)
 * 3. AS_DEDUPLICATE: Merge and remove duplicate OTUs (using SeqKit)
 */

/*
 * Step 1: Pre-processing stats to count total reads used
 */
process AS_PRE_STATS {
    tag "${sample_id}"
    label 'process_low'
    container 'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("read_count_actual.txt"), emit: count
    path "versions_as_pre.yml",                          emit: versions

    script:
    """
    # Count actual reads (denominator for downstream analysis)
    seqkit stats -T ${reads} | tail -n 1 | awk '{print \$4}' > read_count_actual.txt
    
    echo "SeqKit (AS_PRE): \$(seqkit version | awk '{print \$2}')" > versions_as_pre.yml
    """
}

/*
 * Step 2: Core clustering process
 */
process AMPLICON_SORTER {
    tag "${sample_id}"
    label 'process_high'
    // Local SIF container with edlib, biopython, matplotlib, numpy
    container "${params.as_container_path}"

    input:
    tuple val(sample_id), path(reads), val(min_len), val(max_len)
    val max_reads
    path as_script  // Path to local amplicon_sorter.py (2025-10-09)

    output:
    tuple val(sample_id), path("out_dir/*_consensussequences.fasta"), emit: raw_consensus, optional: true
    path "out_dir/results.txt",                                      emit: summary,       optional: true
    path "out_dir/*",                                                emit: all_outputs
    path "versions_as_main.yml",                                     emit: versions

    script:
    """
    # Execute the local clustering script
    python3 ${as_script} \\
        -i ${reads} \\
        -n ${max_reads} \\
        -min ${min_len} \\
        -max ${max_len} \\
        -o out_dir

    echo "Amplicon_Sorter_Script: 2025-10-09 (Local)" > versions_as_main.yml
    """
}

/*
 * Step 3: Post-processing to merge and deduplicate OTUs
 */
process AS_DEDUPLICATE {
    tag "${sample_id}"
    label 'process_low'
    publishDir "${params.outdir}/03_amplicon_sorter/${sample_id}", mode: 'copy'
    container 'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0'

    input:
    tuple val(sample_id), path(raw_fastas)

    output:
    tuple val(sample_id), path("${sample_id}_clustered_consensus.fasta"), emit: consensus
    path "versions_as_post.yml",                                         emit: versions

    script:
    """
    # Merge all cluster consensus files and remove exact duplicates
    if [ -n "${raw_fastas}" ]; then
        cat ${raw_fastas} > merged_tmp.fasta
        seqkit rmdup -s merged_tmp.fasta -o ${sample_id}_clustered_consensus.fasta
    else
        touch ${sample_id}_clustered_consensus.fasta
    fi

    echo "SeqKit (AS_POST): \$(seqkit version | awk '{print \$2}')" > versions_as_post.yml
    """
}

```



# プロセス要件定義: 4.OTU_COUNT_TABLE
- module/otu_count.nf
## 1. 目的
`Amplicon_Sorter` の出力（FASTAヘッダー）からリード数を集計し、クラスタリングされなかった分（Unsorted）を含むカウントテーブルを作成する。

## 2. コアロジック
1.  **外部ライブラリへの依存排除**:
    * Singularity環境での可搬性を高めるため、Biopython等の外部モジュールを使用せず、Python標準ライブラリのみでFASTAヘッダーをパースする。
2.  **カウント抽出**:
    * ヘッダー末尾の `(N)` 形式から各クラスターのリード数を抽出する。
3.  **未分類リード（Unsorted）の自動算出**:
    * `actual_divisor`（サンプリング後の総入力リード数）から全クラスターの合計値を差し引き、どの群にも属さなかったリード数を算出・追加する。
4.  **出力フォーマット**:
    * `sample_id`, `cluster_id`, `read_count` の3カラムからなる TSV 形式。

## 3. 実装上の注意点
* **Singularity互換性**: 標準ライブラリのみを使用することで、`docker://` 経由で取得したイメージが Singularity で実行される際も、追加のレイヤー構築なしで安定動作する。
* **データ整合性**: `AMPLICON_SORTER` から引き継いだ「実際の分母」を用いることで、計算の正確性を担保する。

## 4. 使用ツールとコンテナ
* **Python**: 3.10.x (Standard Library only)
* **Container**: `docker://quay.io/biocontainers/python:3.10.12`

## 5. 実装コード (Nextflow)

```nextflow
/*
 * OTU_COUNT Module
 * * Description:
 * Processes the OTU count matrix and incorporates total read counts per sample
 * to provide a normalized or summarized view of the community composition.
 * * Requirements:
 * - Container with Biopython and Python 3.11+
 */

process OTU_COUNT_TABLE {
    label 'process_medium'
    publishDir "${params.outdir}/04_otu_merge", mode: 'copy'

    // Using the verified Biopython container (Python 3.11 based)
    container 'https://depot.galaxyproject.org/singularity/biopython:1.79--py311h1425ee9_1'

    input:
    path count_matrix      // Generated by OTU_MERGE
    path divisor_files    // List of read_count_actual.txt from AMPLICON_SORTER

    output:
    path "otu_table_final.tsv",     emit: final_table
    path "otu_table_summary.txt",   emit: summary
    path "versions_otu_count.yml",  emit: versions

    script:
    """
    # Use python script to integrate actual read counts and finalize the matrix
    python3 << 'EOF'
import os
import glob

# 1. Load sample divisors (total reads used in Amplicon_Sorter)
sample_totals = {}
divisor_files = "${divisor_files}".split()
for df in divisor_files:
    # Assuming directory structure or naming links to sample_id
    # If the file is just 'read_count_actual.txt', we infer sample from path
    sample_id = os.path.basename(os.path.dirname(df))
    with open(df, 'r') as fh:
        count = int(fh.read().strip())
        sample_totals[sample_id] = count

# 2. Process count matrix and add summary info
with open("${count_matrix}", 'r') as f_in, open("otu_table_final.tsv", 'w') as f_out:
    header = f_in.readline().strip().split('\\t')
    # Add a row or header for total counts if needed, or just pass through
    f_out.write("\\t".join(header) + "\\n")
    
    for line in f_in:
        f_out.write(line)

# 3. Generate a simple summary for provenance
with open("otu_table_summary.txt", 'w') as summary:
    summary.write("Sample\\tTotal_Reads_Used\\n")
    for s, val in sorted(sample_totals.items()):
        summary.write(f"{s}\\t{val}\\n")

EOF

    # Capture versions
    echo "Python (Biopython container): \$(python3 --version)" > versions_otu_count.yml
    echo "Biopython: 1.79" >> versions_otu_count.yml
    """
}

```






# プロセス要件定義: 5.OTU_MERGE
- module/otu_merge.nf
## 1. 目的
全サンプルから得られたコンセンサス配列を集約し、配列が完全に一致（100% Identity）するものを単一の OTU (Operational Taxonomic Unit) として統合する。これにより、サンプルを跨いだ比較を可能にし、後続の BLAST アノテーションの重複計算を排除する。

## 2. コアロジック
1.  **配列ベースの集約**:
    * 各サンプルの `consensus.fasta` から配列文字列を抽出し、ハッシュマップ（辞書型）のキーとして管理する。
    * 全く同じ塩基配列を持つクラスターを、サンプルを問わず一つの OTU ID（例: OTU_0001）に統合する。
2.  **カウントマトリックスの作成**:
    * 行を OTU、列をサンプル名とするマトリックスを生成する。
    * 各セルには、そのサンプルにおいて該当配列（OTU）が何リード存在したかの値を格納する。
3.  **代表配列の出力**:
    * 統合された各 OTU に対して一つの代表配列を選出し、BLAST 入力用の単一 FASTA ファイルを出力する。

## 3. 実装上の注意点
* **メモリ効率**: 全サンプルの配列をメモリ上に保持するため、サンプル数や OTU 数が膨大な場合は、Python のジェネレータや外部ツール（vsearch 等）への切り替えを検討する。今回は Python 標準ライブラリと Pandas での実装を想定。
* **IDの不変性**: `OTU_0001` などの ID は、このプロセス以降、すべての統計・アノテーション結果と紐付く一貫した識別子となる。

## 4. 使用ツールとコンテナ
* **Python**: 3.10.x
* **Container**: `docker://quay.io/biocontainers/python:3.10.12`

## 5. 実装コード (Nextflow)

```nextflow
/*
 * OTU_MERGE Module
 * * Description:
 * Integrates clustered consensus sequences (OTUs) from all samples,
 * performs 100% identity clustering to identify unique global OTUs,
 * and generates a count matrix (OTU table).
 * * Requirements:
 * - Container with Biopython and Python 3.11+
 */

process OTU_MERGE {
    label 'process_medium'
    publishDir "${params.outdir}/04_otu_merge", mode: 'copy'

    // Using the verified Biopython container (Python 3.11 based)
    container 'https://depot.galaxyproject.org/singularity/biopython:1.79--py311h1425ee9_1'

    input:
    tuple val(sample_ids), path(consensus_fastas) // List of all consensus files
    tuple val(sample_ids_count), path(count_files)   // List of read_count_actual.txt files

    output:
    path "integrated_unique_otus.fasta", emit: otu_fasta
    path "otu_count_matrix.tsv",         emit: count_matrix
    path "versions_otu_merge.yml",       emit: versions

    script:
    """
    # 1. Integrate all sequences into a single global FASTA
    # Use python script to parse fasta and track sample origins
    python3 << 'EOF'
import os
from Bio import SeqIO
from collections import defaultdict

otu_counts = defaultdict(lambda: defaultdict(int))
unique_seqs = {}

# Iterate through each sample's consensus file
fastas = "${consensus_fastas}".split()
for f in fastas:
    sample_id = os.path.basename(f).replace("_clustered_consensus.fasta", "")
    for record in SeqIO.parse(f, "fasta"):
        seq = str(record.seq).upper()
        # Use sequence itself as key for 100% identity merging
        if seq not in unique_seqs:
            otu_id = f"OTU_{len(unique_seqs) + 1}"
            unique_seqs[seq] = otu_id
        
        target_otu = unique_seqs[seq]
        # In this logic, each cluster from Amplicon_Sorter is treated as 1 unit or weighted by its size if available
        # Here we increment by 1 for the existence of the cluster in the sample
        otu_counts[target_otu][sample_id] += 1

# Write integrated unique OTUs to FASTA
with open("integrated_unique_otus.fasta", "w") as fa:
    for seq, otu_id in unique_seqs.items():
        fa.write(f">{otu_id}\\n{seq}\\n")

# Write OTU count matrix (TSV)
samples = sorted("${sample_ids}".split())
with open("otu_count_matrix.tsv", "w") as tsv:
    tsv.write("OTU_ID\\t" + "\\t".join(samples) + "\\n")
    for otu_id in sorted(unique_seqs.values(), key=lambda x: int(x.split('_')[1])):
        row = [otu_id]
        for s in samples:
            row.append(str(otu_counts[otu_id][s]))
        tsv.write("\\t".join(row) + "\\n")
EOF

    # 2. Capture versions for provenance
    echo "Python (Biopython container): \$(python3 --version)" > versions_otu_merge.yml
    echo "Biopython: 1.79" >> versions_otu_merge.yml
    """
}

```


# プロセス要件定義: 6.BLAST_ANNOTATE (統合アノテーション版)
- module/blast_annotate.nf
## 1. 目的
`OTU_MERGE` によって全サンプルから集約・抽出された「ユニークな代表配列（OTU）」に対し、一括で系統同定（アノテーション）を行う。サンプルを跨いで重複する配列への検索を一度に集約することで、計算リソースと解析時間を劇的に削減する。

## 2. コアロジック
1.  **一括検索 (Batch Search)**:
    * 個別のサンプルごとではなく、パイプライン全体で生成された単一の代表配列ファイル（`unique_otus.fasta`）をクエリとして実行する。
2.  **厳格なフィルタリングパラメータ**:
    * **`-max_target_seqs 1`**: 各 OTU に対して最良のヒット（Best Hit）のみを報告。
    * **`-qcov_hsp_perc 70`**: クエリの全体長に対して70%以上のアライメント範囲を要求。
    * **`-evalue 1e-10`**: 偽陽性を抑える厳格な有意性閾値。
3.  **マッピングキーの維持**:
    * 出力の第1カラム（`qseqid`）には `OTU_0001` 形式の ID が保持される。これにより、`merged_otu_counts.tsv` との完全な紐付けが可能となる。

## 3. 実装上の注意点
* **実行効率**: サンプル数 $N$ 回ではなく $1$ 回の実行となるため、スレッド数（`-num_threads`）を最大限活用できるようリソースを割り当てる。
* **データベースバインド**: ローカルデータベースのディレクトリ（`db_dir`）を Singularity の実行時オプションで適切にマウント（bind）する必要がある。

## 4. 使用ツールとコンテナ
* **BLAST+**: v2.15.0
* **Container**: `docker://quay.io/biocontainers/blast:2.15.0--pl5321h6f7f691_1`

## 5. 実装コード (Nextflow)

```nextflow
/*
 * BLAST_ANNOTATE: 
 * Perform batch BLAST search for the integrated unique OTUs.
 * Minimizes redundant computations by annotating each unique sequence only once.
 * Software: BLAST+ v2.16.0
 */
process BLAST_ANNOTATE {
    label 'process_high'
    publishDir "${params.outdir}/05_annotation", mode: 'copy'

    // Updated to the verified stable URI
    container 'https://depot.galaxyproject.org/singularity/blast:2.16.0--h66d330f_4'

    input:
    path unique_otus_fasta  // Integrated FASTA file from OTU_MERGE
    path db_dir             // Directory containing BLAST database (mounted via nextflow.config or --blast_db_dir)
    val db_name             // Name of the database (e.g., 'nt')
    val blast_type          // BLAST algorithm (e.g., 'blastn')

    output:
    path "all_otus_blast_results.tsv", emit: blast_results
    path "versions_blast.yml",         emit: versions

    script:
    """
    # Execute batch BLAST search
    # Params: min 70% query coverage, max 1 target, evalue 1e-10
    ${blast_type} \\
        -query ${unique_otus_fasta} \\
        -db ${db_dir}/${db_name} \\
        -max_target_seqs 1 \\
        -qcov_hsp_perc 70 \\
        -evalue 1e-10 \\
        -num_threads ${task.cpus} \\
        -outfmt "6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids" \\
        -out all_otus_blast_results.tsv

    # Capture BLAST version for provenance
    echo "BLAST+: \$(${blast_type} -version | head -n 1 | awk '{print \$2}')" > versions_blast.yml
    """
}

```

# プロセス要件定義: 7.CREATE_BIOM (統合データ & 例外処理版)
- module/create_biom.nf
## 1. 目的
`OTU_MERGE` で作成された全サンプル集計済みのカウントテーブル（`Unsorted` 行を含む）と、`BLAST_ANNOTATE` による代表配列のアノテーション結果を統合する。配列が存在しない `Unsorted` リードについても、統計上の欠損が生じないよう適切に例外処理を行い、単一の BIOM ファイルおよび統合レポート（TSV）を出力する。

## 2. コアロジック
1.  **系統情報の解決 (TaxonKit)**:
    * BLAST ヒットがあった TaxID について、NCBI Taxonomy データベースを参照し 7 階層の系統情報を生成する。
2.  **データ統合 (Left Join & Exception Handling)**:
    * `merged_otu_counts.tsv` をベースに BLAST 結果を結合する。
    * **例外処理**: `Unsorted` 行、および BLAST でヒットがなかった OTU に対しては、系統情報を一律で `Unclassified;Unclassified;...` と補完する。これにより、BIOM ファイルの構造的整合性を保つ。
3.  **BIOM ファイルの生成**:
    * すべてのサンプル列、OTU ID（`OTU_xxxx` および `Unsorted`）、および補完済みの系統情報を一つの BIOM v1 (JSON) ファイルに集約する。

## 3. 実装上の注意点
* **配列の不在への対応**: `Unsorted` 行は配列データを持たないため、BLAST 結果には含まれない。プログラム側でこの欠損を「エラー」ではなく「未分類（Unclassified）」として定義することを徹底する。
* **数値の完全性**: サンプルごとの総リード数が、フィルタリング後のリード総数と一致することを確認する。

## 4. 使用ツールとコンテナ
* **TaxonKit**: v0.14.x
* **Python**: 3.x (pandas, biom-format)
* **Container**: `docker://quay.io/biocontainers/mulled-v2-4299ec503e944d18721c5f8df654f5c9071c3c90:218408f615364e0306e9b4187e14d1867e91409f-0`

## 5. 実装コード (Nextflow)

```nextflow
/*
 * CREATE_BIOM Module: Separated into TaxonKit and Python(Pandas/Biom)
 */

process TAXONKIT_LINEAGE {
    label 'process_low'
    container 'https://depot.galaxyproject.org/singularity/taxonkit:0.18.0--h9ee0642_0'

    input:
    path blast_results
    path tax_db_dir

    output:
    path "reformatted_lineage.txt", emit: lineage
    path "versions_taxonkit.yml",  emit: versions

    script:
    """
    export TAXONKIT_DB="${tax_db_dir}"
    awk -F'\\t' '{print \$NF}' ${blast_results} | \
        taxonkit lineage | taxonkit reformat -f "{k};{p};{c};{o};{f};{g};{s}" > reformatted_lineage.txt
    echo "TaxonKit: \$(taxonkit version)" > versions_taxonkit.yml
    """
}

process BIOM_GENERATE {
    label 'process_medium'
    publishDir "${params.outdir}/01_Final_Outputs", mode: 'copy'
    container 'https://depot.galaxyproject.org/singularity/biom-format:2.1.15'
    // Note: Biom-format image should include pandas. If not, quay.io/biocontainers/pandas:2.2.1 can be used with a combined script.

    input:
    path merged_counts
    path lineage
    path blast_results

    output:
    path "merged_results.biom",   emit: biom
    path "merged_otu_report.tsv", emit: tsv_report

    script:
    """
    # Python script for final matrix joining and BIOM export
    """
}

```

# プロセス要件定義: 8. FINALIZE_RESULTS
- module/finalize_results.nf
## 1. 目的
パイプラインの各工程で生成された主要な成果物を、整理されたディレクトリ構造で最終出力フォルダ（`results`）に配置する。ユーザーが解析結果を俯瞰しやすくし、可視化ツールや報告書作成への移行をスムーズにする。

## 2. 収集対象ファイルと整理構造
最終出力ディレクトリ（`${params.outdir}`）の下に以下の構造で整理します。

1.  **`01_Final_Analysis/`**: 
    * `merged_results.biom`: 全サンプル統合済みの可視化用ファイル。
    * `merged_otu_report.tsv`: アノテーション情報付きの統合カウントテーブル。
2.  **`02_Sequences/`**:
    * `unique_otus.fasta`: OTU IDが付与されたユニークな代表配列（BLASTに使用したもの）。
    * `all_samples_consensus.fasta`: 全サンプルの全クラスターコンセンサス配列。
3.  **`03_Statistics/`**:
    * `read_tracking_summary.tsv`: 各工程でのリード数推移。
    * `phylum_distribution_summary.tsv`: 門レベルの組成集計表。
    * `*.png`: 組成確認用の簡易グラフ。

## 3. 実装上の注意点
* **シンボリックリンクの回避**: `publishDir` の `mode: 'copy'` を使用し、中間ファイルが削除されても成果物が残るようにする。
* **一括集約**: 各プロセスから出力される `path` を `collect()` して受け取り、一回のプロセスで整理を行う。

## 4. 実装コード (Nextflow)

```nextflow
/*
 * FINALIZE_RESULTS:
 * Consolidates all key outputs into a structured directory for the user.
 */
process FINALIZE_RESULTS {
    label 'process_low'
    publishDir "\${params.outdir}", mode: 'copy'

    input:
    path merged_biom         // From CREATE_BIOM
    path merged_report       // From CREATE_BIOM
    path unique_otu_fasta    // From OTU_MERGE
    path all_consensus_fasta // From AMPLICON_SORTER (collected)
    path stats_tables        // From SUMMARY_REPORT (collected)
    path plots               // From SUMMARY_REPORT (collected)

    output:
    path "01_Final_Analysis/*"
    path "02_Sequences/*"
    path "03_Statistics/*"

    script:
    """
    # Create directory structure
    mkdir -p 01_Final_Analysis 02_Sequences 03_Statistics

    # 1. Final Tables & BIOM
    cp \${merged_biom} 01_Final_Analysis/
    cp \${merged_report} 01_Final_Analysis/

    # 2. Sequences
    cp \${unique_otu_fasta} 02_Sequences/
    # 統合コンセンサスを作成（全サンプルのfastaを一つにまとめる）
    cat \${all_consensus_fasta} > 02_Sequences/all_samples_consensus.fasta

    # 3. Statistics & Plots
    cp \${stats_tables} 03_Statistics/
    cp \${plots} 03_Statistics/
    """
}
```


# プロセス要件定義: 9. GENERATE_PROVENANCE (高信頼性版)
- generate_provenance.nf
## 1. 目的
実際に解析に使用された各コンテナー内から取得したバージョン情報と、並列・統合の構造を明示したパイプライン図を統合し、完全な解析証跡レポートを作成する。

## 2. コアロジック
1.  **バージョン集約**: 各プロセス（fastp, blast等）から出力された `versions.yml` を `collect` して読み込み、重複を排除してリスト化する。
2.  **並列化構造の可視化**: Mermaid図において、サンプルごとの並列処理セクションと、全サンプル統合セクションを視覚的に分離する。

## 3. 実装コード (Nextflow)

```nextflow

/*
 * GENERATE_PROVENANCE Module
 * * Description:
 * Aggregates version information from all processes and creates a 
 * Markdown report and a plain text file for audit trails.
 * * Requirements:
 * - Container with Python 3.11+
 */

process GENERATE_PROVENANCE {
    label 'process_low'
    publishDir "${params.outdir}/00_provenance", mode: 'copy'

    // Using the verified Biopython container (Python 3.11 based)
    container 'https://depot.galaxyproject.org/singularity/biopython:1.79--py311h1425ee9_1'

    input:
    path v_files // List of all version files (versions_*.yml) collected via .collect()
    val workflow_info // workflow object from main.nf

    output:
    path "software_versions.md",  emit: md_report
    path "software_versions.txt", emit: txt_report
    path "provenance.yml",        emit: yaml_summary

    script:
    """
    python3 << 'EOF'
    import datetime

    # 1. Collect versions from all input files
    version_data = {}
    files = "${v_files}".split()
    for f in files:
        with open(f, 'r') as fh:
            for line in fh:
                if ':' in line:
                    key, val = line.strip().split(':', 1)
                    version_data[key.strip()] = val.strip()

    # 2. Generate Markdown Report
    with open("software_versions.md", "w") as md:
        md.write("# Software Provenance Report\\n\\n")
        md.write(f"**Analysis Date:** {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\\n")
        md.write(f"**Nextflow Version:** ${workflow_info.nextflow.version}\\n")
        md.write(f"**Command Line:** `${workflow_info.commandLine}`\\n\\n")
        md.write("| Software/Tool | Version |\\n")
        md.write("| :--- | :--- |\\n")
        for tool, ver in sorted(version_data.items()):
            md.write(f"| {tool} | {ver} |\\n")

    # 3. Generate plain text and yaml summary
    with open("software_versions.txt", "w") as txt:
        for tool, ver in sorted(version_data.items()):
            txt.write(f"{tool}: {ver}\\n")

    with open("provenance.yml", "w") as yml:
        yml.write("software_versions:\\n")
        for tool, ver in sorted(version_data.items()):
            yml.write(f"  {tool}: {ver}\\n")
    EOF
    """
}

```

# パイプライン統合定義: main.nf
- **役割**: 全 11 モジュールのオーケストレーション、データの並列分散および集約制御。

## 1. ワークフローの設計思想
1.  **I/O 最適化**: 同一ディレクトリ内の複数サンプル（バーコード違い）に対し、重い前処理（`Porechop`）を一度だけ実行し、その後に各サンプルへリードを分配する設計。
2.  **型安全なパース**: `sampleSheet.csv` から数値型や文字列型を適切にキャストして読み込み、下流プロセスでのエラーを未然に防ぐ。
3.  **完全なトレーサビリティ**: 全プロセスの `versions.yml` を収集し、解析終了時に自動で証跡レポートを作成する。

## 2. 実装コード (Nextflow)

```nextflow
nextflow.enable.dsl = 2

/*
 * 1. Import Refactored Modules
 */
include { PORECHOP_TRIM as PREPROCESS_ADAPTERS; SEQKIT_CLEAN } from './modules/preprocess.nf'
include { DEMULTIPLEX                                       } from './modules/demux.nf'
include { KRAKEN2_CLASSIFY; KRAKENTOOLS_EXTRACT             } from './modules/kraken2_filter.nf'
include { AS_PRE_STATS; AMPLICON_SORTER; AS_DEDUPLICATE      } from './modules/amplicon_sorter.nf'
include { OTU_COUNT_TABLE                                   } from './modules/otu_count.nf'
include { OTU_MERGE                                         } from './modules/otu_merge.nf'
include { BLAST_ANNOTATE                                    } from './modules/blast_annotate.nf'
include { FINALIZE_RESULTS                                  } from './modules/finalize_results.nf'
include { GENERATE_PROVENANCE                               } from './modules/generate_provenance.nf'

workflow {

    // --- 2. Input Channel Creation ---
    
    // Parse sampleSheet.csv and create a metadata map for each sample
    ch_sample_sheet = Channel.fromPath(params.sample_sheet)
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                sample_id:  row.sample,     // Updated to match your CSV column 'sample'
                fastq_dir:  row.fastq_dir,
                min_len:    row.min_len.toInteger(),
                max_len:    row.max_len.toInteger(),
                f_idx:      row.fwd_index ?: "",
                f_prm:      row.fwd_primer,
                r_idx:      row.rev_index ?: "",
                r_prm:      row.rev_primer
            ]
            return meta
        }

    // Grouping by fastq_dir to run PREPROCESS only once per directory
    ch_for_preprocess = ch_sample_sheet
        .map { meta -> [ meta.fastq_dir, file("${meta.fastq_dir}/*.{fastq,fastq.gz,fq,fq.gz}") ] }
        .groupTuple()
        .map { dir, files -> [ dir.split('/')[-1], files.flatten() ] }

    // --- 3. Parallel Processing Phase (Per-Directory / Per-Sample) ---

    // [A] Pre-processing (Adapter trimming via Porechop)
    PREPROCESS_ADAPTERS(ch_for_preprocess)

    // [B] Redistribute trimmed reads back to individual samples
    ch_demux_input = ch_sample_sheet
        .map { meta -> [ meta.fastq_dir.split('/')[-1], meta ] }
        .combine(PREPROCESS_ADAPTERS.out.chopped_fastq, by: 0)
        .map { dir_id, meta, fastq ->
            [ meta.sample_id, fastq, meta.min_len, meta.max_len, meta.f_idx, meta.f_prm, meta.r_idx, meta.r_prm ]
        }

    // [C] High-precision Demultiplexing
    DEMULTIPLEX(ch_demux_input)

    // Length Filtering (SeqKit Clean) - Using metadata from CSV via DEMULTIPLEX out
    SEQKIT_CLEAN(DEMULTIPLEX.out.reads)

    // [D] Host Removal / Target Extraction (Kraken2 Filter Split)
    KRAKEN2_CLASSIFY(
        SEQKIT_CLEAN.out.reads,
        params.kraken2_db
    )
    
    ch_for_extraction = KRAKEN2_CLASSIFY.out.output
        .join(KRAKEN2_CLASSIFY.out.report)
        .join(SEQKIT_CLEAN.out.reads)

    KRAKENTOOLS_EXTRACT(
        ch_for_extraction,
        params.target_taxid
    )

    // [E] Clustering & Variant Identification (Amplicon_Sorter Split)
    // E-1: Count reads before sorting (Denominator for count table)
    AS_PRE_STATS(KRAKENTOOLS_EXTRACT.out.filtered_fastq)

    // E-2: Core Clustering (Requires metadata re-join for min/max lengths)
    ch_as_input = KRAKENTOOLS_EXTRACT.out.filtered_fastq
        .join(ch_sample_sheet.map { m -> [m.sample_id, m.min_len, m.max_len] })
    
    AMPLICON_SORTER(
        ch_as_input,
        params.max_reads,
        file(params.as_script_path)
    )

    // E-3: Post-clustering Deduplication (Using SeqKit)
    AS_DEDUPLICATE(AMPLICON_SORTER.out.raw_consensus)

    // [F] Per-sample Read Counting
    OTU_COUNT_TABLE(
        AS_DEDUPLICATE.out.consensus,
        AS_PRE_STATS.out.count
    )

    // --- 4. Global Consolidation Phase ---

    ch_all_consensus = AS_DEDUPLICATE.out.consensus.map { it[1] }.collect()
    ch_all_counts    = OTU_COUNT_TABLE.out.final_table.collect()

    // [G] Merge OTUs across all samples (100% identity)
    OTU_MERGE(ch_all_consensus, ch_all_counts)

    // [H] Batch BLAST Annotation
    BLAST_ANNOTATE(
        OTU_MERGE.out.otu_fasta,
        params.blast_db_dir,
        params.blast_db_name,
        params.blast_type
    )

    // --- 5. Result Finalization & Provenance ---

    // Organize final results
    FINALIZE_RESULTS(
        OTU_MERGE.out.count_matrix,
        OTU_MERGE.out.otu_fasta,
        ch_all_consensus,
        BLAST_ANNOTATE.out.blast_results
    )

    // Aggregate version info from all utilized containers
    ch_versions = PREPROCESS_ADAPTERS.out.versions
        .mix(SEQKIT_CLEAN.out.versions)
        .mix(DEMULTIPLEX.out.versions)
        .mix(KRAKEN2_CLASSIFY.out.versions)
        .mix(KRAKENTOOLS_EXTRACT.out.versions)
        .mix(AS_PRE_STATS.out.versions)
        .mix(AMPLICON_SORTER.out.versions)
        .mix(AS_DEDUPLICATE.out.versions)
        .mix(OTU_COUNT_TABLE.out.versions)
        .mix(OTU_MERGE.out.versions)
        .mix(BLAST_ANNOTATE.out.versions)
        .collect()

    GENERATE_PROVENANCE(ch_versions, workflow)
}

```






