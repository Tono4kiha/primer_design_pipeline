

# primer_design_pipeline

**English** | **中文**

A K-mer based primer design pipeline for any target genomic region.  
It enumerates candidate primers (18–27 nt) within the target region plus flanking sequences, filters and scores them, and performs genome-wide specificity evaluation using BLAST.

基于 K-mer 策略的引物设计流程，适用于任意目标基因组区域。  
自动枚举目标区域及侧翼序列中的候选引物（18–27 nt），进行过滤、打分，并结合 BLAST 进行全基因组特异性评估。

---

## Features | 功能特性

- Enumerate candidate primers (F/R) as K-mers  
  枚举候选引物（F/R）
- Filter by GC, Tm, homopolymer, low complexity, repeats, hairpin, homodimer, 3′ self-interaction  
  按 GC、Tm、均聚物、低复杂度、重复、发卡、同源二聚体、3′ 自相互作用过滤
- Pair F/R candidates based on genomic coordinates, with amplicon size and interaction filters  
  基于基因组坐标配对 F/R，并过滤扩增子大小及二聚体相互作用
- Pre-BLAST scoring and selection of top pairs  
  BLAST 前打分并筛选最优配对
- Parse genome BLAST results to classify specificity  
  解析全基因组 BLAST 结果并分类特异性
- Output: UNIQUE_PRODUCT / MULTIPLE_PRODUCTS / NO_PRODUCT_IN_REFERENCE  
  输出：唯一产物 / 多重产物 / 参考基因组中无产物

---

## Workflow | 工作流程

| Step | Script | Description 描述 |
|------|--------|------------------|
| 1 | `01_generate_candidates.py` | Enumerate and filter candidate primers 枚举并过滤候选引物 |
| 2 | `02_pair_primers.py` | Pair F/R candidates, filter amplicon and interactions 配对并过滤 |
| 3 | `03_prepare_blast.py` | Select top pairs and generate BLAST input files 选择最优配对并生成 BLAST 输入 |
| 4 | `04_parse_blast.py` | Parse BLAST hits and classify specificity 解析 BLAST 结果并分类 |

---

## Requirements | 环境依赖

- Linux (recommended) | Linux（推荐）
- Python 3.x
- Python package: `primer3-py==2.3.0`
- BLAST+ (`blastn`, `makeblastdb`)

---

## Installation | 安装

请参照官方文档安装所需依赖：

- Python: <https://www.python.org/downloads/>
- primer3-py: <https://primer3-py.readthedocs.io/>
- BLAST+: <https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download>

> **Note:** 安装后请确保 `blastn` 和 `makeblastdb` 已正确加入 `PATH`。

---

## Input Files | 输入文件

### 1. `input/target_region.fa` (required | 必需)

FASTA 文件，只包含一个序列，由以下部分组成：

- 目标区域上游侧翼序列（如 500 bp）
- 目标基因本体（包括 exon / intron）
- 目标区域下游侧翼序列（如 500 bp）

示例：

```fasta
>target_region
ACGTACGTACGT......
```

### 2. `input/regions.tsv` (optional | 可选)

用于标注引物所在区域（exon / intron / flank）。  
坐标格式：**0-based, half-open [start, end)**

示例：

```text
label    start   end
upstream        0       500
exon1   500     650
intron1 650     1200
exon2   1200    1350
intron2 1350    1800
exon3   1800    1950
downstream      1950    2450
```

---

## Usage | 使用方法

```bash
# 1. 生成候选引物
python 01_generate_candidates.py \
    --fasta input/target_region.fa \
    --regions input/regions.tsv \
    --outdir out

# 2. 配对引物
python 02_pair_primers.py \
    --fasta input/target_region.fa \
    --candidates out/candidates_pass.tsv \
    --out out/pairs_pass.tsv

# 3. 选择 top 配对并生成 BLAST 输入
python 03_prepare_blast.py \
    --fasta input/target_region.fa \
    --pairs out/pairs_pass.tsv \
    --top 200 \
    --outdir out/blast

# 4. 构建参考基因组 BLAST 数据库（示例）
makeblastdb \
    -in reference_genome.fa \
    -dbtype nucl \
    -parse_seqids \
    -out ref_db

# 5. 对引物序列进行 BLAST（示例）
blastn \
    -task blastn-short \
    -query out/blast/blast_primers.fa \
    -db ref_db \
    -dust no \
    -soft_masking false \
    -evalue 1000 \
    -max_target_seqs 10000 \
    -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qseq sseq' \
    -out out/blast/primer_hits.tsv

# 6. 解析 BLAST 结果
python 04_parse_blast.py \
    --pairs out/blast/top_pairs.tsv \
    --blast out/blast/primer_hits.tsv \
    --out out/blast/pairs_with_blast.tsv
```

---

## Configuration | 参数配置

所有参数均在 `settings.py` 中修改，关键参数如下：

| Parameter | Value | Description |
|-----------|-------|-------------|
| `PRIMER_MIN_LEN` / `PRIMER_MAX_LEN` | 18 / 27 | 引物长度范围 |
| `PRIMER_MIN_TM` / `PRIMER_OPT_TM` / `PRIMER_MAX_TM` | 57.0 / 60.0 / 63.0 | 引物 Tm 范围及最优值 |
| `PRIMER_MIN_GC` / `PRIMER_MAX_GC` | 40.0 / 60.0 | 引物 GC 含量范围 |
| `PRODUCT_MIN` / `PRODUCT_OPT` / `PRODUCT_MAX` | 50 / 250 / 700 | 扩增子长度范围及最优值 |
| `MAX_DELTA_TM` | 2.5 | F/R Tm 最大差值 |
| `MAX_HAIRPIN_TM` | 47.0 | 发卡结构最大 Tm |
| `MAX_HOMODIMER_TM` | 47.0 | 同源二聚体最大 Tm |
| `MAX_SELF_END_TM` | 47.0 | 3′ 自相互作用最大 Tm |
| `MAX_HETERODIMER_TM` | 47.0 | F/R 异源二聚体最大 Tm |
| `MAX_PAIR_END_TM` | 47.0 | F/R 3′ 相互作用最大 Tm |
| `MV_CONC` / `DV_CONC` / `DNTP_CONC` / `DNA_CONC` | 50.0 / 1.5 / 0.6 / 50.0 | 热力学计算用离子浓度、dNTP 浓度、DNA 浓度 |

> 根据实际 PCR 反应体系调整 `MV_CONC`、`DV_CONC`、`DNTP_CONC`、`DNA_CONC`。

---

## Outputs | 输出文件

| File | Description 描述 |
|------|------------------|
| `out/candidates_all.tsv` | 所有候选引物及其过滤状态、失败原因 |
| `out/candidates_pass.tsv` | 通过所有过滤的候选引物 |
| `out/pairs_pass.tsv` | 通过配对的引物对及其评分 |
| `out/blast/top_pairs.tsv` | 用于 BLAST 的 top 引物对 |
| `out/blast/blast_primers.fa` | 引物序列 FASTA |
| `out/blast/blast_amplicons.fa` | 扩增子序列 FASTA |
| `out/blast/pairs_with_blast.tsv` | 最终分类结果（UNIQUE_PRODUCT / MULTIPLE_PRODUCTS / NO_PRODUCT_IN_REFERENCE） |

---

## Notes | 注意事项

- 推荐在 Linux 环境下运行。
- 请确保 BLAST+ 已正确安装且路径可用。
- 使用 `blastn-short` 模式处理短引物序列（<50 nt）。
- `UNIQUE_PRODUCT` 为最值得优先实验验证的引物对。
- `NO_PRODUCT_IN_REFERENCE` 可能是参考基因组未包含目标位点所致，不能简单视为引物无效。
