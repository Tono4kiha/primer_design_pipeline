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
- Re-score passed pairs with normalized penalty model  
  对通过硬过滤的引物对进行标准化罚分重打分
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
| 2b | `02b_rescore_pairs.py` | Re-score passed pairs with normalized penalties 对通过配对的引物对重新打分 |
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
- BLAST+: <https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/>

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
