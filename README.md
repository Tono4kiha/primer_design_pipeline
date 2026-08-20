
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

整个流程由 5 个 Python 脚本组成，按顺序执行：

### Step 1: `01_generate_candidates.py`

**功能**：  
对目标区域及侧翼序列枚举长度 18–27 nt 的 K-mer，分别作为正向 (F) 和反向 (R) 候选引物。  
对每个候选引物进行以下过滤：

- GC 含量
- Tm
- homopolymer (均聚物)
- 低复杂度
- 串联重复
- 发卡结构 (hairpin)
- 同源二聚体 (homodimer)
- 3′ 端自相互作用 (self 3′ interaction)

热力学计算直接调用 `primer3-py` 2.3.0 的 `calc_tm`、`calc_hairpin`、`calc_homodimer`、`calc_end_stability` API。

**输入**：
- `input/target_region.fa`（必需）
- `input/regions.tsv`（可选）

**输出**：
- `out/candidates_all.tsv`：所有候选引物，包含过滤状态及失败原因
- `out/candidates_pass.tsv`：通过所有过滤的候选引物

---

### Step 2: `02_pair_primers.py`

**功能**：  
根据基因组坐标将 F/R 候选引物配对，只考虑扩增子长度在 50–700 bp 的组合。  
然后过滤以下条件：

- ΔTm (F 与 R 的 Tm 差值)
- heterodimer (F/R 异源二聚体)
- F/R 之间的 3′ 端相互作用

注意：此处没有对所有 F/R 做完整笛卡尔积，而是先利用坐标限制扩增子大小，再计算相对昂贵的热力学参数，效率较高。

**输入**：
- `input/target_region.fa`
- `out/candidates_pass.tsv`

**输出**：
- `out/pairs_pass.tsv`：通过硬过滤的引物对，包含详细信息和初步打分 `score_preblast`

---

### Step 2b: `02b_rescore_pairs.py`

**功能**：  
对已通过硬过滤的引物对进行重新打分。  
将每个软性指标（如 Tm、GC、二级结构、产物大小等）转换为 0–1 之间的惩罚值（0=理想，1=接近极限），按权重求和得到 `score_preblast`。  
最终按分数升序排序，分数越低越好。  
本步骤不使用 BLAST，基因组特异性作为下游的硬性生物学标准单独处理。

**输入**：
- `out/pairs_pass.tsv`

**输出**：
- `out/pairs_rescored.tsv`：重新打分并排序后的引物对，新增各惩罚项、总分及排名

---

### Step 3: `03_prepare_blast.py`

**功能**：  
从大量合格引物对中选出 pre-BLAST score 最好的前 N 对（例如 200 对），生成用于 BLAST 的 FASTA 文件。

**输入**：
- `input/target_region.fa`
- `out/pairs_rescored.tsv`（或 `out/pairs_pass.tsv`，取决于是否执行了 Step 2b）
- `--top` 参数指定选取数量

**输出**：
- `out/blast/top_pairs.tsv`：选出的引物对
- `out/blast/blast_primers.fa`：引物序列 FASTA
- `out/blast/blast_amplicons.fa`：扩增子序列 FASTA

---

### Step 4: `04_parse_blast.py`

**功能**：  
读取 F/R 的全基因组 BLAST 结果。  
不是简单统计 BLAST hits 数量，而是根据 chromosome/scaffold、方向和距离，重新构建基因组上真正可能形成的 PCR 产物。  
最终将每对引物分类为：

- `UNIQUE_PRODUCT`：基因组中仅有一个预期产物（最值得优先验证）
- `MULTIPLE_PRODUCTS`：存在多个潜在产物
- `NO_PRODUCT_IN_REFERENCE`：参考基因组中未找到可形成的产物（可能因参考基因组缺少目标位点）

默认 BLAST 过滤条件：  
- 引物 query 覆盖率 ≥ 80%
- 最大错配数 = 3
- 最大 gap = 0
- 3′ 端窗口 = 5 nt，最大错配 = 1

**输入**：
- `out/blast/top_pairs.tsv`
- `out/blast/primer_hits.tsv`（BLAST 输出）https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html

**输出**：
- `out/blast/pairs_with_blast.tsv`：最终分类结果

---

## Requirements | 环境依赖

- Linux (recommended) | Linux（推荐）
- Python 3.x
- Python package: `primer3-py==2.3.0`
- BLAST+ (`blastn`, `makeblastdb`)

---

## Installation | 安装

请参照官方文档安装所需依赖：

- Python: <https://www.python.org/downloads/>| 中文名 | Name | E-mail | Blog / Other | Role |
- primer3-py: <https://pypi.org/project/primer3-py/>
- Blast <https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html/>

The script sources are complex, and we thank the authors of some scripts in the corresponding directories.
If there are any omissions, please inform us in the issue.

The main creators include:

<table>
  
  <tr>
    <td>中文名</td>
    <td>Name</td>
    <td>E-mail</td>
    <td>Blog / Other </td>
  </tr>
  
  <tr>
    <td>李硕</td>
    <td>Biols0208</td>
    <td>lishuo6008@outlook.com</td>
    <td>https://bioinformls.com  https://www.researchgate.net/profile/Shuo_Li37 </td>
  </tr>

  <tr>
    <td>靳展</td>
    <td>Tono4kiha</td>
    <td>1165777233@qq.com</td>
    <td>https://github.com/Tono4kiha</td>
  </tr>

</table>

