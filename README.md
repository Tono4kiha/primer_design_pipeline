# sdY primer-pair exhaustive screening pipeline

Purpose: enumerate 18-27 nt primer candidates across one Chinook salmon sdY genomic region, filter single oligos, construct 50-700 bp primer pairs, rank them, and then evaluate genome-wide specificity using BLAST hits for both primers.

## 1. Input

`input/sdy_region.fa`: exactly one FASTA record containing upstream flank + sdY gene body + downstream flank.

Optional `regions.tsv`: 0-based, half-open coordinates.

```text
label    start    end
upstream 0        500
exon1    500      650
intron1  650      1200
...
```

## 2. Install

```bash
pip install -r requirements.txt
```

Install NCBI BLAST+ separately if doing local genome BLAST.

## 3. Generate single-primer candidates

```bash
python 01_generate_candidates.py \
  --fasta input/sdy_region.fa \
  --regions input/regions.tsv \
  --outdir out
```

Outputs:
- `out/candidates_all.tsv`
- `out/candidates_pass.tsv`

## 4. Build primer pairs

```bash
python 02_pair_primers.py \
  --fasta input/sdy_region.fa \
  --candidates out/candidates_pass.tsv \
  --out out/pairs_pass.tsv
```

The output contains all pairs passing product-size, delta-Tm, heterodimer, and 3-prime interaction filters, plus a pre-BLAST score.

## 5. Select top pairs and prepare BLAST queries

```bash
python 03_prepare_blast.py \
  --fasta input/sdy_region.fa \
  --pairs out/pairs_pass.tsv \
  --top 200 \
  --outdir out/blast
```

Outputs:
- `out/blast/top_pairs.tsv`
- `out/blast/blast_primers.fa`
- `out/blast/blast_amplicons.fa`

## 6. Build a local Chinook genome BLAST database

```bash
makeblastdb \
  -in chinook_genome.fa \
  -dbtype nucl \
  -parse_seqids \
  -out chinook_genome_db
```

## 7. BLAST all F/R primers

```bash
blastn \
  -task blastn-short \
  -query out/blast/blast_primers.fa \
  -db chinook_genome_db \
  -dust no \
  -soft_masking false \
  -evalue 1000 \
  -max_target_seqs 10000 \
  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qseq sseq' \
  -out out/blast/primer_hits.tsv
```

Optional whole-amplicon BLAST can be run separately, but primer-pair reconstruction is the primary specificity test.

## 8. Reconstruct possible genome-wide PCR products

```bash
python 04_parse_blast.py \
  --pairs out/blast/top_pairs.tsv \
  --blast out/blast/primer_hits.tsv \
  --out out/blast/pairs_with_blast.tsv
```

`blast_status` values:
- `UNIQUE_PRODUCT`: exactly one plausible genome product under the configured BLAST/PCR criteria.
- `MULTIPLE_PRODUCTS`: more than one plausible product; likely off-target risk.
- `NO_PRODUCT_IN_REFERENCE`: no reconstructed product; interpret carefully if the reference assembly does not contain the intended sdY locus.

## 9. Parameters to change first

Edit `settings.py` rather than changing scripts.

Important groups:
- primer length: 18-27 nt
- primer Tm: 57-63 C, optimum 60 C
- GC: 40-60%
- product size: 50-700 bp
- pair delta-Tm: <=2.5 C
- secondary-structure cutoffs: 47 C
- PCR salt/Mg/dNTP/oligo concentrations
- BLAST 3-prime mismatch tolerance
- off-target product reconstruction window

For final experimental design, PCR chemistry parameters should be changed to match the intended reaction mixture.
