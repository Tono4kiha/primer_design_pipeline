#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import heapq
from pathlib import Path

from primer_utils import read_single_fasta, write_fasta_record


def iter_pairs(path):
    with open(path, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            yield row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--pairs', required=True)
    ap.add_argument('--top', type=int, default=200)
    ap.add_argument('--outdir', default='out/blast')
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    top_rows = heapq.nsmallest(
        args.top,
        iter_pairs(args.pairs),
        key=lambda r: float(r['score_preblast'])
    )
    top_rows.sort(key=lambda r: float(r['score_preblast']))
    if not top_rows:
        raise SystemExit('No primer pairs found.')

    _, template = read_single_fasta(args.fasta)

    top_pairs_path = outdir / 'top_pairs.tsv'
    with open(top_pairs_path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=top_rows[0].keys(), delimiter='\t')
        writer.writeheader()
        writer.writerows(top_rows)

    primer_fa = outdir / 'blast_primers.fa'
    amplicon_fa = outdir / 'blast_amplicons.fa'
    map_tsv = outdir / 'blast_query_map.tsv'

    with open(primer_fa, 'w') as pfh, open(amplicon_fa, 'w') as afh, open(map_tsv, 'w', newline='') as mfh:
        map_writer = csv.writer(mfh, delimiter='\t')
        map_writer.writerow(['qseqid', 'pair_id', 'role'])

        for row in top_rows:
            pid = row['pair_id']
            fq = f'{pid}_F'
            rq = f'{pid}_R'
            write_fasta_record(pfh, fq, row['F_seq'])
            write_fasta_record(pfh, rq, row['R_seq'])
            map_writer.writerow([fq, pid, 'F'])
            map_writer.writerow([rq, pid, 'R'])

            start = int(row['product_start'])
            end = int(row['product_end'])
            write_fasta_record(afh, pid, template[start:end])

    print(f'Wrote: {top_pairs_path}')
    print(f'Wrote: {primer_fa}')
    print(f'Wrote: {amplicon_fa}')
    print(f'Wrote: {map_tsv}')
    print('\nRun primer BLAST with outfmt fields expected by 04_parse_blast.py:')
    print(
        'blastn -task blastn-short -query blast_primers.fa -db <chinook_genome_db> '
        '-dust no -soft_masking false -evalue 1000 -max_target_seqs 10000 '
        '-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend '
        'sstart send evalue bitscore qlen qseq sseq" -out primer_hits.tsv'
    )


if __name__ == '__main__':
    main()
