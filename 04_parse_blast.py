#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import settings as S

BLAST_FIELDS = [
    'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
    'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore',
    'qlen', 'qseq', 'sseq'
]


def three_prime_mismatches(hit):
    qstart = int(hit['qstart'])
    qlen = int(hit['qlen'])
    qseq = hit['qseq'].upper()
    sseq = hit['sseq'].upper()

    qpos = qstart - 1
    mismatches = 0
    covered = set()
    cutoff = qlen - S.BLAST_3P_WINDOW + 1

    for qbase, sbase in zip(qseq, sseq):
        if qbase != '-':
            qpos += 1
            if qpos >= cutoff:
                covered.add(qpos)
                if sbase == '-' or qbase != sbase:
                    mismatches += 1

    complete_3p = all(pos in covered for pos in range(cutoff, qlen + 1))
    return mismatches, complete_3p


def accept_hit(hit):
    qlen = int(hit['qlen'])
    aligned_query = sum(1 for b in hit['qseq'] if b != '-')
    coverage = aligned_query / qlen
    mm3, complete_3p = three_prime_mismatches(hit)

    if coverage < S.BLAST_MIN_QUERY_COVERAGE:
        return False, coverage, mm3
    if int(hit['mismatch']) > S.BLAST_MAX_MISMATCH:
        return False, coverage, mm3
    if int(hit['gapopen']) > S.BLAST_MAX_GAPOPEN:
        return False, coverage, mm3
    if not complete_3p:
        return False, coverage, mm3
    if mm3 > S.BLAST_MAX_3P_MISMATCH:
        return False, coverage, mm3
    return True, coverage, mm3


def normalize_hit(hit, coverage, mm3):
    sstart = int(hit['sstart'])
    send = int(hit['send'])
    return {
        'sseqid': hit['sseqid'],
        'strand': '+' if sstart <= send else '-',
        'left': min(sstart, send),
        'right': max(sstart, send),
        'coverage': coverage,
        'mismatch': int(hit['mismatch']),
        'mm3': mm3,
        'bitscore': float(hit['bitscore']),
    }


def reconstruct_products(f_hits, r_hits):
    products = {}
    r_by_subject = defaultdict(list)
    for r in r_hits:
        r_by_subject[r['sseqid']].append(r)

    for f in f_hits:
        for r in r_by_subject.get(f['sseqid'], []):
            # Expected orientation: F query matches + strand on the left,
            # R query matches - strand on the right.
            if f['strand'] == '+' and r['strand'] == '-' and f['right'] < r['left']:
                start = f['left']
                end = r['right']
                size = end - start + 1
                orientation = 'F+_R-'
            # Off-targets can also occur with primer identities swapped in genomic orientation.
            elif r['strand'] == '+' and f['strand'] == '-' and r['right'] < f['left']:
                start = r['left']
                end = f['right']
                size = end - start + 1
                orientation = 'R+_F-'
            else:
                continue

            if S.OFFTARGET_PRODUCT_MIN <= size <= S.OFFTARGET_PRODUCT_MAX:
                key = (f['sseqid'], start, end, orientation)
                products[key] = {
                    'sseqid': f['sseqid'],
                    'start': start,
                    'end': end,
                    'size': size,
                    'orientation': orientation,
                }
    return list(products.values())


def load_pairs(path):
    rows = []
    with open(path, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            rows.append(row)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pairs', required=True, help='top_pairs.tsv from 03_prepare_blast.py')
    ap.add_argument('--blast', required=True, help='primer_hits.tsv from blastn-short')
    ap.add_argument('--out', default='out/blast/pairs_with_blast.tsv')
    args = ap.parse_args()

    pair_hits = defaultdict(lambda: {'F': [], 'R': []})

    with open(args.blast, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t', fieldnames=BLAST_FIELDS)
        for hit in reader:
            qseqid = hit['qseqid']
            if '_' not in qseqid:
                continue
            pair_id, role = qseqid.rsplit('_', 1)
            if role not in ('F', 'R'):
                continue

            ok, coverage, mm3 = accept_hit(hit)
            if not ok:
                continue
            pair_hits[pair_id][role].append(normalize_hit(hit, coverage, mm3))

    rows = load_pairs(args.pairs)
    extra_fields = [
        'blast_F_accepted_hits', 'blast_R_accepted_hits',
        'predicted_genome_amplicons', 'blast_status', 'genome_amplicons_detail'
    ]
    out_fields = list(rows[0].keys()) + extra_fields if rows else extra_fields

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, 'w', newline='') as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=out_fields, delimiter='\t')
        writer.writeheader()

        for row in rows:
            pid = row['pair_id']
            f_hits = pair_hits[pid]['F']
            r_hits = pair_hits[pid]['R']
            products = reconstruct_products(f_hits, r_hits)
            nprod = len(products)

            if nprod == 1:
                status = 'UNIQUE_PRODUCT'
            elif nprod == 0:
                status = 'NO_PRODUCT_IN_REFERENCE'
            else:
                status = 'MULTIPLE_PRODUCTS'

            detail = ';'.join(
                f"{p['sseqid']}:{p['start']}-{p['end']}:{p['size']}:{p['orientation']}"
                for p in sorted(products, key=lambda x: (x['sseqid'], x['start'], x['end']))
            )

            row.update({
                'blast_F_accepted_hits': len(f_hits),
                'blast_R_accepted_hits': len(r_hits),
                'predicted_genome_amplicons': nprod,
                'blast_status': status,
                'genome_amplicons_detail': detail,
            })
            writer.writerow(row)

    print(f'Wrote: {args.out}')
    print('Interpretation: UNIQUE_PRODUCT is strongest; MULTIPLE_PRODUCTS suggests off-target risk;')
    print('NO_PRODUCT_IN_REFERENCE is inconclusive if the reference assembly lacks the target sdY locus.')


if __name__ == '__main__':
    main()
