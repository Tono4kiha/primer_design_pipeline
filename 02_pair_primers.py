#!/usr/bin/env python3
from __future__ import annotations

import argparse
import bisect
import csv
from collections import Counter
from pathlib import Path

import settings as S
from primer_utils import gc_pct, pair_thermo, preblast_score, read_single_fasta

PAIR_FIELDS = [
    'pair_id',
    'F_id', 'R_id', 'F_seq', 'R_seq',
    'F_start', 'F_end', 'R_start', 'R_end',
    'F_region', 'R_region',
    'F_length', 'R_length',
    'F_tm', 'R_tm', 'delta_tm',
    'F_gc', 'R_gc', 'F_end_gc5', 'R_end_gc5',
    'F_hairpin_tm', 'R_hairpin_tm',
    'F_homodimer_tm', 'R_homodimer_tm',
    'F_self_end_tm', 'R_self_end_tm',
    'heterodimer_tm', 'pair_end_tm',
    'product_start', 'product_end', 'product_size', 'product_gc',
    'score_preblast'
]


def load_candidates(path):
    left, right = [], []
    with open(path, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            if row.get('status') != 'PASS':
                continue
            row['start_i'] = int(row['start'])
            row['end_i'] = int(row['end'])
            row['length_i'] = int(row['length'])
            row['tm_f'] = float(row['tm'])
            if row['strand'] == 'F':
                left.append(row)
            elif row['strand'] == 'R':
                right.append(row)
    right.sort(key=lambda x: x['end_i'])
    return left, right


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--candidates', required=True)
    ap.add_argument('--out', default='out/pairs_pass.tsv')
    args = ap.parse_args()

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    _, template = read_single_fasta(args.fasta)
    left, right = load_candidates(args.candidates)
    right_ends = [r['end_i'] for r in right]

    stats = Counter()
    pair_no = 0

    with open(args.out, 'w', newline='') as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=PAIR_FIELDS, delimiter='\t')
        writer.writeheader()

        for f in left:
            min_r_end = f['start_i'] + S.PRODUCT_MIN
            max_r_end = f['start_i'] + S.PRODUCT_MAX
            lo = bisect.bisect_left(right_ends, min_r_end)
            hi = bisect.bisect_right(right_ends, max_r_end)

            for idx in range(lo, hi):
                r = right[idx]
                stats['coordinate_candidates'] += 1

                if r['start_i'] < f['end_i']:
                    stats['overlap'] += 1
                    continue

                delta_tm = abs(f['tm_f'] - r['tm_f'])
                if delta_tm > S.MAX_DELTA_TM:
                    stats['delta_tm'] += 1
                    continue

                product_start = f['start_i']
                product_end = r['end_i']
                product_size = product_end - product_start
                if not (S.PRODUCT_MIN <= product_size <= S.PRODUCT_MAX):
                    stats['product_size'] += 1
                    continue

                heterodimer_tm, pair_end_tm = pair_thermo(f['primer_seq'], r['primer_seq'])
                if heterodimer_tm > S.MAX_HETERODIMER_TM:
                    stats['heterodimer'] += 1
                    continue
                if pair_end_tm > S.MAX_PAIR_END_TM:
                    stats['pair_3prime'] += 1
                    continue

                product_seq = template[product_start:product_end]
                pair_no += 1
                row = {
                    'pair_id': f'PAIR{pair_no:09d}',
                    'F_id': f['candidate_id'],
                    'R_id': r['candidate_id'],
                    'F_seq': f['primer_seq'],
                    'R_seq': r['primer_seq'],
                    'F_start': f['start_i'],
                    'F_end': f['end_i'],
                    'R_start': r['start_i'],
                    'R_end': r['end_i'],
                    'F_region': f['region'],
                    'R_region': r['region'],
                    'F_length': f['length_i'],
                    'R_length': r['length_i'],
                    'F_tm': f['tm'],
                    'R_tm': r['tm'],
                    'delta_tm': round(delta_tm, 3),
                    'F_gc': f['gc'],
                    'R_gc': r['gc'],
                    'F_end_gc5': f['end_gc5'],
                    'R_end_gc5': r['end_gc5'],
                    'F_hairpin_tm': f['hairpin_tm'],
                    'R_hairpin_tm': r['hairpin_tm'],
                    'F_homodimer_tm': f['homodimer_tm'],
                    'R_homodimer_tm': r['homodimer_tm'],
                    'F_self_end_tm': f['self_end_tm'],
                    'R_self_end_tm': r['self_end_tm'],
                    'heterodimer_tm': round(heterodimer_tm, 3),
                    'pair_end_tm': round(pair_end_tm, 3),
                    'product_start': product_start,
                    'product_end': product_end,
                    'product_size': product_size,
                    'product_gc': round(gc_pct(product_seq), 3),
                    'score_preblast': 0.0,
                }
                row['score_preblast'] = round(preblast_score(row), 4)
                writer.writerow(row)
                stats['PASS'] += 1

    print(f'Forward candidates: {len(left)}')
    print(f'Reverse candidates: {len(right)}')
    for key in sorted(stats):
        print(f'{key}: {stats[key]}')
    print(f'Wrote: {args.out}')


if __name__ == '__main__':
    main()
