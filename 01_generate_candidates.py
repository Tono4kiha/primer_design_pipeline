#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path

import settings as S
from primer_utils import (
    cheap_filter,
    end_gc_5,
    gc_pct,
    interval_labels,
    load_regions,
    max_homopolymer,
    max_tandem_copies,
    read_single_fasta,
    revcomp,
    shannon_entropy,
    single_thermo,
)

FIELDS = [
    'candidate_id', 'strand', 'start', 'end', 'length', 'primer_seq',
    'binding_seq_plus', 'region', 'gc', 'tm', 'end_gc5',
    'homopolymer_max', 'entropy', 'tandem_copies_max',
    'hairpin_tm', 'homodimer_tm', 'self_end_tm', 'status', 'fail_reason'
]


def evaluate(candidate_id, strand, start, end, primer_seq, binding_seq, regions):
    row = {
        'candidate_id': candidate_id,
        'strand': strand,
        'start': start,
        'end': end,
        'length': len(primer_seq),
        'primer_seq': primer_seq,
        'binding_seq_plus': binding_seq,
        'region': interval_labels(start, end, regions),
        'gc': round(gc_pct(primer_seq), 3) if all(b in 'ACGT' for b in primer_seq) else 'NA',
        'tm': 'NA',
        'end_gc5': end_gc_5(primer_seq),
        'homopolymer_max': max_homopolymer(primer_seq),
        'entropy': round(shannon_entropy(primer_seq), 4),
        'tandem_copies_max': max_tandem_copies(primer_seq),
        'hairpin_tm': 'NA',
        'homodimer_tm': 'NA',
        'self_end_tm': 'NA',
        'status': 'FAIL',
        'fail_reason': '',
    }

    reasons = cheap_filter(primer_seq)
    if reasons:
        row['fail_reason'] = ';'.join(reasons)
        return row

    tm, hairpin_tm, homodimer_tm, self_end_tm = single_thermo(primer_seq)
    row['tm'] = round(tm, 3)
    row['hairpin_tm'] = round(hairpin_tm, 3)
    row['homodimer_tm'] = round(homodimer_tm, 3)
    row['self_end_tm'] = round(self_end_tm, 3)

    if not (S.PRIMER_MIN_TM <= tm <= S.PRIMER_MAX_TM):
        reasons.append('Tm')
    if hairpin_tm > S.MAX_HAIRPIN_TM:
        reasons.append('hairpin')
    if homodimer_tm > S.MAX_HOMODIMER_TM:
        reasons.append('homodimer')
    if self_end_tm > S.MAX_SELF_END_TM:
        reasons.append('self_3prime')

    if reasons:
        row['fail_reason'] = ';'.join(reasons)
    else:
        row['status'] = 'PASS'
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--fasta', required=True, help='One-record FASTA: sdY gene body + flanks')
    ap.add_argument('--regions', default=None, help='Optional TSV: label,start,end; 0-based half-open')
    ap.add_argument('--outdir', default='out')
    args = ap.parse_args()

    Path(args.outdir).mkdir(parents=True, exist_ok=True)
    _, seq = read_single_fasta(args.fasta)
    regions = load_regions(args.regions)

    all_path = Path(args.outdir) / 'candidates_all.tsv'
    pass_path = Path(args.outdir) / 'candidates_pass.tsv'
    stats = Counter()

    with open(all_path, 'w', newline='') as all_fh, open(pass_path, 'w', newline='') as pass_fh:
        all_writer = csv.DictWriter(all_fh, fieldnames=FIELDS, delimiter='\t')
        pass_writer = csv.DictWriter(pass_fh, fieldnames=FIELDS, delimiter='\t')
        all_writer.writeheader()
        pass_writer.writeheader()

        n = len(seq)
        for k in range(S.PRIMER_MIN_LEN, S.PRIMER_MAX_LEN + 1):
            for start in range(0, n - k + 1):
                end = start + k
                binding = seq[start:end]

                f_seq = binding
                f_id = f'F_{start}_{end}_{k}'
                f_row = evaluate(f_id, 'F', start, end, f_seq, binding, regions)
                all_writer.writerow(f_row)
                stats[f_row['status']] += 1
                if f_row['status'] == 'PASS':
                    pass_writer.writerow(f_row)

                r_seq = revcomp(binding)
                r_id = f'R_{start}_{end}_{k}'
                r_row = evaluate(r_id, 'R', start, end, r_seq, binding, regions)
                all_writer.writerow(r_row)
                stats[r_row['status']] += 1
                if r_row['status'] == 'PASS':
                    pass_writer.writerow(r_row)

    print(f'Sequence length: {len(seq)} bp')
    print(f'PASS candidates: {stats["PASS"]}')
    print(f'FAIL candidates: {stats["FAIL"]}')
    print(f'Wrote: {all_path}')
    print(f'Wrote: {pass_path}')


if __name__ == '__main__':
    main()
