from __future__ import annotations

import csv
import math
from functools import lru_cache
from pathlib import Path

from primer3 import bindings

import settings as S


def read_single_fasta(path: str | Path) -> tuple[str, str]:
    name = None
    chunks = []
    with open(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    raise ValueError('Input FASTA must contain exactly one sequence record.')
                name = line[1:].split()[0]
            else:
                if name is None:
                    raise ValueError('FASTA sequence encountered before header.')
                chunks.append(line.upper())
    if name is None:
        raise ValueError('No FASTA record found.')
    seq = ''.join(chunks)
    return name, seq


def write_fasta_record(fh, name: str, seq: str, width: int = 80) -> None:
    fh.write(f'>{name}\n')
    for i in range(0, len(seq), width):
        fh.write(seq[i:i + width] + '\n')


def revcomp(seq: str) -> str:
    table = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
    return seq.translate(table)[::-1]


def gc_pct(seq: str) -> float:
    return 100.0 * (seq.count('G') + seq.count('C')) / len(seq)


def max_homopolymer(seq: str) -> int:
    best = cur = 1
    for a, b in zip(seq, seq[1:]):
        if a == b:
            cur += 1
            best = max(best, cur)
        else:
            cur = 1
    return best


def shannon_entropy(seq: str) -> float:
    n = len(seq)
    entropy = 0.0
    for base in 'ACGT':
        p = seq.count(base) / n
        if p > 0:
            entropy -= p * math.log2(p)
    return entropy


def max_tandem_copies(seq: str, unit_lengths=(2, 3)) -> int:
    best = 1
    n = len(seq)
    for u in unit_lengths:
        for i in range(0, n - 2 * u + 1):
            unit = seq[i:i + u]
            copies = 1
            j = i + u
            while j + u <= n and seq[j:j + u] == unit:
                copies += 1
                j += u
            best = max(best, copies)
    return best


def end_gc_5(seq: str) -> int:
    tail = seq[-5:]
    return tail.count('G') + tail.count('C')


def load_regions(path: str | Path | None):
    if path is None:
        return []
    regions = []
    with open(path, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        required = {'label', 'start', 'end'}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError('regions TSV requires columns: label, start, end')
        for row in reader:
            regions.append((int(row['start']), int(row['end']), row['label']))
    return regions


def interval_labels(start: int, end: int, regions) -> str:
    labels = [label for a, b, label in regions if start < b and end > a]
    return ';'.join(labels) if labels else 'NA'


def _thermo_kwargs():
    return dict(
        mv_conc=S.MV_CONC,
        dv_conc=S.DV_CONC,
        dntp_conc=S.DNTP_CONC,
        dna_conc=S.DNA_CONC,
    )


@lru_cache(maxsize=200000)
def single_thermo(seq: str) -> tuple[float, float, float, float]:
    kw = _thermo_kwargs()
    tm = bindings.calc_tm(seq, **kw)
    hairpin = bindings.calc_hairpin(
        seq, temp_c=S.THERMO_TEMP_C, max_loop=S.MAX_LOOP, **kw
    )
    homodimer = bindings.calc_homodimer(
        seq, temp_c=S.THERMO_TEMP_C, max_loop=S.MAX_LOOP, **kw
    )
    self_end = bindings.calc_end_stability(
        seq, seq, temp_c=S.THERMO_TEMP_C, max_loop=S.MAX_LOOP, **kw
    )
    return float(tm), float(hairpin.tm), float(homodimer.tm), float(self_end.tm)


@lru_cache(maxsize=500000)
def pair_thermo(seq1: str, seq2: str) -> tuple[float, float]:
    kw = _thermo_kwargs()
    hetero = bindings.calc_heterodimer(
        seq1, seq2, temp_c=S.THERMO_TEMP_C, max_loop=S.MAX_LOOP, **kw
    )
    end12 = bindings.calc_end_stability(
        seq1, seq2, temp_c=S.THERMO_TEMP_C, max_loop=S.MAX_LOOP, **kw
    )
    end21 = bindings.calc_end_stability(
        seq2, seq1, temp_c=S.THERMO_TEMP_C, max_loop=S.MAX_LOOP, **kw
    )
    return float(hetero.tm), max(float(end12.tm), float(end21.tm))


def cheap_filter(seq: str) -> list[str]:
    reasons = []
    if any(b not in 'ACGT' for b in seq):
        reasons.append('non_ACGT')
        return reasons

    gc = gc_pct(seq)
    if gc < S.PRIMER_MIN_GC or gc > S.PRIMER_MAX_GC:
        reasons.append('GC')
    if max_homopolymer(seq) > S.MAX_HOMOPOLYMER:
        reasons.append('homopolymer')
    if shannon_entropy(seq) < S.MIN_SHANNON_ENTROPY:
        reasons.append('low_entropy')
    if max_tandem_copies(seq) > S.MAX_TANDEM_COPIES:
        reasons.append('tandem_repeat')
    return reasons


def preblast_score(row: dict) -> float:
    f_tm = float(row['F_tm'])
    r_tm = float(row['R_tm'])
    f_gc = float(row['F_gc'])
    r_gc = float(row['R_gc'])
    delta_tm = abs(f_tm - r_tm)
    size = int(row['product_size'])

    score = 0.0
    score += abs(f_tm - S.PRIMER_OPT_TM)
    score += abs(r_tm - S.PRIMER_OPT_TM)
    score += 1.5 * delta_tm
    score += 0.15 * abs(f_gc - 50.0)
    score += 0.15 * abs(r_gc - 50.0)
    score += 0.004 * abs(size - S.PRODUCT_OPT)

    # Mild preference for ~2 GC bases in the last five bases; not a hard filter.
    score += 0.5 * abs(int(row['F_end_gc5']) - 2)
    score += 0.5 * abs(int(row['R_end_gc5']) - 2)

    # Penalize secondary structures as their predicted Tm approaches the cutoff.
    for key in ('F_hairpin_tm', 'R_hairpin_tm', 'F_homodimer_tm', 'R_homodimer_tm',
                'heterodimer_tm', 'pair_end_tm'):
        value = float(row[key])
        score += max(0.0, value - 35.0) / 5.0
    return score
