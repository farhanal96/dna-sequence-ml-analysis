"""Feature extraction utilities for DNA sequences."""

import re
import numpy as np
import pandas as pd
from collections import Counter


NUCLEOTIDES = ["A", "T", "G", "C"]
DINUCLEOTIDES = [a + b for a in NUCLEOTIDES for b in NUCLEOTIDES]
TRINUCLEOTIDES = [a + b + c for a in NUCLEOTIDES for b in NUCLEOTIDES for c in NUCLEOTIDES]


def gc_content(seq: str) -> float:
    seq = seq.upper()
    if len(seq) == 0:
        return 0.0
    g = seq.count("G")
    c = seq.count("C")
    return (g + c) / len(seq)


def nucleotide_frequencies(seq: str) -> dict:
    seq = seq.upper()
    n = len(seq) or 1
    return {nt: seq.count(nt) / n for nt in NUCLEOTIDES}


def kmer_frequencies(seq: str, k: int = 2) -> dict:
    seq = seq.upper()
    n = len(seq) - k + 1
    if n <= 0:
        return {km: 0.0 for km in (DINUCLEOTIDES if k == 2 else TRINUCLEOTIDES)}
    kmers = DINUCLEOTIDES if k == 2 else TRINUCLEOTIDES
    counts = Counter(seq[i : i + k] for i in range(n))
    return {km: counts.get(km, 0) / n for km in kmers}


def at_gc_ratio(seq: str) -> float:
    seq = seq.upper()
    gc = seq.count("G") + seq.count("C")
    at = seq.count("A") + seq.count("T")
    return at / gc if gc > 0 else 0.0


def cpg_ratio(seq: str) -> float:
    """Observed / expected CpG ratio — relevant for promoter detection."""
    seq = seq.upper()
    n = len(seq)
    if n < 2:
        return 0.0
    obs = seq.count("CG")
    c = seq.count("C")
    g = seq.count("G")
    expected = (c * g) / n if n > 0 else 1
    return obs / expected if expected > 0 else 0.0


def runs(seq: str) -> dict:
    """Longest homopolymer run per nucleotide (normalised by seq length)."""
    seq = seq.upper()
    n = len(seq) or 1
    result = {}
    for nt in NUCLEOTIDES:
        pattern = nt + "+"
        matches = re.findall(pattern, seq)
        result[f"max_run_{nt}"] = max((len(m) for m in matches), default=0) / n
    return result


def extract_features(seq: str) -> dict:
    """Return a flat feature dict for one sequence."""
    seq = seq.upper().replace(" ", "").replace("\n", "")
    feats = {}
    feats["length"] = len(seq)
    feats["gc_content"] = gc_content(seq)
    feats["at_gc_ratio"] = at_gc_ratio(seq)
    feats["cpg_ratio"] = cpg_ratio(seq)
    feats.update(nucleotide_frequencies(seq))
    feats.update(kmer_frequencies(seq, k=2))
    feats.update(kmer_frequencies(seq, k=3))
    feats.update(runs(seq))
    return feats


def sequences_to_dataframe(records: list) -> pd.DataFrame:
    """
    records: list of dicts with keys 'id', 'description', 'sequence', 'label'
    Returns a feature DataFrame with label column.
    """
    rows = []
    for rec in records:
        feats = extract_features(rec["sequence"])
        feats["id"] = rec["id"]
        feats["label"] = rec["label"]
        rows.append(feats)
    df = pd.DataFrame(rows)
    df = df.set_index("id")
    return df
