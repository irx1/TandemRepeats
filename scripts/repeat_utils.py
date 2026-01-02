import pandas as pd
import re
from collections import Counter
from itertools import groupby
import numpy as np
import gc

def calculate_lps(encoded_seq):
    """
    Longest Pure Segment (LPS): return the single motif block with the largest
    repeat count from an encodedSeq.

    Example:
      "<CGG>9<AGG>1<CGG>9<AGG>1<CGG>9" -> "<CGG>9"

    Notes:
    - "Pure segment" here means one uninterrupted <motif>count block.
    - Ties are broken by choosing the first max block encountered.
    - Returns None if input is None/invalid or contains no <motif>count pairs.
    """
    if not isinstance(encoded_seq, str) or not encoded_seq:
        return None

    pairs = re.findall(r"<([^>]+)>(\d+)", encoded_seq)
    if not pairs:
        return None

    best_motif, best_count = None, -1
    for motif, count_str in pairs:
        c = int(count_str)
        if c > best_count:
            best_motif, best_count = motif, c

    return f"<{best_motif}>{best_count}"



def remove_singletons(encoded_seq, min_count=2, min_motif_len=2, min_bp_prop=0.01, main_kmers=None):
    """
    Simplify encodedSeqs by removing small motifs and merging adjacent repeats.

    Parameters:
    - encoded_seq: input string like '<AAG>10<GAG>1<AAG>15'
    - min_count: minimum count to keep motif
    - min_motif_len: minimum length of motif to keep
    - min_bp_prop: minimum proportion of base pairs (motif_len * count / total bp)
    - main_kmers: optional set of motifs to always keep (e.g., {"AGGAT", "ATGGG"})

    Returns:
    - Cleaned encoded string
    """

    try:
        pattern = r"<([^>]+)>(\d+)"
        pairs = [(motif, int(count)) for motif, count in re.findall(pattern, encoded_seq)]
    
        total_bp = sum(len(m) * c for m, c in pairs)
    
        # Keep if:
        # - Count ≥ min_count
        # - Motif length ≥ min_motif_len
        # - Contributes ≥ min_bp_prop of total bp
        # - OR it's in main_kmers (if provided)
        filtered = [
            (m, c)
            for m, c in pairs
            if (c >= min_count and len(m) >= min_motif_len and (len(m) * c / total_bp) >= min_bp_prop)
            or (main_kmers is not None and m in main_kmers)
        ]
    
        # Merge adjacent identical motifs
        merged = []
        for motif, group in groupby(filtered, key=lambda x: x[0]):
            total_count = sum(c for _, c in group)
            merged.append(f"<{motif}>{total_count}")
    
        return ''.join(merged)
    except:
        return None

def calculate_repeat_length(encoded_seq):

    if not isinstance(encoded_seq, str):
        return 0  
    pattern = r'<([A-Za-z]+)>(\d+)'
    matches = re.findall(pattern, encoded_seq)
    total_length = sum(len(motif) * int(count) for motif, count in matches)
    
    return total_length

def reverse_complement(motif):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(motif))

def simplest_motif(motif):
    
    i = (motif + motif).find(motif, 1, -1)
    if i != -1:
        motif = motif[:i]
    return min(motif[i:] + motif[:i] for i in range(len(motif)))

def reverse_motif(motif):
    
    complement = str.maketrans("ATCG", "TAGC")
    return motif.translate(complement)[::-1]

def reverse_kmers(kmers):
    return [kmer[::-1] for kmer in kmers]

def max_motif(motif):
    
    i = (motif + motif).find(motif, 1, -1)
    if i != -1:
        motif = motif[:i]
    return max(motif[i:] + motif[:i] for i in range(len(motif)))


def extract_motifs(encoded_seqs, min_count=2):
    motif_counter = Counter()
    for seq in encoded_seqs:
        if not isinstance(seq, str):
            return []
        matches = re.findall(r'<([^>]+)>(\d+)', seq)
        for motif, count in matches:
            if int(count) >= min_count:
                motif_counter[motif] += int(count)
    return [motif for motif, count in motif_counter.most_common() if count >= min_count]

def count_kmers(s, min_k=2, max_k=7):
    kmer_dict = {}
    kmer_lengths = list(range(min_k, max_k + 1))

    for k in kmer_lengths:
        for i in range(len(s) - k + 1):
            kmer = s[i:i + k]
            kmer = simplest_motif(kmer)
            kmer_dict[kmer] = kmer_dict.get(kmer, 0) + 1

    return kmer_dict
