from __future__ import annotations

from collections import Counter
from typing import Dict, List, Tuple, Optional, Any
import re


def wrap_sw_align(
    s: str,
    m: str,
    match: int = 2,
    mismatch: int = -7,
    gap: int = -7,
) -> Dict[str, Any]:
    """
    Smith-Waterman local alignment of sequence s against a circular motif m.
    Motif index wraps (j-1) % L, enabling cyclic phase.

    Returns:
      {
        "score": int,
        "region": (start_idx, end_idx) on s for aligned chars (ungapped s coords),
        "aligned_s": str,
        "aligned_m": str,
        "motif_len": int,
        "motif": str
      }
    """
    n, L = len(s), len(m)
    if n == 0 or L == 0:
        return {
            "score": 0,
            "region": (None, None),
            "aligned_s": "",
            "aligned_m": "",
            "motif_len": L,
            "motif": m,
        }

    # H: DP scores; P: pointers (0=None, 1=diag, 2=up, 3=left)
    H = [[0] * L for _ in range(n + 1)]
    P = [[0] * L for _ in range(n + 1)]
    best = (0, 0, 0)  # (score, i, j)

    for i in range(1, n + 1):
        si = s[i - 1]
        for j in range(L):
            pj = (j - 1) % L
            sub = match if si == m[j] else mismatch

            diag = H[i - 1][pj] + sub
            up = H[i - 1][j] + gap
            left = H[i][pj] + gap

            val, ptr = 0, 0
            if diag >= val:
                val, ptr = diag, 1
            if up > val:
                val, ptr = up, 2
            if left > val:
                val, ptr = left, 3
            if val < 0:
                val, ptr = 0, 0

            H[i][j] = val
            P[i][j] = ptr

            if val > best[0]:
                best = (val, i, j)

    score, i, j = best

    # traceback
    a_s: List[str] = []
    a_m: List[str] = []
    used_i: List[int] = []

    while i > 0 and P[i][j] != 0 and H[i][j] > 0:
        ptr = P[i][j]
        if ptr == 1:  # diag
            a_s.append(s[i - 1])
            a_m.append(m[j])
            used_i.append(i - 1)
            i, j = i - 1, (j - 1) % L
        elif ptr == 2:  # up (consume s, gap in motif)
            a_s.append(s[i - 1])
            a_m.append("-")
            used_i.append(i - 1)
            i = i - 1
        else:  # left (gap in s, consume motif)
            a_s.append("-")
            a_m.append(m[j])
            j = (j - 1) % L

    a_s.reverse()
    a_m.reverse()

    region = (min(used_i), max(used_i) + 1) if used_i else (None, None)
    return {
        "score": score,
        "region": region,
        "aligned_s": "".join(a_s),
        "aligned_m": "".join(a_m),
        "motif_len": L,
        "motif": m,
    }


def split_units(aligned_s: str, aligned_m: str, L: int) -> List[Tuple[str, str, int, int]]:
    """
    Cut the gapped alignment into motif-sized units by counting motif advances.
    Returns a list of (s_seg, m_seg, start_col, end_col) in alignment columns.
    """
    units: List[Tuple[str, str, int, int]] = []
    start = 0
    adv = 0
    for k, ch in enumerate(aligned_m):
        if ch != "-":        # motif advanced here (diag or left)
            adv += 1
            if adv % L == 0: # completed one motif cycle
                units.append((aligned_s[start : k + 1], aligned_m[start : k + 1], start, k + 1))
                start = k + 1
    return units


def rle(labels: List[str]) -> List[Tuple[str, int]]:
    """Run-length encode labels -> list of (label, count)."""
    if not labels:
        return []
    out: List[Tuple[str, int]] = []
    cur = labels[0]
    cnt = 1
    for x in labels[1:]:
        if x == cur:
            cnt += 1
        else:
            out.append((cur, cnt))
            cur, cnt = x, 1
    out.append((cur, cnt))
    return out


def label_units_with_catalog(
    units: List[Tuple[str, str, int, int]],
    L: int,
    catalog: Dict[str, str],
) -> Tuple[List[str], Dict[str, str]]:
    """
    Assign each unique (ungapped) unit sequence a global ID like '52mer1', '52mer2', ...
    catalog: dict seq->id, mutated in place. Returns labels (per unit) and catalog.
    """
    prefix = f"{L}mer"

    # find next index for this prefix
    taken_nums: List[int] = []
    for v in catalog.values():
        if v.startswith(prefix):
            suffix = v[len(prefix):]
            if suffix.isdigit():
                taken_nums.append(int(suffix))
    next_num = (max(taken_nums) + 1) if taken_nums else 1

    labels: List[str] = []
    for s_seg, m_seg, *_ in units:
        seq = s_seg.replace("-", "")
        if seq not in catalog:
            catalog[seq] = f"{prefix}{next_num}"
            next_num += 1
        labels.append(catalog[seq])

    return labels, catalog


def encode_vntr_in_place(
    s: str,
    m: str,
    catalog: Optional[Dict[str, str]] = None,
    match: int = 2,
    mismatch: int = -7,
    gap: int = -7,
) -> Tuple[str, Dict[str, str], List[Tuple[str, int]]]:
    """
    Align a VNTR region in s to cyclic motif m and replace the aligned region with tags
    like <52mer1>2<52mer3>1 while keeping DNA flanks.

    Returns:
      represented: str (s with VNTR replaced by tags)
      catalog: dict seq->id (global unit IDs across calls)
      compressed: run-length encoded list of (unit_label, count)
    """
    if catalog is None:
        catalog = {}

    res = wrap_sw_align(s, m, match=match, mismatch=mismatch, gap=gap)
    start, end = res["region"]

    if start is None:
        return s, catalog, []

    units = split_units(res["aligned_s"], res["aligned_m"], res["motif_len"])
    labels, catalog = label_units_with_catalog(units, res["motif_len"], catalog)
    compressed = rle(labels)
    tag = "".join(f"<{lab}>{cnt}" for lab, cnt in compressed)

    represented = s[:start] + tag + s[end:]
    return represented, catalog, compressed


def unit_midline(s_seg: str, m_seg: str) -> str:
    """'|' for match, '.' for mismatch, ' ' for gaps."""
    out: List[str] = []
    for a, b in zip(s_seg, m_seg):
        if a == "-" or b == "-":
            out.append(" ")
        elif a == b:
            out.append("|")
        else:
            out.append(".")
    return "".join(out)


def unit_stats(s_seg: str, m_seg: str) -> Dict[str, int]:
    """
    Per-unit alignment stats.
      matches: aligned matches
      mismatches: aligned mismatches
      ins: insertion vs motif (s has base, motif has '-')
      dels: deletion vs motif (s has '-', motif has base)
      s_len: ungapped length of s_seg
    """
    matches = sum(1 for a, b in zip(s_seg, m_seg) if a != "-" and b != "-" and a == b)
    mism = sum(1 for a, b in zip(s_seg, m_seg) if a != "-" and b != "-" and a != b)
    ins = sum(1 for a, b in zip(s_seg, m_seg) if a != "-" and b == "-")
    dels = sum(1 for a, b in zip(s_seg, m_seg) if a == "-" and b != "-")
    s_len = sum(1 for a in s_seg if a != "-")
    return {"matches": matches, "mismatches": mism, "ins": ins, "dels": dels, "s_len": s_len}

