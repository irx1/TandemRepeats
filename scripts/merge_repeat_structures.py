# Utilities to merge/label alleles (encodedSeq strings) by motif-structure patterns.

import re
import gc
import pandas as pd


def merge_adjacent_motifs(motif_list):
    """
    Merge consecutive identical motifs.

    Input:  list of (motif_token, count) where motif_token already includes angle brackets
            e.g. [("<CTG>", 10), ("<CTG>", 2), ("<CTA>", 1)].
    Output: merged list in same format.
    """
    merged = []
    for motif, count in motif_list:
        if merged and merged[-1][0] == motif:
            merged[-1] = (motif, merged[-1][1] + count)
        else:
            merged.append((motif, count))
    return merged


def _encodedseq_to_pattern(encoded_seq, *, remove_single_letter=True, remove_ints=False, int_length=1, min_motif_len=2):
    """
    Convert an encodedSeq like '<CTG>10<CTA>1<CTG>15' into a motif-pattern key that
    preserves order but drops small artifacts and merges split motifs.

    Returns a string like '<CTG><CTA><CTG>' (counts removed, order preserved).
    """
    matches = re.findall(r"<(.*?)>(\d+)", encoded_seq)
    if not matches:
        return None

    # filter motifs/counts
    filtered = []
    for motif, count_str in matches:
        c = int(count_str)
        if remove_single_letter and len(motif) == 1:
            continue
        if len(motif) < min_motif_len:
            continue
        if remove_ints and c <= int_length:
            continue
        filtered.append((motif, c))

    if not filtered:
        return None

    merged = merge_adjacent_motifs([(f"<{m}>", c) for m, c in filtered])
    return "".join([m for m, _ in merged])


def summarize_structure_haplotypes(encoded_seqs, *, threshold=0.01, remove_single_letter=True,
                                  remove_ints=False, int_length=1, min_motif_len=2):
    """
    Identify major motif-structure patterns among encodedSeqs and summarize count ranges per motif step.

    Parameters
    ----------
    encoded_seqs : list[str] or pd.Series
    threshold : float
        Minimum fraction of sequences required for a pattern to be considered "major".
    remove_single_letter : bool
        Drop motifs of length 1.
    remove_ints : bool
        Drop motif blocks with count <= int_length.
    int_length : int
        Threshold for remove_ints.
    min_motif_len : int
        Minimum motif length to keep (default 2 keeps di/tri/etc).

    Returns
    -------
    pd.DataFrame with columns:
      - motif_pattern
      - hap_freq
      - num_sequences
      - haplotype_representation   (ranges per motif step, like "<CTG>10-30<CTA>1-3")
      - major_allele
      - major_allele_freq
    """
    if isinstance(encoded_seqs, pd.Series):
        encoded_seqs = encoded_seqs.dropna().tolist()
    else:
        encoded_seqs = [x for x in encoded_seqs if isinstance(x, str) and x]

    parsed = []
    for seq in encoded_seqs:
        matches = re.findall(r"<(.*?)>(\d+)", seq)
        if not matches:
            continue

        filtered = []
        for motif, count_str in matches:
            c = int(count_str)
            if remove_single_letter and len(motif) == 1:
                continue
            if len(motif) < min_motif_len:
                continue
            if remove_ints and c <= int_length:
                continue
            filtered.append((motif, c))

        if not filtered:
            continue

        merged = merge_adjacent_motifs([(f"<{m}>", c) for m, c in filtered])
        motif_pattern = "".join([m for m, _ in merged])
        motif_counts = {idx: c for idx, (_, c) in enumerate(merged)}

        parsed.append({"encodedSeq": seq, "motif_pattern": motif_pattern, "motif_counts": motif_counts})

    if not parsed:
        return pd.DataFrame(columns=[
            "motif_pattern", "hap_freq", "num_sequences", "haplotype_representation",
            "major_allele", "major_allele_freq"
        ])

    dfp = pd.DataFrame(parsed)
    pattern_freq = dfp["motif_pattern"].value_counts(normalize=True)
    major_patterns = pattern_freq[pattern_freq >= threshold].index.tolist()

    haplotypes = []
    for pattern in major_patterns:
        subset = dfp[dfp["motif_pattern"] == pattern]
        counts = [row["motif_counts"] for _, row in subset.iterrows()]
        counts_df = pd.DataFrame(counts).fillna(0).astype(int)

        # motif order in this pattern (strip < >)
        motif_order = re.findall(r"<(.*?)>", pattern)

        parts = []
        for idx, motif in enumerate(motif_order):
            if idx not in counts_df.columns:
                continue
            mn = int(counts_df[idx].min())
            mx = int(counts_df[idx].max())
            if mn == mx:
                parts.append(f"<{motif}>{mn}")
            else:
                parts.append(f"<{motif}>{mn}-{mx}")

        haplotype_representation = "".join(parts)

        major_allele = subset["encodedSeq"].value_counts().idxmax()
        major_allele_count = int(subset["encodedSeq"].value_counts().max())
        major_allele_freq = round(major_allele_count / len(subset), 4)

        haplotypes.append({
            "motif_pattern": pattern,
            "hap_freq": len(subset) / len(dfp),
            "num_sequences": len(subset),
            "haplotype_representation": haplotype_representation,
            "major_allele": major_allele,
            "major_allele_freq": major_allele_freq
        })

    return pd.DataFrame(haplotypes).sort_values(["hap_freq", "num_sequences"], ascending=False).reset_index(drop=True)


def merge_alleles_by_structure(df, *, encoded_seq_column="encodedSeq", group_column="TRID",
                              threshold=0.02, remove_single_letter=True, remove_ints=False,
                              int_length=1, min_motif_len=2, output_column="allele_structure"):
    """
    Assign each allele (row) to a major motif-structure haplotype within each group (e.g., TRID).
    Alleles whose structure doesn't match a major pattern are labeled "Other".

    This is what your prior code did (assign common haplotypes), but with a clearer name
    and returning ONLY the new column (pd.Series) so you can attach it to df.

    Returns
    -------
    pd.Series aligned to df.index with values like:
      - "<CTG>10-30<CTA>1-3" (haplotype_representation)
      - "Other"
      - "Unassigned" (if encodedSeq missing)
    """
    out = pd.Series(index=df.index, data="Unassigned", name=output_column)

    working = df.dropna(subset=[encoded_seq_column]).copy()

    for trid, group in working.groupby(group_column):
        encoded_seqs = group[encoded_seq_column].tolist()

        hap_df = summarize_structure_haplotypes(
            encoded_seqs,
            threshold=threshold,
            remove_single_letter=remove_single_letter,
            remove_ints=remove_ints,
            int_length=int_length,
            min_motif_len=min_motif_len
        )

        if hap_df.empty:
            out.loc[group.index] = "Other"
            continue

        # map motif-pattern -> human-readable representation
        pattern_to_rep = dict(zip(hap_df["motif_pattern"], hap_df["haplotype_representation"]))

        def seq_to_rep(seq):
            pat = _encodedseq_to_pattern(
                seq,
                remove_single_letter=remove_single_letter,
                remove_ints=remove_ints,
                int_length=int_length,
                min_motif_len=min_motif_len
            )
            if pat is None:
                return "Other"
            return pattern_to_rep.get(pat, "Other")

        out.loc[group.index] = group[encoded_seq_column].apply(seq_to_rep)

        del encoded_seqs, hap_df, pattern_to_rep
        gc.collect()

    return out
