from repeat_utils import count_kmers, simplest_motif, max_motif, extract_motifs
import pandas as pd
from collections import Counter

def segment_repeat(
    s,
    min_k=2,
    max_k=7,
    reverse=False,
    trimEnds=True,
    setkmers=None,
    include_kmers=None,
    max_unmatched_threshold=None,
):

    encodedSeq = None
    kmers = []
    nMosaicMotifs = None

    try:
        if not isinstance(s, str):
            return pd.Series(
                {"kmers": None, "encodedSeq": None, "nMosaicMotifs": None}
            )

        s2 = s  # [20:-20]  # Trim flanks if desired

        if setkmers is not None:
            kmers = sorted(setkmers, key=len, reverse=True)
        else:
            kmers = count_kmers(s2, min_k=min_k, max_k=max_k)

            has_kmers_4plus = any(len(kmer) >= 4 for kmer in kmers)

            kmers = {
                kmer: count
                for kmer, count in kmers.items()
                if s.find(
                    kmer
                    * (
                        4 if has_kmers_4plus and len(kmer) == 2
                        else 3 if has_kmers_4plus and len(kmer) == 3
                        else 8 if has_kmers_4plus and len(kmer) == 1
                        else 2
                    )
                )
                != -1
            }

            kmers = sorted(set(kmers.keys()), key=len, reverse=True)
            kmers = list(set(simplest_motif(kmer) for kmer in kmers))

            if reverse:
                kmers = list(set(max_motif(kmer) for kmer in kmers))

            if include_kmers is not None:
                kmers = list(set(kmers).union(include_kmers))

            kmers = sorted(kmers, key=len, reverse=True)

#        if not kmers:
#            return pd.Series({"kmers": [], "encodedSeq": None, "nMosaicMotifs": 0})

        all_encoded_sequences = []
        current_encoded_seq = []
        unmatched_seq = ""
        i = 0
        no_match_count = 0

        longestKmerLen = max(len(kmer) for kmer in kmers)

        if max_unmatched_threshold is None:
            max_unmatched_threshold = max(15, max_k + 1)

        while i < len(s):
            found_kmer = False

            for kmer in kmers:
                k = len(kmer)
                if s[i : i + k] == kmer:
                    if unmatched_seq:
                        current_encoded_seq.append(
                            (unmatched_seq[0], len(unmatched_seq))
                            if len(set(unmatched_seq)) == 1
                            else (unmatched_seq, 1)
                        )
                        unmatched_seq = ""

                    count = 0
                    while i + k <= len(s) and s[i : i + k] == kmer:
                        count += 1
                        i += k
                        if any(
                            s[i : i + len(lk)] == lk
                            for lk in kmers
                            if len(lk) > k
                        ):
                            break

                    if count > 1 or (len(kmer) <= 2 and count > 5):
                        no_match_count = 0
                        if current_encoded_seq and current_encoded_seq[-1][0] == kmer:
                            current_encoded_seq[-1] = (
                                kmer,
                                current_encoded_seq[-1][1] + count,
                            )
                        else:
                            current_encoded_seq.append((kmer, count))
                    else:
                        if len(kmer) in (1, 2):
                            unmatched_seq += kmer
                            no_match_count += k
                        else:
                            current_encoded_seq.append((kmer, count))

                    found_kmer = True
                    break

            if not found_kmer:
                unmatched_seq += s[i]
                no_match_count += 1
                i += 1

            if no_match_count >= max_unmatched_threshold:
                if len(set(unmatched_seq)) > 1:
                    all_encoded_sequences.append(current_encoded_seq)
                    current_encoded_seq = []
                    unmatched_seq = ""
                    no_match_count = 0

        if unmatched_seq:
            current_encoded_seq.append(
                (unmatched_seq[0], len(unmatched_seq))
                if len(set(unmatched_seq)) == 1
                else (unmatched_seq, 1)
            )

        if current_encoded_seq:
            all_encoded_sequences.append(current_encoded_seq)

        longest_encoded_seq = max(
            all_encoded_sequences,
            key=lambda seq: max((count for _, count in seq), default=0),
            default=[],
        )

        if longest_encoded_seq:
            encodedSeq = "".join(
                f"<{kmer}>{count}" for kmer, count in longest_encoded_seq
            )

            kmers2 = extract_motifs([encodedSeq], min_count=2)
            kmers2 = [k for k in kmers2 if len(k) >= min_k]

            contains_long_kmers = any(
                len(kmer) > min_k for kmer, _ in longest_encoded_seq
            )

            if trimEnds:
                while longest_encoded_seq and (
                    longest_encoded_seq[0][1] < 2
                    or (
                        contains_long_kmers
                        and len(longest_encoded_seq[0][0]) == 1
                        and longest_encoded_seq[0][1] < 6
                    )
                    or (
                        contains_long_kmers
                        and len(longest_encoded_seq[0][0]) == 2
                        and longest_encoded_seq[0][1] < 5
                    )
                    or (
                        contains_long_kmers
                        and len(longest_encoded_seq[0][0]) == 3
                        and longest_encoded_seq[0][1] < 3
                    )
                ):
                    longest_encoded_seq.pop(0)

                while longest_encoded_seq and (
                    longest_encoded_seq[-1][1] < 2
                    or (
                        contains_long_kmers
                        and len(longest_encoded_seq[-1][0]) == 1
                        and longest_encoded_seq[-1][1] < 6
                    )
                    or (
                        contains_long_kmers
                        and len(longest_encoded_seq[-1][0]) == 2
                        and longest_encoded_seq[-1][1] < 5
                    )
                    or (
                        contains_long_kmers
                        and len(longest_encoded_seq[-1][0]) == 3
                        and longest_encoded_seq[-1][1] < 3
                    )
                ):
                    longest_encoded_seq.pop(-1)

                encodedSeq = "".join(
                    f"<{kmer}>{count}" for kmer, count in longest_encoded_seq
                )

            else:
                while longest_encoded_seq and longest_encoded_seq[0][0] not in kmers2:
                    longest_encoded_seq.pop(0)
                while longest_encoded_seq and longest_encoded_seq[-1][0] not in kmers2:
                    longest_encoded_seq.pop(-1)

                encodedSeq = "".join(
                    f"<{kmer}>{count}" for kmer, count in longest_encoded_seq
                )

        unique_kmers = {
            simplest_motif(kmer)
            for kmer, count in longest_encoded_seq
            if len(kmer) > min_k and count >= 3
        }
        nMosaicMotifs = len(unique_kmers)

    except Exception:
        encodedSeq, kmers, nMosaicMotifs = None, [], None

    return pd.Series(
        {"kmers": kmers, "encodedSeq": encodedSeq, "nMosaicMotifs": nMosaicMotifs}
    )

