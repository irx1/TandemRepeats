#!/usr/bin/env python3

import sys

# allow importing from scripts/ when running from repo root
sys.path.append("scripts")

from vntr_utils import encode_vntr_in_place


def main():
    print("=" * 80)
    print("VNTR decomposition demo (WDR7)")
    print("Demonstrates: motif provided in different reading frame still aligns")
    print("=" * 80)

    # WDR7 reference sequence +- 100bp, from UCSC genome browser
    s = (
        """ATAATACTTCACATTCATGCCGAATTTGAAGAGGTCATAGGGAAAGGATG
	ATCTTGAGCAATGCATTTTAGTGAACATGGGGGTCTTTATCGTTCCACAT
	GTTATGTCCCTTACAGGATTTCACATATCCCTATTCTCAGACAGGAATAG
	GGATATGTGAGCTATGATAGTTATGTCCCTTATAGAATTTCACATATCCC
	TATTCTCAGGCAGGAATAGGGATATGTGAGCTATGATAGTTATGTCCCTT
	ATAGAATTTCACATATCCCTATTCTCAGACAGGAATAGGGATATGTGAAC
	TATGATAGTTATGTCCCTTACAGGATTTCACATATCCCTATTCTCAGACA
	GGAATAGGGATATGTGAGCTATGATAGTTATGTCCCTTATAGAATTTCAC
	ATATCCCTATTCTCAGACAGGAATAGGGATATGTGAGCTATGATAGTTAT
	GTCCCTTACAGGATTTCACATATCCCTATTCTCAGGCAGGAATAGGGATA
	TGTGAACTATGATAGTTATGTCCCTTATAGAATTTCACAAATCCCTATTC
	TCAGGCAGGAAGGCTATTTGACCTGACCAAAGGCATTTATCAAAGTAGGG
	ATGGCTCAGGTTATAAGAGGAACCCAGTTGTAATGGCCTTTTCTAAAATC
	TTTACATCTAG"""
        ).replace("\n", "")
		

    # Motif in one phase
    m1 = "GTTATGTCCCTTATAGAATTTCACATATCCCTATTCTCAGACAGGAATAGGGATATGTGAGCTATGATA"

    # Same motif content, different "reading frame" (split in half then swapped)
    m2 = (
        "TCTCAGACAGGAATAGGGATATGTGAGCTATGATA"
        "GTTATGTCCCTTATAGAATTTCACATATCCCTAT"
    )

    rep1, catalog1, comp1 = encode_vntr_in_place(s, m1, catalog={})
    rep2, catalog2, comp2 = encode_vntr_in_place(s, m2, catalog={})

    print("\n[m1] motif:")
    print(m1)
    print("\n[m1] represented sequence (VNTR replaced by <LmerX> tags):")
    print(rep1)
    print("\n[m1] compressed units (run-length encoding):")
    print(comp1)
    print("\n[m1] unique units:", len(catalog1))

    print("\n" + "-" * 80)

    print("\n[m2] motif (same motif, different phase):")
    print(m2)
    print("\n[m2] represented sequence (VNTR replaced by <LmerX> tags):")
    print(rep2)
    print("\n[m2] compressed units (run-length encoding):")
    print(comp2)
    print("\n[m2] unique units:", len(catalog2))

    print("\nSanity: tags present?")
    print("m1:", "<" in rep1, " | m2:", "<" in rep2)


if __name__ == "__main__":
    main()

