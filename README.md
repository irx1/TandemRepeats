# TandemRepeats

Utilities for decomposing and analyzing tandem repeat sequences from longread data 

# Core Functionality

segment_repeat : unbiased decomposition of a DNA sequence into repeat motifs and counts. Parameters to improve decomposition highlighted in example_decompose.py
	Designed for short tandem repeats (motif sizes <10)

encode_vntr_in_place : decomposition of VNTRs using tandem repeat finder algorithm for alignment. Requires a consensus motif.
	
	See scripts/vntr_decompose_example.py for example of VNTR decomposition 

