# TandemRepeats

Utilities for decomposing and analyzing tandem repeat sequences from long-read data.

## Methods

**segment_repeat**  
Unbiased decomposition of a DNA sequence into repeat motifs and counts.  
Designed for short tandem repeats (motif sizes <10).  
Parameters to improve decomposition are demonstrated in `scripts/example_decompose.py`.

**encode_vntr_in_place**  
Decomposition of VNTRs using motif-based alignment.  
Requires a consensus motif.  
See `scripts/vntr_decompose_example.py` for a worked VNTR example.

