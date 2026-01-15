import os
import pandas as pd
import numpy as np
from joblib import Parallel, delayed
from tqdm import tqdm
from Levenshtein import distance as levenshtein_distance
from RepeatXuTRACTR import segment_repeat
from collections import Counter
from collections import defaultdict
from biotite.sequence.phylo import neighbor_joining as neighjoin
from biotite.sequence.phylo import upgma
import networkx as nx
from ete3 import Tree, TreeStyle, TextFace, NodeStyle, ImgFace ,  RectFace
from IPython.display import Image
from PIL import Image, ImageDraw
from matplotlib.patches import Patch
import matplotlib.image as mpimg
import tempfile  
import matplotlib.pyplot as plt

from phylo_helper_functions import (
    calculate_motif_frequencies_and_filter,
    replace_encoded_seq,
    replace_encoded_seq_simple,
    build_conv,
    create_cost_dict,
    convert_cost_dict,
    min_cost_to_convert_exact
)
from phylo_helper_functions import (
    create_tree_from_distance_matrix,
    render_tree_with_colored_sequences,
    render_tree_to_tempfile,
    calculate_max_codified_length
)



DATA_PATH = os.path.join(
    os.path.dirname(__file__),
    "..",
    "data",
    "edit_distance_test_data.csv"
)

df = pd.read_csv(DATA_PATH)
df['encodedSeq'] = df['Allele'].apply(segment_repeat, setkmers=['AAAGG','AAGGG', 'AAGG', 'AAAG', 'AGGG', 'AAG','A'])['encodedSeq']

motif_frequencies, total_alleles, filtered_motifs = \
    calculate_motif_frequencies_and_filter(df["encodedSeq"], min_count=2)

conv = build_conv(
    filtered_motifs,
)

df["codified_seq"] = df["encodedSeq"].apply(lambda x: replace_encoded_seq(x, conv))
df["codified_seq_simple"] = df['encodedSeq'].apply(lambda x: replace_encoded_seq_simple(x, conv))
df["codified_seq_noO"] = df["codified_seq"].str.replace("O", "", regex=False)

test_locus = df.copy().reset_index(drop=True)

# Option 1: keep only unique encoded sequences
# test_locus = test_locus[~test_locus["codified_seq"].duplicated()].reset_index(drop=True)

# Option 2: unique after removing "Other" symbols
# test_locus = test_locus[~test_locus["codified_seq_noO"].duplicated()].reset_index(drop=True)

# Option 4: take largest N + random sample of the rest
# N_LARGEST = 25
# N_RANDOM  = 100
# tmp = test_locus.assign(allele_len=test_locus["Allele"].str.len())
# top = tmp.nlargest(N_LARGEST, "allele_len")
# rest = tmp.drop(top.index)
# test_locus = pd.concat([top, rest.sample(n=min(N_RANDOM, len(rest)), random_state=42)]) \
#                 .drop(columns=["allele_len"]) \
#                 .reset_index(drop=True)

# Option 5: random subset only
# N_RANDOM = 100
# test_locus = test_locus.sample(n=min(N_RANDOM, len(test_locus)), random_state=42).reset_index(drop=True)


# ----------------------------
# Optional: build cost dict (only needed for DIST_METHOD="cost_dict")
# ----------------------------
### can adjust these numbers to penalize non-repetetive motifs
# base_costs = {"ins_min": 5, "ins_max": 30, "dup_min": 1, "dup_max": 2, "sub": 25}
# motifs_for_costs = [m for m in conv.keys() if m != "Other"]
# for m in motifs_for_costs:
#     motif_frequencies.setdefault(m, 0.0)
# cost_dict = create_cost_dict(motifs_for_costs, motif_frequencies, base_costs)
# converted_cost_dict = convert_cost_dict(cost_dict, conv)


DIST_METHOD = "cost"  # "cost_dict" | "cost" | "levenshtein"

seqs = test_locus["codified_seq"].tolist()
num_sequences = len(seqs)
distance_matrix = np.zeros((num_sequences, num_sequences), dtype=float)

def compute_distance(i, j, method=DIST_METHOD):
    seq1 = seqs[i]
    seq2 = seqs[j]

    if method == "cost_dict":
        # requires converted_cost_dict defined above
        d = min_cost_to_convert_exact(seq1, seq2, converted_cost_dict)
    elif method == "cost":
        d = min_cost_to_convert_exact(seq1, seq2)
    elif method == "levenshtein":
        d = levenshtein_distance(seq1, seq2)
    else:
        raise ValueError(f"Unknown DIST_METHOD: {method}")

    return i, j, d

pairs = [(i, j) for i in range(num_sequences) for j in range(i, num_sequences)]

results = Parallel(n_jobs=-1, backend="loky")(
    delayed(compute_distance)(i, j) for i, j in tqdm(pairs, desc="Computing distances")
)

for i, j, dist in results:
    distance_matrix[i, j] = dist
    if i != j:
        distance_matrix[j, i] = dist


qPalette = [
    "#aec7e8",  # light blue
    "#1f77b4",  # blue
    "#ff7f0e",  # orange
    "#ffbb78",  # light orange
    "#2ca02c",  # green
    "#98df8a",  # light green
    "#17becf",  # cyan
    "#d62728",  # red
    "#ff9896",  # light red
    "#9467bd",  # purple
    "#c5b0d5",  # light purple
    "#8c564b",  # brown
    "#c49c94",  # light brown
    "#e377c2",  # pink
    "#f7b6d2",  # light pink
    "#7f7f7f",  # gray
    "#c7c7c7",  # light gray
    "#bcbd22",  # olive
    "#dbdb8d",  # light olive
    "#9edae5",  # light cyan
]

motif_to_color = {code: color for code, color in zip(conv.values(), qPalette)}
motif_to_color["O"] = "gray"


max_codified_length = calculate_max_codified_length(df['encodedSeq'], conv)
#tree = create_tree_from_distance_matrix(distance_matrix, method='')
tree = create_tree_from_distance_matrix(distance_matrix, method='upgma')
ts = TreeStyle()
#sort_by_height(tree)
#reverse_tree(tree)
#tree.ladderize(direction=0)
ts.branch_vertical_margin = 1
ts.show_branch_length = False
ts.show_scale = False
ts.mode = "r"
ts.margin_left = 10
ts.margin_right = 10
ts.margin_top = 10
ts.margin_bottom = 10
ts.show_leaf_name = False
ts.scale = 1

# Render the tree with colored sequences
tree = render_tree_with_colored_sequences(tree, ts, df, motif_to_color, max_codified_length, conv, simple=True)
tree_img_path = render_tree_to_tempfile(tree, ts)
#tree.render(file_name='%%inline', w=800, h=1200, tree_style=ts)
tree.render(file_name='/Users/isaacxu/TandemRepeats/data/example_tree.png', w=400, h=800, tree_style=ts, dpi=300)

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import os

# invert conv: code -> motif
code_to_motif = {code: motif for motif, code in conv.items()}

handles = [
    Patch(
        facecolor=color,
        edgecolor="black",
        label=code_to_motif.get(code, code)
    )
    for code, color in motif_to_color.items()
    if code != "O"   # optional: skip "Other" if you want
]

fig, ax = plt.subplots(figsize=(4, max(2, 0.3 * len(handles))))
ax.axis("off")
ax.legend(handles=handles, loc="upper left", frameon=False)

OUT_LEGEND = os.path.join(
    os.path.dirname(__file__),
    "..",
    "data",
    "motif_color_legend.png"
)

plt.tight_layout()
plt.savefig(OUT_LEGEND, dpi=300)
plt.close()



