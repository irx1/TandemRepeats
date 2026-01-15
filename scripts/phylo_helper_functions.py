from collections import Counter
import re

def build_conv(filtered_motifs, max_motifs=None, include_motifs=None, overrides=None):
    alpha = list('ABCDEFGHIJKLMNPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz')
    if max_motifs is not None:
        filtered_motifs = filtered_motifs[:max_motifs]

    conv = dict(zip(filtered_motifs, alpha[:len(filtered_motifs)]))

    include_motifs = include_motifs or []
    overrides = overrides or {}

    used = set(conv.values())
    avail = [c for c in alpha if c not in used]

    for m in include_motifs:
        if m not in conv:
            if not avail:
                raise ValueError("Ran out of encoding letters.")
            conv[m] = avail.pop(0)

    for m, letter in overrides.items():
        conv[m] = letter

    conv["Other"] = "O"
    return conv


def calculate_motif_frequencies_and_filter(encoded_seqs, min_count=1):
    """
    Calculate allele-based motif frequencies and return motifs filtered by count.

    Parameters:
    - encoded_seqs (list): List of encoded sequences with motifs.
    - min_count (int): Minimum count of each motif in an allele for it to be considered.

    Returns:
    - dict: Motif frequencies based on allele presence.
    - int: Total number of alleles in the input.
    - list: Filtered list of motifs that meet the minimum count threshold, sorted by allele frequency.
    """
    # Initialize counters
    motif_allele_counter = Counter()
    motif_counter = Counter()
    total_alleles = 0
    
    # Process each allele (encoded sequence)
    for seq in encoded_seqs:
        if not isinstance(seq, str):
            continue  # Skip non-string entries
        total_alleles += 1
        # Extract motifs and their counts
        motif_pattern = re.compile(r'<([^>]+)>(\d+)')
        motifs = motif_pattern.findall(seq)
        
        # Track motifs that meet the minimum count within this allele
        unique_motifs_in_allele = set()
        for motif, count in motifs:
            if int(count) >= min_count:
                #motif = simplest_motif(motif)
                motif_counter[motif] += int(count)
                unique_motifs_in_allele.add(motif)
        # Update counter for each unique motif in this allele
        motif_allele_counter.update(unique_motifs_in_allele)
    
    # Calculate motif frequencies for each motif
    motif_frequencies = {motif: count / total_alleles for motif, count in motif_allele_counter.items()}
    filtered_motifs = [motif for motif, count in motif_counter.most_common() if count >= min_count]
    
    
    return motif_frequencies, total_alleles, filtered_motifs

def replace_encoded_seq(encoded_seq, kmer_dict):
    """
    Replace encoded motifs with single-character symbols.

    Example:
        <CAG>10<CAA>3  ->  AAAAAAAAAABBB
    """
    matches = re.findall(r'<([^>]+)>(\d+)', encoded_seq)
    result = []

    other = kmer_dict.get('Other', '')

    for motif, count in matches:
        symbol = kmer_dict.get(motif, other)
        result.append(symbol * int(count))

    return ''.join(result)


def replace_encoded_seq_simple(
    encoded_seq,
    kmer_dict,
    min_count=1,
    min_motif_len=1
):
    """
    Simplified encoding that drops small / low-support motif runs.

    Keeps a motif block if:
      - count >= min_count
      - and motif length >= min_motif_len
    """
    matches = re.findall(r'<([^>]+)>(\d+)', encoded_seq)
    result = []

    for motif, count in matches:
        count = int(count)

        if count < min_count:
            continue
        if len(motif) < min_motif_len:
            continue

        result.append(
            kmer_dict.get(motif, kmer_dict.get('Other', '')) * count
        )

    return ''.join(result)


def calculate_cost(value, min_value, max_value, min_cost, max_cost):
    """
    Scale a value (e.g., allele frequency) to a cost within a given range.
    
    Parameters:
    - value (float): The input value to scale (e.g., allele frequency).
    - min_value (float): Minimum possible value (e.g., smallest allele frequency).
    - max_value (float): Maximum possible value (e.g., highest allele frequency).
    - min_cost (float): Minimum cost in the output range.
    - max_cost (float): Maximum cost in the output range.
    
    Returns:
    - float: Scaled cost value.
    """
    scaled_cost = max_cost - (value - min_value) / (max_value - min_value) * (max_cost - min_cost)
    return scaled_cost

def create_cost_dict(motifs, motif_frequencies, base_costs=None):
    """
    Creates a cost dictionary for motifs with insertion, deletion, duplication, and contraction
    costs based on allele frequencies. Single-letter motifs are assigned the same costs as 'Other'.

    Parameters:
    - motifs (list): List of motifs.
    - motif_frequencies (dict): Dictionary of motifs and their allele frequencies.
    - base_costs (dict): Base costs for operations, including min and max values for scaling.

    Returns:
    - dict: Cost dictionary for each motif with operation-specific costs.
    """
    # Default base cost parameters if none provided
    if base_costs is None:
        base_costs = {
            'ins_min': 50, 'ins_max': 100,  # Range for insertion and deletion
            'dup_min': 1, 'dup_max': 5,    # Range for duplication and contraction
            'sub': 10                      # Substitution cost scaling factor
        }

    cost_dict = defaultdict(dict)

    # Define "Other" default costs
    other_costs = {
        'ins': 2,
        'del': 2,
        'dup': 1,
        'contract': 1,
        'sub': 2  # Default substitution cost for "O" motif
    }
    cost_dict['O'] = other_costs

    # Find min and max allele frequencies for scaling
    min_freq = min(motif_frequencies.values())
    max_freq = max(motif_frequencies.values())

    for motif in motifs:
        if len(motif) == 1:  # Assign single-letter motifs the same costs as "Other"
            cost_dict[motif] = other_costs
        else:
            allele_frequency = motif_frequencies.get(motif, min_freq)  # Use min_freq if motif not in frequency data

            # Calculate insertion and deletion costs, scaled between 5 and 100
            ins_cost = calculate_cost(allele_frequency, min_freq, max_freq, base_costs['ins_min'], base_costs['ins_max'])
            del_cost = ins_cost  # Use the same scaling for deletion

            # Calculate duplication and contraction costs, scaled between 1 and 5
            dup_cost = calculate_cost(allele_frequency, min_freq, max_freq, base_costs['dup_min'], base_costs['dup_max'])

            # Populate the cost dictionary for multi-letter motifs
            cost_dict[motif]['ins'] = ins_cost
            cost_dict[motif]['del'] = del_cost
            cost_dict[motif]['dup'] = dup_cost
            cost_dict[motif]['contract'] = dup_cost

    # Populate substitution costs based on Levenshtein distance for multi-letter motifs only
    for i, motif1 in enumerate(motifs):
        for motif2 in motifs[i + 1:]:
            if len(motif1) > 1 and len(motif2) > 1:  # Exclude single-letter motifs
                substitution_cost = levenshtein_distance(motif1, motif2) * base_costs['sub']
                cost_dict[motif1][f'sub_{motif2}'] = substitution_cost
                cost_dict[motif2][f'sub_{motif1}'] = substitution_cost  # Symmetric cost

    return cost_dict


def convert_cost_dict(cost_dict, conv):
    """
    Converts the keys in cost_dict to single-character codes based on the conv mapping.

    Parameters:
    - cost_dict (dict): Original cost dictionary with motif names as keys.
    - conv (dict): Dictionary mapping motif names to single-character codes.

    Returns:
    - dict: Updated cost dictionary with single-character codes as keys.
    """
    # Initialize new cost dictionary
    new_cost_dict = {}

    for motif, costs in cost_dict.items():
        # Get the single-character code for the motif
        char_code = conv.get(motif, 'O')  # Default to 'O' if motif not in conv
        new_cost_dict[char_code] = {}
        
        for cost_type, cost_value in costs.items():
            # If cost_type is substitution, convert the target motif to its character code
            if cost_type.startswith('sub_'):
                target_motif = cost_type.split('_')[1]
                target_char = conv.get(target_motif, 'O')
                new_cost_dict[char_code][f'sub_{target_char}'] = cost_value
            else:
                # For non-substitution types, keep the cost type as is
                new_cost_dict[char_code][cost_type] = cost_value

    return new_cost_dict

def min_cost_to_convert_exact(s, t, cost_dict=None, ins_cost=10, del_cost=10, sub_cost=10, dup_cost=1, contract_cost=1):
    """
    Computes the minimum cost to convert string s to string t using a cost dictionary
    or default costs for edit operations including duplication and contraction.
    Now includes a look-ahead for interrupted duplication patterns.

    Parameters:
    - s (str): Source codified sequence (e.g., motifs).
    - t (str): Target codified sequence.
    - cost_dict (dict, optional): Dictionary with motif-specific operation costs.
    - ins_cost, del_cost, sub_cost, dup_cost, contract_cost (int): Operation costs.

    Returns:
    - int: Minimum cost to convert s to t.
    """
    n, m = len(s), len(t)
    dp = [[float('inf')] * (m + 1) for _ in range(n + 1)]

    # Base cases
    dp[0][0] = 0
    for i in range(1, n + 1):
        motif = s[i - 1]
        dp[i][0] = dp[i - 1][0] + (cost_dict[motif]['del'] if cost_dict and motif in cost_dict else del_cost)
    for j in range(1, m + 1):
        motif = t[j - 1]
        dp[0][j] = dp[0][j - 1] + (cost_dict[motif]['ins'] if cost_dict and motif in cost_dict else ins_cost)

    # Fill DP table
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            motif_s = s[i - 1]
            motif_t = t[j - 1]

            # Substitution or match
            if motif_s == motif_t:
                dp[i][j] = dp[i - 1][j - 1]
            else:
                sub_key = f'sub_{motif_t}'
                sub_c = cost_dict[motif_s].get(sub_key, sub_cost) if cost_dict and motif_s in cost_dict else sub_cost
                dp[i][j] = dp[i - 1][j - 1] + sub_c

            # Insertion
            ins_c = cost_dict[motif_t]['ins'] if cost_dict and motif_t in cost_dict else ins_cost
            dp[i][j] = min(dp[i][j], dp[i][j - 1] + ins_c)

            # Deletion
            del_c = cost_dict[motif_s]['del'] if cost_dict and motif_s in cost_dict else del_cost
            dp[i][j] = min(dp[i][j], dp[i - 1][j] + del_c)

            # Duplication (simple)
            if j > 1 and t[j - 1] == t[j - 2]:
                dup_c = cost_dict[motif_t]['dup'] if cost_dict and motif_t in cost_dict else dup_cost
                dp[i][j] = min(dp[i][j], dp[i][j - 1] + dup_c)

            # Advanced duplication with interruption:
            # Pattern: motif motif different motif (e.g., A A B) → duplicate A twice, then substitute middle A to B
            if j > 2 and t[j - 3] == t[j - 1] and t[j - 2] != t[j - 3]:
                motif = t[j - 3]
                dup_c = cost_dict[motif].get('dup', dup_cost) if cost_dict and motif in cost_dict else dup_cost
                sub_c = cost_dict[motif].get(f'sub_{t[j - 2]}', sub_cost) if cost_dict and motif in cost_dict else sub_cost
                total_cost = dp[i][j - 2] + dup_c  + ins_c  # 2 duplications + 1 insertion
                dp[i][j] = min(dp[i][j], total_cost)


            # Contraction
            if i > 1 and s[i - 1] == s[i - 2]:
                contract_c = cost_dict[motif_s]['contract'] if cost_dict and motif_s in cost_dict else contract_cost
                dp[i][j] = min(dp[i][j], dp[i - 1][j] + contract_c)

            # Advanced contraction with interruption (reverse logic of duplication)
            if i > 2 and s[i - 3] == s[i - 1] and s[i - 2] != s[i - 3]:
                motif = s[i - 3]
                contract_c = cost_dict[motif].get('contract', contract_cost) if cost_dict and motif in cost_dict else contract_cost
                sub_c = cost_dict[motif].get(f'sub_{s[i - 2]}', sub_cost) if cost_dict and motif in cost_dict else sub_cost
                total_cost = dp[i - 2][j] + contract_c  + del_c
                dp[i][j] = min(dp[i][j], total_cost)

    return dp[n][m]



def calculate_max_codified_length(encoded_seqs, conv):
    """
    Calculate the maximum codified length from the encoded sequences.
    
    Parameters:
    - encoded_seqs (pd.Series): Series containing encoded sequences like '<AAAAG>10<AAAGG>5'.
    - conv (dict): Dictionary mapping full motifs to their codified single-letter representation.

    Returns:
    - int: Maximum codified length across all sequences.
    """
    max_length = 0
    
    for encoded_seq in encoded_seqs:
        # Parse each motif and its count in the encoded sequence
        motif_matches = re.findall(r'<([^>]+)>(\d+)', encoded_seq)
        total_length = sum(len(motif) * int(count) for motif, count in motif_matches if motif in conv)
        
        # Update max length if the current total length is greater
        max_length = max(max_length, total_length)
    
    return max_length


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

def create_tree_from_distance_matrix(distance_matrix, method='upgma'):
    if method == 'upgma':
        t = upgma(distance_matrix).to_newick(include_distance=True)
    else:
        t = neighjoin(distance_matrix).to_newick(include_distance=True)
    return Tree(t)

# Step 2: Render tree with codified sequence images
def render_tree_with_colored_sequences(tree, ts, t2, motif_to_color, max_codified_length, conv, simple=True):
    nstyle = NodeStyle()
    nstyle["vt_line_width"] = 0
    nstyle["hz_line_width"] = 0
    nstyle["shape"] = "sphere"
    nstyle["size"] = 0
    nstyle["fgcolor"] = "darkred"

    for n in tree.traverse():
        n.set_style(nstyle)

    for leaf in tree.iter_leaves():
        if simple==False:
            codified_seq = t2.loc[int(leaf.name), 'codified_seq']
            #codified_seq = t2.loc[int(leaf.name), 'codified_seq_noO']
        else:
            codified_seq = t2.loc[int(leaf.name), 'codified_seq_simple']
        img = generate_colored_line_codified(codified_seq, motif_to_color, max_codified_length, conv)
        img_path = pil_image_to_temp_file(img)
        face = ImgFace(img_path)
        #leaf.add_face(face, column=0, position="aligned")
        leaf.add_face(face, column=1, position="branch-right")

    return tree

# Step 3: Generate codified sequence images with proportional widths
def generate_colored_line_codified(codified_seq, motif_to_color, max_codified_length, conv, base_box_width=1, box_height=10):

    
    img_width = max_codified_length * base_box_width
    img = Image.new('RGB', (img_width, box_height), 'white')
    draw = ImageDraw.Draw(img)
    x_offset = 0
    for letter in codified_seq:
        full_motif = next((motif for motif, code in conv.items() if code == letter), letter)
        motif_length = len(full_motif)
        box_width = base_box_width * motif_length
        color = motif_to_color.get(letter, "gray")
        draw.rectangle([x_offset, 0, x_offset + box_width, box_height], fill=color)
        x_offset += box_width

    return img

# Helper function to save a PIL image to a temporary file
def pil_image_to_temp_file(img):
    temp_file = tempfile.NamedTemporaryFile(delete=False, suffix=".png")
    img.save(temp_file.name, "PNG")
    return temp_file.name

# Step 4: Render and save the final tree image
def render_tree_to_tempfile(tree, ts):
    temp_file = tempfile.NamedTemporaryFile(delete=False, suffix=".png")
    tree.render(temp_file.name, w=800, h=1000, tree_style=ts)
    return temp_file.name


def reverse_tree(tree):
    """Reverse the order of children for each node in the tree."""
    for node in tree.traverse():
        node.children = node.children[::-1]

def sort_by_height(tree):
    """Sort the children of each node by their height."""
    for node in tree.traverse():
        node.children.sort(key=lambda child: child.get_distance(tree), reverse=True) 
