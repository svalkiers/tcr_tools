import pandas as pd

from os.path import dirname, abspath, join
from modules.tcrdist1.all_genes import all_genes # for recognized genes
from modules.tcrdist1.translation import get_translation
from modules.tcrdist1.tcr_sampler import get_j_cdr3_nucseq

ROOT = dirname(dirname(dirname(abspath(__file__))))
DATA = join(ROOT, '/data')

IMGT = pd.read_csv(join(DATA,'imgt_reference.tsv'), sep='\t')
mapping = pd.read_csv(join(DATA, 'adaptive_imgt_mapping.csv'))

adaptive_to_imgt_human = mapping.loc[mapping['species'] == 'human'].set_index('adaptive')['imgt'].fillna('NA').to_dict()
adaptive_to_imgt_mouse = mapping.loc[mapping['species'] == 'mouse'].set_index('adaptive')['imgt'].fillna('NA').to_dict()

v_fam_freq = pd.read_csv(join(DATA,"adaptive_v_fam_to_imgt_gene.txt"), sep="\t")
adaptive_vfam_mapping = dict(zip(v_fam_freq.adaptive_v_family, v_fam_freq.imgt_v_allele))

aminoacids = 'ACDEFGHIKLMNPQRSTVWY'
_aminoacids_set = set(aminoacids)

def _is_aaseq(seq:str):
    """
    Check if string contains non-amino acid characters.
    Returns True if string only contains standard amino acid characters.
    """
    try:
        return all(c in _aminoacids_set for c in seq)
    except TypeError:
        return False

def _is_cdr3(seq:str):
    """
    Checks if string is a valid CDR3 amino acid sequence,
    according to the following defenitions:
        - First amino acid character is C.
        - Last amino acid is F or W.
        - Sequence exclusively contains valid amino acid characters.
    """
    try:
        return (_is_aaseq(seq)
            and (seq[0] == 'C')
            and (seq[-1] in ['F', 'W', 'C'])
            and (len(seq) <= 30)
            and (len(seq) >= 6))
    # Exclude non-string type input
    except TypeError:
        return Falses

def capture_nucseq(tcr, organism='human'):
    """
    Identifies CDR3 nucleotide sequence from full TCR rearrangement.
    """
    vgene, cdr3, rearrangement, jgene = tcr

    # figure out the cdr3 nucseq
    cdr3_nucseq = None
    jgene_cdr3_nucseq = get_j_cdr3_nucseq(organism,jgene).upper() # added .upper()

    for offset in range(3):
        protseq = get_translation( rearrangement, '+{}'.format(offset+1) )
        for ctrim in range(1,4):
            if cdr3[:-ctrim] in protseq:
                start = offset + 3*(protseq.index(cdr3[:-ctrim]))
                length = 3*(len(cdr3)-ctrim)
                cdr3_nucseq = rearrangement[ start : start + length ]
                if ctrim:
                    cdr3_nucseq += jgene_cdr3_nucseq[-3*ctrim:]
                if cdr3 != get_translation( cdr3_nucseq, '+1' ):
                    cdr3_nucseq = ''
                break
        if cdr3_nucseq=='': # failure signal
            cdr3_nucseq = None
            break

    if cdr3_nucseq is None:
        print('parse cdr3_nucseq failed:',tcr)
        return None
    else:
        return cdr3_nucseq

def select_highest_duplicate_count(group):
    assert "chain" in group.columns, "DataFrame must contain [chain] column."
    if 'A' not in group['chain'].values or 'B' not in group['chain'].values:
        return pd.DataFrame()  # Return an empty DataFrame to ignore this group
    a_max_idx = group[group['chain'] == 'A']['duplicate_count'].idxmax()
    b_max_idx = group[group['chain'] == 'B']['duplicate_count'].idxmax()
    return group.loc[[a_max_idx, b_max_idx]]