import pandas as pd

from parsing_functions import adaptive_to_imgt_human, adaptive_to_imgt_mouse, adaptive_vfam_mapping, IMGT

def read_adaptive(path, organism='human', old=False, recover_unresolved=True, allele_level=False):
    # Read file
    df = pd.read_csv(path, sep="\t", low_memory=False)
    # Remove out-of-frame sequences
    df = df[df.frame_type=="In"]
    # Desired columns
    if old:
        prefix = "cdr3_"
    else:
        prefix = ""
    cols = [
        'templates',
        f'{prefix}rearrangement', 
        f'{prefix}amino_acid',
        'v_family',
        'j_family',
        'v_gene',
        'v_allele',
        'j_gene',
        'j_allele'
        ]
    df = df[cols]
    # Gene parsing
    unresolved = df[df.v_gene=="unresolved"]
    df["v_call"] = df.v_gene.map(adaptive_to_imgt_human)
    df["j_call"] = df.j_gene.map(adaptive_to_imgt_human)
    df = df.dropna(subset=["v_call","j_call"])
    # If True, recover V genes by imputing the most common allele within the V family
    if recover_unresolved:
        unresolved["v_call"] = unresolved.v_family.map(adaptive_vfam_mapping)
        unresolved["j_call"] = unresolved.j_gene.map(adaptive_to_imgt_human)
        print(f"Recovered {unresolved.v_imgt.dropna().shape[0]} V genes")
        df = pd.concat([df,unresolved])
    df[["v_allele","j_allele"]] = df[["v_allele","j_allele"]].fillna(0)
    # # df["v_call"] = df.apply(lambda v: v.v_imgt.split("*")[0] + f"*0{int(v.v_allele)}" if v.v_allele!=0 else v, axis=1)
    # # df["j_call"] = df.apply(lambda j: j.j_imgt.split("*")[0] + f"*0{int(j.j_allele)}" if j.j_allele!=0 else j, axis=1)
    df = df.dropna(subset=["v_imgt"])
    # if allele_level:
        # df["v_call"] = [df.v_imgt.iloc[n].split("*")[0]+f"*0{int(i)}" if int(i)!=0 else df.v_imgt.iloc[n] for n,i in enumerate(df.v_allele)]
        # df["j_call"] = df.j_imgt
    functional = IMGT[IMGT['fct']=='F']
    df = df[df.v_call.isin(functional.imgt_allele_name)]
    df = df[['templates',f'{prefix}rearrangement',f'{prefix}amino_acid','v_call','j_call']]
    df = df[df[f"{prefix}amino_acid"].apply(lambda cdr3: _is_cdr3(cdr3))]
    df = df.dropna(subset=[f'{prefix}rearrangement',f'{prefix}amino_acid','v_call','j_call'])
    df = df.sort_values(by="templates", ascending=False)
    df = df.rename(columns={f"{prefix}rearrangement":"junction",f"{prefix}amino_acid":"junction_aa"})
    df = df.drop_duplicates()
    return df.reset_index(drop=True)