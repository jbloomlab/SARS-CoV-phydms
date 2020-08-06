"""
Utils phydms.
"""

from ete3 import PhyloTree, TreeStyle, NodeStyle
import pandas as pd


def create_nodestyle(color):
    ns = NodeStyle()
    ns["bgcolor"] = color
    return ns


def viztree(fname, clade_dictionary):
    with open(fname) as f:
        treestring = f.read()
    t = PhyloTree(treestring)
    # get node names
    node_names = [x.name for x in t.iter_descendants("postorder") if x.name]

    # color the clades
    for key in clade_dictionary:
        specs = clade_dictionary[key]
        present = all([x in node_names for x in specs[0]])
        if present:
            (t.get_common_ancestor(specs[0])
              .set_style(create_nodestyle(specs[1])))
    return t


def omegabysiteresults(models, prefix, sig_cutoff):
    for m in models:
        df = pd.read_csv(f"{prefix}/{m}_omegabysite.txt",
                         sep='\t',
                         comment='#',)
        df = df[df["Q"] < sig_cutoff]
        print(f"Model {m}")
        print(f"There are {len(df[df.P > 1.0])} sites with sig "
              f"(Q < {sig_cutoff}) evidence for omega > 1.")
        print(f"There are {len(df[df.P < 1.0])} sites with sig "
              f"(Q < {sig_cutoff}) evidence for omega < 1.")
        print(df)
        print()


def relresults(models, prefix, sig_cutoff):
    for m in models:
        if "gammaomega" in m:
            df = pd.read_csv(f"{prefix}/"
                             f"{m}_posteriorprobabilities.csv")
            df = df[df["fdr"] < 0.05]
            print(f"Model {m}")
            print(f"There are {len(df[df['p(omega > 1)'] > 1.0])} sites with "
                  f"sig (fdr < {sig_cutoff}) evidence for omega > 1.")
            if len(df) > 0:
                print(df)
            print()
