import argparse
parser = argparse.ArgumentParser(prog='prepare_pile_in_region.py', description='''
    Merge count by region
''')
parser.add_argument('--pilein')
parser.add_argument('--output')
args = parser.parse_args()

def rename_df(df):
    new_col = {}
    for i in df.columns:
        new_col[i] = 'V' + str(i + 1)
    df = df.rename(columns = new_col)
    return df

def convert_signal_to_score(group, start, end, score):
    out = []
    dup_ind = group.duplicated(subset = start)
    group_unique = group[dup_ind == False]
    for i in range(group_unique.shape[0]):
        line = group_unique.iloc[i]
        n = line[end] - line[start]
        out = out + [ str(line[score]) for j in range(n) ]
    out = ','.join(out)
    return out

import pandas as pd
df = pd.read_table(args.pilein, header = None, compression = 'gzip')
df = rename_df(df)
df = df.groupby(['V5', 'V6', 'V7', 'V8', 'V11']).apply(convert_signal_to_score, 'V2', 'V3', 'V4')
df = pd.DataFrame(df)
df.to_csv(args.output, sep = '\t', compression = 'gzip', header = False)
