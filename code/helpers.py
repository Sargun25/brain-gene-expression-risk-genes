import numpy as np, pandas as pd, re

project = '/content/drive/MyDrive/hd_project'
expression = pd.read_csv(f'{project}/expression_dk_named.csv', index_col=0)
expr_L = expression.loc[[r for r in expression.index if r.startswith(('L_', 'B_'))]]
X = expr_L.values
cols = np.array(expr_L.columns)
alias = {'CIAO2A': 'FAM96A'}

def present(genes):
    out = []
    for g in genes:
        if g in expr_L.columns:
            out.append(g)
        elif g in alias and alias[g] in expr_L.columns:
            out.append(alias[g])
    return sorted(set(out))

def gene_list(df):
    s = df['MAPPED_GENE'].dropna()
    parts = s.apply(lambda x: re.split(r'\s+-\s+|\s*,\s*|\s*;\s*|\s+x\s+', x))
    return sorted({g.strip() for lst in parts for g in lst if g.strip()})

def score(genes, n_perm=10000, seed=0):
    from statsmodels.stats.multitest import multipletests
    sel = present(genes)
    k = len(sel)
    sel_idx = np.where(np.isin(cols, sel))[0]
    pool = np.setdiff1d(np.arange(len(cols)), sel_idx)
    rng = np.random.default_rng(seed)
    obs = X[:, sel_idx].mean(axis=1)
    null = np.empty((n_perm, X.shape[0]))
    for i in range(n_perm):
        null[i] = X[:, rng.choice(pool, size=k, replace=False)].mean(axis=1)
    z = (obs - null.mean(0)) / null.std(0)
    p = ((null >= obs).sum(0) + 1) / (n_perm + 1)
    q = multipletests(p, method='fdr_bh')[1]
    return pd.DataFrame({'z': z, 'p': p, 'q_fdr': q}, index=expr_L.index), k

def score_roi(genes, rois, n_perm=10000, seed=0):
    sel = present(genes)
    k = len(sel)
    sel_idx = np.where(np.isin(cols, sel))[0]
    pool = np.setdiff1d(np.arange(len(cols)), sel_idx)
    ridx = [list(expr_L.index).index(r) for r in rois]
    Xr = X[ridx]
    rng = np.random.default_rng(seed)
    obs = Xr[:, sel_idx].mean()
    null = np.array([Xr[:, rng.choice(pool, size=k, replace=False)].mean() for _ in range(n_perm)])
    return round(float((obs - null.mean()) / null.std()), 3), round(float(((null >= obs).sum() + 1) / (n_perm + 1)), 4)
