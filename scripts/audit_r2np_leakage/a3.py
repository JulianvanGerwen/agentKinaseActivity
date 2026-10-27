import sys, pickle, numpy as np, pandas as pd
sys.path.insert(0,'/workspace/src')
from kinact.explore.lib_variants import DevData, Engine, sums_measured, stat_from_sums, measured_topk
O='/workspace/cache/audit_r2np/'; import os; os.makedirs(O, exist_ok=True)
dd=DevData(); E=Engine(dd)
kin=pd.read_parquet('/workspace/external/ksm_index/kinases.parquet')
U=dd.U.set_index('kinase')
acc2rc=dict(zip(kin.uniprot_acc,kin.residue_class))
rc=np.array([acc2rc.get(U.loc[k,'uniprot'],'ST') for k in dd.kin])
res=np.array([s.split('_')[1][0] for s in dd.sites])
isY=res=='Y'
uni=pd.read_parquet('/workspace/external/ksm_index/universe.parquet',columns=['site_id','n_sources','comp_n_conditions','src_compendium'])
from kinact.library import SiteMapper
sm=SiteMapper(); keys=sm.from_keys(uni.site_id)
uni.index=keys; uni=uni[~uni.index.duplicated()]
nsrc=uni['n_sources'].reindex(dd.sites).fillna(0).to_numpy(float)
ncond=dd.M.sum(1).astype(float)  # DEV conditions measured
rng=np.random.default_rng(0)
def null_D(score_site):
    # kinase-agnostic per residue-class score matrix (K x S) -> class-matched site ranking
    D=np.full((len(dd.kin),len(dd.sites)),np.nan,np.float32)
    for i,c in enumerate(rc):
        ok = isY if c=='Y' else (~isY if c=='ST' else np.ones_like(isY))
        D[i,ok]=score_site[ok]+rng.normal(0,1e-6,ok.sum())  # tiny jitter: break ties differently per kinase
    return D
nulls={'null_resclass_random':null_D(rng.random(len(dd.sites))),
       'null_popularity_nsrc':null_D(nsrc+rng.random(len(dd.sites))*0.01),
       'null_popularity_devcond':null_D(ncond+rng.random(len(dd.sites))*0.01)}
scores={}
for nm,D in nulls.items():
    sel=measured_topk(dd,D,200)
    for k in (50,200):
        scores[f'{nm}@{k}']=stat_from_sums(dd,sums_measured(dd,sel,k,0),'wz',1)
for nm,lib,k,g in [('r2np50m','r2np',50,1),('cant50m','cantley',50,1),('r3_50m','r3',50,1)]:
    scores[nm]=stat_from_sums(dd,sums_measured(dd,E.sel(lib),k,g),'wz',1)
scores['psp']=E.score(dict(stat='wz',lib=dict(type='static',name='psp_all'),min_sub=1))
# within-residue-class AUROC: mask kinases of other class relative to target -> use auroc_scored with mask per class
out=[]
for nm,s in scores.items():
    P=dd.ev.pairs(s); P['cfg']=nm
    S=dd.ev._matrix(s)
    for cl in ['Y','ST']:
        mask=np.broadcast_to((rc==cl)[None,:],S.shape)
        Pm=dd.ev.pairs(s,mask=mask)
        P[f'auroc_within_{cl}']=Pm['auroc_scored'].values
    out.append(P)
P=pd.concat(out); P['trc']=P.kinase.map(dict(zip(dd.kin,rc)))
P.to_parquet(O+'pairs_nulls.parquet')
P=P[P.in_M1sel]
def macro(df,col='auroc'):
    return df.groupby('kinase')[col].mean().mean()
rows=[]
for nm,g in P.groupby('cfg',sort=False):
    r={'cfg':nm,'M1sel':macro(g)}
    for s in 'DHO':
        r[s]=macro(g[g.source==s])
        r[s+'_Ytgt']=macro(g[(g.source==s)&(g.trc=='Y')])
        r[s+'_STtgt']=macro(g[(g.source==s)&(g.trc=='ST')])
        gy=g[(g.source==s)&(g.trc=='Y')]; r[s+'_Ytgt_withinY']=macro(gy,'auroc_within_Y')
        gs=g[(g.source==s)&(g.trc=='ST')]; r[s+'_STtgt_withinST']=macro(gs,'auroc_within_ST')
    rows.append(r)
T=pd.DataFrame(rows).set_index('cfg'); pd.set_option('display.width',250)
print(T.round(3).to_string())
T.to_csv(O+'nulls_table.tsv',sep='\t')
print(P[P.cfg=='psp'].groupby(['source','trc']).agg(n=('kinase','size'),nk=('kinase','nunique')))
