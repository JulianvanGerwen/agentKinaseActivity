import sys, pickle, numpy as np, pandas as pd
sys.path.insert(0,'/workspace/src')
from kinact.explore.lib_variants import DevData, Engine, sums_measured, stat_from_sums, pred_weights
OUT='/workspace/cache/audit_r2np/'; import os; os.makedirs(OUT, exist_ok=True)
dd=DevData(); E=Engine(dd)
cfgs={
 'r2np200g':dict(stat='wz',lib=dict(type='pred',name='r2np',k=200,mode='global',gamma=0),min_sub=1),
 'r2np50m':dict(stat='wz',lib=dict(type='pred',name='r2np',k=50,mode='measured',gamma=1),min_sub=1),
 'r3_200m':dict(stat='wz',lib=dict(type='pred',name='r3',k=200,mode='measured',gamma=4),min_sub=1),
 'cant100m':dict(stat='wz',lib=dict(type='pred',name='cantley',k=100,mode='measured',gamma=0),min_sub=1),
 'cant50m':dict(stat='wz',lib=dict(type='pred',name='cantley',k=50,mode='measured',gamma=1),min_sub=1),
 'psp':dict(stat='wz',lib=dict(type='static',name='psp_all'),min_sub=1),
 'cur':dict(stat='wz',lib=dict(type='static',name='curated_union'),min_sub=1),
}
scores={n:E.score(c) for n,c in cfgs.items()}
# r2np measured top-50 excluding the kinase's own curated substrates (curated_union) -> "non-PSP" part
for nm,lib,k,g in [('r2np50m_noncur','r2np',50,1),('cant50m_noncur','cantley',50,1),('r3_50m_noncur','r3',50,1)]:
    ov=E.overlap('curated_union',lib)
    # need sel with exclusion then take first k valid -> approximate: exclude within top-KMAX and keep first k remaining
    sel=E.sel(lib).copy(); sel[ov]=-1
    # compact
    order=np.argsort(sel<0,axis=2,kind='stable'); sel2=np.take_along_axis(sel,order,2)
    S=sums_measured(dd,sel2,k,g); scores[nm]=stat_from_sums(dd,S,'wz',1)
    # fraction of top-k that were curated
    frac=ov[:,:,:k].mean(); print(nm,'frac top-k curated',frac)
    if nm.startswith('r2np'): pickle.dump(ov,open(OUT+'ov_r2np.pkl','wb'))
res={}
for n,s in scores.items():
    P=dd.ev.pairs(s); P['cfg']=n; res[n]=P
    m=dd.ev.metric(P,'M1-sel'); bs={src:dd.ev.metric(P,f'M1-sel:source={src}')['value'] for src in 'DHO'}
    print(f"{n:16s} M1-sel {m['value']:.3f}  D {bs['D']:.3f} H {bs['H']:.3f} O {bs['O']:.3f}")
pd.concat(res.values()).to_parquet(OUT+'pairs.parquet')
