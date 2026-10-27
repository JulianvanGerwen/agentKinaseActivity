# Response enrichment: fraction of a kinase's measured top-50 sites (non-curated) that are among the condition's
# most-changed 5% sites in the GT-sign direction; target vs non-target (same residue class), per source.
import sys, numpy as np, pandas as pd
sys.path.insert(0,'/workspace/src')
from kinact.explore.lib_variants import DevData, Engine, measured_topk
O='/workspace/cache/audit_r2np/'; import os; os.makedirs(O, exist_ok=True)
dd=DevData(); E=Engine(dd)
kin=pd.read_parquet('/workspace/external/ksm_index/kinases.parquet'); acc2rc=dict(zip(kin.uniprot_acc,kin.residue_class))
U=dd.U.set_index('kinase'); rc=np.array([acc2rc.get(U.loc[k,'uniprot'],'ST') for k in dd.kin])
tg=dd.ev.targets; tg=tg[tg.in_M1sel]
C=len(dd.conds); V=np.where(dd.M,dd.V0,np.nan)
# per condition: down-set (bottom 5%) and up-set (top 5%)
lo=np.nanquantile(V,0.05,axis=0); hi=np.nanquantile(V,0.95,axis=0)
DOWN=(V<=lo[None,:])&dd.M; UP=(V>=hi[None,:])&dd.M
condpos={c:i for i,c in enumerate(dd.conds)}
rows=[]
k=50
for lib in ['r2np','cantley','r3']:
    sel=E.sel(lib)[:,:,:k]
    ov=E.overlap('curated_union',lib)[:,:,:k]
    for _,t in tg.iterrows():
        c=condpos[t.conditionID]; kk=dd.kin.get_loc(t.kinase)
        R=DOWN[:,c] if t.sign<0 else UP[:,c]
        def frac(ki,excl):
            s=sel[c,ki]; v=(s>=0)&(~ov[c,ki] if excl else True)
            return R[s[v]].mean() if v.any() else np.nan
        same=np.flatnonzero((rc==rc[kk])&dd.ev.neg_mask[c])
        for excl in (False,True):
            ft=frac(kk,excl); fn=np.nanmean([frac(j,excl) for j in same])
            rows.append(dict(lib=lib,noncur=excl,source=t.source,kinase=t.kinase,trc=rc[kk],cond=t.conditionID,f_tgt=ft,f_non=fn))
R=pd.DataFrame(rows); R.to_parquet(O+'enrich.parquet')
g=R.groupby(['lib','noncur','source','kinase'])[['f_tgt','f_non']].mean().groupby(['lib','noncur','source']).mean()
g['ratio']=g.f_tgt/g.f_non; print(g.round(3).to_string())
