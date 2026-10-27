import sys, numpy as np, pandas as pd
sys.path.insert(0,'/workspace/src')
from kinact.explore.lib_variants import DevData, pred_dense
from kinact.library import SiteMapper
dd=DevData()
uni=pd.read_parquet('/workspace/external/ksm_index/universe.parquet',columns=['site_id','src_psp','src_cptac','src_ochoa','src_compendium','n_sources'])
uni.index=SiteMapper().from_keys(uni.site_id); uni=uni[~uni.index.duplicated()].reindex(dd.sites)
componly=(uni.src_compendium==True)&~(uni.src_psp==True)&~(uni.src_cptac==True)&~(uni.src_ochoa==True)
co=componly.fillna(False).to_numpy(bool); ns=uni.n_sources.fillna(0).to_numpy()
res=np.array([s.split('_')[1][0] for s in dd.sites])
print('DEV sites',len(co),'compendium-only frac',co.mean())
for lib in ['r2np','r3','cantley']:
    D=pred_dense(dd,lib); Dn=np.where(np.isnan(D),-np.inf,D)
    top=np.argsort(-Dn,axis=1)[:,:200]
    ok=np.isfinite(np.take_along_axis(Dn,top,1))
    f=[];b=[];n=[];nb=[]
    for i in range(D.shape[0]):
        t=top[i][ok[i]]
        if len(t)<50: continue
        cand=np.isfinite(D[i])
        f.append(co[t].mean()); b.append(co[cand].mean()); n.append(ns[t].mean()); nb.append(ns[cand].mean())
    print(lib,'top200 compendium-only frac %.3f vs candidate base %.3f | mean n_sources %.2f vs %.2f'%(np.mean(f),np.mean(b),np.mean(n),np.mean(nb)))
