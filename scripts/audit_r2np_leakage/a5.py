import pandas as pd, numpy as np
O='/workspace/cache/audit_r2np/'; import os; os.makedirs(O, exist_ok=True)
P=pd.read_parquet(O+'pairs_nulls.parquet'); P=P[P.in_M1sel]
c=pd.read_parquet('/workspace/cache/conditions.parquet').set_index('conditionID')
sfk=['ABL1','SRC','FYN','YES1','LYN','LCK','EPHA2','KIT','HCK','FGR','BLK','ABL2','BTK','EPHB4','DDR1','PDGFRB']
q=P[P.kinase.isin(sfk)]
q=q.assign(cond=q.conditionID.map(c.Condition),pmid=q.conditionID.map(c.PMID),cell=q.conditionID.map(c.Sample))
pd.set_option('display.width',250); pd.set_option('display.max_rows',200)
t=q.pivot_table(index=['source','pmid','cond','cell','kinase'],columns='cfg',values='auroc')
print(t[['r2np50m','cant50m','r3_50m','psp','null_resclass_random@50','null_popularity_devcond@50']].round(2).query("source!='D'").to_string())
# D summary by drug (pooled over kinases)
d=q[q.source=='D']; d=d.assign(drug=d.cond.str.split('_').str[0])
print(d.pivot_table(index=['drug','cell'],columns='cfg',values='auroc')[['r2np50m','cant50m','r3_50m','psp','null_resclass_random@50']].round(2))
print(d.pivot_table(index=['drug','cell'],columns='cfg',values='auroc_within_Y')[['r2np50m','cant50m','r3_50m','psp','null_resclass_random@50']].round(2))
