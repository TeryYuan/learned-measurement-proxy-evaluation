"""Weighted block Ridge with training-only dose interactions and dual Cholesky."""
from b1_common import *
from scipy.linalg import cho_factor,cho_solve
GRID=[.001,.1,10.,'DROP']
ROUTES=['D_cond','D_struct','D_meas_only','D_pred','D_meas','D_perm']
def candidates(route):
    if route=='D_cond':return [('DROP','DROP')]
    if route=='D_struct':return [(a,'DROP') for a in GRID]
    if route=='D_meas_only':return [('DROP',a) for a in GRID]
    return list(itertools.product(GRID,GRID))
def keypair(a,b):return f'X={a}|Z={b}'
def tie_key(pair):
    a,b=pair;active=sum(v!='DROP' for v in pair)
    return (active,-(np.inf if b=='DROP' else b),-(np.inf if a=='DROP' else a))
def select(scores):
    best=min(scores.values());ties=[p for p,v in scores.items() if v<=best+1e-10];return min(ties,key=tie_key)
class WeightedScale:
    def fit(self,a,q,center_only=False):
        a=np.asarray(a,dtype=float);self.mean=q@a;self.constant=np.ptp(a,axis=0)==0;self.mean[self.constant]=a[0,self.constant]
        self.var=q@(a-self.mean)**2;self.scale=np.sqrt(self.var);self.scale[self.constant]=1
        self.center_only=center_only
        if center_only:self.scale=np.ones(a.shape[1])
        return self
    def transform(self,a):
        v=(np.asarray(a)-self.mean)/self.scale;v[:,self.constant]=0;return v
class BlockTransform:
    def fit(self,c,x,z,route,q):
        self.cs=WeightedScale().fit(c,q);self.bits=WeightedScale().fit(x[:,:2048],q,True);self.desc=WeightedScale().fit(x[:,2048:],q);self.zs=WeightedScale().fit(z,q)
        self.categories=sorted(set(route));self.energyX=1.;self.energyZ=1.
        b=self.transform(c,x,z,route)
        self.energyX=float(np.sqrt(q@np.sum(b['X']**2,axis=1)));self.energyZ=float(np.sqrt(q@np.sum(b['Z']**2,axis=1)))
        return self
    def transform(self,c,x,z,route):
        cc=self.cs.transform(c);ud=cc[:,1,None];X=np.column_stack([self.bits.transform(x[:,:2048]),self.desc.transform(x[:,2048:])]);Z=self.zs.transform(z)
        xp=np.column_stack([X,ud*X]);zp=np.column_stack([Z,ud*Z]);onehot=(np.asarray(route)[:,None]==np.array(self.categories)[None,:]).astype(float)
        return dict(C=np.column_stack([cc,onehot]),X=xp/self.energyX if self.energyX>0 else np.zeros_like(xp),Z=zp/self.energyZ if self.energyZ>0 else np.zeros_like(zp))
    def support(self,blocks,q):
        rows=[]
        for name,b in blocks.items():
            bc=(b-q@b)*np.sqrt(q)[:,None];gram=bc@bc.T;ev=np.linalg.eigvalsh(gram);tol=max(bc.shape)*np.finfo(float).eps*max(float(ev.max()),0)
            rows.append(dict(block=name,n_columns=b.shape[1],constant_base_columns=int(self.cs.constant.sum()) if name=='C' else int(self.bits.constant.sum()+self.desc.constant.sum()) if name=='X' else int(self.zs.constant.sum()),unscaled_energy=np.nan if name=='C' else self.energyX if name=='X' else self.energyZ,scaled_weighted_energy=float(np.sqrt(q@np.sum(b*b,axis=1))),centered_gram_numeric_rank=int(np.sum(ev>tol)),rank_tolerance=tol,zero_energy=(name=='X' and self.energyX==0) or (name=='Z' and self.energyZ==0)))
        return rows
class DualContext:
    def __init__(self,blocks,y,q):
        self.q=np.asarray(q);self.sqrtq=np.sqrt(q);self.muy=float(q@y);self.y=np.asarray(y);self.t=self.sqrtq*(y-self.muy)
        self.means={k:q@v for k,v in blocks.items()};self.centered={k:v-self.means[k] for k,v in blocks.items()};self.grams={k:(v*self.sqrtq[:,None])@(v*self.sqrtq[:,None]).T for k,v in self.centered.items()}
    def query_cache(self,blocks):return {k:(v-self.means[k])@self.centered[k].T for k,v in blocks.items()}
    def solve(self,lx,lz,querycache,save_coefficients=False):
        penalties={'C':.01}
        if lx!='DROP':penalties['X']=float(lx)
        if lz!='DROP':penalties['Z']=float(lz)
        G=np.eye(len(self.y))
        for k,lam in penalties.items():G+=self.grams[k]/lam
        u=cho_solve(cho_factor(G,lower=True,check_finite=True),self.t,check_finite=True);v=self.sqrtq*u
        pred=np.full(next(iter(querycache.values())).shape[0],self.muy)
        for k,lam in penalties.items():pred+=(querycache[k]@v)/lam
        coeff={k:self.centered[k].T@v/lam for k,lam in penalties.items()} if save_coefficients else None
        if not np.isfinite(pred).all():raise FloatingPointError('non-finite candidate prediction')
        return pred,coeff
class BandedModel:
    def __init__(self,transform,dual,coeff,train_ids,allowed,lx,lz,context,zmean,columns,genesets):
        self.transformer=transform;self.block_means=dual.means;self.intercept=dual.muy;self.coefficients=coeff;self.train_ids=list(train_ids);self.allowed=set(allowed);self.lambdaX=lx;self.lambdaZ=lz;self.context=context;self.measured_training_zmean=zmean;self.structure_columns=columns;self.genesets=genesets
    def predict(self,ids,c,x,z,route):
        guard(ids,self.allowed)
        if set(ids)&set(self.train_ids):raise PermissionError('Outer query overlaps fitted training parents')
        b=self.transformer.transform(c,x,z,route);p=np.full(len(ids),self.intercept)
        for k,v in self.coefficients.items():p+=(b[k]-self.block_means[k])@v
        return p
def permute_vectors(ids,z,outer,inner):
    ids=sorted(ids);seed=2718+100*outer+inner;source=np.array(ids);np.random.Generator(np.random.PCG64(seed)).shuffle(source)
    d=pd.DataFrame(z.loc[source].to_numpy(),index=ids,columns=z.columns)
    return d,[dict(context=f'outer{outer}_'+('refit' if inner==3 else f'inner{inner}'),seed=seed,parent_id=p,assigned_z_parent=s) for p,s in zip(ids,source)]
