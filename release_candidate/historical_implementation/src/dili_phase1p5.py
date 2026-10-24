"""Finite forest-partition residual fusion. No N learning entry point."""
import pathlib,json,datetime,time,warnings
import numpy as np
from scipy.special import expit,logit
from scipy.optimize import minimize
from scipy.linalg import eigh
from dili_phase1_log import append_jsonl
ROOT=pathlib.Path(__file__).resolve().parents[1]
LOG=ROOT/'logs/phase1p5'
TASKS=['BSEPi','MTX_MP','ARE']
ROUTES=['R_score','R_leaf','R_gate','R_gate_perm','R_gate_chem','R_gate_struct']
LAMBDAS=[.001,.01,.1,1.,10.,'NO_CORRECTION']
EPS=1e-6

def event(kind,**kw):
    append_jsonl(LOG/'events.jsonl',{'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'event':kind,**kw})

def guard(ids,allowed,forbidden=(),role='O'):
    if role!='O':raise PermissionError('Phase1.5 learned access is O only; N and dedicated C are closed')
    a=list(map(str,ids));s=set(a)
    if len(s)!=len(a) or not s<=set(allowed) or s&set(forbidden):raise PermissionError('Unauthorized or held-out identity in learned operation')

def leaf_kernel(a,b):
    a=np.asarray(a);b=np.asarray(b)
    if a.ndim!=2 or b.ndim!=2 or a.shape[1]!=b.shape[1] or not a.shape[1]:raise ValueError('Same forest and tree columns required')
    k=np.zeros((len(a),len(b)),np.float64)
    for t in range(a.shape[1]):k+=a[:,t,None]==b[None,:,t]
    return k/a.shape[1]

def tanimoto(a,b):
    a=np.asarray(a,dtype=np.float64);b=np.asarray(b,dtype=np.float64)
    inter=a@b.T;den=a.sum(1)[:,None]+b.sum(1)[None,:]-inter
    return np.divide(inter,den,out=np.ones_like(inter),where=den!=0)

def build_kernels(parts,q,z,struct):
    gate={branch:sum(q[:,k,None]*q[None,:,k]*parts[(task,branch)] for k,task in enumerate(TASKS))/3 for branch in ['real','perm','chem']}
    out={'R_score':z@z.T/3,'R_leaf':sum(parts[(t,'real')] for t in TASKS)/3,'R_gate':gate['real'],
         'R_gate_perm':gate['perm'],'R_gate_chem':gate['chem'],'R_gate_struct':(q@q.T/3)*struct}
    for k,t in enumerate(TASKS):out['R_gate_without_'+t]=gate['real']-q[:,k,None]*q[None,:,k]*parts[(t,'real')]/3
    return out

def spectral(k,rank_cap=64):
    k=np.asarray(k,np.float64)
    if not np.isfinite(k).all() or not np.allclose(k,k.T,atol=1e-12,rtol=1e-12):raise ValueError('Nonfinite/asymmetric training kernel')
    mu,u=eigh((k+k.T)/2,check_finite=True)
    scale=max(1.,float(mu[-1]));tol=1e-10*scale
    if mu[0]<-1e-8*scale:raise ValueError('Materially negative kernel eigenvalue')
    ix=np.flatnonzero(mu>tol)[-rank_cap:][::-1];vals=mu[ix];vecs=u[:,ix]
    basis={'values':vals,'vectors':vecs}
    pos=mu[mu>tol]
    info={'rank':len(ix),'positive_rank':len(pos),'eigen_min':float(mu[0]),'eigen_max':float(mu[-1]),
          'retained_eigensum_fraction':float(vals.sum()/pos.sum()) if len(pos) else 0.,'roundoff_negative_count':int((mu<0).sum())}
    return vecs*np.sqrt(vals),basis,info

def project(kqs,basis):
    return np.asarray(kqs)@basis['vectors']/np.sqrt(basis['values'])

def objective(w,phi,offset,y,lam):
    eta=offset+phi@w;p=expit(eta)
    return float(np.mean(np.logaddexp(0,eta)-y*eta)+lam*np.dot(w,w)/2),phi.T@(p-y)/len(y)+lam*w

def fit_residual(phi,raw_cf,y,lam):
    off=logit(np.clip(raw_cf,EPS,1-EPS));zero=np.zeros(phi.shape[1]);l0=float(np.mean(np.logaddexp(0,off)-y*off))
    if lam=='NO_CORRECTION' or phi.shape[1]==0:
        return zero,{'success':True,'iterations':0,'retried':False,'objective':l0,'objective_zero':l0,'weight_norm':0.,'norm_bound':None,'message':'explicit zero correction' if lam=='NO_CORRECTION' else 'zero rank'}
    result=None;retried=False
    for maxiter in [2000,5000]:
        result=minimize(objective,zero,args=(phi,off,np.asarray(y),float(lam)),method='L-BFGS-B',jac=True,
                        options={'maxiter':maxiter,'gtol':1e-7,'ftol':1e-10})
        if result.success:break
        retried=True
    norm=float(np.linalg.norm(result.x));bound=float(np.sqrt(2*l0/float(lam)))
    valid=bool(result.success and np.isfinite(result.fun) and result.fun<=l0+1e-9 and norm<=bound+1e-7)
    return result.x,{'success':valid,'iterations':int(result.nit),'retried':retried,'objective':float(result.fun),
       'objective_zero':l0,'weight_norm':norm,'norm_bound':bound,'message':str(result.message)}

def residual_predict(raw,phi,w,lam):
    base=np.clip(raw,EPS,1-EPS)
    if lam=='NO_CORRECTION' or not np.any(w):return base.copy(),np.zeros(len(base))
    h=phi@w
    return expit(logit(base)+h),h

def choose_lambda(rows):
    ok=[r for r in rows if r['valid']]
    best=max(r['mean_AP'] for r in ok)
    tied=[r for r in ok if best-r['mean_AP']<=1e-12]
    return sorted(tied,key=lambda r: (0 if r['lambda']=='NO_CORRECTION' else 1,-float(r['lambda']) if r['lambda']!='NO_CORRECTION' else 0))[0]
