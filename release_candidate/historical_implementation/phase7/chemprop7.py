"""Chemprop 2.3.1 D-MPNN components, explicit fixed Adam/BCE loop (no Noam)."""
from common7 import *
import torch,chemprop,joblib,random
from torch import nn
from rdkit import Chem
from chemprop.nn import BondMessagePassing,MeanAggregation,BinaryClassificationFFN
from chemprop.featurizers import SimpleMoleculeMolGraphFeaturizer
from chemprop.data.collate import BatchMolGraph
torch.set_num_threads(4)
torch.set_num_interop_threads(1)
GRAPH_CACHE=None

class GraphHead(nn.Module):
 def __init__(self,n_desc,dropout,d_v,d_e):
  super().__init__();self.mp=BondMessagePassing(d_v=d_v,d_e=d_e,d_h=128,depth=3,dropout=dropout);self.agg=MeanAggregation();self.head=BinaryClassificationFFN(input_dim=128+n_desc,hidden_dim=128,n_layers=2,dropout=dropout)
 def forward(self,bmg,xd):
  # Auxiliary scores are molecular descriptors, concatenated only AFTER aggregation.
  h=self.agg(self.mp(bmg),bmg.batch);return self.head.train_step(torch.cat((h,xd),dim=1)).squeeze(1)

def graphs(smiles):
 global GRAPH_CACHE
 if GRAPH_CACHE is None:
  p=A/'graph_cache.joblib'
  if p.exists():
   a=joblib.load(p);assert np.array_equal(a['smiles'],smiles);GRAPH_CACHE=a['graphs']
  else:
   f=SimpleMoleculeMolGraphFeaturizer();GRAPH_CACHE=[f(Chem.MolFromSmiles(s)) for s in smiles];assert len(GRAPH_CACHE)==750
   for g in GRAPH_CACHE:assert len(g.V)>0 and np.isfinite(g.V).all() and np.isfinite(g.E).all()
   joblib.dump(dict(smiles=smiles,graphs=GRAPH_CACHE,chemprop_version=chemprop.__version__),p,compress=3);event('A_technical','deterministic_graph_cache',identities=len(smiles),model_queries=0)
 return GRAPH_CACHE

def fit_graph(smiles,x,y,train,query,dropout):
 torch.manual_seed(4109);random.seed(4109);torch.use_deterministic_algorithms(True)
 gs=graphs(smiles);mean=x[train].mean(0);scale=x[train].std(0);scale[scale==0]=1;xd=torch.tensor((x-mean)/scale,dtype=torch.float32);ys=torch.tensor(y,dtype=torch.float32)
 net=GraphHead(x.shape[1],dropout,gs[0].V.shape[1],gs[0].E.shape[1]);dims=[net.head.ffn[0][-1].in_features]+[layer[-1].out_features for layer in net.head.ffn];assert dims==[128+x.shape[1],128,128,1]
 opt=torch.optim.Adam(net.parameters(),lr=.001);loss=nn.BCEWithLogitsLoss();rng=torch.Generator(device='cpu').manual_seed(4109);trace=[];steps=0
 for epoch in range(30):
  net.train();order=train[torch.randperm(len(train),generator=rng).numpy()];epoch_loss=0.
  for start in range(0,len(order),64):
   idx=order[start:start+64];batch=BatchMolGraph([gs[i] for i in idx]);opt.zero_grad(set_to_none=True);value=loss(net(batch,xd[idx]),ys[idx]);assert torch.isfinite(value);value.backward();opt.step();epoch_loss+=float(value.detach())*len(idx);steps+=1
  trace.append(epoch_loss/len(train))
 net.eval();pv=[]
 with torch.no_grad():
  for start in range(0,len(query),64):
   idx=query[start:start+64];pv.append(net(BatchMolGraph([gs[i] for i in idx]),xd[idx]).sigmoid().numpy())
 p=np.concatenate(pv);obj=dict(state_dict={k:v.detach().cpu() for k,v in net.state_dict().items()},scaler_mean=mean,scaler_scale=scale,architecture=dict(n_desc=x.shape[1],dropout=dropout,d_v=gs[0].V.shape[1],d_e=gs[0].E.shape[1]),classes_=np.array([0,1]),version=chemprop.__version__,train_idx=train,feature_order='original8 plus optional original3',epochs=30)
 info=dict(epochs=30,optimizer='Adam',optimizer_lr=.001,scheduler=None,batch_size=64,optimization_steps=steps,ffn_dimensions=dims,aggregation='mean',message_hidden=128,message_depth=3,loss='BCEWithLogitsLoss',class_weights=None,batch_normalization=False,training_loss_by_epoch=trace,dropout=dropout,torch_threads=torch.get_num_threads(),torch_interop_threads=torch.get_num_interop_threads(),device='cpu',pretrained=False,query_targets_used=False,scaler_fit_rows=len(train),data_loader_workers=0)
 return obj,p,info
