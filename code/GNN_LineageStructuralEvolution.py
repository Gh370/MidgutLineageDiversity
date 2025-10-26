import os, torch, pandas as pd, numpy as np
from tqdm import tqdm
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader
from torch_geometric.nn import GCNConv, global_mean_pool
from sklearn.metrics.pairwise import cosine_similarity
import torch.nn.functional as F

graph_dir = r"~/.../code/GraphStructs"
all_lineages = [f.split("node_")[1].replace(".csv", "")
                for f in os.listdir(graph_dir) if f.startswith("node_")]
grouped = {d: [lid for lid in all_lineages if f"_{d}d" in lid]
           for d in ["33","43","53","63"]}

def load_graph(lid):
    node_df = pd.read_csv(f"{graph_dir}/node_{lid}.csv")
    edge_df = pd.read_csv(f"{graph_dir}/edge_{lid}.csv")

    gen_leaf = node_df[["generation", "is_leaf"]].values
    const_sz = np.full((gen_leaf.shape[0], 1), np.log(len(node_df)))  
    node_feat = np.concatenate([gen_leaf, const_sz], axis=1)         

    return Data(
        x=torch.tensor(node_feat, dtype=torch.float),
        edge_index=torch.tensor(edge_df.values.T, dtype=torch.long),
        y=torch.tensor([0]),
        lineage_id=lid
    )

# Graph Convolution AutoEncoder 
class GCAE(torch.nn.Module):
    def __init__(self, in_dim=3, hid=32, out_dim=16):
        super().__init__()
        self.enc1 = GCNConv(in_dim, hid)
        self.enc2 = GCNConv(hid, out_dim)
        self.dec  = torch.nn.Linear(out_dim, in_dim)   

    def encode(self, x, edge_index, batch):
        h = self.enc1(x, edge_index).relu()
        z = self.enc2(h, edge_index)
        return global_mean_pool(z, batch)          

    def forward(self, data):
        z = self.encode(data.x, data.edge_index, data.batch)
        x_hat = self.dec(z)[data.batch]            
        loss  = F.mse_loss(x_hat, data.x)
        return loss, z

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model  = GCAE().to(device)
opt    = torch.optim.Adam(model.parameters(), lr=1e-3)

train_graphs = [load_graph(lid) for day in ["43","53","63"] for lid in grouped[day]]
train_loader = DataLoader(train_graphs, batch_size=8, shuffle=True)

for epoch in range(1000):
    model.train(); total = 0
    for batch in train_loader:
        batch = batch.to(device); opt.zero_grad()
        loss, _ = model(batch)
        loss.backward(); opt.step()
        total += loss.item()
    if epoch % 20 == 0:
        print(f"Epoch {epoch:03d} | AE-loss {total:.4f}")

@torch.no_grad()
def encode_set(ids):
    model.eval(); vecs=[]
    for lid in tqdm(ids):
        g = load_graph(lid).to(device)
        g.batch = torch.zeros(g.x.size(0), dtype=torch.long, device=device)
        _, z = model(g)
        vecs.append(z.cpu().numpy())
    return pd.DataFrame(np.vstack(vecs), index=ids)

vecs = {d: encode_set(grouped[d]) for d in ["33","43","53","63"]}

def match(src_df, tgt_df, tgt_day, thr=0.9):
    sims = cosine_similarity(src_df.values, tgt_df.values)
    best = sims.argmax(axis=1)
    score = sims.max(axis=1)
    match = np.array(tgt_df.index)[best]
    match[score < thr] = "Unassigned"
    return pd.DataFrame({"Lineage": src_df.index,
                         f"Matched_{tgt_day}": match,
                         "Similarity": score})

pairs = [("33","43"),
         ("43","53"),
         ("53","63")]

for src, tgt in pairs:
    out = match(vecs[src], vecs[tgt], tgt)
    out.to_csv(f"{graph_dir}/match_{src}_to_{tgt}.csv", index=False)

print("✅ GC-AE embedding & matching finished")