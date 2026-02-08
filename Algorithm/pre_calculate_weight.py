import igraph as ig
import sys
G = ig.Graph.Read_GraphML(sys.argv[1])
G_name = G.vs['name']
pre_weight_dict = {}
for edge in range(len(G.es)):
    now_weight = G.es[edge]['weight']
    if now_weight>1:
        pre_weight_dict[tuple(sorted([G_name[G.es[edge].source],G_name[G.es[edge].target]]))] = now_weight
import pickle
with open('pre_weight_dict.pkl','wb') as f:
    pickle.dump(pre_weight_dict,f)
