import sys
import numpy as np
import igraph as ig
import math
import bisect
import sys
import numpy as np
import igraph as ig
from itertools import combinations
from collections import Counter
from collections import defaultdict, deque

import argparse
parser = argparse.ArgumentParser(description="Process some parameters.")
parser.add_argument('-m','--MboI', required=True, help="MboI restriction site file")
parser.add_argument('-g','--Graph', required=True, help="Graph file")
parser.add_argument('-t','--top_frag_dict', required=True,help="top_frag_dict.npy")
parser.add_argument('-d','--droplets',required=True,help='fragments in each droplet')
parser.add_argument('-x', '--max_size',type=int, default=30, help="max size threshold for connected components")
parser.add_argument('-y', '--fragment_each',type=int, default=13, help="stop random walk if the mean size of d-LHCC exceed this number of fragments for each droplet") 
parser.add_argument('-z', '--expand',type=int, default=5, help="expansion cutoff for neighborhood")
parser.add_argument('-i', '--inter_ratio',type=float, default=0.08, help="inter-chromosomal connections allowed in each droplet")
args = parser.parse_args()
chrindex={}
for i in range(1,23):
    chrindex[str(i)]=i
chrindex['X']=23
chrindex['Y']=24
import pickle
#with open(args.model,'rb') as f:
#    model_poly = pickle.load(f)
with open('pre_weight_dict.pkl','rb') as f:
    pre_weight_dict = pickle.load(f)

def make_MboI(filepath):
    """Load MboI restriction site data"""
    MboI = {}
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            chrom = parts[0]
            if chrom not in chrindex:
                continue
            try:
                MboI[chrom] = np.array([1] + [int(x) for x in parts[1:]])
            except ValueError:
                continue
    return MboI

MboI = make_MboI(args.MboI)
G = ig.Graph.Read_GraphML(args.Graph)
top_frag_dict = np.load(args.top_frag_dict,allow_pickle=True).item()
G_name = G.vs['name']
with open('top_frags.txt','r') as f:
    top_frags = f.readlines()[0]
    top_frags = set(top_frags.split(' ')[:-1])
def get_k(distance):
    if distance<1000:
        return 1
    elif distance>1e5:
        return 0.2
    else:
        return -(np.log10(distance)-3)*0.4+1
def get_edge_distance(edge):
    return calculate_distance_frag(edge[0],edge[1])
def calculate_distance_frag(frag1,frag2):
    ch_loc1 = frag1.split('#')
    ch_loc2 = frag2.split('#')
    if ch_loc1[0]!=ch_loc2[0]:
        return 1000000000
    return abs(int(ch_loc1[1])-int(ch_loc2[1]))
def find_components_union_find(n, edges):
    """
    使用并查集找连通分量
    
    参数:
        n: 节点数量 (节点编号从0到n-1)
        edges: 边列表 [(u, v), ...]
    
    返回:
        components: 列表的列表，每个子列表是一个连通分量
    """
    # 初始化并查集
    parent = list(range(n))
    size = [1] * n
    
    def find(x):
        # 路径压缩
        while parent[x] != x:
            parent[x] = parent[parent[x]]  # 路径压缩优化
            x = parent[x]
        return x
    
    def union(x, y):
        root_x = find(x)
        root_y = find(y)
        if root_x == root_y:
            return False
        
        # 按大小合并（小树合并到大树）
        if size[root_x] < size[root_y]:
            root_x, root_y = root_y, root_x
        
        parent[root_y] = root_x
        size[root_x] += size[root_y]
        return True
    
    # 处理所有边
    for u, v in edges:
        union(u, v)
    
    # 收集连通分量
    from collections import defaultdict
    comp_dict = defaultdict(list)
    
    for i in range(n):
        root = find(i)
        comp_dict[root].append(i)
    
    return list(comp_dict.values())
def get_sum_weight(node,expand=0):
    final_weight = defaultdict(int)
    expand_node = expand_component([node],cut=expand)
    all_neighbor_weight_dict = {node:get_neighbor_weights(G,node,G_name) for node in expand_node}
    #print(all_neighbor_weight_dict)
    for node1 in all_neighbor_weight_dict:
        _,weight_dict = all_neighbor_weight_dict[node1]  #weight_dict: neighbor, weight
        for node2 in weight_dict:
            final_weight[node2]+=weight_dict[node2]
    return final_weight  #neighbor, weight
def get_neighbor_weights(graph,node,G_name,remove_set = set()):
    try:
        neighbors = graph.neighbors(node,mode=ig.ALL)
    except ValueError:
        return set(),{}
    neighbors = [G_name[i] for i in neighbors]
    neighbors = [i for i in neighbors if i not in remove_set]
    all_edges = [tuple(sorted([node,neighbor])) for neighbor in neighbors]
    neighbor_weights = [pre_weight_dict.get(edge,1)*get_k(get_edge_distance(edge)) for edge in all_edges]
    #neighbor_weights = [pre_weight_dict.get(edge,1) for edge in all_edges]
    #eids = graph.get_eids(all_edges)
    #print(eids)
    #neighbor_weights = [graph.es[eid]['weight'] for eid in eids]
    #neighbor_weights = [pre_weight_dict[edge] for edge in all_edges]
    #neighbor_weights = [graph.es[graph.get_eid(node,neighbor)]['weight'] for neighbor in neighbors]
    weight_dict = dict([(neighbors[i],neighbor_weights[i]) for i in range(len(neighbors))])
    return set(neighbors),weight_dict
def make_cluster_frag_dict(random_droplet):
    ch_dict = split_by_chr(random_droplet)
    cluster_frag_dict = {}
    for ch in ch_dict:
        clusters = frags2clusters(ch_dict[ch])
        for cluster in clusters:
            frag0 = f'{ch}#{cluster[0]}'
            cluster_frag_dict[frag0] = [f'{ch}#{i}' for i in cluster]
    return cluster_frag_dict
def make_neighbor_weight_reverse(random_neighbors):
    neighbor_weight_reverse = defaultdict(dict) #key, neighbor, values: {node, weight}
    for node in random_neighbors:
        for neighbor in random_neighbors[node]:
            neighbor_weight_reverse[neighbor][node] = random_neighbors[node][neighbor]
    return neighbor_weight_reverse
def get_comp_node_neighbor_dict(components0,components0_expand,G,G_name,expand):
    comp_node_neighbor_dict = {} #keys: components ,values: dict{keys:node,values:[set(neighbors),weight_dict]}
    #print(len(components0))
    for i in range(len(components0)):
        #print(i,len(components0_expand[i]))
        comp_node_neighbor_dict[i] = {node:get_neighbor_weights(G,node,G_name,remove_set = set(components0_expand[i])) for node in components0_expand[i]}
    return comp_node_neighbor_dict
def calculate_prob_of_each_neighbor(comp_node_neighbor_dict):
    comp_node_prob = {}#comp:{node,prob to pick this node}
    comp_neighbor_prob = defaultdict(dict) #comp:{node:{neighbor:prob to pick this neighbor}}
    cummulative_neighbor_prob = defaultdict(dict)
    for comp in comp_node_neighbor_dict:
        now_comp_prob = {}
        for node in comp_node_neighbor_dict[comp]:
            weight_dict = comp_node_neighbor_dict[comp][node][1]
            total_weight = sum(weight_dict.values())
            now_comp_prob[node] = total_weight
            node_prob_dict = {i:weight_dict[i]/total_weight for i in weight_dict if i not in comp_node_neighbor_dict[comp].keys()}
            comp_neighbor_prob[comp][node] = node_prob_dict
        comp_total_weight = sum(now_comp_prob.values())
        try:
            comp_node_prob[comp] = {node:now_comp_prob[node]/comp_total_weight for node in now_comp_prob}
        except ZeroDivisionError:
            comp_node_prob[comp] = {}
    for comp in comp_node_prob:
        comp_prob = comp_node_prob[comp]
        if comp_prob=={}:
            cummulative_neighbor_prob[comp] = {}
        for node in comp_neighbor_prob[comp]:
            for neighbor in comp_neighbor_prob[comp][node]:
                now_val = comp_node_prob[comp][node]*comp_neighbor_prob[comp][node][neighbor]
                try:
                    cummulative_neighbor_prob[comp][neighbor]+=now_val
                except KeyError:
                    cummulative_neighbor_prob[comp][neighbor]=now_val
    #print(np.max([sum(cummulative_neighbor_prob[comp].values()) for comp in cummulative_neighbor_prob]),np.min([sum(cummulative_neighbor_prob[comp].values()) for comp in cummulative_neighbor_prob]))
    return cummulative_neighbor_prob
def connected_comp_summary(frag_comp_index):
    to_connect = set()
    same_dict = {}  #统计每个comp和哪些连接了（component邻接表）
    for frag in frag_comp_index:
        if len(frag_comp_index[frag])<2:
            continue
        for comp1,comp2 in combinations(frag_comp_index[frag],2):
            to_connect.add(tuple(sorted([comp1,comp2])))
    return to_connect
def make_same_dict(final_components):
    same_dict = defaultdict(set)  #统计每个comp和哪些连接了（component邻接表）
    for comps in final_components:
        for comp1,comp2 in combinations(comps,2):
            same_dict[comp1].add(comp2)
            same_dict[comp2].add(comp1)
    for i in same_dict:
        same_dict[i].add(i)
    return same_dict    
    
#def combine_sequential(common_neighbor_weight):
def expand_component(comp,cut=5):
    comp_set = set(comp)
    for i in comp:
        ch_loc = i.split('#')
        ch = ch_loc[0]
        loc = int(ch_loc[1])
        for j in range(max(0,loc-cut),min(loc+cut+1,len(MboI[ch]))):
            now_frag = ch+'#'+str(j)
            if now_frag not in top_frags:
                comp_set.add(now_frag)
    return list(comp_set)
def get_component_size(components0):
    frag_clusternum_dict = {}
    frag_fragnum_dict = {}
    for comp in components0:
        ch_dict = split_by_chr(comp)
        cluster_num = sum([len(frags2clusters(ch_dict[ch])) for ch in ch_dict])
        for frag in comp:
            frag_clusternum_dict[frag] = cluster_num
            frag_fragnum_dict[frag] = len(comp)
    #return frag_clusternum_dict
    return frag_fragnum_dict
def get_mode(i):
    mode='inter'
    if i[0].split('#')[0]==i[1].split('#')[0]:
        mode='intra'
    return mode
def get_prob_from_weight(i,true_common_neighbor_weight,prob_dict,expand):
    mode = get_mode(i)
    weight = true_common_neighbor_weight[i]
    if weight not in prob_dict[expand][mode]:
        prob=1
    else:
        prob = prob_dict[expand][mode][weight]
    return prob
def get_frag_comp_index(components0):
    frag_comp_index = defaultdict(set)
    for idx,comp in enumerate(components0):
        for frag in comp:
            frag_comp_index[frag].add(idx)
    return frag_comp_index

import bisect
def sample_neighbor(neighbors,cum_probs,r):
    idx = bisect.bisect(cum_probs, r)
    return neighbors[idx]
def simple_connect_components(components, edges_to_connect, max_size, max_connections=None):
    """
    按顺序连接组件，检查大小限制
    
    参数:
        components: 初始组件列表 [[frag1, frag2, ...], ...]
        edges_to_connect: 待连接的边列表 [(comp1_idx, comp2_idx), ...]
        max_size: 每个连通分量的最大fragment数量
        max_connections: 最大连接数（None表示连接所有）
    
    返回:
        final_components: 最终的连通分量列表
    """
    n = len(components)
    parent = list(range(n))  # 并查集父节点
    
    # 预计算每个组件的大小
    comp_sizes = [len(comp) for comp in components]
    
    def find(x):
        """查找根节点"""
        while parent[x] != x:
            parent[x] = parent[parent[x]]  # 路径压缩
            x = parent[x]
        return x
    
    # 按顺序处理每条边
    connections_made = 0
    for comp1, comp2 in edges_to_connect:
        # 检查是否达到最大连接数
        if max_connections is not None and connections_made >= max_connections:
            break
        
        # 查找根节点
        root1 = find(comp1)
        root2 = find(comp2)
        
        # 如果已经在同一分量中，跳过
        if root1 == root2:
            continue
        
        # 计算两个分量的大小
        size1 = sum(comp_sizes[i] for i in range(n) if find(i) == root1)
        size2 = sum(comp_sizes[i] for i in range(n) if find(i) == root2)
        
        # 检查合并后是否超过限制
        if size1 + size2 > 1.5*max_size:
            continue
        if size1 + size2 > max_size and min(size1,size2)>3:
            continue
        
        # 合并
        parent[root2] = root1
        connections_made += 1
    
    # 收集最终结果
    from collections import defaultdict
    groups = defaultdict(list)
    for i in range(n):
        root = find(i)
        groups[root].append(i)
    
    # 转换为fragment列表
    final_components = []
    for comp_indices in groups.values():
        combined = []
        for idx in comp_indices:
            combined.extend(components[idx])
        final_components.append(sorted(combined))
    
    return final_components
def get_chs(component):
    return set([i.split('#')[0] for i in component])
def find_components_and_expand(r,G,expand):
    
    bc_frags_bc = []
    r = r.strip().split(' ')
    bc = r[0]
    r = r[1:]
    for i in r:
        bc_frags_bc.append(i)
    bc_frags_set = set(bc_frags_bc)
    #print(len([i for i in true_common_neighbor_weight if true_common_neighbor_weight[i]>0]),len(true_common_neighbor_weight))
    #print(bc,top_frag_dict[bc])
    #print('start component')
    G_sub = G.induced_subgraph(bc_frags_bc)
    try:
        G_sub.delete_vertices(list(top_frag_dict[bc].keys()))
        G_sub.add_vertices(list(top_frag_dict[bc].keys()))
    except:
        pass
    for i in top_frag_dict[bc].keys():
        try:
            G_sub.add_edges(top_frag_dict[bc][i])
        except:
            pass
    subname = G_sub.vs['name']
    components0 = G_sub.components()               
    components_new = list(components0)
    components0 = [[subname[i] for i in comp] for comp in components0]  ##所有components
    components0_expand = [expand_component(i,expand) for i in components0]
    return components0,components0_expand,components_new,bc,subname,bc_frags_set

def random_walk_connection(components0,components0_expand,components_new,bc,subname,bc_frags_set,fragment_each = 20):
    global random_nums
    intra_pair = set([tuple(sorted([comp1,comp2])) for comp1,comp2 in combinations(range(len(components0)),2) if get_chs(components0[comp1])&get_chs(components0[comp2])!=set()])
    frag_comp_index = get_frag_comp_index(components0)   #assign each node to component
    comp_node_neighbor_dict = get_comp_node_neighbor_dict(components0,components0_expand,G,G_name,args.expand)
    cummulative_neighbor_prob = calculate_prob_of_each_neighbor(comp_node_neighbor_dict)
    final_components_num = len(bc_frags_set)/fragment_each
    max_connections = int(len(components0) - final_components_num)
    count={'intra':0,'inter':0}
    true_connect = 0
    max_connections_mode = {'intra':max_connections,'inter':max_connections*args.inter_ratio}
    to_connect = set()
    same_dict = {i:{i} for i in range(len(components0))}
    cycle = 0
    max_cycle = 20000
    if max_connections<5:
        max_cycle = 100
    cummulative_neighbor_prob_neighbors = {comp:list(cummulative_neighbor_prob[comp].keys()) for comp in cummulative_neighbor_prob}
    cummulative_neighbor_prob_probs = {comp:list(cummulative_neighbor_prob[comp].values()) for comp in cummulative_neighbor_prob}
    comp_cumsums = {comp:np.cumsum(cummulative_neighbor_prob_probs[comp]) for comp in cummulative_neighbor_prob_probs}
    final_components = [[i] for i in range(len(components_new))]
    candidates = defaultdict(int)
    while True:
        cycle+=1
        for comp_large in final_components:
            if sum([len(components0[i]) for i in comp_large])>=args.max_size:
                continue
            try:
                r = random_nums.pop()
            except IndexError:
                random_nums = list(np.random.random(1000000))
            comp = comp_large[int(r*len(comp_large))]
            neighbors = cummulative_neighbor_prob_neighbors[comp]
            if len(neighbors)==0:
                continue
            cum_probs = comp_cumsums[comp]
            try:
                r = random_nums.pop()
            except IndexError:
                random_nums = list(np.random.random(1000000))
            rand_neighbor = sample_neighbor(neighbors,cum_probs,r)
            now_connected_comp = set()
            should_connect=False
            if len(frag_comp_index[rand_neighbor]-{comp})>=1: #这个rand neighbor已经出现过了，它出现的总和加上新的，不超过最大
                for i in frag_comp_index[rand_neighbor]:
                    now_pair = tuple(sorted([i,comp]))
                    if i==comp:
                        continue
                    candidates[now_pair]+=1 
            frag_comp_index[rand_neighbor].add(comp)   #如果这个连接制造出了新的大comp，就不连接。
        candidates_cutoff = [i for i in candidates if candidates[i]>3]
        candidates_filtered = []
        inter_count = 0
        for candidate in candidates_cutoff:
            if candidate in intra_pair:
                candidates_filtered.append(candidate)
            elif inter_count<=max_connections*args.inter_ratio:
                inter_count+=1
                candidates_filtered.append(candidate)
        """
        inter和intra分开
        """
        final_components = find_components_union_find(len(components0), candidates_filtered)
        if len(final_components)<=final_components_num:
            break
        if cycle>max_cycle:
            break
    """
    从candidates从上往下选，每次选的时候check size
    """
    candidates = dict(sorted(candidates.items(),key=lambda x:x[1],reverse=True))
    final_components =  simple_connect_components(components0, candidates.keys(), args.max_size, max_connections)
    return final_components,bc,subname
def write_dLHCC(bc,components_new,fout,subname):
    for comp in components_new:
        if comp==[]:
            continue
        fout.write(bc)
        for i in comp:
            fout.write(' ' + i)
        fout.write('\n')
def process_single_droplet(r,G,top_frag_dict,G_name,fout,max_size,fragment_each,inter_ratio,expand):
    components0,components0_expand,components_new,bc,subname,bc_frags_set = find_components_and_expand(r,G,expand)
    components_new,bc,subname = random_walk_connection(components0,components0_expand,components_new,bc,subname,bc_frags_set,fragment_each = fragment_each)
    write_dLHCC(bc,components_new,fout,subname)
def process_all_droplets(args):
    global random_nums
    random_nums = list(np.random.random(1000000))
    file_suffix = args.droplets.split('_')[-1]
    output_dir = '/'.join(args.droplets.split('/')[:-1])
    fout_name = 'Components_combine_local_expand_optimized_intra_RW_inter_close_'+str(args.max_size)+'_'+str(args.fragment_each)+'_'+str(args.expand)+'_'+str(round(args.inter_ratio,2))+'_'+file_suffix                                #expand 
    with open(args.droplets,'r') as fin,open(fout_name,'w') as fout:
        l=0
        for r in fin:
            l += 1
        #prepare，make induced subgraph, find components
            if l % 50 == 0:
                sys.stdout.write("\r%d" % l)
                sys.stdout.flush()
            process_single_droplet(r,G,top_frag_dict,G_name,fout,args.max_size,args.fragment_each,args.inter_ratio,args.expand)

process_all_droplets(args)
