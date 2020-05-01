"""
This module contains many functional that report on the score gained by different actions.
All these functions are used to decide which action will be used.
Most functions here do not actually perform the action, but just report on the scores gained
by performing the function. If they do perform the action, it is usually only in order to compute
the gained score, and the action is then reverted.
"""
import time
import heapq
import pandas as pd
import networkx as nx
from monet.module import Module
import numpy as np


def best_samples_to_add(module, glob_var):
    """"returns the best neighbors that max the immediate score and the score it would provide
    in case there is no such node will return None and -inf"""
    heap = []
    max_pats_to_add = glob_var.max_pats_per_action
    start_weight = module.get_weight()
    best_sample = (-float("inf"), None)

    for name, pat in glob_var.patients.items():
        if pat.is_in_module():
            continue

        tmp_weight = module.add_patient(pat)
        module.remove_patient(pat)

        if tmp_weight < start_weight:
            continue
        if tmp_weight > best_sample[0]:
            best_sample = tmp_weight, pat
        if len(heap) < max_pats_to_add:
            heapq.heappush(heap, (tmp_weight, name))
        elif tmp_weight > heap[0][0]:
            heapq.heappushpop(heap, (tmp_weight, name))

    for weight, sample_name in heap:
        module.add_patient(glob_var.patients[sample_name])
    total_addition = module.get_weight()
    for weight, sample_name in heap:
        module.remove_patient(glob_var.patients[sample_name])

    if best_sample[0] > total_addition:
        maximal_addition = best_sample[0]
        to_add = [best_sample[1]]
    else:
        maximal_addition = total_addition
        to_add = [glob_var.patients[sample_name] for weight, sample_name in heap]
    if maximal_addition > start_weight and len(to_add) > 0:
        return maximal_addition - start_weight, to_add
    return -float("inf"), None


def which_patient_to_remove(module, glob_var):
    """returns the node that is most beneficial to remove and
    the score it would provide
    in case there is no such node will return None and -inf"""
    if len(module.get_patients()) == glob_var.min_mod_size:
        return -float("inf"), None
    cur_max = module.get_weight()
    start_weight = module.get_weight()
    pat_to_del = None
    for name, pat in module.get_patients().items():
        tmp_weight = module.remove_patient(pat)
        if tmp_weight > cur_max:
            cur_max = tmp_weight
            pat_to_del = pat
        module.add_patient(pat)
    if pat_to_del:
        return cur_max - start_weight, pat_to_del
    return -float("inf"), None

    
def top_samples_to_switch(mod, glob_var):
    max_pats_to_switch = glob_var.max_pats_per_action
    init_num = (len(mod.get_patients_names_as_list()))
    res = []
    init_weight = sum([tmp_mod.get_weight() for mod_name, tmp_mod in glob_var.modules.items()])
    for mod2_name, mod2 in glob_var.modules.items():
        if mod == mod2:
            continue
        weight, samples_list = best_samples_to_switch(mod, mod2, glob_var)
        res += samples_list
        cur_num = (len(mod.get_patients_names_as_list()))
    if len(res) == 0:
        return -float("inf"), None
    res = sorted(res, key=lambda x: x[1], reverse=True)
    used = set()
    used_indices = set()
    for i, (pat, score, mod2) in enumerate(res):
        if pat in used:
            continue
        pat_mod = glob_var.patients[pat].module
        if len(pat_mod.get_patients()) == glob_var.min_mod_size:
            continue
        if len(used) == max_pats_to_switch:
            break
        switch_2_patients(pat, mod, mod2, glob_var) 
        used.add(pat)
        used_indices.add(i)
    curr_weight = sum([tmp_mod.get_weight() for mod_name, tmp_mod in glob_var.modules.items()])
    ret = []
    for i, (pat, score, mod2) in enumerate(res):
        if pat in used and i in used_indices:
            used.remove(pat)
            switch_2_patients(pat, mod, mod2, glob_var)
            ret += [(pat, score, mod2)]
    res = ret
    if len(res) == 0:
        return -float("inf"), None
    if (curr_weight-init_weight) > res[0][1]:
        return curr_weight-init_weight, res
    return res[0][1], [res[0]]


def best_samples_to_switch(mod1, mod2, glob_var):
    mod1_names = mod1.get_patients_names_as_list()
    mod2_names = mod2.get_patients_names_as_list()
    samples_dict = {node:{"in_mod": 0, "diff_mod": 0} for node in (mod1_names+mod2_names)}
    start_weight = mod1.get_weight() + mod2.get_weight()
    for omic_name, omic in mod1.get_omics().items():
        for node_name, degree in omic.graph.subgraph(mod1_names).degree(weight='weight'):
            samples_dict[node_name]["in_mod"] += degree
        for node_name in mod2_names:
            samples_dict[node_name]["diff_mod"] += \
                omic.graph.subgraph(mod1_names + [node_name]).degree(node_name, weight="weight")
    for omic_name, omic in mod2.get_omics().items():
        for node_name, degree in omic.graph.subgraph(mod2_names).degree(weight='weight'):
            samples_dict[node_name]["in_mod"] += degree
        for node_name in mod1_names:
            samples_dict[node_name]["diff_mod"] += \
                omic.graph.subgraph(mod2_names + [node_name]).degree(node_name, weight="weight")
    samples_list = [(node, dict["diff_mod"]-dict["in_mod"], mod2) for node, dict in samples_dict.items()]
    samples_list = sorted(samples_list, key=lambda x: x[1], reverse=True)
    num_switches = 0

    for i in range(min(len(samples_list), glob_var.max_pats_per_action)):
        if samples_list[i][1] <= 0:
            if i == 0:
                return -float("inf"), []
            break
        num_switches += 1
        switch_2_patients(samples_list[i][0], mod1, mod2, glob_var)
    samples_list = samples_list[:num_switches]
    switch_all_weight = mod1.get_weight() + mod2.get_weight() - start_weight
    for pat, diff, fake_mod in samples_list:
        switch_2_patients(pat, mod1, fake_mod, glob_var)

    if switch_all_weight < samples_list[0][1]:
        samples_list = [samples_list[0]]
        switch_all_weight = samples_list[0][1]
    if switch_all_weight > 0:
        return switch_all_weight, samples_list
    return -float("inf"), []


def switch_2_patients(pat, mod1, mod2, glob_var):
    if isinstance(pat, str):
        pat = glob_var.patients[pat]
    if pat.get_module() == mod1:
        mod1.remove_patient(pat)
        mod2.add_patient(pat)
    elif pat.get_module() == mod2:
        mod2.remove_patient(pat)
        mod1.add_patient(pat)
    else:
	# This should never occur!
        import pdb;pdb.set_trace()
    return


def which_omic_to_add_to_module(mod, glob_var):
    """get a module and see what happens when you add different omics"""
    cur_max = mod.get_weight(), None
    start_weight = mod.get_weight()
    for omic in glob_var.omics:
        if omic in mod.get_omics():
            continue
        mod.add_omic(omic, glob_var)
        tmp = mod.get_weight()
        if tmp > cur_max[0]:
            cur_max = tmp, omic
        mod.remove_omic(omic, glob_var)
    if cur_max[1]:
        return cur_max[0] - start_weight, cur_max[1]
    return -float("inf"), None


def which_omic_to_remove_from_module(mod, glob_var):
    cur_max = mod.get_weight(), None
    start_weight = mod.get_weight()
    if len(mod.get_omics()) <= 1:
        return -float("inf"), None
    for omic in mod.get_omics():
        mod.remove_omic(omic, glob_var)
        tmp = mod.get_weight()
        if tmp > cur_max[0]:
            cur_max = tmp, omic
        mod.add_omic(omic, glob_var)
    if cur_max[1]:
        return cur_max[0] - start_weight, cur_max[1]
    return -float("inf"), None


def weight_of_merged_modules(mod1, mod2, glob_var):
    start_weight = mod1.get_weight() + mod2.get_weight()
    patients_list = mod1.get_patients_names_as_list() + mod2.get_patients_names_as_list()
    omics_lists = [set(mod1.get_omics().values()) | set(mod2.get_omics().values()), set(mod1.get_omics().values()),
                   set(mod2.get_omics().values()),
                   set(mod1.get_omics().values()).intersection(set(mod2.get_omics().values()))]
    weights = [0 for i in omics_lists]
    for i in range(len(weights)):
        for omic in omics_lists[i]:
            weights[i] += omic.graph.subgraph(patients_list).size('weight')
    if max(weights) < start_weight:
        return -float("inf"), None
    return max(weights) - start_weight, (weights.index(max(weights)), (mod1, mod2))


def score_of_split_module(mod, glob_var):
    start_weight = mod.get_weight()
    adj = np.zeros((len(mod.get_patients()), len(mod.get_patients())))
    pat_list = mod.get_patients_names_as_list()
    if len(pat_list) == 0:
        raise Exception("no patients in split module {}.".format(mod))
    lst = list(mod.get_omics().items())
    lst.sort(key=lambda x: x[0])
    for name, omic in lst:
        adj += nx.adjacency_matrix(omic.graph.subgraph(pat_list), nodelist=pat_list)
    adj = pd.DataFrame(adj, index=pat_list, columns=pat_list)
    joined_subgraph = nx.from_pandas_adjacency(adj)
    heavy_subweight, heavy_subnodes = find_heaviest_subgraph(joined_subgraph, weight=mod.get_weight(), 
                                      min_size=glob_var.min_mod_size, max_size=len(pat_list) - glob_var.min_mod_size)
    if heavy_subnodes is None:
        return -float("inf"), None
    left_out = set(heavy_subnodes).symmetric_difference(set(pat_list))
    if heavy_subnodes is None or len(left_out) <= 1 or len(heavy_subnodes) <= 1:
        return -float("inf"), None
    other_weight = joined_subgraph.subgraph(left_out).size('weight')
    heavy_weight = joined_subgraph.subgraph(heavy_subnodes).size('weight')
    end_weight = heavy_weight + other_weight
    if start_weight > end_weight:
        return -float("inf"), None
    else:
        return end_weight - start_weight, (mod, heavy_subnodes)


def find_heaviest_subgraph(orig_g, weight=0, min_size=1, max_size=None):
    if max_size is None:
        max_size = len(orig_g.nodes())
    if not weight:
        weight = orig_g.size('weight')
    max_subgraph = -float("inf"), None
    g = orig_g.copy()
    graph_size = g.number_of_nodes()
    nodes_degrees = {node_name: node_degree for node_name, node_degree in g.degree(weight='weight')}
    while graph_size > 1:
        lightest_node = min(list(nodes_degrees.items()), key=lambda x: x[1])
        nodes_degrees.pop(lightest_node[0])
        for u, v in g.edges(lightest_node[0]):
            if u == lightest_node[0]: 
                nodes_degrees[v] -= g.edges[u, v]['weight']
        g.remove_node(lightest_node[0])
        graph_size -= 1
        weight -= lightest_node[1]
        cur_graph_size = len(g.nodes())
        if weight > max_subgraph[0] and cur_graph_size >= min_size and cur_graph_size <= max_size:
            max_subgraph = weight, list(g.nodes())
    return max_subgraph


def weight_of_split_and_add_omic(mod, glob_var):
    pot_omics = set(glob_var.omics.values()).symmetric_difference(set(mod.get_omics().values()))
    pats = mod.get_patients_names_as_list()
    start_weight = mod.get_weight()
    max_improvement = -float("inf"), None
    for omic in pot_omics:
        tmp_graph = omic.graph.subgraph(pats)
        tmp_weight = tmp_graph.size('weight')
        weight, sub_nodes = find_heaviest_subgraph(tmp_graph, tmp_weight,
                                      min_size=glob_var.min_mod_size, max_size=len(pats) - glob_var.min_mod_size)
        if sub_nodes is None:
            continue
        weight2 = 0
        for sub_omic in mod.get_omics().values():
            weight += sub_omic.graph.subgraph(sub_nodes).size('weight')
            weight2 += sub_omic.graph.subgraph(set(pats).symmetric_difference(set(sub_nodes))).size('weight')
        if weight + weight2 - start_weight > max_improvement[0]:
            max_improvement = (weight + weight2 - start_weight), (omic, sub_nodes)
    if max_improvement[0] > 0:
        return max_improvement
    else:
        return -float("inf"), None
        
        
def weight_of_split_and_remove_omic(mod, glob_var):
    # With 1 omic this is equivalent to module splitting.
    pot_omics = set(mod.get_omics().values())
    pats = mod.get_patients_names_as_list()
    start_weight = mod.get_weight()
    max_improvement = -float("inf"), None
    for omic in pot_omics:
        tmp_graph = omic.graph.subgraph(pats)
        tmp_weight = tmp_graph.size('weight')
        weight, sub_nodes = find_heaviest_subgraph(tmp_graph, tmp_weight,
                                      min_size=glob_var.min_mod_size, max_size=len(pats) - glob_var.min_mod_size)
        if sub_nodes is None:
            continue
        weight2 = 0
        #for sub_omic in pot_omics.symmetric_difference(set([omic])):
        for sub_omic in pot_omics:
            weight2 += sub_omic.graph.subgraph(set(pats).symmetric_difference(set(sub_nodes))).size('weight')
        if weight + weight2 - start_weight > max_improvement[0]:
            max_improvement = (weight + weight2 - start_weight), (omic, sub_nodes)
    if max_improvement[0] > 0:
        return max_improvement
    else:
        return -float("inf"), None


def weight_of_new_module(mod, glob_var):
    free_patients = []
    max_weight = -float("inf")
    max_omic = None
    max_patients = None
    for pat_name, pat in glob_var.patients.items():
        if not pat.is_in_module():
            free_patients += [pat_name]
    for omic_name, omic in glob_var.omics.items():
        tmp_weight, tmp_nodes = find_heaviest_subgraph(omic.graph.subgraph(free_patients), min_size=glob_var.min_mod_size)
        if tmp_nodes is None:
            continue
        if tmp_weight > max_weight:
            max_weight = tmp_weight
            max_omic = omic
            max_patients = tmp_nodes
    if max_weight > 0 and len(max_patients) >= 2:
        return max_weight, (max_patients, max_omic)
    else:
        return -float("inf"), None

def weight_of_spreading_module(mod, glob_var):
    if len(mod.patients) > 2 * glob_var.min_mod_size:
        return -float("inf"), None
    mod_weight  = mod.weight
    new_pat_mods = {}
    total_weight = sum([cur_mod.weight for cur_mod in glob_var.modules.values()])
    orig_mod_pats = []
    for pat in set(mod.patients.values()).copy():
        orig_mod_pats.append(pat)
        pat_best = 0, None
        mod.remove_patient(pat)
        for cur_mod_name, cur_mod in glob_var.modules.items():
            if cur_mod == mod:
                continue
            before_mod_weight = cur_mod.weight
            cur_mod.add_patient(pat)
            after_mod_weight = cur_mod.weight
            cur_mod.remove_patient(pat)
            weight_diff = after_mod_weight - before_mod_weight
            if weight_diff > pat_best[0]:
                pat_best = weight_diff, cur_mod
        mod.add_patient(pat)
        new_pat_mods[pat] = pat_best[1]

    for pat in orig_mod_pats:
        pat_new_mod = new_pat_mods[pat]
        mod.remove_patient(pat)
        if pat_new_mod is not None:
            pat_new_mod.add_patient(pat)
    new_total_weight = sum([cur_mod.weight for cur_mod in glob_var.modules.values()])
    weight_diff = new_total_weight - total_weight

    for pat in orig_mod_pats:
        pat_new_mod = new_pat_mods[pat]
        if pat_new_mod is not None:
            pat_new_mod.remove_patient(pat)
        mod.add_patient(pat)
        
    if weight_diff > 0:
        return weight_diff, new_pat_mods
    else:
        return -float("inf"), None

