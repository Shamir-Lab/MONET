import heapq
import math
import networkx as nx
import numpy as np
import pandas as pd
import time
import os
from monet.patient import Patient
from monet.module import Module
from monet.omic import Omic
from monet.globals import Globals
from monet.module_score import *
import random

class Monet:
    """
    Represents a single MONET run.
    """
    # a list of the actions considered by MONET in each iteration. Each action correponds to one function in the list.
    functions = [best_samples_to_add, which_patient_to_remove, which_omic_to_add_to_module,
                 which_omic_to_remove_from_module, score_of_split_module, weight_of_split_and_add_omic, weight_of_split_and_remove_omic,
                 weight_of_new_module, top_samples_to_switch, weight_of_spreading_module, weight_of_merged_modules]
    functions_names = ["best_samples_to_add", "which_patient_to_remove", "which_omic_to_add_to_module",
                       "which_omic_to_remove_from_module", "score_of_split_module", "weight_of_split_and_add_omic", "weight_of_split_and_remove_omic",
                       "weight_of_new_module", "top_samples_to_switch", "weight_of_spreading_module", "weight_of_merged_modules"]

    def build_a_graph_similarity(self, distances):
        res = distances.values
        np.fill_diagonal(res, 0)
        g = nx.from_numpy_matrix(res)
        mapping = {}
        for i in range(len(distances.columns)):
            mapping[i] = distances.columns[i]
        nx.relabel_nodes(g, mapping, False)
        return g, [], [], 0

    def create_env(self, glob_var, data, is_input_raw, percentile_shift, percentile_remove_edge):
        """
	Create all the variables used during MONET's run: 
	modules, omics etc, and associating them with a Global instance.
	"""
        num_omics = len(data.values())
        all_pat_names = set()
        for i in range(num_omics):
            all_pat_names.update(list(data.values())[i].columns.values)
        for pat in all_pat_names:
            glob_var.patients.update({pat:Patient(pat)})
            
        for omic, dat in data.items():
            self.omic = omic
            if is_input_raw:
                assert False, "only similarity input matrices are corrently supported"
            else:
                graph, means, covs, percentile = self.build_a_graph_similarity(dat)
            if percentile_remove_edge is not None:
                all_weights = []
                for edge in graph.edges:
                    all_weights.append(graph.edges[edge]['weight'])
                all_weights_array = np.array(all_weights)
                positive_thresh = np.percentile(all_weights_array[all_weights_array > 0], percentile_remove_edge)
                negative_thresh = np.percentile(all_weights_array[all_weights_array < 0], 100 - percentile_remove_edge)
                all_edges = []
                for edge in graph.edges:
                    all_edges.append(edge)
                for edge in all_edges:
                    cur_weight = graph.edges[edge]['weight']
                    if (cur_weight > 0 and cur_weight < positive_thresh) or (cur_weight < 0 and cur_weight > negative_thresh):
                        graph.remove_edge(edge[0], edge[1])
                
            cur_graph_pats = set(graph.nodes)
            missing_pats = all_pat_names - cur_graph_pats
            for missing_pat in missing_pats:
                graph.add_node(missing_pat)
            glob_var.omics.update({omic: Omic(graph=graph, name=omic)})
            glob_var.gmm_params.update({omic: {'mean': means, 'cov': covs, 'percentile': percentile}})
        return glob_var

    def get_seeds(self, glob_var, num_of_seeds=3, num_of_samples_in_seed=10):
        """
	Create seed modules.
	"""
        lst = list(glob_var.omics.items())
        lst.sort(key=lambda x: x[0])
        omics_list = [omic for name, omic in lst]
        pat_list = list(glob_var.patients.keys())
        adj = np.zeros((len(pat_list), len(pat_list)))
        for name, omic in lst:
            adj += nx.adjacency_matrix(omic.graph.subgraph(pat_list), nodelist=pat_list)
        adj = pd.DataFrame(adj, index=pat_list, columns=pat_list)
        joined_subgraph = nx.from_pandas_adjacency(adj)
        omic_graphs = [joined_subgraph]

        for i in range(num_of_seeds):

            omic_graph = omic_graphs[0]
            cur_nodes = list(sorted(omic_graph.nodes()))
            adj = list(omic_graph.adjacency())

            if len(cur_nodes) == 0:
                break

            rand_pat_index = random.randint(0, len(cur_nodes) - 1)
            rand_pat_name = cur_nodes[rand_pat_index]
            rand_pat_in_adj = [pat[0] for pat in adj].index(rand_pat_name)
            neighbors = [(key, adj[rand_pat_in_adj][1][key]['weight']) for key in adj[rand_pat_in_adj][1]]
            neighbors = sorted(neighbors, key=lambda x: x[1], reverse=True)[:(num_of_samples_in_seed - 1)]
            nodes = {rand_pat_name: glob_var.patients[rand_pat_name]}
            for nei in neighbors:
                if nei[1] > 0 and nei[0] != rand_pat_name:
                    nodes.update({nei[0]:glob_var.patients[nei[0]]})
            mod_weight = omic_graph.subgraph(list(nodes.keys())).size('weight')
            if mod_weight > 0 and len(nodes) > 1 and len(nodes) >= glob_var.min_mod_size:
                Module(glob_var=glob_var, patients=nodes, omics=[omic for omic in omics_list], weight=mod_weight)
                for k in range(len(omic_graphs)):
                    omic_graph = omic_graphs[k]
                    remaining_nodes = list(sorted(set(cur_nodes) - set(nodes.keys())))
                    omic_graphs[k] = omic_graph.subgraph(remaining_nodes)
        return glob_var


    def get_next_step(self, mod, glob_var):
        """
	this function decided what is the next action that will be executed.
	"""
        max_res = (-1, (-float("inf"), None))
        for func_i in range(len(self.functions)):
            if func_i <= 9:  # only one module needed
                tmp = self.functions[func_i](mod, glob_var)
                if tmp[0] > max_res[1][0]:
                    max_res = (func_i, tmp)
            else:
                for mod2 in glob_var.modules.values():
                    if mod2 == mod:
                        continue
                    tmp = self.functions[func_i](mod, mod2, glob_var)
                    if tmp[0] > max_res[1][0]:
                        max_res = (func_i, tmp)
        return max_res


    def exec_next_step(self, mod, max_res, glob_var):
        """
	this function actually performs an action, given that the 
	algorithm already decided what the next action will be.
	"""
        if max_res[1][0] == -float("inf") or max_res[1][0] < 0:
            return glob_var
        func_i = max_res[0]
        glob_var.actions[func_i] += 1
        if func_i == 0:  # add
            for sample in max_res[1][1]:
                mod.add_patient(sample)
        elif func_i == 1:  # remove
            mod.remove_patient(max_res[1][1])
            if len(mod.get_patients()) <= 1:
                glob_var = glob_var.kill_module(mod)
        elif func_i == 2:  # add omic
            mod.add_omic(max_res[1][1], glob_var)
        elif func_i == 3:  # remove omic
            mod.remove_omic(max_res[1][1], glob_var)
        elif func_i == 4:  # split
            glob_var = mod.split_module(max_res[1][1][1], glob_var)
        elif func_i == 5:  # split and add omic
            glob_var = mod.split_and_add_omic(omic=max_res[1][1][0], sub_nodes=max_res[1][1][1], glob_var=glob_var)
        elif func_i == 6:  # split and remove omic
            glob_var = mod.split_and_remove_omic(omic=max_res[1][1][0], sub_nodes=max_res[1][1][1], glob_var=glob_var)
        elif func_i == 7:  # create new module
            new_mod = Module(glob_var)
            new_mod.add_omic(max_res[1][1][1], glob_var)
            for pat in max_res[1][1][0]:
                new_mod.add_patient(glob_var.patients[pat])
        elif func_i == 8:  # transfer
            pats = [(pat, mod2) for pat, weight, mod2 in max_res[1][1]]
            for pat, mod2 in pats:
                switch_2_patients(glob_var.patients[pat], mod, mod2, glob_var)
        elif func_i == 9:  # spread module
            mod.spread_module(max_res[1][1], glob_var)
        elif func_i == 10:  # merge
            glob_var = mod.merge_with_module(max_res[1], glob_var)
        return glob_var

    def create_seeds_from_solution(self, glob_var, init_modules):
        for mod_name, pat_ids in init_modules.items():
            omics = glob_var.omics
            pat_dict = {}
            for pat_id in pat_ids:
                pat_dict[pat_id] = glob_var.patients[pat_id]
            mod_weight = 0
            for omic in omics.values():
                mod_weight += omic.graph.subgraph(list(pat_dict.keys())).size('weight')
            Module(glob_var=glob_var, patients=pat_dict, omics=omics, weight=mod_weight)
        return glob_var

    def main_loop(self, data, is_input_raw=False, init_modules=None, iters = 500, 
                  num_of_seeds=10, num_of_samples_in_seed=10, min_mod_size=10, max_pats_per_action=10, percentile_shift=None, percentile_remove_edge=None):
        """
	:param data: the omic graphs in which MONET will detect modules. data should be a dict mapping omic names to the omic graphs.
	Each omic graph is a pandas DataFrame with shape (num_samples, num_samples), where each entry in the dataframe signifies the similarity
	between a pair of samples. Column and row names should be the sample ids.
	:param is_input_raw: True <=> the input data is the raw feature values (instead of the omic graphs). Currently only omic graphs inputs are supported,
	so this value has to be False.
	:type is_input_raw: bool.
	:param init_modules: an optional module initialization for MONET. A dict mapping between module names to sample ids (as appear in the feature names in "data").
	All modules are initialized to cover all omics. Set to None to use MONET's seed finding algorithm for initialization.
	:param iters: maximal number of iterations.
	:type iters: maximal number of iterations.
	:param num_of_seeds: number of seeds to create in MONET's module initialization algorithm.
	:type num_of_seeds: int.
	:param num_of_samples_in_seed: number of samples to put in each seeds to create in MONET's module initialization algorithm.
	:type num_of_samples_in_seed: int.
	:param min_mod_size: minimal size (number of samples) for a MONET module.
	:type min_mod_size: int.
	:param max_pats_per_action: maximal number of samples in a single MONET action 
	(maximal number of samples added to a module or replaced between modules in a single action).
	:type max_pats_per_action: int.
	:param percentile_shift: not currently being used.
	:param percentile_remove_edge: only edges with weight percentile above (for positive weights) or below (for negative weights) this percentile
	are kept in the graph. For example, percentile_remove_edge=90 keeps only the 10% edges with highest positive weight and lowest negative weight in the graph.
	None keeps all edges in the graph.
	:type percentile_remove_edge: int.
	:return: MONET's detected modules in the data. The returned object is a dict mapping module names to Module objects. Every module instance includes its
	set of patients (under the "patients" attribute) and its set of omics (the "omics" attribute).

	Note that MONET uses both the numpy random seed and the "random" module's random seed. There, to we recommend you set both seeds before executing the algorithm:
	numpy.random.seed(rand_seed)
	random.seed(rand_seed)
	"""
        start_time = time.time()
        glob_var = Globals(len(self.functions))
        glob_var = self.create_env(glob_var, data, is_input_raw, percentile_shift, percentile_remove_edge)
        glob_var.min_mod_size = min_mod_size
        glob_var.max_pats_per_action = max_pats_per_action

        if init_modules is None:
            glob_var = self.get_seeds(glob_var, num_of_seeds=num_of_seeds, num_of_samples_in_seed=num_of_samples_in_seed)
        else:
            glob_var = self.create_seeds_from_solution(glob_var, init_modules)

        for _, some_mod in glob_var.modules.copy().items():
            if len(some_mod.patients) < min_mod_size:
                print('killing a small module before starting')
                glob_var.kill_module(some_mod)

        total_weight = sum(mod.get_weight() for mod in glob_var.modules.values())
        converged_modules = {}
        did_action = False
        mod_index = 0
        iterations = 0

        iteration_times = []
        while iterations < iters:
            iteration_start_time = time.time()
            prev_weight = total_weight

            active_module_names = list(sorted(set(glob_var.modules.keys()) - set(converged_modules.keys())))
            if len(active_module_names) == 0:
                if not did_action:
                    print("converged, total score: {}.".format(total_weight))
                    break
                else:
                    converged_modules = {}
                    did_action = False
                    active_module_names = list(sorted(glob_var.modules.keys()))
            mod_name = random.choice(active_module_names)
            mod = glob_var.modules[mod_name]

            max_res = self.get_next_step(mod, glob_var)
            print(str(mod) + str(max_res))
            glob_var = self.exec_next_step(mod, max_res, glob_var)
            for _, some_mod in glob_var.modules.copy().items():
                if len(some_mod.get_patients()) <= 1 or not some_mod.get_omics():
                    glob_var.kill_module(some_mod)
                    print('removing zombie module')

            total_weight = sum([mod.get_weight() for name, mod in glob_var.modules.items()])
            iterations += 1
            if iterations % 10 == 0:
                print("iteration: " + str(iterations))
                print("num of modules: " + str(len(glob_var.modules)))
                print("total_weight: " + str(total_weight))
                print("actions: " + str(glob_var.actions))

            # Assert module sizes
            for _, some_mod in glob_var.modules.copy().items():
                assert len(some_mod.patients) >= min_mod_size

            if total_weight <= prev_weight or max_res[1][0] == -float("inf"):
                if mod_name in glob_var.modules:
                    converged_modules.update({mod_name: glob_var.modules[mod_name]})
            else: # the score deviates from the score we expected
                if not (abs(total_weight - prev_weight - max_res[1][0]) < 0.01):
		    # This signifies a bug and should never occur:
		    # that the difference in the objective function from the 
		    # previous iteration is different from the difference
		    # the algorithm expected for the function.
                    import pdb;pdb.set_trace()
                did_action = True
                assert abs(total_weight - prev_weight - max_res[1][0]) < 0.01
                did_action = True
            iteration_times.append(time.time() - iteration_start_time)
        for mod_name, mod in glob_var.modules.copy().items():
            if mod.get_size() <= glob_var.min_mod_size and not is_mod_significant(mod, glob_var):
                print("module {} with patients {} on omics {} is not significant.".format((mod_name, mod), mod.get_patients(), mod.get_omics().keys()))
                glob_var.kill_module(mod)
        total_time = time.time() - start_time
        return glob_var, total_time, None, iteration_times, total_weight


def is_mod_significant(mod, glob_var, percentile=95, iterations=500):
    """
    Assess the statisitcal significance of a module by sampling modules or similar size.
    """
    draws = [0 for i in range(iterations)]
    mod_size = len(mod.get_patients())
    if mod_size <= 1:
        return False
    for i in range(iterations):
        samps = random.sample(glob_var.patients.keys(), mod_size)
        lst = list(mod.get_omics().items())
        lst.sort(key=lambda x: x[0])
        for name, omic in lst:
            draws[i] += omic.graph.subgraph(samps).size('weight')
    num_to_beat = np.percentile(draws, percentile)
    if mod.get_weight() > num_to_beat:
        return True
    return False

