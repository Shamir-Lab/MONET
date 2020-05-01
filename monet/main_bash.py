"""
Module used to run MONET from the command line. This is used by the R code that executes MONET, 
but it is not required when running MONET from python.
"""

import sys
import os.path

import random
import numpy as np
import pandas as pd
import networkx as nx
import argparse
from monet import data_prep
from monet.monet import Monet

def export_monet_ret(monet_ret, output_dir_path):
    res, total_time, super_g, iteration_times, total_weight = monet_ret
    exp_string = ''
    mod_string = ''
    patients_with_modules = []
    for mod_id, module in res.modules.items():
        for patient in module.patients:
            patients_with_modules.append(patient)
            exp_string = exp_string + patient + ' ' + str(mod_id) + '\n'
        mod_string = mod_string + str(mod_id) + ' ' + ','.join(list(module.get_omics().keys())) + '\n'

    for omic_id, omic in res.omics.items():
        pd.DataFrame(nx.to_numpy_matrix(omic.graph)).to_csv(os.path.join(output_dir_path, 'graph_' + omic_id))

    pats_without_mods = set(res.patients.keys()) - set(patients_with_modules)
    for patient in pats_without_mods:
        exp_string = exp_string + patient + ' None\n'

    with open(os.path.join(output_dir_path, 'pats'), 'wb') as f:
        f.write(bytes(exp_string, 'ascii'))
    with open(os.path.join(output_dir_path, 'time'), 'wb') as f:
        f.write(bytes(str(total_time), 'ascii'))
    with open(os.path.join(output_dir_path, 'weight'), 'wb') as f:
        f.write(bytes(str(total_weight), 'ascii'))
    with open(os.path.join(output_dir_path, 'mods'), 'wb') as f:
        f.write(bytes(mod_string, 'ascii'))
    with open(os.path.join(output_dir_path, 'iter_times'), 'wb') as f:
        iter_times_string = '\n'.join([str(x) for x in iteration_times])
        f.write(bytes(iter_times_string, 'ascii'))
    
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir_path", help="directory which includes all omic files")
    parser.add_argument("--output_dir_path", help="directory where the output will be saved")
    parser.add_argument("--num_iters", help="number of iterations run by the algorithm", type=int, default=1e3)
    parser.add_argument("--num_seeds", help="number of seeds to use", type=int, default=30)
    parser.add_argument("--min_mod_size", help="minimal module size", type=int, default=10)
    parser.add_argument("--max_pats_per_action", help="maximal number of samples to be added to  a module / switched between modules in a single action", type=int, default=10)
    parser.add_argument("--seed_size", help="number of samples in a seed", type=int, default=15)
    parser.add_argument("--rand_seed", help="value used as seed for randomization", type=int, default=42)
    parser.add_argument("--input_type", help="the type of input - either similarity data (default) or raw", type=str, default='similarity', choices=['similarity', 'raw'])
    parser.add_argument("--percentile_shift", help="the percentile shifted to weight edge (default - no shift. Should be between 0 and 100)", type=float, default=None)
    parser.add_argument("--percentile_remove_edge", help="A percentile, in which case positive and negative edges under this percentile are removed", type=float, default=None)
    args = parser.parse_args()
    dir_path = args.dir_path
    output_dir_path = args.output_dir_path
    num_seeds = int(args.num_seeds)
    seed_size = int(args.seed_size)
    rand_seed = int(args.rand_seed)
    num_iters = int(args.num_iters)
    min_mod_size = int(args.min_mod_size)
    max_pats_per_action = int(args.max_pats_per_action)
    assert not ((args.percentile_shift is not None) and (args.input_type == 'similarity'))
    if args.percentile_shift is not None:
        percentile_shift = float(args.percentile_shift)
    else:
        percentile_shift = None
    is_input_raw = args.input_type  == 'raw'

    np.random.seed(rand_seed)
    random.seed(rand_seed)

    dist_mats, raw_data = data_prep.process_data(path=os.path.join(dir_path, 'dist'), is_input_raw=is_input_raw)
    init_file_name = os.path.join(dir_path, 'init', 'init')
    if os.path.exists(init_file_name):
        npm = np.loadtxt(fname=init_file_name, delimiter=',', dtype=str)
        init_modules = {}
        for i in range(npm.shape[0]):
            if npm[i, 1] not in init_modules:
                init_modules[npm[i, 1]] = []
            init_modules[npm[i, 1]].append(npm[i, 0])
    else:
        init_modules = None

    monet_ret = Monet().main_loop(data=dist_mats, is_input_raw=is_input_raw, init_modules=init_modules, iters=num_iters, min_mod_size=min_mod_size, max_pats_per_action=max_pats_per_action,
                                  num_of_samples_in_seed=seed_size, num_of_seeds=num_seeds, percentile_shift=percentile_shift, percentile_remove_edge=args.percentile_remove_edge)
    export_monet_ret(monet_ret, output_dir_path)

