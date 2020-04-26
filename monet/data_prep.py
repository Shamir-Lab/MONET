"""
Module contains code to read matrices written by the R code that executes MONET.
It is not required if running MONET from python.
"""

from sklearn.metrics.cluster import adjusted_rand_score
import numpy as np
import pandas as pd
from monet import Monet
import os
from scipy import stats
import random
import operator
import networkx as nx
#import matplotlib.pyplot as plt
import time

def expression_matrix_to_similarity(mat):
    """gets an expression matrix (single omic) and returns dist matrix and
    correlation between samples"""
    from scipy.spatial.distance import squareform, pdist
    import math
    
    df = pd.DataFrame.corr(mat)
    np.fill_diagonal(df.values, 0)
    return df


def standardize_data(df):
    """standardize the data"""
    return df.mean(1), df.std(1), df.sub(df.mean(1), axis=0).div(df.std(1), axis=0)


def readfiles(path):
    data = {}
    for filename in os.listdir(path):
        full_filename = os.path.join(path, filename)
        if os.path.isfile(full_filename):
            npm = np.loadtxt(fname=full_filename, delimiter=',', dtype=str)
            npm = pd.DataFrame(data=pd.to_numeric(npm[1:, 1:].flatten(), errors='coerce').reshape(npm.shape[0]-1, npm.shape[1]-1), columns=npm[0, 1:],
                               index=npm[1:, 0])
            data[filename.split(".")[1]] = npm
    # del data[2] # negative values
    return data

def process_data(mats=None, path=r'', is_input_raw=False):
    """gets raw data or a directory to raw data and return a list of similarity matrices"""
    if not mats:
        mats = readfiles(path)
    dist_matrices = {}
    for omic, mat in mats.items():
        if is_input_raw:
            mean1, std1, stand = standardize_data(mat)
            dist_matrices[omic] = expression_matrix_to_similarity(stand)
        else:
            dist_matrices[omic] = mat

    return dist_matrices, mats




