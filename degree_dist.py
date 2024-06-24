import struct
import time
import sys
from os import path as osp
import numpy as np
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator
import seaborn as sns

def timer(func):
    def wrapper(*args, **kwargs):
        print(f'Function {func.__name__} starts.')
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f'Done. It takes {elapsed_time : .3f} seconds.\n')
        return result
    return wrapper

@timer
def from_txt(file_name : str):
    f = open(file = file_name, mode = 'r')
    
    # The maximum node_id used.
    max_node_id = 0
    
    # Read edges and update adjacency_list and degree_list.
    lines = f.readlines()
    row_ind = []
    col_ind = []
    for line_buffer in lines:
        # Skip the comments.
        if line_buffer[0] == '#':
            continue
        # Line format: <from_node> <to_node>
        nodes = line_buffer.split()
        # Read node ids.
        from_node = int(nodes[0])
        to_node = int(nodes[1])
        
        max_node_id = max(max_node_id, from_node)
        max_node_id = max(max_node_id, to_node)
        
        row_ind.append(from_node)
        col_ind.append(to_node)
        
        row_ind.append(to_node)
        col_ind.append(from_node)
        
    f.close()
    
    row_ind = np.array(row_ind)
    col_ind = np.array(col_ind)
    
    data = np.ones_like(row_ind)
    
    adj_sparse = csr_matrix((data, (row_ind, col_ind)))
    
    return adj_sparse

@timer
def from_bin(file_name : str):
    f = open(file = file_name, mode = 'rb')
    
    # Read whole size.
    all_size = struct.unpack('I', f.read(struct.calcsize('I')))[0]
    
    # Read all data into raw_data.
    raw_data = struct.unpack(f'{all_size}I', f.read(struct.calcsize(f'{all_size}I')))
    
    pos = 1
    row_ind = []
    col_ind = []
    
    # Read vectors one by one.
    for i in range(raw_data[0]):
        # Read vector size.
        num_neighbors = raw_data[pos]
        pos += 1
        for j in range(num_neighbors):
            row_ind.append(i)
            col_ind.append(raw_data[pos + j])
            
        pos += num_neighbors
        
    row_ind = np.array(row_ind)
    col_ind = np.array(col_ind)
    
    data = np.ones_like(row_ind)
    
    adj_sparse = csr_matrix((data, (row_ind, col_ind)))
    
    return adj_sparse

@timer
def get_degrees_fractions(degree_array_path):
    degrees = np.load(degree_array_path)
    num_nodes = len(degrees)
    # Calculate frequencies.
    unique, counts = np.unique(degrees, return_counts = True)
    fractions = counts / num_nodes
    return unique, fractions

def plot_degree_dist(unique, fractions, dataset_name):        
    plt.rcParams.update({
            "pgf.texsystem": "pdflatex",
            "font.family": "serif",
            "text.usetex": True,
            "pgf.rcfonts": False,
        })
    
    sns.set_theme(style = "ticks")
    fig = plt.figure(figsize = (4, 4))
    ax = fig.add_subplot(111)
    
    ax.scatter(unique, fractions, marker = '.', s = 50, color = 'red', edgecolor = 'none', alpha = .5)
    
    if dataset_name != 'synthetic':
        ax.set_xscale('log')
        min_ticks = 5
        ax.xaxis.set_major_locator(LogLocator(base = 10.0, numticks = max(10, min_ticks)))
        ax.xaxis.set_minor_locator(LogLocator(base = 10.0, subs = [2, 3, 4, 5, 6, 7, 8, 9], numticks = max(10, min_ticks)))
    ax.set_yscale('log')
    ax.tick_params(axis = 'x', labelsize = 14)
    ax.tick_params(axis = 'y', labelsize = 14)
    ax.set_ylabel('Fraction', fontsize = 14)
    ax.set_xlabel('Degree', fontsize = 14)

    ax.grid(axis = 'y', linestyle = '--', alpha = 0.7)
    ax.set_title(dataset_name if dataset_name != 'synthetic' else r'Erd\H{o}s-R{\'e}nyi', fontsize = 14)
    
    # Make it square.
    ax.set_aspect(1.0 / ax.get_data_ratio(), adjustable = 'box')
    plt.savefig(f'./results/degree-dist/{dataset_name}.pdf', bbox_inches = 'tight')
    plt.show()
    
if __name__ == '__main__':
    arguments = sys.argv
    if len(arguments) < 2:
        print('Lack arguments!')
        raise
    dataset_name = arguments[1]
    compressed_sorted_bin_path = './datasets/' + dataset_name + '-compressed_sorted.bin'
    degree_array_path = './datasets/' + dataset_name + '-degree_list.npy'
    
    if osp.exists(degree_array_path):
        print('Using the cached degree array.')
        unique, fractions = get_degrees_fractions(degree_array_path)
    elif osp.exists(compressed_sorted_bin_path):
        print('Reading graphs.')
        adj_sparse = from_bin(compressed_sorted_bin_path)
        num_nodes = adj_sparse.shape[0]
        # Calculate degrees.
        degrees = np.array(adj_sparse.sum(axis = 1)).flatten()
        # Save file.
        np.save(degree_array_path, degrees)
        # Calculate frequencies.
        unique, counts = np.unique(degrees, return_counts = True)
        fractions = counts / num_nodes
    else:
        print('No graph file, please run main.cpp first to get graph in compressed binary format.')
        raise
        
    plot_degree_dist(unique, fractions, dataset_name)