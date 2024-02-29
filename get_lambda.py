import struct
import time
from os import path as osp
import numpy as np
from scipy.sparse import csr_matrix, diags
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import eigs, eigsh
import sys

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
def get_P_lambda(A : csr_matrix):
	degree_vector = np.asarray(A.sum(axis = 1)).flatten()
	D_inv = diags(1.0 / degree_vector, format = 'csr')
	P = A @ D_inv

	w, v = eigs(P, k = 2, which = 'LM')
	return np.abs(w[1])

@timer
def get_L_lambda_2(A : csr_matrix):
	D = diags(np.asarray(A.sum(axis = 1)).flatten(), format = 'csr')
	L = D - A
	
	w, v = eigsh(L, k = 2, which = 'SA')
	return np.abs(w[1])

@timer
def get_L_lambda_n(A : csr_matrix):
	D = diags(np.asarray(A.sum(axis = 1)).flatten(), format = 'csr')
	L = D - A
	
	w, v = eigsh(L, k = 1, which = 'LA')
	return np.abs(w[0])

if __name__ == '__main__':
	arguments = sys.argv
	if len(arguments) < 2:
		print('Lack arguments!')
		raise
	dataset_name = arguments[1]
	compressed_sorted_bin_path = './datasets/' + dataset_name + '-compressed_sorted.bin'
	lambda_path = './datasets/' + dataset_name + '.lambda'

	if osp.exists(compressed_sorted_bin_path):
		adj_sparse = from_bin(compressed_sorted_bin_path)
	else:
		print('No graph file, please run main.cpp first to get graph in compressed binary format.')
		raise

	# A = get_largest_connected_component(adj_sparse = adj_sparse)
	A = adj_sparse

	lambda_n = get_L_lambda_n(A)
	print('lambda_n = ', lambda_n)
	
	lambda_2 = get_L_lambda_2(A)
	print('lambda_2 = ', lambda_2)
	
	print('eps = 1e-2, L_max = ', np.log(100 * np.power(lambda_n, 0.5) / np.power(lambda_2, 1.5)))
	