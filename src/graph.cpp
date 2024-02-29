#include "../include/graph.h"

Graph::Graph() : num_nodes(0), num_edges(0), d_min(0), d_max(0), max_degree_node(0), transition_matrix(), is_transition_matrix_valid(false)
{
}

bool Graph::operator==(const Graph & other)
{
	return this -> adjacency_list == other.adjacency_list;
}

// Build graph from given file (const char *)
bool Graph::from_txt(const char * file_name)
{
	// Clear.
	this -> num_nodes = size_type(0);
	this -> num_edges = size_type(0);
	this -> d_min = degree_type(0);
	this -> d_max = degree_type(0);
	this -> max_degree_node = size_type(0);
	this -> transition_matrix = Eigen::SparseMatrix<double>();
	this -> is_transition_matrix_valid = false;

	this -> adjacency_list.clear();
	// this -> degree_list.clear();

	this -> adjacency_list.resize(100000000, std::vector<size_type>());
	// this -> degree_list.resize(100000000, degree_type(0));

	// Close syncing with stdio.
	std::ios::sync_with_stdio(false);
	// Open the file as a stream.
	std::ifstream file_stream(file_name);
	// File not found.
	if(!file_stream.is_open())
	{
		file_stream.close();
		std::ios::sync_with_stdio(true);
		return false;
	}

	// The maximum node_id used.
	size_type max_node_id = size_type(0);

	{
		// Create a string stream and read all contents into the stream.
		std::stringstream string_stream;
		string_stream << file_stream.rdbuf();
		// Close file.
		file_stream.close();

		// Read edges and update adjacency_list and degree_list.
		for(std::string line_buffer; std::getline(string_stream, line_buffer); )
		{
			// Skip the comments.
			if(line_buffer[0] == '#')
				continue;
			// Read node ids.
			std::istringstream line_stream(line_buffer);
			size_type from_node = size_type(0), to_node = size_type(0);
			// Line format: <from_node> <to_node>
			line_stream >> from_node >> to_node;

			max_node_id = std::max(max_node_id, from_node);
			max_node_id = std::max(max_node_id, to_node);
			
			// Increase size.
			if(from_node >= this -> adjacency_list.size() 
				|| to_node >= this -> adjacency_list.size()
				// || from_node >= this -> degree_list.size()
				// || to_node >= this -> degree_list.size()
			)
			{
				this -> adjacency_list.resize(2 * std::max(from_node, to_node), std::vector<size_type>());
				// this -> degree_list.resize(2 * std::max(from_node, to_node), degree_type(0));
			}

			this -> adjacency_list[from_node].push_back(to_node);
			this -> adjacency_list[to_node].push_back(from_node);

			this -> num_edges += 1;

			// this -> degree_list[from_node] += 1;
			// this -> degree_list[to_node] += 1;
		}
	}
	
	// Resize to fit.
	this -> adjacency_list.resize(max_node_id + 1);
	// this -> degree_list.resize(max_node_id + 1);

	// Count num_nodes.
	this -> num_nodes = Graph::size_type(0);
	for(const auto & v : this -> adjacency_list)
	{
		if(!v.empty())
			this -> num_nodes++;
	}

	// Get degree parameters.
	this -> d_max = Graph::size_type(0);
	this -> max_degree_node = Graph::size_type(0);
	for(Graph::size_type i = 0; i < this -> adjacency_list.size(); i++)
	{
		Graph::degree_type d = static_cast<Graph::degree_type>(this -> adjacency_list[i].size());
		// Get node with maximum degree.
		if(d > this -> d_max)
		{
			this -> d_max = d;
			this -> max_degree_node = i;
		}
	}
	this -> d_min = d_max;
	for(const auto & v : this -> adjacency_list)
	{
		Graph::degree_type d = static_cast<Graph::degree_type>(v.size());
		if(v.size() > degree_type(0))
			this -> d_min = std::min(d, this -> d_min);
	}
	
	std::ios::sync_with_stdio(true);
	return true;
}

bool Graph::to_bin(const char * file_name)
{
	// Close syncing with stdio.
	std::ios::sync_with_stdio(false);

	std::ofstream file_stream(file_name, std::ios::out | std::ifstream::binary);
	if(!file_stream.is_open())
	{
		file_stream.close();
		// Resyncing with stdio.
		std::ios::sync_with_stdio(true);
		return false;
	}
	// Store total number of numbers (adjacency_list_size itself + adjacency_list.size() + 2 * num_edges) to file.
	// Format: <all_size = 1 + n + 2m> <adj_list_size> <adj_list[0].size()> <adj_list[0][0]> <adj_list[0][1]> ... <adj_list[1].size()> ...
	Graph::size_type all_size = 1 + this -> adjacency_list.size() + 2 * this -> num_edges;
	file_stream.write(reinterpret_cast<const char *>(&all_size), sizeof(Graph::size_type));

	// Store adjacency_list to file.
	// Store size of the outer vector.
	Graph::size_type adjacency_list_size = this -> adjacency_list.size();
	file_stream.write(reinterpret_cast<const char *>(&adjacency_list_size), sizeof(Graph::size_type));    

	// Now write each vector one by one.
	for(auto & v : this -> adjacency_list)
	{
		// Store its size.
		Graph::size_type size = v.size();
		file_stream.write(reinterpret_cast<const char *>(&size), sizeof(Graph::size_type));

		// Store its contents.
		file_stream.write(reinterpret_cast<const char *>(&v[0]), v.size() * sizeof(Graph::size_type));
	}
	file_stream.close();
	// Resyncing with stdio.
	std::ios::sync_with_stdio(true);
	return true; 
}

bool Graph::from_bin(const char * file_name)
{
	// Clear.
	this -> num_nodes = size_type(0);
	this -> num_edges = size_type(0);
	this -> d_min = degree_type(0);
	this -> d_max = degree_type(0);
	this -> max_degree_node = size_type(0);
	this -> transition_matrix = Eigen::SparseMatrix<double>();
	this -> is_transition_matrix_valid = false;

	this -> adjacency_list.clear();
	// this -> degree_list.clear();

	// Close syncing with stdio.
	std::ios::sync_with_stdio(false);

	std::ifstream file_stream(file_name, std::ios::in | std::ifstream::binary);
	if(!file_stream.is_open())
	{
		file_stream.close();
		// Resyncing with stdio.
		std::ios::sync_with_stdio(true);
		return false;
	}

	// Read whole size.
	Graph::size_type all_size = Graph::size_type(0);
	file_stream.read(reinterpret_cast<char *>(&all_size), sizeof(Graph::size_type));

	// Read all data into raw_data.
	std::unique_ptr<Graph::size_type> raw_data(new Graph::size_type[all_size]);

	// Read size of adjacency_list.
	file_stream.read(reinterpret_cast<char *>(&raw_data.get()[0]), sizeof(Graph::size_type) * all_size);

	// Reserve space.
	this -> adjacency_list.resize(raw_data.get()[0], std::vector<Graph::size_type>());
	// this -> degree_list.resize(raw_data[0], Graph::degree_type(0));
	Graph::size_type pos = Graph::size_type(1);

	// Read vectors one by one.
	for(Graph::size_type i = Graph::size_type(0); i < raw_data.get()[0]; i++)
	{
		// Read vector size.
		Graph::size_type num_neighbors = raw_data.get()[pos++];
		this -> num_edges += num_neighbors;

		this -> adjacency_list[i].resize(num_neighbors, Graph::size_type(0));
		// this -> degree_list[i] = reinterpret_cast<Graph::degree_type>(num_neighbors);

		// Read each element.
		Graph::size_type node = Graph::size_type(0);
		for(Graph::size_type j = Graph::size_type(0); j < num_neighbors; j++)
		{
			node = raw_data.get()[pos++];
			this -> adjacency_list[i][j] = node;
		}
	}

	// Release raw_data.
	raw_data = nullptr;

	// Note: each edge is counted twice.
	this -> num_edges /= Graph::size_type(2);

	// Count num_nodes.
	this -> num_nodes = Graph::size_type(0);
	for(const auto & v : this -> adjacency_list)
	{
		if(!v.empty())
			this -> num_nodes++;
	}

	// Get degree parameters.
	this -> d_max = Graph::size_type(0);
	this -> max_degree_node = Graph::size_type(0);
	for(Graph::size_type i = 0; i < this -> adjacency_list.size(); i++)
	{
		Graph::degree_type d = static_cast<Graph::degree_type>(this -> adjacency_list[i].size());
		// Get node with maximum degree.
		if(d > this -> d_max)
		{
			this -> d_max = d;
			this -> max_degree_node = i;
		}
	}
	this -> d_min = d_max;
	for(const auto & v : this -> adjacency_list)
	{
		Graph::degree_type d = static_cast<Graph::degree_type>(v.size());
		if(v.size() > degree_type(0))
			this -> d_min = std::min(d, this -> d_min);
	}

	// Close file.
	file_stream.close();
	// Resyncing with stdio.
	std::ios::sync_with_stdio(true);
	return true; 
}

bool Graph::compress_adjacency_list()
{
	// Clear.
	this -> is_transition_matrix_valid = false;

	// Map which maps old node_id to new_id.
	std::map<Graph::size_type, Graph::size_type> map_to_new;
	Graph::size_type new_id = 0;

	// Get the mapping.
	for(Graph::size_type i = 0; i < this -> adjacency_list.size(); i++)
	{
		if(this -> adjacency_list[i].size() > 0)
		{
			map_to_new[i] = new_id;
			new_id++;
		}
	}
	// Move.
	for(Graph::size_type i = 0; i < this -> adjacency_list.size(); i++)
	{
		if(this -> adjacency_list[i].size() > 0)
		{
			std::for_each(this -> adjacency_list[i].begin(), this -> adjacency_list[i].end(), [&](Graph::size_type & a){a = map_to_new[a];});

			if(map_to_new[i] != i)
			{
				this -> adjacency_list[map_to_new[i]] = (this -> adjacency_list[i]);
				this -> adjacency_list[i].clear();
			}
			// this -> num_edges += this -> adjacency_list[map_to_new[i]].size();
		}
	}

	// Get max degree node.
	this -> d_max = Graph::size_type(0);
	this -> max_degree_node = Graph::size_type(0);
	for(Graph::size_type i = 0; i < this -> adjacency_list.size(); i++)
	{
		Graph::degree_type d = static_cast<Graph::degree_type>(this -> adjacency_list[i].size());
		// Get node with maximum degree.
		if(d > this -> d_max)
		{
			this -> d_max = d;
			this -> max_degree_node = i;
		}
	}

	return true;
}

bool Graph::sort_adjacency_list()
{
	for(auto & v : this -> adjacency_list)
	{
		std::sort(v.begin(), v.end(), [&](const Graph::size_type & v1, const Graph::size_type & v2) -> bool {return this -> adjacency_list[v1].size() < this -> adjacency_list[v2].size();});
	}
	return true;
}

void Graph::print_adjacency_list()
{
	for(Graph::size_type i = 0; i < this -> adjacency_list.size(); i++)
	{
		std::cout << "Node " << i << ": ";
		for(const auto & n : this -> adjacency_list[i])
		{
			std::cout << n << " ";
		}
		std::cout << std::endl;
	}
}

void Graph::print_statistics()
{
	std::cout << "Statistics:" << std::endl;
	std::cout << "\tnum_nodes: " << this -> get_num_nodes() << std::endl;
	std::cout << "\tnum_edges: " << this -> get_num_edges() << std::endl;
	std::cout << "\td_min: " << this -> get_d_min() << std::endl;
	std::cout << "\td_max: " << this -> get_d_max() << std::endl;
	std::cout << "\td_avg: " << this -> get_d_avg() << std::endl;
}

void Graph::init_transition_matrix()
{
	this -> transition_matrix = Eigen::SparseMatrix<double>(this -> adjacency_list.size(), this -> adjacency_list.size());
	std::vector<Eigen::Triplet<double, int64_t>> triplets;
    for(Graph::size_type i = 0; i < this -> adjacency_list.size(); i++)
	{
        for(const Graph::size_type & j : this -> adjacency_list[i])
		{
            triplets.push_back({i, j, 1.0 / this -> adjacency_list[j].size()});
        }
    }

	this -> transition_matrix.setFromTriplets(triplets.begin(), triplets.end());
	this -> is_transition_matrix_valid = true;
}

// Get graph parameters.
Graph::size_type Graph::get_num_nodes() const
{
	return this -> num_nodes;
}

Graph::size_type Graph::get_num_edges() const
{
	return this -> num_edges;
}

Graph::degree_type Graph::get_d_min() const
{
	return this -> d_min;
}

Graph::degree_type Graph::get_d_max() const
{
	return this -> d_max;
}

double Graph::get_d_avg() const
{
	return this -> num_nodes != size_type(0) ? 2.0 * this -> num_edges / this -> num_nodes : degree_type(0);
}

const std::vector< std::vector<Graph::size_type> > & Graph::get_adjacency_list() const
{
	return this -> adjacency_list;
}

const bool Graph::get_is_transition_matrix_valid() const
{
	return this -> is_transition_matrix_valid;
}

const Eigen::SparseMatrix<double, 0, int64_t> & Graph::get_transition_matrix() const
{
	if(this -> is_transition_matrix_valid == false)
		throw;
	return this -> transition_matrix;
}

const Graph::degree_type Graph::get_d(Graph::size_type n) const
{
	return this -> adjacency_list[n].size();
}

const Graph::size_type Graph::get_max_degree_node() const
{
	return this -> max_degree_node;
}