#include "../include/node_sampler.h"

std::vector<Graph::size_type> NodeSampler_Tools::read_node_ids(const char * file_name)
{
	// Close syncing with stdio.
	std::ios::sync_with_stdio(false);
	// Open the file as a stream.
	std::ifstream file_stream(file_name);
	// File not found.
	if(!file_stream.is_open())
	{
		file_stream.close();
		std::ios::sync_with_stdio(true);
		return {};
	}

	std::vector<Graph::size_type> res;
	Graph::size_type node_id = 0;

	while(file_stream >> node_id)
	{
		res.push_back(node_id);
	}
		
	// Close file.
	file_stream.close();
	std::ios::sync_with_stdio(true);
	return res;
}

bool NodeSampler_Tools::save_node_ids(const char * file_name, const std::vector<Graph::size_type> & nodes, bool overwrite)
{
	// Close syncing with stdio.
	std::ios::sync_with_stdio(false);
	// Open the file as a stream.
	std::ofstream file_stream(file_name, overwrite ? std::ios_base::out : std::ios_base::app);
	// File not found.
	if(!file_stream.is_open())
	{
		file_stream.close();
		std::ios::sync_with_stdio(true);
		return false;
	}

	for(auto & node_id : nodes)
	{
		file_stream << node_id << '\n';
	}
		
	// Close file.
	file_stream.close();
	std::ios::sync_with_stdio(true);
	return true;
}

NodeSampler::NodeSampler(const Graph & G)
{
	this -> nodes.reserve(G.get_num_nodes());
	const auto & adj_list = G.get_adjacency_list();

	for(Graph::size_type i = 0; i < adj_list.size(); i++)
	{
		if(adj_list[i].size() > 0)
			this -> nodes.push_back(i);
	}

	std::sort(this -> nodes.begin(), this -> nodes.end(), 
			[&](const Graph::size_type & v1, const Graph::size_type & v2) -> bool {return adj_list[v1].size() > adj_list[v2].size();});

	this -> degrees.resize(G.get_num_nodes());
	std::transform(this -> nodes.begin(), this -> nodes.end(), this -> degrees.begin(), 
		[&](Graph::size_type & idx) -> Graph::degree_type {return adj_list[idx].size();});
}

// Return (top num) nodes with minimum / maximum degree.
std::vector<Graph::size_type> NodeSampler::min_degree(Graph::size_type num)
{
	if(num > this -> nodes.size())
	{
		throw std::out_of_range("num is larger than num_nodes.");
	}

	return std::vector<Graph::size_type>(this -> nodes.end() - num, this -> nodes.end());
}

std::vector<Graph::size_type> NodeSampler::max_degree(Graph::size_type num)
{
	if(num > this -> nodes.size())
	{
		throw std::out_of_range("num is larger than num_nodes.");
	}
	return std::vector<Graph::size_type>(this -> nodes.begin(), this -> nodes.begin() + num);
}

// Return uniformly random sampled nodes (w/, w/o replacement).
std::vector<Graph::size_type> NodeSampler::uniform_with_replacement(Graph::size_type num)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<Graph::size_type> d(0, this -> nodes.size() - 1);
	
	std::vector<Graph::size_type> res;
	for(Graph::size_type i = 0; i < num; i++)
	{
		res.push_back(this -> nodes[d(gen)]);
	}
	return res;
}


std::vector<Graph::size_type> NodeSampler::uniform_without_replacement(Graph::size_type num)
{
	std::vector<Graph::size_type> tmp_nodes(this -> nodes);
	std::shuffle(tmp_nodes.begin(), tmp_nodes.end(), std::mt19937_64{std::random_device{}()});
	return std::vector<Graph::size_type>(tmp_nodes.begin(), tmp_nodes.begin() + num);
}

// Return degree uniformly random (Pr(node i is sampled) ~ d(i) = d(i) / 2m) sampled nodes (w/, w/o replacement).
std::vector<Graph::size_type> NodeSampler::degree_uniform_with_replacement(Graph::size_type num)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::discrete_distribution<Graph::size_type> d(this -> degrees.begin(), this -> degrees.end());
	
	std::vector<Graph::size_type> res;
	for(Graph::size_type i = 0; i < num; i++)
	{
		res.push_back(this -> nodes[d(gen)]);
	}
	return res;
}

std::vector<Graph::size_type> NodeSampler::degree_uniform_without_replacement(Graph::size_type num)
{
	return NodeSampler_Tools::weighted_sampling_without_replacement<Graph::size_type, Graph::degree_type>(this -> nodes, this -> degrees, num);
}