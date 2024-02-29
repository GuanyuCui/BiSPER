#ifndef NODE_SAMPLER_H
#define NODE_SAMPLER_H

#include "graph.h"
#include <random>

namespace NodeSampler_Tools
{
	template <typename T, typename T_weight = double>
	// https://stackoverflow.com/questions/53632441/c-sampling-from-discrete-distribution-without-replacement.
	std::vector<T> weighted_sampling_without_replacement(const std::vector<T> & items, const std::vector<T_weight> & weights, Graph::size_type num)
	{
		if(items.size() != weights.size())
		{
			throw std::invalid_argument("Different number of items and weights.");
		}
		
		// Store index and key.
		std::vector<std::pair<size_t, double>> idx_keys(items.size());
		
		std::random_device rd;
		std::mt19937 gen(rd());
		// Generate random number in [0, 1].
		std::uniform_real_distribution<double> dis(0, std::nextafter(1, std::numeric_limits<double>::max()));
		// For each item, u_i = random(0, 1) and k_i = u^(1 / wi).
		for(Graph::size_type i = 0; i < idx_keys.size(); i++)
		{
			idx_keys[i] = std::make_pair(i, std::pow(dis(gen), 1.0 / weights[i]));
		}
		
		// Build a priority_queue such that index with larger key has higher priority.
		auto lambda_cmp = [](const std::pair<size_t, double> & a, const std::pair<size_t, double> & b) -> bool {return a.second < b.second;};
		std::priority_queue<std::pair<size_t, double>, std::vector<std::pair<size_t, double> >, decltype(lambda_cmp) > q(idx_keys.begin(), idx_keys.end(), lambda_cmp);
		
		// Get top num items.
		std::vector<T> res;
		for(Graph::size_type i = 0; i < num; i++)
		{
			res.push_back(items[q.top().first]);
			q.pop();
		}
		return res;
	}

	std::vector<Graph::size_type> read_node_ids(const char * file_name);
	bool save_node_ids(const char * file_name, const std::vector<Graph::size_type> & nodes, bool overwrite = true);
}

class NodeSampler
{
	public:
		// NodeSampler needs to bind a graph (and sort adjacency list if needed).
		NodeSampler(const Graph & G);
		~NodeSampler() = default;

		// Return (top num) nodes with minimum / maximum degree.
		std::vector<Graph::size_type> min_degree(Graph::size_type num);
		std::vector<Graph::size_type> max_degree(Graph::size_type num);

		// Return uniformly random sampled nodes (w/, w/o replacement).
		std::vector<Graph::size_type> uniform_with_replacement(Graph::size_type num);
		std::vector<Graph::size_type> uniform_without_replacement(Graph::size_type num);

		// Return degree uniformly random (Pr(node i is sampled) ~ d(i) = d(i) / 2m) sampled nodes (w/, w/o replacement).
		std::vector<Graph::size_type> degree_uniform_with_replacement(Graph::size_type num);
		std::vector<Graph::size_type> degree_uniform_without_replacement(Graph::size_type num);

	private:
		// Stores node_id with non-zero degree in degree descending order.
		std::vector<Graph::size_type> nodes;
		// Stores the degrees of the nodes (in corresponding order).
		std::vector<Graph::degree_type> degrees;

};

#endif