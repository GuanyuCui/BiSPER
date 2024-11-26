#include <iostream>
#include <iomanip>
#include <chrono>
#include <set>
#include <sstream>
#include <thread>
#if __has_include(<filesystem>)
	#include <filesystem>
	namespace fs = std::filesystem;
#elif __has_include(<experimental/filesystem>)
	#include <experimental/filesystem> 
	namespace fs = std::experimental::filesystem;
#else
	error "Missing the <filesystem> header."
#endif

#include "../include/graph.h"
#include "../include/node_sampler.h"
#include "../include/spalgorithms.h"
#include "../include/timer.h"

void print_duration(const Timer & timer, const std::string & unit = "ms")
{
	std::cout << "Time taken: " << timer.get_duration(unit) << " " + unit + "." << std::endl;
}

bool write_string_to_file(const char * file_name, const std::string & content, std::ios_base::openmode mode = std::ios_base::out | std::ios_base::app)
{
	// Open the file as a stream.
	std::ofstream file_stream(file_name, mode);
	// File not found.
	if(!file_stream.is_open())
	{
		file_stream.close();
		throw std::invalid_argument("File not found.");
		return false;
	}

	file_stream << content;

	// Close file.
	file_stream.close();	
	return true;
}

std::string doubletostring(double value, int precision)
{
	std::stringstream stream;
	stream << std::fixed << std::setprecision(precision) << value;
	return stream.str();
}

bool build_graph_from_dataset(Graph & G, const std::string & dataset_name)
{
	std::string txt_path = "./datasets/" + dataset_name + ".txt";
	std::string bin_path = "./datasets/" + dataset_name + ".bin";
	std::string compressed_sorted_bin_path = "./datasets/" + dataset_name + "-compressed_sorted.bin";
	// Reading graph.
	if(fs::exists(compressed_sorted_bin_path))
	{
		std::cout << "Reading graph (node_id compressed, sorted adjacency list) from binary file into memory..." << std::endl;
		G.from_bin(compressed_sorted_bin_path.c_str());
		std::cout << "Done." << std::endl;
	}
	else if(fs::exists(bin_path))
	{
		std::cout << "Reading graph from binary file into memory..." << std::endl;
		G.from_bin(bin_path.c_str());
		std::cout << "Done." << std::endl;

		std::cout << "Compressing node_id..." << std::endl;
		G.compress_adjacency_list();
		std::cout << "Done." << std::endl;

		std::cout << "Sorting adjacency list..." << std::endl;
		G.sort_adjacency_list();
		std::cout << "Done." << std::endl;
		
		std::cout << "Writing graph (sorted adjacency list) from memory into binary..." << std::endl;
		G.to_bin(compressed_sorted_bin_path.c_str());
		std::cout << "Done." << std::endl;
	}
	else if(fs::exists(txt_path))
	{
		std::cout << "Reading graph from text file into memory..." << std::endl;
		G.from_txt(txt_path.c_str());
		std::cout << "Done." << std::endl;

		std::cout << "Compressing node_id..." << std::endl;
		G.compress_adjacency_list();
		std::cout << "Done." << std::endl;

		std::cout << "Sorting adjacency list..." << std::endl;
		G.sort_adjacency_list();
		std::cout << "Done." << std::endl;

		std::cout << "Writing graph from memory into binary..." << std::endl;
		G.to_bin(compressed_sorted_bin_path.c_str());
		std::cout << "Done." << std::endl;
	}
	else
	{
		std::cout << "Graph dataset not found!" << std::endl;
		return false;
	}
	return true;
}

double read_lambda(const std::string & dataset_name)
{
	std::string lambda_file = "./datasets/" + dataset_name + ".lambda";
	// Open the file as a stream.
	std::ifstream file_stream(lambda_file.c_str());
	// File not found.
	if(!file_stream.is_open())
	{
		file_stream.close();
		throw std::invalid_argument("File not found. Run get_lambda.py with the binary graph files.");
		return 0.0;
	}

	double lambda = 0.0;

	file_stream >> lambda;

	// Close file.
	file_stream.close();	
	return lambda;
}

std::vector<Graph::size_type> get_query_nodes(const std::string & dataset_name, NodeSampler & sampler, const std::string & method, Graph::size_type num)
{
	if(method != "degree" && method != "uniform" && method == "max_degree" && method == "min_degree")
	{
		throw std::invalid_argument("Parameter 'method' should be one of 'degree', 'uniform', 'max_degree', 'min_degree'.");
	}

	std::string query_nodes_file = "./samples/" + dataset_name + "-" + method + ".query";
	std::vector<Graph::size_type> ret;
	if(fs::exists(query_nodes_file))
	{
		std::cout << "Reading queries from file..." << std::endl;
		ret = NodeSampler_Tools::read_node_ids(query_nodes_file.c_str());
		std::cout << "Done." << std::endl;
	}
	else
	{
		if(method == "degree")
		{
			ret = sampler.degree_uniform_without_replacement(num);
		}
		else if(method == "uniform")
		{
			ret = sampler.uniform_without_replacement(num);
		}
		else if(method == "max_degree")
		{
			ret = sampler.max_degree(num);
		}
		else // if(method == "min_degree")
		{
			ret = sampler.min_degree(num);
		}
		std::cout << "Saving query nodes file..." << std::endl;
		NodeSampler_Tools::save_node_ids(query_nodes_file.c_str(), ret);
		std::cout << "Done." << std::endl;
	}
	return ret;
}

std::vector<double> get_ground_truth(const std::string & dataset_name, Graph & G, const std::string & query_node_method, const std::vector<Graph::size_type> & query_nodes, Algorithms::length_type L_max)
{
	std::string ground_truth_file = "./samples/" + dataset_name + "-" + query_node_method + "-L_max=" + std::to_string(L_max) + ".ground_truth";
	std::vector<double> ret;
	if(fs::exists(ground_truth_file))
	{
		std::cout << "Reading ground-truth from file..." << std::endl;
		
		// Close syncing with stdio.
		std::ios::sync_with_stdio(false);
		// Open the file as a stream.
		std::ifstream file_stream(ground_truth_file);
		// File not found.
		if(!file_stream.is_open())
		{
			file_stream.close();
			std::ios::sync_with_stdio(true);
			return {};
		}

		double value = 0.0;
		while(file_stream >> value)
		{
			ret.push_back(value);
		}
			
		// Close file.
		file_stream.close();
		std::ios::sync_with_stdio(true);
		std::cout << "Done." << std::endl;
	}
	else
	{
		std::cout << "Generating ground-truth..." << std::endl;
		SinglePairAlgorithms alg(G);

		if(dataset_name == "Friendster" && G.get_is_transition_matrix_valid() == false)
		{
			std::cout << "Initializing transition matrix..." << std::endl;
			Timer timer;
			timer.clear();
			timer.start();
			G.init_transition_matrix();
			timer.stop();
			std::cout << "Done." << std::endl;
			print_duration(timer, "ms");
		}

		std::string algorithm_name = (dataset_name == "Friendster" ? "Power-Iteration-Eigen" : "Power-Iteration");
		alg.init_data_stucture(algorithm_name, L_max);

		double total_time = 0.0;
		Timer timer;
		for(size_t i = 0; i < query_nodes.size(); i += 2)
		{
			Graph::size_type s = query_nodes[i];
			Graph::size_type t = query_nodes[i + 1];

			// Ground-truth generation.
			std::cout << "Pair " << i / 2 + 1 << "..." << std::endl;
			// Run algorithm. 
			timer.clear();
			timer.start();
			double r = alg.run(algorithm_name, s, t, L_max, 0.0, 0.0);
			timer.stop();
			std::cout << "Done. " << std::endl;
			total_time += timer.get_duration("ms");
			print_duration(timer, "ms");
			ret.push_back(r);
		}
		double average_time = total_time / (query_nodes.size() / 2);
		std::cout << "All done. Average time: " << average_time << "ms." << std::endl;
		write_string_to_file("./samples/ground_truth_generation_time.txt", dataset_name + ":" + doubletostring(average_time, 16) + '\n');

		std::cout << "Saving ground-truth file..." << std::endl;
		// Close syncing with stdio.
		std::ios::sync_with_stdio(false);
		// Open the file as a stream.
		std::ofstream file_stream(ground_truth_file);
		// File not found.
		if(!file_stream.is_open())
		{
			file_stream.close();
			std::ios::sync_with_stdio(true);
			return ret;
		}

		for(double & value : ret)
		{
			file_stream << std::fixed << std::setprecision(16) << value << '\n';
		}
			
		// Close file.
		file_stream.close();
		std::ios::sync_with_stdio(true);
		std::cout << "Done." << std::endl;
	}
	return ret;
}

std::vector<double> get_ground_truth_parallel(const std::string & dataset_name, Graph & G, const std::string & query_node_method, const std::vector<Graph::size_type> & query_nodes, Algorithms::length_type L_max)
{
	std::string ground_truth_file = "./samples/" + dataset_name + "-" + query_node_method + "-L_max=" + std::to_string(L_max) + ".ground_truth";
	if(fs::exists(ground_truth_file))
	{
		std::cout << "Reading ground-truth from file..." << std::endl;
		
		// Close syncing with stdio.
		std::ios::sync_with_stdio(false);
		// Open the file as a stream.
		std::ifstream file_stream(ground_truth_file);
		// File not found.
		if(!file_stream.is_open())
		{
			file_stream.close();
			std::ios::sync_with_stdio(true);
			return {};
		}
		
		std::vector<double> ret;
		double value = 0.0;
		while(file_stream >> value)
		{
			ret.push_back(value);
		}
			
		// Close file.
		file_stream.close();
		std::ios::sync_with_stdio(true);
		std::cout << "Done." << std::endl;
		return ret;
	}
	else
	{
		std::cout << "Generating ground-truth..." << std::endl;
		// Number of query pairs.
		Graph::size_type num_query_pairs = query_nodes.size() / 2;
		// Number of threads.
		size_t num_threads = std::thread::hardware_concurrency();
		// The return values (ground-truth) of each pair, indexed by the index of the pair (ranged in [0, num_query_pairs - 1]).
		std::vector<double> ret(num_query_pairs, 0.0);
		// Total query time of each thread.
		std::vector<double> query_times(num_threads, 0.0);

		if(dataset_name == "Friendster" && G.get_is_transition_matrix_valid() == false)
		{
			std::cout << "Initializing transition matrix..." << std::endl;
			Timer timer;
			timer.clear();
			timer.start();
			G.init_transition_matrix();
			timer.stop();
			std::cout << "Done." << std::endl;
			print_duration(timer, "ms");
		}

		// The task of each thread (indexed by thread_index, ranged in [0, num_threads - 1]):
		// Calculate the ground-truth of each pair from start_index (included) to end_index (not included),
		// write the results to ret, and the total_time to query_times.
		auto task = [&](size_t thread_index, size_t start_index, size_t end_index)
		{
			SinglePairAlgorithms alg(G);
			std::string algorithm_name = (dataset_name == "Friendster" ? "Power-Iteration-Eigen" : "Power-Iteration");
			alg.init_data_stucture(algorithm_name, L_max);
			
			double total_time = 0.0;
			Timer timer;
			for(size_t index = start_index; index < end_index; index++)
			{
				Graph::size_type s = query_nodes[2 * index];
				Graph::size_type t = query_nodes[2 * index + 1];
				std::cout << "Pair " << index << "..." << std::endl;
				std::string algorithm_name = (dataset_name == "Friendster" ? "Power-Iteration-Eigen" : "Power-Iteration");
				// Run algorithm. 
				timer.clear();
				timer.start();
				double r = alg.run(algorithm_name, s, t, L_max, 0.0, 0.0);
				timer.stop();
				std::cout << "Done. " << std::endl;
				total_time += timer.get_duration("ms");
				print_duration(timer, "ms");
				ret[index] = r;
			}
			query_times[thread_index] = total_time;
		};

		std::vector<std::thread> threads;
		// Partition tasks.
		size_t index_step = std::ceil(double(num_query_pairs) / double(num_threads));
		for(size_t i = 0; i < num_threads; i++)
		{
			size_t start_index = index_step * i;
			size_t end_index = std::min(index_step * (i + 1), size_t(num_query_pairs));
			threads.emplace_back(task, i, start_index, end_index);
			if(end_index == num_query_pairs)
				break;
		}

		for(auto & thread : threads)
		{
			thread.join();
		}

		double average_time = std::accumulate(query_times.begin(), query_times.end(), 0.0) / num_query_pairs;
		std::cout << "All done. Average time:" << average_time << std::endl;
		write_string_to_file("./samples/ground_truth_generation_time.txt", dataset_name + ":" + doubletostring(average_time, 16) + '\n');

		std::cout << "Saving ground-truth file..." << std::endl;
		// Close syncing with stdio.
		std::ios::sync_with_stdio(false);
		// Open the file as a stream.
		std::ofstream file_stream(ground_truth_file);
		// File not found.
		if(!file_stream.is_open())
		{
			file_stream.close();
			std::ios::sync_with_stdio(true);
			return ret;
		}

		for(double & value : ret)
		{
			file_stream << std::fixed << std::setprecision(16) << value << '\n';
		}
			
		// Close file.
		file_stream.close();
		std::ios::sync_with_stdio(true);
		std::cout << "Done." << std::endl;
		return ret;
	}
}

// Save Pr_index
void save_Pr_index(const std::string & dataset_name, const std::vector<std::vector<double>> & data, 
					Algorithms::size_type num_landmarks, Algorithms::size_type num_pr_sample = 10000)
{
	std::string index_file = "./datasets/" + dataset_name + "-num_pr_sample=" + std::to_string(num_pr_sample) + "-num_landmarks=" + std::to_string(num_landmarks) + ".Prindex";
	std::ofstream outfile(index_file, std::ios::binary);
	
	if (!outfile) {
		std::cerr << "Cannot open file!" << std::endl;
		return;
	}

	// Write outer size.
	size_t outer_size = data.size();
	outfile.write(reinterpret_cast<const char*>(&outer_size), sizeof(outer_size));

	// Write inner vectors.
	for (const auto & inner_vec : data) {
		// Write inner vector size.
		size_t inner_size = inner_vec.size();
		outfile.write(reinterpret_cast<const char*>(&inner_size), sizeof(inner_size));

		// Write inner vector content.
		outfile.write(reinterpret_cast<const char*>(inner_vec.data()), inner_size * sizeof(double));
	}

	outfile.close();
}

// get Pr_index from file or build from scratch.
std::vector<std::vector<double>> get_Pr_index(const std::string & dataset_name, SinglePairAlgorithms & alg, 
												std::vector<Graph::size_type> & vl, std::vector<long long> & index_vl, 
												Algorithms::size_type num_landmarks, Algorithms::size_type num_pr_sample = 10000)
{
	std::string index_file = "./datasets/" + dataset_name + "-num_pr_sample=" + std::to_string(num_pr_sample) + "-num_landmarks=" + std::to_string(num_landmarks) + ".Prindex";
	std::vector<std::vector<double>> data;
	if(fs::exists(index_file))
	{
		std::cout << "Reading indices from file..." << std::endl;
		std::ifstream infile(index_file, std::ios::binary);
		
		if(!infile)
		{
			std::cerr << "Cannot open file!" << std::endl;
			return data;
		}

		// Read outer size.
		size_t outer_size;
		infile.read(reinterpret_cast<char*>(&outer_size), sizeof(outer_size));
		data.resize(outer_size);

		// Read in vector.
		for (size_t i = 0; i < outer_size; ++i)
		{
			// Read in vector size.
			size_t inner_size;
			infile.read(reinterpret_cast<char*>(&inner_size), sizeof(inner_size));
			data[i].resize(inner_size);

			// Read content.
			infile.read(reinterpret_cast<char*>(data[i].data()), inner_size * sizeof(double));
		}

		infile.close();
		std::cout << "Done." << std::endl;
		return data;
	}
	else
	{
		std::cout << "Precomputing indices..." << std::endl;
		data = alg.precompute_Pr(vl, index_vl, num_pr_sample);
		std::cout << "Done." << std::endl;
		std::cout << "Saving indices to file..." << std::endl;
		save_Pr_index(dataset_name, data, num_landmarks, num_pr_sample);
		std::cout << "Done." << std::endl;
		return data;
	}

}

template <typename KeyType>
double get_max_error(const std::map<KeyType, double> & ground_truth, const std::map<KeyType, double> & results)
{
	double max_error = 0.0;
	for(const auto & pair : results)
	{
		double error = std::abs(pair.second - ground_truth.at(pair.first));
		max_error = std::max(max_error, error);
	}
	return max_error;
}

void print_usage()
{
	std::cout << "---------------------------------" << std::endl;
	std::cout << "|\t\tSPER\t\t|" << std::endl;
	std::cout << "---------------------------------" << std::endl;

	std::cout << "Arguments:" << std::endl;
	std::cout << "\t--help: print this help message." << std::endl;
	std::cout << "\t--dataset: set the dataset to run algorithms. (default: Facebook)" << std::endl;
	std::cout << "\t--algorithm: set the competitor algorithm. (default: BiSPER)" << std::endl;
	std::cout << "\t--num_query: set the number of query pairs. (default: 100)" << std::endl;
	std::cout << "\t--L_max: set the maximum random walk length in the algorithm, auto if input \"auto\" (default: auto)" << std::endl;
	std::cout << "\t--eps: set the absolute error parameter \\epsilon in the algorithm. (default: 1e-2)" << std::endl;
	std::cout << "\t--p_f: set the failure probability p_f in the algorithm. (default: 1e-2)" << std::endl;
	std::cout << "\t--num_landmarks: set the number of landmarks. (default: 100)" << std::endl;
	std::cout << "\t--num_samples: set the number of sampled random walks. (default: 10000)" << std::endl;
	std::cout << "\t--r_max: set the threshold in landmark push algorithms. (default: 1e-4)" << std::endl;
}

int main(int argc, char * argv[])
{
	std::string dataset_name = "Facebook";
	std::string algorithm_name = "BiSPER";
	Graph::size_type num_query_pairs = size_t(100);
	bool auto_L_max = true;
	Algorithms::length_type L_max = 0;
	double eps = 1e-2;
	double p_f = 1e-2;
	Algorithms::size_type num_landmarks = 100;
	Algorithms::size_type num_samples = 10000;
	double r_max = 1e-4;

	// Iterate through all arguments.
	for(int i = 1; i < argc; i++)
	{
		std::string arg = argv[i];
		// Parse arguments.
		if(arg == "--help")
		{
			print_usage();
			return 0;
		}
		else if(arg == "--dataset")
		{
			if (i + 1 < argc)
			{
				dataset_name = argv[i + 1];
				// Skip the next argument as it has been read.
				i++;
			}
			else
			{
				std::cerr << "--dataset requires a value!" << std::endl;
				return -1;
			}
		}
		else if(arg == "--algorithm")
		{
			if (i + 1 < argc)
			{
				algorithm_name = argv[i + 1];
				i++;
			}
			else
			{
				std::cerr << "--algorithm requires a value!" << std::endl;
				return -1;
			}
		}
		else if(arg == "--num_query")
		{
			if (i + 1 < argc)
			{
				num_query_pairs = std::stoull(argv[i + 1]);
				i++;
			}
			else
			{
				std::cerr << "--num_query requires a value!" << std::endl;
				return -1;
			}
		}
		else if(arg == "--L_max")
		{
			if (i + 1 < argc)
			{
				try
				{
					L_max = std::stoull(argv[i + 1]);
					auto_L_max = false;
				}
				catch(...)
				{
					auto_L_max = true; 
				}
				i++;
			}
			else
			{
				std::cerr << "--L requires a value!" << std::endl;
				return -1;
			}
		}

		else if(arg == "--eps")
		{
			if (i + 1 < argc)
			{
				eps = std::stod(argv[i + 1]);
				i++;
			}
			else
			{
				std::cerr << "--eps requires a value!" << std::endl;
				return -1;
			}
		}
		else if(arg == "--p_f")
		{
			if (i + 1 < argc)
			{
				p_f = std::stod(argv[i + 1]);
				i++;
			}
			else
			{
				std::cerr << "--p_f requires a value!" << std::endl;
				return -1;
			}
		}
		else if(arg == "--num_landmarks")
		{
			if (i + 1 < argc)
			{
				num_landmarks = std::stoull(argv[i + 1]);
				i++;
			}
			else
			{
				std::cerr << "--num_landmarks requires a value!" << std::endl;
				return -1;
			}
		}
		else if(arg == "--num_samples")
		{
			if (i + 1 < argc)
			{
				num_samples = std::stoull(argv[i + 1]);
				i++;
			}
			else
			{
				std::cerr << "--num_samples requires a value!" << std::endl;
				return -1;
			}
		}
		else if(arg == "--r_max")
		{
			if (i + 1 < argc)
			{
				r_max = std::stod(argv[i + 1]);
				i++;
			}
			else
			{
				std::cerr << "--r_max requires a value!" << std::endl;
				return -1;
			}
		}
		else
		{
			std::cerr << "Unknown argument: " << arg << std::endl;
			print_usage();
			return -1;
		}
	}
	
	std::cout << "Dataset: " << dataset_name << std::endl;
	std::cout << "Algorithm: " << algorithm_name << std::endl;
	std::cout << "num_query: " << num_query_pairs << std::endl;
	std::cout << "L_max: " << (auto_L_max ? "auto" : std::to_string(L_max)) << std::endl;
	std::cout << "epsilon: " << eps << std::endl;
	std::cout << "p_f: " << p_f << std::endl;
	std::cout << "num_landmarks: " << num_landmarks << std::endl;
	std::cout << "num_samples: " << num_samples << std::endl;
	std::cout << "r_max: " << r_max << std::endl;
	std::cout << std::endl;

	// The graph.
	Graph G;
	// Timer.
	Timer timer;

	// Read graph.
	timer.clear();
	timer.start();
	build_graph_from_dataset(G, dataset_name);
	timer.stop();
	print_duration(timer, "ms");

	if(algorithm_name == "Power-Iteration-Eigen")
	{
		std::cout << "Initializing transition matrix..." << std::endl;
		timer.clear();
		timer.start();
		G.init_transition_matrix();
		timer.stop();
		std::cout << "Done." << std::endl;
		print_duration(timer, "ms");
	}

	double lambda = 0.90;
	if(auto_L_max)
	{
		// Get lambda.
		lambda = read_lambda(dataset_name);
	}
	
	// Print statistics.
	G.print_statistics();
	std::cout << std::endl;

	// Algorithm.
	SinglePairAlgorithms alg(G);

	// Experiment 0.
	if(num_query_pairs == 1)
	{
		Graph::size_type s = 1;
		Graph::size_type t = 1000;

		std::cout << "Querying pair (1, 1000)..."  << std::endl;

		// L_max.
		Graph::degree_type d_s = G.get_d(s);
		Graph::degree_type d_t = G.get_d(t);

		std::cout << "d(s) = " << d_s << " , d(t) = " << d_t << std::endl;
		if(auto_L_max)
		{
			L_max = std::max(1., std::ceil(std::log(2 * (1. / d_s + 1. / d_t) / (eps * (1. - lambda))) / std::log(1. / lambda)));
		}
		
		std::cout << "Initializing data structures..." << std::endl;
		timer.clear();
		timer.start();
		alg.init_data_stucture(algorithm_name, L_max);
		timer.stop();
		std::cout << "Done." << std::endl;
		print_duration(timer, "ms");
		
		// Run competitor algorithm. 
		std::cout << "Running algorithm for ground truth..." << std::endl;
		timer.clear();
		timer.start();
		double result = alg.run("Power-Iteration", s, t, L_max, eps, p_f);
		timer.stop();
		std::cout << "Done. Result = " << result << "." << std::endl;
		print_duration(timer, "ms");
		std::cout << std::endl;

		// File to write results.
		write_string_to_file("results/Experiment-0.out", 
			"Dataset:" + dataset_name + 
			"\tL_max:" + (auto_L_max ? "auto" : std::to_string(L_max)) +
			"\tER:" + doubletostring(result, 16) + "\n"
		);

		return 0;
	}

	// Node sampler.
	NodeSampler sampler(G);

	// Get landmark nodes.
	std::vector<Graph::size_type> vl;
	std::vector<long long> index_vl(G.get_num_nodes(), 0);
	if(algorithm_name == "RW-vl" || algorithm_name == "Push-vl" || algorithm_name == "Bipush-vl")
	{
		vl = sampler.max_degree(num_landmarks);
		auto generate_index_vl = [&G](std::vector<Graph::size_type> & vl, std::vector<long long> & index_vl) 
		{
			// std::vector<int> index_vl(n, 0);
			for(size_t i = 1; i <= vl.size(); i++)
			{
				index_vl[vl[i - 1]] = i;
			}
			long long u_ind = 0;
			for(size_t i = 0; i < G.get_num_nodes(); i++)
			{
				if(index_vl[i] == 0)
				{
					index_vl[i] = u_ind;
					u_ind -= 1;
				}
			}
		};
		generate_index_vl(vl, index_vl);
	}

	// Get query nodes.
	std::vector<Graph::size_type> query_nodes = get_query_nodes(dataset_name, sampler, "uniform", 2 * num_query_pairs);
	std::cout << std::endl;
	// Get ground-truth values.
	Algorithms::length_type L_max_ground_truth = L_max;
	// Experiment II/III.
	if(auto_L_max)
	{
		L_max_ground_truth = std::max(1., std::ceil(std::log(4 / (1e-17 * (1. - lambda))) / std::log(1. / lambda)));
	}
	std::vector<double> ground_truths = get_ground_truth_parallel(dataset_name, G, "uniform", query_nodes, L_max_ground_truth);
	std::cout << std::endl;
	
	// Vectors that store error and query time of each query node.
	std::vector<double> errors;
	std::vector<double> query_times;

	std::cout << "Initializing data structures..." << std::endl;
	timer.clear();
	timer.start();
	if(auto_L_max)
	{
		Algorithms::length_type L_max_eps = std::max(1., std::ceil(std::log(4 / (eps * (1. - lambda))) / std::log(1. / lambda)));
		if(algorithm_name == "AbWalk" || algorithm_name == "Push" || algorithm_name == "Bipush")
		{
			alg.init_data_stucture_landmark(algorithm_name, G.get_max_degree_node());
		}
		else if(algorithm_name == "RW-vl"|| algorithm_name == "Push-vl" || algorithm_name == "Bipush-vl")
		{
			std::vector<std::vector<double>> Pr = get_Pr_index(dataset_name, alg, vl, index_vl, num_landmarks, 100000);
			alg.init_data_stucture_vl(algorithm_name, vl, index_vl, Pr, num_samples);
		}
		else
		{
			alg.init_data_stucture(algorithm_name, L_max_eps);
		}
	}
	else
	{
		if(algorithm_name == "AbWalk" || algorithm_name == "Push" || algorithm_name == "Bipush")
		{
			alg.init_data_stucture_landmark(algorithm_name, G.get_max_degree_node());
		}
		else if(algorithm_name == "RW-vl"|| algorithm_name == "Push-vl" || algorithm_name == "Bipush-vl")
		{
			std::vector<std::vector<double>> Pr = get_Pr_index(dataset_name, alg, vl, index_vl, num_landmarks, 100000);
			alg.init_data_stucture_vl(algorithm_name, vl, index_vl, Pr, num_samples);
		}
		else
		{
			alg.init_data_stucture(algorithm_name, L_max);
		}
	}
	timer.stop();
	std::cout << "Done." << std::endl;
	print_duration(timer, "ms");

	for(size_t i = 0; i < num_query_pairs; i++)
	{
		Graph::size_type s = query_nodes[2 * i];
		Graph::size_type t = query_nodes[2 * i + 1];

		std::cout << "Querying pair No." << i << ": (s, t) = " << "(" << s << ", " << t << ")" << std::endl;
		// L_max.
		Graph::degree_type d_s = G.get_d(s);
		Graph::degree_type d_t = G.get_d(t);
		if(auto_L_max)
		{
			L_max = std::max(1., std::ceil(std::log(2 * (1. / d_s + 1. / d_t) / (eps * (1. - lambda))) / std::log(1. / lambda)));
		}

		double ground_truth = ground_truths[i];
		
		// Run competitor algorithm. 
		std::cout << "Running competitor algorithm..." << std::endl;
		double result = 0.0;
		timer.clear();
		timer.start();
		if(algorithm_name == "AbWalk" || algorithm_name == "Push" || algorithm_name == "Bipush"
			|| algorithm_name == "RW-vl" || algorithm_name == "Push-vl" || algorithm_name == "Bipush-vl")
		{
			result = alg.run_landmark(algorithm_name, s, t, num_landmarks, num_samples, r_max);
		}
		else
		{
			result = alg.run(algorithm_name, s, t, L_max, eps, p_f);
		}
		timer.stop();
		std::cout << "Competitor algorithm done. Result = " << result << "." << std::endl;
		print_duration(timer, "ms");
		std::cout << std::endl;

		// Calculate error.
		double error = std::abs(ground_truth - result);
		std::cout << "Error = " << error << std::endl << std::endl;
		errors.push_back(error);

		query_times.push_back(timer.get_duration("ms"));
	}

	// File to write results.
	std::string output_file = "results/" + dataset_name + ".out";
	write_string_to_file(output_file.c_str(), 
		"Dataset:" + dataset_name + 
		"\tAlgorithm:" + algorithm_name + 
		"\tnum_query:" + std::to_string(num_query_pairs) + 
		"\tL_max:" + (auto_L_max ? "auto" : std::to_string(L_max)) + 
		"\teps:" + doubletostring(eps, 16) + 
		"\tp_f:" + doubletostring(p_f, 16) +
		"\tnum_samples:" + std::to_string(num_samples) + 
		"\tr_max:" + doubletostring(r_max, 16)
	);

	double max_error = *std::max_element(errors.begin(), errors.end());
	write_string_to_file(output_file.c_str(), "\tmax.error:" + doubletostring(max_error, 16));

	double avg_error = std::accumulate(errors.begin(), errors.end(), 0.0) / errors.size();
	write_string_to_file(output_file.c_str(), "\tavg.error:" + doubletostring(avg_error, 16));

	double avg_query_time = std::accumulate(query_times.begin(), query_times.end(), 0.0) / query_times.size();
	write_string_to_file(output_file.c_str(), "\tavg.time(ms):" + doubletostring(avg_query_time, 16) + "\n");

	return 0;
}