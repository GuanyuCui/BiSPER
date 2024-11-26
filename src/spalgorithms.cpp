#include <array>
#include <map>
#include <stdexcept>
#include <random>
#include <utility>
#include <algorithm>

#include "../include/spalgorithms.h"

// -------------------- SinglePairAlgorithms --------------------
// -------------------- Public --------------------
SinglePairAlgorithms::SinglePairAlgorithms(const Graph & G) : Algorithms(G) {}

void SinglePairAlgorithms::init_data_stucture(const std::string & algorithm_name, Algorithms::length_type L_max)
{
	if(algorithm_name == "Power-Iteration")
	{
		this -> mat_s_n_times_2 = std::vector<std::vector<double>>(G.get_adjacency_list().size(), std::vector<double>(2, 0.0));
		this -> mat_t_n_times_2 = std::vector<std::vector<double>>(G.get_adjacency_list().size(), std::vector<double>(2, 0.0));
	}
	else if(algorithm_name == "AMC")
	{
		this -> vec_1_n = std::vector<double>(G.get_adjacency_list().size(), 0.0);
		this -> vec_2_n = std::vector<double>(G.get_adjacency_list().size(), 0.0);
	}
	else if(algorithm_name == "GEER")
	{
		this -> mat_s_2_times_n = std::vector<std::vector<double>>(2, std::vector<double>(G.get_adjacency_list().size(), 0.0));
		this -> mat_t_2_times_n = std::vector<std::vector<double>>(2, std::vector<double>(G.get_adjacency_list().size(), 0.0));
	}
	else if(algorithm_name == "BiSPER")
	{
		this -> BITs_s = std::vector<SinglePairAlgorithms_Tools::BIT<double>>(G.get_adjacency_list().size(), SinglePairAlgorithms_Tools::BIT<double>(L_max + 2));
		this -> BITs_t = std::vector<SinglePairAlgorithms_Tools::BIT<double>>(G.get_adjacency_list().size(), SinglePairAlgorithms_Tools::BIT<double>(L_max + 2));

		this -> BITs_visited_s = std::vector<bool>(G.get_adjacency_list().size(), false);
		this -> BITs_visited_t = std::vector<bool>(G.get_adjacency_list().size(), false);

		// In case of degeneration.
		this -> init_data_stucture("Power-Iteration", L_max);
	}
}

void SinglePairAlgorithms::init_data_stucture_landmark(const std::string & algorithm_name, Graph::size_type landmark)
{
	if(algorithm_name == "AbWalk")
	{
		// Do nothing.
	}
	else if(algorithm_name == "Bipush")
	{
		this -> vec_bool_n = std::vector<bool>(G.get_adjacency_list().size(), false);

		this -> vec_s_1_n = std::vector<double>(G.get_adjacency_list().size(), 0.0);
		this -> vec_t_1_n = std::vector<double>(G.get_adjacency_list().size(), 0.0);
		this -> vec_s_2_n = std::vector<double>(G.get_adjacency_list().size(), 0.0);
		this -> vec_t_2_n = std::vector<double>(G.get_adjacency_list().size(), 0.0);
	}
	else if(algorithm_name == "Push")
	{
		this -> vec_1_n = std::vector<double>(G.get_adjacency_list().size(), 0.0);
		this -> vec_bool_n = std::vector<bool>(G.get_adjacency_list().size(), false);

		this -> vec_s_1_n = std::vector<double>(G.get_adjacency_list().size(), 0.0);
		this -> vec_t_1_n = std::vector<double>(G.get_adjacency_list().size(), 0.0);
	}
}

void SinglePairAlgorithms::init_data_stucture_vl(const std::string & algorithm_name, std::vector<Graph::size_type> & vl, std::vector<long long> & index_vl, std::vector<std::vector<double>> & Pr, Algorithms::size_type num_samples)
{
	if(algorithm_name == "RW-vl")
	{
		this -> vl = vl;
		this -> index_vl = index_vl;
		//std::cout << "Precomouting SC_Inverse..." << std::endl;
		this -> SCInverse = this -> precompute_SchurComplementInverse(this -> vl, this -> index_vl, Pr);
		//std::cout << "Done!" << std::endl;
		//std::cout << "Precomouting P_SCInverse..." << std::endl;
		this -> Index = this -> combine_P_SCInverse(Pr, this -> SCInverse, this -> index_vl);
		//std::cout << "Done!" << std::endl;
	}
	else if(algorithm_name == "Push-vl" || algorithm_name == "Bipush-vl")
	{
		this -> vec_bool_n = std::vector<bool>(G.get_adjacency_list().size(), false);

		this -> vec_s_1_n = std::vector<double>(G.get_adjacency_list().size(), 0.0);
		this -> vec_t_1_n = std::vector<double>(G.get_adjacency_list().size(), 0.0);
		this -> vec_s_2_n = std::vector<double>(G.get_adjacency_list().size(), 0.0);
		this -> vec_t_2_n = std::vector<double>(G.get_adjacency_list().size(), 0.0);

		this -> vl = vl;
		this -> index_vl = index_vl;
		//std::cout << "Precomouting SC_Inverse..." << std::endl;
		this -> SCInverse = this -> precompute_SchurComplementInverse(this -> vl, this -> index_vl, Pr);
		//std::cout << "Done!" << std::endl;
		//std::cout << "Precomouting P_SCInverse..." << std::endl;
		this -> Index = this -> combine_P_SCInverse(Pr, this -> SCInverse, this -> index_vl);
		//std::cout << "Done!" << std::endl;
	}
}

void SinglePairAlgorithms::reset_zeros(const std::string & algorithm_name)
{
	if(algorithm_name == "Power-Iteration")
	{
		for(auto & vec : mat_s_n_times_2)
		{
			std::fill(vec.begin(), vec.end(), 0.0);
		}
		for(auto & vec : mat_t_n_times_2)
		{
			std::fill(vec.begin(), vec.end(), 0.0);
		}
	}
	else if(algorithm_name == "AMC")
	{
		std::fill(vec_1_n.begin(), vec_1_n.end(), 0.0);
		std::fill(vec_2_n.begin(), vec_2_n.end(), 0.0);
	}
	else if(algorithm_name == "GEER")
	{
		for(auto & vec : mat_s_2_times_n)
		{
			std::fill(vec.begin(), vec.end(), 0.0);
		}
		for(auto & vec : mat_t_2_times_n)
		{
			std::fill(vec.begin(), vec.end(), 0.0);
		}
	}
	else if(algorithm_name == "AbWalk" || algorithm_name == "RW-vl")
	{
		// Do nothing.
	}
	else if(algorithm_name == "Bipush")
	{
		std::fill(vec_bool_n.begin(), vec_bool_n.end(), false);
		std::fill(vec_s_1_n.begin(), vec_s_1_n.end(), 0.0);
		std::fill(vec_t_1_n.begin(), vec_t_1_n.end(), 0.0);
		std::fill(vec_s_2_n.begin(), vec_s_2_n.end(), 0.0);
		std::fill(vec_t_2_n.begin(), vec_t_2_n.end(), 0.0);
	}
	else if(algorithm_name == "Push")
	{
		std::fill(vec_1_n.begin(), vec_1_n.end(), 0.0);
		std::fill(vec_bool_n.begin(), vec_bool_n.end(), false);
		std::fill(vec_s_1_n.begin(), vec_s_1_n.end(), 0.0);
		std::fill(vec_t_1_n.begin(), vec_t_1_n.end(), 0.0);
	}
	else if(algorithm_name == "Push-vl" || algorithm_name == "Bipush-vl")
	{
		std::fill(vec_bool_n.begin(), vec_bool_n.end(), false);
		std::fill(vec_s_1_n.begin(), vec_s_1_n.end(), 0.0);
		std::fill(vec_t_1_n.begin(), vec_t_1_n.end(), 0.0);
		std::fill(vec_s_2_n.begin(), vec_s_2_n.end(), 0.0);
		std::fill(vec_t_2_n.begin(), vec_t_2_n.end(), 0.0);
	}
	else if(algorithm_name == "BiSPER")
	{
		Graph::size_type sz = G.get_adjacency_list().size();
		for(Graph::size_type i = 0; i < sz; i++)
		{
			if(BITs_visited_s[i])
				BITs_s[i].reset_zeros();
		}
		for(Graph::size_type i = 0; i < sz; i++)
		{
			if(BITs_visited_t[i])
				BITs_t[i].reset_zeros();
		}
		std::fill(BITs_visited_s.begin(), BITs_visited_s.end(), false);
		std::fill(BITs_visited_t.begin(), BITs_visited_t.end(), false);
	}
	else
	{
		std::cout << "Algorithm not found!" << std::endl;
		throw;
	}
}

// Interface for running single-pair effective resistance algorithms.
double SinglePairAlgorithms::run(const std::string & algorithm_name, Graph::size_type s, Graph::size_type t, Algorithms::length_type L_max, double eps, double p_f)
{
	if(algorithm_name == "Power-Iteration")
	{
		double R = 0.0;
		R += this -> Power_Iteration_ER(s, t, L_max);
		return R;
	}
	else if(algorithm_name == "Power-Iteration-Eigen")
	{
		double R = 0.0;
		R += this -> Power_Iteration_ER_Eigen(s, t, L_max);
		return R;
	}
	else if(algorithm_name == "AMC")
	{
		double R = 0.0;
		R = this -> AMC(s, t, L_max, eps / 2.0, p_f);
		return R;
	}
	else if(algorithm_name == "GEER")
	{
		double R = 0.0;
		R = this -> GEER(s, t, L_max, eps / 2.0, p_f);
		return R;
	}
	else if(algorithm_name == "BiSPER")
	{
		double R = 0.0;
		R = this -> BiSPER(s, t, L_max, eps / 2.0, p_f);
		return R;
	}
	throw std::invalid_argument("Algorithm not found!");
	return -1;
}

double SinglePairAlgorithms::run_landmark(const std::string & algorithm_name, Graph::size_type s, Graph::size_type t, Algorithms::size_type num_landmarks, Algorithms::size_type num_samples, double r_max)
{
	if(algorithm_name == "AbWalk")
	{
		double R = 0.0;
		R = this -> AbWalk(s, t, G.get_max_degree_node(), num_samples);
		return R;
	}
	else if(algorithm_name == "Bipush")
	{
		double R = 0.0;
		R = this -> Bipush(s, t, G.get_max_degree_node(), num_samples, r_max);
		return R;
	}
	else if(algorithm_name == "Push")
	{
		double R = 0.0;
		R = this -> Push(s, t, G.get_max_degree_node(), r_max);
		return R;
	}
	else if(algorithm_name == "RW-vl")
	{
		double R = 0.0;
		R = this -> RW_vl(s, t, num_landmarks, num_samples);
		return R;
	}
	else if(algorithm_name == "Push-vl")
	{
		double R = 0.0;
		R = this -> Push_vl(s, t, num_landmarks, r_max, num_samples);
		return R;
	}
	else if(algorithm_name == "Bipush-vl")
	{
		double R = 0.0;
		R = this -> Bipush_vl(s, t, num_landmarks, r_max, num_samples);
		return R;
	}
	throw std::invalid_argument("Algorithm not found!");
	return -1;
}

// ------------------------------ Private ------------------------------
// -------------------- Single-pair transition probability tool algorithms --------------------
// Ground-truth generation.
double SinglePairAlgorithms::Power_Iteration_ER(Graph::size_type s, Graph::size_type t, Algorithms::length_type L_max)
{
	// Holding node u such that q_s^{(k)}(u) > 0 and q_s^{(k + 1)}(u) > 0 node_id.
	std::array<std::queue<Graph::size_type>, 2> active_nodes_s;
	// Holding q_s^{(k)}(u) > 0 and q_s^{(k + 1)}(u) > 0 values.
	std::vector<std::vector<double>> & r_s = this -> mat_s_n_times_2;

	// Holding node u such that q_t^{(k)}(u) > 0 and q_t^{(k + 1)}(u) > 0 node_id.
	std::array<std::queue<Graph::size_type>, 2> active_nodes_t;
	// Holding q_t^{(k)}(u) > 0 and q_t^{(k + 1)}(u) > 0 values.
	std::vector<std::vector<double>> & r_t = this -> mat_t_n_times_2;

	this -> reset_zeros("Power-Iteration");

	Graph::degree_type d_s = G.get_d(s);
	Graph::degree_type d_t = G.get_d(t);

	Graph::size_type m = G.get_num_edges();
	const std::vector<std::vector<Graph::size_type>> & adjacency_list = G.get_adjacency_list();
	Graph::size_type adj_size = adjacency_list.size();

	// Node s.
	r_s[s][0] = 1.0;
	active_nodes_s[0].push(s);
	Graph::size_type push_cost_s = d_s;

	// Node t.
	r_t[t][0] = 1.0;
	active_nodes_t[0].push(t);
	Graph::size_type push_cost_t = d_t;

	// \sum_{\ell = 0}^{L_max} q_s^{(\ell)}(s) / d(s) - q_s^{(\ell)}(t) / d(t)
	double sum_qss_div_ds_minus_qst_div_dt = 0.0;
	// \sum_{\ell = 0}^{L_max} q_t^{(\ell)}(t) / d(t) - q_t^{(\ell)}(s) / d(s)
	double sum_qtt_div_dt_minus_qts_div_ds = 0.0;


	// Push l-th layer.
	for(Algorithms::length_type l = 0; l <= L_max; l++)
	{
		if(push_cost_s < m / 2)
		{
			push_cost_s = 0; 
			while(!active_nodes_s[l % 2].empty())
			{
				// Get a node with non-zero q.
				Graph::size_type u = active_nodes_s[l % 2].front();
				active_nodes_s[l % 2].pop();
				double r_u = r_s[u][l % 2];
				Graph::degree_type d_u = G.get_d(u);

				// Push (u, v).
				sum_qss_div_ds_minus_qst_div_dt += ((u == s) ? r_u / d_s : 0.0) - ((u == t) ? r_u / d_t : 0.0);
				for(const Graph::size_type & v : adjacency_list[u])
				{
					bool is_old_residue_small = (r_s[v][(l + 1) % 2] <= 0.0);
					r_s[v][(l + 1) % 2] += r_u / d_u;
					// bool is_new_residue_large = (r_s[v][(l + 1) % 2] > 0.0);

					// Node v is active.
					if(is_old_residue_small /*&& is_new_residue_large*/)
					{
						active_nodes_s[(l + 1) % 2].push(v);
						push_cost_s += G.get_d(v);
					}
				}
				r_s[u][l % 2] = 0.0;
			}
		}
		else
		{
			for(Graph::size_type u = 0; u < adj_size; u++)
			{
				double r_u = r_s[u][l % 2];
				Graph::degree_type d_u = G.get_d(u);
				// Push (u, v).
				sum_qss_div_ds_minus_qst_div_dt += ((u == s) ? r_u / d_s : 0.0) - ((u == t) ? r_u / d_t : 0.0);
				for(const Graph::size_type & v : adjacency_list[u])
				{
					r_s[v][(l + 1) % 2] += r_u / d_u;
				}
				r_s[u][l % 2] = 0.0;
			}
		}
	}

	for(Algorithms::length_type l = 0; l <= L_max; l++)
	{
		if(push_cost_t < m / 2)
		{
			push_cost_t = 0;
			while(!active_nodes_t[l % 2].empty())
			{
				// Get a node with non-zero q.
				Graph::size_type u = active_nodes_t[l % 2].front();
				active_nodes_t[l % 2].pop();
				double r_u = r_t[u][l % 2];
				Graph::degree_type d_u = G.get_d(u);

				// Push (u, v).
				sum_qtt_div_dt_minus_qts_div_ds += ((u == t) ? r_u / d_t : 0.0) - ((u == s) ? r_u / d_s : 0.0);
				for(const Graph::size_type & v : adjacency_list[u])
				{
					bool is_old_residue_small = (r_t[v][(l + 1) % 2] <= 0.0);
					r_t[v][(l + 1) % 2] += r_u / d_u;
					// bool is_new_residue_large = (r_t[v][(l + 1) % 2] > 0.0);

					// Node v is active.
					if(is_old_residue_small /*&& is_new_residue_large*/)
					{
						active_nodes_t[(l + 1) % 2].push(v);
						push_cost_t += G.get_d(v);
					}
				}
				r_t[u][l % 2] = 0.0;
			}
		}
		else
		{
			for(Graph::size_type u = 0; u < adj_size; u++)
			{
				double r_u = r_t[u][l % 2];
				Graph::degree_type d_u = G.get_d(u);
				// Push (u, v).
				sum_qtt_div_dt_minus_qts_div_ds += ((u == t) ? r_u / d_t : 0.0) - ((u == s) ? r_u / d_s : 0.0);
				for(const Graph::size_type & v : adjacency_list[u])
				{
					r_t[v][(l + 1) % 2] += r_u / d_u;
				}
				r_t[u][l % 2] = 0.0;
			}
		}
	}
	return sum_qss_div_ds_minus_qst_div_dt + sum_qtt_div_dt_minus_qts_div_ds;
}

double SinglePairAlgorithms::Power_Iteration_ER_Eigen(Graph::size_type s, Graph::size_type t, Algorithms::length_type L_max)
{
	double R = 0.0;

	Graph::degree_type d_s = G.get_d(s);
	Graph::degree_type d_t = G.get_d(t);

	Eigen::VectorXd v_s = Eigen::VectorXd::Zero(G.get_adjacency_list().size());
	v_s[s] = 1.0;

	Eigen::VectorXd v_t = Eigen::VectorXd::Zero(G.get_adjacency_list().size());
	v_t[t] = 1.0;

	R += (v_s[s] / d_s  - v_s[t] / d_t + v_t[t] / d_t - v_t[s] / d_s);
	for(Algorithms::length_type l = 1; l <= L_max; l++)
	{
		v_s = G.get_transition_matrix() * v_s;
		v_t = G.get_transition_matrix() * v_t;
		R += (v_s[s] / d_s  - v_s[t] / d_t + v_t[t] / d_t - v_t[s] / d_s);
	}
	return R;
}

// From SIGMOD'23 paper "Efficient Estimation of Pairwise Effective Resistance".
double SinglePairAlgorithms::AdaptiveMonteCarlo(Graph::size_type s, Graph::size_type t, Algorithms::length_type len_walk, double eps, double p_f, std::vector<double> & svec, std::vector<double> & tvec, double psi_s_prime, double psi_t_prime, Algorithms::size_type tau)
{
	Graph::degree_type ds = G.get_d(s);
	Graph::degree_type dt = G.get_d(t);
	double epsilon = eps;

	Algorithms::size_type eta_star_all = std::ceil(8.0 * std::log(tau / p_f) * (len_walk + 1) * (len_walk + 1) * std::pow(psi_s_prime / ds + psi_t_prime / dt, 2.0) / epsilon / epsilon);
	Algorithms::size_type eta_all = std::ceil(std::max(eta_star_all * 1.0 / pow(2, tau - 1), 1.0));

	double mu_all = 0.0;

	// Random device as seeding source.
	std::random_device rd;
	// Seeding.
	std::mt19937_64 gen(rd());

	const std::vector<std::vector<Graph::size_type>> & adjacency_list = G.get_adjacency_list();

	for(Algorithms::size_type i = 0; i < tau; i++)
	{
		mu_all = 0.0;
		double sigma_all = 0.0;
		for(Algorithms::size_type k = 0; k < eta_all; k++)
		{
			Graph::size_type cur = s;
			double xss = 0.0;
			double xst = 0.0;
			for(Algorithms::length_type j = 0; j <= len_walk; j++)
			{
				Graph::degree_type deg = G.get_d(cur);
				std::uniform_int_distribution<Graph::degree_type> dis(0, deg - 1);
				Graph::size_type l = dis(gen);
				cur = adjacency_list[cur][l];
				xss += svec[cur] * 1.0;
				xst += tvec[cur] * 1.0;
			}

			cur = t;
			double xtt = 0.0;
			double xts = 0.0;
			for(Algorithms::length_type j = 0; j <= len_walk; j++)
			{
				Graph::degree_type deg = G.get_d(cur);
				std::uniform_int_distribution<Graph::degree_type> dis(0, deg - 1);
				Graph::size_type l = dis(gen);
				cur = adjacency_list[cur][l];
				xtt += tvec[cur] * 1.0;
				xts += svec[cur] * 1.0;
			}

			double mu_update = xss / ds + xtt / dt - xst / dt - xts / ds;
			mu_all += mu_update;
			sigma_all = sigma_all + std::pow(mu_update, 2);
		}

		mu_all = mu_all / eta_all;
		sigma_all = sigma_all / eta_all - std::pow(mu_all, 2);
		// Bug!
		// eta_all = eta_all * 2;

		double deltaprime = p_f / tau;
		double psi_all = 2 * (len_walk + 1) * (psi_s_prime / ds + psi_t_prime / dt);
		double eps_f = std::sqrt(2.0 * sigma_all * std::log(3.0 / deltaprime) / eta_all) + 3.0 * psi_all * std::log(3.0 / deltaprime) / eta_all;

		if(eps_f < epsilon)
		{
			break;
		}
		// Bug fix.
		eta_all = eta_all * 2;
	}

	double rf = mu_all;
	return rf;
}

double SinglePairAlgorithms::AMC(Graph::size_type s, Graph::size_type t, Algorithms::length_type L_max, double eps, double p_f, Algorithms::size_type tau)
{
	std::vector<double> & svec = this -> vec_1_n;
	std::vector<double> & tvec = this -> vec_2_n;

	this -> reset_zeros("AMC");

	svec[s] = 1.0;
	tvec[t] = 1.0;
	double er = svec[s] / G.get_d(s) + tvec[t] / G.get_d(t) - 2.0 * svec[t] / G.get_d(t);
	er += AdaptiveMonteCarlo(s, t, L_max, eps, p_f, svec, tvec, 1.0, 1.0);
	svec[s] = 0.0;
	tvec[t] = 0.0;
	return er;
}

// This function is taken (with some neccessary variable names modified) from Renchi Yang's GitHub Repo.
// However, there are some bugs in it, e.g., only one svec / tvec will mixup residue value in different layers.
double SinglePairAlgorithms::Power_Iteration_AMC_GEER(Graph::size_type s, Graph::size_type t, Algorithms::length_type L_max, double eps, double p_f, std::vector<std::vector<double>> & svecs, std::vector<std::vector<double>> & tvecs, double & smax, double & tmax, Algorithms::length_type & ellb, Algorithms::size_type tau)
{
	svecs[0][s] = 1.0;
	std::set<Graph::size_type> S;
	S.insert(s);

	tvecs[0][t] = 1.0;
	std::set<Graph::size_type> T;
	T.insert(t);

	ellb = 0;
	double pi = 0.0;

	Algorithms::length_type ell = L_max;
	double epsilon = eps;
	Graph::degree_type ds = G.get_d(s);
	Graph::degree_type dt = G.get_d(t);
	double delta = p_f;

	const std::vector<std::vector<Graph::size_type>> & adjacency_list = G.get_adjacency_list();

	for(Algorithms::length_type i = 0; i <= ell; i++)
	{
		smax = 0.0;
		tmax = 0.0;
		std::set<Graph::size_type> tmpS;
		std::set<Graph::size_type> tmpT;

		for(std::set<Graph::size_type>::iterator it = S.begin(); it != S.end(); ++it)
		{
			Graph::size_type v = *it;
			double residue = svecs[i % 2][v];
			if(v == s) 
				pi += residue / (double)G.get_d(v);
			svecs[i % 2][v] = 0.0;
			for(const auto & u: adjacency_list[v])
			{
				tmpS.insert(u);
				double update = residue / (double)G.get_d(u);
				// Bug! Mixup residue values in different layers.
				// svec[u] += update;
				// if(svec[u] > smax)
				// {
				// 	smax = svec[u];
				// }
				svecs[(i + 1) % 2][u] += update;
				if(svecs[(i + 1) % 2][u] > smax)
				{
				 	smax = svecs[(i + 1) % 2][u];
				}
			}
		}
		S = tmpS;

		for(std::set<Graph::size_type>::iterator it = T.begin(); it != T.end(); ++it)
		{
			Graph::size_type v = *it;

			double residue = tvecs[i % 2][v];
			if(v == t) 
				pi += residue / (double)G.get_d(v);
			if(v == s)
				pi -= 2 * residue / (double)G.get_d(t);
			tvecs[i % 2][v] = 0.0;
			for(const auto & u: adjacency_list[v])
			{
				tmpT.insert(u);
				double update = residue / (double)G.get_d(u);
				// Bug! Mixup residue values in different layers.
				// tvec[u] += update;
				// if(tvec[u] > tmax)
				// {
				// 	tmax = tvec[u];
				// }
				tvecs[(i + 1) % 2][u] += update;
				if(tvecs[(i + 1) % 2][u] > tmax)
				{
					tmax = tvecs[(i + 1) % 2][u];
				}
			}
		}
		T = tmpT;

		ellb += 1;

		double picost = 0.0;
		for(std::set<Graph::size_type>::iterator it = S.begin(); it != S.end(); ++it)
		{
			uint v = *it;
			picost += G.get_d(v);
		}

		for(std::set<Graph::size_type>::iterator it = T.begin(); it != T.end(); ++it)
		{
			uint v = *it;
			picost += G.get_d(v);
		}

		Algorithms::length_type lenwalk = ell - ellb;
		smax = 2 * smax;
		tmax = 2 * tmax;
		Algorithms::size_type eta_star_s = Algorithms::size_type(8 * std::log(3 * (tau + 1) / delta) * (lenwalk + 1) * (lenwalk + 1) * std::pow(std::max(smax / ds, tmax / dt), 2) / epsilon / epsilon);
		Algorithms::size_type eta_star_t = Algorithms::size_type(8 * std::log(3 * (tau + 1) / delta) * (lenwalk + 1) * (lenwalk + 1) * tmax * tmax / epsilon / epsilon / dt / dt);
		double mccost = (eta_star_s + eta_star_t);

		if(picost > mccost)
		{
			break;
		}
	}
	return pi;
}

double SinglePairAlgorithms::GEER(Graph::size_type s, Graph::size_type t, Algorithms::length_type L_max, double eps, double p_f, Algorithms::size_type tau)
{
	std::vector<std::vector<double>> & svecs = this -> mat_s_2_times_n;
	std::vector<std::vector<double>> & tvecs = this -> mat_t_2_times_n;

	this -> reset_zeros("GEER");

	Algorithms::length_type ellb = 0;
	double smax = 0.0;
	double tmax = 0.0;

	double rb = Power_Iteration_AMC_GEER(s, t, L_max, eps, p_f, svecs, tvecs, smax, tmax, ellb, tau);

	double rf = 0.0;
	if(L_max > ellb)
	{
		Algorithms::length_type ellf = L_max - ellb;
		rf = AdaptiveMonteCarlo(s, t, ellf, eps, p_f, svecs[ellb % 2], tvecs[ellb % 2], smax, tmax, tau);
	}

	double er = rf + rb;
	return er;
}

// From SIGMOD'23 paper "Efficient Resistance Distance Computation: the Power of Landmark-based Approaches".
double SinglePairAlgorithms::AbWalk(Graph::size_type s, Graph::size_type t, Graph::size_type landmark, Algorithms::size_type num_samples)
{
	this -> reset_zeros("AbWalk");

	// Random device as seeding source.
	std::random_device rd;
	// Seeding.
	std::mt19937 gen(rd());

	if(s == t)
	{
		return 0.0;
	}
	double rst = 0.0;
	if (t == landmark)
	{
		for (Algorithms::size_type i = 0; i < num_samples; ++i)
		{
			double tao_s = 0.0;
			Graph::size_type current = s;
			while (current != landmark)
			{
				if(current == s)
				{
					tao_s += 1.0;
				}
				std::uniform_int_distribution<Graph::degree_type> dis(0, G.get_d(current) - 1);
				current = G.get_adjacency_list()[current][dis(gen)];
			}
			rst += tao_s / (G.get_d(s) * num_samples);
		}
		return rst;
	}
	else if (s == landmark)
	{
		for (Algorithms::size_type i = 0; i < num_samples; ++i)
		{
			double tao_t = 0.0;
			Graph::size_type current = t;
			while (current != landmark)
			{
				if (current == t)
				{
					tao_t += 1.0;
				}
				std::uniform_int_distribution<Graph::degree_type> dis(0, G.get_d(current) - 1);
				current = G.get_adjacency_list()[current][dis(gen)];
			}
			rst += tao_t / (G.get_d(t) * num_samples);
		}
		return rst;
	}
	else
	{
		for (Algorithms::size_type i = 0; i < num_samples; ++i)
		{
			double tao_ss = 0.0;
			double tao_st = 0.0;
			Graph::size_type current = s;
			while (current != landmark)
			{
				if (current == s)
				{
					tao_ss += 1.0;
				}
				if (current == t)
				{
					tao_st += 1.0;
				}
				std::uniform_int_distribution<Graph::degree_type> dis(0, G.get_d(current) - 1);
				current = G.get_adjacency_list()[current][dis(gen)];
			}
			double tao_ts = 0.0;
			double tao_tt = 0.0;
			current = t;
			while (current != landmark)
			{
				if (current == s)
				{
					tao_ts += 1.0;
				}
				if (current == t)
				{
					tao_tt += 1.0;
				}
				std::uniform_int_distribution<Graph::degree_type> dis(0, G.get_d(current) - 1);
				current = G.get_adjacency_list()[current][dis(gen)];;
			}
			rst += (tao_ss / (G.get_d(s) * num_samples)) + (tao_tt / (G.get_d(t) * num_samples)) - (tao_st / (G.get_d(t) * num_samples)) - (tao_ts / (G.get_d(s) * num_samples));
		}
		return rst;
	}
}

void SinglePairAlgorithms::SS_landmark_local_push_res(Graph::size_type s, Graph::size_type landmark, double r_max, std::vector<double> & res, std::vector<double> & result)
{
	std::vector<bool> & isInQueue = this -> vec_bool_n;

	std::queue<Graph::size_type> r_q;

	const std::vector<std::vector<Graph::size_type>> & adjacency_list = G.get_adjacency_list();

	double push_threshold = r_max;

	res[s] = 1.0;
	r_q.push(s);
	isInQueue[s] = true;
	
	while(r_q.size() > 0)
	{
		Graph::size_type current_node = r_q.front();
		r_q.pop();
		isInQueue[current_node] = false;

		result[current_node] += res[current_node];
		double increment = res[current_node];
		res[current_node] = 0;
		for(size_t i = 0; i < adjacency_list[current_node].size(); i++)
		{
			if(adjacency_list[current_node][i] != landmark)
			{
				Graph::size_type updateNode = adjacency_list[current_node][i];
				res[updateNode] += increment / G.get_d(current_node);
				if(!isInQueue[updateNode] && res[updateNode] >= push_threshold)
				{
					r_q.push(updateNode);
					isInQueue[updateNode] = true;
				}
			}
		}
	}
}

double SinglePairAlgorithms::Bipush(Graph::size_type s, Graph::size_type t, Graph::size_type landmark, Algorithms::size_type num_samples, double r_max)
{
	double r = 0.0;
	std::vector<double> & X_s = this -> vec_s_1_n;
	std::vector<double> & X_t = this -> vec_t_1_n;
	std::vector<double> & R_s = this -> vec_s_2_n;
	std::vector<double> & R_t = this -> vec_t_2_n;

	this -> reset_zeros("Bipush");

	// Random device as seeding source.
	std::random_device rd;
	// Seeding.
	std::mt19937 gen(rd());

	if(s == t)
	{
		return 0.0;
	}
	else if(s == landmark)
	{
		SS_landmark_local_push_res(t, landmark, r_max, R_t, X_t);
		r = X_t[t] / G.get_d(t);
		for(Algorithms::size_type i = 0; i < num_samples; i++)
		{
			Graph::size_type v = t;
			for(;;)
			{
				r += R_t[v] / G.get_d(v) / num_samples;
				std::uniform_int_distribution<Graph::degree_type> dis(0, G.get_d(v) - 1);
				v = G.get_adjacency_list()[v][dis(gen)];
				if(v == landmark)
					break;
			}
		}
		return r;
	}
	else if(t == landmark)
	{
		SS_landmark_local_push_res(s, landmark, r_max, R_s, X_s);
		r = X_s[s] / G.get_d(s);
		for(Algorithms::size_type i = 0; i < num_samples; i++)
		{
			Graph::size_type v = s;
			for(;;)
			{
				r += R_s[v] / G.get_d(v) / num_samples;
				std::uniform_int_distribution<Graph::degree_type> dis(0, G.get_d(v) - 1);
				v = G.get_adjacency_list()[v][dis(gen)];
				if(v == landmark)
					break;
			}
		}
		return r;
	}
	SS_landmark_local_push_res(s, landmark, r_max, R_s, X_s);
	SS_landmark_local_push_res(t, landmark, r_max, R_t, X_t);

	r = X_s[s] / G.get_d(s) - X_s[t] / G.get_d(t) + X_t[t] / G.get_d(t) - X_t[s] / G.get_d(s);

	// Run random walks.
	for (Algorithms::size_type i = 0; i < num_samples; i++)
	{
		// Random walk from s.
		Graph::size_type v = s;
		for (;;)
		{
			r = r + R_s[v] / G.get_d(v) / num_samples - R_t[v] / G.get_d(v) / num_samples;
			std::uniform_int_distribution<Graph::degree_type> dis(0, G.get_d(v) - 1);
			v = G.get_adjacency_list()[v][dis(gen)];
			if(v == landmark)
				break;
		}
		// Random walk from t.
		v = t;
		for (;;)
		{
			r = r + R_t[v] / G.get_d(v) / num_samples - R_s[v] / G.get_d(v) / num_samples;
			std::uniform_int_distribution<Graph::degree_type> dis(0, G.get_d(v) - 1);
			v = G.get_adjacency_list()[v][dis(gen)];
			if (v == landmark)
				break;
		}
	}

	return r;
}

void SinglePairAlgorithms::SS_landmark_local_push(Graph::size_type s, Graph::size_type landmark, double r_max, std::vector<double> & result)
{
	std::vector<double> & res = this -> vec_1_n;
	std::vector<bool> & isInQueue = this -> vec_bool_n;
	std::queue<Graph::size_type> r_q;

	const std::vector<std::vector<Graph::size_type>> & adjacency_list = G.get_adjacency_list();

	double push_threshold = r_max;

	res[s] = 1.0;
	r_q.push(s);
	isInQueue[s] = true;
	
	while(r_q.size() > 0)
	{
		Graph::size_type current_node = r_q.front();
		r_q.pop();
		isInQueue[current_node] = false;

		result[current_node] += res[current_node];
		double increment = res[current_node];
		res[current_node] = 0;
		for(size_t i = 0; i < adjacency_list[current_node].size(); i++)
		{
			if(adjacency_list[current_node][i] != landmark)
			{
				Graph::size_type updateNode = adjacency_list[current_node][i];
				res[updateNode] += increment / G.get_d(current_node);
				if(!isInQueue[updateNode] && res[updateNode] >= push_threshold)
				{
					r_q.push(updateNode);
					isInQueue[updateNode] = true;
				}
			}
		}
	}
}

double SinglePairAlgorithms::Push(Graph::size_type s, Graph::size_type t, Graph::size_type landmark, double r_max)
{
	double r = 0.0;
	std::vector<double> & X_s = this -> vec_s_1_n;
	std::vector<double> & X_t = this -> vec_t_1_n;

	this -> reset_zeros("Push");

	if(s == t)
	{
		return r;
	}
	else if(s == landmark)
	{
		this -> SS_landmark_local_push(t, landmark, r_max, X_t);
		r = X_t[t] / G.get_d(t);
		return r;
	}
	else if(t == landmark)
	{
		this -> SS_landmark_local_push(s, landmark, r_max, X_s);
		r = X_s[s] / G.get_d(s);
		return r;
	}
	this -> SS_landmark_local_push(s, landmark, r_max, X_s);
	this -> SS_landmark_local_push(t, landmark, r_max, X_t);
	
	r = X_s[s] / G.get_d(s) - X_s[t] / G.get_d(t) + X_t[t] / G.get_d(t) - X_t[s] / G.get_d(s);

	return r;
}


// From SIGMOD'24 paper "Efficient and Provable Effective Resistance Computation on Large Graphs: an Index-based Approach".
Eigen::MatrixXd SinglePairAlgorithms::vectorOfVectorsToEigenMatrix(const std::vector<std::vector<double>> & data)
{
	size_t rows = data.size();
	size_t cols = data[0].size();

	Eigen::MatrixXd matrix(rows, cols);
	for (size_t i = 0; i < rows; ++i)
	{
		for(size_t j = 0; j < cols; ++j)
		{
			matrix(i, j) = data[i][j];
		}
	}
	return matrix;
}

Eigen::MatrixXd SinglePairAlgorithms::pinv_eigen_based(Eigen::MatrixXd & origin, const double er)
{
	Eigen::JacobiSVD<Eigen::MatrixXd> svd_holder(origin, Eigen::ComputeThinU | Eigen::ComputeThinV);
	Eigen::MatrixXd U = svd_holder.matrixU();
	Eigen::MatrixXd V = svd_holder.matrixV();
	Eigen::MatrixXd D = svd_holder.singularValues();

	// S
	Eigen::MatrixXd S(V.cols(), U.cols());
	S.setZero();

	for (Eigen::Index i = 0; i < D.size(); ++i)
	{
		if (D(i, 0) > er)
		{
			S(i, i) = 1.0 / D(i, 0);
		}
		else
		{
			S(i, i) = 0.0;
		}
	}

	// pinv_matrix = V * S * U^T
	return V * S * U.transpose();
}

Eigen::MatrixXd SinglePairAlgorithms::pseudo_inverse(const std::vector<std::vector<double>> & data)
{
	Eigen::MatrixXd matrix = vectorOfVectorsToEigenMatrix(data);
	Eigen::MatrixXd pseudoInv = pinv_eigen_based(matrix, 1e-10);
	return pseudoInv;
}

std::vector<std::vector<double>> SinglePairAlgorithms::precompute_Pr(std::vector<Graph::size_type> & vl, std::vector<long long> & index_vl, Algorithms::size_type num_samples)
{
	std::vector<std::vector<double>> Pr(G.get_num_nodes() - vl.size(), std::vector<double>(vl.size(), 0.0));
	// Random device as seeding source.
	std::random_device rd;
	// Seeding.
	std::mt19937 gen(rd());

	for (Algorithms::size_type times = 1; times <= num_samples; ++times)
	{
		for (Graph::size_type i = 0; i < G.get_num_nodes(); ++i)
		{
			if(!(index_vl[i] > 0))
			{
				Graph::size_type current = i;
				while (!(index_vl[current] > 0))
				{
					std::uniform_int_distribution<Graph::degree_type> dis(0, G.get_d(current) - 1);
					current = G.get_adjacency_list()[current][dis(gen)];
				}
				long long ij = index_vl[current] - 1;
				long long ii = -index_vl[i];
				Pr[ii][ij] += 1.0 / num_samples;
			}            
		}
	}
	return Pr;
}

Eigen::MatrixXd SinglePairAlgorithms::precompute_SchurComplementInverse(std::vector<Graph::size_type> & vl, std::vector<long long> & index_vl, std::vector<std::vector<double>> & P)
{
	std::vector<std::vector<double>> SC(vl.size(), std::vector<double>(vl.size(), 0.0));

	for (size_t i = 0; i < vl.size(); ++i)
	{
		Graph::size_type v = vl[i];
		for (Graph::size_type neighbor : G.get_adjacency_list()[v])
		{
			if (!(index_vl[neighbor] > 0))
			{
				long long iu = -index_vl[neighbor];
				for (size_t j = 0; j < vl.size(); ++j)
				{
					SC[i][j] -= P[iu][j];
				}
			}
			else
			{
				long long iv = index_vl[neighbor] - 1;
				SC[i][iv] += -1.0;
			}
		}
		SC[i][i] += G.get_d(v);
	}

	return pseudo_inverse(SC);
}

std::vector<std::vector<double>> SinglePairAlgorithms::combine_P_SCInverse(std::vector<std::vector<double>> & P, Eigen::MatrixXd SCInverse, std::vector<long long> & index_vl)
{
	if(!(SCInverse.rows() == P[0].size()))
	{
		throw;
	}
	Eigen::Index selectNum = SCInverse.rows();
	Eigen::Index n = P.size() + selectNum;
	std::vector<std::vector<double>> Index(n, std::vector<double>(selectNum, 0.0));
	Eigen::Index v_offset = 0;
	for(Eigen::Index i = 0; i < n; i++)
	{
		if(index_vl[i] > 0)
		{
			for(Eigen::Index j = 0; j < selectNum; j++)
			{
				Index[i][j] = SCInverse(index_vl[i] - 1, j);
			}
			v_offset += 1;
		}
		else
		{
			Index[i] = P[i - v_offset];
		}
	}
	return Index;
}

void SinglePairAlgorithms::vl_absorbed_push(Graph::size_type s, const std::vector<long long> & index_vl, double r_max, std::vector<double> & r, std::vector<double> & q)
{
	// std::vector<double> r(G.get_num_nodes(), 0.0);
	// std::vector<double> q(G.get_num_nodes(), 0.0);
	r[s] = 1.0;
	std::vector<Graph::size_type> push_queue = {s};
	// std::vector<bool> is_inQueue(G.get_num_nodes(), false);
	std::vector<bool> & is_inQueue = this -> vec_bool_n;
	is_inQueue[s] = true;

	while(!push_queue.empty())
	{
		Graph::size_type u = push_queue.back();
		push_queue.pop_back();
		is_inQueue[u] = false;
		q[u] += r[u];

		for (size_t i = 0; i < G.get_adjacency_list()[u].size(); ++i)
		{
			Graph::size_type neighbor = G.get_adjacency_list()[u][i];
			if (!(index_vl[neighbor] > 0))
			{
				r[neighbor] += r[u] / G.get_d(u);
				if (r[neighbor] >= G.get_d(neighbor) * r_max && !is_inQueue[neighbor])
				{
					push_queue.push_back(neighbor);
					is_inQueue[neighbor] = true;
				}
			}
		}
		r[u] = 0.0;
	}
}


double SinglePairAlgorithms::RW_vl(Graph::size_type s, Graph::size_type t, Graph::size_type num_landmarks, Algorithms::size_type num_samples)
{
	this -> reset_zeros("RW-vl");

	// Random device as seeding source.
	std::random_device rd;
	// Seeding.
	std::mt19937 gen(rd());

	if ((index_vl[s] > 0) && (index_vl[t] > 0))
	{
		long long is = index_vl[s] - 1;
		long long it = index_vl[t] - 1;
		return SCInverse(is, is) + SCInverse(it, it) - 2 * SCInverse(is, it);
	}
	else if (!(index_vl[s] > 0) && (index_vl[t] > 0))
	{
		long long it = index_vl[t] - 1;

		double qs_s = 0;
		for(Algorithms::size_type i = 0; i < num_samples; i++)
		{
			Graph::size_type current = s;
			while(!(index_vl[current] > 0))
			{
				if(current == s)
				{
					qs_s += 1.0 / num_samples;
				}
				std::uniform_int_distribution<Graph::degree_type> dis(0, G.get_d(current) - 1);
				current = G.get_adjacency_list()[current][dis(gen)];
			}
		}

		std::vector<double> ps_sci(num_landmarks, 0.0);
		for(Graph::size_type i = 0;i < num_landmarks; ++i)
		{
			for(Graph::size_type j = 0; j < num_landmarks; j++)
			{
				ps_sci[i] += Index[s][j] * SCInverse(j, i);
			}
		}
		double ps_sci_ps = 0.0, ps_scit = 0.0;
		for (Graph::size_type j = 0; j < num_landmarks; ++j)
		{
			ps_sci_ps += ps_sci[j] * Index[s][j];
			ps_scit += Index[s][j] * SCInverse(j, it);
		}
		return qs_s / G.get_d(s) + ps_sci_ps + SCInverse(it, it) - 2 * ps_scit;
	}
	else if((index_vl[s] > 0) && !(index_vl[t] > 0))
	{
		long long is = index_vl[s] - 1;
		
		double qt_t = 0;
		for(Algorithms::size_type i = 0; i < num_samples; i++)
		{
			Graph::size_type current = t;
			while(!(index_vl[current] > 0))
			{
				if(current == t)
				{
					qt_t += 1.0 / num_samples;
				}
				std::uniform_int_distribution<Graph::degree_type> dis(0, G.get_d(current) - 1);
				current = G.get_adjacency_list()[current][dis(gen)];
			}
		}

		std::vector<double> pt_sci(num_landmarks, 0.0);
		for (Graph::size_type i = 0; i < num_landmarks; ++i)
		{
			for(Graph::size_type j = 0; j < num_landmarks; j++)
			{
				pt_sci[i] += Index[t][j] * SCInverse(j, i);
			}
		}
		double pt_sci_pt = 0.0, pt_scis = 0.0;
		for (Graph::size_type j = 0; j < num_landmarks; ++j)
		{
			pt_sci_pt += pt_sci[j] * Index[t][j];
			pt_scis += Index[t][j] * SCInverse(j, is);
		}
		return qt_t / G.get_d(t) + pt_sci_pt + SCInverse(is, is) - 2 * pt_scis;
	}
	else if (!(index_vl[s] > 0) && !(index_vl[t] > 0))
	{
		double qs_s = 0;
		double qs_t = 0;
		for(Algorithms::size_type i = 0; i < num_samples; i++)
		{
			Graph::size_type current = s;
			while(!(index_vl[current] > 0))
			{
				if(current == s)
				{
					qs_s += 1.0 / num_samples;
				}
				if(current == t)
				{
					qs_t += 1.0 / num_samples;
				}
				std::uniform_int_distribution<Graph::degree_type> dis(0, G.get_d(current) - 1);
				current = G.get_adjacency_list()[current][dis(gen)];
			}
		}

		double qt_t = 0;
		double qt_s = 0;
		for(Algorithms::size_type i = 0; i < num_samples; i++)
		{
			Graph::size_type current = t;
			while(!(index_vl[current] > 0))
			{
				if(current == t)
				{
					qt_t += 1.0 / num_samples;
				}
				if(current == s)
				{
					qt_s += 1.0 / num_samples;
				}
				std::uniform_int_distribution<Graph::degree_type> dis(0, G.get_d(current) - 1);
				current = G.get_adjacency_list()[current][dis(gen)];
			}
		}

		std::vector<double> ps_sci(num_landmarks, 0.0);
		std::vector<double> pt_sci(num_landmarks, 0.0);
		for (Graph::size_type i = 0; i < num_landmarks; ++i)
		{
			for(Graph::size_type j = 0; j < num_landmarks; j++)
			{
				ps_sci[i] += Index[s][j] * SCInverse(j, i);
				pt_sci[i] += Index[t][j] * SCInverse(j, i);
			}
		}
		double ps_sci_ps = 0.0, pt_sci_pt = 0.0, ps_sci_pt = 0.0;
		for (Graph::size_type j = 0; j < num_landmarks; ++j)
		{
			ps_sci_ps += ps_sci[j] * Index[s][j];
			pt_sci_pt += pt_sci[j] * Index[t][j];
			ps_sci_pt += ps_sci[j] * Index[t][j];
		}
		return qs_s / G.get_d(s) + qt_t / G.get_d(t) -  qs_t / G.get_d(t) - qt_s / G.get_d(s) + ps_sci_ps + pt_sci_pt - 2 * ps_sci_pt;
	}
	std::cout << "Unexpected Return!!" << std::endl;
	return 0.0;
}

double SinglePairAlgorithms::Push_vl(Graph::size_type s, Graph::size_type t, Graph::size_type num_landmarks, double r_max, Algorithms::size_type num_samples)
{
	std::vector<double> & qs = this -> vec_s_1_n;
	std::vector<double> & rs = this -> vec_s_2_n;
	std::vector<double> & qt = this -> vec_t_1_n;
	std::vector<double> & rt = this -> vec_t_2_n;

	this -> reset_zeros("Push-vl");

	if ((index_vl[s] > 0) && (index_vl[t] > 0))
	{
		long long is = index_vl[s] - 1;
		long long it = index_vl[t] - 1;
		return SCInverse(is, is) + SCInverse(it, it) - 2 * SCInverse(is, it);
	}
	else if (!(index_vl[s] > 0) && (index_vl[t] > 0))
	{
		long long it = index_vl[t] - 1;
		this -> vl_absorbed_push(s, index_vl, r_max, rs, qs);
		std::vector<double> ps_sci(num_landmarks, 0.0);
		for (Graph::size_type i = 0; i < num_landmarks; ++i)
		{
			for(Graph::size_type j = 0; j < num_landmarks; j++)
			{
				ps_sci[i] += Index[s][j] * SCInverse(j, i);
			}
		}
		double ps_sci_ps = 0.0, ps_scit = 0.0;
		for(Graph::size_type j = 0; j < num_landmarks; ++j)
		{
			ps_sci_ps += ps_sci[j] * Index[s][j];
			ps_scit += Index[s][j] * SCInverse(j, it);
		}
		return qs[s] / G.get_d(s) + ps_sci_ps + SCInverse(it, it) - 2 * ps_scit;
	}
	else if ((index_vl[s] > 0) && !(index_vl[t] > 0))
	{
		long long is = index_vl[s] - 1;
		vl_absorbed_push(t, index_vl, r_max, rt, qt);
		std::vector<double> pt_sci(num_landmarks, 0.0);
		for(Graph::size_type i = 0;i < num_landmarks; ++i)
		{
			for(Graph::size_type j = 0; j < num_landmarks; j++)
			{
				pt_sci[i] += Index[t][j] * SCInverse(j, i);
			}
		}
		double pt_sci_pt = 0.0, pt_scis = 0.0;
		for (Graph::size_type j = 0; j < num_landmarks; ++j)
		{
			pt_sci_pt += pt_sci[j] * Index[t][j];
			pt_scis += Index[t][j] * SCInverse(j, is);
		}
		return qt[t] / G.get_d(t) + pt_sci_pt + SCInverse(is, is) - 2 * pt_scis;
	}
	else if(!(index_vl[s] > 0) && !(index_vl[t] > 0))
	{
		vl_absorbed_push(s, index_vl, r_max, rs, qs);
		vl_absorbed_push(t, index_vl, r_max, rt, qt);
		std::vector<double> ps_sci(num_landmarks, 0.0);
		std::vector<double> pt_sci(num_landmarks, 0.0);
		for(Graph::size_type i = 0; i < num_landmarks; ++i)
		{
			for(Graph::size_type j = 0; j < num_landmarks; j++)
			{
				ps_sci[i] += Index[s][j] * SCInverse(j, i);
				pt_sci[i] += Index[t][j] * SCInverse(j, i);
			}
		}
		double ps_sci_ps = 0.0, pt_sci_pt = 0.0, ps_sci_pt = 0.0;
		for(Graph::size_type j = 0; j < num_landmarks; ++j)
		{
			ps_sci_ps += ps_sci[j] * Index[s][j];
			pt_sci_pt += pt_sci[j] * Index[t][j];
			ps_sci_pt += ps_sci[j] * Index[t][j];
		}
		return qs[s] / G.get_d(s) + qt[t] / G.get_d(t) - qs[t] / G.get_d(t) - qt[s] / G.get_d(s) + ps_sci_ps + pt_sci_pt - 2 * ps_sci_pt;
	}
	std::cout << "Unexpected Return!!" << std::endl;
	return 0.0;
}

double SinglePairAlgorithms::Bipush_vl(Graph::size_type s, Graph::size_type t, Graph::size_type num_landmarks, double r_max, Algorithms::size_type num_samples)
{
	std::vector<double> & qs = this -> vec_s_1_n;
	std::vector<double> & rs = this -> vec_s_2_n;
	std::vector<double> & qt = this -> vec_t_1_n;
	std::vector<double> & rt = this -> vec_t_2_n;

	this -> reset_zeros("Bipush-vl");

	// Random device as seeding source.
	std::random_device rd;
	// Seeding.
	std::mt19937 gen(rd());
	
	if((index_vl[s] > 0) && (index_vl[t] > 0))
	{
		long long is = index_vl[s] - 1;
		long long it = index_vl[t] - 1;
		return SCInverse(is, is) + SCInverse(it, it) - 2 * SCInverse(is, it);
	}
	else if(!(index_vl[s] > 0) && (index_vl[t] > 0))
	{
		long long it = index_vl[t] - 1;
		vl_absorbed_push(s, index_vl, r_max, rs, qs);
		std::vector<double> ps_sci(num_landmarks, 0.0);
		for(Graph::size_type i = 0; i < num_landmarks; ++i)
		{
			for(Graph::size_type j = 0; j < num_landmarks; j++)
			{
				ps_sci[i] += Index[s][j] * SCInverse(j, i);
			}
		}
		double ps_sci_ps = 0.0, ps_scit = 0.0;
		for(Graph::size_type j = 0; j < num_landmarks; ++j)
		{
			ps_sci_ps += ps_sci[j] * Index[s][j];
			ps_scit += Index[s][j] * SCInverse(j, it);
		}
		double rst = qs[s] / G.get_d(s) + ps_sci_ps + SCInverse(it, it) - 2 * ps_scit;

		for(Algorithms::size_type i = 0; i < num_samples; ++i)
		{
			Graph::size_type current = s;
			while (!(index_vl[current] > 0))
			{
				rst += rs[current] / (G.get_d(current) * num_samples);
				std::uniform_int_distribution<Graph::degree_type> dis(0, G.get_d(current) - 1);
				current = G.get_adjacency_list()[current][dis(gen)];
			}
		}
		return rst;
	}
	else if(index_vl[s] > 0 && !(index_vl[t] > 0))
	{
		long long is = index_vl[s] - 1;
		
		SinglePairAlgorithms::vl_absorbed_push(t, index_vl, r_max, rt, qt);
		std::vector<double> pt_sci(num_landmarks, 0.0);
		for(Graph::size_type i = 0; i < num_landmarks; ++i)
		{
			for(Graph::size_type j = 0; j < num_landmarks; j++)
			{
				pt_sci[i] += Index[t][j] * SCInverse(j, i);
			}
		}
		double pt_sci_pt = 0.0, pt_scis = 0.0;
		for(Graph::size_type j = 0; j < num_landmarks; ++j)
		{
			pt_sci_pt += pt_sci[j] * Index[t][j];
			pt_scis += Index[t][j] * SCInverse(j,is);
		}
		double rst = qt[t] / G.get_d(t) + pt_sci_pt + SCInverse(is, is) - 2 * pt_scis;

		for(Algorithms::size_type i = 0; i < num_samples; ++i)
		{
			Graph::size_type current = t;
			while (!(index_vl[current] > 0))
			{
				rst += rt[current] / (G.get_d(current) * num_samples);
				std::uniform_int_distribution<Graph::degree_type> dis(0, G.get_d(current) - 1);
				current = G.get_adjacency_list()[current][dis(gen)];
			}
		}
		return rst;
	}
	else if(!(index_vl[s] > 0) && !(index_vl[t] > 0))
	{
		vl_absorbed_push(s, index_vl, r_max, rs, qs);
		vl_absorbed_push(t, index_vl, r_max, rt, qt);
		std::vector<double> ps_sci(num_landmarks, 0.0);
		std::vector<double> pt_sci(num_landmarks, 0.0);
		for(Graph::size_type i = 0; i < num_landmarks; ++i)
		{
			for(Graph::size_type j = 0; j < num_landmarks; j++)
			{
				ps_sci[i] += Index[s][j] * SCInverse(j, i);
				pt_sci[i] += Index[t][j] * SCInverse(j, i);
			}
		}
		double ps_sci_ps = 0.0, pt_sci_pt = 0.0, ps_sci_pt = 0.0;
		for (Graph::size_type j = 0; j < num_landmarks; ++j)
		{
			ps_sci_ps += ps_sci[j] * Index[s][j];
			pt_sci_pt += pt_sci[j] * Index[t][j];
			ps_sci_pt += ps_sci[j] * Index[t][j];
		}
		double rst = qs[s] / G.get_d(s) + qt[t] / G.get_d(t) - qs[t] / G.get_d(t) - qt[s] / G.get_d(s) + ps_sci_ps + pt_sci_pt - 2 * ps_sci_pt;

		// double current_time = get_current_time_sec_method();
		for (Algorithms::size_type i = 0; i < num_samples; ++i)
		{
			Graph::size_type current = s;
			while (!(index_vl[current] > 0))
			{
				rst += (rs[current] - rt[current]) / (G.get_d(current) * num_samples);
				std::uniform_int_distribution<Graph::degree_type> dis(0, G.get_d(current) - 1);
				current = G.get_adjacency_list()[current][dis(gen)];
			}
		}
		for (Algorithms::size_type i = 0; i < num_samples; ++i)
		{
			Graph::size_type current = t;
			while (!(index_vl[current] > 0))
			{
				rst += (rt[current] - rs[current]) / (G.get_d(current) * num_samples);
				std::uniform_int_distribution<Graph::degree_type> dis(0, G.get_d(current) - 1);
				current = G.get_adjacency_list()[current][dis(gen)];
			}
		}
		return rst;
	}
	std::cout << "unexpected return!!" << std::endl;
	throw;
	return 0.0;
}



// Ours.
double SinglePairAlgorithms::BiSPER_(Graph::size_type s, Graph::size_type t, Algorithms::length_type L_max, double eps, double p_f, double r_max)
{
	Graph::degree_type d_s = G.get_d(s);
	Graph::degree_type d_t = G.get_d(t);
	Graph::degree_type d = std::min(d_s, d_t);

	// Holding r^{(k)}(u) / d(u) > r_max and r^{(k + 1)}(u) / d(u) > r_max node_id.
	std::array<std::queue<Graph::size_type>, 2> active_node_s;
	std::array<std::queue<Graph::size_type>, 2> active_node_t;

	// n BITs, each one holds the residues and their prefix sum \sum_{k = 0}^{\ell} r^{(k)}(u) / d(u).
	std::vector<SinglePairAlgorithms_Tools::BIT<double>> & Q_s = BITs_s;
	std::vector<SinglePairAlgorithms_Tools::BIT<double>> & Q_t = BITs_t;

	// Push phase.

	// \sum_{\ell = 0}^{L_max} q_s^{(\ell)}(s) / d(s) - \sum_{\ell = 0}^{L_max} q_s^{(\ell)}(t) / d(t) 
	double sum_s_qs_div_ds_minus_qt_div_dt = 0.0;
	// \sum_{\ell = 0}^{L_max} q_t^{(\ell)}(t) / d(t) - \sum_{\ell = 0}^{L_max} q_t^{(\ell)}(s) / d(s) 
	double sum_t_qt_div_dt_minus_qs_div_ds = 0.0;
	
	double sum_s_q_all = 0.0;
	double sum_t_q_all = 0.0;

	// r_s[s] = 1.0;
	Q_s[s].increase(0, 1.0 / d_s);
	BITs_visited_s[s] = true;
	if(1.0 > d_s * r_max)
	{
		active_node_s[0].push(s);
	}
	else
	{
		// No push.
		sum_s_qs_div_ds_minus_qt_div_dt = (s == t) ? 0.0 : 1.0 / d_s;
		sum_s_q_all = 1.0;
	}

	// r_t[t] = 1.0;
	Q_t[t].increase(0, 1.0 / d_t);
	BITs_visited_t[t] = true;
	if(1.0 > d_t * r_max)
	{
		active_node_t[0].push(t);
	}
	else
	{
		// No push.
		sum_t_qt_div_dt_minus_qs_div_ds = (s == t) ? 0.0 : 1.0 / d_t;
		sum_t_q_all = 0.0;
	}

	for(Algorithms::size_type l = 0; l <= L_max; l++)
	{
		while(!active_node_s[l % 2].empty())
		{
			// Extract node u.
			Graph::size_type u = active_node_s[l % 2].front();
			active_node_s[l % 2].pop();

			Graph::degree_type d_u = G.get_d(u);
			// r_s[u][l] / d_u
			double r_u_div_d_u = Q_s[u].get_value(l);

			// Push(u, l)
			sum_s_qs_div_ds_minus_qt_div_dt += ((u == s) ? r_u_div_d_u : 0.0) - ((u == t) ? r_u_div_d_u : 0.0);
			sum_s_q_all += r_u_div_d_u * d_u;
			for(const Graph::size_type & v : G.get_adjacency_list()[u])
			{
				Graph::degree_type d_v = G.get_d(v);
				double old_residue_div_d_v = Q_s[v].get_value(l + 1);

				bool is_old_residue_small = (old_residue_div_d_v <= r_max);
				// r_s[v][l + 1] += r_u / d_u;
				Q_s[v].increase(l + 1, r_u_div_d_u / d_v);
				BITs_visited_s[v] = true;
				bool is_new_residue_large = (old_residue_div_d_v + r_u_div_d_u / d_v > r_max);

				// Node v is active.
				if(is_old_residue_small && is_new_residue_large)
				{
					active_node_s[(l + 1) % 2].push(v);
				}
			}
			Q_s[u].increase(l, -r_u_div_d_u);
			// r_s[u][l] = 0.0;
		}
	}

	for(Algorithms::size_type l = 0; l <= L_max; l++)
	{
		while(!active_node_t[l % 2].empty())
		{
			// Extract node u.
			Graph::size_type u = active_node_t[l % 2].front();
			active_node_t[l % 2].pop();

			Graph::degree_type d_u = G.get_d(u);
			// r_t[u][l].
			double r_u_div_d_u = Q_t[u].get_value(l);

			// Push(u, l)
			sum_t_qt_div_dt_minus_qs_div_ds += ((u == t) ? r_u_div_d_u : 0.0) - ((u == s) ? r_u_div_d_u : 0.0);
			sum_t_q_all += r_u_div_d_u * d_u;
			for(const Graph::size_type & v : G.get_adjacency_list()[u])
			{
				Graph::degree_type d_v = G.get_d(v);
				double old_residue_div_d_v = Q_t[v].get_value(l + 1);

				bool is_old_residue_small = (old_residue_div_d_v <= r_max);
				// r_t[v][l + 1] += r_u / d_u;
				Q_t[v].increase(l + 1, r_u_div_d_u / d_v);
				BITs_visited_t[v] = true;
				bool is_new_residue_large = (old_residue_div_d_v + r_u_div_d_u / d_v > r_max);

				// Node v is active.
				if(is_old_residue_small && is_new_residue_large)
				{
					active_node_t[(l + 1) % 2].push(v);
				}
			}
			Q_t[u].increase(l, -r_u_div_d_u);
			// r_t[u][l] = 0.0;
		}
	}
	
	double T_k_bound = std::min(2.0 * (L_max + 1) - (sum_s_q_all + sum_t_q_all), (L_max + 1) * (L_max + 2) * r_max);

	// Monte-Carlo phase.
	// Total number of random walks generated.
	Algorithms::size_type num_rw_upper_bound = 0;
	if(r_max >= 1.0 / d)
		num_rw_upper_bound = Algorithms::size_type(std::max(0.0, std::ceil(8.0 * (L_max + 1) * (L_max + 1) * std::log(2.0 / p_f) / (eps * d * eps * d))));
	else
		num_rw_upper_bound = Algorithms::size_type(std::max(0.0, std::ceil(2.0 * (T_k_bound / eps) * (T_k_bound / eps) * std::log(2.0 / p_f))));

	if(num_rw_upper_bound == 0)
		return sum_s_qs_div_ds_minus_qt_div_dt + sum_t_qt_div_dt_minus_qs_div_ds;
	
	const auto & adjacency_list = G.get_adjacency_list();

	// Random device as seeding source.
	std::random_device rd;
	// Seeding.
	std::mt19937 gen(rd());

	// Sum of random variables T_k.
	double sum_T_k = 0.0;
	// Sum of square of random variables T_k^2.
	double sum_T_k_squared = 0.0;

	for(Algorithms::size_type num_rw = 1; num_rw <= num_rw_upper_bound; num_rw++)
	{
		// Random variable T_k.
		double T_k = 0.0;
		Graph::size_type now_node = t;
		T_k += Q_t[now_node].get_prefix_sum(L_max) - Q_s[now_node].get_prefix_sum(L_max);

		for(Algorithms::length_type l = 1; l <= L_max; l++)
		{
			// Select a random neighbor.
			Graph::degree_type d_node = adjacency_list[now_node].size();
			std::uniform_int_distribution<Graph::degree_type> dis(0, d_node - 1);
			now_node = adjacency_list[now_node][dis(gen)];
			T_k += Q_t[now_node].get_prefix_sum(L_max - l) - Q_s[now_node].get_prefix_sum(L_max - l);
		}

		now_node = s;
		T_k += Q_s[now_node].get_prefix_sum(L_max) - Q_t[now_node].get_prefix_sum(L_max);

		for(Algorithms::length_type l = 1; l <= L_max; l++)
		{
			// Select a random neighbor.
			Graph::degree_type d_node = adjacency_list[now_node].size();
			std::uniform_int_distribution<Graph::degree_type> dis(0, d_node - 1);
			now_node = adjacency_list[now_node][dis(gen)];
			T_k += Q_s[now_node].get_prefix_sum(L_max - l) - Q_t[now_node].get_prefix_sum(L_max - l);
		}

		sum_T_k += T_k;
		sum_T_k_squared += (T_k * T_k);

		double emp_var = sum_T_k_squared / num_rw - (sum_T_k / num_rw) * (sum_T_k / num_rw);
		if(std::sqrt(2.0 * emp_var * std::log(3.0 / p_f) / num_rw) + 6.0 * T_k_bound * std::log(3.0 / p_f) / num_rw <= eps)
		{
			double S = sum_s_qs_div_ds_minus_qt_div_dt + sum_t_qt_div_dt_minus_qs_div_ds + sum_T_k / num_rw;
			return S;
		}
	}

	return sum_s_qs_div_ds_minus_qt_div_dt + sum_t_qt_div_dt_minus_qs_div_ds + sum_T_k / num_rw_upper_bound;
}

double SinglePairAlgorithms::AdaptiveMonteCarlo_BiSPER(Graph::size_type s, Graph::size_type t, Algorithms::length_type L_max, double eps, double p_f)
{
	Graph::degree_type d_s = G.get_d(s);
	Graph::degree_type d_t = G.get_d(t);
	Graph::degree_type d = std::min(d_s, d_t);

	// \sum_{\ell = 0}^{L_max} q_s^{(\ell)}(s) / d(s) - \sum_{\ell = 0}^{L_max} q_s^{(\ell)}(t) / d(t) 
	double sum_s_qs_div_ds_minus_qt_div_dt = (s != t) ? 1.0 / d_s : 0.0;
	// \sum_{\ell = 0}^{L_max} q_t^{(\ell)}(t) / d(t) - \sum_{\ell = 0}^{L_max} q_t^{(\ell)}(s) / d(s) 
	double sum_t_qt_div_dt_minus_qs_div_ds = (s != t) ? 1.0 / d_t : 0.0;
		
	double T_k_bound = std::min(2.0 * L_max, (L_max + 1) * (L_max + 2) * 1.0 / d);

	// Monte-Carlo phase.
	// Total number of random walks generated.
	Algorithms::size_type num_rw_upper_bound = std::ceil(8.0 * (L_max + 1) * (L_max + 1) * std::log(2.0 / p_f) / (eps * d * eps * d));

	if(num_rw_upper_bound == 0)
		return sum_s_qs_div_ds_minus_qt_div_dt + sum_t_qt_div_dt_minus_qs_div_ds;
	
	const auto & adjacency_list = G.get_adjacency_list();

	// Random device as seeding source.
	std::random_device rd;
	// Seeding.
	std::mt19937 gen(rd());

	// Sum of random variables T_k.
	double sum_T_k = 0.0;
	// Sum of square of random variables T_k^2.
	double sum_T_k_squared = 0.0;

	for(Algorithms::size_type num_rw = 1; num_rw <= num_rw_upper_bound; num_rw++)
	{
		// Random variable T_k.
		double T_k = 0.0;
		Graph::size_type now_node = t;
		T_k += (((now_node == t) ? 1.0 / d_t : 0.0) - ((now_node == s) ? 1.0 / d_s : 0.0));

		for(Algorithms::length_type l = 1; l <= L_max; l++)
		{
			// Select a random neighbor.
			Graph::degree_type d_node = adjacency_list[now_node].size();
			std::uniform_int_distribution<Graph::degree_type> dis(0, d_node - 1);
			now_node = adjacency_list[now_node][dis(gen)];
			T_k += (((now_node == t) ? 1.0 / d_t : 0.0) - ((now_node == s) ? 1.0 / d_s : 0.0));
		}

		now_node = s;
		T_k += (((now_node == s) ? 1.0 / d_s : 0.0) - ((now_node == t) ? 1.0 / d_t : 0.0));

		for(Algorithms::length_type l = 1; l <= L_max; l++)
		{
			// Select a random neighbor.
			Graph::degree_type d_node = adjacency_list[now_node].size();
			std::uniform_int_distribution<Graph::degree_type> dis(0, d_node - 1);
			now_node = adjacency_list[now_node][dis(gen)];
			T_k += (((now_node == s) ? 1.0 / d_s : 0.0) - ((now_node == t) ? 1.0 / d_t : 0.0));
		}

		sum_T_k += T_k;
		sum_T_k_squared += (T_k * T_k);

		double emp_var = sum_T_k_squared / num_rw - (sum_T_k / num_rw) * (sum_T_k / num_rw);
		if(std::sqrt(2.0 * emp_var * std::log(3.0 / p_f) / num_rw) + 6.0 * T_k_bound * std::log(3.0 / p_f) / num_rw <= eps)
		{
			double S = sum_s_qs_div_ds_minus_qt_div_dt + sum_t_qt_div_dt_minus_qs_div_ds + sum_T_k / num_rw;
			return S;
		}
	}

	return sum_s_qs_div_ds_minus_qt_div_dt + sum_t_qt_div_dt_minus_qs_div_ds + sum_T_k / num_rw_upper_bound;
}

double SinglePairAlgorithms::BiSPER(Graph::size_type s, Graph::size_type t, Algorithms::length_type L_max, double eps, double p_f)
{
	double r_max = 0.0;
	Graph::size_type m = G.get_num_edges();
	Graph::degree_type d = std::min(G.get_d(s), G.get_d(t));

	// L_max >= m^{1 / 2} * eps * d / (2 * log^{1 / 2}(2 / p_f)) && L_max >= 2 * m^{3 / 4} * eps^{1 / 2} / (3^{3 / 4} * log^{1 / 4}(2 / p_f)))
	if(L_max >= std::sqrt(m) * eps * d / (2.0 * std::sqrt(std::log(2.0 / p_f))) 
		&& L_max >= 2.0 * std::pow(m, 3.0 / 4.0) * std::sqrt(eps) / (std::pow(3.0, 3.0 / 4.0) * std::pow(std::log(2.0 / p_f), 1.0 / 4.0)))
	{
		std::cout << "Warning: degenerated to Power-Iteration." << std::endl;
		r_max = 0.0;
		this -> reset_zeros("Power-Iteration");
		return Power_Iteration_ER(s, t, L_max);
	}
	// d >= 2^{5 / 3} * (L_max + 1)^{1 / 3} * log^{1 / 3}(2 / p_f) / (3^{1 / 2} * eps^{2 / 3}) && d >= 2 * (L_max + 1) * log^{1 / 2}(2 / p_f) / (m^{1 / 2} * eps)
	else if(d >= std::pow(2.0, 5.0 / 3.0) * std::cbrt(L_max + 1) * std::cbrt(std::log(2.0 / p_f)) / (std::sqrt(3.0) * std::pow(eps, 2.0 / 3.0)) 
		&& d >= 2.0 * (L_max + 1) * std::sqrt(std::log(2.0 / p_f)) / (std::sqrt(m) * eps))
	{
		std::cout << "Warning: degenerated to Adaptive-Monte-Carlo." << std::endl;
		r_max = 1.0 / d;
		return AdaptiveMonteCarlo_BiSPER(s, t, L_max, eps, p_f);
	}
	else
	{
		// r_max = eps^{2/3} / (2^{2 / 3} * (L_max + 1)^{2 / 3} * (L_max + 2)^{2 / 3} * log^{1 / 3}(2 / p_f))
		r_max = std::pow(eps, 2.0 / 3.0) / (std::pow(2.0, 2.0 / 3.0) * std::pow(L_max + 1, 2.0 / 3.0) * std::pow(L_max + 2, 2.0 / 3.0) * std::cbrt(std::log(2.0 / p_f)));
	}
	this -> reset_zeros("BiSPER");
	return BiSPER_(s, t, L_max, eps, p_f, r_max);
}
