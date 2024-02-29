#ifndef SINGLE_PAIR_ALGORITHMS_H
#define SINGLE_PAIR_ALGORITHMS_H

#include "algorithms.h"
#include "timer.h"
#include <set>
#include <unordered_set>

namespace SinglePairAlgorithms_Tools
{
	template <typename T>
	class BIT 
	{
		public:
			BIT(size_t size) : size(size), tree(size + 1, T(0)) {}
			void reset_zeros()
			{
				std::fill(tree.begin(), tree.end(), T(0));
			}
			// Update index-th value of contents by the amount of delta.
			void increase(size_t index, T delta)
			{
				// Convert 0-based index to 1-based index for BIT.
				index++;
				while(index <= size)
				{
					tree[index] += delta;
					// Move to the next node in the tree.
					index += (index & -index);
				}
			}
			// Return: sum of values in [0, index].
			T get_prefix_sum(size_t index) const 
			{
				index++;
				T result = T(0);
				while(index > 0)
				{
					result += tree[index];
					// Move to the parent node in the tree.
					index -= (index & -index);
				}
				return result;
			}
			T get_value(size_t index) const 
			{
				if(index == 0)
					return this -> get_prefix_sum(0);
				return this -> get_prefix_sum(index) - get_prefix_sum(index - 1);
			}
		private:
			size_t size;
			std::vector<T> tree;
	};
}

class SinglePairAlgorithms : public Algorithms
{
	public:
		SinglePairAlgorithms(const Graph & G);
		virtual ~SinglePairAlgorithms() = default;
		void init_data_stucture(const std::string & algorithm_name, Algorithms::length_type L_max);
		void reset_zeros(const std::string & algorithm_name);
		// -------------------- Single-pair --------------------
		// Interface for running single-pair effective resistance algorithms.
		// For AMC, GEER, and our BiSPER
		double run(const std::string & algorithm_name, Graph::size_type s, Graph::size_type t, Algorithms::length_type L_max, double eps, double p_f);
		double run_landmark(const std::string & algorithm_name, Graph::size_type s, Graph::size_type t, Algorithms::size_type num_samples, double r_max);

	private:
		// -------------------- Single-pair transition probability tool algorithms --------------------
		// Ground-truth generation.
		double Power_Iteration_ER(Graph::size_type s, Graph::size_type t, Algorithms::length_type L_max);
		double Power_Iteration_ER_Eigen(Graph::size_type s, Graph::size_type t, Algorithms::length_type L_max);

		// From SIGMOD'23 paper "Efficient Estimation of Pairwise Effective Resistance".
		double AdaptiveMonteCarlo(Graph::size_type s, Graph::size_type t, Algorithms::length_type len_walk, double eps, double p_f, std::vector<double> & svec, std::vector<double> & tvec, double psi_s_prime, double psi_t_prime, Algorithms::size_type tau = 5);
		double AMC(Graph::size_type s, Graph::size_type t, Algorithms::length_type L_max, double eps, double p_f, Algorithms::size_type tau = 5);
		double Power_Iteration_AMC_GEER(Graph::size_type s, Graph::size_type t, Algorithms::length_type L_max, double eps, double p_f, std::vector<std::vector<double>> & svecs, std::vector<std::vector<double>> & tvecs, double & smax, double & tmax, Algorithms::length_type & ellb, Algorithms::size_type tau = 5);
		double GEER(Graph::size_type s, Graph::size_type t, Algorithms::length_type L_max, double eps, double p_f, Algorithms::size_type tau = 5);

		// From SIGMOD'23 paper "Efficient Resistance Distance Computation: the Power of Landmark-based Approaches".
		void SS_landmark_local_push_res(Graph::size_type s, Graph::size_type landmark, double r_max, std::vector<double> & res, std::vector<double> & result);
		double Bipush(Graph::size_type s, Graph::size_type t, Graph::size_type landmark, Algorithms::size_type num_samples, double r_max);
		void SS_landmark_local_push(Graph::size_type s, Graph::size_type landmark, double r_max, std::vector<double> & result);
		double Push(Graph::size_type s, Graph::size_type t, Graph::size_type landmark, double r_max);

		// Ours.
		double BiSPER_(Graph::size_type s, Graph::size_type t, Algorithms::length_type L_max, double eps, double p_f, double r_max);
		double AdaptiveMonteCarlo_BiSPER(Graph::size_type s, Graph::size_type t, Algorithms::length_type L_max, double eps, double p_f);
		double BiSPER(Graph::size_type s, Graph::size_type t, Algorithms::length_type L_max, double eps, double p_f);

		// Power-Iteration
		std::vector<std::vector<double>> mat_s_n_times_2;
		std::vector<std::vector<double>> mat_t_n_times_2;

		// GEER
		std::vector<std::vector<double>> mat_s_2_times_n;
		std::vector<std::vector<double>> mat_t_2_times_n;

		// AMC, SS_landmark_local_push_res, SS_landmark_local_push
		std::vector<double> vec_1_n;
		std::vector<double> vec_2_n;
		std::vector<bool> vec_bool_n;

		// Bipush, Push
		std::vector<double> vec_s_1_n;
		std::vector<double> vec_t_1_n;
		std::vector<double> vec_s_2_n;
		std::vector<double> vec_t_2_n;

		// BiSPER
		std::vector<SinglePairAlgorithms_Tools::BIT<double>> BITs_s;
		std::vector<SinglePairAlgorithms_Tools::BIT<double>> BITs_t;

		std::vector<bool> BITs_visited_s;
		std::vector<bool> BITs_visited_t;
};

#endif