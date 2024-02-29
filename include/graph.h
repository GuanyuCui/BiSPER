#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <map>
#include <queue>
#include <algorithm>

#include "eigen-3.4.0/Eigen/Sparse"

class Graph
{
    public:
        using size_type = uint32_t;
        using degree_type = uint32_t;

        // Constructor & destructor.
        Graph();
        ~Graph() = default;

        // Comparator.
        bool operator==(const Graph & other);

        // Reading & writing.
        bool from_txt(const char * file_name);
        bool to_bin(const char * file_name);
        bool from_bin(const char * file_name);
        
        // Compress adjacency list, such that there is no node_id with zero degree.
        bool compress_adjacency_list();
        // Sort the adjacency list, such that neighbors are sorted in the degree ascending order.
        bool sort_adjacency_list();
        // Print adjacency list.
        void print_adjacency_list();
        // Print statistics.
        void print_statistics();
        // Initialize the transion matrix.
        void init_transition_matrix();

        // Get parameters.
        Graph::size_type get_num_nodes() const;
        Graph::size_type get_num_edges() const;
        
        Graph::degree_type get_d_min() const;
        Graph::degree_type get_d_max() const;
        double get_d_avg() const;

        // Get adjacency_list.
        const std::vector< std::vector<size_type> > & get_adjacency_list() const;
        // Get transition_matrix.
        const bool get_is_transition_matrix_valid() const;
        const Eigen::SparseMatrix<double, 0, int64_t> & get_transition_matrix() const;

        // Get degree.
        const Graph::degree_type get_d(Graph::size_type n) const;

        // Get node with maximum degree.
        const Graph::size_type get_max_degree_node() const;
    private:
        // Adjacency list.
        // Format: <node_id> -> <neighbors>: vector<node_id>
        std::vector< std::vector<size_type> > adjacency_list;
        // Format: <node_id> -> <degree>: degree_type
        // std::vector<degree_type> degree_list;

        // Parameters.
        Graph::size_type num_nodes;
        Graph::size_type num_edges;

        Graph::degree_type d_min;
        Graph::degree_type d_max;

        Graph::size_type max_degree_node;

        Eigen::SparseMatrix<double, 0, int64_t> transition_matrix;
        bool is_transition_matrix_valid;
};

#endif