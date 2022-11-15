#ifndef CLUSTERING_UTIL_H
#define CLUSTERING_UTIL_H

#include <map>
#include <set>
#include <vector>
#include <armadillo>

typedef struct LinkConstraint {

    // must-link
    std::map<int, std::set<int>> ml_graph;
    // cannot-link
    std::map<int, std::set<int>> cl_graph;

} LinkConstraint;

bool sort_by_value(const std::pair<int, double> &a, const std::pair<int, double> &b);
LinkConstraint transitive_closure(std::vector<std::pair<int, int>> &ml, std::vector<std::pair<int, int>> &cl, int n);
void display_graph(std::map<int, std::set<int>> &map);
void print_pairs(std::vector<std::pair<int, int>> &cl_vector);
void add_both(std::map<int, std::set<int>> &graph, int i, int j);
void dfs(int i, std::map<int, std::set<int>> &graph, std::vector<bool> &visited, std::vector<int> &component);
std::map<int, std::set<int>> get_ml_map(int n, std::vector<std::pair<int, int>> &ml);
std::map<int, std::set<int>> build_must_link_map(std::map<int, std::set<int>> &old_ml_map, int i, int j);
std::vector<std::pair<int, int>> build_global_cannot_link_pairs(std::map<int, std::set<int>> &ml_map, std::vector<std::pair<int, int>> &local_cl_pairs);
std::set<int> update_cannot_link(std::vector<std::pair<int, int>> &local_cl_pairs, int i, int j);
std::vector<std::pair<int, int>> build_global_must_link_pairs(std::map<int, std::set<int>> &ml_map);

std::pair<int, int> find_branch_norm(arma::mat &Z);
void print_header_sdp(std::ostream &log_file);
void print_log_sdp(std::ostream &log_file, int n, int n_bin, int node_parent, int node, double lb_parent, double lb,
                   double time, int cp_iter, int cp_flag, int n_ineq, int n_sdp_fixing, double time_sdp_fixing, int n_fixed,
                   double ub, double gub, int i, double node_gap, double gap, int open, bool update);
#endif //CLUSTERING_UTIL_H
