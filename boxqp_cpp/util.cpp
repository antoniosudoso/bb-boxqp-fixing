#include "util.h"
#include <iomanip>

// get ml map from ml constraints
std::map<int, std::set<int>> get_ml_map(int n, std::vector<std::pair<int, int>> &ml) {

    std::map<int, std::set<int>> ml_graph;
    std::map<int, std::set<int>> ml_map;

    for (int i = 0; i < n; i++) {
        ml_graph.insert(std::pair<int, std::set<int>> (i, {}));
    }

    for (auto &pair_ml : ml) {
        add_both(ml_graph, pair_ml.first, pair_ml.second);
    }

    int components_counter = 0;
    std::vector<bool> visited(n, false);
    for (int i = 0; i < n; i++) {
        if (!visited[i]) {
            std::vector<int> component;
            dfs(i, ml_graph, visited, component);
            std::set<int> component_set(component.begin(), component.end());
            ml_map[components_counter] = component_set;
            components_counter++;
        }
    }

    return ml_map;
}


// build the must link map history
std::map<int, std::set<int>> build_must_link_map(std::map<int, std::set<int>> &old_ml_map, int i, int j) {

    std::map<int, std::set<int>> ml_map;

    int n = old_ml_map.size();

    // init new map
    for (int k = 0; k < n - 1; k++){
        ml_map.insert(std::pair<int, std::set<int>>(k, {}));
    }

    for (int k = 0; k < j; k++) {
        // copy keys from old map to new map
        for (auto &elem : old_ml_map[k]) {
            ml_map[k].insert(elem);
        }
    }

    for (auto &elem : old_ml_map[j]) {
        ml_map[i].insert(elem);
    }

    // shift
    for (int k = j + 1; k < n; k++) {
        for (auto &elem : old_ml_map[k]) {
            ml_map[k - 1].insert(elem);
        }
    }

    return ml_map;
}

// build global cannot link pairs from local cannot link and map history
std::vector<std::pair<int, int>> build_global_cannot_link_pairs(std::map<int, std::set<int>> &ml_map, std::vector<std::pair<int, int>> &local_cl_pairs) {

    std::vector<std::pair<int, int>> global_cl_pairs;

    for (auto &elem : local_cl_pairs) {

        int key_first = elem.first;
        int key_second = elem.second;
        std::set<int> set_first = ml_map[key_first];
        std::set<int> set_second = ml_map[key_second];
        for (auto &set_elem_1 : set_first) {
            for (auto &set_elem_2 : set_second) {
                global_cl_pairs.emplace_back(set_elem_1, set_elem_2);
            }
        }
    }

    return global_cl_pairs;
}

// update indices in the local cannot link
std::set<int> update_cannot_link(std::vector<std::pair<int, int>> &local_cl_pairs, int i, int j) {
    int size = local_cl_pairs.size();

    for (int k = 0; k < size; k++) {

        if (local_cl_pairs[k].first == j)
            local_cl_pairs[k].first = i;
        if (local_cl_pairs[k].first > j)
            local_cl_pairs[k].first -= 1;

        if (local_cl_pairs[k].second == j)
            local_cl_pairs[k].second = i;
        if (local_cl_pairs[k].second > j)
            local_cl_pairs[k].second -= 1;

        // check right order
        if (local_cl_pairs[k].first >= local_cl_pairs[k].second) {
            int dmy = local_cl_pairs[k].second;
            local_cl_pairs[k].second = local_cl_pairs[k].first;
            local_cl_pairs[k].first = dmy;
        }

    }

    std::set<int> dup_indices = {};
    for (int t = 0; t < size; t++) {
        for (int s = t + 1; s < size; s++) {
            std::pair<int, int> p1 = local_cl_pairs[t];
            std::pair<int, int> p2 = local_cl_pairs[s];
            if (p1.first == p2.first && p1.second == p2.second) {
                dup_indices.insert(s);
            }
        }
    }

    if (!dup_indices.empty()) {
        // delete indices (indices are sorted in the set)
        std::set<int>::reverse_iterator rit;
        for (rit = dup_indices.rbegin(); rit != dup_indices.rend(); rit++)
            local_cl_pairs.erase(local_cl_pairs.begin() + *rit);
    }

    return dup_indices;

}

// build global must link from history map
std::vector<std::pair<int, int>> build_global_must_link_pairs(std::map<int, std::set<int>> &ml_map) {

    std::vector<std::pair<int, int>> ml_pairs;

    for (auto &elem_map : ml_map) {
        std::set<int> set_i = elem_map.second;
        int size = set_i.size();
        if (size > 1) {
            bool is_first = true;
            int first;
            for (auto &elem_set : set_i) {
                if (is_first) {
                    first = elem_set;
                    is_first = false;
                    continue;
                }
                ml_pairs.emplace_back(first, elem_set);
            }
        }
    }

    return ml_pairs;

}

void print_pairs(std::vector<std::pair<int, int>> &cl_vector) {
    for (auto &elem : cl_vector) {
        std::cout << "(" << elem.first << " " << elem.second << ")" << " ";
    }
    std::cout << "\n";
}

// sort the vector elements by second element of pairs
bool sort_by_value(const std::pair<int, double> &a, const std::pair<int, double> &b) {
    return (a.second < b.second);
}

void add_both(std::map<int, std::set<int>> &graph, int i, int j) {
    graph[i].insert(j);
    graph[j].insert(i);
}

void dfs(int i, std::map<int, std::set<int>> &graph, std::vector<bool> &visited, std::vector<int> &component) {
    visited[i] = true;
    for (auto &j : graph[i]) {
        if (!visited[j]) {
            dfs(j, graph, visited, component);
        }
    }
    component.push_back(i);
}


LinkConstraint transitive_closure(std::vector<std::pair<int, int>> &ml, std::vector<std::pair<int, int>> &cl, int n) {

    std::map<int, std::set<int>> ml_graph;
    std::map<int, std::set<int>> cl_graph;

    for (int i = 0; i < n; i++) {
        ml_graph.insert(std::pair<int, std::set<int>> (i, {}));
        cl_graph.insert(std::pair<int, std::set<int>> (i, {}));
    }

    for (auto &pair_ml : ml) {
        add_both(ml_graph, pair_ml.first, pair_ml.second);
    }

    std::vector<bool> visited(n, false);
    for (int i = 0; i < n; i++) {
        if (!visited[i]) {
            std::vector<int> component;
            dfs(i, ml_graph, visited, component);
            for (auto &x1 : component) {
                for (auto &x2 : component) {
                    if (x1 != x2) {
                        ml_graph[x1].insert(x2);
                    }
                }
            }
        }
    }

    for (auto &pair_cl : cl) {
        int i = pair_cl.first;
        int j = pair_cl.second;
        add_both(cl_graph, i, j);
        for (auto &y : ml_graph[j]) {
            add_both(cl_graph, i, y);
        }
        for (auto &x : ml_graph[i]) {
            add_both(cl_graph, x, j);
            for (auto &y : ml_graph[j]) {
                add_both(cl_graph, x, y);
            }
        }
    }

    for (int i = 0; i < n; i++) {
        for (auto &j : ml_graph[i]) {
            std::set<int> set_j = cl_graph[i];
            if (j != i) {
                if (set_j.count(j)) {
                    std::fprintf(stderr, "Inconsistent constraints between %d and %d", i, j);
                    std::cerr << "ML\n";
                    for (auto &pair : ml){
                        std::cerr << pair.first << " - " << pair.second << "\n";
                    }
                    std::cerr << "CL\n";
                    for (auto &pair : cl){
                        std::cerr << pair.first << " - " << pair.second << "\n";
                    }
                }
            }
        }
    }

    return LinkConstraint{ml_graph, cl_graph};

}

void display_graph(std::map<int, std::set<int>> &map) {

    for (auto &map_elem : map) {
        int key = map_elem.first;
        std::set<int> value = map_elem.second;
        if (value.empty())
            continue;
        std::printf("%d: ", key);
        std::printf("{");
        for (auto &set_elem : value) {
            std::printf(" %d ", set_elem);
        }
        std::printf("}\n");
    }
}


bool comparator_find_branch(std::pair<std::pair<int, int>, double> &a, std::pair<std::pair<int, int>, double> &b) {
    return (a.second < b.second);
}

// branching decision
std::pair<int, int> find_branch(arma::mat &Z) {

    std::vector<std::pair<std::pair<int, int>, double>> branch_data;

    int n = Z.n_rows;
    for (int i = 0; i < n ; i++) {
        for (int j = i + 1; j < n; j++) {
            double min1 = std::min(Z(i, j), Z(i, i) - Z(i, j));
            std::pair<int, int> i_j = std::pair<int, int>(i, j);
            std::pair<std::pair<int, int>, double> data1(i_j, min1);
            double min2 = std::min(Z(i, j), Z(j, j) - Z(i, j));
            std::pair<int, int> j_i = std::pair<int, int>(j, i);
            std::pair<std::pair<int, int>, double> data2(j_i, min2);
            branch_data.push_back(data1);
            branch_data.push_back(data2);
        }
    }

    auto max_elem = std::max_element(branch_data.begin(), branch_data.end(), comparator_find_branch);
    int max_i = max_elem.base()->first.first;
    int max_j = max_elem.base()->first.second;
    double max_val = max_elem.base()->second;

    std::pair<int, int> max_pair(max_i, max_j);

    if (max_val < 1e-5) {
        // Z optimal
        return {-1, -1};
    }

    if (max_j < max_i)
        max_pair = std::pair<int, int>(max_j, max_i);

    return max_pair;
}

// branching decision with norm
std::pair<int, int> find_branch_norm(arma::mat &Z) {

    std::vector<std::pair<std::pair<int, int>, double>> branch_data;

    int n = Z.n_rows;
    for (int i = 0; i < n ; i++) {
        for (int j = i + 1; j < n; j++) {
            double norm = std::pow(arma::norm(Z.row(i) - Z.row(j), 2), 2);
            double first = Z(i, j);
            double min1 = std::min(first, norm);
            std::pair<int, int> i_j = std::pair<int, int>(i, j);
            std::pair<std::pair<int, int>, double> data1(i_j, min1);
            branch_data.push_back(data1);
        }
    }

    auto max_elem = std::max_element(branch_data.begin(), branch_data.end(), comparator_find_branch);
    int max_i = max_elem.base()->first.first;
    int max_j = max_elem.base()->first.second;
    double max_val = max_elem.base()->second;

    std::pair<int, int> max_pair(max_i, max_j);

    if (max_val < 1e-5) {
        // Z optimal
        return {-1, -1};
    }

    if (max_j < max_i)
        max_pair = std::pair<int, int>(max_j, max_i);

    return max_pair;
}

void print_header_sdp(std::ostream &log_file) {

    log_file << "\n" << "|" <<
              std::setw(5) << "N" << "|" <<
              std::setw(5) << "BIN" << "|" <<
              std::setw(8) << "PARENT" << "|" <<
              std::setw(8) << "NODE" << "|" <<
              std::setw(12) << "LB_PAR" << "|" <<
              std::setw(12) << "LB" << "|" <<
              std::setw(10) << "TIME (s)" << "|" <<
              std::setw(8) << "CP_ITER" << "|" <<
              std::setw(8) << "CP_FLAG" << "|" <<
              std::setw(10) << "CP_INEQ" << "|" <<
              std::setw(8) << "SDP_FIX" << "|" <<
              std::setw(9) << "TIME_FIX" << "|" <<
              std::setw(8) << "N_FIX" << "|" <<
              std::setw(12) << "UB" << "|" <<
              std::setw(12) << "GUB" << "|" <<
              std::setw(6) << "I_FIX" << "|" <<
              std::setw(13) << "NODE_GAP" << "|" <<
              std::setw(13) << "GAP" << "|" <<
              std::setw(6) << "OPEN" << "|"
              << std::endl;

}

void print_log_sdp(std::ostream &log_file, int n, int n_bin, int node_parent, int node, double lb_parent, double lb,
                   double time, int cp_iter, int cp_flag, int n_ineq, int n_sdp_fixing, double time_sdp_fixing, int n_fixed,
                   double ub, double gub, int i, double node_gap, double gap, int open, bool update) {

    if (!update) {

        log_file << "|" <<
                  std::setw(5) << n << "|" <<
                  std::setw(5) << n_bin << "|" <<
                  std::setw(8) << node_parent << "|" <<
                  std::setw(8) << node << "|" <<
                  std::setw(12) << lb_parent << "|" <<
                  std::setw(12) << lb << "|" <<
                  std::setw(10) << time << "|" <<
                  std::setw(8) << cp_iter << "|" <<
                  std::setw(8) << cp_flag << "|" <<
                  std::setw(10) << n_ineq << "|" <<
                  std::setw(8) << n_sdp_fixing << "|" <<
                  std::setw(9) << time_sdp_fixing << "|" <<
                  std::setw(8) << n_fixed << "|" <<
                  std::setw(12) << ub << "|" <<
                  std::setw(12) << gub << "|" <<
                  std::setw(6) << i << "|" <<
                  std::setw(13) << node_gap << "|" <<
                  std::setw(13) << gap << "|" <<
                  std::setw(6) << open << "|"
                  << std::endl;

    } else {

        log_file << "|" <<
                  std::setw(5) << n << "|" <<
                  std::setw(5) << n_bin << "|" <<
                  std::setw(8) << node_parent << "|" <<
                  std::setw(8) << node << "|" <<
                  std::setw(12) << lb_parent << "|" <<
                  std::setw(12) << lb << "|" <<
                  std::setw(10) << time << "|" <<
                  std::setw(8) << cp_iter << "|" <<
                  std::setw(8) << cp_flag << "|" <<
                  std::setw(10) << n_ineq << "|" <<
                  std::setw(8) << n_sdp_fixing << "|" <<
                  std::setw(9) << time_sdp_fixing << "|" <<
                  std::setw(8) << n_fixed << "|" <<
                  std::setw(12) << ub << "|" <<
                  std::setw(11) << gub << "*|" <<
                  std::setw(6) << i << "|" <<
                  std::setw(13) << node_gap << "|" <<
                  std::setw(13) << gap << "|" <<
                  std::setw(6) << open << "|"
                  << std::endl;

    }
}