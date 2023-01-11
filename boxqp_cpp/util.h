#ifndef CLUSTERING_UTIL_H
#define CLUSTERING_UTIL_H

#include <map>
#include <set>
#include <vector>
#include <armadillo>
#include "config_params.h"

void save_x_to_file(Config *config, arma::vec &x);
void print_header_sdp(std::stringstream &log_string);
void print_log_sdp(std::stringstream &log_string, int n, int n_bin, int node_parent, int node, double lb_parent, double lb,
                   double time, int cp_iter, int cp_flag, int n_ineq, int n_sdp_fixing, double time_sdp_fixing, int n_fixed,
                   double ub, double gub, int i, double node_gap, double gap, int open, bool update);

#endif //CLUSTERING_UTIL_H
