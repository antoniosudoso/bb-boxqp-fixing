#ifndef CLUSTERING_NODE_H
#define CLUSTERING_NODE_H

#include <armadillo>
#include <map>
#include <set>
#include <vector>

class Node {

public:

    // variables to fix to 0 or 1
    std::vector<std::pair<int, int>> x_fix;
    // variables inside the box
    std::vector<int> x_in;

    // lower bound
    double lb;
    // upper bound
    double ub;
    // node id
    int id;

};


class SDPNode : public Node {


public:

    std::vector<std::vector<arma::sp_mat>> B_vector;
    arma::vec l_vector;

};

typedef struct NodeData {

    SDPNode *node;
    int i;
    int value;

} NodeData;

typedef struct JobData {

    int type;
    NodeData *node_data;

} JobData;

typedef struct SDPResult {

    int n;
    int n_bin;
    std::vector<std::pair<int, int>> init_x_fix;
    int n_sdp_fixing;
    double time_sdp_fixing;
    int n_fixed;
    double lb;
    double ub;
    arma::vec x_sdp;
    arma::vec x_heuristic;
    int n_ineq;
    int cp_iter;
    int cp_flag;
    int branching_index;
    int branching_type;
    std::vector<std::vector<arma::sp_mat>> B_vector;
    arma::vec l_vector;

} SDPResult;


#endif //CLUSTERING_NODE_H
