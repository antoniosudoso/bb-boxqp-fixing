#ifndef CLUSTERING_CONFIG_PARAMS_H
#define CLUSTERING_CONFIG_PARAMS_H

#define BEST_FIRST 0
#define DEPTH_FIRST 1
#define BREADTH_FIRST 2

#define ZERO 0
#define ONE 1
#define IN 2

#define BINARY 2
#define TERNARY 3

extern std::ofstream log_file;

typedef struct Config {

    const char *data_path;
    const char *log_path;
    const char *result_path;

    // branch and bound setting
    double branch_and_bound_tol;
    int branch_and_bound_parallel;
    int branch_and_bound_max_nodes;
    int branch_and_bound_visiting_strategy;
    int branch_and_bound_fixing;
    double branch_and_bound_fixing_tol;

    // matlab setting
    int matlab_session_threads_root;
    int matlab_session_threads_child;
    const char *matlab_source_folder;

    // sdp setting
    double sdp_solver_tol;
    int sdp_solver_verbose;

    // additional solvers
    const char *snopt_folder;
    const char *snopt_license;
    const char *gurobi_folder;
    const char *sdpnal_folder;

    // cutting plane setting
    int cp_max_iter;
    double cp_tol;
    int cp_max_ineq;
    double cp_perc_ineq;
    double cp_eps_ineq;
    double cp_eps_active;

} Config;


#endif //CLUSTERING_CONFIG_PARAMS_H
