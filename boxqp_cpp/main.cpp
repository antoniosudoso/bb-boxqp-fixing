#include <iostream>
#include <map>
#include <algorithm>
#include <armadillo>
#include "sdp_branch_and_bound.h"
#include "config_params.h"

std::ofstream log_file;
std::stringstream log_string;

std::map<std::string, std::string> read_params(std::string &config_file) {

    std::map<std::string, std::string> config_map = {};

    std::ifstream cFile (config_file);
    if (cFile.is_open()) {
        std::string line;
        while (getline(cFile, line)){
            line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
            if(line[0] == '#' || line.empty())
                continue;
            auto delimiterPos = line.find('=');
            auto key = line.substr(0, delimiterPos);
            auto value = line.substr(delimiterPos + 1);
            config_map.insert(std::pair<std::string, std::string>(key, value));
        }

    }
    else {
        std::cerr << "Couldn't open config file for reading.\n";
    }

    return config_map;
}

// read data Q and c
InputData *read_data(const char *filename) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << strerror(errno) << "\n";
        exit(EXIT_FAILURE);
    }
    // read the header n
    int n;
    file >> n;
    arma::vec c(n);
    for (int i = 0; i < n; i++) {
        file >> c(i);
    }
    arma::mat Q(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            file >> Q(i, j);
        }
    }
    auto *input_data = new InputData();
    input_data->Q = Q;
    input_data->c = c;
    return input_data;
}

void run(int argc, char **argv) {

    if (argc < 4) {
        std::cerr << "Usage: <config.txt> <DATA_PATH> <LOG_PATH> <RESULT_PATH>" << std::endl;
        exit(EXIT_FAILURE);
    }

    auto *config = new Config();

    std::string config_file = argv[1];
    config->data_path = argv[2];
    config->log_path = argv[3];
    log_file.open(config->log_path);
    config->result_path = argv[4];

    std::map<std::string, std::string> config_map = read_params(config_file);

    InputData *input_data = read_data(config->data_path);

    // branch and bound setting
    config->branch_and_bound_tol = std::stod(config_map["BRANCH_AND_BOUND_TOL"]);
    config->branch_and_bound_parallel = std::stoi(config_map["BRANCH_AND_BOUND_PARALLEL"]);
    config->branch_and_bound_max_nodes = std::stoi(config_map["BRANCH_AND_BOUND_MAX_NODES"]);
    config->branch_and_bound_visiting_strategy = std::stoi(config_map["BRANCH_AND_BOUND_VISITING_STRATEGY"]);
    config->branch_and_bound_fixing = std::stoi(config_map["BRANCH_AND_BOUND_FIXING"]);
    config->branch_and_bound_fixing_tol = std::stod(config_map["BRANCH_AND_BOUND_FIXING_TOL"]);

    // matlab
    config->matlab_session_threads_root = std::stoi(config_map["MATLAB_SESSION_THREADS_ROOT"]);
    config->matlab_session_threads_child = std::stoi(config_map["MATLAB_SESSION_THREADS_CHILD"]);
    config->matlab_source_folder = config_map["MATLAB_SOURCE_FOLDER"].c_str();

    // sdp solver
    config->sdp_solver_tol = std::stod(config_map["SDP_SOLVER_TOL"]);
    config->sdp_solver_verbose = std::stoi(config_map["SDP_SOLVER_VERBOSE"]);

    // cutting plane
    config->cp_max_iter = std::stoi(config_map["CP_MAX_ITER"]);
    config->cp_tol = std::stod(config_map["CP_TOL"]);
    config->cp_max_ineq = std::stoi(config_map["CP_MAX_INEQ"]);
    config->cp_perc_ineq = std::stod(config_map["CP_PERC_INEQ"]);
    config->cp_eps_ineq = std::stod(config_map["CP_EPS_INEQ"]);
    config->cp_eps_active = std::stod(config_map["CP_EPS_ACTIVE"]);

    // additional solvers
    config->gurobi_folder = config_map["GUROBI_FOLDER"].c_str();
    config->snopt_folder = config_map["SNOPT_FOLDER"].c_str();
    config->snopt_license = config_map["SNOPT_LICENSE"].c_str();
    config->sdpnal_folder = config_map["SDPNAL_FOLDER"].c_str();

    log_file << "\n" << "DATA_PATH: " << config->data_path <<"\n";
    log_file << "LOG_PATH: " << config->log_path << "\n\n";

    log_file << "BRANCH_AND_BOUND_TOL: " << config->branch_and_bound_tol << "\n";
    log_file << "BRANCH_AND_BOUND_PARALLEL: " << config->branch_and_bound_parallel << "\n";
    log_file << "BRANCH_AND_BOUND_MAX_NODES: " <<  config->branch_and_bound_max_nodes << "\n";
    log_file << "BRANCH_AND_BOUND_VISITING_STRATEGY: " << config->branch_and_bound_visiting_strategy << "\n";
    log_file << "BRANCH_AND_BOUND_FIXING: " << config->branch_and_bound_fixing << "\n";
    log_file << "BRANCH_AND_BOUND_FIXING_TOL: " << config->branch_and_bound_fixing_tol << "\n\n";

    log_file << "MATLAB_SESSION_THREADS_ROOT: " << config->matlab_session_threads_root << "\n";
    log_file << "MATLAB_SESSION_THREADS_CHILD: " << config->matlab_session_threads_child << "\n";
    log_file << "MATLAB_SOURCE_FOLDER: " << config->matlab_source_folder << "\n\n";

    log_file << "SDP_SOLVER_TOL: " << config->sdp_solver_tol << "\n";
    log_file << "SDP_SOLVER_VERBOSE: " << config->sdp_solver_verbose << "\n\n";

    log_file << "CP_MAX_ITER: " << config->cp_max_iter << "\n";
    log_file << "CP_TOL: " << config->cp_tol << "\n";
    log_file << "CP_MAX_INEQ: " << config->cp_max_ineq << "\n";
    log_file << "CP_PERC_INEQ: " << config->cp_perc_ineq << "\n";
    log_file << "CP_EPS_INEQ: " << config->cp_eps_ineq << "\n";
    log_file << "CP_EPS_ACTIVE: " << config->cp_eps_active << "\n\n";

    log_file << "GUROBI_FOLDER: " << config->gurobi_folder << "\n";
    log_file << "SNOPT_FOLDER: " << config->snopt_folder << "\n";
    log_file << "SNOPT_LICENSE: " << config->snopt_license << "\n";
    log_file << "SDPNAL_FOLDER: " << config->sdpnal_folder << "\n";

    log_file.flush();
    sdp_branch_and_bound(config, input_data);
    log_file.close();

    delete (config);

}

int main(int argc, char **argv) {

    // The B&B algorithm for solving the nonconvex qudratic programming problem:
    //
    //    min      1/2*x'*Q*x + c'*x
    //    s.t.     0 <=  x <= 1
    //
    //    Q, c: data of the QP formulation above

    run(argc, argv);

    return EXIT_SUCCESS;
}
